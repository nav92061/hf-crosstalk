#!/usr/bin/env python
"""External accuracy check (pre-specified plan PRESPEC_external_accuracy.md, 2026-10-10).

For each LCTSC case: download CT + RTSTRUCT from TCIA (NBIA API v1), run
TotalSegmentator 2.18.0 (task 'total', default non-fast 1.5 mm models, multilabel
output; heart = label 51) on the DICOM directory as in the main study, rasterise
the RTSTRUCT ROI named 'Heart' (case-insensitive exact match) onto the CT grid,
and compute Dice, HD95 and volumes. One CSV row is appended per case immediately
after it finishes; images are deleted after each case.

Usage: python lctsc_external_accuracy_job.py <case_list.csv> <out_csv> [--limit N] [--done prev.csv]
"""
import argparse, csv, io, os, shutil, subprocess, sys, time, traceback, zipfile
from pathlib import Path

import numpy as np
import pydicom
import requests
import SimpleITK as sitk
from scipy import ndimage
from skimage.draw import polygon as sk_polygon

NBIA = "https://services.cancerimagingarchive.net/nbia-api/services/v1/getImage"
HEART_LABEL = 51  # class_map['total'][51] == 'heart' in TotalSegmentator 2.18.0
FIELDS = ["PatientID", "status", "reason", "ct_series_uid", "rtstruct_series_uid",
          "roi_name_matched", "n_slices", "spacing_x_mm", "spacing_y_mm", "slice_spacing_mm",
          "n_contour_slices", "contour_z_max_offset_mm", "for_uid_match", "ref_series_match",
          "dice", "hd95_mm", "hd95_ts_to_ref_mm", "hd95_ref_to_ts_mm",
          "vol_ref_ml", "vol_ts_ml", "vol_diff_ml", "vol_diff_pct",
          "ref_sup_z_mm", "ref_inf_z_mm", "ts_sup_z_mm", "ts_inf_z_mm",
          "ts_runtime_s", "ts_exit_code"]


def download_series(uid, dest: Path, retries=3):
    dest.mkdir(parents=True, exist_ok=True)
    for k in range(retries):
        try:
            r = requests.get(NBIA, params={"SeriesInstanceUID": uid}, timeout=600)
            r.raise_for_status()
            with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                for n in z.namelist():
                    if n.lower().endswith(".dcm"):
                        (dest / Path(n).name).write_bytes(z.read(n))
            files = list(dest.glob("*.dcm"))
            if not files:
                raise RuntimeError("no DICOM files in archive")
            return files
        except Exception:
            if k == retries - 1:
                raise
            time.sleep(10 * (k + 1))


def read_ct(ct_dir: Path):
    reader = sitk.ImageSeriesReader()
    sids = reader.GetGDCMSeriesIDs(str(ct_dir))
    if len(sids) != 1:
        raise RuntimeError(f"expected 1 CT series in dir, found {len(sids)}")
    files = reader.GetGDCMSeriesFileNames(str(ct_dir), sids[0])
    reader.SetFileNames(files)
    img = reader.Execute()
    ds0 = pydicom.dcmread(files[0], stop_before_pixels=True)
    return img, ds0


def rasterise_heart(rs_path: Path, ct: sitk.Image, ct_ds0):
    rs = pydicom.dcmread(str(rs_path))
    rois = {int(r.ROINumber): str(r.ROIName) for r in rs.StructureSetROISequence}
    matches = [n for n, name in rois.items() if name.strip().lower() == "heart"]
    info = {"roi_name_matched": "", "for_uid_match": "", "ref_series_match": "",
            "n_contour_slices": 0, "contour_z_max_offset_mm": np.nan}
    if len(matches) != 1:
        return None, info, f"Heart ROI matches={len(matches)}; ROIs={sorted(rois.values())}"
    roi_no = matches[0]
    info["roi_name_matched"] = rois[roi_no]
    # alignment checks: frame of reference and referenced CT series
    for_uids = {str(r.FrameOfReferenceUID) for r in rs.ReferencedFrameOfReferenceSequence}
    info["for_uid_match"] = str(ct_ds0.FrameOfReferenceUID) in for_uids
    ref_series = set()
    for fr in rs.ReferencedFrameOfReferenceSequence:
        for st in getattr(fr, "RTReferencedStudySequence", []):
            for se in getattr(st, "RTReferencedSeriesSequence", []):
                ref_series.add(str(se.SeriesInstanceUID))
    info["ref_series_match"] = (str(ct_ds0.SeriesInstanceUID) in ref_series) if ref_series else "not_listed"
    contours = None
    for rc in rs.ROIContourSequence:
        if int(rc.ReferencedROINumber) == roi_no:
            contours = getattr(rc, "ContourSequence", None)
    if not contours:
        return None, info, "Heart ROI has no contour data"
    size = ct.GetSize()  # (x, y, z)
    mask = np.zeros((size[2], size[1], size[0]), dtype=bool)  # z, y, x
    zs, offsets = set(), []
    for c in contours:
        pts = np.asarray(c.ContourData, dtype=float).reshape(-1, 3)
        idx = np.array([ct.TransformPhysicalPointToContinuousIndex(tuple(p)) for p in pts])
        zi = idx[:, 2]
        k = int(round(float(np.mean(zi))))
        offsets.append(float(np.max(np.abs(zi - k))) * ct.GetSpacing()[2])
        if k < 0 or k >= size[2]:
            return None, info, f"contour slice index {k} outside CT (0..{size[2]-1})"
        rr, cc = sk_polygon(idx[:, 1], idx[:, 0], shape=(size[1], size[0]))
        sl = np.zeros((size[1], size[0]), dtype=bool)
        sl[rr, cc] = True
        mask[k] ^= sl  # even-odd rule handles holes / multiple contours
        zs.add(k)
    info["n_contour_slices"] = len(zs)
    info["contour_z_max_offset_mm"] = round(max(offsets), 3)
    if not info["for_uid_match"]:
        return None, info, "FrameOfReferenceUID mismatch between RTSTRUCT and CT"
    if max(offsets) > 0.5 * ct.GetSpacing()[2]:
        return None, info, f"contour points off CT slice planes by up to {max(offsets):.2f} mm"
    return mask, info, ""


def surface_distances(a, b, spacing_zyx):
    """Directed distances (mm) from surface voxels of a to the surface of b."""
    sa = a & ~ndimage.binary_erosion(a, structure=ndimage.generate_binary_structure(3, 1))
    sb = b & ~ndimage.binary_erosion(b, structure=ndimage.generate_binary_structure(3, 1))
    dt_b = ndimage.distance_transform_edt(~sb, sampling=spacing_zyx)
    return dt_b[sa]


def z_extent_mm(mask, ct):
    zs = np.where(mask.any(axis=(1, 2)))[0]
    if len(zs) == 0:
        return np.nan, np.nan
    zw = [ct.TransformIndexToPhysicalPoint((0, 0, int(k)))[2] for k in (zs.min(), zs.max())]
    return round(max(zw), 2), round(min(zw), 2)


def process(row, work: Path):
    pid = row["PatientID"]
    out = {f: "" for f in FIELDS}
    out.update(PatientID=pid, ct_series_uid=row["CT_SeriesUID"], rtstruct_series_uid=row["RTSTRUCT_SeriesUID"])
    cdir, rdir, tdir = work / "ct", work / "rs", work / "ts"
    for d in (cdir, rdir, tdir):
        shutil.rmtree(d, ignore_errors=True)
    try:
        download_series(row["CT_SeriesUID"], cdir)
        rs_files = download_series(row["RTSTRUCT_SeriesUID"], rdir)
    except Exception as e:
        out.update(status="failed", reason=f"download error: {e!r}"[:300]); return out
    try:
        ct, ds0 = read_ct(cdir)
        sp = ct.GetSpacing()
        out.update(n_slices=ct.GetSize()[2], spacing_x_mm=round(sp[0], 4), spacing_y_mm=round(sp[1], 4),
                   slice_spacing_mm=round(sp[2], 4))
        ref, info, why = rasterise_heart(rs_files[0], ct, ds0)
        out.update(info)
        if ref is None:
            out.update(status="excluded_reference", reason=why); return out
        # TotalSegmentator on the DICOM directory (as in main study; default, not --fast)
        seg_path = tdir / "seg.nii.gz"
        tdir.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        p = subprocess.run(["TotalSegmentator", "-i", str(cdir), "-o", str(seg_path), "--ml",
                            "--task", "total"], capture_output=True, text=True)
        out["ts_runtime_s"] = round(time.time() - t0, 1)
        out["ts_exit_code"] = p.returncode
        if p.returncode != 0 or not seg_path.exists():
            out.update(status="failed", reason=("TotalSegmentator error: " + p.stderr[-300:]).replace("\n", " "))
            return out
        seg = sitk.ReadImage(str(seg_path))
        # bring TS labels onto the exact CT grid (physical-space, nearest neighbour)
        seg_on_ct = sitk.Resample(seg, ct, sitk.Transform(), sitk.sitkNearestNeighbor, 0, seg.GetPixelID())
        ts = sitk.GetArrayFromImage(seg_on_ct) == HEART_LABEL
        vox_ml = sp[0] * sp[1] * sp[2] / 1000.0
        vr, vt = ref.sum() * vox_ml, ts.sum() * vox_ml
        out.update(vol_ref_ml=round(vr, 2), vol_ts_ml=round(vt, 2), vol_diff_ml=round(vt - vr, 2),
                   vol_diff_pct=round(100 * (vt - vr) / vr, 2) if vr > 0 else "")
        out["ref_sup_z_mm"], out["ref_inf_z_mm"] = z_extent_mm(ref, ct)
        out["ts_sup_z_mm"], out["ts_inf_z_mm"] = z_extent_mm(ts, ct)
        if ts.sum() == 0:
            out.update(status="failed", reason="TotalSegmentator heart mask empty", dice=0.0); return out
        inter = np.logical_and(ref, ts).sum()
        out["dice"] = round(2 * inter / (ref.sum() + ts.sum()), 4)
        spz = (sp[2], sp[1], sp[0])
        d_tr = surface_distances(ts, ref, spz)
        d_rt = surface_distances(ref, ts, spz)
        out["hd95_mm"] = round(float(np.percentile(np.concatenate([d_tr, d_rt]), 95)), 2)
        out["hd95_ts_to_ref_mm"] = round(float(np.percentile(d_tr, 95)), 2)
        out["hd95_ref_to_ts_mm"] = round(float(np.percentile(d_rt, 95)), 2)
        out["status"] = "ok"
        return out
    except Exception as e:
        out.update(status="failed", reason=("exception: " + repr(e) + " " + traceback.format_exc()[-200:]).replace("\n", " ")[:400])
        return out
    finally:
        for d in (cdir, rdir, tdir):
            shutil.rmtree(d, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_list"); ap.add_argument("out_csv")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--done", default="")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.case_list)))
    if a.limit:
        rows = rows[:a.limit]
    out_csv = Path(a.out_csv); out_csv.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if a.done and Path(a.done).exists():
        prev = list(csv.DictReader(open(a.done)))
        done = {r["PatientID"] for r in prev}
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, FIELDS); w.writeheader(); [w.writerow({k: r.get(k, "") for k in FIELDS}) for r in prev]
    elif not out_csv.exists():
        with open(out_csv, "w", newline="") as f:
            csv.DictWriter(f, FIELDS).writeheader()
    work = Path("/tmp/lctsc_work")
    for i, row in enumerate(rows):
        if row["PatientID"] in done:
            continue
        t0 = time.time()
        res = process(row, work)
        with open(out_csv, "a", newline="") as f:
            csv.DictWriter(f, FIELDS).writerow(res)
            f.flush(); os.fsync(f.fileno())
        print(f"[{i+1}/{len(rows)}] {row['PatientID']} {res['status']} dice={res['dice']} "
              f"hd95={res['hd95_mm']} {time.time()-t0:.0f}s {res['reason'][:120]}", flush=True)


if __name__ == "__main__":
    main()

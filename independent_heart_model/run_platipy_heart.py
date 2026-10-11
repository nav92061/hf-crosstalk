#!/usr/bin/env python3
"""PlatiPy open whole-heart model (nnU-Net Task400_OPEN_HEART_1FOLD, 3d_lowres) on TCIA series.
One series at a time: download DICOM, convert, segment with PlatiPy's default whole-heart
settings, compute heart volume on the native grid, delete the images. Writes a row per series
to out/volumes_platipy_heart.csv after each series. Plan: PRESPEC_independent_heart_model.md.

Usage: python run_platipy_heart.py series.csv geometry.csv [max_series]
"""
import sys, os, io, json, time, shutil, hashlib, zipfile, traceback, platform
os.environ.setdefault("RESULTS_FOLDER", "/opt/nnunet_models")
os.environ.setdefault("nnUNet_raw_data_base", "/tmp/nnraw")
os.environ.setdefault("nnUNet_preprocessed", "/tmp/nnpre")
import numpy as np, pandas as pd, requests, SimpleITK as sitk, torch
import importlib.metadata as md
from platipy.imaging.projects.nnunet.run import run_segmentation
from platipy.imaging.projects.cardiac.run import HYBRID_SETTINGS_DEFAULTS

BASE = "https://services.cancerimagingarchive.net/nbia-api/services/v1"
ZIP = "/opt/nnunet_zip/Task400.zip"; ZIP_MD5 = "feb63660bce90c2b6ecf561c8463b6a9"
OUT = os.environ.get("PLATIPY_OUT", "out"); os.makedirs(OUT, exist_ok=True)
ser = pd.read_csv(sys.argv[1]).series_uid.tolist()
geo = pd.read_csv(sys.argv[2]).set_index("series_uid")
if len(sys.argv) > 3: ser = ser[: int(sys.argv[3])]
res_path = f"{OUT}/volumes_platipy_heart.csv"
SETTINGS = dict(HYBRID_SETTINGS_DEFAULTS["nnunet_settings"])  # PlatiPy's own whole-heart settings, unmodified

# ---- provenance and check 1: weights MD5
md5 = hashlib.md5(open(ZIP, "rb").read()).hexdigest()
ckpts = sorted(os.path.join(dp, f) for dp, _, fs in os.walk(os.environ["RESULTS_FOLDER"]) for f in fs if f.endswith(".model"))
prov = {"weights_zip_md5": md5, "weights_zip_md5_expected": ZIP_MD5,
        "checkpoint_files": {os.path.relpath(p, os.environ["RESULTS_FOLDER"]): hashlib.md5(open(p, "rb").read()).hexdigest() for p in ckpts},
        "nnunet_settings": SETTINGS,
        "platipy": md.version("platipy"), "nnunet": md.version("nnunet"),
        "torch": torch.__version__, "cuda": torch.version.cuda,
        "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "python": platform.python_version(), "sitk": sitk.Version_VersionString(),
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
json.dump(prov, open(f"{OUT}/provenance.json", "w"), indent=1, default=str)
assert md5 == ZIP_MD5, "weights MD5 mismatch - halting"
assert ckpts, "no installed nnU-Net checkpoints - halting"
assert torch.cuda.is_available(), "no GPU visible - halting"

S = requests.Session()
def fetch_series(uid, dest):
    for attempt in range(4):
        try:
            r = S.get(f"{BASE}/getImage", params={"SeriesInstanceUID": uid}, timeout=600)
            r.raise_for_status()
            os.makedirs(dest, exist_ok=True)
            zipfile.ZipFile(io.BytesIO(r.content)).extractall(dest)
            return attempt
        except Exception:
            if attempt == 3: raise
            time.sleep(10 * (attempt + 1))

def read_series(dcm_dir):
    rd = sitk.ImageSeriesReader()
    ids = rd.GetGDCMSeriesIDs(dcm_dir)
    files = max((rd.GetGDCMSeriesFileNames(dcm_dir, i) for i in ids), key=len)
    rd.SetFileNames(files)
    return rd.Execute(), len(files)

done = set(pd.read_csv(res_path).series_uid) if os.path.exists(res_path) else set()
log = open(f"{OUT}/run_log.txt", "a")
for k, uid in enumerate(ser, 1):
    if uid in done: continue
    row = {"series_uid": uid, "error": None}
    t0 = time.time(); work = f"/tmp/s_{k}"
    try:
        row["download_retries"] = fetch_series(uid, f"{work}/dcm")
        img, nfiles = read_series(f"{work}/dcm")
        sp = img.GetSpacing(); sz = img.GetSize(); dirm = np.array(img.GetDirection()).reshape(3, 3)
        row.update(n_files=nfiles, nx=sz[0], ny=sz[1], nz=sz[2], sx=sp[0], sy=sp[1], sz_mm=sp[2])
        # check 3: extent along patient z vs ImagePositionPatient extent
        corners = [img.TransformIndexToPhysicalPoint((i, j, kk)) for i in (0, sz[0]-1) for j in (0, sz[1]-1) for kk in (0, sz[2]-1)]
        zs = [c[2] for c in corners]
        ext = max(zs) - min(zs) + abs(dirm[2] @ np.array(sp))
        ref = geo.extent_ipp.get(uid, np.nan)
        row["extent_image"] = ext; row["extent_ipp"] = ref
        tol = max(abs(geo.voxel_z.get(uid, sp[2])), 0.5)
        row["extent_check_pass"] = bool(np.isfinite(ref) and abs(ext - ref) <= tol)
        t1 = time.time(); res = run_segmentation(img, SETTINGS); row["infer_s"] = round(time.time() - t1, 1)
        row["n_labels"] = len(res)                                           # check 4
        assert len(res) == 1, f"expected one label, got {sorted(res)}"
        m = sitk.GetArrayFromImage(res["Struct_0"]).astype(bool)            # z,y,x
        assert m.shape == (sz[2], sz[1], sz[0]), f"mask grid {m.shape} != image {sz[::-1]}"
        n = int(m.sum()); vv = float(np.prod(sp))
        row["heart_nvox"] = n; row["heart_vol_mm3"] = n * vv
        row["heart_empty"] = n == 0
        row["heart_edge"] = bool(n and (m[0].any() or m[-1].any()))        # check 5
    except Exception as e:
        row["error"] = f"{type(e).__name__}: {str(e)[:300]}"
        traceback.print_exc(file=log)
    finally:
        shutil.rmtree(work, ignore_errors=True)
        torch.cuda.empty_cache()
    row["total_s"] = round(time.time() - t0, 1)
    pd.DataFrame([row]).to_csv(res_path, mode="a", header=not os.path.exists(res_path), index=False)
    print(f"[{k}/{len(ser)}] {uid[-12:]} {row['total_s']}s err={row['error']}", flush=True)
    log.write(json.dumps({"k": k, "uid": uid, "error": row["error"], "total_s": row["total_s"]}) + "\n"); log.flush()

prov["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
json.dump(prov, open(f"{OUT}/provenance.json", "w"), indent=1, default=str)
print("DONE")

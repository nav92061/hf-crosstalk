"""Phase 0: pull DICOM headers for one representative slice per series.

Uses TCIA getSOPInstanceUIDs + getSingleImage. Never downloads a full series.
"""
import urllib.request, urllib.parse, json, io, os, time
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import pydicom

BASE = "https://services.cancerimagingarchive.net/nbia-api/services/v1"

TAGS = [
    "Manufacturer", "ManufacturerModelName", "SoftwareVersions", "StationName",
    "KVP", "XRayTubeCurrent", "Exposure", "ExposureTime", "ExposureModulationType",
    "SliceThickness", "SpacingBetweenSlices", "PixelSpacing",
    "ConvolutionKernel", "ReconstructionDiameter", "DataCollectionDiameter",
    "ContrastBolusAgent", "ContrastBolusRoute", "ContrastBolusVolume",
    "StudyDate", "SeriesDate", "AcquisitionDate", "ImageType",
    "Rows", "Columns", "PatientPosition", "BodyPartExamined",
    "SeriesDescription", "ProtocolName", "StudyDescription",
    "AcquisitionNumber", "SliceLocation", "RescaleSlope", "RescaleIntercept",
    "FilterType", "SingleCollimationWidth", "TotalCollimationWidth",
    "SpiralPitchFactor", "CTDIvol", "GantryDetectorTilt", "TableHeight",
    "PatientSex", "PatientAge",
]


def _get(url, timeout=300):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as f:
                return f.read()
        except Exception as e:
            if attempt == 3:
                raise
            time.sleep(2 * (attempt + 1))


def header_for_series(series_uid):
    """Fetch a middle-ish slice and return its header as a dict."""
    b = _get(f"{BASE}/getSOPInstanceUIDs?" + urllib.parse.urlencode(
        {"SeriesInstanceUID": series_uid}))
    sops = [d["SOPInstanceUID"] for d in json.loads(b)]
    if not sops:
        return {"series_uid": series_uid, "error": "no_sops"}
    sop = sops[len(sops) // 2]  # middle slice: avoids scout/first-slice oddities
    raw = _get(f"{BASE}/getSingleImage?" + urllib.parse.urlencode(
        {"SeriesInstanceUID": series_uid, "SOPInstanceUID": sop}))
    ds = pydicom.dcmread(io.BytesIO(raw), stop_before_pixels=True, force=True)
    out = {"series_uid": series_uid, "n_sops": len(sops),
           "probe_sop": sop, "bytes": len(raw), "error": ""}
    for t in TAGS:
        v = getattr(ds, t, None)
        if v is None:
            out[t] = None
        elif isinstance(v, pydicom.multival.MultiValue):
            out[t] = "\\".join(str(x) for x in v)
        else:
            out[t] = str(v)
    return out


def run(series_uids, out_csv, workers=6):
    done = {}
    if os.path.exists(out_csv):
        prev = pd.read_csv(out_csv)
        done = {r["series_uid"]: r.to_dict() for _, r in prev.iterrows()}
    todo = [u for u in series_uids if u not in done]
    print(f"{len(done)} cached, {len(todo)} to fetch", flush=True)
    results = list(done.values())
    n = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(header_for_series, u): u for u in todo}
        for fut in as_completed(futs):
            u = futs[fut]
            try:
                results.append(fut.result())
            except Exception as e:
                results.append({"series_uid": u, "error": f"{type(e).__name__}: {e}"})
            n += 1
            if n % 50 == 0:
                print(f"  {n}/{len(todo)}", flush=True)
                pd.DataFrame(results).to_csv(out_csv, index=False)
    df = pd.DataFrame(results)
    df.to_csv(out_csv, index=False)
    print(f"wrote {out_csv}: {len(df)} rows, "
          f"{(df.error.astype(str).str.len() > 0).sum()} errors", flush=True)
    return df


if __name__ == "__main__":
    coh = pd.read_csv("cohort_series.csv")
    run(coh.series_uid.tolist(), "dicom_headers_raw.csv")

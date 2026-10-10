
"""Fetch ImagePositionPatient for every slice of each series by streaming only the
DICOM header (first ~16-48 KB) of each instance from the TCIA NBIA API."""
import requests, io, os, sys, time, json, pydicom, pandas as pd
from concurrent.futures import ThreadPoolExecutor
BASE="https://services.cancerimagingarchive.net/nbia-api/services/v1"
S=requests.Session()
def head(series, sop):
    for att in range(4):
        try:
            with S.get(f"{BASE}/getSingleImage", params={"SeriesInstanceUID":series,"SOPInstanceUID":sop}, stream=True, timeout=90) as r:
                r.raise_for_status(); buf=b""
                for ch in r.iter_content(8192):
                    buf+=ch
                    if len(buf)>=16384:
                        try:
                            ds=pydicom.dcmread(io.BytesIO(buf), stop_before_pixels=True, force=True)
                            if "ImagePositionPatient" in ds: break
                        except Exception: pass
                    if len(buf)>65536: break
            ds=pydicom.dcmread(io.BytesIO(buf), stop_before_pixels=True, force=True)
            ipp=ds.get("ImagePositionPatient"); iop=ds.get("ImageOrientationPatient")
            it=ds.get("ImageType")
            return dict(series_uid=series, sop=sop, inst=ds.get("InstanceNumber"),
                        x=float(ipp[0]) if ipp else None, y=float(ipp[1]) if ipp else None, z=float(ipp[2]) if ipp else None,
                        iop="\\".join(str(v) for v in iop) if iop else None, thick=ds.get("SliceThickness"),
                        image_type="\\".join(it) if it else None, err="")
        except Exception as e:
            if att==3: return dict(series_uid=series, sop=sop, err=f"{type(e).__name__}: {e}")
            time.sleep(2*(att+1))
def series_rows(u):
    r=S.get(f"{BASE}/getSOPInstanceUIDs", params={"SeriesInstanceUID":u}, timeout=90); sops=[d["SOPInstanceUID"] for d in r.json()]
    with ThreadPoolExecutor(8) as ex: return list(ex.map(lambda s: head(u,s), sops))
uids=pd.read_csv("series.csv").series_uid.tolist()
os.makedirs("per_series",exist_ok=True)
t0=time.time()
for i,u in enumerate(uids):
    f=f"per_series/{u}.csv"
    if os.path.exists(f): continue
    rows=series_rows(u); pd.DataFrame(rows).to_csv(f+".tmp",index=False); os.replace(f+".tmp",f)
    if i%10==0: print(f"{i+1}/{len(uids)} series, {time.time()-t0:.0f}s", flush=True)
print("DONE", time.time()-t0, flush=True)

import urllib.request, urllib.parse, json, io, time
import pandas as pd
from pydicom import dcmread
BASE="https://services.cancerimagingarchive.net/nbia-api/services/v1"
def api(ep,**p):
    u=f"{BASE}/{ep}?"+urllib.parse.urlencode(p)
    with urllib.request.urlopen(u,timeout=300) as f: return json.loads(f.read())
def mid_slice(uid):
    s=api("getSOPInstanceUIDs",SeriesInstanceUID=uid)
    ids=[x["SOPInstanceUID"] for x in s]
    if not ids: return None
    sop=ids[len(ids)//2]
    u=f"{BASE}/getSingleImage?"+urllib.parse.urlencode({"SeriesInstanceUID":uid,"SOPInstanceUID":sop})
    with urllib.request.urlopen(u,timeout=300) as f: raw=f.read()
    return dcmread(io.BytesIO(raw), stop_before_pixels=True, force=True)

pp=pd.read_csv("recon_pairs_detail.csv")
uids=sorted(set(pp.ua)|set(pp.ub))
out=[]
for i,u in enumerate(uids,1):
    try:
        d=mid_slice(u)
        out.append(dict(series_uid=u,
            FrameOfReferenceUID=getattr(d,"FrameOfReferenceUID",None),
            AcquisitionTime=getattr(d,"AcquisitionTime",None),
            ContentTime=getattr(d,"ContentTime",None),
            SeriesTime=getattr(d,"SeriesTime",None),
            ImagePositionPatient=str(getattr(d,"ImagePositionPatient",None)),
            err=""))
    except Exception as e:
        out.append(dict(series_uid=u,err=str(e)[:80]))
    if i%40==0: print(f"{i}/{len(uids)}",flush=True)
pd.DataFrame(out).to_csv("recon_forid.csv",index=False)
print("done",len(out))


import io, os, sys, time, json, requests, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
from pydicom import dcmread
from scipy import ndimage
BASE='https://services.cancerimagingarchive.net/nbia-api/services/v1/'
OUT='sharpness_raw.csv'; LOG='progress.log'
ser=pd.read_csv('series87.csv').series_uid.tolist()
G=pd.read_csv('series_geometry_ipp.csv').set_index('series_uid')
done=set(pd.read_csv(OUT).series_uid) if os.path.exists(OUT) else set()
S=requests.Session()
def get(params, path, stream=False, nbytes=None):
    for k in range(5):
        try:
            r=S.get(BASE+path,params=params,timeout=(20,60),stream=stream)
            if r.status_code==200:
                if nbytes:
                    b=r.raw.read(nbytes); r.close(); return b
                return r.content if path=='getSingleImage' else r.json()
        except Exception as e: err=e
        time.sleep(2*(k+1))
    raise RuntimeError(f'failed {path} {params}')
def hdr(s, sop):
    for nb in (16384, 65536, None):
        b=get({'SeriesInstanceUID':s,'SOPInstanceUID':sop},'getSingleImage',stream=True,nbytes=nb) if nb else get({'SeriesInstanceUID':s,'SOPInstanceUID':sop},'getSingleImage')
        try:
            d=dcmread(io.BytesIO(b),stop_before_pixels=True,force=True,specific_tags=['ImagePositionPatient','ImageOrientationPatient'])
            ipp=np.array(d.ImagePositionPatient,float); iop=np.array(d.ImageOrientationPatient,float)
            return sop, ipp, iop
        except Exception: continue
    raise RuntimeError('hdr '+sop)
def disk(r):
    y,x=np.ogrid[-r:r+1,-r:r+1]; return x*x+y*y<=r*r
SE=disk(10)
def measure(b):
    d=dcmread(io.BytesIO(b),force=True)
    hu=d.pixel_array.astype(np.float64)*float(getattr(d,'RescaleSlope',1))+float(getattr(d,'RescaleIntercept',0))
    psr,psc=[float(v) for v in d.PixelSpacing]
    m=hu>-500
    lab,n=ndimage.label(m)
    if n==0: return None
    sizes=ndimage.sum(m,lab,range(1,n+1)); body=lab==(np.argmax(sizes)+1)
    er=ndimage.binary_erosion(body,structure=SE)
    cy,cx=ndimage.center_of_mass(er)
    nr=int(round(128/psr)); nc=int(round(128/psc))
    r0=int(round(cy-nr/2)); c0=int(round(cx-nc/2))
    r0=min(max(r0,0),hu.shape[0]-nr); c0=min(max(c0,0),hu.shape[1]-nc)
    roi=hu[r0:r0+nr,c0:c0+nc]
    roi=roi-roi.mean()
    win=np.outer(np.hanning(nr),np.hanning(nc))
    F=np.abs(np.fft.fft2(roi*win))**2
    fr=np.fft.fftfreq(nr,d=psr); fc=np.fft.fftfreq(nc,d=psc)
    R=np.sqrt(fr[:,None]**2+fc[None,:]**2)
    F[0,0]=0
    hf=F[R>0.4].sum()/F.sum()
    resid=hu-ndimage.median_filter(hu,size=5)
    nm=er&(hu>=-100)&(hu<=200)
    noise=float(np.std(resid[nm])) if nm.sum()>100 else np.nan
    roi_in_body=float(er[r0:r0+nr,c0:c0+nc].mean())
    return dict(hf=hf,noise=noise,psr=psr,psc=psc,nyq=0.5/max(psr,psc),roi_body_frac=roi_in_body,
                n_noise_px=int(nm.sum()),kernel=str(getattr(d,'ConvolutionKernel','')),rows=hu.shape[0])
ex=ThreadPoolExecutor(8)
t0=time.time()
for i,s in enumerate(ser):
    if s in done: continue
    sops=[x['SOPInstanceUID'] for x in get({'SeriesInstanceUID':s},'getSOPInstanceUIDs')]
    futs=[ex.submit(hdr,s,q) for q in sops]; H=[f.result(timeout=300) for f in futs]
    nrm=np.cross(H[0][2][:3],H[0][2][3:])
    pos=np.array([h[1]@nrm for h in H])
    order=np.argsort(pos); N=len(order)
    rows=[]
    for p in (0.40,0.45,0.50,0.55,0.60):
        k=int(round(p*(N-1))); sop=H[order[k]][0]
        img=get({'SeriesInstanceUID':s,'SOPInstanceUID':sop},'getSingleImage')
        m=measure(img); del img
        m.update(series_uid=s,pct=p,slice_rank=k,n_inst=N,pos=float(pos[order[k]]),sop=sop,
                 normal='/'.join(f'{v:.2f}' for v in nrm),pos_min=float(pos.min()),pos_max=float(pos.max()))
        rows.append(m)
    pd.DataFrame(rows).to_csv(OUT,mode='a',header=not os.path.exists(OUT),index=False)
    with open(LOG,'a') as f: f.write(f'{i+1}/{len(ser)} {s[-12:]} N={N} t={time.time()-t0:.0f}s\n')
print('done', time.time()-t0)

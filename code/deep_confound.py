# -*- coding: utf-8 -*-
"""Q2 的决定性检查: K1/2 与响应强度的正相关, 是真的, 还是"未饱和"造成的机械假象?

混杂假设: 弱细胞的响应曲线在测试浓度范围内没饱和 -> 表观 Dmax 被压低
          -> 一半的门槛也跟着压低 -> 表观 K1/2 偏小。
检验: 用每细胞的饱和指数 sat = D_last / D_max;
     若弱细胞 sat 明显更低, 则该相关是截断假象;
     若在 sat 高的细胞里相关仍在, 才是真的结构。
"""
import os, re, io, json, contextlib
import numpy as np
import scipy.io as sio
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root
DATA = os.path.join(ROOT, "data_external", "moore2024_fret")
SRC  = os.path.join(ROOT, "code", "code54_K12_joint_fit_rupture_quantification.py")
src = open(SRC, encoding="utf-8").read()
src = re.sub(r'ROOT\s*=\s*"[^"]*"', 'ROOT=%r' % DATA, src)
src = re.sub(r'CSV\s*=\s*"[^"]*"',  'CSV=%r' % os.path.join(ROOT,"results","P7_Moore2024_FCD-Weber","code49b_unit_table.csv"), src)
ns = {"__name__":"code54"}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src,SRC,"exec"), ns)
FILES = ns["FILES"]

CURVES=[]
for name,bg in FILES.items():
    p=os.path.join(DATA,name+".mat")
    if not os.path.isfile(p): continue
    try: rd=sio.loadmat(p)["reorgData"]["resp_data"][0,0]
    except Exception: continue
    B=100.0 if name=="230831_FOV2" else float(bg)
    for ci in range(rd.shape[1]):
        try:
            a=rd["a"][0,ci].astype(float); s=rd["s"][0,ci].astype(float)
        except Exception: continue
        if a.shape!=(35,20) or s.shape!=(35,20): continue
        lev=defaultdict(list)
        for row in range(35):
            sr=s[row]; sv=float(sr.max())
            stim=np.where(sr>=sv-1e-9)[0]
            pre=np.arange(0,stim[0]) if len(stim) else None
            if pre is None or len(pre)<4 or len(stim)<6: continue
            da=float(np.median(a[row,pre[-4:]]))-float(np.median(a[row,stim[-6:]]))
            if sv-B>0: lev[round(sv,4)].append(da)
        good={k:float(np.median(v)) for k,v in lev.items() if len(v)>=4}
        if len(good)<3: continue
        ks=sorted(good)
        CURVES.append((B,np.array(ks,float),np.array([good[k] for k in ks],float)))

def k12(D,T):
    if D.max()<=0: return np.nan
    thr=0.5*D.max()
    if D.min()>=thr: return np.nan
    i=int(np.argmax(D>=thr))
    x0,x1=np.log10(T[i-1]),np.log10(T[i]); y0,y1=D[i-1],D[i]
    return 10**(x0+(thr-y0)*(x1-x0)/(y1-y0))

def spear(x,y):
    rx=np.argsort(np.argsort(x)); ry=np.argsort(np.argsort(y))
    return np.corrcoef(rx,ry)[0,1]

print("=== 饱和指数 sat = D_last / D_max ===")
print(" %6s %5s %9s %9s %11s" % ("B","n","med sat","p10","corr(sat,Dmax)"))
by_B=defaultdict(list)
for B,T,D in CURVES:
    by_B[B].append((D.max(), D[-1]/D.max() if D.max()>0 else np.nan, k12(D,T)))
for B in sorted(by_B):
    d=np.array(by_B[B])
    print(" %6.2f %5d %9.3f %9.3f %+11.3f"
          % (B, len(d), np.nanmedian(d[:,1]), np.nanpercentile(d[:,1],10),
             spear(d[:,0], d[:,1])))

print("\n=== 在饱和细胞里重测 K1/2 与强度的关系 ===")
print("  (只保留 sat >= 0.95 的细胞)")
print(" %6s %7s %9s %11s" % ("B","n_kept","med sat","rho(K12,Dmax)"))
allr=[]
for B in sorted(by_B):
    d=[x for x in by_B[B] if np.isfinite(x[2]) and x[1]>=0.95]
    if len(d)<10: continue
    d=np.array(d)
    r=spear(d[:,2], d[:,0]); allr.append(r)
    print(" %6.2f %7d %9.3f %+11.3f" % (B, len(d), np.median(d[:,1]), r))
print("  mean rho (saturated only) = %+.3f" % (np.mean(allr) if allr else float('nan')))

print("\n=== 饱和细胞中 弱/强 两组的 K1/2 差 ===")
print(" %6s %6s %10s %10s %9s" % ("B","n","weak","strong","dev(dex)"))
for B in sorted(by_B):
    d=[x for x in by_B[B] if np.isfinite(x[2]) and x[1]>=0.95]
    if len(d)<20: continue
    d=np.array(d); m=np.median(d[:,0])
    wk=d[d[:,0]<m,2]; st=d[d[:,0]>=m,2]
    if len(wk)>3 and len(st)>3:
        print(" %6.2f %6d %10.3f %10.3f %+9.3f"
              % (B, len(d), 10**np.median(np.log10(wk)), 10**np.median(np.log10(st)),
                 np.log10(10**np.median(np.log10(wk))/10**np.median(np.log10(st)))))

print("\n=== B=0 的 1.30 dex 离散是不是被少数细胞拉出来的 ===")
v=np.log10([x[2] for x in by_B[0] if np.isfinite(x[2])])
print("  n=%d  p5=%.2f p50=%.2f p95=%.2f  sd=%.2f"
      % (len(v), np.percentile(v,5), np.median(v), np.percentile(v,95), np.std(v)))
print("  去掉 |z|>3 后 sd=%.2f" % np.std(v[np.abs((v-np.median(v))/np.std(v))<=3]))

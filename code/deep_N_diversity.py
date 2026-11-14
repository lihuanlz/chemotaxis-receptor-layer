# -*- coding: utf-8 -*-
"""
paper_chemotaxis · 深挖 C
========================
两个问题:
  Q1  N=39 是数据定的, 还是 Ki-Ka-N-amax 脊上的一个点?  -> N 的 profile likelihood
  Q2  弱响应细胞是噪声, 还是结构?                        -> 每细胞响应分布 + K1/2 离散度
"""
import os, re, io, json, contextlib
import numpy as np
from scipy.optimize import least_squares
import scipy.io as sio
from collections import defaultdict
import csv as _csv

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root
SRC  = os.path.join(ROOT, "code", "code54_K12_joint_fit_rupture_quantification.py")
DATA = os.path.join(ROOT, "data_external", "moore2024_fret")
CSV  = os.path.join(ROOT, "results", "P7_Moore2024_FCD-Weber", "code49b_unit_table.csv")

src = open(SRC, encoding="utf-8").read()
src = re.sub(r'ROOT\s*=\s*"[^"]*"', 'ROOT=%r' % DATA, src)
src = re.sub(r'CSV\s*=\s*"[^"]*"',  'CSV=%r'  % CSV,  src)
ns = {"__name__": "code54"}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src, SRC, "exec"), ns)
FILES, gL, TGT, amp_model = ns["FILES"], ns["gL"], ns["TGT"], ns["amp_model"]

# ---------------- 抽曲线 ----------------
CURVES = []
TOT = defaultdict(int)
for name, bg in FILES.items():
    p = os.path.join(DATA, name + ".mat")
    if not os.path.isfile(p): continue
    try: rd = sio.loadmat(p)["reorgData"]["resp_data"][0, 0]
    except Exception: continue
    B = 100.0 if name == "230831_FOV2" else float(bg)
    for ci in range(rd.shape[1]):
        try:
            a = rd["a"][0, ci].astype(float); s = rd["s"][0, ci].astype(float)
        except Exception: continue
        if a.shape != (35,20) or s.shape != (35,20): continue
        lev = defaultdict(list)
        for row in range(35):
            sr = s[row]; sv = float(sr.max())
            stim = np.where(sr >= sv-1e-9)[0]
            pre = np.arange(0, stim[0]) if len(stim) else None
            if pre is None or len(pre) < 4 or len(stim) < 6: continue
            da = float(np.median(a[row, pre[-4:]])) - float(np.median(a[row, stim[-6:]]))
            if sv - B > 0: lev[round(sv,4)].append(da)
        good = {k: float(np.median(v)) for k,v in lev.items() if len(v) >= 4}
        if len(good) < 3: continue
        TOT[B] += 1
        ks = sorted(good)
        CURVES.append((B, np.array(ks,float), np.array([good[k] for k in ks],float)))
print("cells:", len(CURVES))

def k12_rel(B, T, D):
    if D.max() <= 0: return np.nan
    thr = 0.5*D.max()
    if D.min() >= thr: return np.nan
    i = int(np.argmax(D >= thr))
    x0,x1 = np.log10(T[i-1]), np.log10(T[i]); y0,y1 = D[i-1], D[i]
    return 10**(x0 + (thr-y0)*(x1-x0)/(y1-y0))

K12 = defaultdict(list)
DMAX = defaultdict(list)
for B,T,D in CURVES:
    k = k12_rel(B,T,D)
    if np.isfinite(k): K12[B].append(k)
    DMAX[B].append(D.max())

Bs = np.array(sorted(K12), float)
obs = np.array([10**np.nanmedian(np.log10(np.array(K12[b]))) for b in Bs])

rows = list(_csv.DictReader(open(CSV, encoding="utf-8")))
B_E = np.array([float(r["B_uM"]) for r in rows])
F_E = np.array([float(r["F_uM"]) for r in rows])
R_E = np.array([float(r["R_a"]) for r in rows])

def K12_total(B, Ki, Ka, N):
    Q = np.exp(gL(B,Ki,Ka) + TGT/N)
    den = 1.0/Ki - Q/Ka
    if den <= 0: return np.nan
    x = (Q-1.0)/den
    return x if x > 0 else np.nan

def fit_fixed_N(N, lam, n_starts=8, seed=3):
    """固定 N, 只拟合 Ki, Ka, amax。"""
    def resid(th):
        Ki, Ka, amax = np.exp(th)
        ra = (amp_model(B_E,F_E,Ki,Ka,N,amax) - R_E)/np.std(R_E)
        rk = []
        for B,o in zip(Bs,obs):
            x = K12_total(B,Ki,Ka,N)
            rk.append(np.log10(x/o) if np.isfinite(x) and x>0 else 5.0)
        return np.concatenate([ra, np.sqrt(lam)*np.array(rk)/0.5])
    best=None; rng=np.random.default_rng(seed); x0=np.log([10.,200.,1.4])
    for k in range(n_starts):
        xs = x0 if k==0 else x0 + rng.normal(0,0.7,3)
        try:
            sol = least_squares(resid, xs, bounds=(np.log([0.05,5.0,0.3]), np.log([300.,2e4,3.0])),
                                max_nfev=5000)
            ss = np.sum(resid(sol.x)**2)
            if best is None or ss<best[0]: best=(ss,sol.x)
        except Exception: pass
    Ki,Ka,amax = np.exp(best[1])
    pr = amp_model(B_E,F_E,Ki,Ka,N,amax)
    r2 = 1 - ((pr-R_E)**2).sum()/((R_E-R_E.mean())**2).sum()
    dev=[abs(np.log10(K12_total(B,Ki,Ka,N)/o)) for B,o in zip(Bs,obs)]
    return dict(N=N, r2=float(r2), err=float(np.mean(dev)), Ki=float(Ki), Ka=float(Ka), amax=float(amax))

print("\n=== Q1  N 的 profile likelihood (固定 N, 重拟合 Ki,Ka,amax) ===")
print(" %5s | %-31s | %-31s" % ("N", "amplitude only (lambda=0)", "joint (lambda=1)"))
print(" %5s | %7s %7s %6s %6s %6s | %7s %7s" % ("","R2","err","Ki","Ka","amax","R2","err"))
prof=[]
for N in [3,5,8,10,12,15,20,25,30,39,50,60]:
    a = fit_fixed_N(N, 0.0); b = fit_fixed_N(N, 1.0)
    prof.append((a,b))
    print(" %5d | %7.3f %7.3f %6.1f %6.0f %6.2f | %7.3f %7.3f"
          % (N, a["r2"], a["err"], a["Ki"], a["Ka"], a["amax"], b["r2"], b["err"]))

print("\n=== Q2a  每细胞最大响应分布 (da_max) ===")
print(" %6s %5s %8s %8s %8s %8s %8s" % ("B","n","med","p10","p90","frac<0.1","frac<0.2"))
for b in Bs:
    v = np.array(DMAX[b])
    print(" %6.2f %5d %8.3f %8.3f %8.3f %8.1f%% %8.1f%%"
          % (b, len(v), np.median(v), np.percentile(v,10), np.percentile(v,90),
             100*np.mean(v<0.1), 100*np.mean(v<0.2)))

print("\n=== Q2b  K1/2 的细胞间离散度 (相对定义, 总配体轴) ===")
print(" %6s %5s %9s %8s %8s %8s %9s" % ("B","n","median","IQRlo","IQRhi","sd_dex","CV_dex"))
sd_all=[]
for b in Bs:
    v = np.log10(np.array(K12[b]))
    sd_all.append(np.std(v))
    print(" %6.2f %5d %9.3f %8.3f %8.3f %8.3f %9.3f"
          % (b, len(v), 10**np.median(v), 10**np.percentile(v,25), 10**np.percentile(v,75),
             np.std(v), np.std(v)/abs(np.median(v))))
print("  pooled within-background sd (dex): %.3f   [range %.3f-%.3f]"
      % (np.mean(sd_all), min(sd_all), max(sd_all)))

print("\n=== Q2c  K1/2 是否依赖细胞的响应强度? (每背景内 Spearman) ===")
def spear(x,y):
    rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y))
    return np.corrcoef(rx,ry)[0,1]
print(" %6s %5s %9s" % ("B","n","rho(K12, da_max)"))
rhos=[]
for b in Bs:
    ks=[]; ds=[]
    for BB,T,D in CURVES:
        if BB!=b: continue
        k=k12_rel(BB,T,D)
        if np.isfinite(k): ks.append(k); ds.append(D.max())
    r = spear(np.array(ks), np.array(ds)); rhos.append(r)
    print(" %6.2f %5d %+9.3f" % (b, len(ks), r))
print("  mean rho = %+.3f" % np.mean(rhos))

print("\n=== Q2d  弱响应 vs 强响应 两组的中点 ===")
thr = np.median([np.median(np.array(DMAX[b])) for b in Bs])
print("  (split at da_max = %.3f, the median of per-background medians)" % thr)
print(" %6s %10s %10s %10s" % ("B","weak med","strong med","dev(dex)"))
for b in Bs:
    wk=[]; st=[]
    for BB,T,D in CURVES:
        if BB!=b: continue
        k=k12_rel(BB,T,D)
        if not np.isfinite(k): continue
        (wk if D.max()<thr else st).append(k)
    if wk and st:
        print(" %6.2f %10.3f %10.3f %+10.3f"
              % (b, 10**np.median(np.log10(wk)), 10**np.median(np.log10(st)),
                 np.log10(10**np.median(np.log10(wk))/10**np.median(np.log10(st)))))

json.dump({"profile":[[a,b] for a,b in prof]}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"deep.json"),"w"), indent=1)

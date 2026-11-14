# -*- coding: utf-8 -*-
"""细扫 lambda: 判定"没有点进入接受域"是否只是网格太粗。"""
import os, re, io, json, contextlib
import numpy as np
from scipy.optimize import least_squares

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
import scipy.io as sio
from collections import defaultdict
import csv as _csv

def curves():
    out = []
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
            ks = sorted(good)
            out.append((B, np.array(ks,float), np.array([good[k] for k in ks],float)))
    return out

CURVES = curves()

def table(mode, frac=0.5):
    k12 = defaultdict(list)
    for B,T,D in CURVES:
        if mode=="abs":
            if D.max() < 0.5: continue
            thr = 0.5
        else:
            if D.max() <= 0: continue
            thr = frac*D.max()
        if D.min() >= thr: continue
        i = int(np.argmax(D >= thr))
        x0,x1 = np.log10(T[i-1]), np.log10(T[i]); y0,y1 = D[i-1], D[i]
        k12[B].append(10**(x0 + (thr-y0)*(x1-x0)/(y1-y0)))
    Bs = np.array([b for b in sorted(k12) if k12[b]], float)
    obs = np.array([10**np.nanmedian(np.log10(np.array(k12[b]))) for b in Bs])
    return Bs, obs

rows = list(_csv.DictReader(open(CSV, encoding="utf-8")))
B_E = np.array([float(r["B_uM"]) for r in rows]); F_E = np.array([float(r["F_uM"]) for r in rows])
R_E = np.array([float(r["R_a"]) for r in rows])

def K12_total(B, Ki, Ka, N):
    Q = np.exp(gL(B,Ki,Ka) + TGT/N)
    den = 1.0/Ki - Q/Ka
    if den <= 0: return np.nan
    x = (Q-1.0)/den
    return x if x > 0 else np.nan

LB = np.log([0.05,5.0,2.0,0.3]); UB = np.log([300.0,2e4,60.0,3.0])
def fit(lam, Bs, obs, axis, n_starts=10, seed=1):
    def resid(th):
        Ki,Ka,N,amax = np.exp(th)
        ra = (amp_model(B_E,F_E,Ki,Ka,N,amax) - R_E)/np.std(R_E)
        rk = []
        for B,o in zip(Bs,obs):
            x = K12_total(B,Ki,Ka,N)
            if not np.isfinite(x): rk.append(5.0); continue
            pred = x if axis=="total" else x-B
            if pred <= 0: rk.append(5.0); continue
            rk.append(np.log10(pred/o))
        return np.concatenate([ra, np.sqrt(lam)*np.array(rk)/0.5])
    best=None; rng=np.random.default_rng(seed); x0=np.log([2.,200.,6.,1.])
    for k in range(n_starts):
        xs = x0 if k==0 else x0 + rng.normal(0,0.8,4)
        try:
            sol = least_squares(resid, xs, bounds=(LB,UB), max_nfev=6000)
            ss = np.sum(resid(sol.x)**2)
            if best is None or ss<best[0]: best=(ss,sol.x)
        except Exception: pass
    Ki,Ka,N,amax = np.exp(best[1])
    pr = amp_model(B_E,F_E,Ki,Ka,N,amax)
    r2 = 1 - ((pr-R_E)**2).sum()/((R_E-R_E.mean())**2).sum()
    dev=[]
    for B,o in zip(Bs,obs):
        x = K12_total(B,Ki,Ka,N); pred = x if axis=="total" else x-B
        dev.append(abs(np.log10(pred/o)) if pred>0 else 9.9)
    return dict(lam=lam,r2=float(r2),err=float(np.mean(dev)),Ki=float(Ki),Ka=float(Ka),N=float(N),amax=float(amax))

Bs_abs, obs_abs = table("abs")
Bs_rel, obs_rel = table("rel")

def fine(tag, Bs, obs, axis, lams):
    print("\n=== %s ===" % tag)
    print("  %6s %7s %7s  %-38s %s" % ("lambda","R^2","err","params","accept?"))
    best = None
    for l in lams:
        r = fit(l, Bs, obs, axis)
        ok = (r["r2"]>0.9 and r["err"]<0.1)
        # 距离接受域的归一化距离, 用于找最近点
        d = max((0.9-r["r2"])/0.9, (r["err"]-0.1)/0.1, 0)
        if best is None or d < best[0]: best=(d,r,ok)
        print("  %6.2f %7.3f %7.3f  Ki=%-6.1f Ka=%-6.0f N=%-5.1f amax=%.2f  %s"
              % (r["lam"], r["r2"], r["err"], r["Ki"], r["Ka"], r["N"], r["amax"],
                 "YES" if ok else ""))
    print("  -> closest to acceptance: lambda=%.2f  R2=%.3f err=%.3f"
          % (best[1]["lam"], best[1]["r2"], best[1]["err"]))
    return best

fine("B protocol: axis fixed, ABSOLUTE 0.5 threshold",
     Bs_abs, obs_abs, "total", [2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,7.0])
fine("D protocol: axis fixed, RELATIVE half-max",
     Bs_rel, obs_rel, "total", [0.0,0.5,1.0,2.0,5.0,10.0])

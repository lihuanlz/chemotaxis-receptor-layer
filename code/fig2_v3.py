# -*- coding: utf-8 -*-
"""fig2_v3: (a) N profile likelihood (b) per-background K1/2 spread (c) sensitivity vs gain."""
import os, re, io, contextlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.io as sio
from collections import defaultdict
from scipy.optimize import least_squares
import csv as _csv

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root
DATA = os.path.join(ROOT, "data_external", "moore2024_fret")
SRC  = os.path.join(ROOT, "code", "code54_K12_joint_fit_rupture_quantification.py")
CSV  = os.path.join(ROOT, "results", "P7_Moore2024_FCD-Weber", "code49b_unit_table.csv")
OUT  = os.path.dirname(os.path.abspath(__file__))

src = open(SRC, encoding="utf-8").read()
src = re.sub(r'ROOT\s*=\s*"[^"]*"', 'ROOT=%r' % DATA, src)
src = re.sub(r'CSV\s*=\s*"[^"]*"',  'CSV=%r'  % CSV,  src)
ns = {"__name__": "code54"}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src, SRC, "exec"), ns)
FILES, gL, TGT, amp_model = ns["FILES"], ns["gL"], ns["TGT"], ns["amp_model"]

CURVES = []
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
        CURVES.append((B, np.array(ks,float), np.array([good[k] for k in ks],float)))

def k12(T, D):
    """半最大插值。守卫: 插值分数必须落在 [0,1] 附近, 否则该细胞的中点在两级之间不可定义。"""
    if D.max() <= 0: return np.nan
    thr = 0.5*D.max()
    if D.min() >= thr: return np.nan
    i = int(np.argmax(D >= thr))
    dy = D[i] - D[i-1]
    if dy <= 1e-6: return np.nan
    f = (thr - D[i-1]) / dy
    if not (-0.05 <= f <= 1.05): return np.nan
    x0,x1 = np.log10(T[i-1]), np.log10(T[i])
    return 10**(x0 + f*(x1-x0))

K12 = defaultdict(list); DMAX = defaultdict(list)
for B,T,D in CURVES:
    k = k12(T,D)
    if np.isfinite(k): K12[B].append(k)
    DMAX[B].append(D.max())
Bs = np.array(sorted(K12), float)
obs = np.array([10**np.nanmedian(np.log10(np.array(K12[b]))) for b in Bs])

rows = list(_csv.DictReader(open(CSV, encoding="utf-8")))
B_E = np.array([float(r["B_uM"]) for r in rows]); F_E = np.array([float(r["F_uM"]) for r in rows])
R_E = np.array([float(r["R_a"]) for r in rows])

def K12_total(B, Ki, Ka, N):
    Q = np.exp(gL(B,Ki,Ka) + TGT/N)
    den = 1.0/Ki - Q/Ka
    if den <= 0: return np.nan
    x = (Q-1.0)/den
    return x if x > 0 else np.nan

def fit_N(N, lam=0.0, n_starts=6, seed=3):
    def resid(th):
        Ki,Ka,amax = np.exp(th)
        ra = (amp_model(B_E,F_E,Ki,Ka,N,amax)-R_E)/np.std(R_E)
        rk=[np.log10(K12_total(B,Ki,Ka,N)/o) for B,o in zip(Bs,obs)]
        return np.concatenate([ra, np.sqrt(lam)*np.array(rk)/0.5])
    best=None; rng=np.random.default_rng(seed); x0=np.log([10.,200.,1.4])
    for k in range(n_starts):
        xs = x0 if k==0 else x0+rng.normal(0,0.7,3)
        try:
            sol=least_squares(resid,xs,bounds=(np.log([0.05,5.,0.3]),np.log([300.,2e4,3.])),max_nfev=5000)
            ss=np.sum(resid(sol.x)**2)
            if best is None or ss<best[0]: best=(ss,sol.x)
        except Exception: pass
    Ki,Ka,amax=np.exp(best[1])
    pr=amp_model(B_E,F_E,Ki,Ka,N,amax)
    r2=1-((pr-R_E)**2).sum()/((R_E-R_E.mean())**2).sum()
    err=float(np.mean([abs(np.log10(K12_total(B,Ki,Ka,N)/o)) for B,o in zip(Bs,obs)]))
    return dict(N=N, r2=float(r2), err=err, Ka=float(Ka), Ki=float(Ki))

Ns=[3,5,8,10,12,15,20,25,30,39,50,60]
P=[fit_N(n) for n in Ns]
print("profile:", [(p["N"], round(p["r2"],3), round(p["err"],3), int(p["Ka"])) for p in P])

# ---- 每细胞 (背景中心化 log K1/2 vs da_max)，只用饱和细胞 ----
pts=[]
for B,T,D in CURVES:
    if D.max()<=0 or D[-1]/D.max()<0.95: continue
    k=k12(T,D)
    if np.isfinite(k): pts.append((B, D.max(), np.log10(k)))
pts=np.array(pts)
rel=np.zeros(len(pts))
for b in np.unique(pts[:,0]):
    m=pts[:,0]==b
    rel[m]=pts[m,2]-np.median(pts[m,2])

fig = plt.figure(figsize=(11, 3.6))
gs = fig.add_gridspec(1, 3, wspace=0.42)

# (a) N profile
ax = fig.add_subplot(gs[0,0])
ax.plot([p["N"] for p in P],[p["r2"] for p in P],'o-',color="#1f6fb4",ms=4,lw=1,label="$R^2$ (amplitude)")
ax.axhline(0.9,color='k',ls=':',lw=.8)
ax2 = ax.twinx()
ax2.plot([p["N"] for p in P],[p["Ka"] for p in P],'s--',color="#c94f2b",ms=3.5,lw=1,label="$K_A$")
ax2.set_yscale('log'); ax2.set_ylabel("$K_A$", color="#c94f2b", fontsize=8, labelpad=2)
ax2.tick_params(axis='y', labelcolor="#c94f2b", labelsize=7)
ax2.axhline(2e4, color="#c94f2b", ls=':', lw=.6, alpha=.6)
ax.set_xlabel("cooperativity $N$ (fixed)")
ax.set_ylabel("amplitude $R^2$", color="#1f6fb4")
ax.set_ylim(0.7,1.02); ax.tick_params(labelsize=7)
ax.set_title("a  $N$ is not identifiable", loc='left', fontsize=9)
ax.grid(alpha=.25,lw=.5)

# (b) per-background K1/2 spread
ax = fig.add_subplot(gs[0,1])
data=[np.array(K12[b]) for b in Bs]
bp=ax.boxplot(data, positions=np.arange(len(Bs)), widths=.55, showfliers=False,
              medianprops=dict(color="#c94f2b",lw=1.2), boxprops=dict(color="#1f6fb4"),
              whiskerprops=dict(color="#1f6fb4"), capprops=dict(color="#1f6fb4"))
ax.plot(np.arange(len(Bs)), obs, 'o', ms=3.5, color="#222", label="median")
ax.set_yscale('log'); ax.set_xticks(np.arange(len(Bs)))
ax.set_xticklabels(["%g"%b for b in Bs], fontsize=6.5)
ax.set_xlabel("background $B$ ($\\mu$M)"); ax.set_ylabel("$K_{1/2}$ ($\\mu$M)")
ax.set_title("b  spread is state dependent", loc='left', fontsize=9)
ax.legend(fontsize=6.5, frameon=False); ax.tick_params(labelsize=7); ax.grid(alpha=.25,lw=.5)

# (c) sensitivity vs gain
ax = fig.add_subplot(gs[0,2])
ax.scatter(pts[:,1], rel, s=7, alpha=.35, color="#1f6fb4", edgecolors='none')
bins=np.quantile(pts[:,1], np.linspace(0,1,7))
cent=[]; med=[]
for i in range(len(bins)-1):
    m=(pts[:,1]>=bins[i])&(pts[:,1]<=bins[i+1])
    if m.sum()>5: cent.append(np.median(pts[m,1])); med.append(np.median(rel[m]))
ax.plot(cent, med, 'o-', color="#c94f2b", ms=4, lw=1.2, label="binned median")
from scipy.stats import spearmanr
rho,p = spearmanr(pts[:,1], rel)
ax.axhline(0,color='k',lw=.5)
ax.set_xlabel("per-cell maximal response  $\\Delta a_{\\max}$")
ax.set_ylabel("$\\log_{10} K_{1/2}$ (centred per background)")
ax.set_title("c  sensitivity co-varies with gain", loc='left', fontsize=9)
ax.text(0.03,0.95,"Spearman $\\rho$ = +%.3f\n$p$ = %.0e,  n = %d" % (rho,p,len(pts)),
        transform=ax.transAxes, fontsize=7, va='top')
ax.legend(fontsize=6.5, frameon=False, loc='lower right')
ax.tick_params(labelsize=7); ax.grid(alpha=.25,lw=.5)

fig.savefig(os.path.join(OUT,"fig2_v3.png"), dpi=200, bbox_inches='tight')
print("wrote fig2_v3.png  | rho=%.3f p=%.2e n=%d" % (rho,p,len(pts)))



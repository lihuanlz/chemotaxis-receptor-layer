# -*- coding: utf-8 -*-
"""fig1_v2: 修正后的图 —— 两种中点定义、两条 Pareto 前沿、逐背景残差、删失比例。"""
import os, re, io, json, contextlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict
import scipy.io as sio
import csv as _csv

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root
SRC  = os.path.join(ROOT, "code", "code54_K12_joint_fit_rupture_quantification.py")
DATA = os.path.join(ROOT, "data_external", "moore2024_fret")
CSV  = os.path.join(ROOT, "results", "P7_Moore2024_FCD-Weber", "code49b_unit_table.csv")
OUT  = os.path.dirname(os.path.abspath(__file__))

src = open(SRC, encoding="utf-8").read()
src = re.sub(r'ROOT\s*=\s*"[^"]*"', 'ROOT=%r' % DATA, src)
src = re.sub(r'CSV\s*=\s*"[^"]*"',  'CSV=%r'  % CSV,  src)
ns = {"__name__": "code54"}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src, SRC, "exec"), ns)
FILES, gL, TGT = ns["FILES"], ns["gL"], ns["TGT"]

# ---------- 抽曲线 ----------
def load_curves():
    out, tot = [], defaultdict(int)
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
            tot[B] += 1
            ks = sorted(good)
            out.append((B, np.array(ks,float), np.array([good[k] for k in ks],float)))
    return out, tot

CURVES, TOT = load_curves()

def table(mode):
    k12, cens = defaultdict(list), defaultdict(int)
    for B, T, D in CURVES:
        if mode == "abs":
            if D.max() < 0.5:
                cens[B] += 1; continue
            thr = 0.5
        else:
            if D.max() <= 0:
                cens[B] += 1; continue
            thr = 0.5 * D.max()
        if D.min() >= thr: continue
        i = int(np.argmax(D >= thr))
        x0, x1 = np.log10(T[i-1]), np.log10(T[i]); y0, y1 = D[i-1], D[i]
        k12[B].append(10 ** (x0 + (thr-y0)*(x1-x0)/(y1-y0)))
    Bs = np.array([b for b in sorted(k12) if k12[b]], float)
    obs = np.array([10**np.nanmedian(np.log10(np.array(k12[b]))) for b in Bs])
    q = []
    for b in Bs:
        lv = np.log10(np.array(k12[b])); q.append(np.exp(np.percentile(lv, [25,75])))
    return Bs, obs, np.array(q), cens

Bs_a, obs_a, qa, cens_a = table("abs")
Bs_r, obs_r, qr, cens_r = table("rel")

def K12_total(B, Ki, Ka, N):
    Q = np.exp(gL(B,Ki,Ka) + TGT/N)
    den = 1.0/Ki - Q/Ka
    if den <= 0: return np.nan
    x = (Q-1.0)/den
    return x if x > 0 else np.nan

P0 = (30.7, 133.0, 39.0)
pred = np.array([K12_total(b, *P0) for b in Bs_r])

# 两条 Weber 线（log-log 线性拟合）
pa = np.polyfit(np.log10(np.maximum(Bs_a, 1e-3)), np.log10(obs_a), 1)
pr = np.polyfit(np.log10(np.maximum(Bs_r, 1e-3)), np.log10(obs_r), 1)

res = json.load(open(os.path.join(OUT, "recheck_results.json"), encoding="utf-8"))

fig = plt.figure(figsize=(11, 8.2))
gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.26)

# (a) 两种中点统计
ax = fig.add_subplot(gs[0, 0])
ea = [np.maximum(0.0, obs_a-qa[:,0]), np.maximum(0.0, qa[:,1]-obs_a)]
er = [np.maximum(0.0, obs_r-qr[:,0]), np.maximum(0.0, qr[:,1]-obs_r)]
ax.errorbar(np.maximum(Bs_a, 1e-4), obs_a, yerr=ea, fmt='o', ms=5,
            color="#1f6fb4", label="absolute 0.5 threshold", capsize=2, lw=1)
ax.errorbar(np.maximum(Bs_r, 1e-4)*1.08, obs_r, yerr=er, fmt='s', ms=5,
            color="#c94f2b", label="half of own maximum", capsize=2, lw=1)
# 仿射 Weber 线 K = a + s*B
fa = np.polyfit(Bs_a, obs_a, 1)
fr = np.polyfit(Bs_r, obs_r, 1)
bb = np.logspace(-2, 2, 200)
ax.plot(bb, fa[0]*bb + fa[1], '-', color="#1f6fb4", lw=1, alpha=.7)
ax.plot(bb, fr[0]*bb + fr[1], '-', color="#c94f2b", lw=1, alpha=.7)
dev_a = np.max(np.abs(np.log10((fa[0]*Bs_a+fa[1])/obs_a)))
dev_r = np.max(np.abs(np.log10((fr[0]*Bs_r+fr[1])/obs_r)))
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(5e-3, 3e2); ax.set_ylim(0.5, 3e2)
ax.set_xlabel(r"background MeAsp  $B$ ($\mu$M)")
ax.set_ylabel(r"midpoint  $K_{1/2}$ ($\mu$M)")
ax.set_title("a  two midpoint definitions, same cells", loc='left', fontsize=10)
ax.text(0.03, 0.95,
        "abs:  $K=%.2f+%.2f\\,B$,  max dev %.2f dex\nrel:  $K=%.2f+%.2f\\,B$,  max dev %.2f dex"
        % (fa[1], fa[0], dev_a, fr[1], fr[0], dev_r),
        transform=ax.transAxes, fontsize=7, va='top')
ax.legend(fontsize=7.5, frameon=False, loc='lower right')
ax.grid(alpha=.25, lw=.5)

# (b) 修正前后的 Pareto 前沿
ax = fig.add_subplot(gs[0, 1])
for key, lab, col in (("A", "as published (F axis, absolute)", "#888888"),
                      ("B", "axis corrected only", "#1f6fb4"),
                      ("D", "axis + definition corrected", "#c94f2b")):
    r = res[key]
    ax.plot([x["err"] for x in r], [x["r2"] for x in r], 'o-', ms=3.5, lw=1, color=col, label=lab)
ax.axvrect = None
ax.axvline(0.1, color='k', ls=':', lw=.8)
ax.axhline(0.9, color='k', ls=':', lw=.8)
ax.fill_between([0, 0.1], 0.9, 1.0, color='green', alpha=.08)
ax.set_xlim(0, 0.7); ax.set_ylim(0.3, 1.02)
ax.set_xlabel(r"mean midpoint error (dex)")
ax.set_ylabel(r"amplitude $R^2$")
ax.set_title("b  Pareto front under three protocols", loc='left', fontsize=10)
ax.legend(fontsize=7.5, frameon=False, loc='lower right')
ax.grid(alpha=.25, lw=.5)

# (c) 逐背景残差（振幅最优参数）
ax = fig.add_subplot(gs[1, 0])
ra = np.array([np.log10(K12_total(b, *P0)/o) for b, o in zip(Bs_a, obs_a)])
rr = np.array([np.log10(K12_total(b, *P0)/o) for b, o in zip(Bs_r, obs_r)])
w = 0.35
ax.bar(np.arange(len(Bs_a))-w/2, ra, w, color="#1f6fb4", label="absolute 0.5 (mean |dev| 0.258)")
ax.bar(np.arange(len(Bs_r))+w/2, rr, w, color="#c94f2b", label="half of own max (mean |dev| 0.023)")
ax.axhline(0, color='k', lw=.6)
ax.set_xticks(np.arange(len(Bs_a)))
ax.set_xticklabels(["%g" % b for b in Bs_a], fontsize=7)
ax.set_xlabel(r"background $B$ ($\mu$M)")
ax.set_ylabel(r"$\log_{10}$(model / measured)")
ax.set_title("c  residuals at the amplitude-optimal point", loc='left', fontsize=10)
ax.legend(fontsize=7.5, frameon=False, loc='lower right')
ax.grid(alpha=.25, lw=.5, axis='y')

# (d) 删失比例
ax = fig.add_subplot(gs[1, 1])
fa = [100*cens_a.get(b,0)/TOT[b] for b in Bs_a]
fr = [100*cens_r.get(b,0)/TOT[b] for b in Bs_r]
ax.plot(np.arange(len(Bs_a)), fa, 'o-', color="#1f6fb4", ms=4, lw=1, label="absolute 0.5")
ax.plot(np.arange(len(Bs_a)), fr, 's-', color="#c94f2b", ms=4, lw=1, label="half of own max")
ax.set_xticks(np.arange(len(Bs_a))); ax.set_xticklabels(["%g" % b for b in Bs_a], fontsize=7)
ax.set_xlabel(r"background $B$ ($\mu$M)"); ax.set_ylabel("cells censored (%)")
ax.set_title("d  censoring is a property of the threshold", loc='left', fontsize=10)
ax.legend(fontsize=7.5, frameon=False)
ax.grid(alpha=.25, lw=.5)

fig.savefig(os.path.join(OUT, "fig1_v2.png"), dpi=200, bbox_inches='tight')
print("wrote fig1_v2.png")
print("absolute-threshold Weber fit : log10 K = %.4f*log10 B_eff + ...  (slope %.3f)" % (pa[0], pa[0]))
print("relative-threshold Weber fit : slope %.3f" % pr[0])
json.dump({"obs_abs": obs_a.tolist(), "obs_rel": obs_r.tolist(),
           "Bs": Bs_a.tolist(),
           "cens_abs": {str(b): cens_a.get(b,0) for b in Bs_a},
           "cens_rel": {str(b): cens_r.get(b,0) for b in Bs_a},
           "tot": {str(b): TOT[b] for b in Bs_a}},
          open(os.path.join(OUT, "fig1_v2_data.json"), "w"), indent=1)

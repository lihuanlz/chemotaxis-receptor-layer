# make_fig1.py -- regenerate fig1_chemotaxis.png (6 panels) from archived numbers.
# Data sources (verified against _rerun/logs, 2026-10-15 audit):
#   code54  : whole-population Pareto front (7 lambda values), medians/IQRs, censoring counts
#   code62  : same-source front + whole-population control (lambda = 0, 1, 30)
#   code55/56/57/60 : repair-path scores (#2, #4, #6, #7)
#   code70  : prospective arms V1 / V2 / V2_alt
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

BLUE = "#1f77b4"; ORANGE = "#e08214"; GRAY = "#7f7f7f"
LBLUE = "#56B4E9"; LORANGE = "#E69F00"; GREEN = "#2e8b57"

fig, axes = plt.subplots(2, 3, figsize=(12.4, 7.1))

# ---------------- panel a: Pareto front ----------------
axa = axes[0][0]
lam = [0, 0.1, 0.3, 1, 3, 10, 30]
front = [(0.980, 0.622), (0.977, 0.582), (0.960, 0.503), (0.831, 0.224),
         (0.737, 0.157), (0.578, 0.097), (0.511, 0.074)]
ss = [(0.971, 0.573), (0.868, 0.205), (0.599, 0.071)]          # same-source, lam=0,1,30
wp = [(0.976, 0.624), (0.848, 0.261), (0.525, 0.074)]          # whole-population control (code 62)
axa.plot([p[0] for p in front], [p[1] for p in front], "-o", color=BLUE, ms=7, lw=1.6)
axa.plot([p[0] for p in ss], [p[1] for p in ss], "--s", color=ORANGE, ms=6, lw=1.4,
         mfc="none", label="same-source (code 62)")
axa.plot([p[0] for p in wp], [p[1] for p in wp], ":d", color=GRAY, ms=6, lw=1.4,
         label="whole-population control")
lab_off = {0: (0, 9, "center"), 0.1: (10, -13, "left"), 0.3: (10, -13, "left"),
           1: (-12, 6, "right"), 3: (-12, 6, "right"), 10: (0, 8, "center"),
           30: (-6, -14, "right")}
for l, (x, y) in zip(lam, front):
    dx, dy, ha = lab_off[l]
    axa.annotate(f"lam={l:g}", (x, y), textcoords="offset points",
                 xytext=(dx, dy), ha=ha, fontsize=7, color=BLUE)
axa.add_patch(Rectangle((0.9, 0.0), 0.1, 0.1, fc="#d5e8d4", ec="none", zorder=0))
axa.text(0.965, 0.118, "acceptance region\n$R^2$>0.9, err<0.1 dex", fontsize=7,
         color=GREEN, ha="center", va="bottom")
axa.set_xlim(0.45, 1.03); axa.set_ylim(0.0, 0.70)
axa.set_xlabel("amplitude $R^2$", fontsize=9)
axa.set_ylabel(r"$K_{1/2}$ mean $|\Delta\log_{10}|$ (dex)", fontsize=9)
axa.set_title("Pareto front: no point enters the acceptance region", fontsize=10)
axa.legend(fontsize=7.5, loc="upper left", frameon=False)
axa.tick_params(labelsize=8)
axa.text(-0.16, 1.04, "a", transform=axa.transAxes, fontsize=13, fontweight="bold")

# ---------------- panel b: eight repair paths ----------------
axb = axes[0][1]
axb.plot(0.980, 0.622, "o", color=BLUE, ms=11)
axb.annotate("standard (lam=0)", (0.980, 0.622), textcoords="offset points",
             xytext=(-9, 8), ha="right", fontsize=7.5, color=BLUE)
axb.plot(0.511, 0.074, "o", color=BLUE, ms=11)
axb.annotate("standard (lam=30)", (0.511, 0.074), textcoords="offset points",
             xytext=(0, -15), ha="center", fontsize=7.5, color=BLUE)
rep = [(0.950, 0.257, "#4", (9, -2)), (0.953, 0.182, "#2", None),
       (0.946, 0.179, "#7", None), (0.678, 0.052, "#6", (0, 9))]
for x, y, tag, off in rep:
    axb.plot(x, y, "o", color="#cc79a7", ms=11)
axb.annotate("#4", (0.950, 0.257), textcoords="offset points", xytext=(9, -2),
             ha="left", fontsize=8, color="#cc79a7")
axb.annotate("#2,#7", (0.9495, 0.1805), textcoords="offset points", xytext=(0, 11),
             ha="center", fontsize=8, color="#cc79a7")
axb.annotate("#6", (0.678, 0.052), textcoords="offset points", xytext=(0, 9),
             ha="center", fontsize=8, color="#cc79a7")
txt = ("#1 two subpopulations: directional failure (alpha -> 0.05-0.08)\n"
       "#3 protocol mismatch: excluded as primary cause\n"
       "#5 TCS depletion: front digit-for-digit unchanged\n"
       "#8 asymmetric N: n(B) 1->17 within single-N envelope")
axb.text(0.467, 0.455, txt, fontsize=7.2, color="0.25", va="top")
axb.add_patch(Rectangle((0.9, 0.0), 0.1, 0.1, fc="#d5e8d4", ec="none", zorder=0))
axb.text(0.965, 0.115, "acceptance\nregion", fontsize=7, color=GREEN,
         ha="center", va="bottom")
axb.set_xlim(0.45, 1.03); axb.set_ylim(0.0, 0.70)
axb.set_xlabel("amplitude $R^2$", fontsize=9)
axb.set_ylabel(r"$K_{1/2}$ error (dex)", fontsize=9)
axb.set_title("eight repair paths: all outside the acceptance region", fontsize=10)
axb.tick_params(labelsize=8)
axb.text(-0.16, 1.04, "b", transform=axb.transAxes, fontsize=13, fontweight="bold")

# ---------------- panel c: Weber line ----------------
axc = axes[0][2]
B = [0, 0.01, 0.1, 0.3, 1, 10, 100]
med = [2.03, 2.90, 2.17, 2.58, 3.49, 13.92, 119.56]
lo = [1.58, 2.08, 1.77, 1.79, 2.82, 13.21, 112.41]
hi = [3.10, 3.06, 2.84, 3.64, 4.16, 15.90, 132.87]
x0 = 4e-3  # off-scale position for B=0 on the log axis
xs = [x0] + B[1:]
import numpy as np
xg = np.logspace(math.log10(x0), math.log10(300), 200)
axc.plot(xg, 1.17 * (1.95 + xg), "--", color=ORANGE, lw=1.6,
         label="Weber line 1.17(1.95 + B)")
axc.errorbar(xs, med, yerr=[[m - l for m, l in zip(med, lo)],
                            [h - m for m, h in zip(med, hi)]],
             fmt="o", color=BLUE, ms=7, capsize=3, label="measured $K_{1/2}$ (IQR)")
axc.set_xscale("log"); axc.set_yscale("log")
axc.set_xlim(3e-3, 3e2); axc.set_ylim(1.5, 300)
axc.set_xticks([x0, 0.01, 0.1, 1, 10, 100])
axc.set_xticklabels(["0", "0.01", "0.1", "1", "10", "100"])
axc.set_xlabel("background $B$ (uM MeAsp)", fontsize=9)
axc.set_ylabel("$K_{1/2}$ (uM)", fontsize=9)
axc.set_title("Weber line: max deviation 0.10 dex, mean 0.03 dex", fontsize=10)
axc.legend(fontsize=7.5, loc="lower right", frameon=False)
axc.tick_params(labelsize=8)
axc.text(-0.16, 1.04, "c", transform=axc.transAxes, fontsize=13, fontweight="bold")
res = [math.log10(m / (1.17 * (1.95 + b))) for m, b in zip(med, B)]
axi = axc.inset_axes([0.06, 0.55, 0.44, 0.40])
axi.bar(range(7), res, width=0.6, color=LBLUE)
axi.axhline(0, color="0.3", lw=0.8)
axi.set_ylim(-0.12, 0.12)
axi.set_ylabel("residual (dex)", fontsize=7)
axi.set_xlabel("B (uM)", fontsize=7, labelpad=1)
axi.set_xticks(range(7))
axi.set_xticklabels(["0", ".01", ".1", ".3", "1", "10", "100"], fontsize=6)
axi.tick_params(labelsize=6)

# ---------------- panel d: right-censoring ----------------
axd = axes[1][0]
est = [89, 32, 22, 32, 25, 56, 80]
cen = [179, 69, 35, 38, 40, 46, 159]
tot = [e + c for e, c in zip(est, cen)]
pct = [int(100.0 * c / t + 0.5) for c, t in zip(cen, tot)]
xlab = ["0", "0.01", "0.1", "0.3", "1", "10", "100"]
axd.bar(range(7), est, width=0.6, color=BLUE, label="estimable")
axd.bar(range(7), cen, width=0.6, bottom=est, color="#c8c8c8", label="right-censored")
for i, (t, p) in enumerate(zip(tot, pct)):
    axd.text(i, t + 5, f"{p}%", ha="center", fontsize=7.5, color="0.25")
axd.set_xticks(range(7)); axd.set_xticklabels(xlab)
axd.set_ylim(0, 300)
axd.set_xlabel("background $B$ (uM)", fontsize=9)
axd.set_ylabel("cells", fontsize=9)
axd.set_title("right-censoring: 45-68% of cells per background", fontsize=10)
axd.legend(fontsize=7.5, loc="upper center", frameon=False)
axd.tick_params(labelsize=8)
axd.text(-0.16, 1.04, "d", transform=axd.transAxes, fontsize=13, fontweight="bold")

# ---------------- panel e: prospective arms ----------------
axe = axes[1][1]
xe = np.linspace(0, 520, 100)
axe.fill_between(xe, xe * 10**-0.15, xe * 10**0.15, color="#d5e8d4", alpha=0.6,
                 zorder=0)
axe.plot(xe, xe, "-", color="0.45", lw=1.2)
arms = [(388, 401.3, BLUE, "V1: pred 388, meas 401.3"),
        (122, 133.2, "#d55e00", "V2: pred 122, meas 133.2"),
        (122, 130.8, GRAY, "V2_alt: pred 122, meas 130.8")]
for x, y, c, lb in arms:
    axe.plot(x, y, "o", color=c, ms=9, label=lb)
axe.set_xlim(0, 520); axe.set_ylim(0, 450)
axe.set_xlabel("predicted $K_{1/2}$ (uM, effective axis)", fontsize=9)
axe.set_ylabel("measured $K_{1/2}$ (uM)", fontsize=9)
axe.set_title("prospective arms (code 70): P-2 hit 2/2", fontsize=10)
axe.legend(fontsize=7.5, loc="upper left", frameon=False)
axe.text(505, 60, "+/-0.15 dex frozen tolerance", fontsize=7.5, color=GREEN,
         ha="right")
axe.tick_params(labelsize=8)
axe.text(-0.16, 1.04, "e", transform=axe.transAxes, fontsize=13, fontweight="bold")

# ---------------- panel f: amplitude & structural clauses ----------------
axf = axes[1][2]
grp = ["V1", "V2", "V2_alt (control)"]
amp = [0.105, 0.126, 0.069]
strc = [-0.005, 0.177, 0.202]
xp = np.arange(3)
b1 = axf.bar(xp - 0.17, amp, width=0.32, color=LORANGE,
             label="P-1b max amplitude dev")
b2 = axf.bar(xp + 0.17, strc, width=0.32, color=LBLUE,
             label="P-1a $R^2$(logF) $- R^2$(logr)")
l1, = axf.plot([-0.5, 2.5], [0.3, 0.3], "--", color="#0072b2", lw=1.3,
               label="P-1a decision line +0.3")
l2, = axf.plot([-0.5, 2.5], [0.08, 0.08], "--", color=LORANGE, lw=1.3,
               label="P-1b tolerance 0.08")
for x, v in zip(xp - 0.17, amp):
    axf.text(x, v + 0.006, f"{v:.3f}", ha="center", fontsize=7.5)
for x, v in zip(xp + 0.17, strc):
    axf.text(x, v + (0.006 if v > 0 else -0.024), f"{v:+.3f}", ha="center",
             fontsize=7.5)
axf.set_xticks(xp); axf.set_xticklabels(grp, fontsize=8)
axf.set_xlim(-0.5, 2.5); axf.set_ylim(-0.05, 0.345)
axf.set_ylabel("deviation", fontsize=9)
axf.set_title("P-1b missed at low F; P-1a underpowered", fontsize=10, pad=46)
axf.legend(handles=[b1, b2, l1, l2], fontsize=7.0, loc="lower center",
           bbox_to_anchor=(0.5, 1.005), ncol=2, frameon=False,
           columnspacing=0.9, handlelength=1.6)
axf.tick_params(labelsize=8)
axf.text(-0.16, 1.04, "f", transform=axf.transAxes, fontsize=13, fontweight="bold")

fig.tight_layout()
fig.savefig("fig1_chemotaxis.png", dpi=121)
print("saved fig1_chemotaxis.png")

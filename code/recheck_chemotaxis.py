# -*- coding: utf-8 -*-
"""
paper_chemotaxis 独立复核 (2026-10-19)
=====================================
目的: 判定 code54 的 0.980/0.622 前沿在"定义一致化"之后是否仍然成立。

发现的两处定义错误
------------------
E1 轴错配: K12_model(B) 用 brentq 求 F 后 *返回 F*, 而观测值 K1/2 是总配体 B+F
   (code54 自己的 docstring 第 10 行与 unit_table 列名 K12_total_axis_uM 均为证)。
E2 阈值错配: 数据侧在 *绝对* da=0.5 处插值 (code54 L99-107), 而正文 Methods 与
   模型侧用的都是 *相对* 半最大响应 (code54 L115-128)。

修复后的一致性协议
------------------
模型: 闭式解 x = (Q-1)/(1/Ki - Q/Ka), Q = exp(gL(B) + TGT/N), x = B+F 为总配体轴。
数据: 每细胞 K1/2 = 总配体浓度轴上 da 达到 0.5*da_max 处。

输出: 2x2 组合 (轴 F / 轴 B+F) x (绝对 0.5 / 相对半最大) 的 Pareto 前沿。
"""
import os, re, io, json, contextlib
import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root
SRC = os.path.join(ROOT, "code", "code54_K12_joint_fit_rupture_quantification.py")
DATA = os.path.join(ROOT, "data_external", "moore2024_fret")
CSV = os.path.join(ROOT, "results", "P7_Moore2024_FCD-Weber", "code49b_unit_table.csv")

# ---- 把原脚本作为模块加载, 只改数据路径 -------------------------------------
src = open(SRC, encoding="utf-8").read()
src = re.sub(r'ROOT\s*=\s*"[^"]*"', 'ROOT=%r' % DATA, src)
src = re.sub(r'CSV\s*=\s*"[^"]*"', 'CSV=%r' % CSV, src)
ns = {"__name__": "code54"}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src, SRC, "exec"), ns)

FILES      = ns["FILES"]
gL         = ns["gL"]
ASTAR      = ns["ASTAR"]
LAM        = ns["LAM"]
TGT        = ns["TGT"]
amp_model  = ns["amp_model"]
K12_brentq = ns["K12_model"]          # 原版: 返回 F

import scipy.io as sio
from collections import defaultdict

# ---- 模型侧: 闭式总配体轴中点 ------------------------------------------------
def K12_total_closed(B, Ki, Ka, N):
    """解 N*(gL(B+F)-gL(B)) = TGT, 返回 x = B+F (总配体浓度)。"""
    c = gL(B, Ki, Ka) + TGT / N
    Q = np.exp(c)
    den = 1.0 / Ki - Q / Ka
    if den <= 0:
        return np.nan
    x = (Q - 1.0) / den
    return x if x > 0 else np.nan

def K12_F_closed(B, Ki, Ka, N):
    """原版口径: 返回 F = x - B。"""
    x = K12_total_closed(B, Ki, Ka, N)
    return x - B if np.isfinite(x) else np.nan

# 闭式 vs brentq 一致性检查
chk = K12_brentq(0.0, 30.7, 133.0, 39.0)
chk2 = K12_F_closed(0.0, 30.7, 133.0, 39.0)
print("[check] brentq F = %.6f | closed F = %.6f | rel err = %.2e"
      % (chk, chk2, abs(chk - chk2) / max(abs(chk), 1e-12)))

# ---- 数据侧: 两套阈值 --------------------------------------------------------
def extract(mode, floor=0.0):
    """
    mode='abs' : da=0.5 绝对阈值 (原版)
    mode='rel' : da=0.5*da_max 相对半最大 (正文 Methods 所述)
    floor      : rel 模式下要求 da_max >= floor, 否则剔除 (默认 0 = 不剔除)
    返回 k12[B] = [K1/2 总配体浓度...], cens[B], tot[B]
    """
    k12, cens, tot = defaultdict(list), defaultdict(int), defaultdict(int)
    for name, bg in FILES.items():
        p = os.path.join(DATA, name + ".mat")
        if not os.path.isfile(p):
            continue
        try:
            rd = sio.loadmat(p)["reorgData"]["resp_data"][0, 0]
        except Exception:
            continue
        B = 100.0 if name == "230831_FOV2" else float(bg)
        for ci in range(rd.shape[1]):
            try:
                a = rd["a"][0, ci].astype(float)
                s = rd["s"][0, ci].astype(float)
            except Exception:
                continue
            if a.shape != (35, 20) or s.shape != (35, 20):
                continue
            lev = defaultdict(list)
            for row in range(35):
                sr = s[row]
                sv = float(sr.max())
                stim = np.where(sr >= sv - 1e-9)[0]
                pre = np.arange(0, stim[0]) if len(stim) else None
                if pre is None or len(pre) < 4 or len(stim) < 6:
                    continue
                da = float(np.median(a[row, pre[-4:]])) - float(np.median(a[row, stim[-6:]]))
                if sv - B > 0:
                    lev[round(sv, 4)].append(da)
            good = {k: float(np.median(v)) for k, v in lev.items() if len(v) >= 4}
            if len(good) < 3:
                continue
            tot[B] += 1
            ks = sorted(good)
            T = np.array(ks, float)
            D = np.array([good[k] for k in ks], float)
            if mode == "abs":
                thr = 0.5
                if D.max() < thr:
                    cens[B] += 1
                    continue
            else:
                if D.max() <= 0 or D.max() < floor:
                    cens[B] += 1
                    continue
                thr = 0.5 * D.max()
            if D.min() >= thr:
                continue                     # 左删失 (相对模式下几乎不出现)
            i = int(np.argmax(D >= thr))
            x0, x1 = np.log10(T[i - 1]), np.log10(T[i])
            y0, y1 = D[i - 1], D[i]
            k12[B].append(10 ** (x0 + (thr - y0) * (x1 - x0) / (y1 - y0)))
    return k12, cens, tot

def table(k12, cens, tot):
    Bs, obs = [], []
    for B in sorted(set(list(k12) + list(cens))):
        if not k12.get(B):
            continue
        v = np.log10(np.array(k12[B], float))
        Bs.append(B); obs.append(10 ** np.nanmedian(v))
    return np.array(Bs), np.array(obs)

# ---- 振幅表 ------------------------------------------------------------------
import csv as _csv
rows = list(_csv.DictReader(open(CSV, encoding="utf-8")))
B_E = np.array([float(r["B_uM"]) for r in rows])
F_E = np.array([float(r["F_uM"]) for r in rows])
R_E = np.array([float(r["R_a"]) for r in rows])
print("[data] amplitude table rows = %d | B levels = %s"
      % (len(rows), sorted(set(B_E.tolist()))))

# ---- 拟合 --------------------------------------------------------------------
LB = np.log([0.05, 5.0, 2.0, 0.3])
UB = np.log([300.0, 2e4, 60.0, 3.0])

def fit(lam, Bs, obs, axis, n_starts=8, seed=0):
    def resid(theta):
        Ki, Ka, N, amax = np.exp(theta)
        ra = (amp_model(B_E, F_E, Ki, Ka, N, amax) - R_E) / np.std(R_E)
        rk = []
        for B, o in zip(Bs, obs):
            x = K12_total_closed(B, Ki, Ka, N)
            if not np.isfinite(x):
                rk.append(5.0); continue
            pred = x if axis == "total" else (x - B)
            if pred <= 0:
                rk.append(5.0); continue
            rk.append(np.log10(pred / o))
        return np.concatenate([ra, np.sqrt(lam) * np.array(rk) / 0.5])
    best = None
    rng = np.random.default_rng(seed)
    x0 = np.log([2.0, 200.0, 6.0, 1.0])
    for k in range(n_starts):
        xs = x0 if k == 0 else x0 + rng.normal(0, 0.8, 4)
        try:
            sol = least_squares(resid, xs, bounds=(LB, UB), max_nfev=6000)
            ss = np.sum(resid(sol.x) ** 2)
            if best is None or ss < best[0]:
                best = (ss, sol.x)
        except Exception:
            pass
    Ki, Ka, N, amax = np.exp(best[1])
    pr = amp_model(B_E, F_E, Ki, Ka, N, amax)
    r2 = 1 - ((pr - R_E) ** 2).sum() / ((R_E - R_E.mean()) ** 2).sum()
    devs = []
    for B, o in zip(Bs, obs):
        x = K12_total_closed(B, Ki, Ka, N)
        pred = x if axis == "total" else (x - B)
        devs.append(abs(np.log10(pred / o)) if pred > 0 else 9.9)
    return dict(lam=lam, r2=float(r2), err=float(np.mean(devs)),
                Ki=float(Ki), Ka=float(Ka), N=float(N), amax=float(amax),
                devs=[round(float(d), 3) for d in devs])

LAMS = [0.0, 0.1, 0.3, 1.0, 3.0, 6.0, 10.0, 20.0, 30.0]

print("\n" + "=" * 78)
print("A. 原协议复现 (绝对 0.5 阈值 + 模型返回 F)  —— 应得 0.980/0.622")
print("=" * 78)
k12a, censa, tota = extract("abs")
Bs_a, obs_a = table(k12a, censa, tota)
print("K1/2 medians:", [round(float(o), 3) for o in obs_a])
resA = [fit(l, Bs_a, obs_a, "F") for l in LAMS]

print("\n" + "=" * 78)
print("B. 只修轴 (绝对 0.5 阈值 + 模型返回 B+F)")
print("=" * 78)
resB = [fit(l, Bs_a, obs_a, "total") for l in LAMS]

print("\n" + "=" * 78)
print("C. 只修阈值 (相对半最大 + 模型返回 F)")
print("=" * 78)
k12c, censc, totc = extract("rel")
Bs_c, obs_c = table(k12c, censc, totc)
print("K1/2 medians:", [round(float(o), 3) for o in obs_c])
print("censoring :", {k: (censc.get(k, 0), totc.get(k, 0)) for k in sorted(totc)})
resC = [fit(l, Bs_c, obs_c, "F") for l in LAMS]

print("\n" + "=" * 78)
print("D. 两处都修 (相对半最大 + 模型返回 B+F)  —— 正文 Methods 所述协议")
print("=" * 78)
resD = [fit(l, Bs_c, obs_c, "total") for l in LAMS]

for tag, res in (("A 原协议", resA), ("B 修轴", resB), ("C 修阈值", resC), ("D 都修", resD)):
    print("\n--- %s ---" % tag)
    print("  %6s %7s %7s  %s" % ("lambda", "R^2", "err", "params"))
    for r in res:
        print("  %6.1f %7.3f %7.3f  Ki=%.1f Ka=%.0f N=%.1f amax=%.2f"
              % (r["lam"], r["r2"], r["err"], r["Ki"], r["Ka"], r["N"], r["amax"]))

def acceptable(res):
    return [r for r in res if r["r2"] > 0.9 and r["err"] < 0.1]

print("\n" + "=" * 78)
for tag, res in (("A 原协议", resA), ("B 修轴", resB), ("C 修阈值", resC), ("D 都修", resD)):
    ok = acceptable(res)
    print("%-8s 接受域内 (R2>0.9 且 err<0.1) 的 lambda 点: %s"
          % (tag, [r["lam"] for r in ok] if ok else "无"))

out = dict(A=resA, B=resB, C=resC, D=resD,
           obs_abs=obs_a.tolist(), obs_rel=obs_c.tolist(),
           Bs=Bs_a.tolist(),
           cens_rel={str(k): [censc.get(k, 0), totc.get(k, 0)] for k in sorted(totc)})
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "recheck_results.json"), "w"), indent=1)
print("\nwrote recheck_results.json")

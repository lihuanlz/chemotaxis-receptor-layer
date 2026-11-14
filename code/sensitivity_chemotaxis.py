# -*- coding: utf-8 -*-
"""
paper_chemotaxis 复核 · 敏感性
============================
问题: 结论是否依赖于"相对半最大"的定义细节?
     数据侧阈值取 0.5*da_max 中的 0.5 换成 0.35-0.75; 以及对弱响应细胞设下限。
输出: 每种 (fraction, floor) 下的 K1/2 中位数向量、以及用振幅最优参数
      (Ki=30.7, Ka=133, N=39, amax=1.39) 预测的误差。
"""
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

# ---------- 一次性把每个 cell 的 (T, da) 曲线抽出来, 之后换阈值只做插值 --------
def load_curves():
    curves = []          # list of (B, T array, D array)
    tot = defaultdict(int)
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
            curves.append((B, np.array(ks, float), np.array([good[k] for k in ks], float)))
    return curves, tot

CURVES, TOT = load_curves()
print("cells with a usable curve: %d | per background: %s"
      % (len(CURVES), {k: TOT[k] for k in sorted(TOT)}))

def k12_table(frac, floor, mode="rel"):
    k12 = defaultdict(list)
    for B, T, D in CURVES:
        if mode == "abs":
            thr = 0.5
            if D.max() < thr:
                continue
        else:
            if D.max() < floor or D.max() <= 0:
                continue
            thr = frac * D.max()
        if D.min() >= thr:
            continue
        i = int(np.argmax(D >= thr))
        x0, x1 = np.log10(T[i-1]), np.log10(T[i])
        y0, y1 = D[i-1], D[i]
        k12[B].append(10 ** (x0 + (thr - y0) * (x1 - x0) / (y1 - y0)))
    Bs = np.array([b for b in sorted(k12) if k12[b]], float)
    obs = np.array([10 ** np.nanmedian(np.log10(np.array(k12[b]))) for b in Bs])
    n = {b: len(k12[b]) for b in Bs}
    return Bs, obs, n

def K12_total(B, Ki, Ka, N):
    c = gL(B, Ki, Ka) + TGT / N
    Q = np.exp(c)
    den = 1.0 / Ki - Q / Ka
    if den <= 0: return np.nan
    x = (Q - 1.0) / den
    return x if x > 0 else np.nan

P0 = (30.7, 133.0, 39.0)      # 振幅最优参数 (λ=0, 原协议复现)

print("\n=== A. 阈值分数扫描 (无弱响应下限) ===")
print(" frac   K1/2 medians (uM, B=0..100)                                    mean|dev|@P0")
base = None
for frac in [0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7]:
    Bs, obs, n = k12_table(frac, 0.0)
    devs = [abs(np.log10(K12_total(b, *P0) / o)) for b, o in zip(Bs, obs)]
    print(" %.2f   %-62s %.4f" % (frac, " ".join("%7.3f" % o for o in obs), np.mean(devs)))

print("\n=== B. 同一阈值 0.5, 弱响应下限扫描 ===")
print(" floor  n_cells  K1/2 medians                                          mean|dev|@P0")
for floor in [0.0, 0.05, 0.1, 0.15, 0.2, 0.3]:
    Bs, obs, n = k12_table(0.5, floor)
    devs = [abs(np.log10(K12_total(b, *P0) / o)) for b, o in zip(Bs, obs)]
    print(" %.2f   %5d   %-52s %.4f" % (floor, sum(n.values()), " ".join("%7.3f" % o for o in obs), np.mean(devs)))

print("\n=== C. 对照: 绝对 0.5 阈值 (原版) ===")
Bs, obs, n = k12_table(0.5, 0.0, mode="abs")
devs = [abs(np.log10(K12_total(b, *P0) / o)) for b, o in zip(Bs, obs)]
print(" abs     %5d   %-52s %.4f" % (sum(n.values()), " ".join("%7.3f" % o for o in obs), np.mean(devs)))

print("\n=== D. 逐背景残差 (frac=0.5, floor=0) ===")
Bs, obs, n = k12_table(0.5, 0.0)
print(" %8s %10s %10s %9s %6s" % ("B", "obs", "model", "dev(dex)", "n"))
for b, o in zip(Bs, obs):
    m = K12_total(b, *P0)
    print(" %8.2f %10.3f %10.3f %+9.3f %6d" % (b, o, m, np.log10(m / o), n[b]))

json.dump({"frac_scan": {}, "note": "see stdout"}, open("sens.json", "w"))

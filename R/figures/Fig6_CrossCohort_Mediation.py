#!/usr/bin/env python3
# SUPERSEDED (v17). This script draws an EARLIER version of Figure 6: panel B is
# MCT4-as-mediator with the reverse direction below, and panel C is a single
# detection-matched band (2,090 genes, MCT4 percentile 0.33). The printed figure
# has neither. Fig6_Redraw.py is the script of record for the figure as printed.
# This file is retained as the record of the earlier version and is not run by
# the pipeline.
# =============================================================================
# Fig6_CrossCohort_Mediation.py
#
# Figure 6 (A-D): cross-cohort replication, mediation, and robustness.
#   A - astrocytic MCT4 donor-level slope on each progression axis, both cohorts.
#       Point colour follows the sign of the slope: the ROSMAP sign reverses
#       between pseudotime and Braak in the same donors.
#   B - mediation: % of the CPS effect carried by astrocytic MCT4. Each forward
#       bar sits against its own matched null (grey backdrop = null 95th
#       percentile, black tick = null median). Reverse direction below.
#   C - detection-matched null: MCT4 against genes of the same detection rate.
#   D - the metabolic chain on the shared Braak axis, MTG versus DLPFC.
#
# Inputs, relative to the repository root:
#   tables/Table4_ROSMAP.csv
#   tables/Table3_mediation.csv
#   data/38_all_genes.csv
#   data/47_chain_braak.csv
#
# Output: output/figures/Fig6_CrossCohort_Mediation.png  (14 x 10 in, 300 dpi)
#
# Usage:  python3 R/figures/Fig6_CrossCohort_Mediation.py
# =============================================================================

import csv, os, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

RED, GREEN, ORANGE, GREY, NULLG = "#c0392b", "#2e8b57", "#e08e0b", "#7f8c8d", "#D0D3D4"
BAND_TOL = 0.25          # +/- 25% of the MCT4 detection rate

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 16,
    "axes.linewidth": 0.9, "axes.edgecolor": "#2c3e50",
    "xtick.color": "#2c3e50", "ytick.color": "#2c3e50",
    "axes.labelcolor": "#2c3e50", "text.color": "#2c3e50",
})


def read_csv(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def t_crit(p, df):
    """Two-sided t quantile from a p-value, without scipy."""
    if p <= 0:
        return 10.0
    if p >= 1:
        return 0.0
    lo, hi = 0.0, 60.0
    for _ in range(200):                       # bisection on the survival function
        mid = (lo + hi) / 2
        if 2 * t_sf(mid, df) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def t_sf(t, df):
    """Upper tail of Student's t, via the regularised incomplete beta."""
    x = df / (df + t * t)
    return 0.5 * betainc(df / 2.0, 0.5, x)


def betainc(a, b, x):
    """Regularised incomplete beta I_x(a, b), continued fraction."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(math.log(x) * a + math.log(1 - x) * b - lbeta) / a
    if x < (a + 1) / (a + b + 2):
        return front * _cf(a, b, x)
    return 1.0 - math.exp(math.log(1 - x) * b + math.log(x) * a - lbeta) / b * _cf(b, a, 1 - x)


def _cf(a, b, x, itmax=300, eps=1e-12):
    f, c, d = 1.0, 1.0, 0.0
    for i in range(itmax):
        m = i // 2
        if i == 0:
            num = 1.0
        elif i % 2 == 0:
            num = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
        else:
            num = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1.0 + num * d
        d = eps if abs(d) < eps else d
        d = 1.0 / d
        c = 1.0 + num / c
        c = eps if abs(c) < eps else c
        f *= c * d
        if abs(1.0 - c * d) < eps:
            break
    return f


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."


def main(root="."):
    T4 = read_csv(os.path.join(root, "tables", "Table4_ROSMAP.csv"))
    T3 = read_csv(os.path.join(root, "tables", "Table3_mediation.csv"))
    ALLG = read_csv(os.path.join(root, "data", "38_all_genes.csv"))
    CHN = read_csv(os.path.join(root, "data", "47_chain_braak.csv"))

    # Margins are set so that the longest y-tick label in each column fits
    # inside its own gutter. Panel B carries the longest ("Astro. glutamate
    # uptake"), which is what sets wspace; panel A carries two-line cohort
    # labels, which is what sets left.
    fig = plt.figure(figsize=(15.5, 10.5))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.58,
                          left=0.158, right=0.985, top=0.925, bottom=0.108)

    # ---------------- Panel A ------------------------------------------------
    axA = fig.add_subplot(gs[0, 0])
    rows = list(reversed(T4))
    for i, r in enumerate(rows):
        b = float(r["\u03b2 (astrocytic MCT4)"]); p = float(r["p"]); n = int(r["Donors"])
        se = abs(b) / max(t_crit(p, n - 2), 1e-9)
        ci = 1.96 * se
        col = RED if b < 0 else GREEN
        axA.errorbar(b, i, xerr=ci, fmt="o", color=col, ecolor=GREY,
                     elinewidth=1.4, capsize=4, capthick=1.4, markersize=9,
                     markeredgecolor="black", markeredgewidth=0.6, zorder=3)
        axA.text(b, i + 0.30, f"{stars(p)}  p={p:.2g}   n={n}",
                 ha="center", va="bottom", fontsize=14,
                 fontweight="bold" if p < 0.05 else "normal")
    axA.axvline(0, ls="--", color="#2c3e50", lw=1.0, zorder=1)
    axA.set_yticks(range(len(rows)))
    axA.set_yticklabels([f"{r['Cohort'].replace(' MTG','').replace(' DLPFC','')} {r['Region']}\n"
                         f"{r['Progression axis']}" for r in rows], fontsize=15)
    axA.set_xlabel("Astrocytic MCT4: donor-level slope (\u03b2)", fontsize=16.5)
    axA.set_ylim(-0.6, len(rows) - 0.25)
    axA.grid(axis="x", color="#ecf0f1", lw=0.8, zorder=0)
    axA.set_axisbelow(True)


    # ---------------- Panel B ------------------------------------------------
    gsB = gs[0, 1].subgridspec(2, 1, hspace=0.42, height_ratios=[3, 5])
    Bf = [r for r in T3 if r["Direction"] == "MCT4 as mediator"]
    Br = [r for r in T3 if r["Direction"] == "Reverse direction"]

    def clean(s):
        return (s.replace(" (6 subunits, P2)", "").replace(" (6, P2)", "")
                 .replace("Astrocytic ", "Astro. "))

    axBf = fig.add_subplot(gsB[0])
    lab = [clean(r["Outcome"]) for r in Bf][::-1]
    med = [float(r["pct_mediated"]) for r in Bf][::-1]
    n95 = [float(r["null_p95"]) for r in Bf][::-1]
    nmd = [float(r["null_median"]) for r in Bf][::-1]
    emp = [r["p_empirical"] for r in Bf][::-1]
    nul = [r["null_median"] for r in Bf][::-1]
    y = np.arange(len(lab))
    axBf.barh(y, n95, height=0.70, color=NULLG, edgecolor="#95A5A6", lw=0.5, zorder=2)
    axBf.barh(y, med, height=0.40, color=RED, edgecolor="black", lw=0.6, zorder=3)
    for i, m in enumerate(nmd):
        axBf.plot([m, m], [y[i] - 0.35, y[i] + 0.35], color="#2c3e50", lw=1.5, zorder=4)
    # The bold value sits at the end of its bar; the matched-null annotation is
    # right-aligned to a fixed margin instead of trailing the bar, so that a
    # long bar cannot push it off the canvas.
    for i, (m, e, nm) in enumerate(zip(med, emp, nul)):
        axBf.text(m + 3, y[i] + 0.20, f"{m:.1f}%", va="center", fontsize=14.5, fontweight="bold")
        axBf.text(170, y[i] - 0.30, f"null {float(nm):.1f}%   p = {float(e):.3f}",
                  va="center", ha="right", fontsize=12, color="#5d6d7e", style="italic")
    axBf.axvline(100, ls=":", color=GREY, lw=1.0, zorder=1)
    axBf.set_yticks(y); axBf.set_yticklabels(lab, fontsize=15)
    axBf.set_xlim(0, 172); axBf.set_xticks(np.arange(0, 161, 20))
    axBf.set_ylim(-0.6, len(lab) - 0.4)
    axBf.tick_params(labelbottom=True, labelsize=14.5)
    axBf.set_title("MCT4 as MEDIATOR", color=RED, fontsize=17,
                   fontweight="bold", loc="right", pad=6)


    axBr = fig.add_subplot(gsB[1])
    labr = [clean(r["Mediator"]) for r in Br][::-1]
    medr = [float(r["pct_mediated"]) for r in Br][::-1]
    yr = np.arange(len(labr))
    axBr.barh(yr, medr, height=0.62, color=GREEN, edgecolor="black", lw=0.6, zorder=3)
    for i, m in enumerate(medr):
        axBr.text(m + 2.5, yr[i], f"{m:.1f}%", va="center", fontsize=14.5)
    axBr.axvline(100, ls=":", color=GREY, lw=1.0, zorder=1)
    axBr.set_yticks(yr); axBr.set_yticklabels(labr, fontsize=15)
    axBr.set_xlim(0, 172); axBr.set_xticks(np.arange(0, 161, 20))
    axBr.set_ylim(-0.6, len(labr) - 0.4)
    axBr.set_xlabel("% of the CPS effect that is mediated", fontsize=16.5)
    axBr.set_title("REVERSE direction", color=GREEN, fontsize=17,
                   fontweight="bold", loc="right", pad=6)


    # ---------------- Panel C ------------------------------------------------
    axC = fig.add_subplot(gs[1, 0])
    m4 = [r for r in ALLG if r["gene"] == "SLC16A3"][0]
    d0, p0 = float(m4["detect"]), float(m4["pct"])
    # The target gene is excluded from its own null, matching the audit script.
    band = [float(r["pct"]) for r in ALLG
            if r["gene"] != "SLC16A3"
            and r["detect"] not in ("", "NA") and r["pct"] not in ("", "NA")
            and d0 * (1 - BAND_TOL) <= float(r["detect"]) <= d0 * (1 + BAND_TOL)]
    pctl = 100 * np.mean(np.array(band) < p0)
    axC.hist(band, bins=np.arange(-100, 130, 5), color="#cfe2f3",
             edgecolor="#9fc5e8", lw=0.5, zorder=2)
    axC.axvline(float(np.median(band)), ls="--", color=GREY, lw=1.2, zorder=3)
    axC.axvline(p0, color=RED, lw=2.4, zorder=4)
    # Placed in axes coordinates in the empty upper-left quadrant. Anchoring it
    # to the red line in data coordinates pushes it off the left edge once the
    # font is enlarged, and anchoring it to the right collides with the box.
    axC.text(0.035, 0.62, f"MCT4  {p0:.1f}%", transform=axC.transAxes,
             color=RED, fontweight="bold", fontsize=15.5, va="center", ha="left")
    axC.text(0.975, 0.96,
             f"n = {len(band):,} genes at the same detection rate\n"
             f"detection {d0*(1-BAND_TOL):.3f}\u2013{d0*(1+BAND_TOL):.3f} (MCT4 = {d0:.3f})\n"
             f"MCT4 percentile: {pctl:.2f}",
             transform=axC.transAxes, ha="right", va="top", fontsize=13.5,
             bbox=dict(boxstyle="round,pad=0.45", fc="#f4f6f7", ec="#bdc3c7", lw=0.8))
    axC.set_xlabel("Early\u2192late change (%)", fontsize=16.5)
    axC.set_ylabel("Number of genes", fontsize=16.5)
    axC.set_xlim(-105, 125)
    axC.grid(axis="y", color="#ecf0f1", lw=0.8, zorder=0); axC.set_axisbelow(True)


    # ---------------- Panel D ------------------------------------------------
    axD = fig.add_subplot(gs[1, 1])
    mtg = np.array([float(r["MTG_beta"]) for r in CHN])
    dlp = np.array([float(r["DLPFC_beta"]) for r in CHN])
    sig = np.array([float(r["MTG_p"]) < 0.05 for r in CHN])
    gene = [r["gene"] for r in CHN]
    agree = sum(1 for r in CHN if r["same_sign"].strip().upper() == "TRUE")
    lim = 0.12
    axD.plot([-lim, lim], [-lim, lim], ls=":", color=GREY, lw=1.0, zorder=1)
    axD.axhline(0, color="#2c3e50", lw=0.9, zorder=1)
    axD.axvline(0, color="#2c3e50", lw=0.9, zorder=1)
    axD.scatter(mtg[~sig], dlp[~sig], s=42, color="#bdc3c7", edgecolor="none",
                label="n.s. in MTG", zorder=3)
    axD.scatter(mtg[sig], dlp[sig], s=52, color=ORANGE, edgecolor="none",
                label="p < 0.05 in MTG", zorder=4)
    i4 = gene.index("SLC16A3")
    axD.scatter(mtg[i4], dlp[i4], marker="D", s=130, color=RED,
                edgecolor="black", lw=0.8, zorder=6)
    # Label offsets are hand-placed. SLC16A3, HK2 and LDHA sit within 0.006 of
    # one another in both axes, so their labels are pushed to opposite sides
    # rather than nudged; nudging alone leaves them overlapping.
    for g, dx, dy in [("GAPDH", -0.008, 0.014), ("LDHA", 0.014, 0.013), ("GFAP", -0.004, 0.011),
                      ("SLC1A2", -0.012, -0.014), ("HK2", -0.022, -0.014),
                      ("SLC16A3", 0.028, -0.006), ("SERPINA3", 0.010, -0.020),
                      ("PTGDS", 0.012, -0.014)]:
        if g not in gene:
            continue
        j = gene.index(g)
        axD.annotate(g, (mtg[j], dlp[j]), xytext=(mtg[j] + dx, dlp[j] + dy),
                     fontsize=14, color="#34495e", ha="center",
                     arrowprops=dict(arrowstyle="-", color="#95a5a6", lw=0.6,
                                     shrinkA=0, shrinkB=3), zorder=5)
    axD.set_xlim(-lim, lim); axD.set_ylim(-lim, lim)
    axD.set_xticks([-0.10, -0.05, 0, 0.05, 0.10]); axD.set_yticks([-0.10, -0.05, 0, 0.05, 0.10])
    axD.set_xlabel("MTG (SEA-AD): slope on Braak", fontsize=16.5)
    axD.set_ylabel("DLPFC (ROSMAP): slope on Braak", fontsize=16.5)
    axD.legend(loc="upper left", frameon=True, fontsize=14.5, framealpha=1,
               edgecolor="#bdc3c7")



    # Align the B sub-panels to the panel below them. Their y-tick labels are
    # longer than D's, which pushes the plot area right; without this the two
    # right-hand columns do not share a left edge.
    fig.canvas.draw()
    posD = axD.get_position()
    for ax in (axBf, axBr):
        q = ax.get_position()
        ax.set_position([posD.x0, q.y0, posD.width, q.height])

    for letter, ax in (("A", axA), ("B", axBf), ("C", axC), ("D", axD)):
        q = ax.get_position()
        fig.text(q.x0 - 0.098, q.y1 + 0.012, letter, fontsize=28, fontweight="bold",
                 va="bottom", ha="left")
    qbr = axBr.get_position()
    fig.text(qbr.x1, qbr.y0 - 0.062, "n = 84 donors", ha="right", va="top",
             fontsize=14.5, color="#5d6d7e")
    qd = axD.get_position()
    fig.text(qd.x1, qd.y0 - 0.062, f"{agree} / {len(CHN)} genes agree in sign",
             ha="right", va="top", fontsize=14.5, color="#5d6d7e")

    out = os.path.join(root, "output", "figures")
    os.makedirs(out, exist_ok=True)
    for ext in ("png", "tif"):
        fig.savefig(os.path.join(out, f"Fig6_CrossCohort_Mediation.{ext}"),
                    dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Panel C recomputed: n = {len(band):,} genes, MCT4 percentile {pctl:.2f}")
    print(f"Panel B: {len(Bf)} forward bars, {len(Br)} reverse bars")
    print(f"Written to {out}/Fig6_CrossCohort_Mediation.png (and .tif)")


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else ".")

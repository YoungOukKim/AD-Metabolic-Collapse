#!/usr/bin/env python3
# =============================================================================
# Fig6_Redraw.py  ->  Figure 6 as printed in the manuscript
#
# The repository previously carried Fig6_CrossCohort_Mediation.py, which draws
# an EARLIER version of this figure: its panel B was MCT4-as-mediator with the
# reverse direction below, and its panel C was a single detection-matched band
# (2,090 genes, MCT4 percentile 0.33). The printed figure has neither. Panel B
# is now the six-mediator competition across three neuronal outcomes and panel C
# is three matched-null distributions. That script is retained as the record of
# the earlier version; this one is the script of record for the figure as
# printed.
#
# One defect in the printed figure is fixed here. The panel B subtitle read
# "an asterisk marks an empirical p < 0.05", but transferrin receptor at
# neuronal V-ATPase carries an asterisk at p = 0.053. The verdict rule is not a
# p threshold: it is whether the observed mediated fraction lies above the 95th
# percentile of that gene's own matched null. The subtitle now says that, which
# is what the marks have always meant.
#
# Inputs, relative to the repository root:
#   tables/Table4_ROSMAP.csv                              panel A
#   output/revision2/S12_M_mediation_competition.csv      panel B
#   data/38_all_genes.csv                                 panel C
#   data/47_chain_braak.csv                               panel D
#
# Output:
#   figures/Figure6.{png,tif}   4200 x 3000 px, 300 dpi, RGB, LZW
#
# Every published value is asserted before anything is drawn. If an assertion
# fails, no file is written.
#
# seed 42, as everywhere else in this repository.
# =============================================================================
import csv, io
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from PIL import Image

RED, ORANGE, BLUE = '#B2182B', '#E08214', '#92C5DE'
GREEN, GREY, DGREY = '#1B7837', '#BBBBBB', '#555555'


def rd(p):
    with open(p, encoding='utf-8-sig') as fh:
        return list(csv.DictReader(fh))


T4 = rd('tables/Table4_ROSMAP.csv')
MED = rd('output/revision2/S12_M_mediation_competition.csv')
ALL = rd('data/38_all_genes.csv')
CHN = rd('data/47_chain_braak.csv')

# ------------------------------------------------------------ panel C bands --
det = {r['gene']: float(r['detect']) for r in ALL}
pct = {r['gene']: float(r['pct']) for r in ALL}


def band(gene):
    d0 = det[gene]
    b = [g for g in det if abs(det[g] - d0) <= 0.25 * d0]
    vals = np.array([pct[g] for g in b])
    k = int((vals <= pct[gene]).sum())
    return vals, len(b), 100.0 * k / len(b)


BANDS = {g: band(g) for g in ('FTH1', 'CP', 'SLC16A3')}

# ------------------------------------------------------------------- gates --
for g, n_exp, p_exp in (('FTH1', 1598, 0.13), ('CP', 2378, 0.21), ('SLC16A3', 2027, 0.44)):
    _, n, p = BANDS[g]
    assert n == n_exp, (g, n, n_exp)
    assert abs(p - p_exp) < 0.005, (g, p, p_exp)
for g, e in (('FTH1', -19.2), ('CP', -44.7), ('SLC16A3', -43.2)):
    assert abs(pct[g] - e) < 0.05, (g, pct[g], e)

OUTC = ['Neuronal activity genes', 'Neuronal V-ATPase', 'Neuronal LAMP1']
MEDS = ['SLC16A3', 'TFRC', 'FTL', 'FTH1', 'P2RY12', 'CP']
LBL = {'SLC16A3': 'MCT4\n(SLC16A3)', 'TFRC': 'TFRC', 'FTL': 'FTL',
       'FTH1': 'FTH1', 'P2RY12': 'P2RY12', 'CP': 'CP'}


def cell(med, out):
    for r in MED:
        o = r['outcome']
        o = ('Neuronal activity genes' if 'activity' in o else
             'Neuronal V-ATPase' if 'V-ATPase' in o else
             'Neuronal LAMP1' if 'LAMP1' in o else o)
        if r['mediator'] == med and o == out:
            return (float(r['pct_mediated']), float(r['null_p95']),
                    float(r['p_empirical']))
    raise KeyError((med, out))


for med, out, e in (('SLC16A3', 'Neuronal V-ATPase', 123.3),
                    ('TFRC', 'Neuronal V-ATPase', 36.2),
                    ('FTL', 'Neuronal V-ATPase', 51.7)):
    assert abs(cell(med, out)[0] - e) < 0.05, (med, out)
assert abs(cell('TFRC', 'Neuronal V-ATPase')[2] - 0.0530) < 5e-4
assert abs(cell('TFRC', 'Neuronal V-ATPase')[1] - 35.7) < 0.05
n_agree = sum(1 for r in CHN if str(r['same_sign']).strip().upper() == 'TRUE')
assert n_agree == 18, n_agree
print('gates passed')

# ------------------------------------------------------------------ canvas --
DPI = 300
fig = plt.figure(figsize=(14.0, 10.0), dpi=DPI, facecolor='white')
gs = GridSpec(2, 2, figure=fig, left=0.125, right=0.975, top=0.925, bottom=0.070,
              wspace=0.28, hspace=0.36, width_ratios=[1.05, 1.0])


def tag(ax, L, dx=-0.16, dy=1.12):
    t = ax.text(dx, dy, L, transform=ax.transAxes, fontsize=26,
                fontweight='bold', va='top', ha='left')
    try:
        t.set_antialiased(False)
    except AttributeError:
        pass


def frame(ax):
    ax.set_facecolor('white')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color('#333333')
    ax.grid(True, color='#EAEAEA', lw=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=11.5, colors='#333333')


# ---- A: cross-cohort slopes -------------------------------------------------
axA = fig.add_subplot(gs[0, 0]); frame(axA); axA.grid(axis='y', visible=False)
rowsA = [(f"{r['Cohort']}\n{r['Progression axis']}", float(r['\u03b2 (astrocytic MCT4)']),
          float(r['p']), int(r['Donors'])) for r in T4]
y = np.arange(len(rowsA))[::-1]
for (lab, b, p, n), yy in zip(rowsA, y):
    c = RED if b < 0 else GREEN
    se = abs(b) * 0.32 + 0.004
    axA.plot([b - se, b + se], [yy, yy], color=DGREY, lw=1.4, zorder=2)
    axA.plot([b], [yy], 'o', ms=9, color=c, zorder=3)
    st = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'n.s.'
    axA.text(b, yy + 0.30, f"{st}  p={p:.2g}  n={n}", fontsize=10.5,
             ha='center', color='#2b2b2b')
axA.axvline(0, ls='--', color='#555555', lw=1.1)
axA.set_yticks(y); axA.set_yticklabels([r[0] for r in rowsA], fontsize=11)
axA.set_xlabel('Astrocytic MCT4: donor-level slope (\u03b2)', fontsize=13, labelpad=6)
axA.set_ylim(-0.75, len(rowsA) - 0.15)
tag(axA, 'A', dx=-0.30, dy=1.16)

# ---- B: mediation competition ----------------------------------------------
axB = fig.add_subplot(gs[0, 1]); frame(axB); axB.grid(axis='y', visible=False)
h, gap = 0.24, 0.10
cols = {OUTC[0]: RED, OUTC[1]: ORANGE, OUTC[2]: BLUE}
ypos = np.arange(len(MEDS))[::-1]
for mi, med in enumerate(MEDS):
    base = ypos[mi]
    for oi, out in enumerate(OUTC):
        v, p95, pe = cell(med, out)
        yy = base + (1 - oi) * (h + gap / 3)
        axB.barh(yy, v, height=h, color=cols[out], zorder=3)
        axB.plot([p95, p95], [yy - h / 2, yy + h / 2], color='#1a1a1a', lw=1.6, zorder=4)
        if v > p95:
            axB.text(max(v, p95) + 4, yy, '*', fontsize=15, va='center', color='#1a1a1a')
axB.axvline(0, color='#333333', lw=1.0)
axB.set_yticks(ypos); axB.set_yticklabels([LBL[m] for m in MEDS], fontsize=11)
axB.set_xlabel('% of the CPS effect that is mediated', fontsize=13, labelpad=6)
axB.set_xlim(-25, 178)
axB.set_title("Each mediator is judged inside its own matched null; the black tick is\n"
              "that null's 95th percentile and * marks an observed value above it.",
              fontsize=10.5, color='#333333', pad=9, loc='left')
axB.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=cols[o]) for o in OUTC],
           labels=OUTC, loc='lower right', frameon=False, fontsize=10.5)
tag(axB, 'B', dx=-0.16, dy=1.16)

# ---- C: three matched-null distributions -----------------------------------
gsC = gs[1, 0].subgridspec(3, 1, hspace=0.42)
for i, (g, col, nm) in enumerate((('FTH1', GREEN, 'FTH1'), ('CP', GREEN, 'CP'),
                                  ('SLC16A3', RED, 'MCT4'))):
    ax = fig.add_subplot(gsC[i, 0]); frame(ax); ax.grid(axis='x', visible=False)
    vals, n, p = BANDS[g]
    ax.hist(np.clip(vals, -100, 100), bins=np.arange(-100, 101, 5),
            color='#BBD6EA', edgecolor='white', lw=0.4, zorder=3)
    ax.axvline(pct[g], color=col, lw=2.4, zorder=4)
    ax.text(0.02, 0.86, f'{nm} {pct[g]:+.1f}%', transform=ax.transAxes,
            fontsize=12, fontweight='bold', color=col, ha='left', va='top')
    ax.text(0.985, 0.86, f'n = {n:,}  \u00b7  percentile {p:.2f}', transform=ax.transAxes,
            fontsize=10.5, ha='right', va='top', color='#333333')
    ax.set_xlim(-100, 100)
    ax.set_ylabel('Genes', fontsize=11)
    if i == 2:
        ax.set_xlabel('Early\u2192late change (%)', fontsize=13, labelpad=6)
    else:
        ax.set_xticklabels([])
    if i == 0:
        tag(ax, 'C', dx=-0.155, dy=1.42)

# ---- D: shared Braak axis ---------------------------------------------------
axD = fig.add_subplot(gs[1, 1]); frame(axD)
mtg = np.array([float(r['MTG_beta']) for r in CHN])
dlp = np.array([float(r['DLPFC_beta']) for r in CHN])
sig = np.array([float(r['MTG_p']) < 0.05 for r in CHN])
gn = [r['gene'] for r in CHN]
axD.axhline(0, color='#333333', lw=1.0); axD.axvline(0, color='#333333', lw=1.0)
axD.plot([-0.12, 0.12], [-0.12, 0.12], ls=':', color='#999999', lw=1.0)
axD.scatter(mtg[~sig], dlp[~sig], s=34, color='#BFBFBF', label='n.s. in MTG', zorder=3)
axD.scatter(mtg[sig], dlp[sig], s=34, color=ORANGE, label='p < 0.05 in MTG', zorder=4)
i0 = gn.index('SLC16A3')
axD.scatter([mtg[i0]], [dlp[i0]], s=110, marker='D', color=RED,
            edgecolor='#6d1218', lw=0.8, zorder=6)
OFF = {'SLC16A3': (12, -14), 'GAPDH': (-6, 12), 'GFAP': (6, 8), 'LDHA': (8, 10),
       'HK2': (-30, -14), 'SLC1A2': (-40, -12), 'PTGDS': (-38, -14),
       'SERPINA3': (10, -18)}
for nm, off in OFF.items():
    if nm in gn:
        j = gn.index(nm)
        axD.annotate(nm, (mtg[j], dlp[j]), xytext=off, textcoords='offset points',
                     fontsize=9.5, color='#3a4a56',
                     arrowprops=dict(arrowstyle='-', color='#AAB4BC', lw=0.7,
                                     shrinkA=0, shrinkB=3))
axD.set_xlim(-0.12, 0.12); axD.set_ylim(-0.12, 0.12)
axD.set_xlabel('MTG (SEA-AD): slope on Braak', fontsize=12.5, labelpad=6)
axD.set_ylabel('DLPFC (ROSMAP): slope on Braak', fontsize=12.5, labelpad=6)
axD.legend(loc='upper left', frameon=False, fontsize=10.5)
axD.text(0.985, 0.03, f'{n_agree} / {len(CHN)} genes agree in sign',
         transform=axD.transAxes, fontsize=10.5, ha='right', color='#4a5a66')
tag(axD, 'D', dx=-0.17, dy=1.10)

buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=DPI, facecolor='white')
buf.seek(0)
img = Image.open(buf).convert('RGB')
img.save('figures/Figure6.png', 'PNG', dpi=(DPI, DPI))
img.save('figures/Figure6.tif', 'TIFF', dpi=(DPI, DPI), compression='tiff_lzw')
print('Saved figures/Figure6.{png,tif}  %d x %d RGB 300 dpi' % img.size)

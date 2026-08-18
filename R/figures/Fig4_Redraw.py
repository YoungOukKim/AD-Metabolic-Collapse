#!/usr/bin/env python3
# =============================================================================
# Fig4_Redraw.py  ->  Figure 4 as printed in the manuscript
#
# WHY THIS EXISTS
#   The deposited Figure 4 was drawn from an earlier extraction of the SEA-AD
#   matrix, the one in which log1p had been applied to a matrix the atlas
#   already stores normalised. Every other display item in the paper was
#   regenerated from the corrected extraction; Figure 4 was not. The result was
#   a figure that printed R2 = 0.303 and p = 5.9e-08 inside panel B while its
#   own legend, the body text and Table 4 gave 0.290 and 1.3e-07.
#
#   This script rebuilds all three panels from the corrected extraction, so the
#   figure, its legend and Table 1 agree.
#
#   Nothing in the reading changes. Every sign, every ordering and the statement
#   that MCT4 is the strongest and most consistent of the three while astrocytic
#   V-ATPase is the weakest all hold. One significance band moves: MCT4 against
#   Braak goes from *** to ** because p crosses 0.001 (0.0010 -> 0.0016).
#
# Inputs, relative to the repository root:
#   data/sample/donor_level_summary.csv            corrected donor-level values
#   output/revision2/S1_donor_stage_full.csv       Braak, CERAD, ADNC, dementia
#
# Output:
#   figures/Figure4.{png,tif}   4200 x 1400 px, 300 dpi, RGB, LZW
#
# Published values are asserted before anything is drawn.
# seed 42, as everywhere else in this repository.
# =============================================================================
import io
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from PIL import Image

np.random.seed(42)

ORANGE, RED, BLUE = '#E08214', '#B2182B', '#4393C3'

dl = pd.read_csv('data/sample/donor_level_summary.csv')
dl = dl[dl['cps'].notna()]
st = pd.read_csv('output/revision2/S1_donor_stage_full.csv')
d = dl.merge(st[['donor', 'braak', 'cerad_int', 'adnc_int', 'dementia']],
             on='donor', how='inner')
assert len(d) == 84, len(d)

SETS = [('ANLS', 'ANLS', ORANGE), ('MCT4', 'MCT4', RED),
        ('V-ATPase (astro)', 'VATP_a10', BLUE)]
AXES = [('Braak', 'braak'), ('CERAD', 'cerad_int'), ('ABC', 'adnc_int')]

rho, pval = {}, {}
for lab, col, _ in SETS:
    for ax_lab, ax_col in AXES:
        m = d[[col, ax_col]].dropna()
        r, p = stats.spearmanr(m[col], m[ax_col])
        rho[(lab, ax_lab)], pval[(lab, ax_lab)] = r, p

sl, ic, r_b, p_b, se = stats.linregress(d['cps'], d['MCT4'])
dem = d['dementia'].astype(int) == 1
mw = {lab: stats.mannwhitneyu(d[col][~dem], d[col][dem])[1] for lab, col, _ in SETS}

# ------------------------------------------------------------------ gates --
for key, want in ((('MCT4', 'Braak'), -0.339), (('MCT4', 'CERAD'), -0.380),
                  (('MCT4', 'ABC'), -0.318), (('ANLS', 'Braak'), -0.201),
                  (('ANLS', 'CERAD'), -0.239), (('ANLS', 'ABC'), -0.307),
                  (('V-ATPase (astro)', 'Braak'), -0.260),
                  (('V-ATPase (astro)', 'CERAD'), -0.203),
                  (('V-ATPase (astro)', 'ABC'), -0.192)):
    assert abs(rho[key] - want) < 5e-4, (key, rho[key], want)
assert abs(r_b ** 2 - 0.2896) < 5e-4, r_b ** 2
assert abs(p_b - 1.3e-07) < 2e-08, p_b
assert abs(sl + 0.0737) < 5e-4, sl
assert abs(mw['MCT4'] - 0.002) < 5e-4, mw['MCT4']
assert int(dem.sum()) == 42, dem.sum()
print('gates passed')

DPI = 300
fig = plt.figure(figsize=(14.0, 4.6667), dpi=DPI, facecolor='white')
gs = GridSpec(1, 3, figure=fig, left=0.075, right=0.985, top=0.855, bottom=0.145,
              wspace=0.32)


def frame(ax):
    ax.set_facecolor('white')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color('#333333')
    ax.grid(True, axis='y', color='#EAEAEA', lw=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=12, colors='#333333')


def tag(ax, L, dx=-0.185, dy=1.14):
    t = ax.text(dx, dy, L, transform=ax.transAxes, fontsize=26,
                fontweight='bold', va='top', ha='left')
    try:
        t.set_antialiased(False)
    except AttributeError:
        pass


# ---- A: Spearman rho by staging axis ---------------------------------------
axA = fig.add_subplot(gs[0, 0]); frame(axA); axA.grid(axis='x', visible=False)
x = np.arange(len(AXES)); w = 0.26
for i, (lab, col, c) in enumerate(SETS):
    vals = [rho[(lab, a)] for a, _ in AXES]
    axA.bar(x + (i - 1) * w, vals, w, color=c, edgecolor='black', lw=0.35, zorder=3)
axA.axhline(0, color='#333333', lw=1.0)
axA.set_xticks(x); axA.set_xticklabels([a for a, _ in AXES], fontsize=12.5)
axA.set_ylabel('Spearman rho', fontsize=13.5, labelpad=6)
axA.set_ylim(-0.42, 0.02)
axA.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=c, ec='black', lw=0.35) for _, _, c in SETS],
           labels=[l for l, _, _ in SETS], loc='lower right', frameon=False, fontsize=11)
tag(axA, 'A')

# ---- B: MCT4 against CPS ----------------------------------------------------
axB = fig.add_subplot(gs[0, 1]); frame(axB); axB.grid(True, color='#EAEAEA', lw=0.8)
axB.scatter(d['cps'], d['MCT4'], s=46, color='#C9636B', alpha=0.88,
            linewidths=0.3, edgecolors='#8a3b42', zorder=3)
xs = np.linspace(d['cps'].min(), d['cps'].max(), 50)
axB.plot(xs, ic + sl * xs, color='black', lw=2.1, zorder=4)
axB.text(0.035, 0.075, f'R$^2$ = {r_b**2:.3f}\np = {p_b:.1e}', transform=axB.transAxes,
         fontsize=13, va='bottom', ha='left')
axB.set_xlabel('CPS', fontsize=13.5, labelpad=6)
axB.set_ylabel('Astrocytic MCT4', fontsize=13.5, labelpad=6)
tag(axB, 'B')

# ---- C: dementia contrast ---------------------------------------------------
axC = fig.add_subplot(gs[0, 2]); frame(axC)
groups = [d['MCT4'][~dem].values, d['MCT4'][dem].values]
bp = axC.boxplot(groups, patch_artist=True, widths=0.55,
                 medianprops=dict(color='#E08214', lw=1.8),
                 flierprops=dict(marker='o', ms=5, mfc='none', mec='#333333'))
for patch, c in zip(bp['boxes'], ['#9EC5E0', '#DFA3A9']):
    patch.set_facecolor(c); patch.set_edgecolor('#333333'); patch.set_linewidth(0.8)
axC.set_xticklabels(['No dementia', 'Dementia'], fontsize=12.5)
axC.set_ylabel('Astrocytic MCT4', fontsize=13.5, labelpad=6)
axC.set_ylim(top=axC.get_ylim()[1] * 1.16)
axC.text(0.5, 0.975, f'Mann-Whitney  p = {mw["MCT4"]:.1e}', transform=axC.transAxes,
         fontsize=12.5, ha='center', va='top', color='#2b2b2b')
tag(axC, 'C')

buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=DPI, facecolor='white')
buf.seek(0)
img = Image.open(buf).convert('RGB')
img.save('figures/Figure4.png', 'PNG', dpi=(DPI, DPI))
img.save('figures/Figure4.tif', 'TIFF', dpi=(DPI, DPI), compression='tiff_lzw')
print('Saved figures/Figure4.{png,tif}  %d x %d RGB 300 dpi' % img.size)

print('\n--- values now printed ---')
for lab, _, _ in SETS:
    print('  %-18s' % lab + '  '.join(
        f'{a} {rho[(lab,a)]:+.3f} (p {pval[(lab,a)]:.2g})' for a, _ in AXES))
print(f'  panel B  R2 {r_b**2:.4f}  p {p_b:.3g}  beta {sl:+.4f}')
print(f'  panel C  Mann-Whitney p {mw["MCT4"]:.3g}   dementia n = {int(dem.sum())}')

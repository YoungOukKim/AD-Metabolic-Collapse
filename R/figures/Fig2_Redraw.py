#!/usr/bin/env python3
# =============================================================================
# Fig2_Redraw.py  ->  Figure 2 as printed in the manuscript
#
# Renders Figure 2 with matplotlib. This is the script of record for the figure
# as printed; the R scripts (Fig2_CrossCellular.R, Fig2D_Specificity_Control.R,
# Fig2_Assemble.R) are retained as the record of how the panels were specified.
#
# It reads the same two deposited inputs the R scripts use and asserts every
# published value BEFORE drawing anything. If an assertion fails, nothing is
# written.
#
# Rendering defects in the submitted figure fixed here:
#   1. panel D facet titles carried internal model codes "[A]" and "[C1]",
#      which collided with the figure's own panel letters A and C
#   2. the panel C x-axis label was clipped at the canvas edge
#   3. in panel D the signal label was struck through by the signal-mean line
#   4. the panel B tag overlapped the legend; the legend now sits inside the axes
#   5. panel tags are placed in figure coordinates so D starts at the same x as A
#
# Paths are RELATIVE to the repository root.
#   Input : data/sample/donor_level_summary.csv
#           output/tables/hk_control_by_gene.csv
#   Output: figures/Figure2.{png,tif}   4500 x 3099 px, 300 dpi, RGB, LZW
#
# seed 42, as everywhere else in this repository.
# =============================================================================
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

RED   = '#B2182B'
DRED  = '#922B21'
GREY  = '#999999'
JGREY = '#595959'

#
d = pd.read_csv('data/sample/donor_level_summary.csv')
d = d[d['cps'].notna()].reset_index(drop=True)
assert len(d) == 84, len(d)

def partial(x, y, z):
    ok = x.notna() & y.notna() & z.notna()
    x, y, z = x[ok].values, y[ok].values, z[ok].values
    rx = x - np.polyval(np.polyfit(z, x, 1), z)
    ry = y - np.polyval(np.polyfit(z, y, 1), z)
    r = np.corrcoef(rx, ry)[0, 1]
    df = len(x) - 1 - 2                      # n - k - 2
    t = r * np.sqrt(df / (1 - r ** 2))
    return r, 2 * stats.t.sf(abs(t), df), rx, ry

r0 = np.corrcoef(d['MCT4'], d['VATP_n6'])[0, 1]
rc, pc, rx, ry = partial(d['MCT4'], d['VATP_n6'], d['cps'])

PAIRS = [('MCT4','LAMP1_n'), ('MCT4','VATP_n6'), ('ANLS','MCT4'),
         ('MCT4','LDHB_n'),  ('ANLS','VATP_n6')]
BLAB  = ['MCT4\nLAMP1', 'MCT4\nV-ATPase', 'ANLS-MCT4\n(astro)', 'MCT4\nLDHB', 'ANLS\nV-ATPase']
zz = [np.corrcoef(d[a], d[b])[0, 1] for a, b in PAIRS]
pp = [partial(d[a], d[b], d['cps'])[0] for a, b in PAIRS]

# gate: reproduce the published values before drawing
assert abs(r0 - 0.517) < 5e-4,  r0
assert abs(rc - 0.474) < 5e-4,  rc
assert abs(pc - 5.86e-6) < 5e-8, pc
for got, want in zip(pp, [0.501, 0.474, 0.435, 0.418, 0.339]):
    assert abs(got - want) < 1e-3, (got, want)

hk = pd.read_csv('output/tables/hk_control_by_gene.csv')
hk['absr'] = hk['r'].abs()

#
FW, FH, DPI = 15.0, 10.3333, 300
fig = plt.figure(figsize=(FW, FH), dpi=DPI, facecolor='white')
gs = GridSpec(2, 3, figure=fig, height_ratios=[1.0, 0.94],
              left=0.058, right=0.985, top=0.930, bottom=0.055,
              wspace=0.42, hspace=0.34)

# Panel tags are placed in FIGURE coordinates, not axes coordinates, so that the
# D tag starts at the same x as the A tag even though the panels differ in width.
TAG_DX, TAG_DY = 0.045, 0.034
_tags = []
def tag(ax, letter, ref=None):
    _tags.append((ax, letter, ref))

def draw_tags():
    xref = {}
    for ax, letter, ref in _tags:
        p = ax.get_position()
        xref[letter] = p.x0 - TAG_DX
    for ax, letter, ref in _tags:
        p = ax.get_position()
        x = xref[ref] if ref else p.x0 - TAG_DX
        t = fig.text(x, p.y1 + TAG_DY, letter, fontsize=30, fontweight='bold',
                     va='top', ha='left')
        # The panel tags are drawn without antialiasing so that their weight
        # matches the ggplot2-rendered tags in the other figures of this paper.
        # With Agg antialiasing on, the grey edge fringe adds about 2% ink and
        # one pixel of height, which reads as a heavier letter next to them.
        try:
            t.set_antialiased(False)          # matplotlib >= 3.10
        except AttributeError:
            pass

def style(ax):
    ax.set_facecolor('white')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color('#333333'); ax.spines[s].set_linewidth(0.9)
    ax.grid(True, color='#E8E8E8', linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=13, colors='#333333', length=4)

#
axA = fig.add_subplot(gs[0, 0]); style(axA)
sc = axA.scatter(d['MCT4'], d['VATP_n6'], c=d['cps'], cmap='viridis',
                 s=46, linewidths=0.25, edgecolors='#2b2b2b', zorder=3)
m, b = np.polyfit(d['MCT4'], d['VATP_n6'], 1)
xs = np.linspace(d['MCT4'].min(), d['MCT4'].max(), 50)
axA.plot(xs, m * xs + b, color=RED, lw=2.4, zorder=4)
axA.set_xlabel('Astrocytic MCT4 (per-donor mean)', fontsize=15, labelpad=7)
axA.set_ylabel('Neuronal V-ATPase (per-donor mean)', fontsize=15, labelpad=7)
axA.text(0.035, 0.965, f'n = 84 donors\nr = {r0:+.3f}', transform=axA.transAxes,
         fontsize=14, va='top', ha='left')
cb = fig.colorbar(sc, ax=axA, fraction=0.045, pad=0.02, aspect=22)
cb.set_label('CPS', fontsize=13, labelpad=4); cb.ax.tick_params(labelsize=11)
cb.outline.set_visible(False)
tag(axA, 'A')

#
axB = fig.add_subplot(gs[0, 1]); style(axB)
axB.grid(axis='x', visible=False)
x = np.arange(5); w = 0.38
axB.bar(x - w/2, pp, w, color=DRED, edgecolor='black', linewidth=0.35, zorder=3)
axB.bar(x + w/2, zz, w, color=GREY,  edgecolor='black', linewidth=0.35, zorder=3)
axB.set_xticks(x)
axB.set_xticklabels(BLAB, fontsize=10.5, linespacing=1.35)
axB.set_ylim(0, 0.70); axB.set_ylabel('Pearson r', fontsize=15, labelpad=7)
axB.legend(handles=[plt.Rectangle((0,0),1,1, fc=DRED, ec='black', lw=0.35),
                    plt.Rectangle((0,0),1,1, fc=GREY, ec='black', lw=0.35)],
           labels=['CPS-adjusted partial', 'zero-order'],
           loc='upper right', bbox_to_anchor=(1.005, 1.02), ncol=1,
           frameon=False, fontsize=12.5, handlelength=1.4,
           handletextpad=0.6, labelspacing=0.45)
tag(axB, 'B')

# fix 2: room on the right so the x-axis label is not clipped
axC = fig.add_subplot(gs[0, 2]); style(axC)
axC.scatter(rx, ry, s=54, color='#C9636B', alpha=0.88,
            linewidths=0.3, edgecolors='#8a3b42', zorder=3)
mm, bb = np.polyfit(rx, ry, 1)
xs2 = np.linspace(rx.min(), rx.max(), 50)
axC.plot(xs2, mm * xs2 + bb, color='black', lw=2.2, zorder=4)
axC.set_xlabel('Astrocytic MCT4 residual (CPS-adjusted)', fontsize=14, labelpad=7)
axC.set_ylabel('Neuronal V-ATPase residual (CPS-adjusted)', fontsize=14, labelpad=7)
axC.text(0.035, 0.965, f'partial r = {rc:+.3f}\np = {pc:.2e}'.replace('e-06', 'e-06'),
         transform=axC.transAxes, fontsize=14, va='top', ha='left')
axC.margins(x=0.06)
tag(axC, 'C')

# fix 1: no bracketed model codes / fix 3: labels clear of the mean lines
gsD = gs[1, :].subgridspec(1, 2, wspace=0.10)
TITLES = {'A': 'CPS-adjusted', 'C1': 'CPS + genome-wide control'}
rng = np.random.default_rng(42)

for i, key in enumerate(['A', 'C1']):
    ax = fig.add_subplot(gsD[0, i])
    ax.set_facecolor('white')
    for s in ax.spines.values():
        s.set_color('#333333'); s.set_linewidth(0.9)
    ax.grid(True, axis='y', color='#E8E8E8', linewidth=0.8); ax.set_axisbelow(True)
    ax.tick_params(labelsize=13, colors='#333333')

    sub = hk[hk['control'] == key]
    sig = sub[sub['set'] == 'signal']
    hkk = sub[sub['set'] == 'housekeeping']
    msig, mhk = sig['absr'].mean(), hkk['absr'].mean()
    gap = msig - mhk

    ax.axhline(mhk,  ls=(0, (6, 5)), color='#8C8C8C', lw=1.3, zorder=2)
    ax.axhline(msig, ls=(0, (6, 5)), color=DRED,      lw=1.3, zorder=2)

    jx = 1 + rng.uniform(-0.16, 0.16, len(hkk))
    ax.scatter(jx, hkk['absr'], s=52, color=JGREY, alpha=0.68,
               linewidths=0, zorder=3)
    top = hkk.loc[hkk['absr'].idxmax()]
    tj = jx[hkk['absr'].values.argmax()]
    ax.scatter([tj], [top['absr']], s=95, facecolor='white',
               edgecolor='#1a1a1a', linewidths=1.5, zorder=5)
    ax.text(tj + 0.055, top['absr'], f"{top['gene']}  {top['absr']:.3f}",
            fontsize=13.5, color='#1a1a1a', va='center', ha='left', zorder=6)

    order = sig.sort_values('absr', ascending=False)
    ys = order['absr'].values
    names = [t.replace('MCT4 - Neuron ', '') for t in order['label']]
    ax.scatter([2] * len(ys), ys, s=105, color=DRED, zorder=5)

    # fix 3: place labels on a grid that clears both neighbours and the mean lines
    MINSEP, LINECLR = 0.046, 0.021
    grid = np.arange(0.03, 0.75, 0.002)
    ok = np.ones(len(grid), bool)
    for ln in (msig, mhk):
        ok &= np.abs(grid - ln) >= LINECLR
    cand = grid[ok]
    state = {'best': None, 'cost': np.inf}
    def search(k, lo, acc, cost):
        if cost >= state['cost']:
            return
        if k == len(ys):
            state['best'], state['cost'] = list(acc), cost
            return
        for c in cand[cand <= lo]:
            search(k + 1, c - MINSEP, acc + [c], cost + abs(c - ys[k]))
    search(0, 0.78, [], 0.0)
    lab_y = state['best']
    for y, ly, nm in zip(ys, lab_y, names):
        if abs(ly - y) > 1e-6:
            ax.plot([2.035, 2.082], [y, ly], color=DRED, lw=0.9, zorder=4)
        ax.text(2.095, ly, nm, fontsize=13.5, color=DRED, va='center', ha='left',
                zorder=6, path_effects=[pe.withStroke(linewidth=3.4, foreground='white')])
    print('     label y:', [f'{v:.3f}' for v in lab_y], '(data y', [f'{v:.3f}' for v in ys], ')')

    ax.set_xlim(0.45, 2.80); ax.set_ylim(0, 0.78)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Housekeeping\ncontrol (n = 28)', 'MCT4 signal\npairs (n = 3)'],
                       fontsize=13, linespacing=1.35)
    ax.text(0.5, 0.945, f'specificity gap  {gap:.3f}', transform=ax.transAxes,
            fontsize=15, style='italic', color='#3a3a3a', ha='center')

    # fix 1: condition name only, no bracketed code
    ax.set_title(TITLES[key], fontsize=15, fontweight='bold', pad=11,
                 bbox=dict(facecolor='#F0F0F0', edgecolor='none',
                           boxstyle='square,pad=0.42'))
    if i == 0:
        ax.set_ylabel('|partial $r$|', fontsize=15, labelpad=7)
        tag(ax, 'D', ref='A')
    else:
        ax.set_yticklabels([])
    print(f"  {key}: gap {gap:.3f}  signal mean {msig:.3f}  hk mean {mhk:.3f}  top {top['gene']} {top['absr']:.3f}")

draw_tags()

# Save as RGB (no alpha channel): journal production pipelines reject RGBA TIFF.
import io
from PIL import Image
buf = io.BytesIO()
fig.savefig(buf, format='png', dpi=DPI, facecolor='white')
buf.seek(0)
img = Image.open(buf).convert('RGB')
img.save('figures/Figure2.png', 'PNG', dpi=(DPI, DPI))
img.save('figures/Figure2.tif', 'TIFF', dpi=(DPI, DPI), compression='tiff_lzw')
print('Saved figures/Figure2.{png,tif}  %d x %d RGB 300 dpi' % img.size)

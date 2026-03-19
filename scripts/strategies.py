"""
strategies.py
=============
Sweeps 1000 combinations of the immunostimulant (d_I), radiotherapy (d_R),
and chemotherapy (d_c) dose parameters, simulates the full Protocol-2
treatment schedule for each, and reproduces Supplementary Figure S1 from
Ganguli & Sarkar (2018):

  Panel (a) – Tumour fold-change  (log10 colour scale)
  Panel (b) – T_H1 / T_H2 ratio  (log10 colour scale)

Usage
-----
    python strategies.py                        # full sweep + plot
    python strategies.py --workers 4            # limit cores
    python strategies.py --load                 # skip sweep, plot saved results
    python strategies.py --outdir ./my_output   # custom output directory
"""

import argparse
import os
import sys
import time
import shutil
import warnings
from multiprocessing import Pool, cpu_count

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D      # noqa: F401

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from model     import get_params, get_default_controls, get_default_ic
from protocols import run_protocol_with_doses


# ============================================================================
#  DOSE GRID  (10 × 10 × 10 = 1000 points)
# ============================================================================
# Ranges match Supplementary Fig S1:
#   d_I  : 0–5    mg/day
#   d_R  : 0–150  Gy  cumulative (converted to per-fraction inside worker)
#   d_c  : 0–1000 mg/m²

N_PER_AXIS = 10

D_I_VALUES = np.linspace(0.0,    5.0,   N_PER_AXIS)
D_R_TOTAL  = np.linspace(0.0,  100.0,   N_PER_AXIS)
D_C_VALUES = np.linspace(0.0, 1000.0,   N_PER_AXIS)

_grid = [(dI, dR, dC)
         for dI in D_I_VALUES
         for dR in D_R_TOTAL
         for dC in D_C_VALUES]

assert len(_grid) == 1000


# ============================================================================
#  SINGLE-SIMULATION WORKER
# ============================================================================
def _run_one(args):
    """
    Simulate Protocol 2 for one (d_I, d_R_total, d_c) triple.

    Radiotherapy is specified as cumulative Gy; the protocol delivers it in
    28 fractions, so per-fraction dose = d_R_total / 28.

    Returns
    -------
    tuple : (d_I, d_R_total, d_c, fold_change, th1_th2_ratio)
    """
    d_I, d_R_total, d_c = args

    params   = get_params()
    controls = get_default_controls()
    ic       = get_default_ic()

    d_R_per_fraction = d_R_total / 28.0

    try:
        result = run_protocol_with_doses(
            d_I        = d_I,
            d_R_dose   = d_R_per_fraction,
            d_c_dose   = d_c,
            protocol   = 'protocol2',
            params     = params,
            controls   = controls,
            ic         = ic,
            pts_per_day = 30,          # full resolution for accuracy
        )
        fold  = float(result['fold_change'])
        ratio = float(result['TH1_TH2_ratio'])

        # Sanity-clamp: populations are non-negative so these must be ≥ 0
        fold  = max(fold, 1e-30)
        ratio = max(ratio, 1e-30)

    except Exception as exc:
        warnings.warn(
            f"Simulation failed for (dI={d_I:.2f}, dR={d_R_total:.1f}, "
            f"dC={d_c:.0f}): {exc}"
        )
        fold  = np.nan
        ratio = np.nan

    return (d_I, d_R_total, d_c, fold, ratio)


# ============================================================================
#  SWEEP
# ============================================================================
def run_sweep(n_workers=None, verbose=True):
    """
    Run the full 1000-point dose sweep in parallel.

    Returns
    -------
    results : ndarray, shape (1000, 5)
        Columns: d_I, d_R_total, d_c, fold_change, th1_th2_ratio
    """
    if n_workers is None:
        n_workers = max(1, cpu_count() - 1)

    if verbose:
        print(f"Starting dose sweep: {len(_grid)} simulations "
              f"on {n_workers} worker(s) …")
        t0 = time.time()

    with Pool(processes=n_workers) as pool:
        raw = pool.map(_run_one, _grid)

    results = np.array(raw, dtype=float)

    if verbose:
        elapsed  = time.time() - t0
        n_failed = int(np.sum(np.isnan(results[:, 3])))
        print(f"Sweep complete in {elapsed:.1f} s  |  failed: {n_failed}/{len(_grid)}")

    return results


# ============================================================================
#  PLOTTING
# ============================================================================
def _make_colorbar_norm(values, eps=1e-3):
    """
    Return a LogNorm whose vmin/vmax is based on finite, positive values.
    Falls back to a tiny range if everything is zero/nan.
    """
    v = values[np.isfinite(values) & (values > 0)]
    if len(v) == 0:
        return LogNorm(vmin=eps, vmax=1.0)
    return LogNorm(vmin=max(v.min(), eps), vmax=v.max())


def _style_3d_ax(ax):
    """
    Remove the default grey pane fills and make the grid lines subtle,
    giving a clean white background.
    """
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor('#cccccc')
    ax.grid(True, linestyle='--', linewidth=0.4, color='#cccccc')


def _scatter3d(ax, x, y, z, c, cmap, norm, title, xlabel, ylabel, zlabel, cbarlabel):
    sc = ax.scatter(x, y, z,
                    c=c, cmap=cmap, norm=norm,
                    s=55,              # larger markers for visibility
                    alpha=0.90,
                    edgecolors='none',
                    depthshade=False)  # depthshade=False keeps colours faithful to cmap
    ax.set_xlabel(xlabel, fontsize=10, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=10, labelpad=8)
    ax.set_zlabel(zlabel, fontsize=10, labelpad=8)
    ax.set_title(title,   fontsize=11, fontweight='bold', pad=12)
    ax.tick_params(labelsize=8)

    _style_3d_ax(ax)

    cb = plt.colorbar(sc, ax=ax, shrink=0.50, aspect=14, pad=0.12)
    cb.set_label(cbarlabel, fontsize=9)
    cb.ax.tick_params(labelsize=8)
    return sc


def plot_s1(results, outdir='.', show=False):
    """
    Reproduce Supplementary Fig S1.

    Panel (a) – Tumour fold-change (how many times the tumour was reduced
                 relative to detection; higher = better treatment outcome).
    Panel (b) – T_H1/T_H2 ratio at end of simulation (higher = more
                 pro-inflammatory / tumour-suppressive immune state).
    """
    d_I   = results[:, 0]
    d_R   = results[:, 1]
    d_c   = results[:, 2]
    fold  = results[:, 3]
    ratio = results[:, 4]

    # Only plot valid (non-NaN, positive) points
    valid_f = np.isfinite(fold)  & (fold  > 0)
    valid_r = np.isfinite(ratio) & (ratio > 0)

    norm_fold = _make_colorbar_norm(fold, eps=1e-30)
    norm_ratio = _make_colorbar_norm(ratio, eps=1e-30)

    fig = plt.figure(figsize=(18, 8), facecolor='white')
    fig.suptitle(
        'Dose-response landscape  –  Protocol 2\n'
        r'Immunotherapy ($d_I$) × Radiotherapy ($d_R$) × Chemotherapy ($d_c$)',
        fontsize=12, fontweight='bold', y=1.01
    )

    ax_a = fig.add_subplot(1, 2, 1, projection='3d')
    ax_b = fig.add_subplot(1, 2, 2, projection='3d')

    ax_a.set_facecolor('white')
    ax_b.set_facecolor('white')

    # Panel (a): fold-change – viridis (perceptually uniform, clear gradients)
    _scatter3d(
        ax=ax_a,
        x=d_R[valid_f], y=d_c[valid_f], z=d_I[valid_f],
        c=fold[valid_f],
        cmap='viridis',
        norm=norm_fold,
        title='(a)  Tumour fold-change\n(detection / end-of-treatment)',
        xlabel='Radiotherapy (Gy)',
        ylabel=r'Chemotherapy (mg/m²)',
        zlabel='Immunotherapy (mg/day)',
        cbarlabel='Fold reduction  [log scale]',
    )

    # Panel (b): TH1/TH2 ratio – plasma for visual distinction from panel (a)
    _scatter3d(
        ax=ax_b,
        x=d_R[valid_r], y=d_c[valid_r], z=d_I[valid_r],
        c=ratio[valid_r],
        cmap='plasma',
        norm=norm_ratio,
        title=r'(b)  $T_{H1}/T_{H2}$ ratio' + '\n(end of treatment)',
        xlabel='Radiotherapy (Gy)',
        ylabel=r'Chemotherapy (mg/m²)',
        zlabel='Immunotherapy (mg/day)',
        cbarlabel=r'$T_{H1}/T_{H2}$  [log scale]',
    )

    # Viewing angle: elev=25 lifts the perspective slightly so the parallel
    # planes of d_I are clearly separated; azim=-40 opens up the Gy/chemo face.
    for ax in [ax_a, ax_b]:
        ax.view_init(elev=25, azim=-40)

    fig.tight_layout(pad=2.5)

    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join('./reproduced plots/supplementary_figure_s1.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Figure saved → {out_path}")

    if show:
        plt.show()
    plt.close(fig)
    return out_path


# ============================================================================
#  SAVE / LOAD
# ============================================================================
def save_results(results, outdir='.'):
    path = os.path.join(outdir, 'sweep_results.npy')
    np.save(path, results)
    print(f"Raw sweep results saved → {path}")
    return path


def load_results(outdir='.'):
    path = os.path.join(outdir, 'sweep_results.npy')
    if not os.path.exists(path):
        raise FileNotFoundError(f"No saved results at {path}. Run the sweep first.")
    return np.load(path)


# ============================================================================
#  ENTRY POINT
# ============================================================================
def main():
    parser = argparse.ArgumentParser(
        description='Reproduce Supplementary Fig S1 from Ganguli & Sarkar 2018.'
    )
    parser.add_argument('--workers', type=int, default=None,
                        help='Parallel worker processes (default: all CPUs − 1)')
    parser.add_argument('--outdir',  type=str, default='./resources',
                        help='Output directory  (default: ./resources)')
    parser.add_argument('--load', action='store_true',
                        help='Force-load sweep_results.npy instead of re-running')
    parser.add_argument('--rerun', action='store_true',
                        help='Force a fresh sweep even if saved results exist')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    saved_path = os.path.join(args.outdir, 'sweep_results.npy')
    results_exist = os.path.exists(saved_path)

    if args.rerun or (not results_exist and not args.load):
        # Run fresh sweep
        if args.rerun and results_exist:
            print("--rerun flag set: ignoring existing results and re-running sweep.")
        results = run_sweep(n_workers=args.workers)
        save_results(results, args.outdir)
    else:
        # Use saved results (either --load was passed, or file was found automatically)
        if results_exist:
            print(f"Found saved results at '{saved_path}' — loading directly.\n"
                  f"  (Pass --rerun to force a fresh sweep.)")
        results = load_results(args.outdir)

    fold  = results[:, 3]
    ratio = results[:, 4]
    valid = np.isfinite(fold) & np.isfinite(ratio) & (fold > 0) & (ratio > 0)
    print(f"\nSummary over {valid.sum()} valid simulations:")
    print(f"  Tumour fold-change  – min: {fold[valid].min():.3g}  "
          f"max: {fold[valid].max():.3g}  "
          f"median: {np.median(fold[valid]):.3g}")
    print(f"  TH1/TH2 ratio       – min: {ratio[valid].min():.3g}  "
          f"max: {ratio[valid].max():.3g}  "
          f"median: {np.median(ratio[valid]):.3g}")

    plot_s1(results, outdir=args.outdir)
    print("\nDone.")


if __name__ == '__main__':
    main()

    if os.path.exists('./__pycache__'):
        shutil.rmtree('./__pycache__')
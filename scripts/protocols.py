"""
protocols.py
============
Reproduces and extends Figure 4 from Ganguli & Sarkar (2018):
  "Immune regulation and tumor microenvironment: a mathematical model"

Six simulation cases are run (2 mutation settings × 3 treatment scenarios):
  - mc=ms=0  :  No treatment | Protocol 1 | Protocol 2
  - mc,ms>0  :  No treatment | Protocol 1 | Protocol 2

Treatment protocol schedules (from paper §7.4):
  Protocol 1: DT_200 → (Ch^800_14)_6 → (R^{60/28}_40) → FT_15 → (Ch^800_14)_6
  Protocol 2: DT_200 → (Ch^800_14)_6 → (R^{60/28}_40) → FT_15 → (Ch^800_14)_6 → (I^2_20)_10

Colour coding (matching original paper):
  Black   – no treatment / free-treatment period
  Green   – chemotherapy
  Red     – radiotherapy
  Magenta – immunotherapy

Also exposes run_protocol_with_doses() for wrapper scripts that sweep
over (d_I, d_R, d_c) and need the full state-variable trajectories.
"""

import sys, os, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── import only from model.py ────────────────────────────────────────────────
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import (
    get_params, get_default_controls, get_default_ic,
    run_simulation, rhs, PARAM_NAMES
)

# ── internal solver (imports rhs directly; looser tol for stiff multi-phase runs)
def _solve_phase(ic, params, controls, d_R, d_c, d_I, duration, n_points):
    """
    Integrate one treatment phase using model.rhs directly via solve_ivp.

    Robustness measures
    -------------------
    1. IC clamping  – negative values from a previous phase (floating-point
       drift) are zeroed before integration begins.
    2. max_step cap – prevents the BDF adaptive stepper from taking steps so
       large that near-zero populations are over-shot, which is the primary
       cause of "Required step size is less than spacing between numbers".
    3. Tolerance fallback ladder – if the first attempt fails, two
       progressively looser tolerance pairs are tried before giving up.
    """
    # --- 1. Clamp ICs: populations cannot be negative --------------------
    ic = np.maximum(np.asarray(ic, dtype=float), 0.0)

    t_span = [0.0, float(duration)]

    # Tolerance ladder: tight → medium → loose
    tol_ladder = [
        dict(rtol=1e-6,  atol=1e-8,  max_step=1.0),
        dict(rtol=1e-4,  atol=1e-6,  max_step=0.5),
        dict(rtol=1e-3,  atol=1e-5,  max_step=0.1),
    ]

    sol = None
    last_msg = ''
    for tols in tol_ladder:
        sol = solve_ivp(
            rhs,
            t_span,
            ic,
            args=(params, controls, d_R, d_c, d_I),
            method='BDF',
            dense_output=True,
            **tols
        )
        if sol.success:
            break
        last_msg = sol.message

    if not sol.success:
        raise RuntimeError(f"Phase solver failed: {last_msg}")

    t_eval = np.linspace(0.0, sol.t[-1], n_points)
    y_eval = sol.sol(t_eval).T          # shape (n_points, 13)
    # Clamp output too – dense interpolation can produce tiny negatives
    y_eval = np.maximum(y_eval, 0.0)
    return t_eval, y_eval

# ============================================================================
#  PROTOCOL CONSTANTS  (paper §7.4)
# ============================================================================
DT      = 200          # detection time – no treatment pre-phase (days)

# Chemotherapy
D_C     = 800          # dose  (mg m⁻²)
T_C     = 14           # duration of one cycle (days)
N_C     = 6            # number of cycles per chemo block

# Radiotherapy
D_R     = 60.0 / 28.0  # dose per fraction ≈ 2.14 Gy  (60 Gy in 28 fractions)
T_R     = 40           # total radiotherapy duration (days)

# Free-treatment gap
FT      = 15           # duration (days)

# Immunotherapy  (Protocol 2 only)
D_I     = 2.0          # dose  (mg day⁻¹)
T_I_ON  = 20           # days ON per cycle
T_I_OFF = 1            # rest days between cycles
N_I     = 10           # number of cycles

# Total durations (convenient reference)
_P1_DUR = DT + N_C * T_C + T_R + FT + N_C * T_C                 # 423 days
_P2_DUR = _P1_DUR + N_I * T_I_ON + (N_I - 1) * T_I_OFF          # 632 days


# ============================================================================
#  PHASE-SCHEDULE BUILDERS
# ============================================================================
def _phase(duration, d_R=0.0, d_c=0.0, d_I=0.0, color='black', label='No treatment'):
    return dict(duration=duration, d_R=d_R, d_c=d_c, d_I=d_I, color=color, label=label)


def build_notreatment_phases(total_duration=_P2_DUR):
    """Single free-running phase (no treatment) up to total_duration."""
    return [_phase(total_duration)]


def build_protocol1_phases(d_c=D_C, d_R=D_R):
    """
    DT_200 → (Ch^d_c_{T_C})_{N_C} → (R^d_R_{T_R}) → FT_{FT} → (Ch^d_c_{T_C})_{N_C}

    Parameters
    ----------
    d_c : float   chemotherapy drug concentration (mg m⁻²)
    d_R : float   radiotherapy dose per fraction (Gy)
    """
    phases = [_phase(DT)]                                          # pre-treatment

    for _ in range(N_C):                                           # 1st chemo block
        phases.append(_phase(T_C, d_c=d_c, color='green', label='Chemotherapy'))

    phases.append(_phase(T_R, d_R=d_R, color='red', label='Radiotherapy'))
    phases.append(_phase(FT))                                      # free period

    for _ in range(N_C):                                           # 2nd chemo block
        phases.append(_phase(T_C, d_c=d_c, color='green', label='Chemotherapy'))

    return phases


def build_protocol2_phases(d_c=D_C, d_R=D_R, d_I=D_I):
    """
    Protocol 1 followed by (I^d_I_{T_I_ON})_{N_I}

    Parameters
    ----------
    d_c : float   chemotherapy dose (mg m⁻²)
    d_R : float   radiotherapy dose per fraction (Gy)
    d_I : float   immunotherapy dose (mg day⁻¹)
    """
    phases = build_protocol1_phases(d_c=d_c, d_R=d_R)

    for i in range(N_I):                                           # immuno cycles
        phases.append(_phase(T_I_ON, d_I=d_I, color='magenta', label='Immunotherapy'))
        if i < N_I - 1:                                            # rest between cycles
            phases.append(_phase(T_I_OFF))

    return phases


# ============================================================================
#  PHASED ODE RUNNER
# ============================================================================
def run_phased_simulation(phases, params, controls, ic, pts_per_day=30):
    """
    Chain multiple ODE segments with different treatment doses.

    Each phase ends at the state that seeds the next phase. Time is accumulated
    so the returned arrays form one continuous trajectory.

    Parameters
    ----------
    phases      : list of phase dicts (from build_*_phases helpers)
    params      : ndarray  – model parameters
    controls    : ndarray  – [alpha, beta, f_c, M_c, k_S, M_Tc, M_TH1]
    ic          : ndarray  – initial conditions (length 13)
    pts_per_day : int      – output resolution (points per day)

    Returns
    -------
    t_all    : ndarray, shape (N,)
    y_all    : ndarray, shape (N, 13)   – all 13 state variables
    col_all  : list[str], length N      – colour tag for each time point
    """
    t_all, y_all, col_all = [], [], []
    t_offset  = 0.0
    current_ic = ic.copy()

    for ph in phases:
        dur   = ph['duration']
        if dur <= 0:
            continue
        n_pts = max(int(dur * pts_per_day), 10)

        # Use _solve_phase (wraps model.rhs directly) for robustness across
        # all mutation conditions, including mC=mS=0 where SR=CR=0 always.
        t_phase_raw, y_phase = _solve_phase(
            current_ic, params, controls,
            ph['d_R'], ph['d_c'], ph['d_I'],
            dur, n_pts
        )
        t_phase = t_phase_raw + t_offset

        t_all.append(t_phase)
        y_all.append(y_phase)
        col_all.extend([ph['color']] * len(t_phase))

        t_offset   += dur
        current_ic  = y_phase[-1].copy()

    return np.concatenate(t_all), np.concatenate(y_all, axis=0), col_all


# ============================================================================
#  PUBLIC WRAPPER  – call this from external sweep scripts
# ============================================================================
STATE_NAMES = ['S','SR','C','CR','M1','M2','TH1','TH2','TC','Treg','IL10','IFNgamma','IL2']

def run_protocol_with_doses(
        d_I      = D_I,
        d_R_dose = D_R,
        d_c_dose = D_C,
        mC       = None,
        mS       = None,
        protocol = 'protocol2',
        params   = None,
        controls = None,
        ic       = None,
        pts_per_day = 30
):
    """
    Run a treatment protocol with arbitrary dose parameters.
    Designed for wrapper scripts that sweep over (d_I, d_R_dose, d_c_dose).

    Parameters
    ----------
    d_I       : float   immunotherapy dose  (mg day⁻¹)
    d_R_dose  : float   radiotherapy dose per fraction  (Gy)
    d_c_dose  : float   chemotherapy concentration  (mg m⁻²)
    mC        : float | None   cancer-cell mutation rate  (None → keep params default)
    mS        : float | None   stem-cell mutation rate    (None → keep params default)
    protocol  : str     'none' | 'protocol1' | 'protocol2'
    params    : ndarray | None  – full parameter vector  (None → defaults)
    controls  : ndarray | None  – control vector         (None → defaults)
    ic        : ndarray | None  – initial conditions     (None → defaults)
    pts_per_day : int   output resolution

    Returns
    -------
    dict with keys
        't'              : ndarray  – time axis
        'y'              : ndarray shape (N, 13) – all state variables
        'state_names'    : list[str]
        'colors'         : list[str] – colour per time point
        'final_state'    : ndarray length 13
        'TH1_TH2_ratio'  : float  at end of simulation
        'fold_change'    : float  tumor mass at DT / tumor mass at end
        'phases'         : list of phase dicts used
    """
    if params   is None: params   = get_params()
    if controls is None: controls = get_default_controls()
    if ic       is None: ic       = get_default_ic()

    params = params.copy()
    if mC is not None: params[PARAM_NAMES['mC']] = mC
    if mS is not None: params[PARAM_NAMES['mS']] = mS

    if   protocol == 'none':
        total = _P2_DUR
        phases = build_notreatment_phases(total)
    elif protocol == 'protocol1':
        phases = build_protocol1_phases(d_c=d_c_dose, d_R=d_R_dose)
    elif protocol == 'protocol2':
        phases = build_protocol2_phases(d_c=d_c_dose, d_R=d_R_dose, d_I=d_I)
    else:
        raise ValueError(f"Unknown protocol: {protocol!r}. Choose 'none','protocol1','protocol2'.")

    t_arr, y_arr, colors = run_phased_simulation(phases, params, controls, ic, pts_per_day)

    # ── efficacy metrics ─────────────────────────────────────────────────────
    dt_idx = np.searchsorted(t_arr, DT)
    tumor_detection = float(np.sum(y_arr[dt_idx, :4]))          # S+SR+C+CR at DT
    tumor_end       = float(np.sum(y_arr[-1,    :4]))
    fold_change     = tumor_detection / max(tumor_end, 1e-30)   # >1 means reduction

    TH1_TH2_ratio = float(y_arr[-1, 6]) / max(float(y_arr[-1, 7]), 1e-30)

    return dict(
        t             = t_arr,
        y             = y_arr,
        state_names   = STATE_NAMES,
        colors        = colors,
        final_state   = y_arr[-1].copy(),
        TH1_TH2_ratio = TH1_TH2_ratio,
        fold_change   = fold_change,
        phases        = phases
    )


# ============================================================================
#  INTERNAL – run all 6 figure cases
# ============================================================================
def _run_all_cases(pts_per_day=30):
    """Run the 6 cases needed for the extended Figure 4."""
    params_base = get_params()
    controls    = get_default_controls()
    ic          = get_default_ic()

    # mc=ms=0
    params_no_mut = params_base.copy()
    params_no_mut[PARAM_NAMES['mC']] = 0.0
    params_no_mut[PARAM_NAMES['mS']] = 0.0

    # mc,ms>0  (use model defaults: mC=0.01, mS=4e-7)
    params_mut = params_base.copy()

    cases = {}

    for tag, params in [('nomut', params_no_mut), ('mut', params_mut)]:
        mut_label = r'$m_C=m_S=0$' if tag == 'nomut' else r'$m_C,m_S>0$'
        print(f"\n── {mut_label} ──")

        for proto_key, proto_label, phase_fn in [
            ('notx', 'No Treatment',  lambda: build_notreatment_phases(_P2_DUR)),
            ('p1',   'Protocol 1',    build_protocol1_phases),
            ('p2',   'Protocol 2',    build_protocol2_phases),
        ]:
            print(f"   {proto_label} ...", end=' ', flush=True)
            phases = phase_fn()
            t, y, c = run_phased_simulation(phases, params, controls, ic, pts_per_day)
            cases[f'{tag}_{proto_key}'] = dict(t=t, y=y, colors=c,
                                               mut_label=mut_label,
                                               proto_label=proto_label)
            print(f"done  (t_end={t[-1]:.0f} d)")

    return cases


# ============================================================================
#  PLOTTING HELPERS
# ============================================================================
def _plot_colored_line(ax, t, y, colors, lw=1.2):
    """
    Draw y vs t as a single line whose colour changes segment-by-segment
    according to treatment phase, with no gaps at transitions.
    """
    if len(t) == 0:
        return

    # Collect contiguous colour runs
    runs = []
    start = 0
    cur_col = colors[0]
    for i in range(1, len(colors)):
        if colors[i] != cur_col:
            runs.append((start, i, cur_col))
            start   = i
            cur_col = colors[i]
    runs.append((start, len(colors), cur_col))

    for i, (s, e, col) in enumerate(runs):
        # Extend each segment by one point into the next to avoid white gaps
        end = min(e + 1, len(t))
        ax.plot(t[s:end], y[s:end], color=col, linewidth=lw)


def _make_figure(cases, save_path):
    """Build and save the 4-row × 6-column figure."""

    col_order = [
        ('nomut_notx', r'$m_C=m_S=0$'   + '\nNo Treatment'),
        ('nomut_p1',   r'$m_C=m_S=0$'   + '\nProtocol 1'),
        ('nomut_p2',   r'$m_C=m_S=0$'   + '\nProtocol 2'),
        ('mut_notx',   r'$m_C,m_S>0$'   + '\nNo Treatment'),
        ('mut_p1',     r'$m_C,m_S>0$'   + '\nProtocol 1'),
        ('mut_p2',     r'$m_C,m_S>0$'   + '\nProtocol 2'),
    ]

    row_info = [
        (0, 'Stem Cells $(S)$'),
        (1, 'Stem Resistant $(S_R)$'),
        (2, 'Cancer Cells $(C)$'),
        (3, 'Cancer Resistant $(C_R)$'),
    ]

    n_rows, n_cols = 4, 6
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(24, 13),
                             constrained_layout=False)
    fig.subplots_adjust(left=0.06, right=0.99, top=0.91, bottom=0.10,
                        hspace=0.45, wspace=0.38)

    fig.suptitle(
        'Tumor Dynamics Under Treatment Protocols  '
        r'($m_C=m_S=0$ vs $m_C,m_S>0$  ·  No Tx / Protocol 1 / Protocol 2)',
        fontsize=12, fontweight='bold', y=0.97
    )

    for col_idx, (case_key, col_title) in enumerate(col_order):
        data = cases[case_key]
        t      = data['t']
        y      = data['y']
        colors = data['colors']

        for row_idx, (sv_idx, row_label) in enumerate(row_info):
            ax = axes[row_idx, col_idx]
            _plot_colored_line(ax, t, y[:, sv_idx], colors)

            # Titles & labels
            if row_idx == 0:
                ax.set_title(col_title, fontsize=8.5, pad=5)
            if col_idx == 0:
                ax.set_ylabel(row_label, fontsize=8)
            if row_idx == n_rows - 1:
                ax.set_xlabel('Time (days)', fontsize=8)

            ax.tick_params(labelsize=7)
            ax.set_xlim(left=0, right=t[-1])
            # ax.set_xlim(0, 800) # SETTING THIS TO BE A FIXED THING FOR CONVENIENCE
            ax.ticklabel_format(style='sci', axis='y', scilimits=(-2, 4))
            ax.yaxis.get_offset_text().set_fontsize(6)

            # Shade background lightly to distinguish mutation columns
            if col_idx >= 3:
                ax.set_facecolor('#f7f7ff')

    # Legend
    legend_handles = [
        mpatches.Patch(color='black',   label='No treatment / free period'),
        mpatches.Patch(color='green',   label='Chemotherapy'),
        mpatches.Patch(color='red',     label='Radiotherapy'),
        mpatches.Patch(color='magenta', label='Immunotherapy'),
    ]
    fig.legend(handles=legend_handles, loc='lower center', ncol=4,
               fontsize=10, bbox_to_anchor=(0.5, 0.01),
               frameon=True, edgecolor='grey')

    # Dividing line between mc=ms=0 and mc,ms>0 columns
    line_x = (axes[0, 2].get_position().x1 + axes[0, 3].get_position().x0) / 2
    fig.add_artist(plt.Line2D([line_x, line_x], [0.08, 0.95],
                              transform=fig.transFigure,
                              color='grey', linewidth=1.2, linestyle='--'))

    # Column-group labels
    fig.text(0.26, 0.94, r'$m_C = m_S = 0$  (no resistant cells)',
             ha='center', fontsize=9, style='italic', color='#333333')
    fig.text(0.74, 0.94, r'$m_C, m_S > 0$  (resistant cells present)',
             ha='center', fontsize=9, style='italic', color='#333333')

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\nFigure saved → {save_path}")


# ============================================================================
#  ENTRY POINT
# ============================================================================
if __name__ == '__main__':
    print("=" * 60)
    print("protocols.py  –  Ganguli & Sarkar 2018, Fig. 4 (extended)")
    print("=" * 60)
    print(f"Protocol 1 total duration : {_P1_DUR} days")
    print(f"Protocol 2 total duration : {_P2_DUR} days")

    cases = _run_all_cases(pts_per_day=25)

    out = './reproduced plots/figure_4.png'
    _make_figure(cases, out)
    print("All done.")

    if os.path.exists("./__pycache__"):
        shutil.rmtree("./__pycache__")
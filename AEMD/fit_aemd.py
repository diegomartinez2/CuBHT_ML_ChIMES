#!/usr/bin/env python3
"""
fit_aemd.py
Fit AEMD <ΔT(t)> decay to extract thermal diffusivity α and thermal conductivity κ,
and save a plot overlaying data + fitted curve.

Inputs (from the LAMMPS workflow):
- geom.info          : last non-comment line must be "N_atoms  Volume_A^3  Lz_A"
- aemd_deltaT.dat    : columns "step  time(ps)  Th  Tc  dT"  (time column optional)

Model (step-like initial condition, PBC):
<ΔT(t)> = (8 ΔT0 / π^2) * Σ_{n=0..∞} 1/(2n+1)^2 * exp(-((2n+1)π/Lz)^2 * α * t)

Notes:
- α returned in Å^2/ps, converted to m^2/s via 1 Å^2/ps = 1e-8 m^2/s
- κ computed using classical volumetric heat capacity Cvol ≈ 3 n kB
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np


def read_geom(fname: str = "geom.info") -> tuple[int, float, float]:
    """
    Expects:
      # natoms vol_A3 lz_A
      N  V_A3  Lz_A
    Uses the last non-comment, non-empty line.
    """
    p = Path(fname)
    if not p.exists():
        raise FileNotFoundError(f"Geometry file not found: {fname}")

    lines = []
    for ln in p.read_text().splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        lines.append(s)

    if not lines:
        raise RuntimeError(f"No usable data found in {fname}")

    parts = lines[-1].split()
    if len(parts) < 3:
        raise RuntimeError(f"Expected at least 3 columns in {fname}, got: {lines[-1]}")

    N = int(float(parts[0]))
    V_A3 = float(parts[1])
    Lz_A = float(parts[2])
    if N <= 0 or V_A3 <= 0.0 or Lz_A <= 0.0:
        raise ValueError(f"Non-physical geom values: N={N}, V_A3={V_A3}, Lz_A={Lz_A}")

    return N, V_A3, Lz_A


def read_aemd_deltaT(fname: str = "aemd_deltaT.dat") -> tuple[np.ndarray, np.ndarray | None, np.ndarray, np.ndarray, np.ndarray]:
    """
    Supports two common formats:

    (A) With time column (recommended):
      step  time(ps)  Th  Tc  dT

    (B) Without time column:
      step  Th  Tc  dT

    Returns: step, time_or_None, Th, Tc, dT
    """
    p = Path(fname)
    if not p.exists():
        raise FileNotFoundError(f"ΔT file not found: {fname}")

    data = np.genfromtxt(fname, comments="#")
    if data.size == 0:
        raise RuntimeError(f"No data read from {fname}")
    if data.ndim == 1:
        data = data[None, :]

    ncol = data.shape[1]
    if ncol == 5:
        step = data[:, 0]
        time_ps = data[:, 1]
        Th = data[:, 2]
        Tc = data[:, 3]
        dT = data[:, 4]
        return step, time_ps, Th, Tc, dT
    elif ncol >= 4:
        step = data[:, 0]
        Th = data[:, 1]
        Tc = data[:, 2]
        dT = data[:, 3]
        return step, None, Th, Tc, dT
    else:
        raise RuntimeError(f"Expected 4 or 5 columns in {fname}, got shape {data.shape}")


def model_deltaT(
    t_ps: np.ndarray,
    alpha_A2_per_ps: float,
    dT0: float,
    Lz_A: float,
    nterms: int = 20,
) -> np.ndarray:
    if nterms < 1:
        raise ValueError("nterms must be >= 1")
    odd = (2 * np.arange(nterms) + 1).astype(float)  # 1, 3, 5, ...
    coef = 1.0 / (odd**2)
    lam2 = (odd * np.pi / Lz_A) ** 2  # 1/Å^2
    expo = np.exp(-lam2[:, None] * alpha_A2_per_ps * t_ps[None, :])
    series = np.sum(coef[:, None] * expo, axis=0)
    return (8.0 * dT0 / (np.pi**2)) * series


def fit_alpha(
    t_ps: np.ndarray,
    dT: np.ndarray,
    Lz_A: float,
    nterms: int = 20,
) -> tuple[float, float, float]:
    """
    Lightweight fit (no SciPy):
    - grid search in log-space for alpha
    - for each alpha, solve best dT0 by linear least squares
    Returns: (alpha_A2_per_ps, dT0_fit, mse)
    """
    t_ps = np.asarray(t_ps, dtype=float)
    dT = np.asarray(dT, dtype=float)

    if len(t_ps) < 10:
        raise RuntimeError("Not enough points to fit (need at least ~10).")

    dT_abs0 = float(np.abs(dT[0]))
    if not np.isfinite(dT_abs0) or dT_abs0 == 0.0:
        raise RuntimeError("ΔT(0) is zero or non-finite; cannot fit.")

    # crude guess for tau using 1/e crossing
    target = dT_abs0 / np.e
    idx = np.where(np.abs(dT) < target)[0]
    tau_ps = float(t_ps[idx[0]]) if len(idx) > 0 else float(t_ps[-1])
    tau_ps = max(tau_ps, 1e-6)

    alpha_guess = (Lz_A**2) / (np.pi**2 * tau_ps)  # Å^2/ps

    alphas = alpha_guess * 10 ** np.linspace(-2.0, 2.0, 121)

    best = None  # (mse, alpha, dT0)
    for a in alphas:
        base = model_deltaT(t_ps, a, 1.0, Lz_A, nterms=nterms)
        denom = float(np.dot(base, base))
        if denom <= 1e-30:
            continue
        dT0 = float(np.dot(base, dT) / denom)
        pred = dT0 * base
        mse = float(np.mean((pred - dT) ** 2))
        if (best is None) or (mse < best[0]):
            best = (mse, a, dT0)

    if best is None:
        raise RuntimeError("Fit failed during coarse search.")

    # refine locally around best alpha
    _, a0, _ = best
    for _ in range(6):
        alphas = a0 * 10 ** np.linspace(-0.25, 0.25, 41)
        best_local = None
        for a in alphas:
            base = model_deltaT(t_ps, a, 1.0, Lz_A, nterms=nterms)
            denom = float(np.dot(base, base))
            if denom <= 1e-30:
                continue
            dT0 = float(np.dot(base, dT) / denom)
            pred = dT0 * base
            mse = float(np.mean((pred - dT) ** 2))
            if (best_local is None) or (mse < best_local[0]):
                best_local = (mse, a, dT0)
        if best_local is None:
            break
        best = best_local
        a0 = best[1]

    mse, alpha_fit, dT0_fit = best
    return float(alpha_fit), float(dT0_fit), float(mse)


def save_plot(
    t_ps_all: np.ndarray,
    dT_all: np.ndarray,
    t_ps_fit: np.ndarray,
    dT_fit: np.ndarray,
    alpha_A2ps: float,
    dT0_fit: float,
    Lz_A: float,
    nterms: int,
    plotfile: str,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        raise RuntimeError(
            "matplotlib is required to save the plot. Install it (e.g. pip install matplotlib)."
        ) from e

    # Smooth curve over the fit window
    tmin = float(np.min(t_ps_fit))
    tmax = float(np.max(t_ps_fit))
    t_smooth = np.linspace(tmin, tmax, 600)
    dT_smooth = model_deltaT(t_smooth, alpha_A2ps, dT0_fit, Lz_A, nterms=nterms)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # all data (faint)
    ax.plot(t_ps_all, dT_all, ".", markersize=3, alpha=0.35, label="ΔT data (all)")

    # fit window data (emphasize)
    ax.plot(t_ps_fit, dT_fit, ".", markersize=4, alpha=0.9, label="ΔT data (fit window)")

    # fitted curve
    ax.plot(t_smooth, dT_smooth, "-", linewidth=2.0, label="Fit (series model)")

    ax.set_xlabel("time (ps)")
    ax.set_ylabel("ΔT (K)")
    ax.set_title("AEMD: ΔT(t) data and fit")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(plotfile, dpi=200)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Fit AEMD <ΔT(t)> to extract thermal diffusivity and thermal conductivity, and save a ΔT(t) plot."
    )
    ap.add_argument("--geom", default="geom.info", help="Geometry file written by LAMMPS (geom.info).")
    ap.add_argument("--data", default="aemd_deltaT.dat", help="ΔT file (columns: step [time] Th Tc dT).")
    ap.add_argument("--dt_fs", type=float, default=1.0, help="LAMMPS timestep in fs (used only if file has no time column).")
    ap.add_argument("--nterms", type=int, default=20, help="Number of terms in the series (e.g., 10, 20, 40, 80).")
    ap.add_argument("--tmin_ps", type=float, default=0.0, help="Ignore data earlier than this time (ps).")
    ap.add_argument("--tmax_ps", type=float, default=None, help="Ignore data later than this time (ps).")
    ap.add_argument("--plotfile", default="aemd_fit.png", help="Output image filename (PNG recommended).")
    args = ap.parse_args()

    if args.dt_fs <= 0:
        raise ValueError("--dt_fs must be > 0")
    if args.nterms < 1:
        raise ValueError("--nterms must be >= 1")

    N, V_A3, Lz_A = read_geom(args.geom)
    step, time_ps, Th, Tc, dT = read_aemd_deltaT(args.data)

    # build time axis
    if time_ps is None:
        dt_ps = args.dt_fs * 1.0e-3
        t_ps = step * dt_ps
    else:
        t_ps = time_ps

    # filter
    mask = np.isfinite(dT) & np.isfinite(t_ps)
    mask &= (t_ps >= args.tmin_ps)
    if args.tmax_ps is not None:
        mask &= (t_ps <= args.tmax_ps)

    t_ps_fit = t_ps[mask]
    dT_fit = dT[mask]

    if len(t_ps_fit) < 10:
        raise RuntimeError("Not enough points after applying tmin/tmax filtering.")

    # Ensure increasing time for the fit (not strictly necessary, but safer)
    order = np.argsort(t_ps_fit)
    t_ps_fit = t_ps_fit[order]
    dT_fit = dT_fit[order]

    alpha_A2ps, dT0_fit, mse = fit_alpha(t_ps_fit, dT_fit, Lz_A, nterms=args.nterms)

    # Convert alpha to SI: (Å^2/ps) -> (m^2/s)
    alpha_SI = alpha_A2ps * 1.0e-8

    # Classical volumetric heat capacity: Cvol = 3 * n * kB
    kB = 1.380649e-23  # J/K
    V_m3 = V_A3 * 1.0e-30
    n = N / V_m3  # 1/m^3
    Cvol = 3.0 * n * kB  # J/(m^3 K)

    kappa_WmK = alpha_SI * Cvol

    # Save plot overlaying data + fit
    save_plot(
        t_ps_all=t_ps,
        dT_all=dT,
        t_ps_fit=t_ps_fit,
        dT_fit=dT_fit,
        alpha_A2ps=alpha_A2ps,
        dT0_fit=dT0_fit,
        Lz_A=Lz_A,
        nterms=args.nterms,
        plotfile=args.plotfile,
    )

    print("=== AEMD fit results ===")
    print(f"N atoms             : {N}")
    print(f"Volume (Å^3)        : {V_A3:.6e}")
    print(f"Lz (Å)              : {Lz_A:.6f}")
    print(f"nterms              : {args.nterms}")
    print(f"Fit window (ps)     : [{float(np.min(t_ps_fit)):.6g}, {float(np.max(t_ps_fit)):.6g}]  (N={len(t_ps_fit)})")
    print(f"ΔT0 (fit, K)        : {dT0_fit:.6f}")
    print(f"alpha (Å^2/ps)      : {alpha_A2ps:.6e}")
    print(f"alpha (m^2/s)       : {alpha_SI:.6e}")
    print(f"Cvol (J/m^3/K)      : {Cvol:.6e}  (classical 3 n kB)")
    print(f"kappa (W/m/K)       : {kappa_WmK:.6f}")
    print(f"MSE (K^2)           : {mse:.6e}")
    print(f"Saved plot          : {args.plotfile}")


if __name__ == "__main__":
    main()

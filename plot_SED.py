#!/usr/bin/env python3
"""Plot SEDs from merged photometry files without pyhdust.

Expected input format per line in ``{star}/{star}_photometry_manual.dat``:
    wavelength_um  flux_jy  flux_err_jy  label  [optional extra tokens...]

For radio measurements, prefer direct survey labels in column 4
(e.g., ``laboca``, ``vla``, ``jvla``). Legacy ``radio <dataset>``
and ``radio_upper <dataset>`` rows are also supported.

Trailing tokens after ``label`` are ignored except for legacy
``radio`` rows, where the next token is treated as the dataset label.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.container import ErrorbarContainer
from matplotlib.legend_handler import HandlerErrorbar

try:
    from auxiliary_vieira_et_al_2015_eq20 import calc_n
except Exception:
    calc_n = None


#
# USER CONFIG (edit these in-script if you prefer not passing CLI args)
#
# starlist = ('HD44637', '10CMa', '20Vul', '25Vul', '25Peg', '120Tau', 'epsPsA', 'HR2249',\
#             'omeOri', 'upsCyg', 'zetCrv')
starlist = ('betCMi',)
photometry_file = "{star}/{star}_photometry.dat"
plot_file = "{star}/{star}_SED.png"
iue_avg_file = "{star}/{star}_iue_average_jy.dat"
hpol_avg_file = "{star}/{star}_hpol_average_jy.dat"
ebv = 0.1
rv = 3.1

ir_fit_wavelength_min_um = 8.0
ir_fit_wavelength_max_um = 100.0

radio_fit_wavelength_min_cm = 0.08
radio_fit_wavelength_max_cm = 15

# Dotted-line extrapolation range for IR fit [um].
# Set to None to auto-use fitted-data limits.
ir_fit_plot_wavelength_min_um = ir_fit_wavelength_min_um
ir_fit_plot_wavelength_max_um = ir_fit_wavelength_max_um

# Dotted-line extrapolation range for Radio fit [um].
# Set to None to auto-use fitted-data limits.
radio_fit_plot_wavelength_min_um = 500.0
radio_fit_plot_wavelength_max_um = radio_fit_wavelength_max_cm * 1.0e4

label_order = ["simbad", "dougherty+1991", "iras", "wise", "akari", "laboca", "wendker+2000", "waters+1991", "karma", "vla", "jvla"]

# label -> (legend_label, marker, color)
styles = {
    "simbad": ("Simbad", "o", "tab:blue"),
    "dougherty+1991": ("Dougherty+1991", "^", "tab:green"),
    "iras": ("IRAS", "x", "tab:red"),
    "wise": ("WISE", "s", "tab:orange"),
    "akari": ("AKARI", "v", "tab:purple"),
    "laboca": ("APEX/LABOCA", "D", "tab:brown"),
    "wendker+2000": ("30-m Pico Veleta", "P", "tab:pink"),
    "waters+1991": ("JCMT", "X", "tab:cyan"),
    "karma": ("CARMA", "*", "tab:olive"),
    "vla": ("VLA", "<", "tab:gray"),
    "jvla": ("JVLA", ">", "k"),
}


def ccm89_ab(x_inv_micron: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return CCM89 coefficients a(x), b(x) for x=1/lambda[um].

    For far-IR (x < 0.3; lambda > 3.33 um), extinction is set to 0 here.
    """
    x = np.asarray(x_inv_micron, dtype=float)
    a = np.zeros_like(x)
    b = np.zeros_like(x)

    # IR: 0.3 <= x < 1.1
    m_ir = (x >= 0.3) & (x < 1.1)
    if np.any(m_ir):
        xx = x[m_ir]
        a[m_ir] = 0.574 * xx ** 1.61
        b[m_ir] = -0.527 * xx ** 1.61

    # Optical/NIR: 1.1 <= x <= 3.3
    m_opt = (x >= 1.1) & (x <= 3.3)
    if np.any(m_opt):
        y = x[m_opt] - 1.82
        a[m_opt] = (
            1
            + 0.17699 * y
            - 0.50447 * y**2
            - 0.02427 * y**3
            + 0.72085 * y**4
            + 0.01979 * y**5
            - 0.77530 * y**6
            + 0.32999 * y**7
        )
        b[m_opt] = (
            1.41338 * y
            + 2.28305 * y**2
            + 1.07233 * y**3
            - 5.38434 * y**4
            - 0.62251 * y**5
            + 5.30260 * y**6
            - 2.09002 * y**7
        )

    # UV: 3.3 < x <= 8.0
    m_uv = (x > 3.3) & (x <= 8.0)
    if np.any(m_uv):
        xx = x[m_uv]
        fa = np.zeros_like(xx)
        fb = np.zeros_like(xx)
        m_fuv = xx >= 5.9
        if np.any(m_fuv):
            yy = xx[m_fuv] - 5.9
            fa[m_fuv] = -0.04473 * yy**2 - 0.009779 * yy**3
            fb[m_fuv] = 0.2130 * yy**2 + 0.1207 * yy**3

        a[m_uv] = 1.752 - 0.316 * xx - 0.104 / ((xx - 4.67) ** 2 + 0.341) + fa
        b[m_uv] = -3.090 + 1.825 * xx + 1.206 / ((xx - 4.62) ** 2 + 0.263) + fb

    # For x > 8.0 (far-UV), clip to x=8.0 approximation.
    m_hi = x > 8.0
    if np.any(m_hi):
        aa, bb = ccm89_ab(np.full(np.count_nonzero(m_hi), 8.0))
        a[m_hi], b[m_hi] = aa, bb

    return a, b


def deredden_fnu(wavelength_um: np.ndarray, flux_fnu: np.ndarray, ebv: float, rv: float) -> np.ndarray:
    """Deredden flux density F_nu using CCM89 extinction law."""
    if ebv <= 0:
        return flux_fnu

    wl = np.asarray(wavelength_um, dtype=float)
    flux = np.asarray(flux_fnu, dtype=float)
    valid = wl > 0
    out = flux.copy()

    x = np.zeros_like(wl)
    x[valid] = 1.0 / wl[valid]
    a, b = ccm89_ab(x)
    a_lambda = ebv * (a * rv + b)
    out[valid] = flux[valid] * (10.0 ** (0.4 * a_lambda[valid]))
    return out


def load_photometry(path: Path) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Load detections and upper limits from a photometry file.

    Returns
    -------
    detections, upper_limits : dict[label] -> ndarray shape (N, 3)
        Columns are wavelength_um, flux_jy, err_jy.
    """
    detections_raw: Dict[str, List[Tuple[float, float, float]]] = {}
    upper_raw: Dict[str, List[Tuple[float, float, float]]] = {}

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 4:
                continue

            try:
                wl_um = float(parts[0])
                flux_jy = float(parts[1])
                err_jy = float(parts[2])
            except ValueError:
                continue

            label = parts[3].strip().lower()
            if wl_um <= 0 or flux_jy <= 0:
                continue

            target = detections_raw
            if label.endswith("_upper"):
                label = label[: -len("_upper")]
                target = upper_raw

            # Support compact radio notation:
            # wl flux err radio <dataset> [optional notes...]
            if label == "radio" and len(parts) >= 5:
                radio_label = parts[4].strip().lower().strip("()[]{}.,;:")
                if radio_label:
                    label = radio_label

            target.setdefault(label, []).append((wl_um, flux_jy, err_jy))

    def finalize(raw: Dict[str, List[Tuple[float, float, float]]]) -> Dict[str, np.ndarray]:
        out: Dict[str, np.ndarray] = {}
        for label, rows in raw.items():
            arr = np.array(rows, dtype=float)
            order = np.argsort(arr[:, 0])
            out[label] = arr[order]
        return out

    return finalize(detections_raw), finalize(upper_raw)


def load_spectrum_jy(path: Path) -> np.ndarray:
    """Load a two-column spectrum file: wavelength_um, flux_jy."""
    rows: List[Tuple[float, float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                wl_um = float(parts[0])
                flux_jy = float(parts[1])
            except ValueError:
                continue
            if wl_um > 0 and flux_jy > 0:
                rows.append((wl_um, flux_jy))

    if not rows:
        return np.empty((0, 2), dtype=float)

    arr = np.array(rows, dtype=float)
    order = np.argsort(arr[:, 0])
    return arr[order]


def ordered_labels(*groups: Dict[str, np.ndarray]) -> List[str]:
    labels = set()
    for group in groups:
        labels.update(group.keys())

    ordered = [lbl for lbl in label_order if lbl in labels]
    ordered += sorted(lbl for lbl in labels if lbl not in label_order)
    return ordered


def style_for_label(label: str, idx: int) -> Tuple[str, str, str]:
    if label in styles:
        return styles[label]

    cmap = plt.get_cmap("tab20")
    color = cmap(idx % cmap.N)
    marker_cycle = ["o", "^", "s", "v", "D", "P", "X", "<", ">"]
    marker = marker_cycle[idx % len(marker_cycle)]
    return label, marker, color


def fit_power_law_loglog(
    wavelength_um: np.ndarray,
    flux_fnu: np.ndarray,
    err_fnu: np.ndarray,
) -> Tuple[float, float]:
    """Fit y = A * x^k in log-log space, weighted when possible.

    Returns
    -------
    amplitude, slope
    """
    x = np.asarray(wavelength_um, dtype=float)
    y = np.asarray(flux_fnu, dtype=float)
    e = np.asarray(err_fnu, dtype=float)

    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (x > 0) & (y > 0)
    x, y, e = x[valid], y[valid], e[valid]
    if x.size < 2:
        raise ValueError("Need at least two valid points for power-law fit.")

    lx = np.log10(x)
    ly = np.log10(y)

    sigma_ly = np.full_like(ly, np.nan)
    m = e > 0
    sigma_ly[m] = e[m] / (y[m] * np.log(10.0))

    has_weights = np.isfinite(sigma_ly) & (sigma_ly > 0)
    if np.count_nonzero(has_weights) >= 2:
        lx_fit = lx[has_weights]
        ly_fit = ly[has_weights]
        w = 1.0 / sigma_ly[has_weights] ** 2
    else:
        lx_fit = lx
        ly_fit = ly
        w = np.ones_like(lx_fit)

    A = np.column_stack([np.ones_like(lx_fit), lx_fit])
    Aw = A * np.sqrt(w)[:, None]
    bw = ly_fit * np.sqrt(w)
    intercept, slope = np.linalg.lstsq(Aw, bw, rcond=None)[0]

    amplitude = 10.0 ** intercept
    return amplitude, slope


def collect_ir_points(detections: Dict[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Collect points for IR slope fit, approximating prior behavior.

    Uses:
    - last two Dougherty points (if available),
    - all WISE points,
    - all IRAS points,
    - all AKARI points.
    """
    chunks = []

    if "dougherty+1991" in detections:
        d = detections["dougherty+1991"]
        chunks.append(d[-2:] if len(d) >= 2 else d)
    if "wise" in detections:
        chunks.append(detections["wise"])
    if "iras" in detections:
        chunks.append(detections["iras"])
    if "akari" in detections:
        chunks.append(detections["akari"])

    if not chunks:
        return np.array([]), np.array([]), np.array([])

    arr = np.vstack(chunks)
    order = np.argsort(arr[:, 0])
    arr = arr[order]

    wav = arr[:, 0]
    flux_fnu = arr[:, 1] * 1e-23
    err_fnu = arr[:, 2] * 1e-23
    return wav, flux_fnu, err_fnu


def collect_radio_points(detections: Dict[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Collect candidate points for radio slope fit from all detections."""
    chunks = [arr for arr in detections.values() if len(arr) > 0]
    if not chunks:
        return np.array([]), np.array([]), np.array([])

    arr = np.vstack(chunks)
    order = np.argsort(arr[:, 0])
    arr = arr[order]

    wav = arr[:, 0]
    flux_fnu = arr[:, 1] * 1e-23
    err_fnu = arr[:, 2] * 1e-23
    return wav, flux_fnu, err_fnu


def plot_single_star(
    star: str,
    phot_path: Path,
    out_path: Path,
    ebv: float,
    rv: float,
    fit_ir: bool,
    fit_radio: bool,
    radio_fit_min_cm: float | None,
    radio_fit_max_cm: float | None,
) -> None:
    detections, upper_limits = load_photometry(phot_path)
    show_observed = ebv > 0

    if not detections and not upper_limits:
        print(f"{star}: no plottable data in {phot_path}")
        return

    fig, ax = plt.subplots(figsize=(8, 4), dpi=300)
    ax.set_xlabel(r"$\lambda\, [\mu\mathrm{m}]$", fontsize=13)
    ax.set_ylabel(r"$F_{\nu}\, [\mathrm{erg\, s^{-1}\, cm^{-2}\, Hz^{-1}}]$", fontsize=13)
    ax.set_xscale("log")
    ax.set_yscale("log")
    # ax.grid(False, which="both", ls=":", lw=0.6, alpha=0.6)
    ax.tick_params(axis="both", which="major", labelsize=11)

    def add_fit_interval_lines(min_um: float | None, max_um: float | None, color: str) -> None:
        """Draw faint vertical lines that mark fitting-interval boundaries."""
        line_kwargs = dict(color=color, ls=":", lw=0.8, alpha=0.22, zorder=0)
        if min_um is not None and np.isfinite(min_um) and min_um > 0:
            ax.axvline(min_um, **line_kwargs)
        if max_um is not None and np.isfinite(max_um) and max_um > 0:
            if min_um is None or not np.isclose(min_um, max_um):
                ax.axvline(max_um, **line_kwargs)

    add_fit_interval_lines(ir_fit_wavelength_min_um, ir_fit_wavelength_max_um, "0.35")
    radio_min_um = None if radio_fit_min_cm is None else radio_fit_min_cm * 1.0e4
    radio_max_um = None if radio_fit_max_cm is None else radio_fit_max_cm * 1.0e4
    add_fit_interval_lines(radio_min_um, radio_max_um, "0.15")

    labels = ordered_labels(detections, upper_limits)

    all_x = []
    all_y = []
    for idx, label in enumerate(labels):
        legend_label, marker, color = style_for_label(label, idx)

        if label in detections:
            arr = detections[label]
            x = arr[:, 0]
            y_obs = arr[:, 1] * 1e-23
            yerr_obs = arr[:, 2] * 1e-23

            y = deredden_fnu(x, y_obs, ebv=ebv, rv=rv)
            yerr = deredden_fnu(x, yerr_obs, ebv=ebv, rv=rv)

            if show_observed:
                ax.errorbar(
                    x,
                    y_obs,
                    yerr=yerr_obs,
                    ls="None",
                    marker=marker,
                    markersize=5,
                    markerfacecolor="none",
                    markeredgecolor=color,
                    color=color,
                    alpha=0.25,
                )

            ax.errorbar(
                x,
                y,
                yerr=yerr,
                ls="None",
                marker=marker,
                markersize=7,
                markerfacecolor="none",
                markeredgecolor=color,
                color=color,
                label=legend_label,
            )
            all_x.append(x)
            all_y.append(y)

        if label in upper_limits:
            arr = upper_limits[label]
            x = arr[:, 0]
            y_obs = arr[:, 1] * 1e-23
            yerr_obs = arr[:, 2] * 1e-23

            y = deredden_fnu(x, y_obs, ebv=ebv, rv=rv)
            yerr = deredden_fnu(x, yerr_obs, ebv=ebv, rv=rv)

            if show_observed:
                ax.errorbar(
                    x,
                    y_obs,
                    yerr=yerr_obs,
                    ls="None",
                    marker=marker,
                    markersize=5,
                    color=color,
                    uplims=True,
                    alpha=0.25,
                )

            ax.errorbar(
                x,
                y,
                yerr=yerr,
                ls="None",
                marker=marker,
                markersize=6,
                color=color,
                uplims=True,
                alpha=0.55,
            )
            all_x.append(x)
            all_y.append(y)

    spectrum_paths = [
        ("IUE", Path(iue_avg_file.format(star=star)), "black"),
        ("HPOL", Path(hpol_avg_file.format(star=star)), "tab:red"),
    ]
    for legend_label, spectrum_path, color in spectrum_paths:
        if not spectrum_path.is_file():
            continue
        spec_arr = load_spectrum_jy(spectrum_path)
        if spec_arr.size == 0:
            print(f"{star}: empty/invalid spectrum file {spectrum_path}")
            continue

        x = spec_arr[:, 0]
        y_obs = spec_arr[:, 1] * 1e-23
        y = deredden_fnu(x, y_obs, ebv=ebv, rv=rv)

        if show_observed:
            ax.plot(x, y_obs, color=color, lw=0.8, alpha=0.35, ls="--")

        ax.plot(x, y, color=color, lw=0.8, alpha=0.9, label=legend_label)
        all_x.append(x)
        all_y.append(y)

    ir_fit_text = None
    if fit_ir:
        wav_ir, flux_ir, err_ir = collect_ir_points(detections)
        mask_ir = np.isfinite(wav_ir)
        if ir_fit_wavelength_min_um is not None:
            mask_ir &= wav_ir >= ir_fit_wavelength_min_um
        if ir_fit_wavelength_max_um is not None:
            mask_ir &= wav_ir <= ir_fit_wavelength_max_um
        wav_ir = wav_ir[mask_ir]
        flux_ir = flux_ir[mask_ir]
        err_ir = err_ir[mask_ir]
        if wav_ir.size >= 2:
            flux_ir = deredden_fnu(wav_ir, flux_ir, ebv=ebv, rv=rv)
            err_ir = deredden_fnu(wav_ir, err_ir, ebv=ebv, rv=rv)
            try:
                amplitude, slope = fit_power_law_loglog(wav_ir, flux_ir, err_ir)
                xfit_min = max(0.2, np.nanmin(wav_ir) * 0.9)
                xfit_max = np.nanmax(wav_ir) * 1.1
                if ir_fit_plot_wavelength_min_um is not None:
                    xfit_min = max(0.2, ir_fit_plot_wavelength_min_um)
                if ir_fit_plot_wavelength_max_um is not None:
                    xfit_max = ir_fit_plot_wavelength_max_um
                if xfit_max <= xfit_min:
                    raise ValueError(
                        f"Invalid IR extrapolation range: min={xfit_min} um, max={xfit_max} um"
                    )
                xfit = np.geomspace(xfit_min, xfit_max, 300)
                yfit = amplitude * xfit**slope
                ax.plot(xfit, yfit, color="0.35", ls=":", lw=1.2, label="_nolegend_")

                kappa_display = abs(-slope)
                txt = f"$\kappa$ (IR) = {kappa_display:.2f}"
                if calc_n is not None:
                    try:
                        n_val = float(np.asarray(calc_n(-slope + 2.0, 15))[2])
                        txt += f", n ≃ {abs(n_val):.2f}"
                    except Exception:
                        pass
                ir_fit_text = txt
                print(f"{star}: IR fit kappa={kappa_display:.3f}")
            except Exception as exc:
                print(f"{star}: IR fit skipped ({exc})")

    radio_fit_text = None
    if fit_radio:
        wav_radio, flux_radio, err_radio = collect_radio_points(detections)
        mask_radio = np.isfinite(wav_radio)
        if radio_fit_min_cm is not None:
            mask_radio &= wav_radio >= (radio_fit_min_cm * 1.0e4)
        if radio_fit_max_cm is not None:
            mask_radio &= wav_radio <= (radio_fit_max_cm * 1.0e4)
        wav_radio = wav_radio[mask_radio]
        flux_radio = flux_radio[mask_radio]
        err_radio = err_radio[mask_radio]
        if wav_radio.size >= 2:
            flux_radio = deredden_fnu(wav_radio, flux_radio, ebv=ebv, rv=rv)
            err_radio = deredden_fnu(wav_radio, err_radio, ebv=ebv, rv=rv)
            try:
                amplitude, slope = fit_power_law_loglog(wav_radio, flux_radio, err_radio)
                xfit_min = max(0.2, np.nanmin(wav_radio) * 0.9)
                xfit_max = np.nanmax(wav_radio) * 1.1
                if radio_fit_plot_wavelength_min_um is not None:
                    xfit_min = max(0.2, radio_fit_plot_wavelength_min_um)
                if radio_fit_plot_wavelength_max_um is not None:
                    xfit_max = radio_fit_plot_wavelength_max_um
                if xfit_max <= xfit_min:
                    raise ValueError(
                        f"Invalid Radio extrapolation range: min={xfit_min} um, max={xfit_max} um"
                    )
                xfit = np.geomspace(xfit_min, xfit_max, 300)
                yfit = amplitude * xfit**slope
                ax.plot(xfit, yfit, color="0.15", ls=":", lw=1.2, label="_nolegend_")

                kappa_display = abs(-slope)
                txt = f"$\kappa$ (Radio) = {kappa_display:.2f}"
                if calc_n is not None:
                    try:
                        n_val = float(np.asarray(calc_n(-slope + 2.0, 15))[2])
                        txt += f", n ≃ {abs(n_val):.2f}"
                    except Exception:
                        pass
                radio_fit_text = txt
                print(f"{star}: radio fit kappa={kappa_display:.3f}")
            except Exception as exc:
                print(f"{star}: radio fit skipped ({exc})")

    formula_text = r"$F_{\nu} \propto \lambda^{-\kappa},\ \rho \propto \rho_0^{-n}$"
    ax.text(0.02, 0.13, formula_text, transform=ax.transAxes, fontsize=8)
    if ir_fit_text is not None:
        ax.text(0.02, 0.07, ir_fit_text, transform=ax.transAxes, fontsize=8)
    if radio_fit_text is not None:
        radio_text_y = 0.02 if ir_fit_text is not None else 0.05
        ax.text(0.02, radio_text_y, radio_fit_text, transform=ax.transAxes, fontsize=8)
    if ebv > 0:
        ax.text(
            0.02,
            0.17,
            f"E(B-V) = {ebv:.2f}",
            transform=ax.transAxes,
            fontsize=8,
            color='k',
            fontstyle="italic",
        )

    if all_x and all_y:
        x_all = np.concatenate(all_x)
        y_all = np.concatenate(all_y)
        x_min = np.nanmin(x_all[x_all > 0])
        x_max = np.nanmax(x_all[x_all > 0])
        y_min = np.nanmin(y_all[y_all > 0])
        y_max = np.nanmax(y_all[y_all > 0])

        interval_bounds = [
            ir_fit_wavelength_min_um,
            ir_fit_wavelength_max_um,
            radio_min_um,
            radio_max_um,
        ]
        interval_bounds = [
            val for val in interval_bounds if val is not None and np.isfinite(val) and val > 0
        ]
        if interval_bounds:
            x_min = min(x_min, min(interval_bounds))
            x_max = max(x_max, max(interval_bounds))

        ax.set_xlim(max(0.08, x_min / 2.0), x_max * 2.0)
        ax.set_ylim(y_min / 2.5, y_max * 2.5)

    ax.legend(
        fontsize=9,
        ncol=2,
        numpoints=1,
        handler_map={ErrorbarContainer: HandlerErrorbar(xerr_size=0.0, yerr_size=0.0)},
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, format=out_path.suffix.lstrip(".") or "png")
    up_path = out_path.parent.parent / out_path.name
    up_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(up_path, format=up_path.suffix.lstrip(".") or "png")
    plt.close(fig)

    print(f"{star}: saved {out_path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot SEDs from merged photometry files.")
    parser.add_argument(
        "stars",
        nargs="*",
        default=list(starlist),
        help="Star names (directory/file key). Example: HD44637",
    )
    parser.add_argument(
        "--phot-template",
        default=photometry_file,
        help="Input photometry template (uses {star}).",
    )
    parser.add_argument(
        "--out-template",
        default=plot_file,
        help="Output figure template (uses {star}).",
    )
    parser.add_argument("--ebv", type=float, default=ebv, help="E(B-V) for dereddening.")
    parser.add_argument("--rv", type=float, default=rv, help="R_V for CCM89 dereddening.")
    parser.add_argument(
        "--no-fit-ir",
        action="store_true",
        help="Disable IR power-law fit.",
    )
    parser.add_argument(
        "--no-fit-radio",
        action="store_true",
        help="Disable radio power-law fit.",
    )
    parser.add_argument(
        "--radio-fit-min-cm",
        type=float,
        default=radio_fit_wavelength_min_cm,
        help="Minimum wavelength [cm] for radio fit.",
    )
    parser.add_argument(
        "--radio-fit-max-cm",
        type=float,
        default=radio_fit_wavelength_max_cm,
        help="Maximum wavelength [cm] for radio fit.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    for star in args.stars:
        phot_path = Path(args.phot_template.format(star=star))
        out_path = Path(args.out_template.format(star=star))

        if not phot_path.is_file():
            print(f"{star}: missing {phot_path}")
            continue

        plot_single_star(
            star=star,
            phot_path=phot_path,
            out_path=out_path,
            ebv=args.ebv,
            rv=args.rv,
            fit_ir=not args.no_fit_ir,
            fit_radio=not args.no_fit_radio,
            radio_fit_min_cm=args.radio_fit_min_cm,
            radio_fit_max_cm=args.radio_fit_max_cm,
        )


if __name__ == "__main__":
    main()

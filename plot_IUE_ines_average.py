# IMPORTANT:
# This script expects IUE spectra downloaded from the INES database in an IUE/ subfolder.

from glob import glob
import os

from astropy.io import fits
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from PyAstronomy import pyasl


starlist = ('HD44637', '10CMa', '20Vul', '25Vul', '25Peg', '120Tau', 'epsPsA', 'HR2249',\
            'omeOri', 'upsCyg', 'zetCrv')
# starlist = ('betCMi',)
EBV = 0.00
RV = 2.9
NBINS = 170

# PyAstronomy uses np.NAN internally; keep compatibility with newer NumPy.
if not hasattr(np, 'NAN'):
    np.NAN = np.nan


def collect_iue_files(star_path):
    swp_files = sorted(glob('{}/IUE/*SWP*.FITS'.format(star_path)))
    lwp_files = sorted(glob('{}/IUE/*LW*.FITS'.format(star_path)))
    return swp_files, lwp_files


def flux_flambda_to_jy(wav_angstrom, flux_flambda):
    return flux_flambda * 3.33564095e4 * wav_angstrom**2


def setup_plot():
    fig = plt.figure(figsize=(4, 2.4), dpi=300)
    gs = gridspec.GridSpec(1, 1, height_ratios=[1], width_ratios=[1])
    gs.update(left=0.25, right=0.9, bottom=0.20, top=0.9)
    ax = plt.subplot(gs[0, 0])
    ax.set_xlabel('$\lambda$ [ $\AA$]')
    ax.set_ylabel('$ F_{\lambda}\, \mathrm{[erg\, s^{-1}\, cm^{-2}\, {\AA}^{-1}]}$')
    ax.set_yscale('log')
    return fig, ax


def _weighted_average_flux(fluxes, sigmas):
    weights = np.zeros_like(sigmas, dtype=float)
    positive_sigma = sigmas > 0
    weights[positive_sigma] = 1.0 / (sigmas[positive_sigma]**2)
    weight_sum = np.sum(weights, axis=0)

    flux_aver = np.empty(fluxes.shape[1], dtype=float)
    sigma_aver = np.empty(fluxes.shape[1], dtype=float)

    valid = weight_sum > 0
    flux_aver[valid] = np.sum(weights[:, valid] * fluxes[:, valid], axis=0) / weight_sum[valid]
    sigma_aver[valid] = 1.0 / np.sqrt(weight_sum[valid])

    if np.any(~valid):
        flux_aver[~valid] = np.mean(fluxes[:, ~valid], axis=0)
        if fluxes.shape[0] == 1:
            sigma_aver[~valid] = np.abs(sigmas[0, ~valid])
        else:
            sigma_aver[~valid] = np.std(fluxes[:, ~valid], axis=0)

    return flux_aver, sigma_aver


def _average_iue_files(files, wav_mask, ax=None):
    if len(files) == 0:
        return np.array([]), np.array([]), np.array([])

    wav_ref = None
    flux_rows = []
    sigma_rows = []

    for file in files:
        with fits.open(file) as hdul:
            data = hdul[1].data

        wav = data['WAVELENGTH']
        mask = wav_mask(wav)
        wav_use = wav[mask]
        flux_use = data['FLUX'][mask]
        sigma_use = data['SIGMA'][mask]

        if wav_use.size == 0:
            continue

        if wav_ref is None:
            wav_ref = wav_use
            flux_rows.append(flux_use)
            sigma_rows.append(sigma_use)
            if ax is not None:
                ax.plot(wav_ref, flux_use, lw=0.5, color='grey', alpha=0.6)
        else:
            flux_interp = np.interp(wav_ref, wav_use, flux_use)
            sigma_interp = np.interp(wav_ref, wav_use, sigma_use)
            flux_rows.append(flux_interp)
            sigma_rows.append(sigma_interp)
            if ax is not None:
                ax.plot(wav_ref, flux_interp, lw=0.5, color='grey', alpha=0.6)

    if wav_ref is None or len(flux_rows) == 0:
        return np.array([]), np.array([]), np.array([])

    fluxes = np.asarray(flux_rows, dtype=float)
    sigmas = np.asarray(sigma_rows, dtype=float)
    flux_aver, sigma_aver = _weighted_average_flux(fluxes, sigmas)
    return wav_ref, flux_aver, sigma_aver


def average_iue_swp(files, ax=None):
    return _average_iue_files(files, lambda wav: wav > 1150.0, ax=ax)


def average_iue_lwp(files, ax=None):
    return _average_iue_files(files, lambda wav: (wav > 1980.085) & (wav < 3000.0), ax=ax)


def rebin_spectrum(wav, flux, sigma):
    try:
        rebin, _ = pyasl.binningx0dt(wav, flux, yerr=sigma, nbins=NBINS, x0=wav[0])
        return rebin[:, 0], rebin[:, 1], rebin[:, 2]
    except Exception as exc:
        print('[IUE][WARN] rebin failed ({}), using native sampling'.format(exc))
        return wav, flux, sigma


def save_average_table(star_path, wav_angstrom, flux_flambda):
    os.makedirs(star_path, exist_ok=True)
    dat_outpath = '{}/{}_iue_average_jy.dat'.format(star_path, star_path)
    dat_arr = np.column_stack([wav_angstrom * 1e-4, flux_flambda_to_jy(wav_angstrom, flux_flambda)])
    np.savetxt(
        dat_outpath,
        dat_arr,
        fmt='%.6f %.6e',
        header='wavelength_um flux_jy',
    )
    print('[IUE][OK] wrote {}'.format(dat_outpath))


def process_star(star):
    print('[IUE][STAR] {}'.format(star))
    _, ax = setup_plot()

    swp_files, lwp_files = collect_iue_files(star)
    if len(swp_files) == 0 and len(lwp_files) == 0:
        raise RuntimeError('Missing IUE spectra for {}'.format(star))

    chunks = []
    if len(swp_files) > 0:
        wav1, flux1, sigma1 = average_iue_swp(swp_files, ax=ax)
        if wav1.size > 0:
            chunks.append((wav1, flux1, sigma1))

    if len(lwp_files) > 0:
        wav2, flux2, sigma2 = average_iue_lwp(lwp_files, ax=ax)
        if wav2.size > 0:
            chunks.append((wav2, flux2, sigma2))

    if len(chunks) == 0:
        raise RuntimeError('No usable IUE spectra for {}'.format(star))

    wav = np.concatenate([chunk[0] for chunk in chunks])
    flux = np.concatenate([chunk[1] for chunk in chunks])
    sigma = np.concatenate([chunk[2] for chunk in chunks])

    sort_idx = np.argsort(wav)
    wav = wav[sort_idx]
    flux = flux[sort_idx]
    sigma = sigma[sort_idx]

    wav_rbn, flux_rbn, sigma_rbn = rebin_spectrum(wav, flux, sigma)
    flux_dered = pyasl.unred(wav_rbn, flux_rbn, EBV, R_V=RV)

    yerr = np.where(np.isfinite(sigma_rbn) & (sigma_rbn > 0), sigma_rbn, 0.0)
    ax.errorbar(
        wav_rbn,
        flux_rbn,
        yerr,
        color='grey',
        ecolor='grey',
        lw=1.0,
        elinewidth=0.1,
    )
    ax.errorbar(
        wav_rbn,
        flux_dered,
        yerr,
        color='black',
        ecolor='grey',
        elinewidth=0.1,
    )

    os.makedirs(star, exist_ok=True)
    plt.savefig('{0}/{0}_IUE_ines_average.png'.format(star), format='png')
    save_average_table(star, wav_rbn, flux_rbn)
    # plt.show()
    plt.close()


def main():
    for star in starlist:
        try:
            process_star(star)
        except Exception as exc:
            print('[IUE][ERROR] {}: {}'.format(star, exc))


if __name__ == '__main__':
    main()

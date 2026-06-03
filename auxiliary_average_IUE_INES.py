from astropy.io import fits
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from PyAstronomy import pyasl
import astropy.units as unit


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


def _average_iue_files(files, wav_mask, plots=False):
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
            if plots:
                plt.plot(wav_ref, flux_use, lw=0.5, alpha=0.6)
        else:
            flux_interp = np.interp(wav_ref, wav_use, flux_use)
            sigma_interp = np.interp(wav_ref, wav_use, sigma_use)
            flux_rows.append(flux_interp)
            sigma_rows.append(sigma_interp)
            if plots:
                plt.plot(wav_ref, flux_interp, lw=0.5, alpha=0.6)

    if wav_ref is None or len(flux_rows) == 0:
        return np.array([]), np.array([]), np.array([])

    fluxes = np.asarray(flux_rows, dtype=float)
    sigmas = np.asarray(sigma_rows, dtype=float)
    flux_aver, sigma_aver = _weighted_average_flux(fluxes, sigmas)
    return wav_ref, flux_aver, sigma_aver


def average_iue_swp(swp_files, plots=False):
    return _average_iue_files(swp_files, lambda wav: wav > 1150.0, plots=plots)


def average_iue_lwp(lwp_files, plots=False):
    return _average_iue_files(
        lwp_files,
        lambda wav: (wav > 1980.085) & (wav < 3000.0),
        plots=plots,
    )

###########################################################

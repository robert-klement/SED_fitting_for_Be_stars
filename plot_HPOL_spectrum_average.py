# IMPORTANT:
# This script expects HPOL spectra FITS files (downloaded from MAST) in an HPOL/ or hpol/ subfolder.

import datetime as _dt
from glob import glob
import os

from astropy.io import fits
from astropy.time import Time
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

starlist = ('HD44637', '10CMa', '20Vul', '25Vul', '25Peg', '120Tau', 'epsPsA', 'HR2249',\
            'omeOri', 'upsCyg', 'zetCrv')
# starlist = ('omeOri',)


def _build_wavelength_axis(header, npts):
    crval1 = header.get('CRVAL1')
    cdelt1 = header.get('CDELT1', header.get('CD1_1'))
    crpix1 = float(header.get('CRPIX1', 1.0))

    if crval1 is None or cdelt1 is None:
        raise KeyError('Missing CRVAL1/CDELT1 (or CD1_1) in FITS header.')

    pix = np.arange(npts, dtype=float) + 1.0
    return crval1 + (pix - crpix1) * cdelt1


def _parse_date_obs_to_mjd(date_obs):
    if not date_obs:
        return None

    text = str(date_obs).strip()
    if not text:
        return None

    for fmt in ('isot', 'iso'):
        try:
            return float(Time(text, format=fmt, scale='utc').mjd)
        except Exception:
            pass

    if 'T' in text:
        dpart, tpart = text.split('T', 1)
    else:
        dpart, tpart = text, '00:00:00'

    try:
        if '-' in dpart:
            y, m, d = dpart.split('-')
            dt_obj = _dt.datetime(
                int(y),
                int(m),
                int(d),
                int(tpart.split(':')[0]),
                int(tpart.split(':')[1]),
                int(float(tpart.split(':')[2])),
            )
            return float(Time(dt_obj, scale='utc').mjd)
    except Exception:
        pass

    try:
        if '/' in dpart:
            d, m, y = dpart.split('/')
            dt_obj = _dt.datetime(
                int(y),
                int(m),
                int(d),
                int(tpart.split(':')[0]),
                int(tpart.split(':')[1]),
                int(float(tpart.split(':')[2])),
            )
            return float(Time(dt_obj, scale='utc').mjd)
    except Exception:
        pass

    return None


def _extract_mjd_from_header(header, fallback_mjd):
    if 'MJD-OBS' in header:
        return float(header['MJD-OBS'])
    if 'MJD' in header:
        return float(header['MJD'])
    if 'JD' in header:
        return float(header['JD']) - 2400000.5
    if 'HJD' in header:
        return float(header['HJD'])

    for key in ('DATE-OBS', 'FRAME'):
        if key in header:
            mjd = _parse_date_obs_to_mjd(header[key])
            if mjd is not None:
                return mjd

    return fallback_mjd


def loadfits(fitsfile, data_header=0, ARCES_reduced_by_Peter=False):
    del ARCES_reduced_by_Peter
    mjd_jd2000 = 51544.5

    with fits.open(fitsfile) as hdul:
        data = np.asarray(hdul[data_header].data).squeeze()
        if data.ndim != 1:
            raise ValueError('Expected 1D spectrum in {}, got shape {}.'.format(fitsfile, data.shape))

        header = hdul[0].header
        flux = data
        wl = _build_wavelength_axis(header, flux.size)

        if 'WLSHIFT' in header:
            wl = wl + float(header['WLSHIFT'])

        mjd = _extract_mjd_from_header(header, fallback_mjd=mjd_jd2000)
        dateobs = header.get('DATE-OBS', header.get('FRAME', ''))
        datereduc = header.get('IRAF-TLM', header.get('DATE', ''))

    return wl, flux, mjd, dateobs, datereduc, fitsfile


def collect_hpol_files(star_path, pattern):
    files = glob(star_path + '/hpol/' + pattern)
    files += glob(star_path + '/HPOL/' + pattern)
    return sorted(files)


def flux_flambda_to_jy(wav_angstrom, flux_flambda):
    return flux_flambda * 3.33564095e4 * wav_angstrom**2


def setup_plot():
    fig = plt.figure(figsize=(4, 2.4), dpi=300)
    gs = gridspec.GridSpec(1, 1, height_ratios=[1], width_ratios=[1])
    gs.update(left=0.2, right=0.95, bottom=0.20, top=0.9)
    ax = plt.subplot(gs[0, 0])
    ax.set_xlabel('$\lambda$ [ $\AA$]')
    ax.set_ylabel('$ F_{\lambda}\, \mathrm{[erg\, s^{-1}\, cm^{-2}\, {\AA}^{-1}]}$')
    return fig, ax


def collect_hpol_groups(star_path):
    files_ccd_b = collect_hpol_files(star_path, '*ccd*b_hw.fits')
    files_ccd_r = collect_hpol_files(star_path, '*ccd*r_hw.fits')
    files_ret = collect_hpol_files(star_path, '*ret*.fits')

    files = [files_ccd_b, files_ccd_r, files_ret]
    file_groups = ['ccd_b', 'ccd_r', 'ret']

    if np.sum([len(group) for group in files]) == 0:
        files = [collect_hpol_files(star_path, '*.fits') + collect_hpol_files(star_path, '*.FITS')]
        file_groups = ['all']

    return files, file_groups


def average_hpol_group(files, ax):
    wav_ref = None
    flux_rows = []

    for file in files:
        wav, flux, _, _, _, _ = loadfits(file)
        mask = flux != 0.0
        wav = wav[mask]
        flux = flux[mask]

        if wav.size == 0:
            continue

        if wav_ref is None:
            wav_ref = wav
            flux_rows.append(flux)
        else:
            flux_rows.append(np.interp(wav_ref, wav, flux))

        ax.plot(wav, flux, lw=0.5, color='grey', alpha=0.5)

    if wav_ref is None or len(flux_rows) == 0:
        return np.array([]), np.array([])

    flux_aver = np.average(np.asarray(flux_rows, dtype=float), axis=0)
    return wav_ref, flux_aver


def save_average_table(star_path, wav_angstrom, flux_flambda):
    os.makedirs(star_path, exist_ok=True)
    dat_outpath = '{}/{}_hpol_average_jy.dat'.format(star_path, star_path)
    dat_arr = np.column_stack([wav_angstrom * 1e-4, flux_flambda_to_jy(wav_angstrom, flux_flambda)])
    np.savetxt(
        dat_outpath,
        dat_arr,
        fmt='%.6f %.6e',
        header='wavelength_um flux_jy',
    )
    print('[HPOL][OK] wrote {}'.format(dat_outpath))


def process_star(star):
    print('[HPOL][STAR] {}'.format(star))
    _, ax = setup_plot()
    files, file_groups = collect_hpol_groups(star)

    avg_chunks = []
    for idx, file_group in enumerate(files):
        if len(file_group) == 0:
            print('[HPOL][INFO] {}: no files found'.format(file_groups[idx]))
            continue

        wav_ref, flux_aver = average_hpol_group(file_group, ax)
        if wav_ref.size == 0:
            print('[HPOL][INFO] {}: no usable spectra'.format(file_groups[idx]))
            continue

        ax.plot(wav_ref, flux_aver, lw=0.5, color='black')
        avg_chunks.append(np.column_stack([wav_ref, flux_aver]))

    if len(avg_chunks) == 0:
        raise RuntimeError('[HPOL] No spectra found in {}/hpol or {}/HPOL'.format(star, star))

    avg_data = np.vstack(avg_chunks)
    sort_idx = np.argsort(avg_data[:, 0])
    wav_avg = avg_data[sort_idx, 0]
    flux_avg = avg_data[sort_idx, 1]

    os.makedirs(star, exist_ok=True)
    plt.savefig(star + '/' + star + '_HPOL.png', format='png', dpi=300)
    save_average_table(star, wav_avg, flux_avg)
    # plt.show()
    plt.close()


def main():
    for star in starlist:
        try:
            process_star(star)
        except Exception as exc:
            print('[HPOL][ERROR] {}: {}'.format(star, exc))


if __name__ == '__main__':
    main()

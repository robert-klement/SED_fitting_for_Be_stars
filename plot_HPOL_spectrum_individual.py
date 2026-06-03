import datetime as _dt
from glob import glob

from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

star = 'omeOri'


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


def loadfits(fitsfile, data_header=0):
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

files = sorted(glob(star+'/hpol/*.fits'))

for i in range(len(files)):

    fig = plt.figure(figsize=(4, 2.4), dpi=300)
    gs = gridspec.GridSpec(1, 1, height_ratios=[1], width_ratios=[1])
    gs.update(left=0.2, right=0.95, bottom=0.20, top=0.9)
    ax1 = plt.subplot(gs[0, 0])  # Area for the first plot

    ax1.set_ylim(5e-13, 4e-9)
    ax1.set_yscale('log')
    ax1.set_xlabel('$\lambda$ [ $\AA$]')
    ax1.set_ylabel(
        '$ F_{\lambda}\, \mathrm{[erg\, s^{-1}\, cm^{-2}\, {\AA}^{-1}]}$')


    wl, flux, mjd, dateobs, datereduc, fitsfile = loadfits(files[i])
    wl = wl[(flux != 0.)]
    flux = flux[(flux != 0.)]

    print(fitsfile)

    # if int(dateobs[0:4]) < 1995:
    #     #R = wl[-1] / (wl[-1] - wl[-2])
    #     wav_bin = 25
    # else:
    #     wav_bin = 10

    # ax1.set_title(dateobs)
    # ax1.text(np.mean(wl)*1.1, np.max(flux), '$\lambda$ bin = {} $\AA$'.format(wav_bin))

    ax1.plot(wl, flux, lw=0.5, color='black')


    plt.savefig(files[i]+'.png', format='png', dpi=300)

# plt.savefig(star+'/'+star+'_HPOL.png', format='png', dpi=300)

plt.show()
plt.close()

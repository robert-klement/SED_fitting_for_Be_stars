import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pyhdust as phd
import pyhdust.spectools as spt
import pyhdust.phc as phc
from glob import glob
from PyAstronomy import pyasl
from shutil import copyfile
from astropy.modeling import models, fitting
import os
import astropy.io.fits as fits
import sys
import astropy.units as u
os.environ['PYSYN_CDBS'] = "/home/rklement/Documents/Kurucz_models/"
import pysynphot as S

sys.path.append('/home/rklement/Documents/python/python_tools')
from spec_tools import loadfits
from auxiliary_vieira_et_al_2015_eq20 import calc_n
from auxiliary_average_IUE_INES import *



def readphot(label):
    """Read photometry values (wavelength, flux, error), deredden the flux, and return the values."""
    index = (photometry_labels == label).nonzero()
    wav = photometry_values[index, 0]
    """Flux in erg s-1 cm-2 A-1"""
    # flux = (photometry_values[index, 1] * 3e-13) / (photometry_values[index, 0]**2)
    # err = (photometry_values[index, 2] * 3e-13) / (photometry_values[index, 0]**2)
    """Flux in erg s-1 cm-2 Hz-1"""
    flux = photometry_values[index, 1] * 1e-23
    err = photometry_values[index, 2] * 1e-23
    flux_jy = photometry_values[index, 1]
    err_jy = photometry_values[index, 2]
    wav = wav[0, :]
    flux = flux[0, :]
    err = err[0, :]
    # flux = pyasl.unred(wav*1e4,flux,ebv=ebv)
    # flux_phot_interp = np.interp(np.log10(wav), np.log10(lbd_hdust_phot), np.log10(flux_phot[0,:]))
    return wav, flux, err, flux_jy, err_jy


def readupperlim(label):
    """Read upper limit values (wavelength, flux, error), deredden the flux, and return the values."""
    index = (photometry_labels == label).nonzero()
    wav = photometry_values[index, 0]
    """Flux in erg s-1 cm-2 A-1"""
    # flux = (photometry_values[index, 1] * 3e-13) / (photometry_values[index, 0]**2)
    # err = (photometry_values[index, 2] * 3e-13) / (photometry_values[index, 0]**2)
    """Flux in erg s-1 cm-2 Hz-1"""
    flux = photometry_values[index, 1] * 1e-23
    err = photometry_values[index, 2] * 1e-23
    wav = wav[0, :]
    flux = flux[0, :]
    err = err[0, :]
    # flux = pyasl.unred(wav*1e4,flux,ebv=ebv)
    # flux_phot_interp = np.interp(np.log10(wav), np.log10(lbd_hdust_phot), np.log10(flux_phot[0,:]))
    return wav, flux, err

#################################################################################š

starlist = ['HR2142']

for j, star in enumerate(starlist):
    print(star)

    fig = plt.figure(figsize=(8,4), dpi=300)
    gs = gridspec.GridSpec(1, 1, height_ratios=[1], width_ratios=[1])
    gs.update(left=0.12, right=0.96, bottom=0.16, top=0.92)
    ax1 = plt.subplot(gs[0, 0])  # Area for the first plot

    ax1.set_xlabel('$\lambda\, \mathrm{[\mu m]}$', fontsize=12)
    ax1.set_ylabel(
        '$ F_{\\nu}\, \mathrm{[erg\, s^{-1}\, cm^{-2}\, \mathrm{Hz}^{-1}]}$', fontsize=12)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(1.0e-1, 0.3e2)
    ax1.set_ylim(1e-25, 1e-20)
    ax1.tick_params(axis='both', which='major', labelsize=11)


    # IUE
    swp_files = sorted(glob('{}/IUE/*SWP*.FITS'.format(star)))
    lwp_files = sorted(glob('{}/IUE/*LW*.FITS'.format(star)))

    wav1, flux1, sigma1 = average_iue_swp(swp_files, plots=True)
    wav2, flux2, sigma2 = average_iue_lwp(lwp_files, plots=True)

    wav, flux, sigma = np.concatenate((wav1, wav2)), np.concatenate(
        (flux1, flux2)), np.concatenate((sigma1, sigma2))

    rebin, wav_bin = pyasl.binningx0dt(wav, flux, yerr=sigma, nbins=170, x0=wav[0])
    wav_rbn, flux_rbn, sigma_rbn = rebin[::, 0], rebin[::, 1], rebin[::, 2]
    wl, flux = wav_rbn, flux_rbn

    flux_jy = []
    for i in range(len(flux)):
        flux_jy.append((3.33564095e4 * flux[i] * ((wl[i])**2)))
    # ax1.plot(wl*1e-4, np.array(flux_jy)*1e-23, color='grey', lw=0.8)

    ax1.plot(wl*1e-4, np.array(flux_jy)*1e-23,
            color='black', lw=1.0, label='IUE')

    # # HPOL

    # files_ccd_b = sorted(glob(star+'/HPOL/*ccd*b_hw.fits'))
    # files_ccd_r = sorted(glob(star+'/HPOL/*ccd*r_hw.fits'))
    # files_ret = sorted(glob(star+'/HPOL/*ret*.fits'))

    # files_hpol = [files_ccd_b, files_ccd_r, files_ret]

    # for f, filetype in enumerate(files_hpol):
    #     for i, file in enumerate(filetype):

    #         wav, flux, mjd, dateobs, datereduc, fitsfile = spt.loadfits(file)
    #         wav = wav[(flux != 0.)]
    #         flux = flux[(flux != 0.)]

    #         if i == 0:
    #             wav_ref = wav
    #             fluxes = np.empty((len(filetype), len(wav_ref)))
    #             fluxes[i, :] = flux
    #         else:
    #             fluxes[i, :] = np.interp(wav_ref, wav, flux)

    #         flux_aver = np.empty(len(wav_ref))
    #         for i in range(len(wav_ref)):
    #             flux_aver[i] = np.average(fluxes[:, i])

    #         # ax1.plot(wav*1e-4, flux, lw=0.5, color='grey', alpha=0.5)
    #     if f == 0:
    #         ax1.plot(wav_ref*1e-4, flux_aver*(((wav_ref*1e-4)**2)/3e10), lw=1.0, color='red', label='HPOL CCD')
    #     elif f == 1:
    #         ax1.plot(wav_ref*1e-4, flux_aver*(((wav_ref*1e-4)**2)/3e10), lw=1.0, color='red')
    #     elif f == 2:
    #         ax1.plot(wav_ref*1e-4, flux_aver*(((wav_ref*1e-4)**2)/3e10), lw=1.0, ls='--', color='blue', label='HPOL Reticon')

    # photometric data
    photometry_values = np.loadtxt(
        '{0}/{0}_photometry_manual.dat'.format(star), usecols=[0, 1, 2])
    photometry_labels = np.loadtxt(
        '{0}/{0}_photometry_manual.dat'.format(star), usecols=[3], dtype=np.str)

    phot_labels = ['simbad', 'dougherty+1991', 'iras', 'wise', 'akari']
    legend_labels = ['Ducati+2002', 'Dougherty+1991', 'IRAS', 'WISE', 'Akari']

    # markers

    markers = ['o','^','x','s','v']
    # markers = ['' for x in range(len(phot_labels))]
    # for i in range(0, len(phot_labels)):
    #     markers[i] = str(phc.cycles(i=i, ctype='mk'))
    msize = 8


    # colors
    cmap = 'gist_rainbow'
    cmap2 = 'nipy_spectral'
    colors = phc.gradColor(np.arange(6),
                           vmin=0, vmax=6 - 1, cmapn=cmap2)
    # colors = ('red','blue','green','cyan')

    laboca_wl, laboca_flux, laboca_err = [], [], []
    jcmt_wl, jcmt_flux, jcmt_err = [], [], []
    iram_wl, iram_flux, iram_err = [], [], []
    karma_wl, karma_flux, karma_err = [], [], []

    for i in range(len(phot_labels)):
        photometry = readphot(phot_labels[i])

        leg1 = ax1.errorbar(photometry[0], photometry[1], yerr=photometry[2], ls='None',
                            markersize=msize, marker=markers[i], color=colors[i], markeredgecolor=colors[i], markerfacecolor='None', label=legend_labels[i])
        # leg1 = ax1.errorbar(photometry[0], photometry[1], ls='None', markersize=msize,
        #                     marker=markers[i], color=colors[i], markeredgecolor='none', label=legend_labels[i])

        ax1.legend(fontsize='8', ncol=2, numpoints=1, loc='upper right')

        if phot_labels[i] == 'simbad':
            simbad_wl, simbad_flux = photometry[0], photometry[1]
            # plaw_phot = models.PowerLaw1D(amplitude=simbad_flux[-3], x_0=simbad_wl[-3], alpha=2)
            # x = np.arange(3, 1e5, 10)
            # ax1.plot(x, plaw_phot(x), color='r', ls=':')

        if phot_labels[i] == 'dougherty+1991':
            dougherty_wl, dougherty_flux, dougherty_err = photometry[0], photometry[1], photometry[2]

        if phot_labels[i] == 'iras':
            if len(photometry[0]) > 0:
                iras_wl, iras_flux, iras_err = photometry[0], photometry[1], photometry[2]

                if len(photometry[0][:]) > 1:
                    plaw_init = models.PowerLaw1D(
                        amplitude=iras_flux[0], x_0=iras_wl[0], alpha=1)
                    plaw_init.x_0.fixed = True
                    fit = fitting.LevMarLSQFitter()
                    plaw_fn0 = fit(plaw_init, iras_wl,
                                   iras_flux, weights=1 / iras_err)

                    kappa_iras = plaw_fn0.alpha.value

                    n_iras = calc_n(kappa_iras + 2.0, 15)

                    # ax1.text(4e3, 3e-22, 'IRAS: $\\kappa$ = {:.2f}, $n$ = {:.2f}'.format(
                    #     kappa_iras, n_iras[2]), fontsize=14)
                    x = np.arange(3, 1e5, 10)
                    # ax1.plot(x, plaw_fn0(x),
                    #          color=colors[i], ls='--', alpha=0.5)
            else:
                iras_wl, iras_flux, iras_err = [], [], []
                kappa_iras = 0.

        if phot_labels[i] == 'wise':
            if len(photometry[0]) > 0:
                wise_wl, wise_flux, wise_err = photometry[0], photometry[1], photometry[2]

                if len(photometry[0][:]) > 1:
                    plaw_init = models.PowerLaw1D(
                        amplitude=wise_flux[0], x_0=wise_wl[0], alpha=1)
                    plaw_init.x_0.fixed = True
                    fit = fitting.LevMarLSQFitter()
                    plaw_fn0 = fit(plaw_init, wise_wl,
                                   wise_flux, weights=1 / wise_err)

                    kappa_wise = plaw_fn0.alpha.value

                    n_wise = calc_n(kappa_wise + 2.0, 15)

                    # ax1.text(4e3, 5e-23, 'WISE: $\\kappa$ = {:.2f}, $n$ = {:.2f}'.format(
                    #     kappa_wise, n_wise[2]), fontsize=14)
                    x = np.arange(3, 80, 10)
                    # ax1.plot(x, plaw_fn0(x),
                    #          color=colors[i], ls='--', alpha=0.5)
            else:
                wise_wl, wise_flux, wise_err = [], [], []
                kappa_wise = 0.

        if phot_labels[i] == 'akari':
            if len(photometry[0]) > 0:
                akari_wl, akari_flux, akari_err = photometry[0], photometry[1], photometry[2]

                if len(photometry[0][:]) > 1:
                    plaw_init = models.PowerLaw1D(
                        amplitude=akari_flux[0], x_0=akari_wl[0], alpha=1)
                    plaw_init.x_0.fixed = True
                    fit = fitting.LevMarLSQFitter()
                    plaw_fn0 = fit(plaw_init, akari_wl,
                                   akari_flux, weights=1 / akari_err)

                    kappa_akari = plaw_fn0.alpha.value

                    n_akari = calc_n(kappa_akari + 2.0, 15)

                    # ax1.text(4e3, 1.2e-22, 'AKARI: $\\kappa$ = {:.2f}, $n$ = {:.2f}'.format(
                    #     kappa_akari, n_akari[2]), fontsize=14)
                    x = np.arange(3, 1e5, 10)
                    # ax1.plot(x, plaw_fn0(x),
                    #          color=colors[i], ls='--', alpha=0.5)
            else:
                akari_wl, akari_flux, akari_err = [], [], []
                kappa_akari = 0.

        if (phot_labels[i] == 'laboca'):
            if len(photometry[0]) > 0:
                laboca_wl, laboca_flux, laboca_err = photometry[0], photometry[1], photometry[2]

        if (phot_labels[i] == 'waters+1991'):
            if len(photometry[0]) > 0:
                jcmt_wl, jcmt_flux, jcmt_err = photometry[0], photometry[1], photometry[2]

        if (phot_labels[i] == 'wendker+2000'):
            if len(photometry[0]) > 0:
                iram_wl, iram_flux, iram_err = photometry[0], photometry[1], photometry[2]

        if (phot_labels[i] == 'karma'):
            if len(photometry[0]) > 0:
                karma_wl, karma_flux, karma_err = photometry[0], photometry[1], photometry[2]

        if phot_labels[i] == 'taylor+1990':
            if len(photometry[0]) > 0:
                vla_wl, vla_flux, vla_err = photometry[0], photometry[1], photometry[2]

                if len(photometry[0][:]) > 1:
                    plaw_init = models.PowerLaw1D(
                        amplitude=photometry[1][0], x_0=photometry[0][0], alpha=1)
                    plaw_init.x_0.fixed = True
                    fit = fitting.LevMarLSQFitter()
                    plaw_fn = fit(
                        plaw_init, photometry[0], photometry[1], weights=1 / vla_err)
                    ax1.text(
                        4e3, 1.5e-23, 'VLA: $\\kappa$ = {:.2f}'.format(plaw_fn.alpha.value), fontsize=14)
                    x = np.arange(3, 1e5, 10)
                    ax1.plot(x, plaw_fn(x),
                             color=colors[i], ls='--', alpha=0.5)
            else:
                vla_wl, vla_flux, vla_err = [], [], []

        if phot_labels[i] == 'jvla':
            if len(photometry[0]) > 0:
                jvla_wl, jvla_flux, jvla_err = photometry[0], photometry[1], photometry[2]

                if len(photometry[0][:]) > 1:
                    plaw_init = models.PowerLaw1D(
                        amplitude=photometry[1][0], x_0=photometry[0][0], alpha=1)
                    plaw_init.x_0.fixed = True
                    fit = fitting.LevMarLSQFitter()
                    plaw_fn = fit(
                        plaw_init, photometry[0], photometry[1], weights=1 / jvla_err)
                    ax1.text(
                        4e3, 6e-24, 'JVLA: $\\kappa$ = {:.2f}'.format(plaw_fn.alpha.value), fontsize=14)
                    x = np.arange(3, 1e5, 10)
                    ax1.plot(x, plaw_fn(x),
                             color=colors[i], ls='--', alpha=0.5)
            else:
                jvla_wl, jvla_flux, jvla_err = [], [], []

    # FIT ALL IR DATA
    wav_IR = np.concatenate((dougherty_wl[-2:], wise_wl, akari_wl))
    # print(dougherty_wl[-2:])
    flux_IR = np.concatenate((dougherty_flux[-2:], wise_flux, akari_flux))
    err_IR = np.concatenate((dougherty_err[-2:], wise_err, akari_err))

    flux_IR_Flambda = flux_IR * (2.997924e14 / (wav_IR)**2)
    err_IR_Flambda = err_IR * (2.997924e14 / (wav_IR)**2)

    if len(wav_IR) > 1:
        # F_nu
        plaw_init = models.PowerLaw1D(
            amplitude=flux_IR[0], x_0=wav_IR[0], alpha=1)
        plaw_init.x_0.fixed = True
        fit = fitting.LevMarLSQFitter()
        plaw_fn0 = fit(plaw_init, wav_IR, flux_IR, weights=1 / err_IR)

        kappa_IR = plaw_fn0.alpha.value
        print('kappa_IR = ', kappa_IR)

        # linfit_init = models.Linear1D(slope=-3, intercept=-15)
        # fit = fitting.LevMarLSQFitter()
        # linfit_fn0 = fit(linfit_init, np.log(wav_IR), np.log(flux_IR), weights=1 / err_IR)

        # kappa_IR_log = linfit_fn0.slope.value

        # F_lambda
        # plaw_init = models.PowerLaw1D(amplitude=flux_IR_Flambda[0], x_0=wav_IR[0], alpha=1)
        # plaw_init.x_0.fixed = True
        # fit = fitting.LevMarLSQFitter()
        # plaw_fn = fit(plaw_init, wav_IR, flux_IR_Flambda, weights=1 / err_IR)

        # alpha_IR = plaw_fn.alpha.value

        # linfit_init = models.Linear1D(slope=-3, intercept=-15)
        # fit = fitting.LevMarLSQFitter()
        # linfit_fn = fit(linfit_init, np.log(wav_IR), np.log(flux_IR_Flambda), weights=-1 / err_IR_Flambda)

        # alpha_IR_log = linfit_fn.slope.value

        n_IR = calc_n(kappa_IR + 2.0, 15)
        print('n_IR = ', n_IR[2])

        x = np.array([3, 28])
        # plt.text(4e3, 1e-24, 'IR missions: $\\kappa$ = {:.2f}, $n$ = {:.2f}'.format(
        #     kappa_IR, n_IR[2]), fontsize=14, weight='bold', color='green')
        # plt.plot(x, np.exp(linfit_fn0.intercept) * x**(kappa_IR_log), color='green', ls=':')
        plt.plot(x, plaw_fn0(x), color='grey', ls=':', lw=1.0)

        # rms = np.array([0.008, 0.006])
        # for k in range(1, 3):
        #ax1.text(4e3, 1e-25 - (k - 1) * 9e-26, '{:.2f} cm - {:.4f} mJy - {:.2f}'.format(x[k] * 1e-4, plaw_fn0(x[k]) * 1e26, (plaw_fn0(x[k]) * 1e26) / rms[k - 1]))


    distance = 403

    ###
    sp = S.Icat('ck04models', 22500, 0.0, 3.35) # HR 2142
    sp = sp * S.Extinction(0.14, 'mwavg')
    sp.convert('fnu')

    factor = (((7.06*u.Rsun) / (distance*u.pc))**2).decompose()
    flux_nu_Be = sp.flux * factor

    ax1.plot(sp.wave*1e-4, flux_nu_Be, color='grey', lw=0.5)

    # factor = (((7.5*u.Rsun) / (distance*u.pc))**2).decompose()
    # flux_nu_Be = sp.flux * factor

    # ax1.plot(sp.wave*1e-4, flux_nu_Be, color='grey', lw=0.5)

    # factor = (((8.0*u.Rsun) / (distance*u.pc))**2).decompose()
    # flux_nu_Be = sp.flux * factor

    # ax1.plot(sp.wave*1e-4, flux_nu_Be, color='grey', lw=0.5)

    # ###
    # sp = S.Icat('ck04models', 35000, 0.0, 5.00) # sdO
    # sp = sp * S.Extinction(0.08, 'mwavg')
    # sp.convert('fnu')

    # factor = (((0.62*u.Rsun) / (distance*u.pc))**2).decompose()
    # flux_nu = sp.flux * factor

    # ax1.plot(sp.wave*1e-4, flux_nu, color='purple', lw=0.5)

    ###
    sp = S.Icat('k93models', 43000, 0.0, 5.0) # sdO
    sp = sp * S.Extinction(0.14, 'mwavg')
    sp.convert('fnu')

    factor = (((0.25*u.Rsun) / (distance*u.pc))**2).decompose()
    flux_nu = sp.flux * factor

    ax1.plot(sp.wave*1e-4, flux_nu, color='purple', lw=0.5)

    ###
    sp = S.Icat('k93models', 43000, 0.0, 5.0) # sdO
    sp = sp * S.Extinction(0.14, 'mwavg')
    sp.convert('fnu')

    factor = (((0.6*u.Rsun) / (distance*u.pc))**2).decompose()
    flux_nu = sp.flux * factor

    ax1.plot(sp.wave*1e-4, flux_nu, color='purple', lw=0.5)



    # ax1.axvspan(xmin=0.1150, xmax=0.1950, color='gray',alpha=0.2)
    # ax1.axvspan(xmin=1.63-0.307/2, xmax=1.63+0.307/2, color='gray',alpha=0.2)
    # ax1.axvspan(xmin=1.95, xmax=2.4, color='gray',alpha=0.2)

    ax1.axhline(y=2.6552e-24, xmin=0.02, xmax=0.03, color='r',ls='-', lw=0.6)
    ax1.axhline(y=1.5341e-23, xmin=0.02, xmax=0.03, color='r',ls='-', lw=0.6)


    ax1.axhline(y=1.1340e-24, xmin=0.475, xmax=0.495, color='k',ls='--', lw=0.6)


    # plt.savefig('{}/{}_SED.png'.format(star, starlist[j]), format='png')
    # plt.savefig('{}/{}_SED.pdf'.format(star, starlist[j]), format='pdf')

    plt.show()
    plt.close()

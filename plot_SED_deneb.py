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
sys.path.append('/home/robert/Documents/python/python_tools')
from spec_tools import loadfits
from auxiliary_vieira_et_al_2015_eq20 import calc_n
from auxiliary_average_IUE_INES import *
import astropy.units as u
os.environ['PYSYN_CDBS'] = "/home/robert/Documents/Kurucz_models/"
import pysynphot as S



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

starlist = ['deneb']
# fnames = ['HR2142']

for j, star in enumerate(starlist):
    print(star)

    fig = plt.figure(figsize=(8,4), dpi=300)
    gs = gridspec.GridSpec(1, 1, height_ratios=[1], width_ratios=[1])
    gs.update(left=0.12, right=0.96, bottom=0.16, top=0.92)
    ax1 = plt.subplot(gs[0, 0])  # Area for the first plot

    ax1.set_xlabel('$\lambda\, \mathrm{[\mu m]}$', fontsize=15)
    ax1.set_ylabel(
        '$ F_{\\nu}\, \mathrm{[erg\, s^{-1}\, cm^{-2}\, \mathrm{Hz}^{-1}]}$', fontsize=15)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(1.1e-1, 100000)
    ax1.set_ylim(5e-29, 1e-19)
    ax1.tick_params(axis='both', which='major', labelsize=11)


    # IUE
    swp_files = sorted(glob('{}/IUE/SWP*.FITS'.format(star)))
    lwp_files = sorted(glob('{}/IUE/LW*.FITS'.format(star)))

    wav1, flux1, sigma1 = average_iue_swp(swp_files, plots=False)
    wav2, flux2, sigma2 = average_iue_lwp(lwp_files, plots=False)

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


    # photometric data
    photometry_values = np.loadtxt(
        '{0}/{0}_photometry_manual.dat'.format(star), usecols=[0, 1, 2])
    photometry_labels = np.loadtxt(
        '{0}/{0}_photometry_manual.dat'.format(star), usecols=[3], dtype=np.str)

    phot_labels = ['ducati+2002', 'cobe', 'irtf', 'msx', 'isophot', 'hht', 'scuba', 'vla']
    legend_labels = ['Ducati+2002', 'COBE/DIRBE', 'IRTF', 'MSX', 'ISOPHOT', 'HHT', 'SCUBA', 'VLA']


    # markers

    markers = ['o','^','x','s','v','<','>','1']
    # markers = ['' for x in range(len(phot_labels))]
    # for i in range(0, len(phot_labels)):
    #     markers[i] = str(phc.cycles(i=i, ctype='mk'))
    msize = 8


    # colors
    cmap = 'gist_rainbow'
    cmap2 = 'nipy_spectral'
    colors = phc.gradColor(np.arange(8),
                           vmin=0, vmax=8 - 1, cmapn=cmap2)
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

        ax1.legend(fontsize='9', ncol=2, numpoints=1)

        uplim_labels = phot_labels

    for i in range(len(phot_labels)):
        uplim_labels[i] = uplim_labels[i] + '_upper'

    for i in range(len(uplim_labels)):
        photometry = readupperlim(uplim_labels[i])

        ax1.errorbar(photometry[0], photometry[1], yerr=0, ls='None', markersize=msize, marker=markers[i], color=colors[i], alpha=0.3)
        ax1.errorbar(photometry[0], photometry[1], yerr=photometry[2], ls='None', markersize=msize, marker=markers[i], color=colors[i], uplims=True)  # , alpha=0.3)


    sp = S.Icat('ck04models', 8525, -0.20, 1.5) # kap Dra
    sp = sp * S.Extinction(0.04, 'mwavg')
    sp.convert('fnu')

    factor = (((203*u.Rsun) / (802*u.pc))**2).decompose()
    flux_nu = sp.flux * factor

    ax1.plot(sp.wave*1e-4, flux_nu, color='purple', lw=0.5)


    plt.savefig('{}/{}_SED.png'.format(star, starlist[j]), format='png')
    # plt.savefig('{}/{}_plotphot.pdf'.format(star, starlist[j]), format='pdf')


    plt.close()

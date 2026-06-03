from astropy.io import fits
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from PyAstronomy import pyasl

star = 'omeOri'

# PyAstronomy uses np.NAN internally; keep compatibility with newer NumPy.
if not hasattr(np, 'NAN'):
    np.NAN = np.nan

bin_size = 0.2 # A, high disp.
bin_size = 6  # A, low disp.

# swp_files = sorted(glob('{}/IUE_INES_Rebin_SWP/*SWP*.FITS'.format(star)))
# lwp_files = sorted(glob('{}/IUE_INES_Rebin_LW/*LW*.FITS'.format(star)))
# print('SWP files: {}, LWP/LWR files: {}\n-----------'.format(len(swp_files), len(lwp_files)))
# files = swp_files + lwp_files

files = sorted(glob('{}/IUE/*.FITS'.format(star)))

for i in range(len(files)):

    print(os.path.split(files[i])[1])

    with fits.open(files[i]) as fits_file:
        data = fits_file[1].data

    wav = data['WAVELENGTH']
    flux = data['FLUX']
    sigma = data['SIGMA']

    fig = plt.figure(figsize=(4, 2.4), dpi=300)
    gs = gridspec.GridSpec(1, 1, height_ratios=[1], width_ratios=[1])
    gs.update(left=0.25, right=0.9, bottom=0.20, top=0.9)
    ax = plt.subplot(gs[0, 0])  # Area for the first plot

    rebin = int((wav[-1]-wav[0])/bin_size)

    if rebin>0:
        binned, dt = pyasl.binningx0dt(wav, flux, yerr=sigma, nbins=rebin, x0=wav[0])
        ax.errorbar(binned[::,0], binned[::,1], binned[::,2], lw=0.5, label=os.path.split(files[i])[1])
    else:
        ax.errorbar(wav, flux, sigma, lw=0.5, label=os.path.split(files[i])[1])

    ax.set_xlabel('$\lambda$ [ $\AA$]')
    ax.set_ylabel(
        '$ F_{\lambda}\, \mathrm{[erg\, s^{-1}\, cm^{-2}\, {\AA}^{-1}]}$')

    ax.legend(fontsize=8)

    ax.set_yscale('log')
    ax.set_ylim(1e-11,1.5e-9)

    plt.savefig('{}.png'.format(files[i]), format='png')

    # plt.show()
    plt.close()

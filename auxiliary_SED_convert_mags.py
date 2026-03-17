import os

BANDPASSES_DIR = os.path.join(os.getcwd(), 'bandpasses')
JOHNSON_FILTER_INDEX = {
    'U': 0,
    'B': 1,
    'V': 2,
    'R': 3,
    'I': 4,
    'J': 5,
    'H': 6,
    'K': 7,
    'L': 8,
    'M': 9,
    'N': 10
}
WISE_FILTER_INDEX = {'W1': 0, 'W2': 1, 'W3': 2, 'W4': 3}


def _load_zero_flux_table(filename):
    filters_all = []
    wavs = []
    zero_fluxes = []

    with open(os.path.join(BANDPASSES_DIR, filename), 'r') as zfluxes:
        for line in zfluxes.readlines():
            line = line.strip()
            columns = line.split()
            filters_all.append(columns[0])
            wavs.append(columns[1])
            zero_fluxes.append(columns[3])

    return filters_all, wavs, zero_fluxes


def _convert_johnson(star, input_suffix, output_suffix, source_label):
    filters_all, wavs, zero_fluxes = _load_zero_flux_table('johnson_zero_fluxes.dat')

    filters = []
    mags = []
    magerrs = []

    with open(star + '/' + star + '_' + input_suffix + '.dat', 'r') as jphot:
        for line in jphot.readlines():
            line = line.strip()
            columns = line.split()
            filters.append(columns[0])
            mags.append(columns[1])
            magerrs.append(columns[2])

    wav = []
    flux = []
    err = []

    for i in range(len(filters)):
        if filters[i] in JOHNSON_FILTER_INDEX:
            idx = JOHNSON_FILTER_INDEX[filters[i]]
            wav.append(float(wavs[idx]))
            flux.append(10**(-0.4 * float(mags[i])) * float(zero_fluxes[idx]))
            if magerrs[i] == '--':
                err.append(0.1 * flux[i])
            else:
                fluxerr1 = 10**(-0.4 * (float(mags[i]) + float(magerrs[i]))) * float(zero_fluxes[0])
                fluxerr2 = 10**(-0.4 * (float(mags[i]) - float(magerrs[i]))) * float(zero_fluxes[0])
                err.append(abs(fluxerr1 - fluxerr2) / 2.)

    for i in range(len(flux)):
        flux[i] = round(flux[i], 2)
        if i < 5:
            err[i] = round(flux[i] * 0.1, 2)
        else:
            err[i] = round(err[i], 2)

    outpath = star + '/' + star + '_' + output_suffix + '.dat'
    exists = os.path.isfile(outpath)
    if exists:
        os.remove(outpath)

    with open(outpath, 'a') as photfile:
        for i in range(len(filters)):
            photfile.write(('%s %s %s %s\n' % (wav[i], flux[i], err[i], source_label)))

def convert_johnson_simbad(star):
    _convert_johnson(star, 'simbad_johnson', 'simbad_johnson_jy', 'simbad')


def convert_johnson_dougherty(star):
    _convert_johnson(star, 'dougherty_johnson', 'dougherty_johnson_jy', 'dougherty+1991')


def convert_wise(star):
    filters_all, wavs, zero_fluxes = _load_zero_flux_table('wise_zero_fluxes.dat')

    filters = []
    mags = []
    magerrs = []

    raw_flags = []

    with open(star + '/' + star + '_vizier_wise_mags.dat', 'r') as jphot:
        for line in jphot.readlines():
            line = line.strip()
            columns = line.split()
            filters.append(columns[0])
            mags.append(columns[1])
            magerrs.append(columns[2])
            raw_flags.append(columns[3])

    wav = []
    flux = []
    err = []
    flags = []

    for i in range(len(filters)):
        if filters[i] in WISE_FILTER_INDEX:
            idx = WISE_FILTER_INDEX[filters[i]]
            wav.append(float(wavs[idx]))
            flux.append(10**(-0.4 * float(mags[i])) * float(zero_fluxes[idx]))
            fluxerr1 = 10**(-0.4 * (float(mags[i]) + float(magerrs[i]))) * float(zero_fluxes[idx])
            fluxerr2 = 10**(-0.4 * (float(mags[i]) - float(magerrs[i]))) * float(zero_fluxes[idx])
            err.append(abs(fluxerr1 - fluxerr2) / 2.)
            flags.append(raw_flags[i])

    return wav, flux, err, flags

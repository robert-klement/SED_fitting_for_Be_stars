from astroquery.simbad import Simbad
import numpy as np
import os
import re

from auxiliary_SED_convert_mags import convert_johnson_simbad


starlist = ('HD44637', '10CMa', '20Vul', '25Vul', '25Peg', '120Tau', 'epsPsA', 'HR2249',\
            'omeOri', 'upsCyg', 'zetCrv')
# starlist = ('betCMi',)

print('[SIMBAD] Processing {} stars'.format(len(starlist)))


def simbad_name_candidates(star):
    """Return plausible SIMBAD identifiers for compact internal names."""
    candidates = [star]

    # Common compact style in this repo: omeOri, epsPsA, 10CMa -> ome Ori, eps PsA, 10 CMa
    split_idx = None
    for i in range(1, len(star)):
        if star[i].isupper() and (star[i - 1].islower() or star[i - 1].isdigit()):
            split_idx = i
            break
    spaced = star[:split_idx] + ' ' + star[split_idx:] if split_idx is not None else star
    if spaced != star:
        candidates.append(spaced)

    # e.g., ANCol -> AN Col
    m_caps = re.match(r'^([A-Z]{2,})([A-Z][a-z].*)$', star)
    if m_caps:
        candidates.append('{} {}'.format(m_caps.group(1), m_caps.group(2)))

    if '_' in star:
        candidates.append(star.replace('_', ' '))

    # Preserve order while removing duplicates.
    deduped = []
    seen = set()
    for name in candidates:
        if name not in seen:
            deduped.append(name)
            seen.add(name)
    return deduped


def query_object_with_aliases(custom_simbad, star, **kwargs):
    for candidate in simbad_name_candidates(star):
        try:
            result = custom_simbad.query_object(candidate, **kwargs)
        except Exception:
            continue
        if result is not None and len(result) > 0:
            return result, candidate

    return None, star

for star in starlist:
    print('[SIMBAD][STAR] {}'.format(star))

    if not os.path.isdir(star):
        os.makedirs(star)

    # SIMBAD Johnson magnitudes (new TAP interface)
    johnson_filters = ('U', 'B', 'V', 'R', 'I', 'J', 'H', 'K')

    custom_simbad = Simbad()
    custom_simbad.reset_votable_fields()
    custom_simbad.add_votable_fields(*johnson_filters)
    fluxes, resolved_name = query_object_with_aliases(custom_simbad, star)
    if resolved_name != star:
        print('[SIMBAD][INFO] {} resolved as {}'.format(star, resolved_name))

    # flux_error(<filter>) no longer exists; use flux table and extract flux_err by filter
    custom_simbad = Simbad()
    custom_simbad.reset_votable_fields()
    custom_simbad.add_votable_fields('flux')
    errors, _ = query_object_with_aliases(
        custom_simbad, resolved_name, criteria="filter IN ('U', 'B', 'V', 'R', 'I', 'J', 'H', 'K')"
    )

    outpath = star + '/' + star + '_simbad_johnson.dat'
    if os.path.isfile(outpath):
        os.remove(outpath)

    err_by_filter = {}
    if errors is not None and len(errors) > 0:
        err_col = None
        if 'flux_err' in errors.colnames:
            err_col = 'flux_err'
        elif 'flux.error' in errors.colnames:
            err_col = 'flux.error'

        filter_col = 'flux.filter' if 'flux.filter' in errors.colnames else None

        if err_col is not None and filter_col is not None:
            for row in errors:
                filt = str(row[filter_col]).strip()
                err_val = row[err_col]
                if filt in johnson_filters and str(err_val) != '--':
                    # Keep the smallest reported error if multiple measurements exist.
                    if filt not in err_by_filter or float(err_val) < float(err_by_filter[filt]):
                        err_by_filter[filt] = err_val

    mags_written = 0
    has_flux_rows = (fluxes is not None) and (len(fluxes) > 0)
    with open(outpath, 'a') as photfile:
        if has_flux_rows:
            for filt in johnson_filters:
                if filt in fluxes.colnames:
                    mag = fluxes[filt][0]
                    if str(mag) != '--':
                        mags_written += 1
                        err = err_by_filter.get(filt, '--')
                        photfile.write('%s %s %s\n' % (filt, mag, err))
    print('[SIMBAD][INFO] Johnson mags: {} -> {}'.format(mags_written, outpath))

    jy_outpath = star + '/' + star + '_simbad_johnson_jy.dat'
    if mags_written == 0:
        print('[SIMBAD][INFO] no SIMBAD Johnson photometry for {}, skipping Jy conversion'.format(star))
        if os.path.isfile(jy_outpath):
            os.remove(jy_outpath)
        continue

    convert_johnson_simbad(star)

    # Keep these reads to preserve current script behavior/workflow.
    lbd = np.loadtxt(jy_outpath, usecols=(0))
    flux = np.loadtxt(jy_outpath, usecols=(1))
    err = np.loadtxt(jy_outpath, usecols=(2))
    npts = np.atleast_1d(lbd).size
    print('[SIMBAD][OK] Converted to Jy: {} point(s) -> {}/{}_simbad_johnson_jy.dat'.format(
        npts, star, star
    ))

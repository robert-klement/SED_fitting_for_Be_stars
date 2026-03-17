from astroquery.vizier import Vizier
import numpy as np
import os
import re

from auxiliary_SED_convert_mags import convert_wise
from auxiliary_SED_color_correction_IR_missions import color_corr


starlist = ('HD44637', '10CMa', '20Vul', '25Vul', '25Peg', '120Tau', 'epsPsA', 'HR2249',\
            'omeOri', 'upsCyg', 'zetCrv')
# starlist = ('betCMi',)

print('[WISE] Processing {} stars'.format(len(starlist)))


def vizier_name_candidates(star):
    candidates = [star]
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

    deduped = []
    seen = set()
    for name in candidates:
        if name not in seen:
            deduped.append(name)
            seen.add(name)
    return deduped


def query_vizier_with_aliases(vizier, method_name, star, **kwargs):
    method = getattr(vizier, method_name)
    last_result = []
    for candidate in vizier_name_candidates(star):
        try:
            result = method(candidate, **kwargs)
        except Exception:
            continue
        if len(result) > 0:
            return result, candidate
        last_result = result
    return last_result, star

for star in starlist:
    print('[WISE][STAR] {}'.format(star))
    resolved_alias = None

    if os.path.isfile('{0}/{0}_vizier_wise_jy_cc0.dat'.format(star)):
        os.remove('{0}/{0}_vizier_wise_jy_cc0.dat'.format(star))

    if not os.path.isdir('{}'.format(star)):
        os.makedirs('{}'.format(star))

    # WISE
    v = Vizier(columns=['W1mag', 'W2mag', 'W3mag', 'W4mag'])
    mags, used_name = query_vizier_with_aliases(
        v, 'query_region', star, radius='0d0m15s', catalog=['II/328/allwise']
    )
    if used_name != star and resolved_alias is None:
        resolved_alias = used_name
        print('[WISE][INFO] {} resolved as {}'.format(star, resolved_alias))
    if len(mags) == 0:
        mags = []
    else:
        mags = mags[0]

    v = Vizier(columns=['e_W1mag', 'e_W2mag', 'e_W3mag', 'e_W4mag'])
    errors, used_name = query_vizier_with_aliases(
        v, 'query_region', star, radius='0d0m15s', catalog=['II/328/allwise']
    )
    if used_name != star and resolved_alias is None:
        resolved_alias = used_name
        print('[WISE][INFO] {} resolved as {}'.format(star, resolved_alias))
    if len(errors) == 0:
        errors = []
    else:
        errors = errors[0]

    v = Vizier(columns=['ccf', 'ex', 'var', 'qph'])
    flags, used_name = query_vizier_with_aliases(
        v, 'query_region', star, radius='0d0m15s', catalog=['II/328/allwise']
    )
    if used_name != star and resolved_alias is None:
        resolved_alias = used_name
        print('[WISE][INFO] {} resolved as {}'.format(star, resolved_alias))
    if len(flags) == 0:
        flags = []
    else:
        flags = flags[0]

    if len(mags) != 0:
        mags_outpath = '{0}/{0}_vizier_wise_mags.dat'.format(star)
        if os.path.isfile(mags_outpath):
            os.remove(mags_outpath)

        ccf = flags['ccf'][0]
        ex = flags['ex'][0]
        var = flags['var'][0]
        qph = flags['qph'][0]

        mags_points_written = 0
        with open(mags_outpath, 'a') as photfile:
            for i in range(len(mags.colnames)):
                mag = mags[mags.colnames[i]][0]
                err = errors[errors.colnames[i]][0]

                if err != '--' and ccf[i] == '0' and qph[i] in ['A', 'B']:
                    mags_points_written += 1
                    photfile.write('{} {} {} ccf:{}_ex:{}_var:{}_qph:{}\n'.format(
                        mags.colnames[i][:2], mag, err, ccf, ex, var, qph
                    ))
                    if var[i] not in ['n', '0', '1', '2', '3', '4', '5', '6', '7']:
                        print('[WISE][WARN] {} at {} flagged as variable - {}'.format(
                            star, mags.colnames[i][:2], var[i]
                        ))
                    if ex != 0:
                        print('[WISE][WARN] {} flagged as extended in one or more bands - {}'.format(star, ex))

        print('[WISE][INFO] magnitude points: {} -> {}'.format(mags_points_written, mags_outpath))

        lbd, flux, err, band_flags = convert_wise(star)
        lbd = np.array(lbd)
        flux = np.array(flux)
        err = np.array(err)
        flux = flux * 3e-9 / lbd**2

        flux_corr = color_corr(lbd, flux, 'wise')
        flux_corr = flux_corr * lbd**2 / 3e-9

        cc_outpath = '{0}/{0}_vizier_wise_jy_cc.dat'.format(star)
        if os.path.isfile(cc_outpath):
            os.remove(cc_outpath)

        with open(cc_outpath, 'a') as photfile:
            for i in range(lbd.size):
                photfile.write('{} {:.3f} {:.3f} {}   {}\n'.format(
                    lbd[i], flux_corr[i], err[i], 'wise', band_flags[i]
                ))
        print('[WISE][OK] {} point(s) -> {}'.format(lbd.size, cc_outpath))
    else:
        print('[WISE][INFO] no WISE points')

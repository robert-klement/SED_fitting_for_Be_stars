from astroquery.vizier import Vizier
import numpy as np
import os
import re

from auxiliary_SED_color_correction_IR_missions import color_corr


# starlist = ('HD44637', '10CMa', '20Vul', '25Vul', '25Peg', '120Tau', 'epsPsA', 'HR2249',\
#             'omeOri', 'upsCyg', 'zetCrv')
starlist = ('betCMi',)


print('[AKARI] Processing {} stars'.format(len(starlist)))


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
    print('[AKARI][STAR] {}'.format(star))
    resolved_alias = None

    if not os.path.isdir('{}'.format(star)):
        os.makedirs('{}'.format(star))

    # AKARI IRC
    v = Vizier(columns=['S09', 'S18'])
    fluxes, used_name = query_vizier_with_aliases(
        v, 'query_region', star, radius='0d0m10s', catalog=['II/297/irc']
    )
    if used_name != star and resolved_alias is None:
        resolved_alias = used_name
        print('[AKARI][INFO] {} resolved as {}'.format(star, resolved_alias))
    if len(fluxes) == 0:
        fluxes = []
    else:
        fluxes = fluxes[0]

    v = Vizier(columns=['e_S09', 'e_S18'])
    errors, used_name = query_vizier_with_aliases(
        v, 'query_region', star, radius='0d0m10s', catalog=['II/297/irc']
    )
    if used_name != star and resolved_alias is None:
        resolved_alias = used_name
        print('[AKARI][INFO] {} resolved as {}'.format(star, resolved_alias))
    if len(errors) == 0:
        errors = []
    else:
        errors = errors[0]

    v = Vizier(columns=['X09', 'X18'])
    flags, used_name = query_vizier_with_aliases(
        v, 'query_region', star, radius='0d0m10s', catalog=['II/297/irc']
    )
    if used_name != star and resolved_alias is None:
        resolved_alias = used_name
        print('[AKARI][INFO] {} resolved as {}'.format(star, resolved_alias))
    if len(flags) == 0:
        flags = []
    else:
        flags = flags[0]

    # AKARI FIS
    v = Vizier(columns=['S65', 'S90', 'S140', 'S160'])
    fluxes_fis, used_name = query_vizier_with_aliases(
        v, 'query_region', star, radius='0d0m20s', catalog=['II/298/fis']
    )
    if used_name != star and resolved_alias is None:
        resolved_alias = used_name
        print('[AKARI][INFO] {} resolved as {}'.format(star, resolved_alias))
    if len(fluxes_fis) == 0:
        fluxes_fis = []
    else:
        fluxes_fis = fluxes_fis[0]

    v = Vizier(columns=['e_S65', 'e_S90', 'e_S140', 'e_S160'])
    errors_fis, used_name = query_vizier_with_aliases(
        v, 'query_region', star, radius='0d0m20s', catalog=['II/298/fis']
    )
    if used_name != star and resolved_alias is None:
        resolved_alias = used_name
        print('[AKARI][INFO] {} resolved as {}'.format(star, resolved_alias))
    if len(errors_fis) == 0:
        errors_fis = []
    else:
        errors_fis = errors_fis[0]

    v = Vizier(columns=['q_S65', 'q_S90', 'f_S90', 'q_S140', 'q_S160'])
    flags_fis, used_name = query_vizier_with_aliases(v, 'query_object', star, catalog=['II/298/fis'])
    if used_name != star and resolved_alias is None:
        resolved_alias = used_name
        print('[AKARI][INFO] {} resolved as {}'.format(star, resolved_alias))
    if len(flags_fis) == 0:
        flags_fis = []
    else:
        flags_fis = flags_fis[0]

    if len(errors) != 0:
        raw_outpath = '{0}/{0}_vizier_akari_jy.dat'.format(star)
        if os.path.isfile(raw_outpath):
            os.remove(raw_outpath)

        irc_points_written = 0
        with open(raw_outpath, 'a') as photfile:
            for i in range(len(fluxes.colnames)):
                flux_irc = fluxes[fluxes.colnames[i]][0]
                err_irc = errors[errors.colnames[i]][0]
                x_flags = [flags['X09'][0], flags['X18'][0]]
                if flux_irc != '--' and x_flags[i] == 0:
                    irc_points_written += 1
                    photfile.write('{} {} {} extflags: {}\n'.format(
                        fluxes.colnames[i][-2:], flux_irc, err_irc, x_flags[i]
                    ))

        fis_points_written = 0
        if len(errors_fis) != 0:
            with open(raw_outpath, 'a') as photfile:
                for i in range(len(fluxes_fis.colnames)):
                    flux_fis = fluxes_fis[fluxes_fis.colnames[i]][0]
                    err_fis = errors_fis[errors_fis.colnames[i]][0]
                    qualities = [
                        flags_fis['q_S65'][0],
                        flags_fis['q_S90'][0],
                        flags_fis['q_S140'][0],
                        flags_fis['q_S160'][0],
                    ]
                    if flux_fis != '--' and qualities[i] == 3:
                        fis_points_written += 1
                        photfile.write('{} {:.3f} {:.3f} qualities: {}\n'.format(
                            fluxes_fis.colnames[i][-2:], flux_fis, err_fis, qualities[i]
                        ))
                        if flags_fis['f_S90'][0] != '0':
                            print('[AKARI][WARN] non-0 f_S90 flag: {}'.format(flags_fis['f_S90'][0]))

        print('[AKARI][INFO] raw points: IRC={} FIS={} -> {}'.format(
            irc_points_written, fis_points_written, raw_outpath
        ))

        lbd = np.loadtxt(raw_outpath, usecols=(0))
        flux = np.loadtxt(raw_outpath, usecols=(1))
        err = np.loadtxt(raw_outpath, usecols=(2))
        out_flags = np.loadtxt(raw_outpath, usecols=(4))
        flux = flux * 3e-9 / lbd**2

        cc_outpath = '{0}/{0}_vizier_akari_jy_cc.dat'.format(star)
        if os.path.isfile(cc_outpath):
            os.remove(cc_outpath)

        if not lbd.any():
            print('[AKARI][WARN] extended flag on all AKARI points, removing {}'.format(raw_outpath))
            os.remove(raw_outpath)
        else:
            flux_corr = color_corr(lbd, flux, 'akari')
            flux_corr = flux_corr * lbd**2 / 3e-9
            with open(cc_outpath, 'a') as photfile:
                if lbd.size == 1:
                    photfile.write('{} {:.3f} {:.3f} akari   flags: {:.0f}\n'.format(
                        lbd, flux_corr[0], err, out_flags
                    ))
                else:
                    for i in range(lbd.size):
                        photfile.write('{} {:.3f} {:.3f} akari   flags: {:.0f}\n'.format(
                            lbd[i], flux_corr[i], err[i], out_flags[i]
                        ))
            print('[AKARI][OK] {} point(s) -> {}'.format(lbd.size, cc_outpath))
    else:
        print('[AKARI][INFO] no AKARI points')

from astroquery.vizier import Vizier
import numpy as np
import os
import re

from auxiliary_SED_color_correction_IR_missions import color_corr


starlist = ('HD44637', '10CMa', '20Vul', '25Vul', '25Peg', '120Tau', 'epsPsA', 'HR2249',\
            'omeOri', 'upsCyg', 'zetCrv')
# starlist = ('betCMi',)


print('[IRAS] Processing {} stars'.format(len(starlist)))


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
    print('[IRAS][STAR] {}'.format(star))
    resolved_alias = None

    if not os.path.isdir('{}'.format(star)):
        os.makedirs('{}'.format(star))

    # IRAS
    v = Vizier(columns=['Fnu_12', 'Fnu_25', 'Fnu_60'])
    fluxes, used_name = query_vizier_with_aliases(v, 'query_object', star, catalog=['II/125/main'])
    if used_name != star and resolved_alias is None:
        resolved_alias = used_name
        print('[IRAS][INFO] {} resolved as {}'.format(star, resolved_alias))
    if len(fluxes) == 0:
        fluxes = []
    else:
        fluxes = fluxes[0]

    v = Vizier(columns=['e_Fnu_12', 'e_Fnu_25', 'e_Fnu_60'])
    errors, used_name = query_vizier_with_aliases(v, 'query_object', star, catalog=['II/125/main'])
    if used_name != star and resolved_alias is None:
        resolved_alias = used_name
        print('[IRAS][INFO] {} resolved as {}'.format(star, resolved_alias))
    if len(errors) == 0:
        errors = []
    else:
        errors = errors[0]

    v = Vizier(columns=['q_Fnu_12', 'q_Fnu_25', 'q_Fnu_60', 'TSNR_12', 'TSNR_25', 'TSNR_60',
                        'CC_12', 'CC_25', 'CC_60', 'Confuse', 'SES1_12', 'SES1_25', 'SES1_60',
                        'SES2_12', 'SES2_25', 'SES2_60'])
    flags, used_name = query_vizier_with_aliases(v, 'query_object', star, catalog=['II/125/main'])
    if used_name != star and resolved_alias is None:
        resolved_alias = used_name
        print('[IRAS][INFO] {} resolved as {}'.format(star, resolved_alias))

    v = Vizier(columns=['HSDFlag'])
    hsdflag, used_name = query_vizier_with_aliases(v, 'query_object', star, catalog=['II/125/main'])
    if used_name != star and resolved_alias is None:
        resolved_alias = used_name
        print('[IRAS][INFO] {} resolved as {}'.format(star, resolved_alias))

    if len(flags) == 0:
        flags = []
        hsdflag = []
    else:
        flags = flags[0]
        hsdflag = hsdflag[0]
        hsdflag = hsdflag['HSDFlag'][0]

        q_fnu_12 = flags['q_Fnu_12'][0]
        q_fnu_25 = flags['q_Fnu_25'][0]
        q_fnu_60 = flags['q_Fnu_60'][0]

        tsnr_12 = flags['TSNR_12'][0]
        tsnr_25 = flags['TSNR_25'][0]
        tsnr_60 = flags['TSNR_60'][0]

        cc_12 = flags['CC_12'][0]
        cc_25 = flags['CC_25'][0]
        if cc_25 == '':
            cc_25 = 'X'
        cc_60 = flags['CC_60'][0]
        if cc_60 == '':
            cc_60 = 'X'

        confuse = flags['Confuse'][0]

        ses1_12 = flags['SES1_12'][0]
        ses1_25 = flags['SES1_25'][0]
        ses1_60 = flags['SES1_60'][0]

        ses2_12 = flags['SES2_12'][0]
        ses2_25 = flags['SES2_25'][0]
        ses2_60 = flags['SES2_60'][0]

    raw_outpath = '{0}/{0}_vizier_iras_jy.dat'.format(star)
    if len(fluxes) != 0:
        if os.path.isfile(raw_outpath):
            os.remove(raw_outpath)

        points_written = 0
        with open(raw_outpath, 'a') as photfile:
            for i in range(len(fluxes.colnames)):
                flux = fluxes[fluxes.colnames[i]][0]
                err = errors[errors.colnames[i]][0] * 0.01 * flux
                if err != '--':
                    points_written += 1
                    photfile.write(
                        '{} {:.5f} {:.5f} flags: q:{}{}{}_snr:{:00005}_{:00005}_{:00005}_'
                        'cc:{:1}_{:1}_{:1}_ses1:{:1}_{:1}_{:1}_ses2:{:1}_{:1}_{:1}_conf:{}_hsd:{}\n'.format(
                            fluxes.colnames[i][-2:], flux, err,
                            q_fnu_12, q_fnu_25, q_fnu_60,
                            tsnr_12, tsnr_25, tsnr_60,
                            cc_12, cc_25, cc_60,
                            ses1_12, ses1_25, ses1_60,
                            ses2_12, ses2_25, ses2_60,
                            confuse, hsdflag,
                        )
                    )

        print('[IRAS][OK] {} point(s) -> {}'.format(points_written, raw_outpath))
    else:
        print('[IRAS][INFO] no IRAS points')

    if os.path.isfile(raw_outpath):
        lbd = []
        flux = []
        err = []
        flags_out = []

        with open(raw_outpath, 'r') as f:
            for line in f:
                line = line.strip()
                columns = line.split()

                flags = columns[4]
                qualities = (flags[2:5])
                snr_12 = float(flags[10:15])
                snr_25 = float(flags[16:21])
                snr_60 = float(flags[22:27])
                cc_12, cc_25, cc_60 = flags[31], flags[33], flags[35]
                cc_allowed = ['A', 'B', 'C', 'D']
                ses1_12, ses1_25, ses1_60 = flags[42], flags[44], flags[46]
                ses2_12, ses2_25, ses2_60 = flags[53], flags[55], flags[57]
                conf = flags[64]
                hsd = flags[70]

                if line[0:2] == '12':
                    if (qualities[0] != '1' and cc_12 in cc_allowed and ses1_12 in ['0', '1'] and ses2_12 == '0' and snr_12 > 50 and conf not in ['1', '3', '5', '7', '9', 'B', 'D', 'F']):
                        lbd.append(float(columns[0]))
                        flux.append(float(columns[1]))
                        err.append(float(columns[2]))
                        flags_out.append(flags)
                        if snr_12 < 60:
                            print('{} WARNING: snr12 = {}'.format(star, snr_12 / 10))
                        if conf in ['1', '3', '5', '7', '9', 'B', 'D', 'F']:
                            print('{} WARNING: 12um point confused - {}'.format(star, conf))
                        if hsd in ['1', '3', '5', '7', '9', 'B', 'D', 'F']:
                            print('{} WARNING: 12um in HSD region'.format(star))

                if (line[0:2]) == '25':
                    if (qualities[1] != '1' and cc_25 in cc_allowed and ses1_25 in ['0', '1'] and ses2_25 == '0' and snr_25 > 50 and conf not in ['2', '3', '6', '7', 'A', 'B', 'F']):
                        lbd.append(float(columns[0]))
                        flux.append(float(columns[1]))
                        err.append(float(columns[2]))
                        flags_out.append(flags)
                        if snr_25 < 60:
                            print('{} WARNING: snr25 = {}'.format(star, snr_25 / 10))
                        if conf in ['2', '3', '6', '7', 'A', 'B', 'F']:
                            print('{} WARNING: 25um point confused - {}'.format(star, conf))
                        if hsd in ['2', '3', '6', '7', 'A', 'B', 'F']:
                            print('{} WARNING: 25um in HSD region'.format(star))

                if (line[0:2]) == '60':
                    if (qualities[2] != '1' and cc_60 in cc_allowed and ses1_60 in ['0', '1'] and ses2_60 == '0' and snr_60 > 50 and conf not in ['4', '5', '6', '7', 'D', 'F']):
                        lbd.append(float(columns[0]))
                        flux.append(float(columns[1]))
                        err.append(float(columns[2]))
                        flags_out.append(flags)
                        if snr_60 < 60:
                            print('{} WARNING: snr60 = {}'.format(star, snr_60 / 10))
                        if conf in ['4', '5', '6', '7', 'D', 'F']:
                            print('{} WARNING: 60um point confused - {}'.format(star, conf))
                        if hsd in ['4', '5', '6', '7', 'D', 'F']:
                            print('{} WARNING: 60um point in HSD region'.format(star))

        lbd = np.array(lbd)
        flux = np.array(flux)
        err = np.array(err)

        if lbd.size == 0:
            print('{} INFO: no IRAS points passed quality filters'.format(star))
            continue

        flux = flux * 3e-9 / lbd**2
        flux_corr = color_corr(lbd, flux, 'iras')
        flux_corr = flux_corr * lbd**2 / 3e-9

        cc_outpath = '{0}/{0}_vizier_iras_jy_cc.dat'.format(star)
        if os.path.isfile(cc_outpath):
            os.remove(cc_outpath)

        with open(cc_outpath, 'a') as photfile:
            if lbd.size == 1:
                photfile.write('{} {:.2f} {:.2f} iras {}\n'.format(lbd[0], flux_corr[0], err[0], flags_out[0]))
            else:
                for i in range(lbd.size):
                    photfile.write('{} {:.2f} {:.2f} iras {}\n'.format(lbd[i], flux_corr[i], err[i], flags_out[i]))
    else:
        print('{} INFO: missing {}/{}_vizier_iras_jy.dat'.format(star, star, star))

import os


# USER CONFIG
starlist = ('HD44637', '10CMa', '20Vul', '25Vul', '25Peg', '120Tau', 'epsPsA', 'HR2249',\
            'omeOri', 'upsCyg', 'zetCrv')
# starlist = ('betCMi',)


def main():
    print('[MERGE] Processing {} stars'.format(len(starlist)))

    for star in starlist:
        print('[MERGE][STAR] {}'.format(star))
        star_dir = '{}'.format(star)
        outpath = '{0}/{0}_photometry.dat'.format(star)

        if not os.path.isdir(star_dir):
            os.makedirs(star_dir)

        if os.path.isfile(outpath):
            os.remove(outpath)

        lines_written = 0

        # SIMBAD
        simbad_path = '{0}/{0}_simbad_johnson_jy.dat'.format(star)
        if os.path.isfile(simbad_path):
            with open(simbad_path, 'r') as simbad:
                simbad_text = simbad.read()
                with open(outpath, 'w') as photfile:
                    photfile.write(simbad_text)
            lines_written += len([line for line in simbad_text.splitlines() if line.strip()])
        else:
            print('[MERGE][INFO] {} - no SIMBAD file ({})'.format(star, simbad_path))
            # Create an empty output file so remaining catalogs can still be appended.
            open(outpath, 'w').close()

        # DOUGHERTY
        dougherty_path = '{0}/{0}_dougherty_johnson_jy.dat'.format(star)
        if os.path.isfile(dougherty_path):
            with open(dougherty_path, 'r') as dougherty:
                dougherty_text = dougherty.read()
                with open(outpath, 'a') as photfile:
                    photfile.write(dougherty_text)
            lines_written += len([line for line in dougherty_text.splitlines() if line.strip()])
        else:
            print('[MERGE][INFO] {} - no Dougherty'.format(star))

        # IRAS
        iras_path = '{0}/{0}_vizier_iras_jy_cc.dat'.format(star)
        if os.path.isfile(iras_path):
            with open(iras_path, 'r') as iras:
                iras_text = iras.read()
                with open(outpath, 'a') as photfile:
                    photfile.write(iras_text)
            lines_written += len([line for line in iras_text.splitlines() if line.strip()])
        else:
            print('[MERGE][INFO] {} - no IRAS'.format(star))

        # AKARI
        akari_path = '{0}/{0}_vizier_akari_jy_cc.dat'.format(star)
        if os.path.isfile(akari_path):
            with open(akari_path, 'r') as akari:
                akari_text = akari.read()
                with open(outpath, 'a') as photfile:
                    photfile.write(akari_text)
            lines_written += len([line for line in akari_text.splitlines() if line.strip()])
        else:
            print('[MERGE][INFO] {} - no AKARI'.format(star))

        # WISE
        wise_path = '{0}/{0}_vizier_wise_jy_cc.dat'.format(star)
        if os.path.isfile(wise_path):
            with open(wise_path, 'r') as wise:
                wise_text = wise.read()
                if wise_text:
                    with open(outpath, 'a') as photfile:
                        photfile.write(wise_text)
                    lines_written += len([line for line in wise_text.splitlines() if line.strip()])
                else:
                    print('[MERGE][INFO] {} - empty WISE photfile'.format(star))
        else:
            print('[MERGE][INFO] {} - no WISE'.format(star))

        # RADIO
        radio_path = '{0}/{0}_radio.dat'.format(star)
        if os.path.isfile(radio_path):
            with open(radio_path, 'r') as radio:
                radio_text = radio.read()
                with open(outpath, 'a') as photfile:
                    photfile.write(radio_text)
            lines_written += len([line for line in radio_text.splitlines() if line.strip()])
        else:
            print('[MERGE][INFO] {} - no radio'.format(star))

        print('[MERGE][OK] {} line(s) -> {}'.format(lines_written, outpath))


if __name__ == '__main__':
    main()

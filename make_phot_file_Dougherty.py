# IMPORTANT:
# This script requires a manually created <star>/<star>_dougherty_johnson.dat file.
# The input photometry should be taken from Dougherty et al. (1991), their Table 2.

import os

from auxiliary_SED_convert_mags import convert_johnson_dougherty

# starlist = ('HD44637', '10CMa', '20Vul', '25Vul', '25Peg', '120Tau', 'epsPsA', 'HR2249',\
#             'omeOri', 'upsCyg', 'zetCrv')
starlist = ('betCMi',)


print('[DOUGHERTY] Processing {} stars'.format(len(starlist)))

converted = 0
skipped = 0

for star in starlist:
    print('[DOUGHERTY][STAR] {}'.format(star))

    inpath = '{0}/{0}_dougherty_johnson.dat'.format(star)
    if not os.path.isfile(inpath):
        print('[DOUGHERTY][INFO] missing {}, skipping'.format(inpath))
        skipped += 1
        continue

    convert_johnson_dougherty(star)
    converted += 1
    print('[DOUGHERTY][OK] Wrote {}/{}_dougherty_johnson_jy.dat'.format(star, star))

print('[DOUGHERTY][DONE] converted={} skipped={}'.format(converted, skipped))

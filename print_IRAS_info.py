from astroquery.vizier import Vizier
import os

# starlist = ('HD44637', '10CMa', '20Vul', '25Vul', '25Peg', '120Tau', 'epsPsA', 'HR2249',\
#             'omeOri', 'upsCyg', 'zetCrv')
starlist = ('omeOri',)

for star in starlist:
    v = Vizier(columns=['Fnu_12', 'q_Fnu_12', 'TSNR_12', 'CC_12', 'Confuse'])

    info = v.query_region(star, radius="0d0m20s", catalog=['II/125/main'])

    if len(info) == 0:
        info = []

    if len(info) != 0:
        with open('iras_info.dat', 'a') as f:
            print('{}: {}\n'.format(star, info[0]))
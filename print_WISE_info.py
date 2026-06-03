from astroquery.vizier import Vizier
import os

# starlist = ('HD44637', '10CMa', '20Vul', '25Vul', '25Peg', '120Tau', 'epsPsA', 'HR2249',\
#             'omeOri', 'upsCyg', 'zetCrv')
starlist = ('omeOri',)

exists = os.path.isfile('wise_info.dat')
if exists:
    os.remove('wise_info.dat')

for star in starlist:
    v = Vizier(columns=['ccf', 'ex', 'var', 'qph'])

    info = v.query_region(star, radius="0d0m5s", catalog=['II/328/allwise'])

    if len(info) == 0:
        info = []

    if len(info) != 0:
        # with open('wise_info.dat', 'a') as f:
        print('{}: {}\n'.format(star, info[0]))
            # f.write('{}:\n {}\n\n'.format(star, info[0]))

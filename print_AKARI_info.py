from astroquery.vizier import Vizier
import os

# starlist = ('HD44637', '10CMa', '20Vul', '25Vul', '25Peg', '120Tau', 'epsPsA', 'HR2249',\
#             'omeOri', 'upsCyg', 'zetCrv')
starlist = ('omeOri',)

exists = os.path.isfile('akari_info.dat')
if exists:
    os.remove('akari_info.dat')

for star in starlist:
    v = Vizier(columns=['q_S09', 'Nd09', 'X09', 'q_S18', 'Nd18', 'X18'])

    info = v.query_region(star, radius="0d0m10s", catalog=['II/297/irc'])

    if len(info) == 0:
        info = []

    if len(info) != 0:
        print('{}: {}\n'.format(star, info[0]))
    # with open('akari_info.dat', 'a') as f:
    #     print('{}: {}\n'.format(star, info[0]))
    #     f.write('{}: {}\n'.format(star, info[0]))

for star in starlist:
    v = Vizier(columns=['S65', 'e_S65', 'S90', 'e_S90', 'S140', 'e_S140', 'S160', 'e_S160'])
    # v = Vizier(columns=['q_S65'])#, 'q_S90', 'q_S140','q_S160'])

    info = v.query_object(star, catalog=['II/298/fis'])

    if len(info) == 0:
        info = []

    if len(info) != 0:
        print('{}: {}\n'.format(star, info[0]))
    # with open('akari_info.dat', 'a') as f:
    #     print('{}: {}\n'.format(star, info[0]))
    #     f.write('{}: {}\n'.format(star, info[0]))
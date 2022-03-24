# script to anlayze/plot outputs from recycling scan

import numpy as np
import neutrals_analysis as neut
import matplotlib.pyplot as plt

SOLPSWORK = '/nobackup1/millerma/solps-iter/runs'
shot = 1070614013

scan_abs = np.linspace(0,9,10)
scan_pressure = []
scan_nnsep = []

for scan in scan_abs:

	experiment = 'pump{}'.format(int(scan))
	
	try:	
		# get cryo pressure
		p_cryo = neut.calc_neutral_pressure(SOLPSWORK,shot,experiment,1,'F_CRYO')
		scan_pressure.append(p_cryo)

		# get nn_sep
		nn_sep = neut.get_neutral_density(SOLPSWORK,shot,experiment,1,'sep')
		scan_nnsep.append(nn_sep)

	except Exception as e:
		print(e)
		#print('run didnt finish')

scan_R = 1 - np.array(scan_abs)*.1

fig, ax = plt.subplots(2)
ax[0].plot(scan_R, scan_pressure, 'o')
ax[1].plot(scan_R, scan_nnsep, 'o')

ax[0].set_ylabel('Baffle pressure (mTorr)')
ax[1].set_ylabel('$n_{D,sep} (cm^{-3})$')

ax[0].set_xlabel('Recycling coefficient')
plt.show()

					
	


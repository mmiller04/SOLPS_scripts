# script to anlayze/plot outputs from recycling scan

import numpy as np
import neutrals_analysis as neut
import matplotlib.pyplot as plt
import sys

SOLPSWORK = '/nobackup1/millerma/solps-iter/runs'
shot = 1070614013
attempt = sys.argv[1]

scan_abs = np.linspace(0,9,10)
scan_pF = []
scan_pG = []
scan_nnsep = []
scan_pflux = []
scan_mflux = []

for scan in scan_abs:

	experiment = 'pump{}'.format(int(scan))
	
	#try:	
		# get cryo/omp pressure
	p_cryo = neut.calc_neutral_pressure(SOLPSWORK,shot,experiment,attempt,'F_CRYO')
	p_omp = neut.calc_neutral_pressure(SOLPSWORK,shot,experiment,attempt,'G-SIDE_RAT')

	scan_pF.append(p_cryo)
	scan_pG.append(p_omp)

	# get nn_sep
	nn_sep = neut.get_neutral_density(SOLPSWORK,shot,experiment,attempt,'sep')
	scan_nnsep.append(nn_sep)

	# get pumped flux
	pumped_flux = neut.calc_pumped_flux(SOLPSWORK,shot,experiment,attempt)
	scan_pflux.append(pumped_flux)
		
		# get mom flux
	mom_flux = neut.calc_mom_flux(SOLPSWORK,shot,experiment,attempt)
	scan_mflux.append(mom_flux)

	#except Exception as e:
	#	print(e)
		#print('run didnt finish')

scan_R = 1 - np.array(scan_abs)*.1

figa, axa = plt.subplots(3)
axa[0].plot(scan_R, scan_pF, '-o')
axa[1].plot(scan_R, scan_pflux, '-o')
axa[2].plot(scan_R, scan_mflux, '-o')

axa[0].set_ylabel('$P_{EIRENE} Baffle (mTorr)$')
axa[1].set_ylabel('pumped_flux $(s^{-1})$')
axa[2].set_ylabel('momentum flux density $(kgm^{-2}s^{-1})$')

axa[2].set_xlabel('Recycling coefficient')

figb, axb = plt.subplots(2)
axb[0].plot(scan_R, scan_pG, '-o')
axb[1].plot(scan_R, scan_nnsep, '-o')

axb[0].set_ylabel('$P_{EIRENE} OMP (mTorr)$')
axb[1].set_ylabel('$n_{D}(R_{max}) (cm^{-3})$')

axb[1].set_xlabel('Recycling coefficient')
plt.show()

					
	


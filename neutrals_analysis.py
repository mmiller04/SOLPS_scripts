# set of routines to get neutral info from SOLPS output and analyze it

import numpy as np
import matplotlib.pyplot as plt

import aurora

from scipy.interpolate import interp1d

def calc_neutral_pressure(SOLPSWORK, shot, experiment, attempt, loc, plot=False):

	# initialize aurora solps object
	b2path = '{}/{}/{}/attempt{}'.format(SOLPSWORK, shot, experiment, attempt)

	so = aurora.solps_case(
		b2fstate_path='{}/b2fstate'.format(b2path),
		b2fgmtry_path='{}/b2fgmtry'.format(b2path),
		geqdsk='/nobackup1/millerma/solps-iter/runs/1070614013/pump_on/baserun/g1070614013.00793_610'
	)

	# pull energy density from fort46 file
	edena = so.fort46['edena']

	if loc == 'F_CRYO':
	
		inds = np.where(so.ynodes >= so.ynodes.max()-1e-3)[0] # get indices that correspond to vertices lying at the top of VV
		tris = so.triangles.astype(int) # triangles indices saved as floats for some reason

		xcent, pcent = [], []

		for i in inds:
			
			num = np.where(tris == i)[0] 
			corner = np.where(tris == i)[0]

			# most vertices are shared so need to iterate through shared vertex coords
			for j in range(len(num)):
			
				xcent.append(so.xnodes[tris[num[j]]].mean())
				pcent.append(edena[num[j]])

		# cryopump baffle coords
		R_cryo = [0.574, 0.641] # these correspond to SOLPS grid but may not correspond to real space
		
		# conver to arrays for masking
		xcent, pcent = np.array(xcent), np.array(pcent)

		mask = np.logical_and(xcent > R_cryo[0], xcent < R_cryo[1])

		avg_edena_pump = np.mean(pcent[mask])*1e6 # convert to eV/m^3
		avg_edena_pump_mtorr = avg_edena_pump*1.602e-19/133*1e3 # convert to Pa --> Torr --> mTorr
		pressure_out = 2/3*avg_edena_pump_mtorr # pressure in mTorr


		if plot:

			fig, ax = plt.subplots()

			ax.plot(xcent[mask], pcent[mask], 'x', c='r')
			ax.plot(xcent, pcent, 'o')
			ax.legend(['Pump baffle'])

			ax.set_xlabel('Radius (m)')
			ax.set_ylabel('Pressure (mTorr)')


	return pressure_out
	

def get_neutral_density(SOLPSWORK, shot, experiment, attempt, loc='sep', plot=False):

	# initialize aurora solps object
	b2path = '{}/{}/{}/attempt{}'.format(SOLPSWORK, shot, experiment, attempt)

	so = aurora.solps_case(
		b2fstate_path='{}/b2fstate'.format(b2path),
		b2fgmtry_path='{}/b2fgmtry'.format(b2path),
		geqdsk='/nobackup1/millerma/solps-iter/runs/1070614013/pump_on/baserun/g1070614013.00793_610'
	)

	# pull neutral atomic density from fort44 file
	dab2 = so.fort44['dab2']
	
	# get radial coordinates in flux space
	R_LFS = so.data('cr')[:,38]
	Z_LFS = so.data('cz')[:,38]

	rhop_LFS = aurora.coords.get_rhop_RZ(R_LFS,Z_LFS,so.geqdsk)

	if loc == 'sep':
		psin = 1

	f = interp1d(rhop_LFS, dab2[38,:,0])
	#nn_out = f(psin)
	nn_out = dab2[38,-1,0]

	return nn_out


def get_neutral_density(SOLPSWORK, shot, experiment, attempt, loc='sep', plot=False):

	


	




# balance.nc gives source
# fort.44/fort.46 gives n_dens, energy_dens (pressure) on B2/EIRENE grids

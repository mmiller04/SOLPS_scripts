# set of routines to get neutral info from SOLPS output and analyze it

import numpy as np
import matplotlib.pyplot as plt

import aurora
import netCDF4 as nc
import glob

from scipy.interpolate import interp1d


def calc_neutral_pressure(SOLPSWORK, shot, experiment, name, loc, plot=False):

	# initialize aurora solps object
	b2path = '{}/{}/{}/{}'.format(SOLPSWORK, shot, experiment, name)

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

	
	if loc == 'G-SIDE_RAT':

		inds = np.where(so.xnodes >= so.xnodes.max() - 1e-3)[0] # get indices that correspond to vertices lying at the right of VV
		tris = so.triangles.astype(int) # triangles indices saved as floats for some reason

		xcent, pcent = [], []

		for i in inds:

			num = np.where(tris == i)[0]
			corner = np.where(tris == i)[0]

			# most vertices are shared so need to iterate through shared vertex coords
			for j in range(len(num)):
		
				xcent.append(so.xnodes[tris[num[j]]].mean())
				pcent.append(edena[num[j]])

	
		avg_edena_omp = np.mean(pcent)*1e6 # convert to eV/m^3
		avg_edena_omp_mtorr = avg_edena_omp*1.602e-19/133*1e3 # convert to Pa --> Torr --> mTorr
		pressure_out = 2/3*avg_edena_omp_mtorr # pressure in mTorr

	return pressure_out


def calc_mom_flux(SOLPSWORK, shot, experiment, name, plot=False):

	# initialize aurora solps object
	b2path = '{}/{}/{}/{}'.format(SOLPSWORK, shot, experiment, name)

	so = aurora.solps_case(
		b2fstate_path='{}/b2fstate'.format(b2path),
		b2fgmtry_path='{}/b2fgmtry'.format(b2path),
		geqdsk='/nobackup1/millerma/solps-iter/runs/1070614013/pump_on/baserun/g1070614013.00793_610'
	)

	# pull energy density from fort46 file
	vydena= so.fort46['vydena']
	
	inds = np.where(so.ynodes >= so.ynodes.max()-1e-3)[0] # get indices that correspond to vertices lying at the top of VV
	tris = so.triangles.astype(int) # triangles indices saved as floats for some reason

	xcent, mcent = [], []

	for i in inds:
			
		num = np.where(tris == i)[0] 
		corner = np.where(tris == i)[0]

		# most vertices are shared so need to iterate through shared vertex coords
		for j in range(len(num)):
			
			xcent.append(so.xnodes[tris[num[j]]].mean())
			mcent.append(vydena[num[j]])

	# cryopump baffle coords
	R_cryo = [0.574, 0.641] # these correspond to SOLPS grid but may not correspond to real space
		
	# conver to arrays for masking
	xcent, mcent = np.array(xcent), np.array(mcent)

	mask = np.logical_and(xcent > R_cryo[0], xcent < R_cryo[1])

	avg_mflux_up = np.mean(mcent[mask])/1e3*1e4 # convert to kg/m^2/s^-1


	if plot:

		fig, ax = plt.subplots()

		ax.plot(xcent[mask], mcent[mask], 'x', c='r')
		ax.plot(xcent, mcent, 'o')
		ax.legend(['Pump baffle'])

		ax.set_xlabel('Radius (m)')
		ax.set_ylabel('Pressure (mTorr)')

	return avg_mflux_up
	

def calc_pumped_flux(SOLPSWORK, shot, experiment, name, plot=False):
	
	# initialize aurora solps object
	b2path = '{}/{}/{}/{}'.format(SOLPSWORK, shot, experiment, name)

	so = aurora.solps_case(
		b2fstate_path='{}/b2fstate'.format(b2path),
		b2fgmtry_path='{}/b2fgmtry'.format(b2path),
		geqdsk='/nobackup1/millerma/solps-iter/runs/1070614013/pump_on/baserun/g1070614013.00793_610'
	)

	pump_inds = [82, 83, 96, 97]

	pumped_flux = np.mean(so.fort44['wlpump(A)'][0])

	return pumped_flux


def get_neutral_density(SOLPSWORK, shot, experiment, name, loc='sep', plot=False):

	# initialize aurora solps object
	b2path = '{}/{}/{}/{}'.format(SOLPSWORK, shot, experiment, name)

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


#def get_neutral_density(SOLPSWORK, shot, experiment, attempt, loc='sep', plot=False):

def undo_dumb_stuff(dumb_stuff,dumb_helper):

	psin = np.array(dumb_stuff)+1
	rhop = np.sqrt(psin)

	Rmid = aurora.rad_coord_transform(rhop,'rhop','Rmid',dumb_helper)
	Rsep = aurora.rad_coord_transform(1,'rhop','Rmid',dumb_helper)

	undumb_stuff = Rmid-Rsep

	return undumb_stuff


def load_aurora_socase(SOLPSWORK, shot, experiment, name):

	run_path = '{}/{}/{}/{}'.format(SOLPSWORK,shot,experiment,name)
	base_path = '{}/{}/{}/baserun'.format(SOLPSWORK,shot,experiment)

	so = aurora.solps_case(
		b2fstate_path="{}/b2fstate".format(run_path),
		b2fgmtry_path="{}/b2fgmtry".format(run_path),
		geqdsk = glob.glob("{}/g{}.*".format(base_path,shot))[0]
	)

	return so


def get_2dcontour_quantity(run_path, quantity, so):

	toplot = None

	if quantity == 'nn':
		
		toplot = so.fort46['pdena'][:,0]*1e6		

	elif quantity == 'source':

		fn = '{}/balance.nc'.format(run_path)
		ds = nc.Dataset(fn)
	
		sna = ds['eirene_mc_papl_sna_bal']
		vol = ds['vol'] # cell volumes
		crx = ds['crx'] # x-coords
		cry = ds['cry'] # y-coords

		sna_sum = np.sum(sna, axis=0) # sum over EIRENE strata
		sna_Dplus_vol = sna_sum[1]/vol # get source per vol
		
		sna_b2grid = sna_Dplus_vol[1:-1,1:-1] # cut off first and last elements to match B2 grid from aurora
	
		toplot = sna_b2grid

	return toplot

		


# balance.nc gives source
# fort.44/fort.46 gives n_dens, energy_dens (pressure) on B2/EIRENE grids






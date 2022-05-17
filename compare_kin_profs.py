import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import sys
from scipy.interpolate import interp1d

params = ['ne', 'te']
shot_on = 1070614013
shot_off = 1070614016
#experiment = 'pump_off'
#shot = 1120917011
experiment = sys.argv[1]
name = sys.argv[2]

SOLPSWORK = '/nobackup1/users/millerma/solps-iter/runs'

exp_on_file = '{}/exp_data/lyman_data_{}.pkl'.format(SOLPSWORK, shot_on)
exp_off_file = '{}/exp_data/lyman_data_{}.pkl'.format(SOLPSWORK, shot_off)
with open(exp_on_file, 'rb') as fon:
	exp_on_obj = pkl.load(fon)
with open(exp_off_file, 'rb') as foff:
	exp_off_obj = pkl.load(foff)

# put into dictionary
exp = {'on':{}, 'off':{}}
exp['on'], exp['off'] = {'ne':{}, 'te':{}}, {'ne':{}, 'te':{}}

ne_on, ne_off = exp['on']['ne'], exp['off']['ne']
te_on, te_off = exp['on']['te'], exp['off']['te']

#psin_ne, ne, ne_unc, psin_Te, Te, Te_unc = exp_obj
psin_ne_on, psin_ne_off = exp_on_obj[0]**2, exp_off_obj[0]**2
R_ne_on, R_ne_off = exp_on_obj[2], exp_off_obj[2]
ne_data_on, ne_data_off = exp_on_obj[7], exp_off_obj[7]
ne_unc_on, ne_unc_off = exp_on_obj[8], exp_off_obj[8]

psin_Te_on, psin_Te_off = exp_on_obj[0]**2, exp_off_obj[0]**2
R_Te_on, R_Te_off = exp_on_obj[2], exp_off_obj[2]
Te_data_on, Te_data_off = exp_on_obj[9], exp_off_obj[9]
Te_unc_on, Te_unc_off = exp_on_obj[10], exp_off_obj[10]

# find R_sep by interpolating onto psi grid
R_sep_on = interp1d(psin_ne_on, R_ne_on)(1)
R_sep_off = interp1d(psin_ne_off, R_ne_off)(1)

#ne_on['X'], ne_off['X'] = psin_ne_on, psin_ne_off
ne_on['X'], ne_off['X'] = R_ne_on-R_sep_on, R_ne_off-R_sep_off
ne_on['Y'], ne_off['Y'] = ne_data_on*1e6, ne_data_off*1e6 # convert to m^{-3}
ne_on['Y_unc'], ne_off['Y_unc'] = ne_unc_on*1e6, ne_unc_off*1e6 # convert to m^{-3}

#te_on['X'], te_off['X'] = psin_Te_on, psin_Te_off
te_on['X'], te_off['X'] = R_Te_on-R_sep_on, R_Te_off-R_sep_off
te_on['Y'], te_off['Y'] = Te_data_on, Te_data_off
te_on['Y_unc'], te_off['Y_unc'] = Te_unc_on, Te_unc_off # in eV


SOLPS = {'ne':{}, 'te':{}}

b2path = '{}/{}/{}/{}'.format(SOLPSWORK, shot_on, experiment, name)

# use aurora to get ne/te from solps
import aurora
so = aurora.solps_case(
	b2fstate_path='{}/b2fstate'.format(b2path),
	b2fgmtry_path='{}/b2fgmtry'.format(b2path),
	geqdsk='/nobackup1/millerma/solps-iter/runs/1070614013/pump_on/baserun/g1070614013.00793_610'
)

# midplane index - in theory can be read from b2mn.dat but might be switched
jxa = 38
R_sep_SOLPS = aurora.rad_coord_transform(1,'rhop','Rmid',so.geqdsk)

for kw in params:

	SOLPS[kw]['X'] = so.data('cr')[:,jxa] - R_sep_SOLPS # in units of R - R_sep
	SOLPS[kw]['Y'] = so.data(kw)[:,jxa]


tparams = ['dna0','dpa0','hcib','hce0','vla0_x','vla0_y','vsa0','sig0','alf0']
transport = {'dna0':{}, 'hcib':{}, 'hce0':{}}

transport_file = '{}/{}/{}/{}/b2.transport.inputfile'.format(SOLPSWORK, shot_on, experiment, name)
tlist = []
with open(transport_file,'r') as f:
	for line in f.readlines():
		tlist.append(line.split('data'))

transport['dna0']['X'], transport['dna0']['Y'] = [], []
transport['hcib']['X'], transport['hcib']['Y'] = [], []
transport['hce0']['X'], transport['hce0']['Y'] = [], []

ind = 0
for line in tlist:

	if line[0].strip() == 'n':
		ind = line[1][4]

	if line[0].strip() == 't':
		xcoord_line = line[1].split('=')
		ycoord_line = line[2].split('=')
		
		xcoord = xcoord_line[1].split(',')		
		ycoord = ycoord_line[1].split(',')		

		transport[str(tparams[int(ind)-1])]['X'].append(float(xcoord[0].strip()))
		transport[str(tparams[int(ind)-1])]['Y'].append(float(ycoord[0].strip()))

#from neutrals_analysis import undo_dumb_stuff as uds
#import aurora

#geqdsk = '/nobackup1/millerma/solps-iter/runs/1070614013/basecase/baserun/g1070614013.00793_610'
#b2fstate = '/nobackup1/millerma/solps-iter/runs/1070614013/basecase/attempt5/b2fstate'
#b2fgmtry = '/nobackup1/millerma/solps-iter/runs/1070614013/basecase/attempt5/b2fgmtry'
#so = aurora.solps_case(
#	b2fstate_path=b2fstate,
#	b2fgmtry_path=b2fgmtry,
#	geqdsk=geqdsk
#)

#newD = uds(transport['dna0']['X'],so.geqdsk)
#newXi = uds(transport['hcib']['X'],so.geqdsk)
#newXe = uds(transport['hce0']['X'],so.geqdsk)

fig, ax = plt.subplots(2,2, sharex=True)
#ax[0,0].errorbar(exp['on']['ne']['X'], exp['on']['ne']['Y'], exp['on']['ne']['Y_unc'])
#ax[0,0].errorbar(exp['off']['ne']['X'], exp['off']['ne']['Y'], exp['off']['ne']['Y_unc'])
ax[0,0].plot(SOLPS['ne']['X'], SOLPS['ne']['Y'])
ax[0,0].plot(exp['on']['ne']['X'], exp['on']['ne']['Y'])
#ax[0,0].plot(exp['off']['ne']['X'], exp['off']['ne']['Y'],'.')
ax[0,0].legend(['SOLPS','1070614013 (on)','1070614016 (off)'])
ax[0,0].set_ylabel('$n_{e} (m^{-3})$', fontsize=14)
ax[0,0].tick_params(axis='y', labelsize=14)
ax[0,0].tick_params(axis='x', labelsize=14)

ax[1,0].plot(transport['dna0']['X'], transport['dna0']['Y'],'-o')
ax[1,0].set_ylabel('$D (m^{2}s^{-1})$', fontsize=14)
ax[1,0].set_xlabel('$R - R_{sep} (m)$', fontsize=14)
ax[1,0].tick_params(axis='y', labelsize=14)
ax[1,0].tick_params(axis='x', labelsize=14)

#ax[0,1].errorbar(exp['on']['te']['X'], exp['on']['te']['Y'], exp['on']['te']['Y_unc'])
#ax[0,1].errorbar(exp['off']['te']['X'], exp['off']['te']['Y'], exp['off']['te']['Y_unc'])
ax[0,1].plot(SOLPS['te']['X'], SOLPS['te']['Y'])
ax[0,1].plot(exp['on']['te']['X'], exp['on']['te']['Y'])
#ax[0,1].plot(exp['off']['te']['X'], exp['off']['te']['Y'],'.')
ax[0,1].set_ylabel('$T_{e} (eV)$', fontsize=14)
ax[0,1].tick_params(axis='y', labelsize=14)
ax[0,1].tick_params(axis='x', labelsize=14)

ax[1,1].plot(transport['hcib']['X'], transport['hcib']['Y'],'-o')
ax[1,1].plot(transport['hce0']['X'], transport['hce0']['Y'],'-o')
ax[1,1].set_ylabel('$\chi_{i,e} (m^{2}s^{-1})$', fontsize=14)
ax[1,1].set_xlabel('$R - R_{sep} (m)$', fontsize=14)
ax[1,1].tick_params(axis='y', labelsize=14)
ax[1,1].tick_params(axis='x', labelsize=14)

ax[0,1].set_xlim(SOLPS[kw]['X'][0],SOLPS[kw]['X'][-1]+1e-2)

plt.show()










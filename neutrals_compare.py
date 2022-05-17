import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import sys
import glob

params = ['ne', 'te']
shot_on = 1070614013
shot_off = 1070614016
#experiment = 'pump_on'
#experiment = 'pump_on_MC'
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
exp['on'], exp['off'] = {'nn':{}, 'S_ion':{}, 'flux':{}}, {'nn':{}, 'S_ion':{}, 'flux':{}}

nn_on, nn_off = exp['on']['nn'], exp['off']['nn']
S_ion_on, S_ion_off = exp['on']['S_ion'], exp['off']['S_ion']
flux_on, flux_off = exp['on']['flux'], exp['off']['flux']

psin_nn_on, psin_nn_off = exp_on_obj[0]**2, exp_off_obj[0]**2
R_nn_on, R_nn_off = exp_on_obj[2], exp_off_obj[2]
nn_data_on, nn_data_off = exp_on_obj[3], exp_off_obj[3]
nn_unc_on, nn_unc_off = exp_on_obj[4], exp_off_obj[4]

psin_S_ion_on, psin_S_ion_off = exp_on_obj[0]**2, exp_off_obj[0]**2
R_S_ion_on, R_S_ion_off = exp_on_obj[2], exp_off_obj[2]
S_ion_data_on, S_ion_data_off = exp_on_obj[5], exp_off_obj[5]
S_ion_unc_on, S_ion_unc_off = exp_on_obj[6], exp_off_obj[6]

#nn_on['X'], nn_off['X'] = psin_nn_on, psin_nn_off
nn_on['X'], nn_off['X'] = R_nn_on, R_nn_off
nn_on['Y'], nn_off['Y'] = nn_data_on*1e6, nn_data_off*1e6
nn_on['Y_unc'], nn_off['Y_unc'] = nn_unc_on*1e6, nn_unc_off*1e6

#S_ion_on['X'], S_ion_off['X'] = psin_S_ion_on, psin_S_ion_off
S_ion_on['X'], S_ion_off['X'] = R_S_ion_on, R_S_ion_off
S_ion_on['Y'], S_ion_off['Y'] = S_ion_data_on*1e6, S_ion_data_off*1e6
S_ion_on['Y_unc'], S_ion_off['Y_unc'] = S_ion_unc_on*1e6, S_ion_unc_off*1e6

# try to compute flux from exp data
'''
#remove nans first only from flux arrays
mask_on = ~np.isnan(exp['on']['S_ion']['Y'])
mask_off = ~np.isnan(exp['off']['S_ion']['Y'])

R_exp_on = exp['on']['S_ion']['X'][mask_on]
R_exp_off = exp['off']['S_ion']['X'][mask_off]

exp_flux_on = np.zeros_like(R_exp_on)
exp_flux_off = np.zeros_like(R_exp_off)


cum_int = 0
for i in range(1,len(exp['on']['S_ion']['Y'][mask_on])):

	avg = (exp['on']['S_ion']['Y'][mask_on][i] + exp['on']['S_ion']['Y'][mask_on][i-1])/2
	dx = R_exp_on[i] - R_exp_on[i-1]
	cum_int += avg*dx
	exp_flux_on[i] = cum_int

cum_int = 0
for i in range(1,len(exp['off']['S_ion']['Y'][mask_off])):

	avg = (exp['off']['S_ion']['Y'][mask_off][i] + exp['off']['S_ion']['Y'][mask_off][i-1])/2
	dx = R_exp_off[i] - R_exp_off[i-1]
	cum_int += avg*dx
	exp_flux_off[i] = cum_int

	
exp['on']['flux']['X'] = R_exp_on
exp['off']['flux']['X'] = R_exp_off

exp['on']['flux']['Y'] = exp_flux_on
exp['off']['flux']['Y'] = exp_flux_off
'''

flux_on_file = '{}/exp_data/flux_{}.pkl'.format(SOLPSWORK, shot_on)
flux_off_file = '{}/exp_data/flux_{}.pkl'.format(SOLPSWORK, shot_off)
with open(flux_on_file, 'rb') as fon:
	flux_on_obj = pkl.load(fon)
with open(flux_off_file, 'rb') as foff:
	flux_off_obj = pkl.load(foff)

exp['on']['flux']['X'] = flux_on_obj[0]
exp['off']['flux']['X'] = flux_off_obj[0]


exp['on']['flux']['Y'] = flux_on_obj[1]*1e6
exp['off']['flux']['Y'] = flux_off_obj[1]*1e6


## SOLPS stuff

SOLPS = {'nn':{}, 'S_ion':{}, 'flux':{}, 'fna':{}}
SOLPS['nn']['B2'] = {}
SOLPS['nn']['EIR'] = {}

import aurora

b2path = '{}/{}/{}/{}'.format(SOLPSWORK, shot_on, experiment, name)

so = aurora.solps_case(
	b2fstate_path='{}/b2fstate'.format(b2path),
        b2fgmtry_path='{}/b2fgmtry'.format(b2path),
#	geqdsk=glob.glob('{}/{}/baserun/g{}.00793_610'.format(SOLPSWORK,shot,shot))
	geqdsk='/nobackup1/millerma/solps-iter/runs/1070614013/pump_on/baserun/g1070614013.00793_610'
)

out = so.get_radial_prof(so.fort44['dab2'])
#SOLPS['nn']['X'] = out[2]**2                   
#SOLPS['nn']['Y'] = out[3]                   

R_LFS = so.data('cr')[:,38]
rhop_LFS = aurora.coords.get_rhop_RZ(R_LFS, np.zeros_like(R_LFS), so.geqdsk)
#SOLPS['nn']['B2']['X'] = rhop_LFS**2
SOLPS['nn']['B2']['X'] = R_LFS
SOLPS['nn']['B2']['Y'] = so.fort44['dab2'][38]

Z_LFS = so.data('cz')[:,38]

# get flux from SOLPS
SOLPS['flux']['B2'] = {}
SOLPS['flux']['B2']['X'] = R_LFS
SOLPS['flux']['B2']['Y'] = -so.fort44['rfluxa'][38]
#div_flux = np.zeros_like(rad_flux[:,0])
#div_flux[1:] = np.diff(rad_flux[:,0])/np.diff(R_LFS)
#div_flux[0] = div_flux[1]

#SOLPS['S_ion']['X'] = rhop_LFS**2
#SOLPS['S_ion']['Y'] = -rad_flux


# get source from SOLPS
import netCDF4 as nc
fn = '{}/balance.nc'.format(b2path)
ds = nc.Dataset(fn)

sna = ds['eirene_mc_papl_sna_bal']
vol = ds['vol'] # cell volumes
crx = ds['crx'] # x-coords
cry = ds['cry'] # y-coords

sna_sum = np.sum(sna,axis=0) # sum over EIRENE strata
sna_Dplus_vol = sna_sum[1]/vol # get source per vol

sna_OMP = sna_Dplus_vol[:,39]

R_OMP = np.mean(crx[:,:,39], axis=0) # average the 4 corners
Z_OMP = np.mean(cry[:,:,39], axis=0)
rhop_SOLPS = aurora.coords.get_rhop_RZ(R_OMP, Z_OMP, so.geqdsk)

#SOLPS['S_ion']['X'] = rhop_SOLPS**2
SOLPS['S_ion']['X'] = R_OMP 
SOLPS['S_ion']['Y'] = sna_OMP

# cell faces for flux calculation
hx = np.array(ds['hx'])
hy = np.array(ds['hy'])
hz = np.array(ds['hz'])

rad_SA = vol/hy

rad_SA=1

fna = ds['fna_tot'][1] # gets flux for ions
fnay_m2 = fna[1]/rad_SA # gets y-directed (radial) flux in m^{-2}s^{-1}

fnay_OMP = fnay_m2[:,39]

SOLPS['fna']['B2'] = {}
SOLPS['fna']['B2']['X'] = R_OMP
SOLPS['fna']['B2']['Y'] = fnay_OMP


### try to get triangles info ###
edena = so.fort46['edena']
momy = so.fort46['vydena']
inds = np.where(so.ynodes >= so.ynodes.max()-0.001)[0] # what is the index corresponding to the vertex (x,y)
tris = so.triangles.astype(int)
xtri = []
xcentroid = []
ytri = []
ptri = []
mtri = []
for i in inds:
	num = np.where(tris==i)[0] # what is the number of triangle that has that vertex
	corner = np.where(tris==i)[1]
	for j in range(len(num)):
		xcentroid.append(so.xnodes[tris[num[j]]].mean())
		ptri.append(edena[num[j]])
		mtri.append(momy[num[j]])

# plot energy density at the top
#fig,ax = plt.subplots()
#ax.plot(xcentroid,ptri,'o')
#plt.show()

# get pressure in cryopump baffle
R_Lend = 0.574 # left end of baffle
R_Rend = 0.641 # right end of baffle
xcentroid = np.array(xcentroid)
ptri = np.array(ptri)
mtri = np.array(mtri)
mask = np.logical_and(xcentroid > R_Lend, xcentroid < R_Rend)

avg_edena_pump = np.mean(ptri[mask]) * 1e6 # convert to ev/m^3
avg_edena_pump_mtorr = avg_edena_pump*1.602e-19/133*1e3 # convert to Pa --> Torr --> mTorr
Pavg_pump = 2/3*avg_edena_pump_mtorr # pressure in mTorr

print('pump pressure (mTorr): {}'.format(Pavg_pump))

#fig,ax = plt.subplots()
#ax.plot(xcentroid,ptri,'o')
#ax.plot(xcentroid[mask],ptri[mask],'x',c='r')

#fig,ax = plt.subplots()
#ax.plot(xcentroid,mtri,'o')
#ax.plot(xcentroid[mask],mtri[mask],'x',c='r')
#plt.show()


# get full nn from eirene

xcentt = []
ycentt = []
for t in tris:
	xcentt.append(so.xnodes[t].mean())
	ycentt.append(so.ynodes[t].mean())

xcentt = np.array(xcentt)
ycentt = np.array(ycentt)
mask_OMP_tris = np.logical_and(np.abs(ycentt) < 0.01, xcentt > 0.7)

R_EIRENE_OMP = xcentt[mask_OMP_tris]
Z_EIRENE_OMP = ycentt[mask_OMP_tris]
nn_EIRENE_OMP = so.fort46['pdena'][mask_OMP_tris,0]

rhop_EIRENE_LFS = aurora.coords.get_rhop_RZ(R_EIRENE_OMP,Z_EIRENE_OMP,so.geqdsk)
#SOLPS['nn']['EIR']['X'] = rhop_EIRENE_LFS**2
SOLPS['nn']['EIR']['X'] = R_EIRENE_OMP
SOLPS['nn']['EIR']['Y'] = nn_EIRENE_OMP*1e6


###

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

		transport[str(tparams[int(ind)-1])]['X'].append(float(xcoord[0].strip())+1)
		transport[str(tparams[int(ind)-1])]['Y'].append(float(ycoord[0].strip()))

R_sep_SOLPS = aurora.rad_coord_transform(1,'rhop','Rmid',so.geqdsk)

# stuff for error bar plotting
ff = 1./np.log(10.)

fig, ax = plt.subplots(2, sharex=True)
ax[0].errorbar(exp['on']['nn']['X'], np.log(exp['on']['nn']['Y']), ff*exp['on']['nn']['Y_unc']/exp['on']['nn']['Y'])
ax[0].errorbar(exp['off']['nn']['X'], np.log(exp['off']['nn']['Y']), ff*exp['off']['nn']['Y_unc']/exp['off']['nn']['Y'])
#ax[0].semilogy(exp['on']['nn']['X'], exp['on']['nn']['Y'])
#ax[0].semilogy(exp['off']['nn']['X'], exp['off']['nn']['Y'])
ax[0].plot(SOLPS['nn']['B2']['X'], np.log(SOLPS['nn']['B2']['Y']),'.')
ax[0].plot(SOLPS['nn']['EIR']['X'], np.log(SOLPS['nn']['EIR']['Y']),'.')
ax[0].legend(['B2', 'EIRENE', '1070614013 (on)'])
ax[0].axvline(R_sep_SOLPS,linestyle='--',color='gray')
#ax[0,0].plot(SOLPS['nn']['X'], SOLPS['nn']['Y'])
ax[0].set_ylabel('$log(n_{D}) (m^{-3})$', fontsize=14)
ax[0].tick_params(axis='y', labelsize=14)

#ax[1,0].plot(transport['dna0']['X'], transport['dna0']['Y'],'-o')
#ax[1,0].legend(['$D$'])

ax[1].errorbar(exp['on']['S_ion']['X'], np.log(exp['on']['S_ion']['Y']), ff*exp['on']['S_ion']['Y_unc']/exp['on']['S_ion']['Y'], fmt='.')
ax[1].errorbar(exp['off']['S_ion']['X'], np.log(exp['off']['S_ion']['Y']), ff*exp['off']['S_ion']['Y_unc']/exp['off']['S_ion']['Y'], fmt='.')
#ax[1].semilogy(exp['on']['S_ion']['X'], exp['on']['S_ion']['Y'])
#ax[1].semilogy(exp['off']['S_ion']['X'], exp['off']['S_ion']['Y'])
ax[1].plot(SOLPS['S_ion']['X'], np.log(SOLPS['S_ion']['Y']),'.')
ax[1].axvline(R_sep_SOLPS,linestyle='--',color='gray')
#ax[0,1].plot(SOLPS['te']['X'], SOLPS['te']['Y'])
#ax[0,1].set_ylabel('$S_{ion}$')
ax[1].set_ylabel('$log(S_{D^{+}}) (m^{-3}s^{-1})$', fontsize=14)
ax[1].set_xlabel('$R (m)$', fontsize=14)
ax[1].tick_params(axis='y', labelsize=14)
ax[1].tick_params(axis='x', labelsize=14)

ax[0].set_xlim([0.85,0.92])

ax[0].set_ylim([31,40])
ax[1].set_ylim([47,52])

fig, ax = plt.subplots()
#ax.plot(exp['on']['flux']['X'], exp['on']['flux']['Y'])
#ax.plot(exp['off']['flux']['X'], exp['off']['flux']['Y'])
ax.plot(SOLPS['flux']['B2']['X'], SOLPS['flux']['B2']['Y'],'.')
ax.plot(SOLPS['fna']['B2']['X'], SOLPS['fna']['B2']['Y'],'.')
ax.axvline(R_sep_SOLPS,linestyle='--',color='gray')
#ax.legend(['1070614013 (on)', '1070614016 (off)','fort44', 'balance'])
ax.legend(['fort44', 'balance'])
ax.set_ylabel('$\Gamma_{D^{+}}$', fontsize=14)
#ax[1].plot(transport['hcib']['X'], transport['hcib']['Y'],'-o')
#ax[1].plot(transport['hce0']['X'], transport['hce0']['Y'],'-o')
#ax[1].legend(['$\chi_{i}$','$\chi_{e}$'])

#ax[0,1].set_xlim([SOLPS[kw]['X'][0],SOLPS[kw]['X'][-1]+0.05])

#plt.show()











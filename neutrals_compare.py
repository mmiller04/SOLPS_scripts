import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import sys
import glob

params = ['ne', 'te']
shot = 1070614013
#experiment = 'pump_on'
#experiment = 'pump_on_MC'
#experiment = 'pump_off'
#shot = 1120917011
experiment = sys.argv[1]
attempt = sys.argv[2]

SOLPSWORK = '/nobackup1/users/millerma/solps-iter/runs'

exp_file = '{}/exp_data/lyman_data_{}.pkl'.format(SOLPSWORK, shot)
with open(exp_file, 'rb') as f:
	exp_obj = pkl.load(f)

# put into dictionary
exp = {'nn':{}, 'S_ion':{}}
nn = exp['nn']
S_ion = exp['S_ion']

#psin_ne, ne, ne_unc, psin_Te, Te, Te_unc = exp_obj
psin_nn = exp_obj[0]**2
nn_data = exp_obj[3]
nn_unc = exp_obj[4]
psin_S_ion = exp_obj[0]**2
S_ion_data = exp_obj[5]
S_ion_unc = exp_obj[6]

nn['X'] = psin_nn
nn['Y'], nn['Y_unc'] = nn_data*1e6, nn_unc*1e6
S_ion['X'] = psin_S_ion
S_ion['Y'], S_ion['Y_unc'] = S_ion_data, S_ion_unc


SOLPS = {'nn':{}, 'S_ion':{}}
SOLPS['nn']['B2'] = {}
SOLPS['nn']['EIR'] = {}


import aurora

b2path = '{}/{}/{}/attempt{}'.format(SOLPSWORK, shot, experiment, attempt)

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

# get flux from SOLPS
rad_flux = so.fort44['rfluxa'][38]
#div_flux = np.zeros_like(rad_flux[:,0])
#div_flux[1:] = np.diff(rad_flux[:,0])/np.diff(R_LFS)
#div_flux[0] = div_flux[1]

#SOLPS['S_ion']['X'] = rhop_LFS**2
#SOLPS['S_ion']['Y'] = -rad_flux

# try to compute flux from exp data
R_exp = exp_obj[1]
exp_flux = np.zeros_like(S_ion_data)
cum_int = 0
for i in range(1,len(S_ion_data)):
	cum_int+=(S_ion_data[i-1]+S_ion_data[i])/2*(R_exp[i] - R_exp[i-1])
	exp_flux[i] = cum_int

#exp['S_ion']['Y'] = exp_flux*1e6
exp['S_ion']['Y'] = S_ion_data*1e6

# get source from SOLPS
import netCDF4 as nc
fn = '{}/balance.nc'.format(b2path)
ds = nc.Dataset(fn)

sna = ds['eirene_mc_papl_sna_bal']
snm = ds['eirene_mc_pmpl_sna_bal']
vol = ds['vol'] # cell volumes
crx = ds['crx'] # x-coords
cry = ds['cry'] # y-coords

sna_sum = np.sum(sna,axis=0) # sum over EIRENE strata
sna_Dplus_vol = sna_sum[1]/vol # get source per vol

sna_OMP = sna_Dplus_vol[:,39]

R_OMP = np.mean(crx[:,:,39], axis=0) # average the 4 corners
Z_OMP = np.mean(cry[:,:,39], axis=0)
rhop_SOLPS = aurora.coords.get_rhop_RZ(R_OMP, Z_OMP, so.geqdsk)

SOLPS['S_ion']['X'] = rhop_SOLPS**2
SOLPS['S_ion']['Y'] = sna_OMP

### try to get triangles info ###
edena = so.fort46['edena']
inds = np.where(so.ynodes >= so.ynodes.max()-0.001)[0] # what is the index corresponding to the vertex (x,y)
tris = so.triangles.astype(int)
xtri = []
xcentroid = []
ytri = []
ptri = []
for i in inds:
	num = np.where(tris==i)[0] # what is the number of triangle that has that vertex
	corner = np.where(tris==i)[1]
	for j in range(len(num)):
		xcentroid.append(so.xnodes[tris[num[j]]].mean())
		ptri.append(so.fort46['edena'][num[j]])

# plot energy density at the top
#fig,ax = plt.subplots()
#ax.plot(xcentroid,ptri,'o')
#plt.show()

# get pressure in cryopump baffle
R_Lend = 0.574 # left end of baffle
R_Rend = 0.641 # right end of baffle
xcentroid = np.array(xcentroid)
ptri = np.array(ptri)
mask = np.logical_and(xcentroid > R_Lend, xcentroid < R_Rend)

avg_edena_pump = np.mean(ptri[mask]) * 1e6 # convert to ev/m^3
avg_edena_pump_mtorr = avg_edena_pump*1.602e-19/133*1e3 # convert to Pa --> Torr --> mTorr
Pavg_pump = 2/3*avg_edena_pump_mtorr # pressure in mTorr

print('pump pressure (mTorr): {}'.format(Pavg_pump))


fig,ax = plt.subplots()
ax.plot(xcentroid,ptri,'o')
ax.plot(xcentroid[mask],ptri[mask],'x',c='r')
plt.show()


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

transport_file = '{}/{}/{}/attempt{}/b2.transport.inputfile'.format(SOLPSWORK, shot, experiment, attempt)
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

fig, ax = plt.subplots(2,2, sharex=True)
#ax[0,0].errorbar(exp['nn']['X'], exp['nn']['Y'], exp['nn']['Y_unc'], fmt='.')
#ax[0,0].semilogy(exp['nn']['X'], exp['nn']['Y'], '.')
ax[0,0].semilogy(SOLPS['nn']['B2']['X'], SOLPS['nn']['B2']['Y'],'.')
ax[0,0].semilogy(SOLPS['nn']['EIR']['X'], SOLPS['nn']['EIR']['Y'],'.')
#ax[0,0].plot(SOLPS['nn']['X'], SOLPS['nn']['Y'])
ax[0,0].legend(['exp (Ly-alpha)','SOLPS'])
ax[0,0].set_ylabel('$n_{n}$')

ax[1,0].plot(transport['dna0']['X'], transport['dna0']['Y'],'-o')
ax[1,0].legend(['$D$'])

#ax[0,1].errorbar(exp['te']['X'], exp['te']['Y'], exp['te']['Y_unc'], fmt='.')
ax[0,1].semilogy(exp['S_ion']['X'], exp['S_ion']['Y'], '.')
ax[0,1].semilogy(SOLPS['S_ion']['X'], SOLPS['S_ion']['Y'],'.')
#ax[0,1].plot(SOLPS['te']['X'], SOLPS['te']['Y'])
#ax[0,1].set_ylabel('$S_{ion}$')
ax[0,1].set_ylabel('$\S_{D^{+}}$')

ax[1,1].plot(transport['hcib']['X'], transport['hcib']['Y'],'-o')
ax[1,1].plot(transport['hce0']['X'], transport['hce0']['Y'],'-o')
ax[1,1].legend(['$\chi_{i}$','$\chi_{e}$'])

#ax[0,1].set_xlim([SOLPS[kw]['X'][0],SOLPS[kw]['X'][-1]+0.05])

#plt.show()










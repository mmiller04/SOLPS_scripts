import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import sys

params = ['ne', 'te']
shot = 1070614013
#experiment = 'pump_off'
#shot = 1120917011
experiment = sys.argv[1]
attempt = sys.argv[2]

SOLPSWORK = '/nobackup1/users/millerma/solps-iter/runs'

exp_file = '{}/exp_data/lyman_data_{}.pkl'.format(SOLPSWORK, shot)
with open(exp_file, 'rb') as f:
	exp_obj = pkl.load(f)

# put into dictionary
exp = {'ne':{}, 'te':{}}
ne = exp['ne']
te = exp['te']

#psin_ne, ne, ne_unc, psin_Te, Te, Te_unc = exp_obj
psin_ne = exp_obj[0]**2
ne_data = exp_obj[7]
ne_unc = exp_obj[8]
psin_Te = exp_obj[0]**2
Te_data = exp_obj[9]
Te_unc = exp_obj[10]

ne['X'] = psin_ne
ne['Y'], ne['Y_unc'] = ne_data*1e6, ne_unc*1e6 # convert to m^{-3}
te['X'] = psin_Te
te['Y'], te['Y_unc'] = Te_data, Te_unc # convert to eV

SOLPS = {'ne':{}, 'te':{}}

for kw in params:

	SOLPS_file = '{}/{}/{}/attempt{}/{}3da.last10'.format(SOLPSWORK, shot, experiment, attempt, kw)
	fid_SOLPS = open(SOLPS_file, 'r')
	SOLPS_obj = np.genfromtxt(fid_SOLPS)
	
	SOLPS[kw]['X'] = SOLPS_obj[:,0] + 1 # in units of psi - psi_sep
	SOLPS[kw]['Y'] = SOLPS_obj[:,1]

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
ax[0,0].errorbar(exp['ne']['X'], exp['ne']['Y'], exp['ne']['Y_unc'], fmt='.')
ax[0,0].plot(SOLPS['ne']['X'], SOLPS['ne']['Y'])
ax[0,0].legend(['SOLPS','exp (TS)'])
ax[0,0].set_ylabel('ne')

ax[1,0].plot(transport['dna0']['X'], transport['dna0']['Y'],'-o')
ax[1,0].legend(['$D$'])

ax[0,1].errorbar(exp['te']['X'], exp['te']['Y'], exp['te']['Y_unc'], fmt='.')
ax[0,1].plot(SOLPS['te']['X'], SOLPS['te']['Y'])
ax[0,1].set_ylabel('te')

ax[1,1].plot(transport['hcib']['X'], transport['hcib']['Y'],'-o')
ax[1,1].plot(transport['hce0']['X'], transport['hce0']['Y'],'-o')
ax[1,1].legend(['$\chi_{i}$','$\chi_{e}$'])

ax[0,1].set_xlim([SOLPS[kw]['X'][0],SOLPS[kw]['X'][-1]+0.05])

plt.show()










import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl

params = ['ne', 'te']
shot = 1120917011
attempt = 3

SOLPSWORK = '/nobackup1/users/millerma/solps-iter/runs'

exp_file = '{}/exp_data/lyman_data_{}.pkl'.format(SOLPSWORK, shot)
with open(exp_file, 'rb') as f:
	exp_obj = pkl.load(f)

# put into dictionary
exp = {'X':{}, 'Y':{}}
X = exp['X']
Y = exp['Y']

psin_ne, ne, ne_unc, psin_Te, Te, Te_unc = exp_obj

X['ne'] = psin_ne
X['te'] = psin_Te
Y['ne'], Y['ne_unc'] = ne*1e20, ne_unc*1e20 # convert to m^{-3}
Y['te'], Y['te_unc'] = Te*1e3, Te_unc*1e3 # convert to eV


SOLPS = {'X':{}, 'Y':{}}
X = SOLPS['X']
Y = SOLPS['Y']

for kw in params:

	SOLPS_file = '{}/{}/attempt{}/{}3da.last10'.format(SOLPSWORK, shot, attempt, kw)
	fid_SOLPS = open(SOLPS_file, 'r')
	SOLPS_obj = np.genfromtxt(fid_SOLPS)

	X[kw] = SOLPS_obj[:,0] + 1 # in units of psi - psi_sep
	Y[kw] = SOLPS_obj[:,1]

fig, ax = plt.subplots(2,1, sharex=True)
#ax[0].errorbar(exp['X']['ne'], exp['Y']['ne'], exp['Y']['ne_unc'], fmt='.')
ax[0].plot(exp['X']['ne'], exp['Y']['ne'], 'o')
ax[0].plot(SOLPS['X']['ne'], SOLPS['Y']['ne'])
ax[0].set_ylabel('ne')

#ax[1].errorbar(exp['X']['te'], exp['Y']['te'], exp['Y']['te_unc'], fmt='.')
ax[1].plot(exp['X']['te'], exp['Y']['te'], 'o')
ax[1].plot(SOLPS['X']['te'], SOLPS['Y']['te'])
ax[1].set_ylabel('te')

ax[1].set_xlim([SOLPS['X'][kw][0],SOLPS['X'][kw][-1]])

plt.show()










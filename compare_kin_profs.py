import numpy as np
import matplotlib.pyplot
import pickle as pkl

params = ['ne', 'te']
shot = 1120917011
attempt = 3

SOLPSWORK = '/nobackup1/users/millerma/solps-iter/runs'

exp_file = '{}/exp_data/lyman_data_{}.pkl'.format(SOLPSWORK, shot)
fid_exp = open(exp_file, 'rb')
exp_obj = pkl.load(fid_exp)

# put into dictionary
exp = {'X':{}, 'Y':{}}
X = exp['X']
Y = exp['Y']

psin_ne, ne, ne_unc, psin_Te, Te, Te_unc = exp_obj

X['ne'] = psin_ne
X['te'] = psin_te
Y['ne'], Y['ne_unc'] = ne, ne_unc
Y['te'], Y['te_unc'] = te, te_unc


SOLPS = {'X':{}, 'Y':{}}
X = exp['X']
Y = exp['Y']

for kw in params:

	SOLPS_om_file = '{}/{}/attempt{}/{}3da.last10'.format(SOLPSWORK, shot, attempt, kw)
	fid_SOLPS = open(SOLPS_file, 'r')
	SOLPS_obj = np.genfromtxt(fid_SOLPS)

	X[kw] = SOLPS_obj[:,0]
	Y[kw] = SOLPS_obj[:,1]

fig, ax = plt.subplots(2,1, sharex=True)
ax[0].errorbar(exp['X']['ne'], exp['Y']['ne'], exp['Y']['ne_unc'], fmt='.')
ax[0].plot(SOLPS['X']['ne'], SOLPS['Y']['ne'])
ax[0].set_ylabel('ne')

ax[1].errorbar(exp['X']['te'], exp['Y']['te'], exp['Y']['te_unc'], fmt='.')
ax[1].plot(SOLPS['X']['te'], SOLPS['Y']['te'])
ax[1].set_ylabel('te')












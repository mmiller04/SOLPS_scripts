import pickle as pkl

SOLPSWORK = '/nobackup1/users/millerma/solps-iter/runs'

transport_file = '{}/exp_data/transport_{}.pkl'.format(SOLPSWORK, shot)
with open(transport_file, 'rb') as f:
	tr_obj = pk.load(f)

rhop, R, flux_ion, D_eff, v_eff = tr_obj


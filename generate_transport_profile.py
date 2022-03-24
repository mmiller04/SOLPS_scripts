import pickle as pkl
import numpy as np
from scipy.interpolate import interp1d

shot = 1120917011
attempt = 4

SOLPSWORK = '/nobackup1/millerma/solps-iter/runs'

transport_file = '{}/exp_data/transport_{}.pkl'.format(SOLPSWORK, shot)
with open(transport_file, 'rb') as f:
	tr_obj = pkl.load(f)

rhop, R, flux_ion, D_eff, v_eff = tr_obj

# check how many grid points are in y-direction (ny) - limit to number of points to specify
carre_file = '{}/{}/baserun/carre.dat'.format(SOLPSWORK, shot)
with open(carre_file, 'r') as f:
	big_list = []
	for line in f.readlines():
		big_list.append(line.split('='))
	ny = 0
	for split_line in big_list:
		ny += int(split_line[1].strip()) if split_line[0].strip() == ('npr(2)' or 'npr(3)') else 0

R_sep = interp1d(rhop, R, kind='cubic')(1)

R_fine = R - R_sep
R_SOLPS = np.linspace(R_fine[0], R_fine[-1], ny)

D_eff_SOLPS = interp1d(R_fine, D_eff, kind='cubic')(R_SOLPS)
v_eff_SOLPS = interp1d(R_fine, v_eff, kind='cubic')(R_SOLPS)

# from manual, try multiplying D*1.4 for chi
chii_SOLPS = D_eff_SOLPS*1.4
chie_SOLPS = D_eff_SOLPS*1.4

# define which coeficients to prescribe
tc = {} # transport coefficient dictionary
tc['inds'], tc['data'] = {}, {} # sub directories

tc['inds']['dna0'] = 1 # particle density-driven diffusivity \D_{perp}^{na}
tc['inds']['dpa0'] = 2 # particle pressure-driven diffusivity \D_{perp}^{pa}
tc['inds']['hcib'] = 3 # ion thermal anomalous diffusivity \chi_{i}
tc['inds']['hce0'] = 4 # electron thermal anomalous diffusivity \chi_{e}
tc['inds']['vla0_x'] = 5 # x-component of the anomalous "pinch" velocity
tc['inds']['vla0_y'] = 6 # y-component of the anomalous "pinch" velocity
tc['inds']['vsa0'] = 7 # anomalous viscosity
tc['inds']['sig0'] = 8 # anomalous radial electrical conductivity
tc['inds']['alf0'] = 9 # anomalous radial thermo-electric coefficient

coefs = ['dna0', 'hcib', 'hce0', 'vla0_x']
tc['data']['dna0'] = D_eff_SOLPS
tc['data']['hcib'] = chii_SOLPS
tc['data']['hce0'] = chie_SOLPS
tc['data']['vla0_x'] = v_eff_SOLPS

# write to a file in new attempt directory
outfile = '{}/{}/attempt{}/b2.transport.inputfile'.format(SOLPSWORK, shot, attempt)
with open(outfile, 'w') as f:
	f.write(' &TRANSPORT\n')

	for c in coefs:

		f.write(' ndata(1, {}, 1 )= {} ,\n'.format(tc['inds'][c], ny))

		for yy in range(ny):

			f.write(' tdata(1, {}, {}, 1 )= {:.2E} , tdata(2, {}, {}, 1 )= {:.2f} ,\n'\
				.format(yy, tc['inds'][c], R_SOLPS[yy], yy, tc['inds'][c], tc['data'][c][yy]))

	f.write(' no_pflux=.true.\n')
	f.write(' /')








	



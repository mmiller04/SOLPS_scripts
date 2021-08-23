import pickle as pkl
from scipy.interpolate import interp1d

shot = 1120917011
attempt = 4

SOLPSWORK = '/nobackup1/users/millerma/solps-iter/runs'

transport_file = '{}/exp_data/transport_{}.pkl'.format(SOLPSWORK, shot)
with open(transport_file, 'rb') as f:
	tr_obj = pk.load(f)

rhop, R, flux_ion, D_eff, v_eff = tr_obj

# check how many grid points are in y-direction (ny) - limit to number of points to specify
carre_file = '{}/{}/baserun/carre.dat'
with open(carre_file, 'r') as f:
	big_list = []
	for line in f.readlines():
		big_list.append(line.split('='))
	ny = 0
	for split_line in big_list:
		ny += int(split_line[1].strip()) if split_line[0].strip() == ('npr(2)' or 'npr(3)') else 0

R_sep = interp1d(R, rhop, kind='cubic')(1)

R_fine = R - R_sep
R_SOLPS = np.linspace(R_fine[0], R_fine[-1], ny)

D_eff_SOLPS = np.interp1d(R, D_eff, kind='cubic')(R_SOLPS)
v_eff_SOLPS = np.interp1d(R, v_eff, kind='cubic')(R_SOLPS)

# define which coeficients to prescribe
tc = {} # transport coefficient dictionary
tc['dna0'] = 1 # particle density-driven diffusivity \D_{perp}^{na}
tc['dpa0'] = 2 # particle pressure-driven diffusivity \D_{perp}^{pa}
tc['hcib'] = 3 # ion thermal anomalous diffusivity \chi_{i}
tc['hce0'] = 4 # electron thermal anomalous diffusivity \chi_{e}
tc['vla0_x'] = 5 # x-component of the anomalous "pinch" velocity
tc['vla0_y'] = 6 # y-component of the anomalous "pinch" velocity
tc['vsa0'] = 7 # anomalous viscosity
tc['sig0'] = 8 # anomalous radial electrical conductivity
tc['alf0'] = 9 # anomalous radial thermo-electric coefficient

coefs = ['dna0, vla0_x']

# write to a file in new attempt directory
outfile = '{}/{}/attempt{}/b2.transport.inputfile'.format(SOLPSWORK, shot, attempt)
with open(outfile, 'wb') as f:
	f.write(' &TRANSPORT\n')

	for c in coefs:

		f.write(' ndata(1, {}, 1 )= {} ,\n'.format(tc[c], ny))

		for yy in range(ny):

			f.write(' tdata(1, {}, {}, 1 )= {:.2E} , tdata(2, {}, {}, 1 )= {:.2f} ,\n'\
				.format(yy, tc[c], R_SOLPS[yy], yy, tc[c], D_eff_SOLPS[yy]))

	f.write(' no_pflux=.true.\n')
	f.write(' /')








	



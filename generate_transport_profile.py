import pickle as pkl

shot = 1120917011
attempt = 4

SOLPSWORK = '/nobackup1/users/millerma/solps-iter/runs'

transport_file = '{}/exp_data/transport_{}.pkl'.format(SOLPSWORK, shot)
with open(transport_file, 'rb') as f:
	tr_obj = pk.load(f)

rhop, R, flux_ion, D_eff, v_eff = tr_obj

# write to a file in new attempt directory
outfile = '{}/{}/attempt{}/b2.transport.inputfile'.format(SOLPSWORK, shot, attempt)

tc = {}
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
	



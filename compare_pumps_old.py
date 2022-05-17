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
attempt = sys.argv[1]

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


## SOLPS stuff

import aurora

b2path0 = '{}/{}/pump{}/attempt{}'.format(SOLPSWORK, shot_on, 0, attempt)
#b2path4 = '{}/{}/pump{}/attempt{}'.format(SOLPSWORK, shot_on, 4, attempt)
b2path9 = '{}/{}/pump{}/attempt{}'.format(SOLPSWORK, shot_on, 8, attempt)


so0 = aurora.solps_case(
	b2fstate_path='{}/b2fstate'.format(b2path0),
        b2fgmtry_path='{}/b2fgmtry'.format(b2path0),
	geqdsk='/nobackup1/millerma/solps-iter/runs/1070614013/pump_on/baserun/g1070614013.00793_610'
)

#so4 = aurora.solps_case(
#	b2fstate_path='{}/b2fstate'.format(b2path4),
#        b2fgmtry_path='{}/b2fgmtry'.format(b2path4),
#	geqdsk='/nobackup1/millerma/solps-iter/runs/1070614013/pump_on/baserun/g1070614013.00793_610'
#)

so9 = aurora.solps_case(
	b2fstate_path='{}/b2fstate'.format(b2path9),
        b2fgmtry_path='{}/b2fgmtry'.format(b2path9),
	geqdsk='/nobackup1/millerma/solps-iter/runs/1070614013/pump_on/baserun/g1070614013.00793_610'
)


R0_LFS = so0.data('cr')[:,38]
flux0 = -so0.fort44['rfluxa'][38]
nn0 = so0.fort44['dab2'][38]

#R4_LFS = so4.data('cr')[:,38]
#flux4 = -so4.fort44['rfluxa'][38]
#nn4 = so4.fort44['dab2'][38]

R9_LFS = so9.data('cr')[:,38]
flux9 = -so9.fort44['rfluxa'][38]
nn9 = so9.fort44['dab2'][38]


# get source from SOLPS
import netCDF4 as nc

fn0 = '{}/balance.nc'.format(b2path0)
ds0 = nc.Dataset(fn0)

#fn4 = '{}/balance.nc'.format(b2path4)
#ds4 = nc.Dataset(fn4)

fn9 = '{}/balance.nc'.format(b2path9)
ds9 = nc.Dataset(fn9)

vol = ds0['vol'] # cell volumes
crx = ds0['crx'] # cell volumes

sna0 = ds0['eirene_mc_papl_sna_bal']
sna0_sum = np.sum(sna0,axis=0) # sum over EIRENE strata
sna0_Dplus_vol = sna0_sum[1]/vol # get source per vol

sna0_OMP = sna0_Dplus_vol[:,39]
R0_OMP = np.mean(crx[:,:,39], axis=0) # average the 4 corners


#sna4 = ds4['eirene_mc_papl_sna_bal']
#sna4_sum = np.sum(sna4,axis=0) # sum over EIRENE strata
#sna4_Dplus_vol = sna4_sum[1]/vol # get source per vol

#sna4_OMP = sna4_Dplus_vol[:,39]
#R4_OMP = np.mean(crx[:,:,39], axis=0) # average the 4 corners


sna9 = ds9['eirene_mc_papl_sna_bal']
sna9_sum = np.sum(sna9,axis=0) # sum over EIRENE strata
sna9_Dplus_vol = sna9_sum[1]/vol # get source per vol

sna9_OMP = sna9_Dplus_vol[:,39]
R9_OMP = np.mean(crx[:,:,39], axis=0) # average the 4 corners


R_sep_SOLPS = aurora.rad_coord_transform(1,'rhop','Rmid',so0.geqdsk)

# stuff for error bar plotting
ff = 1./np.log(10.)

fig, ax = plt.subplots()
#ax.errorbar(exp['on']['nn']['X'], np.log(exp['on']['nn']['Y']), ff*exp['on']['nn']['Y_unc']/exp['on']['nn']['Y'])
#ax.errorbar(exp['off']['nn']['X'], np.log(exp['off']['nn']['Y']), ff*exp['off']['nn']['Y_unc']/exp['off']['nn']['Y'])
ax.semilogy(exp['on']['nn']['X'], exp['on']['nn']['Y'])
ax.semilogy(exp['off']['nn']['X'], exp['off']['nn']['Y'])
ax.semilogy(R9_LFS,nn9,'.')
ax.semilogy(R0_LFS,nn0,'.')
#ax.semilogy(R4_LFS,nn4,'.')
ax.legend(['1070614013 (on)', '1070614016 (off)', 'R=0.1', 'R=0.99'])
#ax.legend(['pump0', 'pump4', 'pump9', '1070614013 (on)', '1070614016 (off)'])
ax.axvline(R_sep_SOLPS,linestyle='--',color='gray')
#ax[0,0].plot(SOLPS['nn']['X'], SOLPS['nn']['Y'])
ax.set_ylabel('$log(n_{D}) (m^{-3})$',fontsize=14)
ax.set_xlabel('$R (m)$', fontsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.tick_params(axis='x', labelsize=14)
#ax.set_ylim([31,40])

#ax[1,0].plot(transport['dna0']['X'], transport['dna0']['Y'],'-o')
#ax[1,0].legend(['$D$'])

fig, ax = plt.subplots()
#ax.errorbar(exp['on']['S_ion']['X'], np.log(exp['on']['S_ion']['Y']), ff*exp['on']['S_ion']['Y_unc']/exp['on']['S_ion']['Y'], fmt='.')
#ax.errorbar(exp['off']['S_ion']['X'], np.log(exp['off']['S_ion']['Y']), ff*exp['off']['S_ion']['Y_unc']/exp['off']['S_ion']['Y'], fmt='.')
ax.semilogy(exp['on']['S_ion']['X'], exp['on']['S_ion']['Y'])
ax.semilogy(exp['off']['S_ion']['X'], exp['off']['S_ion']['Y'])
ax.semilogy(R9_OMP,sna9_OMP,'.')
ax.semilogy(R0_OMP,sna0_OMP,'.')
#ax.semilogy(R4_OMP,sna4_OMP,'.')
ax.legend(['1070614013 (on)', '1070614016 (off)', 'R=0.1', 'R=0.99'])
#ax.legend(['pump0', 'pump4', 'pump9', '1070614013 (on)', '1070614016 (off)'])
ax.axvline(R_sep_SOLPS,linestyle='--',color='gray')
#ax[0,1].plot(SOLPS['te']['X'], SOLPS['te']['Y'])
#ax[0,1].set_ylabel('$S_{ion}$')
ax.set_ylabel('$log(S_{D^{+}}) (m^{-3}s^{-1})$', fontsize=14)
ax.set_xlabel('$R (m)$', fontsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.tick_params(axis='x', labelsize=14)
#ax.set_ylim([47,52])
#ax.set_xlim([0.85,0.92])


fig, ax = plt.subplots()
#ax.plot(exp['on']['flux']['X'], exp['on']['flux']['Y'])
#ax.plot(exp['off']['flux']['X'], exp['off']['flux']['Y'])
ax.semilogy(R9_LFS, flux9,'.')
ax.semilogy(R0_LFS, flux0,'.')
#ax.plot(R4_LFS,np.log(flux4),'.')
#ax.legend(['pump0', 'pump4', 'pump9'])
ax.legend(['R=0.1', 'R=0.99'])
ax.axvline(R_sep_SOLPS,linestyle='--',color='gray')
#ax.legend(['1070614013 (on)', '1070614016 (off)','fort44', 'balance'])
ax.set_ylabel('$-\Gamma_{D} (m^{-2}s^{-1})$', fontsize=14)
ax.set_xlabel('R (m)', fontsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.tick_params(axis='x', labelsize=14)
#ax[1].plot(transport['hcib']['X'], transport['hcib']['Y'],'-o')
#ax[1].plot(transport['hce0']['X'], transport['hce0']['Y'],'-o')
#ax[1].legend(['$\chi_{i}$','$\chi_{e}$'])

#ax[0,1].set_xlim([SOLPS[kw]['X'][0],SOLPS[kw]['X'][-1]+0.05])

#plt.show()











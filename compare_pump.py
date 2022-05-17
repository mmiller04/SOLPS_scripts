# script to compare pump on and off

from omfit_classes import omfit_eqdsk
import sys, os
import aurora

import numpy as np
import matplotlib.pyplot as plt

import neutrals_analysis as neut

shot = 1070614013
exp_on = 'transp1_on'
exp_off = 'transp1_off'
name = 'test1'

# get run paths

SOLPSWORK = '/nobackup1/millerma/solps-iter/runs'
run_on = '{}/{}/{}/{}'.format(SOLPSWORK,shot,exp_on,'test_abs1')
run_off = '{}/{}/{}/{}'.format(SOLPSWORK,shot,exp_off,name)

so_on = neut.load_aurora_socase(SOLPSWORK, shot, exp_on, name)
so_off = neut.load_aurora_socase(SOLPSWORK, shot, exp_off, name)


# plot 2d contours

# source
fig, axs = plt.subplots(1, 2, figsize=(10,6), sharex=True)
ax = axs.flatten()

source_on = neut.get_2dcontour_quantity(run_on,'source',so_on)
source_off = neut.get_2dcontour_quantity(run_off,'source',so_off)

source_min = np.minimum(source_on.min(), source_off.min())
source_max = np.maximum(source_on.max(), source_off.max())

so_on.plot2d_b2(source_on, ax=ax[0], scale='log', 
		lb=source_min, ub=source_max,
		label=r'ON')
so_off.plot2d_b2(source_off, ax=ax[1], scale='log',
		lb=source_min, ub=source_max,
		label=r'OFF')

axs[0].set_ylabel('Z (m)', fontsize=14)
axs[1].set_ylabel('Z (m)', fontsize=14)
axs[0].set_xlabel('R (m)', fontsize=14)
axs[1].set_xlabel('R (m)', fontsize=14)
axs[0].tick_params(axis='y', labelsize=14)
axs[0].tick_params(axis='x', labelsize=14)
axs[1].tick_params(axis='y', labelsize=14)
axs[1].tick_params(axis='x', labelsize=14)

fig.suptitle('$S_{D} (m^{-3}s^{-1})$')

# nn
fig, axs = plt.subplots(1, 2, figsize=(10,6), sharex=True)
ax = axs.flatten()

nn_on = neut.get_2dcontour_quantity(run_on,'nn',so_on)
nn_off = neut.get_2dcontour_quantity(run_off,'nn',so_off)

nn_min = np.maximum(np.minimum(nn_on.min(), nn_off.min()),1e10) # set floor as 1e10
nn_max = np.maximum(nn_on.max(), nn_off.max())


so_on.plot2d_eirene(nn_on, ax=ax[0], scale='log',
			lb=nn_min, ub=nn_max,
			label=r'ON')
so_off.plot2d_eirene(nn_off, ax=ax[1], scale='log',
			lb=nn_min, ub=nn_max,
			label=r'OFF')

axs[0].set_ylabel('Z (m)', fontsize=14)
axs[1].set_ylabel('Z (m)', fontsize=14)
axs[0].set_xlabel('R (m)', fontsize=14)
axs[1].set_xlabel('R (m)', fontsize=14)
axs[0].tick_params(axis='y', labelsize=14)
axs[0].tick_params(axis='x', labelsize=14)
axs[1].tick_params(axis='y', labelsize=14)
axs[1].tick_params(axis='x', labelsize=14)

fig.suptitle('$n_{D} (m^{-3})$')

# calculate pressure at the top of pump

# get neutral density profiles




# script to compare pump on and off

from omfit_classes import omfit_eqdsk
import sys, os
import aurora
import matplotlib.pyplot as plt

import neutrals_analysis as neut

shot = 1070614013
exp_on = transp5_on
exp_off = transp5_off
name = test1

# get run paths

SOLPSWORK = '/nobackup1/millerma/solps-iter/runs'
run_on = '{}/{}/{}/{}'.format(SOLPSWORK,shot,exp_on,name)
run_off = '{}/{}/{}/{}'.format(SOLPSWORK,shot,exp_off,name)


# plot 2d contours

# neutral density
fig, axs = plt.subplots(1, 2, figsize=(10,6), sharex=True)
ax = axs.flatten()

fig, ax = neut.plot2d_contour(run_on,'nn','B2',fig=fig,ax=ax)
cbar = plt.colorbar.ColorBase(ax)



# calculate pressure at the top of pump

# get neutral density profiles




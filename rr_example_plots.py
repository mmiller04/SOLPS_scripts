import sys
sys.path.append('/home/millerma/SOLPS_Scripts/SOLPS_Scripts/')
from VesselPlotterMAM import SOLPSPLOT

import numpy as np
import matplotlib.pyplot as plt

PS = ['.','.','.','.','.','-']
control012 = SOLPSPLOT('1070614013',['54'],PsinOffset=-0.005,AVG=False,PlotScheme=PS)

fig1,ax1 = plt.subplots()
control012.RadProf('Ne',Markers=False,Ax=ax1,PlotScheme=[],GRID=True)

fig2,ax2 = plt.subplots()
control012.RadProf('Te',Markers=False,Ax=ax2,PlotScheme=[],GRID=True)

fig3,ax3 = plt.subplots()
control012.RadProf('NeuDen',LOG10=2,Markers=False,Ax=ax3,PlotScheme=PS,GRID=True)

fig4,ax4 = plt.subplots()
control012.PolPlot('NeuDen',LOG10=2,Markers=False,Ax=ax4,PlotScheme=PS,GRID=True)

fig5,ax5 = plt.subplots()
control012.PolPlot('IonPol',LOG10=0,SURF=36,Markers=False,Ax=ax4,PlotScheme=PS,GRID=True)

control012.Contour('NeuDen',LOG10=2)

control012.Contour('IonPol',LOG10=1)

plt.show()

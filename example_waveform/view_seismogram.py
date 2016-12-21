# -*- coding: utf-8 -*-
#==============================================================================
# view_waveform.py
#------------------------------------------------------------------------------
# view_waveform
#
# Created by
#       Keurfon Luu <keurfon.luu@mines-paristech.fr>
#       MINES ParisTech - Centre de Géosciences
#       PSL - Research University
#
# Last updated
#       2016-12-21 15:10
#==============================================================================

# Required modules
#==================
import numpy as np
import matplotlib.pyplot as plt

# Parameters
#============
comp1_file = "./comp1.bin"
comp2_file = "./comp2.bin"
comp3_file = "./comp3.bin"
nrcv = 37

# Import shot gathers
#=====================
comp1 = np.fromfile(comp1_file, dtype = "float32")
comp2 = np.fromfile(comp2_file, dtype = "float32")
comp3 = np.fromfile(comp3_file, dtype = "float32")

comp1 = np.reshape(comp1, (len(comp1)//nrcv, nrcv), order = "F")
comp2 = np.reshape(comp2, (len(comp2)//nrcv, nrcv), order = "F")
comp3 = np.reshape(comp3, (len(comp3)//nrcv, nrcv), order = "F")

# Plot shot gathers
#===================
vmin = np.min(comp1)
vmax = np.max(comp1)

fig = plt.figure(figsize = (20, 8), facecolor = "white")
ax1 = fig.add_subplot(1, 3, 1)
ax2 = fig.add_subplot(1, 3, 2)
ax3 = fig.add_subplot(1, 3, 3)

ax1.contourf(comp1, 200, cmap = "Greys", vmin = vmin, vmax = vmax)
ax1.set_title("Component 1", loc = "right")
ax1.invert_yaxis()

ax2.contourf(comp2, 200, cmap = "Greys", vmin = vmin, vmax = vmax)
ax2.set_title("Component 2", loc = "right")
ax2.invert_yaxis()

ax3.contourf(comp3, 200, cmap = "Greys", vmin = vmin, vmax = vmax)
ax3.set_title("Component 3", loc = "right")
ax3.invert_yaxis()
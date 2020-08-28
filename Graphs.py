#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: spencer
"""

import pickle
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

from vortexAvalanches import *

matplotlib.rcParams.update({'font.size': 24})

with open('newEx3.pkl', 'rb') as input:
    vState = pickle.load(input)

with open('newExSCM3.pkl', 'rb') as input:
    SCM = pickle.load(input)

fig1=plt.figure()
fig1.set_size_inches(SCM.length,SCM.height)
ax=fig1.gca()
ax.set_xlim([0,SCM.length])
ax.set_ylim([0,SCM.height])

ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)

ax.plot(SCM.positions.T[0],SCM.positions.T[1],'x',ms=2)

for v in vState.line_history:
    A=np.asarray(v)
    
    lX,lY=A.T
    ax.plot(lX,lY,'.',color='.8',ms=.5)#linewidth=1)
    
for v in vState.vortex_list:
    A=np.asarray(v.line)
    lX,lY=A.T
    ax.plot(lX,lY,'.',color='k',ms=.5)

fig1.savefig('img.png', transparent=True)
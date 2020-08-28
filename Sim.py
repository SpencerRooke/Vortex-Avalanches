#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 14:46:21 2020

@author: spencer
"""
import pickle
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from progressbar import ProgressBar

import vortexAvalanches as vAval

matplotlib.rcParams.update({'font.size': 24})

sample = vAval.SCMaterial2d(X=12,Y=6,site_density=3, pin_radius=.2,maxpinStrength=10,f_not=.2,xBias=0)
defectNum = len(sample.positions)

fig1=plt.figure()
fig1.set_size_inches(sample.length,sample.height)
ax=fig1.gca()
ax.set_xlim([0,sample.length])
ax.set_ylim([0,sample.height])

ax.plot(sample.positions.T[0],sample.positions.T[1],'x')    
    
#vortexState = vList([])
#vortexState.addVortex(len(sample.positions),sample)
vortexState = vAval.Saturate(sample,.95)

X,Y=vAval.vlist_to_plot(vortexState)        
ax.plot(X,Y,'o')

pbar = ProgressBar()
for i in pbar(range(2000)):
    if(i%5 == 1 and (i/5) < .5*defectNum):vortexState.addVortex(1,sample)
    vortexState.updateForce(sample,dt=.01)

with open('newEx3.pkl', 'wb') as output:
    pickle.dump(vortexState, output, pickle.HIGHEST_PROTOCOL)
    
with open('newExSCM3.pkl', 'wb') as output:
    pickle.dump(sample, output, pickle.HIGHEST_PROTOCOL)

#----
fig1=plt.figure()
fig1.set_size_inches(12,6)
ax=fig1.gca()
ax.set_xlim([0,12])
ax.set_ylim([0,6])

ax.plot(sample.positions.T[0],sample.positions.T[1],'x')

for v in vortexState.line_history:
    A=np.asarray(v)
    
    lX,lY=A.T
    ax.plot(lX,lY,'.',color='.8',ms=.5)#linewidth=1)
    
for v in vortexState.vortex_list:
    A=np.asarray(v.line)
    lX,lY=A.T
    ax.plot(lX,lY,'.',color='k',ms=.5)
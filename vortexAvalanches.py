#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: spencer
"""
import numpy as np
from scipy import special
import matplotlib

matplotlib.rcParams.update({'font.size': 24})

def dist(p1,p2,circumference): #distance around the cylinder
    dist1 = np.sqrt((p1[0]-p2[0])**2 + ((p1[1]-p2[1]))**2)
    dist2 = np.sqrt((p1[0]-p2[0])**2 + (circumference - np.abs(p1[1]-p2[1]))**2)
    return [dist1,dist2]
    
def rHat(p_i,p_j,circ):#on tube
    distance = dist (p_i,p_j,circ)
    hat1 = np.asarray((p_i-p_j) / distance[0])
    sign = 2*((hat1[1]>0)-1/2)
    hat2 = np.asarray((p_i-(p_j+np.asarray([0,sign*circ])) ) / distance[1])
    return [hat1,hat2]

#2d On a Tube
class SCMaterial2d:

        def __init__(self, X=24, Y=24, site_density=1, pin_radius=.15, maxpinStrength=1, f_not=1, xBias=0): #lamda=1   
            #constant inits
            self.xi = pin_radius
            self.f_not = f_not #Φ_0^2/8π^2λ^3, can leave as 1 for now
            f_p = self.f_not*maxpinStrength
            #sample inits
            self.n = site_density
            self.length = X
            self.height = Y # = circumference of the tube
            
            #defect inits
            numPoints = int(self.length*self.height*self.n)
            randPos = np.random.random((numPoints,2))
            randPos = randPos * [self.length,self.height]
            self.positions = randPos #position array of site defects
            self.pinStrengths = np.random.uniform(f_p/5,f_p,numPoints) #At each point is assigned a pin strength
            
            #Bias
            self.bias = xBias *np.array([1,0])
            
class vortex:
    
    def __init__(self, x,y,velocity=[0,0]):
        self.x = x
        self.y = y
        self.velocity=np.asarray(velocity)
        
        self.line = [[x,y]] #Vortex travel history
      
class vList:
    
    def __init__(self, vortex_list):
        self.vortex_list = vortex_list #vortices are added or deleted as they enter/leave the sample  
        self.line_history = [] #Saves historical deleted lines
        
    def updateForce(self,Sample,dt=1, eta=1, reAdd=True, save_condition=10):
        pin_pos = Sample.positions
        pin_Strength = Sample.pinStrengths
        pinNumber = len(pin_pos)
    
        Xi_p = Sample.xi
        f_not = Sample.f_not
        
        tubeCirc = Sample.height
        tubeLength = Sample.length

        velocityArr=[]
        for vortex_i in self.vortex_list:
            v_i = np.array([vortex_i.x,vortex_i.y])
            currentForce = 0
            
            #force from other vortices,works
            for vortex_j in self.vortex_list:
                if(vortex_i!=vortex_j):
                    v_j = np.array([vortex_j.x,vortex_j.y])
                    dist_ij = dist(v_i,v_j, tubeCirc)
                    dHat = rHat(v_i,v_j, tubeCirc)
                    f_ij = special.kn(1,dist_ij[0])*dHat[0] + special.kn(1,dist_ij[1])*dHat[1]
                    currentForce += f_ij
            currentForce = currentForce*f_not
            
            #force from pinning sites, maybe works?
            for k in range(pinNumber):
                dist_ik = dist(v_i, pin_pos[k], tubeCirc)
                f_p = -1*pin_Strength[k]
                #-----
                minArg = np.argmin(dist_ik)
                minDist = dist_ik[minArg]
                check = (Xi_p>minDist)
                #check = np.heaviside((Xi_p-dist_ik), 1) #dont count a defect site at which a vortice is stuck 
                currentForce += (f_p/Xi_p) * minDist * check * rHat(v_i, pin_pos[k],tubeCirc)[minArg]
            
            #force from the left, may want to change its form
            #currentForce += Sample.bias*(1/(vortex_i.x**2+1))
            
            velocityArr = velocityArr+[currentForce/eta]
        #Update position based on last velocity, then update velocity
        numVortices = len(self.vortex_list) 
        #keep=[]
        i = 0
        while (i < numVortices):
            dx, dy = self.vortex_list[i].velocity*dt
            self.vortex_list[i].x += dx
            self.vortex_list[i].y = (self.vortex_list[i].y + dy)%tubeCirc
            self.vortex_list[i].velocity = np.asarray(velocityArr[i])
            self.vortex_list[i].line = self.vortex_list[i].line + [[self.vortex_list[i].x,self.vortex_list[i].y]]
            
            #Delete Vortices that fall off the edge, re-add them to the left edge
            if (self.vortex_list[i].x<0 or self.vortex_list[i].x>tubeLength):
                curLine = self.vortex_list[i].line
                if (len(curLine)>save_condition):self.line_history = self.line_history + [curLine]
                
                if(reAdd):
                    self.vortex_list[i] = vortex((1/50)*Sample.length * np.random.random(),
                                    Sample.height * np.random.random() )
                    i+=1
                else:
                    del(self.vortex_list[i])
                    numVortices -= 1
            else: i+=1


    def addVortex(self,num,Sample):#Adds vortices randomly along the left edge of the sample
        for i in range(num):
            randomY = Sample.height * np.random.random() 
            randomX = (1/50)*Sample.length * np.random.random()
            newVortex = vortex(randomX,randomY)
            self.vortex_list.append(newVortex)
        
    def getAllLines(self):
        cur_Vortices = self.vortex_list
        current_lines = []
        for v in cur_Vortices:
            current_lines = current_lines + [v.line]
        return current_lines, self.line_history

def Saturate(Material,percent = 1):
    vortList = []
    rad = Material.xi/3
    for r in Material.positions:
        if (np.random.rand() < percent):
            newVortex = vortex(r[0]+np.random.rand()*rad-rad/2,r[1]+np.random.rand()*rad-rad/2)
            vortList = vortList + [newVortex]
    return vList(vortList)

def vlist_to_plot(vList):
    X, Y = [],[]
    for vortex in vList.vortex_list:
        X = X + [vortex.x]
        Y = Y + [vortex.y]
    return X, Y

pass 
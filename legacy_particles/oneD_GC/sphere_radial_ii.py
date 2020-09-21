# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 13:22:33 2015

@author: cra
"""

from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import particle_helpers as ph

from scipy.optimize import minimize_scalar


# solve sphere problem in 1D to relate surface potential 
# to surface charge

n0 = 3.
er = 2.

A1 = 5.
I1 = 7.8747 # scaled surface charge

R_balloon = 3. # 
R_mesh = 0.05 # mesh size as a function of debye length

k0 = sqrt(n0/er)

a = A1/(sqrt(2)*k0)
sigma = er*k0*sqrt(2)*I1

params={'a':a,'sigma':sigma,'n0':n0,'er':er}

print ph.calc_sph_surf_pot(params,doPlot=True)

#--- search opt
def f(s):
    params['sigma']=er*k0*sqrt(2)*s
    return (ph.calc_sph_surf_pot(params,doPlot=False)-2)**2
    
res = minimize_scalar(f)    

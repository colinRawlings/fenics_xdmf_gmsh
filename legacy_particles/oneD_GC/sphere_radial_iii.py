# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 13:22:33 2015

@author: cra
"""

from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import particle_helpers as ph



reload(ph)

# solve sphere problem in 1D to relate surface potential 
# to surface charge

n0 = 3.
er = 2.

A1 = 5.
k0 = sqrt(n0/er)

a = A1/(sqrt(2)*k0)

params={'a':a,'n0':n0,'er':er}

sigma = ph.calc_sph_charge(3,params,log=False)
print sigma/(er*k0*sqrt(2))

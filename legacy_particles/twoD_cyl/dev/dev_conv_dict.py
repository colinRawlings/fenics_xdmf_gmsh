# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 15:23:48 2016

@author: cra

practice dict conversion to dolfin.Constant

"""

import particle_helpers as ph
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import fenics_helpers as fh
import matplotlib
import pickle
import os
import time
import os.path
import subprocess as sb
import dolfin as fn



Npts = 5

lD_nm = 9.6 # Debye length

Tdeg = 21

w_nm = 200 # field size
a_small = 1 # 
Rs_mesh = 0.3 # cell size as a function of 1/k0 at the charge sheets


pH_r = np.linspace(3,11,50)
psi_mV={'Si': 10,'sph': 1e-6,'poly': 10}
geom_nm={'w':w_nm,'L':10*lD_nm,'h':5*lD_nm,'a':a_small,'sb':60,'st':60,'dsr':Rs_mesh*lD_nm,'dss':a_small} # cylindrical symmetry
er_dry = {'sph': 1, 'poly':1,'Si': 1} # n.b. fluid dielectric constant will be calc'd by conv_to_ratchet_ND_pars
const={'e': 1.602e-19,'e0': 8.85e-12,'kB': 1.38e-23,\
           'Na': 6.022e23}


geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)


x_ND = np.linspace(0,geom_ND['L'],Npts)


#--- surface models
surf={'pK': 7.5, 'Gamma': 9e18,'C_stern': 2.5} # Silica
surf_ND = ph.calc_surf_ND_pars(surf,geom_ND,params_ND)


def dict_vals_to_Const(a_dict):
    ''' recursively convert the val of all keys in a dict which are floats
    or ints into fn.Constant (s)
    '''
    
    new_dict = {}
    for key,val in a_dict.iteritems():
        if (isinstance(val,float) or isinstance(val,int)):
            new_dict[key] = fn.Constant(val)
        elif isinstance(val,dict):
            new_dict[key] = dict_vals_to_Const(val)
        else:
            new_dict[key]=val
                       
    return new_dict



params_NDd = dict_vals_to_Const(params_ND)

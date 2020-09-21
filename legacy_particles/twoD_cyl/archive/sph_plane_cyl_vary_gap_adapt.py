# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 18:16:20 2015

@author: cra

calculates F for a range of particle positions in the gap


"""

from dolfin import *
import particle_helpers as ph
import numpy as np
import matplotlib.pyplot as plt
import fenics_helpers as fh
import matplotlib
import pickle
import dolfin as fn

reload(ph)
reload(fh)

def x2h(x,geom_nm):
    return x+geom_nm['L']*0.5-geom_nm['a']


plt.close('all')
fh.run_cmd([r'killall',r'paraview'])

parameters.form_compiler.quadrature_degree = 4
parameters.allow_extrapolation = True
parameters["refinement_algorithm"] = "plaza_with_parent_facets"

##============== variables
fName_tail = r'/home/cra/dolfin/particles/twoD_cyl/mesh/cyl_sph_plane_tail_02.geo'
fName_log = r'/home/cra/dolfin/particles/twoD_cyl/logFile.txt'
folderResults = r'/home/cra/dolfin/particles/twoD_cyl/'

Rs_mesh = 1 # cell size as a function of 1/k0 at the charge sheets
Nh = 20 # number of different particle positions to try

R_h_min = 0.2 # min separation between surface and plane

#Lr_nm = np.logspace(np.log10(100),np.log10(300),5)
Lr_nm = np.linspace(80,600,10)


lD_nm = 30 # Debye length
Tdeg = 21
w = 200 # field size
a_nm = 30.

psi_mV={'Si':67,'sph':58,'poly':30}
geom_nm={'w':w,'L':160,'h':50,'a':a_nm,'sb':60,'st':30,'dsr':Rs_mesh*lD_nm,'dss':a_nm/5} # cylindrical symmetry
er_dry = {'sph': 1, 'poly':1,'Si': 1} # n.b. fluid dielectric constant will be calc'd by conv_to_ratchet_ND_pars


#-- object numbering
SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':1,'Si':2}}


#--- required variables for calcs
Fr_int_sph_pl = np.zeros(Nh)
Fr_int_pl_pl = np.zeros(len(Lr_nm))
psi_mV_pl_pl = psi_mV.copy()
psi_mV_pl_pl['sph'] = 1e-6 # switch off sph
psi_mV_sph_pl = psi_mV.copy()
psi_mV_sph_pl['poly'] = 1e-6 # switch off polymer charge
psi_mV_sph = psi_mV_sph_pl.copy()
psi_mV_sph['Si'] = 1e-6 # leave only sphere charged


geom_nm['L'] = np.max(Lr_nm+2*lD_nm)
geom_nm_pl_pl = geom_nm.copy()
geom_nm_pl_pl['a'] = 1 # make very small!

geom_ND,params_ND_sph_pl=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV_sph_pl)
geom_ND_pl_pl,params_ND_pl_pl=ph.conv_to_ratchet_ND_pars(geom_nm_pl_pl,lD_nm,er_dry,Tdeg,psi_mV_pl_pl)
geom_ND,params_ND_sph=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV_sph)
geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)

hr_sph_pl_nm = np.linspace(lD_nm*R_h_min,np.max(Lr_nm)-2*geom_nm['a']-lD_nm*R_h_min,Nh)


#================ calcs
q = len(Lr_nm)-1
N = 20

geom_nm['L'] = Lr_nm[q]
geom_ND['L'] = Lr_nm[q]*1e-9/geom_ND['l0_m']   
geom_ND['h'] = x2h(0,geom_nm)*1e-9/geom_ND['l0_m']    

Fr_ad = np.zeros(N)
f1 = np.zeros(N)

#for n in range(N):
    
Fr_ad[n],Fc_int,u_,mesh = ph.solve_ratchet_problem(geom_ND,params_ND,SDs,fName_tail,\
   doPlots=True,cylSim=True,N_adapt = N,adapt_TOL=0.01)
f = fn.dot(fn.grad(u_.leaf_node()),fn.grad(u_.leaf_node()))*fn.dx

print('n: %d, F: %.5f'%(n,Fr_ad[n]))


#for q in range(len(Lr_nm)):
#    geom_ND['L'] = Lr_nm[q]*1e-9/geom_ND['l0_m']   
#    h_min = lD_nm*R_h_min # don't let the sphere actually touch the bounding wall
#    hr_nm.append(np.linspace(h_min,Lr_nm[q]-2*geom_nm['a']-h_min,Nh))    
#    
#    for p in range(len(hr_nm[q])):
#        showPlots = (Lr_nm[q]==np.max(Lr_nm) and hr_nm[q][p]==np.max(hr_nm[q]))\
#                    or (Lr_nm[q]==np.min(Lr_nm) and hr_nm[q][p]==np.max(hr_nm[q]))
#                    
#        print('Progress: %d/%d'%(p+q*Nh,Nh*len(Lr_nm)))
#        geom_ND['h'] = hr_nm[q][p]*1e-9/geom_ND['l0_m']    
#        Fr_int_gap[q,p],Fc_int,u_ = ph.solve_ratchet_problem(geom_ND,params_ND,SDs,fName_tail,doPlots=showPlots,cylSim=True)


###============== analysis
## see anl_`this script`
#fName_save = folderResults+r'/vary_L_h.pickle'
#
#R = dict()
#for fld in ['geom_ND','params_ND','geom_nm','Lr_nm','hr_nm','Fr_int_gap',\
#            'Fr_int_sph_pl','Fr_sph','hr_sph_pl_nm','Fr_int_pl_pl']:
#    R[fld]=eval(fld)
#    
#with open(fName_save, 'w') as f:
#    pickle.dump(R, f)    
            
            







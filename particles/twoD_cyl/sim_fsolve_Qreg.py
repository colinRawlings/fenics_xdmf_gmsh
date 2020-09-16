# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 14:28:53 2016

@author: cra

solve several Grahame problems in parallel

"""
from six.moves import range

exec(compile(open('/home/cra/dolfin/useful_f/init_cra.py').read(), '/home/cra/dolfin/useful_f/init_cra.py', 'exec'))

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

def calc_err(x_ND,sigma_ND,pH,surf):
    ''' for now assume Grahame holds locally'''    
    
    err = list()
    for p in range(len(x_ND)):
        err.append(ph.calc_phi_d_Grahame(sigma_ND[p],params_ND)-ph.calc_phi_d_SiOH_Beh(sigma_ND[p],pH,surf))
    
    return tuple(err)
    
#--- solve
s_ND = np.zeros((len(pH_r),len(x_ND)))    
    
for p in range(len(pH_r)):
       
       s0 = -1*np.ones(np.shape(x_ND))
       RES = sp.optimize.fsolve(lambda s: calc_err(x_ND,s,pH_r[p],surf_ND),s0,\
           full_output = True,factor=1e-2)    
       
       s_ND[p,:] = RES[0]


f,ax = plt.subplots(nrows=1,ncols=1)

for p in range(len(x_ND)):
    col = fh.chot(p,len(x_ND)) 
    ax.plot(pH_r,s_ND[:,p],'-',color=col)







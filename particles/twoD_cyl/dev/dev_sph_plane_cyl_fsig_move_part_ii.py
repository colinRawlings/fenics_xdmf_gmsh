# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 18:16:20 2015

@author: cra

calculates F for a range of particle positions in the gap

the charge on the gold particle is fixed while the silicon charge follows the 
model of Behrens

"""

#--- to help rope
import numpy as np
import particle_helpers as ph
import fenics_helpers as fh


execfile('/home/cra/dolfin/useful_f/init_cra.py')

fName_tail = r'/home/cra/dolfin/particles/twoD_cyl/mesh/cyl_sph_plane_tail_02.geo'
fName_log = r'/home/cra/dolfin/particles/twoD_cyl/logFile.txt'
folderResults = r'/home/cra/dolfin/particles/twoD_cyl/'

Rs_mesh = 0.2 # cell size as a function of 1/k0 at the charge sheets

R_h_min = 0.1 # min separation between surface and plane
Nh = 10
NL = 3

Tdeg = 21
lD_nm = ph.n0_M_to_lD_nm(0.001,Tdeg) # Debye length
w = 10*lD_nm # field size
a_nm = 30.
Lr_nm = np.linspace(5*lD_nm+2*a_nm,2*R_h_min*lD_nm+lD_nm+2*a_nm,NL)

psi_mV={'Si':67,'sph':-60,'poly':-60}
geom_nm={'w':w,'L':10*lD_nm,'h':R_h_min*lD_nm,'a':a_nm,'sb':2*lD_nm,'st':2*lD_nm,\
        'dsr':Rs_mesh*lD_nm,'dss': a_nm/10} # cylindrical symmetry
er_dry = {'sph': 1e3, 'poly':1,'Si': 1} # n.b. fluid dielectric constant will be calc'd by conv_to_ratchet_ND_pars

surf={'pK': 7.5, 'Gamma': 9e18,'C_stern': 2.5} 
pH = 6.0

#--- get initial pars
geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)
surf_ND = ph.calc_surf_ND_pars(surf,Tdeg)
const = ph.load_const()
sigma_ND_inf,phi_ND_inf = ph.solve_isol_pl_Qreg(pH,surf_ND,lD_nm,Tdeg)        
params_ND['sigma_ND_inf'] = sigma_ND_inf
params_ND['phi_ND_inf'] = phi_ND_inf

SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':1,'Si':2}}

Nx = 1*geom_nm['w']/lD_nm
x_ND = np.linspace(0,geom_ND['w'],Nx)
    
def calc_error(x_ND,sigma_ND,pH,surf_ND):
    err = ph.calc_phi_d_PB_fsig(x_ND,sigma_ND,geom_ND,params_ND,SDs,fName_tail\
            ,set_sigma_poly=False) \
            - ph.calc_phi_d_SiOH_Beh(sigma_ND,pH,surf_ND)
    return err
    
    
#========= calc    
#--- setup: solve
s0 = sigma_ND_inf*np.ones(len(x_ND))


#sigma_NDr = np.zeros((len(Lr_nm),len(x_ND)))
Fchr = np.zeros((len(Lr_nm),Nh))
Ffr =  np.zeros((len(Lr_nm),Nh))

h_min_nm = lD_nm*R_h_min # don't let the sphere actually touch the bounding wall
hr_nm = list()    
sigma_NDrr=list()

for p in range(len(Lr_nm)):
    geom_ND['L'] = Lr_nm[p]*1e-9/geom_ND['l0_m']
    hr_nm.append(np.linspace(Lr_nm[p]-2*geom_nm['a']-h_min_nm,h_min_nm,Nh)) # start from a large separation where s0 sigma_ND_inf should be a good guess
    
    sigma_NDrr.append( float('nan')*np.ones((Nh,len(x_ND))) )  # init to failure value
    s0 = sigma_ND_inf*np.ones(len(x_ND)) # reset to inf value
    
    for q in range(Nh):
        print('Starting: L: (%d/%d), h: (%d/%d)'%(p+1,len(Lr_nm),q+1,Nh))
        
        geom_ND['h'] = hr_nm[p][q]*1e-9/geom_ND['l0_m']  
        
        x,infodict,ier,mesg = sp.optimize.fsolve(lambda s: calc_error(x_ND,s,pH,surf_ND),s0,\
            full_output = True,factor=1,xtol=1e-2)     
        
        if ier==1:    #success
            sigma_NDrr[p][q,:] = x
            s0 = x
            
            #--- repeat soln to recover F's  
            fh.run_cmd([r'killall',r'paraview'])
            phi_d,Ff,Fch,u_,mesh = ph.calc_phi_d_PB_fsig(x_ND,x,geom_ND,params_ND,\
                SDs,fName_tail,doPlots=True,full_output=True,set_sigma_poly=False)
            Fchr[p,q] = Fch        
            Ffr[p,q] = Ff
        else:
            print('failed')
        
    


#====== plotting
#see anl_${this script}

#======== save
fName_save = folderResults+r'move_part_fsig.pickle'

params_ND['sigma_Si']=float('nan') # can't pickle a dolfin expression

R = dict()
for fld in ['x_ND','geom_ND','params_ND','geom_nm','lD_nm','Lr_nm',\
            'sigma_NDrr','surf_ND','hr_nm','Fchr','Ffr']:
    R[fld]=eval(fld)
    
with open(fName_save, 'w') as f:
    pickle.dump(R, f,protocol=-1)    





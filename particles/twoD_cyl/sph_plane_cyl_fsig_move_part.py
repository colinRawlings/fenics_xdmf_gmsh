# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 18:16:20 2015

@author: cra

calculates F for a range of particle positions in the gap

the charge on the gold particle is fixed while the silicon charge follows the 
model of Behrens

"""

#--- to help rope
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import particle_helpers as ph
import fenics_helpers as fh
from six.moves import range


exec(compile(open('/home/cra/dolfin/useful_f/init_cra.py').read(), '/home/cra/dolfin/useful_f/init_cra.py', 'exec'))

fName_tail = r'/home/cra/dolfin/particles/twoD_cyl/mesh/cyl_sph_plane_tail_02.geo'
fName_log = r'/home/cra/dolfin/particles/twoD_cyl/logFile.txt'
folderResults = r'/home/cra/dolfin/particles/twoD_cyl/'

Rs_mesh = 0.2 # cell size as a function of 1/k0 at the charge sheets

R_h_min = 0.1 # min separation between surface and plane
Rx_xLTa = 0.25 # point spacing a function of the Debye length for x<a (sphere radius)
Rx_xGTa = 1.5 # point spacing a function of the Debye length for x>a (sphere radius)

Nh = 20
NL = 10

Tdeg = 21
lD_nm = ph.n0_M_to_lD_nm(0.0001,Tdeg) # Debye length
w = 10*lD_nm # field size
a_nm = 30.
Lr_nm = np.linspace(5*lD_nm+2*a_nm,2*R_h_min*lD_nm+lD_nm+2*a_nm,NL)
#Lr_nm = np.asarray([2*R_h_min*lD_nm+lD_nm+2*a_nm])

psi_mV={'Si':67,'sph':-60,'poly':-67}
geom_nm={'w':w,'L':10*lD_nm,'h':R_h_min*lD_nm,'a':a_nm,'sb':2*lD_nm,'st':2*lD_nm,\
        'dsr':Rs_mesh*lD_nm,'dss': a_nm/10} # cylindrical symmetry
er_dry = {'sph': 1e3, 'poly':4,'Si': 4} # n.b. fluid dielectric constant will be calc'd by conv_to_ratchet_ND_pars

surf={'pK': 7.5, 'Gamma': 9e18,'C_stern': 2.5} 
pH = 5.815 # 5.815: achieve a surface potential of -67 mV

#--- get initial pars
geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)
surf_ND = ph.calc_surf_ND_pars(surf,Tdeg)
const = ph.load_const()
sigma_ND_inf,phi_ND_inf,phi_d_inf_mV = ph.solve_isol_pl_Qreg(pH,surf_ND,lD_nm,Tdeg)        

print(phi_d_inf_mV)

params_ND['sigma_ND_inf'] = sigma_ND_inf
params_ND['phi_ND_inf'] = phi_ND_inf

SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':1,'Si':2}}

x_ND,Nx = ph.calc_fsig_grid_pts(geom_ND,1.5,Rx_xLTa,Rx_xGTa,smoothing=3e-3,doPlots=False)
        
f,axKvg = fh.cfigw()        
def calc_error(x_ND,sigma_ND,geom_ND,params_ND,SDs,fName_tail,pH,surf_ND,progStr):
    
    err = ph.calc_phi_d_PB_fsig(x_ND,sigma_ND,geom_ND,params_ND,SDs,fName_tail\
            ,set_sigma_poly=False) \
            - ph.calc_phi_d_SiOH_Beh(sigma_ND,pH,surf_ND)

    axKvg[0].cla()
    axKvg[0].plot(x_ND,sigma_ND,'.-k')
    fh.cxlbl(axKvg[0],'x','\sigma')
    fh.ctitle(axKvg[0],progStr)    
    
    axKvg[1].cla()
    axKvg[1].plot(x_ND,err,'.-k')
    fh.cxlbl(axKvg[1],'x','\Delta \phi')
    fh.ctitle(axKvg[1],progStr)

    return err
    
    
#========= calc  
    
#--- init 
F={}    
for sDef in ['s_const','s_reg']:
    F[sDef] = dict()    
    for term in ['surf','vol']:
        F[sDef][term] = np.zeros((len(Lr_nm),Nh))
        
h_min_nm = lD_nm*R_h_min # don't let the sphere actually touch the bounding wall
hr_nm = list()    
sigma_NDrr=list()

#--- the calc
for p in range(len(Lr_nm)):
    geom_ND['L'] = Lr_nm[p]*1e-9/geom_ND['l0_m']
    hr_nm.append(np.linspace(Lr_nm[p]-2*geom_nm['a']-h_min_nm,h_min_nm,Nh)) # start from a large separation where s0 sigma_ND_inf should be a good guess
    
    sigma_NDrr.append( float('nan')*np.ones((Nh,len(x_ND))) )  # init to failure value
    s0 = sigma_ND_inf*np.ones(len(x_ND)) # reset to inf value
    
    for q in range(Nh):
        progStr = 'L: (%d/%d), h: (%d/%d)'%(p+1,len(Lr_nm),q+1,Nh)
        
        geom_ND['h'] = hr_nm[p][q]*1e-9/geom_ND['l0_m']  
        s_it_num = 0        
        
        x,infodict,ier,mesg = sp.optimize.fsolve(lambda s: calc_error(x_ND,s,geom_ND,\
            params_ND,SDs,fName_tail,pH,surf_ND,progStr),s0,full_output = True,factor=1,xtol=1e-2)     
        
        if ier==1:    #success
            sigma_NDrr[p][q,:] = x
            s0 = x
            
            #--- repeat soln to recover F's  
            res = fh.run_cmd([r'killall',r'paraview'])
            phi_d,F['s_reg']['vol'][p,q],F['s_reg']['surf'][p,q],u_,mesh = ph.calc_phi_d_PB_fsig(x_ND,x,geom_ND,params_ND,\
                SDs,fName_tail,doPlots=False,full_output=True,set_sigma_poly=False,\
                N_adapt=6,adapt_TOL=3e-3,disp_log=True)
            
            phi_d,F['s_const']['vol'][p,q],F['s_const']['surf'][p,q],u_,mesh = ph.calc_phi_d_PB_fsig(x_ND,sigma_ND_inf*np.ones(len(x_ND)),\
                geom_ND,params_ND,SDs,fName_tail,doPlots=False,full_output=True,\
                set_sigma_poly=False,N_adapt=6,adapt_TOL=3e-3,disp_log=True)
            
        else:
            print('failed')
        
        

#====== plotting
#see anl_${this script}

#======== save
fName_save = folderResults+r'move_part_fsig'+fh.tstamp()+'.pickle'

params_ND['sigma_Si']=float('nan') # can't pickle a dolfin expression

R = dict()
for fld in ['x_ND','geom_ND','params_ND','geom_nm','lD_nm','Lr_nm',\
            'sigma_NDrr','surf_ND','hr_nm','Tdeg','pH','F',\
            'sigma_ND_inf','phi_ND_inf','phi_d_inf_mV']:
    R[fld]=eval(fld)
    
with open(fName_save, 'w') as f:
    pickle.dump(R, f,protocol=-1)    

print(('results written to: '+fName_save))




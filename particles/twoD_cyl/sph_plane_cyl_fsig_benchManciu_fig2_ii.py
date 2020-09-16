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
figFolder = r'/home/cra/Documents/pres/sph_cyl/figs/'

Rs_mesh = 0.3 # cell size as a function of 1/k0 at the charge sheets

R_h_min = 0.1 # min separation between surface and plane
Nh = 20
NL = 8

Tdeg = 21
lD_nm = ph.n0_M_to_lD_nm(0.01,Tdeg) # Manciu Debye length
w = 10*lD_nm # field size
a_small = 0.01*lD_nm
Lr_nm = np.linspace(5*lD_nm,2,8)

psi_mV={'Si':67,'sph':1e-6,'poly':-1e-6}
geom_nm={'w':w,'L':10*lD_nm,'h':0.1*lD_nm,'a':a_small,'sb':2*lD_nm,'st':2*lD_nm,\
        'dsr':Rs_mesh*lD_nm,'dss': a_small} # cylindrical symmetry
er_dry = {'sph': 80, 'poly':1,'Si': 1} # n.b. fluid dielectric constant will be calc'd by conv_to_ratchet_ND_pars

#-- Manciu coefficients
surf={'pK': 7.5, 'Gamma': 3.125e18,'C_stern': 1e3} 
pH = 6.5 

#--- get initial pars
geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)
surf_ND = ph.calc_surf_ND_pars(surf,Tdeg)
const = ph.load_const()
sigma_ND_inf,phi_ND_inf,phi_d_mV = ph.solve_isol_pl_Qreg(pH,surf_ND,lD_nm,Tdeg)        
params_ND['sigma_ND_inf'] = sigma_ND_inf
params_ND['phi_ND_inf'] = phi_ND_inf



#-- object numbering
SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':1,'Si':2}}


Nx = 1*geom_nm['L']/lD_nm

#--- construct spline rep
x_ND = np.linspace(0,geom_ND['w'],Nx)
xc = np.linspace(0,geom_ND['w'],Nx*10)
sigma_ND = 10*np.sin(x_ND*10)

#--- Manciu soln
scaleFac = 8e-4/50.5 # (J/m**2)/mm
Lr_man_nm = np.asarray([1,2,4,6,8])
F_man_ND =  (np.asarray([32,23.5,15,10.5,8])-5.5)*scaleFac*(np.pi*(geom_nm['w']*1e-9)**2)/(const['kB']*(Tdeg+273))



    
def calc_error(x_ND,sigma_ND,pH,surf_ND):
    err = ph.calc_phi_d_PB_fsig(x_ND,sigma_ND,geom_ND,params_ND,SDs,fName_tail\
            ,set_sigma_poly=True) \
            - ph.calc_phi_d_SiOH_Beh(sigma_ND,pH,surf_ND)
    print(('RMS err (%d/%d): %.5e'%(p+1,len(Lr_nm),np.dot(err,err))))
    return err
    
    
#========= calc    
s0 = sigma_ND_inf*np.ones(len(x_ND))


sigma_NDr = np.zeros((len(Lr_nm),len(x_ND)))
Fchr = np.zeros(len(Lr_nm))
Ffr =  np.zeros(len(Lr_nm))

for p in range(len(Lr_nm)):
    
    geom_ND['L'] = Lr_nm[p]*1e-9/geom_ND['l0_m']
    x,infodict,ier,mesg = sp.optimize.fsolve(lambda s: calc_error(x_ND,s,pH,surf_ND),s0,\
        full_output = True,factor=1,xtol=1e-2)     
    
    sigma_NDr[p,:] = float('nan')*np.ones(np.shape(x_ND))
    if ier==1:    
        sigma_NDr[p,:] = x
        s0 = x
        
        #--- repeat soln to calc F     
        phi_d,Ff,Fch,u_,mesh = ph.calc_phi_d_PB_fsig(x_ND,sigma_NDr[p,:]\
            ,geom_ND,params_ND,SDs,fName_tail,full_output=True,set_sigma_poly=True)
        Fchr[p] = Fch        
        Ffr[p] = Ff
    else:
        print('failed')
        
    


#====== plotting
f0,axUP = plt.subplots(nrows=1,ncols=2)
fh.make_cfigw()
fh.place_fig(f0)

x_lD = x_ND*(geom_ND['l0_m']/(lD_nm*1e-9))
Lr_lD = Lr_nm/lD_nm
xc_lD = xc*(geom_ND['l0_m']/(lD_nm*1e-9))

Iinf = np.argmax(Lr_nm)
F_inf = Fchr[Iinf]+Ffr[Iinf]


for p in range(len(Lr_nm)):
    sci = fh.spline_con_gradient(x_ND,sigma_NDr[p,:],[],[],k=3)    
    
    axUP[0].plot(x_lD,sigma_NDr[p,:]/surf_ND['Gamma'],'+',color=fh.chot(p,len(Lr_nm)))
    axUP[0].plot(xc_lD,sci(xc)/surf_ND['Gamma'],'-',\
        color=fh.chot(p,len(Lr_nm)),label='L: %.1f$\lambda_D$'%(Lr_nm[p]/lD_nm))
#
axUP[0].plot(x_lD,sigma_ND_inf*np.ones(len(x_ND))/surf_ND['Gamma']\
            ,'--k',label=r'$L \rightarrow \infty$')
axUP[0].legend()
axUP[0].set_xlabel('$x/\lambda_D$')
axUP[0].set_ylabel('$\sigma / (e \Gamma)$')    
    
#-- cf at edge
axUP[1].plot(Lr_lD,2*Fchr+Ffr-F_inf,'ok',label='fenics') # mult by 2 as two surfaces and integral only conducted over one
axUP[1].plot(Lr_man_nm/lD_nm,F_man_ND,'+-r',label='Manciu')
axUP[1].set_xlabel('$L/\lambda_D$')
axUP[1].set_ylabel('$(F(L) - F(\infty))/ k_B T$')
axUP[1].legend()    
    
plt.draw()

f,axU1 = fh.pfig()

axU1.plot(Lr_lD,2*Fchr+Ffr-F_inf,'ok',label='fenics') # mult by 2 as two surfaces and integral only conducted over one
axU1.plot(Lr_man_nm/lD_nm,F_man_ND,'+-r',label='Manciu')
axU1.legend(frameon=False,loc='best')  
fh.cxlbl(axU1,'L/\lambda_D','(F(L) - F(\infty))/ k_B T')
f.savefig(figFolder+'bench_Manciu_F.pdf')


##======== save
#fName_save = folderResults+r'move_part_in_Qreg_gap.pickle'
#
#params_ND['sigma_Si']=float('nan') # can't pickle a dolfin expression
#
#R = dict()
#for fld in ['geom_ND','params_ND','geom_nm','lD_nm','Lr_nm',\
#            'sigma_NDr','surf_ND','hr_nm','Fchr','Ffr']:
#    R[fld]=eval(fld)
#    
#with open(fName_save, 'w') as f:
#    pickle.dump(R, f,protocol=-1)    





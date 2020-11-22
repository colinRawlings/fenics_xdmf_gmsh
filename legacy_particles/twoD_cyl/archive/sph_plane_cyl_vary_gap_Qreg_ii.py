# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 18:16:20 2015

@author: cra

calculates F for a range of particle positions in the gap


"""

execfile('/home/cra/dolfin/useful_f/init_cra.py')


##============== variables
fName_tail = r'/home/cra/dolfin/particles/twoD_cyl/mesh/cyl_sph_plane_tail_02.geo'
fName_log = r'/home/cra/dolfin/particles/twoD_cyl/logFile.txt'
folderResults = r'/home/cra/dolfin/particles/twoD_cyl/'

Rs_mesh = 0.3 # cell size as a function of 1/k0 at the charge sheets
Nh = 20 # number of different particle positions to try
N_adapt_max = 20
adapt_TOL = 0.001


pH = 10
lD_nm = 9.6 # Debye length

Tdeg = 21
w_nm = 200 # field size
a_small = 1 # 

psi_mVr = np.linspace(10,90,5)

psi_mV={'Si':psi_mVr[0],'sph':1e-6,'poly':psi_mVr[0]}
geom_nm={'w':w_nm,'L':10*lD_nm,'h':5*lD_nm,'a':a_small,'sb':60,'st':60,'dsr':Rs_mesh*lD_nm,'dss':a_small} # cylindrical symmetry
er_dry = {'sph': 1, 'poly':1,'Si': 1} # n.b. fluid dielectric constant will be calc'd by conv_to_ratchet_ND_pars
const={'e': 1.602e-19,'e0': 8.85e-12,'kB': 1.38e-23,\
           'Na': 6.022e23}

#--- surface models
#surf={'pK': 7.5,'alpha': 1, 'Gamma': 9e18,'C_stern': 2.5} # Silica
surf={'pK': 4.9,'alpha': 1, 'Gamma': 0.25e18,'C_stern': 100} # Carboxyl


SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':1,'Si':2}} # domain numbering

Tkelv = Tdeg+273
geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)

#--- calc ND surf pars
surf_ND = dict()
for fld in ['pK','alpha']:
    surf_ND[fld] = surf[fld]
surf_ND['Gamma']= surf['Gamma']*geom_ND['l0_m']**2
surf_ND['C_stern'] = surf['C_stern']*geom_ND['l0_m']**2*const['kB']*Tkelv/const['e']**2

#--- required variables for calcs
Fr_int = np.zeros(len(psi_mVr))
sigma_r_ND = np.zeros(len(psi_mVr))
phi_r_ND = np.zeros(len(psi_mVr))


def calc_phi_d_PB(sigma_ND):
    ''' 
    calc the relationship between the Stern potential and the surface charge dictated
    by the PB relationship
    '''
    
    params_ND['sigma_Si'] = sigma_ND
    params_ND['sigma_poly'] = sigma_ND
    
    Fr_int[p],Fc_int,u_,mesh = ph.solve_ratchet_problem(geom_ND,params_ND,\
        SDs,fName_tail,doPlots=False,cylSim=True,N_adapt = 0,adapt_TOL=adapt_TOL)
        
    return u_(0,geom_ND['L'])
    
def calc_phi_d_Grahame(sigma_ND):
    '''
        Invert the Grahame equation (ND) to obtain the solution for flat 
        isolated surfaces
    '''    
    return 2*np.arcsinh(sigma_ND/np.sqrt(8*params_ND['n0']*params_ND['er_fluid']))
    
    
def calc_phi_d_SiOH_Beh(sigma_ND,pH):
    ''' 
    calc the relationship between the Stern potential and the surface charge dictated
    by the reaction coefficients for the dissociation of the Silanol groups
    '''
    
    #--- check no positive values for sigma_ND
    if np.sum(np.asarray(sigma_ND>0,dtype=int))> 0:
        raise ValueError('sigma_ND cannot exceed 0')
    
    return -sigma_ND/surf_ND['C_stern']+(surf_ND['pK']-pH)*np.log(10)\
            + np.log(-sigma_ND/(surf_ND['Gamma']+sigma_ND))+np.log(surf_ND['alpha'])
    
def calc_err(sigma_ND,pH):
    
    err = calc_phi_d_Grahame(sigma_ND)-calc_phi_d_SiOH_Beh(sigma_ND,pH)
    
    return err

#--- calc fig 2
pH_r = np.linspace(3,11,1e2)
sigma_ND_r = np.zeros(np.shape(pH_r))

for p in range(len(pH_r)):
    RES = sp.optimize.minimize(lambda s: calc_err(s,pH_r[p])**2,-1,\
                bounds=[-surf_ND['Gamma'],0])
    
    sigma_ND_r[p] = float('nan')    
    if RES['success']:
        sigma_ND_r[p] = RES['x'][0]
    else:
        print('calc for pH: %.1f failed'%(pH_r[p]))        

##--- display
sri = np.linspace(-0.5,-surf_ND['Gamma']*0.9,1e3)
phi_d_Gi = calc_phi_d_Grahame(sri)
phi_d_Beh = calc_phi_d_SiOH_Beh(sri,np.max(pH_r))

f0,ax = plt.subplots(nrows=1,ncols=2)
fh.make_cfigw()

ax[0].plot(phi_d_Gi,-sri/surf_ND['Gamma'],'-r',label='Grahame')
ax[0].plot(phi_d_Beh,-sri/surf_ND['Gamma'],'-k',label='Behrens')
ax[0].legend()
ax[0].set_xlabel('$\phi$')
ax[0].set_ylabel('-$\sigma/ \Gamma$')
ax[0].set_title('pH: %.1f'%(np.max(pH_r)))


ax[1].plot(pH_r,-sigma_ND_r/surf_ND['Gamma'],'-k')
ax[1].set_xlabel('pH')
ax[1].set_ylabel('-$\sigma/ \Gamma$')
fh.place_fig(f0)




fh.place_fig(f0)
    








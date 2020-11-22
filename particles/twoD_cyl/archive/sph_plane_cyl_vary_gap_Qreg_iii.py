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


#surf={'pK': 7.5, 'Gamma': 9e18,'C_stern': 2.5}
surf={'pK': 7.5, 'Gamma': 9e18,'C_stern': 2.5}
#surf={'pK': 4.9, 'Gamma': 0.25e18,'C_stern': 100}


SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':1,'Si':2}} # domain numbering

Tkelv = Tdeg+273
geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)

#--- calc ND surf pars
surf_ND['pK'] = surf['pK']    
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
    
    #    print(sigma_ND)    
    
    #--- check no positive values for sigma_ND
    if np.sum(np.asarray(sigma_ND>0,dtype=int))> 0:
        raise ValueError('sigma_ND cannot exceed 0')
    
    return -sigma_ND/surf_ND['C_stern']+(surf_ND['pK']-pH)*np.log(10)\
            + np.log(-sigma_ND/(surf_ND['Gamma']+sigma_ND))
    
def calc_err(sigma_ND,pH):
    
    err = calc_phi_d_Grahame(sigma_ND)-calc_phi_d_SiOH_Beh(sigma_ND,pH)
    
    return err


#--- pre-calc init range

def calc_sigma_SiOH_Beh(phi_d,pH):
    return sp.optimize.fsolve(lambda s: calc_phi_d_SiOH_Beh(s,pH)-phi_d,-1,factor=1e-3)
    
#--- test inverse
phi_dr = np.linspace(0,-25)
sigma_Beh_r = np.zeros(np.shape(phi_dr))

f0,ax = plt.subplots(nrows=1,ncols=2)
fh.make_cfigw()

for p in range(len(phi_dr)):
    sigma_Beh_r[p] = calc_sigma_SiOH_Beh(phi_dr[p],6)
    ax[0].plot(phi_dr[p],sigma_Beh_r[p],'.k')


#print(calc_sigma_SiOH_Beh(-1,6))

    

##--- calc fig 2
#pH_r = np.linspace(3,11,1e3)
#sigma_ND_r = np.zeros(np.shape(pH_r))
#
#for p in range(len(pH_r)):
#    RES = sp.optimize.fsolve(calc_err,-1,args=(pH_r[p],),\
#                    factor=1e-3,full_output=True)
#    sigma_ND_r[p] = RES[0][0]
#
###--- display
##f0,ax = plt.subplots(nrows=1,ncols=2)
##fh.make_cfigw()
##
#
#
#
#sri = np.linspace(-0.5,-150e3,1e3)
#phi_d_Gi = calc_phi_d_Grahame(sri)
#phi_d_Beh = calc_phi_d_SiOH_Beh(sri,pH)
#
#f0,ax = plt.subplots(nrows=1,ncols=2)
#fh.make_cfigw()
#
#ax[0].plot(sri,phi_d_Gi,'-r',label='Grahame')
#ax[0].plot(sri,phi_d_Beh,'-k',label='Behrens')
#ax[0].legend()
#ax[0].set_xlabel('$\sigma$')
#ax[0].set_ylabel('$\phi$')
#
#
#ax[1].plot(pH_r,sigma_ND_r,'-k')
#ax[1].set_xlabel('pH')
#ax[1].set_ylabel('$\sigma$')
#fh.place_fig(f0)
#
#
#
#
#fh.place_fig(f0)
    


#for p in range(len(psi_mVr)):
#    for fld in ['Si','poly']:
#        psi_mV[fld] = psi_mVr[p]
#
#    geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)
#    Fr_int[p],Fc_int,u_,mesh = ph.solve_ratchet_problem(geom_ND,params_ND,\
#        SDs,fName_tail,doPlots=False,cylSim=True,N_adapt = 0,adapt_TOL=adapt_TOL)
#        
#    #--- extract data
#    phi_r_ND[p] = u_(0,geom_ND['L'])
#    sigma_r_ND[p] = params_ND['sigma_Si']


##--- plot
#f0,ax = plt.subplots(nrows=1,ncols=2)
#
#psi_NDi = np.linspace(0,np.max(psi_mVr)*1e-3*const['e']/(const['kB']*Tkelv))
#
#
#
#ax[0].plot(phi_r_ND,sigma_r_ND,'ok')
#ax[0].plot(psi_NDi,sigma_G(psi_NDi),'-k')
#ax[0].set_xlabel('$\phi$')
#ax[0].set_ylabel('$\sigma$')
#
#fh.place_fig(f0)









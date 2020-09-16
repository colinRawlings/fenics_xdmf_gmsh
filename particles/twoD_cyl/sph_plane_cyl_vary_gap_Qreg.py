# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 18:16:20 2015

@author: cra

calculates F for a range of particle positions in the gap


"""
from six.moves import range

exec(compile(open('/home/cra/dolfin/useful_f/init_cra.py').read(), '/home/cra/dolfin/useful_f/init_cra.py', 'exec'))


##============== variables
fName_tail = r'/home/cra/dolfin/particles/twoD_cyl/mesh/cyl_sph_plane_tail_02.geo'
fName_log = r'/home/cra/dolfin/particles/twoD_cyl/logFile.txt'
folderResults = r'/home/cra/dolfin/particles/twoD_cyl/'

Rs_mesh = 0.3 # cell size as a function of 1/k0 at the charge sheets
Nh = 20 # number of different particle positions to try
N_adapt_max = 20
adapt_TOL = 0.001


lD_nm = 30 # Debye length
Tdeg = 21
w_nm = 200 # field size
a_small = 1 # 

psi_mVr = np.linspace(10,90,5)

psi_mV={'Si':psi_mVr[0],'sph':1e-6,'poly':psi_mVr[0]}
geom_nm={'w':w_nm,'L':10*lD_nm,'h':5*lD_nm,'a':a_small,'sb':60,'st':60,'dsr':Rs_mesh*lD_nm,'dss':a_small} # cylindrical symmetry
er_dry = {'sph': 1, 'poly':1,'Si': 1} # n.b. fluid dielectric constant will be calc'd by conv_to_ratchet_ND_pars

surf={'pK': 7.5,'pH': 6,'alpha': 6e23*1e-3, 'Gamma': 8e-18,'C_stern': 2.9}


#-- object numbering
SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':1,'Si':2}}
const={'e': 1.602e-19,'e0': 8.85e-12,'kB': 1.38e-23,\
           'Na': 6.022e23}
Tkelv = Tdeg+273


#--- calc ND surf pars
surf_ND = dict()
for fld in ['pK','pH','alpha']:
    surf_ND[fld] = surf[fld]
surf_ND['Gamma']= surf_ND['Gamma']*geom_ND['l0_m']**2/const['e']


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
    
    
def calc_phi_d_Beh(sigma_ND):
    ''' 
    calc the relationship between the Stern potential and the surface charge dictated
    by the reaction coefficients for the dissociation of the Silanol groups
    '''
    
    return 1
    
    


for p in range(len(psi_mVr)):
    for fld in ['Si','poly']:
        psi_mV[fld] = psi_mVr[p]

    geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)
    Fr_int[p],Fc_int,u_,mesh = ph.solve_ratchet_problem(geom_ND,params_ND,\
        SDs,fName_tail,doPlots=False,cylSim=True,N_adapt = 0,adapt_TOL=adapt_TOL)
        
    #--- extract data
    phi_r_ND[p] = u_(0,geom_ND['L'])
    sigma_r_ND[p] = params_ND['sigma_Si']


#--- plot
f0,ax = plt.subplots(nrows=1,ncols=2)

psi_NDi = np.linspace(0,np.max(psi_mVr)*1e-3*const['e']/(const['kB']*Tkelv))
sigma_G = lambda psi_ND: np.sqrt(8*params_ND['n0']*params_ND['er_fluid'])*np.sinh(0.5*psi_ND) # the Grahame equation


ax[0].plot(phi_r_ND,sigma_r_ND,'ok')
ax[0].plot(psi_NDi,sigma_G(psi_NDi),'-k')
ax[0].set_xlabel('$\phi$')
ax[0].set_ylabel('$\sigma$')

fh.place_fig(f0)









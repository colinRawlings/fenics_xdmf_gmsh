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

use_min = False

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
surf=dict()
surf['Silanol']={'pK': 7.5, 'Gamma': 9e18,'C_stern': 2.5} # Silica
surf['Carboxyl']={'pK': 4.9, 'Gamma': 0.25e18,'C_stern': 100} # Carboxyl


SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':1,'Si':2}} # domain numbering

geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)

#--- calc ND surf pars
surf_ND=dict()
for mat in ['Silanol','Carboxyl']:
    surf_ND[mat] = ph.calc_surf_ND_pars(surf[mat],geom_ND,params_ND)    


#--- required variables for calcs
Fr_int = np.zeros(len(psi_mVr))
sigma_r_ND = np.zeros(len(psi_mVr))
phi_r_ND = np.zeros(len(psi_mVr))


#---- fig2 data
pH_beh = np.arange(4,11,1)
scaleFac = -50e-3/55.5
sigma_beh={'Silanol': np.asarray([0,2,5,11,24,47,float('nan')])*scaleFac,\
            'Carboxyl': np.asarray([2,5.5,12,22,36,43,44.5])*scaleFac}


##============== define equations
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
    
    
def calc_err(sigma_ND,pH,surf):
    
    err = ph.calc_phi_d_Grahame(sigma_ND,params_ND)-ph.calc_phi_d_SiOH_Beh(sigma_ND,pH,surf)
    
    return err

#--- calc fig 2
pH_r = np.linspace(3,11,50)

sigma_ND_r = dict()
N_it = 0 # count total evaluations required

for mat in ['Silanol','Carboxyl']:
    sigma_ND_r[mat] = np.zeros(np.shape(pH_r))    
    s0 = -1
        
    for p in range(len(pH_r)):
        
        if use_min:
            RES = sp.optimize.minimize(lambda s: calc_err(s,pH_r[p],surf_ND[mat])**2,s0,\
                    bounds=[-surf_ND[mat]['Gamma'],0])
            N_it += RES['nfev']
        
            sigma_ND_r[mat][p] = float('nan')    
            if RES['success']:
                sigma_ND_r[mat][p] = RES['x'][0]
                s0 = RES['x'][0]
            else:
                print('calc for pH: %.1f failed'%(pH_r[p]))        

        else: # using fsolve
            RES = sp.optimize.fsolve(lambda s: calc_err(s,pH_r[p],surf_ND[mat]),s0,\
                    full_output = True,factor=1e-2)     
                    
            sigma_ND_r[mat][p] = RES[0][0]
            s0 = RES[0][0]
            
            N_it += RES[1]['nfev']
           

print('Finished %d calcs in %d iterations'%(2*len(pH_r),N_it))


##--- display
f0,ax = plt.subplots(nrows=1,ncols=2)
fh.make_cfigw()

q=0
for mat in ['Silanol','Carboxyl']:
    col = fh.chot(q,2)     
    
    ax[0].plot(pH_r,sigma_ND_r[mat],'-',label='fenics: '+mat,color=col)
    
    ax[1].plot(pH_r,sigma_ND_r[mat]*const['e']/geom_ND['l0_m']**2,'-',\
                label='fenics: '+mat,color=col)
    ax[1].plot(pH_beh,sigma_beh[mat],'o',color=col)
    q+=1


ax[0].set_xlabel('pH')
ax[0].set_ylabel('$\sigma l_0^2/e$ ')
ax[0].legend()

ax[1].set_xlabel('pH')
ax[1].set_ylabel('$\sigma$ (C/m$^2$)')
ax[1].legend()
ax[1].set_ylim((-55e-3,0))
ax[1].set_title('$n_0$: %.1e M ($\lambda_D$: %.1f nm)'%(ph.lD_nm_to_n0_M(lD_nm,Tdeg),lD_nm))


fh.place_fig(f0)
    








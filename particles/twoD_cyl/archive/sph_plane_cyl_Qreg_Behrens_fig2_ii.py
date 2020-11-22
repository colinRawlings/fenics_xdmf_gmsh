# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 18:16:20 2015

@author: cra

calculates F for a range of particle positions in the gap


"""

execfile('/home/cra/dolfin/useful_f/init_cra.py')

Tdeg = 21
lD_nm = 9.6 # Debye length


#--- surface models
surf=dict()
surf['Silanol']={'pK': 7.5, 'Gamma': 9e18,'C_stern': 2.5} # Silica
surf['Carboxyl']={'pK': 4.9, 'Gamma': 0.25e18,'C_stern': 100} # Carboxyl


SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':1,'Si':2}} # domain numbering


#--- calc ND surf pars
surf_ND=dict()
for mat in ['Silanol','Carboxyl']:
    surf_ND[mat] = ph.calc_surf_ND_pars(surf[mat],Tdeg)    


#---- fig2 data
pH_beh = np.arange(4,11,1)
scaleFac = -50e-3/55.5
sigma_beh={'Silanol': np.asarray([0,2,5,11,24,47,float('nan')])*scaleFac,\
            'Carboxyl': np.asarray([2,5.5,12,22,36,43,44.5])*scaleFac}
const = ph.load_const()
l0_m = ph.calc_l0_m(Tdeg)
    
    
def calc_err(sigma_ND,pH,surf):
    
    err = ph.calc_phi_d_Grahame(sigma_ND,lD_nm,Tdeg,ph.calc_er_fluid(Tdeg))-ph.calc_phi_d_SiOH_Beh(sigma_ND,pH,surf)
    
    return err

#--- calc fig 2
pH_r = np.linspace(3,11,50)

sigma_ND_r = dict()
N_it = 0 # count total evaluations required

for mat in ['Silanol','Carboxyl']:
    sigma_ND_r[mat] = np.zeros(np.shape(pH_r))    
    s0 = -1
        
    for p in range(len(pH_r)):
        
        x,infodict,ier,mesg = sp.optimize.fsolve(lambda s: calc_err(s,pH_r[p],surf_ND[mat]),s0,\
                full_output = True,factor=1e-2)     
                
        #--- unpack
        sigma_ND_r[mat][p] = x[0]
        s0 = x[0]
        
        N_it += infodict['nfev']
           

print('Finished %d calcs in %d iterations'%(2*len(pH_r),N_it))


##--- display
f0,ax = plt.subplots(nrows=1,ncols=2)
fh.make_cfigw()

q=0
for mat in ['Silanol','Carboxyl']:
    col = fh.chot(q,2)     
    
    ax[0].plot(pH_r,sigma_ND_r[mat],'-',label='fenics: '+mat,color=col)
    
    ax[1].plot(pH_r,sigma_ND_r[mat]*const['e']/l0_m**2,'-',\
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
    








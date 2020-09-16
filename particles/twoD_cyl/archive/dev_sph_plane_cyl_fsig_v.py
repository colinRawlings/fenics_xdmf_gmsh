# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 18:16:20 2015

@author: cra

calculates F for a range of particle positions in the gap


"""

execfile('/home/cra/dolfin/useful_f/init_cra.py')

fName_tail = r'/home/cra/dolfin/particles/twoD_cyl/mesh/cyl_sph_plane_tail_02.geo'
fName_log = r'/home/cra/dolfin/particles/twoD_cyl/logFile.txt'


Rs_mesh = 0.3 # cell size as a function of 1/k0 at the charge sheets


lD_nm = 9.6 # Debye length
Tdeg = 21
w = 200 # field size
a_small = 1

psi_mV={'Si':67,'sph':1e-6,'poly':1e-6}
geom_nm={'w':w,'L':160,'h':50,'a':a_small,'sb':60,'st':30,'dsr':Rs_mesh*lD_nm,'dss': a_small} # cylindrical symmetry
er_dry = {'sph': 1, 'poly':1,'Si': 1} # n.b. fluid dielectric constant will be calc'd by conv_to_ratchet_ND_pars
surf={'pK': 7.5, 'Gamma': 9e18,'C_stern': 2.5}


#--- calc hr
h_min = lD_nm # don't let the sphere actually touch the bounding wall
        
geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)
surf_ND = ph.calc_surf_ND_pars(surf,geom_ND,params_ND)
const = ph.load_const()

#-- object numbering
SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':1,'Si':2}}

#---- fig2 data
pH_beh = np.arange(4,9,1)
scaleFac = -50e-3/55.5
sigma_beh={'Silanol': np.asarray([0,2,5,11,24,47,float('nan')])*scaleFac,\
            'Carboxyl': np.asarray([2,5.5,12,22,36,43,44.5])*scaleFac}

Nx = 1*geom_nm['L']/lD_nm

#--- construct spline rep
x_ND = np.linspace(0,geom_ND['w'],Nx)
xc = np.linspace(0,geom_ND['w'],Nx*10)
sigma_ND = 10*np.sin(x_ND*10)

class pyExpression(Expression):
    def __init__(self,sci):
        self.sci = sci
    def eval(self,value,x):
        value[0]=self.sci(x[0])

def calc_phi_d_PB(x_ND,sigma_ND):
    ''' 
    calc the relationship between the Stern potential and the surface charge dictated
    by the PB relationship
    '''
    
    #--- fit the spline
    sci = fh.spline_con_gradient(x_ND,sigma_ND,[],[],k=3)
    params_ND['sigma_Si'] = avg(pyExpression(sci))
    
    #--- get phi_d
    Fr_int,Fc_int,u_,mesh = ph.solve_ratchet_problem(geom_ND,params_ND,SDs,fName_tail,\
        doPlots=False,cylSim=True)
    
    phi_d = np.zeros(len(x_ND))
    for p in range(len(x_ND)):
        phi_d[p] = u_(x_ND[p],geom_ND['L'])
        
    return phi_d
    
def calc_error(x_ND,sigma_ND,pH,surf_ND):
    err = calc_phi_d_PB(x_ND,sigma_ND) \
            - ph.calc_phi_d_SiOH_Beh(sigma_ND,pH,surf_ND)
    print('!!!!!!!!!!!RMS err (%d/%d): %.1e'%(p,len(pH_beh),np.dot(err,err)))
    return err
    
#--- setup: solve
s0 = -1*np.ones(len(x_ND))
f0,axUP = plt.subplots(nrows=1,ncols=1)
fh.place_fig(f0)



##--- benchmark vs Grahame
#s0 = np.linspace(1000,7000,len(pH_beh))
#s0c = np.linspace(1000,7000,len(pH_beh)*10)
#
sigma_NDr = np.zeros((len(pH_beh),len(x_ND)))
#phi_NDr = sigma_NDr.copy()
#
#for p in range(len(s0)):
#    s0_x = np.ones(np.shape(x_ND))*s0[p]
#    phi_NDr[p,:] = calc_phi_d_PB(x_ND,s0_x)
#    
#
#axUP.plot(s0,phi_NDr[:,0],'ok',label='fenics')
#axUP.plot(s0c,ph.calc_phi_d_Grahame(s0c,params_ND),'-r',label='Grahame')
#axUP.set_xlabel('$\overline\sigma$')
#axUP.set_ylabel('$\overline\phi_d$')
#axUP.legend()


for p in range(len(pH_beh)):
    RES = sp.optimize.fsolve(lambda s: calc_error(x_ND,s,pH_beh[p],surf_ND),s0,\
        full_output = True,factor=1,xtol=1e-2)     
    sigma_NDr[p,:] = RES[0]
    sci = fh.spline_con_gradient(x_ND,sigma_NDr[p,:],[],[],k=3)    
    #        s0 = sigma_NDr[p,:]    
    
    axUP.plot(x_ND,sigma_NDr[p,:],'+',color=fh.chot(p,len(pH_beh)))
    axUP.plot(xc,sci(xc),'-',color=fh.chot(p,len(pH_beh)),label='fenics, pH: %.1f'%(pH_beh[p]))
    
    axUP.plot(x_ND,sigma_beh['Silanol'][p]*np.ones(len(x_ND))/(const['e']/geom_ND['l0_m']**2)\
        ,'o',color=fh.chot(p,len(pH_beh)))    

axUP.legend()
axUP.set_xlabel('$x/l_0$')
axUP.set_ylabel('$\overline \sigma_ND$')    
    
    
plt.draw()






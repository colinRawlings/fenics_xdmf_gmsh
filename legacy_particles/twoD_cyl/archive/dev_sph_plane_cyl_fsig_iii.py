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
Nx = 10

lD_nm = 30 # Debye length
Tdeg = 21
w = 200 # field size
a_small = 1

psi_mV={'Si':67,'sph':1e-6,'poly':1e-6}
geom_nm={'w':w,'L':160,'h':50,'a':a_small,'sb':60,'st':30,'dsr':Rs_mesh*lD_nm,'dss': a_small} # cylindrical symmetry
er_dry = {'sph': 1, 'poly':1,'Si': 1} # n.b. fluid dielectric constant will be calc'd by conv_to_ratchet_ND_pars

#--- calc hr
h_min = lD_nm # don't let the sphere actually touch the bounding wall
        
geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)

#-- object numbering
SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':1,'Si':2}}


#--- construct spline rep
x_ND = np.linspace(0,geom_ND['w'],Nx)
xc = np.linspace(0,geom_ND['w'],Nx*10)
sigma_ND = 10*np.sin(x_ND*10)

class pyExpression(Expression):
    def __init__(self,sci):
        self.sci = sci
    def eval(self,value,x):
        value[0]=self.sci(x[0])

def calc_phi_d_PB(x_ND,sigma_ND,pyE):
    ''' 
    calc the relationship between the Stern potential and the surface charge dictated
    by the PB relationship
    '''
    
    #--- fit the spline
    sci = fh.spline_con_gradient(x_ND,sigma_ND,[0,geom_ND['w']],[0,0],k=3)
    params_ND['sigma_Si'] = avg(pyExpression(sci))
    
    #--- get phi_d
    Fr_int,Fc_int,u_,mesh = ph.solve_ratchet_problem(geom_ND,params_ND,SDs,fName_tail,\
        doPlots=True,cylSim=True)
    
    phi_d = np.zeros(len(x_ND))
    for p in range(len(x_ND)):
        phi_d[p] = u_(x_ND[p],geom_ND['L'])
        
    return phi_d

pyE = pyExpression(None)


phi_d = calc_phi_d_PB(x_ND,sigma_ND,pyE)



f0,ax = plt.subplots(nrows=1,ncols=2)
ax[0].plot(x_ND,sigma_ND,'ok',label='sampled sin')
ax[0].legend()

ax[1].plot(x_ND,phi_d,'ok')


fh.place_fig(f0)



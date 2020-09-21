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
Nx = 5

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
sc = 10*np.sin(x_ND*30)

sci = fh.spline_con_gradient(x_ND,sc,[0,geom_ND['w']],[0,0],k=len(x_ND)-1)

f0,ax = plt.subplots(nrows=1,ncols=1)
ax.plot(x_ND,sc,'ok',label='sampled sin')
ax.plot(xc,sci(xc),'-k',label='spline')
ax.legend()

#--- for now simply collect sci from global
class pyExpression(Expression):
    def __init__(self,sci):
        self.sci = sci
    def eval(self,value,x):
        value[0]=self.sci(x[0])

fh.place_fig(f0)

pyE = pyExpression(sci)

params_ND['sigma_Si'] = avg(pyExpression(sci))

Fr_int,Fc_int,u_,mesh = ph.solve_ratchet_problem(geom_ND,params_ND,SDs,fName_tail,\
    doPlots=True,cylSim=True)


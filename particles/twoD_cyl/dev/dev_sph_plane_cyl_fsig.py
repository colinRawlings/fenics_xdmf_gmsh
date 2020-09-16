# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 11:15:06 2016

@author: cra

pass interp function for use specifying a BC


"""

execfile('/home/cra/dolfin/useful_f/init_cra.py')



fName_tail = r'/home/cra/dolfin/particles/twoD_cyl/mesh/cyl_sph_plane_tail_01.geo'
fName_log = r'/home/cra/dolfin/particles/twoD_cyl/logFile.txt'


Rs_mesh = 0.3 # cell size as a function of 1/k0 at the charge sheets

Nh = 10 # number of different particle positions to try
lD_nm = 30 # Debye length
Tdeg = 21
w = 200 # field size
a_small = 1

psi_mV={'Si':67,'sph':1e-3,'poly':67}
geom_nm={'w':w,'L':160,'h':50,'a':a_small,'sb':60,'st':60,'dsr':Rs_mesh*lD_nm,'dss':a_small}
er_dry = {'sph': 1, 'poly':1,'Si': 1} # n.b. fluid dielectric constant will be calc'd by conv_to_ratchet_ND_pars
  
       
geom_ND,params_ND=ph.conv_to_ratchet_ND_pars(geom_nm,lD_nm,er_dry,Tdeg,psi_mV)



#-- object numbering
SDs={'vol':{'sph':3,'poly':2,'Si':1,'fluid':0},\
     'surf':{'edge':4,'sph':3,'poly':2,'Si':1}}

#-- construct basis
mesh,sd_MF,bnd_MF=fh.params_to_mesh(geom_ND,fName_tail,'out.geo')
#fh.paraview_mesh(mesh,sd_MF=sd_MF)

V = FunctionSpace(mesh,'CG',2)
u_ = Function(V)
u = TrialFunction(V)
v = TestFunction(V)

#-- define equation
ds = Measure('dS')[bnd_MF]
dx = Measure('dx')[sd_MF]
dx_dry = dx(SDs['vol']['sph'])+dx(SDs['vol']['Si'])+dx(SDs['vol']['poly'])
dx_wet = dx(SDs['vol']['fluid'])  

r,z = SpatialCoordinate(mesh) 

sig_pl1 = Constant(0)

wavelength = 20

class pyExpression(Expression):
    def eval(self,value,x):
        value[0]=np.sin(x[0]*wavelength)


wavelength = 10


#sig_pl2 = Expression('sin(10*x[0])')
#sig_pl2 = Constant(params_ND['sigma_poly'])
sig_pl2 = pyExpression()


F = Constant(params_ND['sigma_sph'])*v('+')*ds(SDs['surf']['sph'])\
    + sig_pl1*v('+')*ds(SDs['surf']['Si'])\
    + avg(sig_pl2)*v('+')*ds(SDs['surf']['poly'])\
    + v/(r+1e-6)*Dx(u,0)*Constant(params_ND['er_fluid'])*dx_wet\
    + v/(r+1e-6)*Dx(u,0)*Constant(params_ND['er_dry']['Si'])*dx_dry\
    - Constant(2*params_ND['n0'])*sinh(u)*v*dx_wet\
    - dot(Constant(params_ND['er_fluid'])*grad(u),grad(v))*dx_wet\
    - dot(Constant(params_ND['er_dry']['Si'])*grad(u),grad(v))*dx_dry
    
F = action(F,u_)
J = derivative(F,u_,u)

#- n.b. with |J| added 
u_el = Constant(params_ND['er_fluid']*0.5)*dot(grad(u_),grad(u_))*2*pi*r*dx_wet+\
       Constant(params_ND['er_dry']['Si']*0.5)*dot(grad(u_),grad(u_))*2*pi*r*dx_dry
Tds = 2*Constant(params_ND['n0'])*(cosh(u_)-1.0-u_*sinh(u_))*2*pi*r*dx_wet

f = u_el - Tds

problem = NonlinearVariationalProblem(F,u_,[],J)

solver = NonlinearVariationalSolver(problem)
solver.parameters['newton_solver']['linear_solver'] ='gmres'
solver.parameters['newton_solver']['error_on_nonconvergence'] = False
solver.solve()

F_int=assemble(f) 

#fh.paraview_fld(u_) 


#--- make plots
Np = 100
f0,axs = plt.subplots(nrows=1,ncols=2)
fh.make_cfigw()
fh.place_fig(f0)

fh.paraview_fld(u_)

zXX = np.linspace(-geom_ND['sb'],geom_ND['L']+geom_ND['st'],Np)
rXX = np.linspace(0,geom_ND['w'],Np)

uXXz = np.zeros(Np)
uXXr = np.zeros(Np)

for q in range(Np):
    uXXz[q] = u_(0,zXX[q])
    uXXr[q] = u_(rXX[q],geom_ND['L']-geom_ND['h']-geom_ND['a'])

#-- plots
axs[0].plot(zXX*geom_ND['l0_m']*1e9,uXXz,'-k')
axs[0].set_xlabel('$z$ (nm)')
axs[0].set_ylabel('$\phi$')

axs[1].plot(rXX*geom_ND['l0_m']*1e9,uXXr,'-k')
axs[1].set_xlabel('$r$ (nm)')
axs[1].set_ylabel('$\phi$')

plt.draw()


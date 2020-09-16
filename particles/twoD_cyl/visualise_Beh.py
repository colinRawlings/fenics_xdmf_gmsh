# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 14:58:13 2016

@author: cra

visualise generalised charge-potential relation when the capacitance
of the Stern layer can be neglected.

"""

from __future__ import print_function
exec(compile(open('/home/cra/dolfin/useful_f/init_cra.py').read(), '/home/cra/dolfin/useful_f/init_cra.py', 'exec'))

Tkelv = 293

x = np.linspace(-10,10,1e2)

f = lambda x: 1./(1+np.exp(-x))

pH = 6

f0,axs=plt.subplots(nrows=1,ncols=1)
fh.place_fig(f0)

axs.plot(x,f(x),'-k')
axs.set_xlabel('$\phi e / k_B T +(\mathrm{pK}-\mathrm{pH})\log 10$')
axs.set_ylabel('$-\sigma / e \Gamma $')


#--- Stern term
surf['Silanol']={'pK': 7.5,'alpha': 1, 'Gamma': 9e18,'C_stern': 2.5} # Silica
surf['Carboxyl']={'pK': 4.9,'alpha': 1, 'Gamma': 0.25e18,'C_stern': 100} # Carboxyl
const={'e': 1.602e-19,'e0': 8.85e-12,'kB': 1.38e-23,\
           'Na': 6.022e23}

phi_act=dict()
phi_Stern = dict()

for mat in ['Silanol','Carboxyl']:
    phi_act[mat] = (surf[mat]['pK']-pH)*np.log(10)
    phi_Stern[mat] = const['e']*surf[mat]['Gamma']/surf[mat]['C_stern']*const['e']/(Tkelv*const['kB'])
    
    print(('%s: pH: %.1f, phi_act: %.3f, phi_Stern: %.3f'%(mat,pH,phi_act[mat],phi_Stern[mat])))


#-- log term
r = np.linspace(1e-6,1,1e3)
phi_log = lambda r: -np.log(1./r-1)

f0,axs=plt.subplots(nrows=1,ncols=1)
fh.place_fig(f0)


axs.plot(r,phi_log(r),'-k',label='$\ln$')
axs.plot(r,-r*phi_Stern['Silanol'],'-r',label='Stern: Silanol')
axs.plot(r,-r*phi_Stern['Carboxyl'],'-b',label='Stern: Carboxyl')
axs.set_xlabel('$-\sigma/ e \Gamma$')
axs.set_ylabel('$e\phi/k_B T$')
axs.legend()


# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 11:19:15 2016

@author: cra

check Debye to n0 relationship

"""

execfile('/home/cra/dolfin/useful_f/init_cra.py')


lD_nm_r = np.asarray([9.6,30,275])
n0_M_r = np.asarray([1e-3,1e-4,1.2e-6])

Tdeg = 20



           
           
#-- forward           
for lD in lD_nm_r:
    n0_M = ph.lD_nm_to_n0_M(lD,Tdeg)
    print('lD: %.1f nm, n0: %.1e M'%(lD,n0_M))
    
    
#-- backward           
for n0 in n0_M_r:
    lD_nm = ph.n0_M_to_lD_nm(n0,Tdeg)
    print('lD: %.1f nm, n0: %.1e M'%(lD_nm,n0))
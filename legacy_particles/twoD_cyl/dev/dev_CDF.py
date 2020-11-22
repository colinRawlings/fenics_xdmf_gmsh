# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 11:34:00 2016

@author: cra
"""

#--- import explicitly to help rope
import numpy as np
import fenics_helpers as fh
import particle_helpers as ph
import scipy.interpolate as si
import scipy.optimize as so

execfile('/home/cra/dolfin/useful_f/init_cra.py')


#--- develop the CDF analysis of particle position
q = len(hr_nm)-1

x = hr_nm[q]+geom_nm['a']
F = Fr_int_gap[q,:]

fit = ph.polyfit_minima(hr_nm[q]+geom_nm['a'],Fr_int_gap[q,:],kT_lim_init,kT_roots)


#--- 
#f,a = fh.cfigw()    
#a[0].plot(x,F,'ok')
#a[0].plot(fit['xc'],fit['yc'])


#def Fc(x_eval,x,F,fit):
x_eval = np.linspace(np.min(x),np.max(x),1e3)

Fout = np.zeros(np.shape(x_eval))    

Ineg = x_eval<np.min(fit['xc'])
Imin = (x_eval>=np.min(fit['xc'])) & (x_eval<np.max(fit['xc']))
Ipos = x_eval>=np.max(fit['xc'])

Fout[Imin] = fit['fP'](x_eval[Imin])
Fout[~Imin] = si.interp1d(x,F,kind='cubic')(x_eval[~Imin])





#--- now conv to int    
thresh = 0.05

P_unnorm = np.exp(-Fout+np.min(Fout))
P_norm = P_unnorm/np.sum(P_unnorm)
CDF = np.cumsum(P_norm)
CDFf = si.interp1d(x_eval,CDF,kind='cubic')

#--- find threshold values
x_neg = so.fsolve(lambda x: CDFf(x)-thresh,fit['x_min'])
x_pos = so.fsolve(lambda x: CDFf(x)-(1-thresh),fit['x_min'])


f,a = fh.cfigw(ncols=3)    
a[0].plot(x,F,'ok')
a[0].plot(fit['xc'],fit['yc'])
a[0].plot(x_eval,Fout,'-r')

a[1].plot(x_eval,P_norm,'-r')
a[2].plot(x_eval,CDFf(x_eval),'-r')
a[2].plot([x_neg,x_neg],[0,1],'--b')
a[2].plot([x_pos,x_pos],[0,1],'--b')

#--- build results
RES={}
RES['x_lo'] = x_neg
RES['x_hi'] = x_pos
RES['x_min'] = fit['x_min']
RES['xc'] = x_eval
RES['Fc'] = Fout
RES['P'] = P_norm
RES['CDF'] = CDF




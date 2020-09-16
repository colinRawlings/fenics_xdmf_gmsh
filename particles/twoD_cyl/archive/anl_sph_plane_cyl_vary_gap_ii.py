# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 18:10:44 2016

@author: cra
"""


from scipy.optimize import curve_fit
import fenics_helpers as fh
import particle_helpers as ph
import pickle
import matplotlib.pyplot as plt
import numpy as np
import fenics_helpers as fh
import matplotlib

reload(fh)
reload(ph)

#---- unpack data
fName_save = r'/home/cra/dolfin/particles/twoD_cyl/vary_L_h.pickle'
flds = ['geom_ND','params_ND','geom_nm','Lr_nm','hr_nm','Fr_int_gap',\
            'Fr_int_sph_pl','Fr_sph','hr_sph_pl_nm']
Rl = pickle.load(open(fName_save,'rb'))

for key,val in Rl.iteritems():
    exec(key+r'=val' )
lD_nm = 30
#----- 

plt.close('all')
f0,axs = plt.subplots(nrows=1,ncols=2)
fh.make_cfigw()


##--- analyse exponential fit to data
#--- fit to all data points
def f_exp(x,A0,l,c):
    return A0*np.exp(-x/l)+c

c_init = np.min(Fr_int_sph_pl)
l_init = lD_nm
A0_init = np.max(Fr_int_sph_pl)-np.min(Fr_int_sph_pl)

hi = np.linspace(np.min(hr_sph_pl_nm),np.max(hr_sph_pl_nm))

popt,pcov = curve_fit(f_exp,hr_sph_pl_nm,Fr_int_sph_pl,p0=[A0_init,l_init,c_init])

axs[0].plot(hr_sph_pl_nm/lD_nm,Fr_int_sph_pl,'o',label='sim')
axs[0].plot(hi/lD_nm,f_exp(hi,popt[0],popt[1],popt[2]),'-',label='exp fit')

#--- predict F energy based on single wall interactions
def h2(h1,geom_nm):
    return geom_nm['L']-2*geom_nm['a']-h1

#Fr_sph_pl = sci.interp1d(hr_sph_pl_nm,Fr_int_sph_pl,kind='quadratic',bounds_error=False)
def Fr_sph_pl(h):
    return f_exp(h,popt[0],popt[1],popt[2])
    
hi = np.linspace(np.min(hr_sph_pl_nm),np.max(hr_sph_pl_nm))
#Fr0 = 2*np.min(Fr_int_sph_pl)-Fr_sp

def Fr_pl_sph_pl(Fr_sph_pl,h,geom_nm):
    return Fr_sph_pl(h)+Fr_sph_pl(h2(h,geom_nm))-Fr_sph

def xc(h,geom_nm):
    return h-geom_nm['L']*0.5+geom_nm['a']

def hrc(h,geom_nm):
    x = xc(h,geom_nm)
    return np.linspace(-np.max(np.abs(x)),np.max(np.abs(x)))+geom_nm['L']*0.5-geom_nm['a']

def xcc(h,geom_nm):
    return xc(hrc(h,geom_nm),geom_nm)

cm = matplotlib.cm.hot

for q in range(len(Lr_nm)):

    geom_nm['L']=Lr_nm[q]
    col = cm(1.*q/(1.*len(Lr_nm)))    
    
    #--- constant offset (all curves)
    axs[1].plot(xc(hr_nm[q],geom_nm),Fr_int_gap[q,:],'o',color=col,\
            label='L-2a: %.1f$\lambda_D$'%((Lr_nm[q]-2*geom_nm['a'])/lD_nm))
    axs[1].plot(xcc(hr_nm[q],geom_nm),Fr_pl_sph_pl(Fr_sph_pl,hrc(hr_nm[q],geom_nm),geom_nm),'-',\
                color=col)
                #    axs[1].plot(np.ones(2)*geom_nm['L']*0.5-geom_nm['a'],[0,np.max(Fr_int_gap[q,:])],':',color=col)
                #    axs[1].plot(-np.ones(2)*geom_nm['L']*0.5+geom_nm['a'],[0,np.max(Fr_int_gap[q,:])],':',color=col)
    
#    #--- constant offset (all curves)
#    Ffit = Fr_pl_sph_pl(Fr_sph_pl,hi,geom_nm)
#    
#    axs[2].plot(hr_nm[q]-Lr_nm[q]*0.5+geom_nm['a'],Fr_int_gap[q,:]-np.min(Fr_int_gap[q,:]),'o',\
#        color=cm(1.*q/(1.*len(Lr_nm))),label='L: %d nm'%(Lr_nm[q]))
#    axs[2].plot(hi-Lr_nm[q]*0.5+geom_nm['a'],Ffit-np.min(Ffit[~np.isnan(Ffit)]),'-',\
#                color=col)
#                #    axs[2].plot(np.ones(2)*geom_nm['L']*0.5-geom_nm['a'],[0,np.max(Fr_int_gap[q,:]-Fr0)],':',color=col)
#                #    axs[2].plot(-np.ones(2)*geom_nm['L']*0.5+geom_nm['a'],[0,np.max(Fr_int_gap[q,:]-Fr0)],':',color=col)
    
t_str = r'$\bar\phi_{sph}$: %.1f, $\bar\phi_{pl}$: %.1f'%(params_ND['psi_sph'],params_ND['psi_Si'])

for ax in axs:
    ax.set_xlabel('$x / \lambda_D$')
    ax.set_ylabel('$F/(k_B T)$')
    ax.set_title(t_str)
    ax.legend()


fh.place_fig(f0)

plt.draw()
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 18:10:44 2016

@author: cra
"""


from scipy.optimize import curve_fit
import fenics_helpers as fh

reload(fh)

plt.close('all')
f0,axs = plt.subplots(nrows=1,ncols=3)
fh.make_cfigw()


##--- analyse exponential fit to data
#--- fit to all data points
def f_exp(x,A0,l,c):
    return A0*np.exp(-x/l)+c

c_init = np.min(Fr_int_sph_pl)
l_init = lD_nm
A0_init = np.max(Fr_int_sph_pl)-np.min(Fr_int_sph_pl)

popt,pcov = curve_fit(f_exp,hr_sph_pl_nm,Fr_int_sph_pl,p0=[A0_init,l_init,c_init])

axs[0].plot(hr_sph_pl_nm,Fr_int_sph_pl,'o',label='sim')
axs[0].plot(hi,f_exp(hi,popt[0],popt[1],popt[2]),'-',label='exp fit')


for q in range(len(Lr_nm)):

    geom_nm['L']=Lr_nm[q]
    col = cm(1.*q/(1.*len(Lr_nm)))    
    
    
    #--- constant offset (all curves)
    axs[1].plot(hr_nm[q]-Lr_nm[q]*0.5+geom_nm['a'],Fr_int_gap[q,:]-Fr0,'o',\
        color=cm(1.*q/(1.*len(Lr_nm))),label='L: %d nm'%(Lr_nm[q]))
    axs[1].plot(hi-Lr_nm[q]*0.5+geom_nm['a'],Fr_pl_sph_pl(Fr_sph_pl,hi,geom_nm)-Fr0-Fr_sph,'-',\
                color=cm(1.*q/(1.*len(Lr_nm))))
    axs[1].plot(np.ones(2)*geom_nm['L']*0.5-geom_nm['a'],[0,np.max(Fr_int_gap[q,:]-Fr0)],':',color=col)
    axs[1].plot(-np.ones(2)*geom_nm['L']*0.5+geom_nm['a'],[0,np.max(Fr_int_gap[q,:]-Fr0)],':',color=col)
    
    #--- constant offset (all curves)
    Ffit = Fr_pl_sph_pl(Fr_sph_pl,hi,geom_nm)
    
    axs[2].plot(hr_nm[q]-Lr_nm[q]*0.5+geom_nm['a'],Fr_int_gap[q,:]-np.min(Fr_int_gap[q,:]),'o',\
        color=cm(1.*q/(1.*len(Lr_nm))),label='L: %d nm'%(Lr_nm[q]))
    axs[2].plot(hi-Lr_nm[q]*0.5+geom_nm['a'],Ffit-np.min(Ffit[~np.isnan(Ffit)]),'-',\
                color=col)
                #    axs[2].plot(np.ones(2)*geom_nm['L']*0.5-geom_nm['a'],[0,np.max(Fr_int_gap[q,:]-Fr0)],':',color=col)
                #    axs[2].plot(-np.ones(2)*geom_nm['L']*0.5+geom_nm['a'],[0,np.max(Fr_int_gap[q,:]-Fr0)],':',color=col)
    

for ax in axs:
    ax.set_xlabel('$x$ (nm)')
    ax.set_ylabel('$F/(k_B T)$')
    fh.outer_legend(ax)




fh.place_fig(f0)

plt.draw()
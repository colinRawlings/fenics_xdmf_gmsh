# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 18:10:44 2016

@author: cra
"""
from __future__ import absolute_import
import six
from six.moves import range
exec(compile(open('/home/cra/dolfin/useful_f/init_cra.py').read(), '/home/cra/dolfin/useful_f/init_cra.py', 'exec'))

from scipy.optimize import curve_fit
import fenics_helpers as fh
import particle_helpers as ph
import pickle
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import fenics_helpers as fh

reload(fh)
reload(ph)

figFolder = '/home/cra/Documents/pres/sph_cyl/figs/'

#---- unpack data
#fName_save = r'/home/cra/dolfin/particles/twoD_cyl/vary_L_h_many_pts_symQpl_er1.pickle' # <-- init (Nov rework)
fName_save = r'/home/cra/dolfin/particles/twoD_cyl/20161111-162938vary_L_h.pickle'

flds = ['geom_ND','params_ND','geom_nm','Lr_nm','hr_nm','Fr_int_gap',\
            'Fr_int_sph_pl','Fr_sph','hr_sph_pl_nm','Fr_int_pl_pl']
Rl = pickle.load(open(fName_save,'rb'))

for key,val in six.iteritems(Rl):
    exec(key+r'=val' )
lD_nm = geom_ND['lD_nm']
#----- 

plt.close('all')
f0,axs = plt.subplots(nrows=2,ncols=2)
axs=np.reshape(axs,np.prod(np.shape(axs))) 
fh.make_cfigb()
t_str = r'$a/\lambda_D$: %.1f, $\bar\phi_{sph}$: %.1f, $\bar\phi_{pl}$: %.1f'\
        %(geom_nm['a']/lD_nm,params_ND['psi_sph'],params_ND['psi_Si'])
cm = matplotlib.cm.hot

##--- anl pl-pl data (constant charge)
Li_nm = np.linspace(np.min(Lr_nm),np.max(Lr_nm))


def Fr_VO_pl_pl(L_nm):
    n0_m = params_ND['n0']/geom_ND['l0_m']**3 # scaled electrolyte concentration
    A_m = np.pi*(geom_ND['w']*geom_ND['l0_m'])**2
    return A_m*4*params_ND['psi_Si']**2*lD_nm*1e-9*n0_m*np.ones(np.shape(L_nm))*np.exp(-L_nm/lD_nm)/(1.+np.exp(-L_nm/lD_nm))

f,axs0 = fh.pfig()

axs0.loglog(Lr_nm/lD_nm,Fr_int_pl_pl-np.min(Fr_int_pl_pl),'ok',label='sim')
axs0.loglog(Li_nm/lD_nm,Fr_VO_pl_pl(Li_nm),'-r',label='Verwey')
axs0.legend(loc='best',frameon=False)
axs0.set_title(r'plane-plane (in sim. area), $e \phi/kT$: %.1f'%(params_ND['psi_Si']))
fh.less_yticks(axs0)
#fh.less_yticks(axs0)
fh.cxlbl(axs0,'L/\lambda_D','(F_{pl-pl}(L) - F_{pl-pl}(\infty)) / (k_B T)')
f.savefig(figFolder+'bench_pl_pl.pdf')



##--- analyse exponential fit to data
#--- fit to all data points
def f_exp(x,A0,l,c):
    return A0*np.exp(-x/l)+c

c_init = np.min(Fr_int_sph_pl)
l_init = lD_nm
A0_init = np.max(Fr_int_sph_pl)-np.min(Fr_int_sph_pl)

hi = np.linspace(np.min(hr_sph_pl_nm),np.max(hr_sph_pl_nm))

popt,pcov = curve_fit(f_exp,hr_sph_pl_nm,Fr_int_sph_pl,p0=[A0_init,l_init,c_init])

f,axs1=fh.pfig()

axs1.semilogx(hr_sph_pl_nm/lD_nm,Fr_int_sph_pl,'ok',label='sim')
axs1.semilogx(hi/lD_nm,f_exp(hi,popt[0],popt[1],popt[2]),'-r',label='fit: $F_0+Ae^{-\kappa h}$')
axs1.autoscale(tight=True)
axs1.legend(loc=1,frameon=False) 
fh.cxlbl(axs1,'h / \lambda_D','F_{sph-pl}/(k_B T)')
axs1.set_title(t_str)    
f.savefig(figFolder+'Fsph_pl_fit.pdf')

#--- predict F energy based on single wall interactions
def h2(h1,geom_nm):
    return geom_nm['L']-2*geom_nm['a']-h1

#Fr_sph_pl = sci.interp1d(hr_sph_pl_nm,Fr_int_sph_pl,kind='quadratic',bounds_error=False)
def Fr_sph_pl(h):
    return f_exp(h,popt[0],popt[1],popt[2])
    
hi = np.linspace(np.min(hr_sph_pl_nm),np.max(hr_sph_pl_nm))

def Fr_pl_sph_pl(Fr_sph_pl,h,geom_nm):
    return Fr_sph_pl(h)+Fr_sph_pl(h2(h,geom_nm))-Fr_sph+Fr_VO_pl_pl(geom_nm['L'])

def x2h(x,geom_nm):
    return x+geom_nm['L']*0.5-geom_nm['a']

def xc(h,geom_nm):
    return h-geom_nm['L']*0.5+geom_nm['a']

def hrc(h,geom_nm):
    x = xc(h,geom_nm)
    return np.linspace(-np.max(np.abs(x)),np.max(np.abs(x)))+geom_nm['L']*0.5-geom_nm['a']

def xcc(h,geom_nm):
    return xc(hrc(h,geom_nm),geom_nm)

Fr0 = np.min(Fr_int_gap)

#--- examine model for Fmin
Fmin_fem = np.zeros(len(Lr_nm))
Fmin_LSA = np.zeros(len(Lr_nm))

for q in range(len(Lr_nm)):
    geom_nm['L']=Lr_nm[q]
    Fmin_fem[q] = np.min(Fr_int_gap[q,:])
    Fmin_LSA[q] = Fr_pl_sph_pl(Fr_sph_pl,x2h(0,geom_nm),geom_nm)


#--- calc integrals over particle position
Ni = 1e4 
q_offs = 3
intPp = np.zeros(len(Lr_nm))
Pshallow_Fint = np.zeros(len(Lr_nm)-q_offs)
Pshallow_Fmin = np.zeros(len(Lr_nm)-q_offs)

#-- calc integrals
for q in range(len(Lr_nm)):
    hc_nm = np.linspace(np.min(hr_nm[q]),np.max(hr_nm[q]),1e4)
    Fc = np.interp(hc_nm,hr_nm[q],Fr_int_gap[q,:])

    Ph = (hc_nm[1]-hc_nm[0])*np.exp(-(Fc-Fr0)) # first order-ish accurate integration
    intPp[q] = sum(Ph)

#--- calc probability ratios
for q in range(len(Lr_nm)-q_offs):
    Pshallow_Fint[q] = 1./(1+(np.exp(-Fr_VO_pl_pl(Lr_nm[q]))*intPp[q+q_offs])\
                        /(np.exp(-Fr_VO_pl_pl(Lr_nm[q+q_offs]))*intPp[q]) )
    Pshallow_Fmin[q] = np.exp(-np.min(Fr_int_gap[q,:])+Fr0-Fr_VO_pl_pl(Lr_nm[q+q_offs]))\
                /(np.exp(-np.min(Fr_int_gap[q,:])+Fr0-Fr_VO_pl_pl(Lr_nm[q+q_offs]))\
                +np.exp(-np.min(Fr_int_gap[q+q_offs,:])+Fr0-Fr_VO_pl_pl(Lr_nm[q])))


#========= plot (some of the) curves
N_plot = 8
if len(Lr_nm) > N_plot:
    I_plot = np.linspace(0,len(Lr_nm)-1,N_plot).astype(np.int64)
else:
    I_plot = list(range(len(Lr_nm)))
        
        
f2,axs2 = fh.pfig()        

for q in I_plot:

    geom_nm['L']=Lr_nm[q]
    col = fh.chot(q,len(Lr_nm))  
    
    lF,lA = 'L-2a: %.1f$\lambda_D$'%((Lr_nm[q]-2*geom_nm['a'])/lD_nm),''
    lF,lA = '',''
    if q==0:
        lF,lA = 'fenics','sum'
    
    #--- constant offset (all curves)
    axs2.plot(xc(hr_nm[q],geom_nm)/lD_nm,Fr_int_gap[q,:],'o',color=col,\
            label=lF)
    axs2.plot(xcc(hr_nm[q],geom_nm)/lD_nm,Fr_pl_sph_pl(Fr_sph_pl,hrc(hr_nm[q],geom_nm),geom_nm),'-',\
                color=col,label=lA)


axs2.set_title('%.1f$<(L-2a)/\lambda_D<$%.1f'%(np.min(Lr_nm-2*geom_nm['a'])/lD_nm,np.max(Lr_nm-2*geom_nm['a'])/lD_nm))
plt.draw()
axs2.legend(loc='best',frameon=False)
fh.cxlbl(axs2,'x/ \lambda_D','F_{pl-sph-pl} / (k_B T) ')
f2.savefig(figFolder+'move_part_test_LSA.pdf')


#axs2.set_title(r'sphere-plane-sphere: '+t_str+'\n'+\
#    'line: $F_{sph-pl}(h)+F_{sph-pl}(L-2a-h)-F_{sph}+F_{pl-pl}$')    
#axs2.legend(bbox_to_anchor=(1.0, 1))

#f3,axs3 = fh.pfig()
#
#axs3.semilogx(Lr_nm/lD_nm,Fmin_fem,'ok',label='sim')
#axs3.semilogx(Lr_nm/lD_nm,Fmin_LSA,'-r',label='LSA')
#axs3.set_xlabel('$L/\lambda_D$')
#axs3.set_ylabel('$\min F_{pl-sph-pl} / (k_B T)$')
#axs3.autoscale(tight=True)
#axs3.legend()
#
#plt.draw()

#--- plot

#--- occupancy probability fig
#fP,axsP = plt.subplots(nrows=2,ncols=1,sharex=True)
#
#fh.make_cfigt()
f1,a1 = fh.pfig()
f2,a2=fh.pfig()

axsP = [a1,a2]


Lr_P = (Lr_nm[q_offs:len(Lr_nm)]-2*geom_nm['a'])/lD_nm
dL_rel = (Lr_nm[q_offs]-Lr_nm[0])/Lr_nm[q_offs]

axsP[0].semilogx(Lr_P,Pshallow_Fint,'-k',label='$f(\int \ e^{-F} \ dh)$ ')
axsP[0].semilogx(Lr_P,Pshallow_Fmin,'-r',label='$f(e^{-F_{min}})$')
#axsP[0].semilogx(Lr_P,(1-dL_rel)/(2-dL_rel)*np.ones(np.shape(Lr_P)),'--k',\
#            label=r'$ \displaystyle \lim_{L \rightarrow \infty} \int \ e^{-F} \ dh$ ')
#axsP[0].semilogx(Lr_P,0.5*np.ones(np.shape(Lr_P)),'--r',\
#            label=r'0.5')
axsP[0].semilogx(Lr_P,(1-dL_rel)/(2-dL_rel)*np.ones(np.shape(Lr_P)),'--k',\
            label=r'')
axsP[0].semilogx(Lr_P,0.5*np.ones(np.shape(Lr_P)),'--r',\
            label=r'')

dL_str = '$\Delta L/L_{deep}$: %.2f'%(dL_rel)   
axsP[0].set_title(dL_str)
axsP[0].legend(loc=4,frameon=False)
fh.cxlbl(axsP[0],'(L_{shallow}-2a)/\lambda_D','P(\mathrm{shallow})')
f1.savefig(figFolder+'p_as_f_L_shallow.pdf')

I = (Pshallow_Fint>0.0001)


axsP[1].semilogx(Lr_P[I],Pshallow_Fint[I]/Pshallow_Fmin[I],'-k')
axsP[1].semilogx(Lr_P,2*(1-dL_rel)/(2-dL_rel)*np.ones(np.shape(Lr_P)),'--k',\
        label=r'$\displaystyle \lim_{L\rightarrow \infty}$')
axsP[1].legend(loc=2,frameon=False)
axsP[1].set_title(dL_str)
fh.cxlbl(axsP[1],'(L_{shallow}-2a)/\lambda_D',r'P(\mathrm{shallow}): \frac{f(\int \ e^{-F} \ dh)}{f(e^{-F_{min}})}')        
f2.savefig(figFolder+'p_as_f_L_shallow_ratio.pdf')
#axsP[1].set_xlabel('$(L_{shallow}-2a)/\lambda_D$')    
#axsP[1].set_ylabel('



#fh.place_fig(f0)
#fh.place_fig(fP)

plt.draw()
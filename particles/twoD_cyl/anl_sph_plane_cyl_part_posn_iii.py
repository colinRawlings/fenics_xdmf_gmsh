# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 18:10:44 2016

@author: cra

analysis of center position for particle as a f of gap distance

"""

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import fenics_helpers as fh
import particle_helpers as ph
import six
from six.moves import range

exec(compile(open('/home/cra/dolfin/useful_f/init_cra.py').read(), '/home/cra/dolfin/useful_f/init_cra.py', 'exec'))
from scipy.optimize import fsolve

#---- unpack data
#fName_save = r'/home/cra/dolfin/particles/twoD_cyl/vary_L_h_symQ_er1_symh.pickle'
#fName_save = r'/home/cra/dolfin/particles/twoD_cyl/vary_L_h_AsymQ_er1_symh.pickle'
#fName_save = r'/home/cra/dolfin/particles/twoD_cyl/vary_L_h_adapt_02.pickle'
fName_save = r'/home/cra/dolfin/particles/twoD_cyl/20160620-123304vary_L_h_asymQ_realEr.pickle'
#fName_save = r'/home/cra/dolfin/particles/twoD_cyl/20160620-134726vary_L_h_symQ_realEr.pickle'


figFolder = '/home/cra/Documents/pres/sph_cyl/figs/'

Rl = pickle.load(open(fName_save,'rb'))

for key,val in six.iteritems(Rl):
    exec(key+r'=val' )
lD_nm = geom_ND['lD_nm']
#----- 
t_str = r'$a/\lambda_D$: %.1f, $\bar\phi_{sph}$: %.1f, $\bar\phi_{pl1}$: %.1f, $\bar\phi_{pl2}$: %.1f'\
        %(geom_nm['a']/lD_nm,params_ND['psi_sph'],params_ND['psi_Si'],params_ND['psi_poly'])
        
        #        %(geom_nm['a']/lD_nm,params_ND['psi_sph'],params_ND['psi_Si'],params_ND['psi_poly'])#+\
        #'$\bar\phi_{pl2}$: %.1f'%(params_ND['psi_poly'])


#---------- helper functions
def x2h(x,geom_nm):
    return x+geom_nm['L']*0.5-geom_nm['a']

def xc(h,geom_nm):
    return h-geom_nm['L']*0.5+geom_nm['a']

def hrc(h,geom_nm):
    x = xc(h,geom_nm)
    return np.linspace(-np.max(np.abs(x)),np.max(np.abs(x)))+geom_nm['L']*0.5-geom_nm['a']

def xcc(h,geom_nm):
    return xc(hrc(h,geom_nm),geom_nm)
#/-----------------

plt.close('all')
f0,axs = plt.subplots(nrows=1,ncols=2)
axs=np.reshape(axs,np.prod(np.shape(axs))) 
fh.make_cfigw()


#==== calc curves
Ni = 1e3
poly_order = 4
kT_lim_init = 1.5 # when making fit
kT_roots = 0.5 # when identifying center
y_min = np.zeros(len(Lr_nm))
y_hiT = y_min.copy()
y_loT = y_min.copy()
yci_r = list()
Fci_r = list()
P_r = list()




#--- do fit
for q in range(len(Lr_nm)):
    
    #    fit = ph.polyfit_minima(hr_nm[q]+geom_nm['a'],Fr_int_gap[q,:],kT_lim_init,kT_roots)
    RES = ph.eval_x_thresh_x_vs_F(hr_nm[q]+geom_nm['a'],Fr_int_gap[q,:])    
    
    y_loT[q],y_hiT[q],y_min[q] = RES['x_lo'],RES['x_hi'],RES['x_min']
    yci_r.append(RES['xc'])
    Fci_r.append(RES['Fc'])        
    
    print(('%d/%d'%(q+1,len(Lr_nm))))
    
    

#========= plot (some of the) curves
N_plot = 20
if len(Lr_nm) > N_plot:
    I_plot = np.linspace(0,len(Lr_nm)-1,N_plot).astype(np.int64)
else:
    I_plot = list(range(len(Lr_nm)))
    
    
for q in I_plot:

    geom_nm['L']=Lr_nm[q]
    col = fh.chot(q,len(Lr_nm))  
    
    axs[0].plot(xc(hr_nm[q],geom_nm)/lD_nm,Fr_int_gap[q,:],'o',color=col,\
            label='L-2a: %.1f$\lambda_D$'%((Lr_nm[q]-2*geom_nm['a'])/lD_nm))
    axs[0].plot(xc(yci_r[q]-geom_nm['a'],geom_nm)/lD_nm,Fci_r[q],'-',color=col)

    

axs[0].set_xlabel(r'$x/ \lambda_D$')
axs[0].set_ylabel(r'$F_{pl-sph-pl} / (k_B T) $')

axs[0].set_title(r'sphere-plane-sphere: '+t_str+'\n'+\
    'line: $F_{sph-pl}(h)+F_{sph-pl}(L-2a-h)-F_{sph}+F_{pl-pl}$')    
axs[0].legend(bbox_to_anchor=(1.0, 1))

axs[1].plot(Lr_nm/lD_nm,y_min/lD_nm,'-b',label='argmin $F$')    
axs[1].plot([0,np.max(Lr_nm/lD_nm)],[0,np.max(Lr_nm/lD_nm)],'-k',label='$L$')    
axs[1].plot(Lr_nm/lD_nm,0.5*Lr_nm/lD_nm,':r',label='$L/2$')
axs[1].plot(Lr_nm/lD_nm,y_hiT/lD_nm,'--b',label='$y: F- \min F = %.1f k_B T$'%(kT_roots))    
axs[1].plot(Lr_nm/lD_nm,y_loT/lD_nm,'--b')    

axs[1].set_xlabel('$L/\lambda_D$')
axs[1].set_ylabel('$y/\lambda_D$')
axs[1].set_ylim([0, np.max(Lr_nm/lD_nm)])
axs[1].set_xlim([0, np.max(Lr_nm/lD_nm)])
axs[1].legend(loc=2)


#--- final BS fig
f2,axs2=fh.pfig()
axs2.plot(Lr_nm/lD_nm,(y_min-Lr_nm*0.5)/lD_nm,'-k',label='min')
axs2.plot(Lr_nm/lD_nm,(y_hiT-Lr_nm*0.5)/lD_nm,'--k',label='5\%-95\%')
axs2.plot(Lr_nm/lD_nm,(y_loT-Lr_nm*0.5)/lD_nm,'--k',label='')
axs2.legend(loc=2,frameon=False)
fh.cxlbl(axs2,r'L/\lambda_D',r'x/\lambda_D')

#f2.savefig(figFolder+'gap_posn_symQ_real_er.pdf')
#f2.savefig(figFolder+'gap_posn_asymQ_real_er_zoom.pdf')

fh.place_fig(f0)
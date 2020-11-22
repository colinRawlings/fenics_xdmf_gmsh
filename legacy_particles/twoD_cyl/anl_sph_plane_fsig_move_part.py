# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 12:59:55 2016

@author: cra
"""

#--- import explicitly to help rope
import numpy as np
import fenics_helpers as fh
import particle_helpers as ph

execfile('/home/cra/dolfin/useful_f/init_cra.py')

#--- multipage pdf output
figFolder = '/home/cra/Documents/pres/sph_cyl/figs/'

#fName = fName_save # last results
#fName = '/home/cra/dolfin/particles/twoD_cyl/move_part_fsig20160616-181529.pickle'
#fName = '/home/cra/dolfin/particles/twoD_cyl/move_part_fsig20160617-152626_erSph1.pickle'
fName = '/home/cra/dolfin/particles/twoD_cyl/move_part_fsig20160620-200933_sym_Q_inf_std_er.pickle'

#--- unpack variables
R = pickle.load(open(fName,'rb'))
for key,val in R.iteritems():
    exec(key+r'=val' )
    

kT_lim_init = 1.5 # when making fit
kT_roots = 2 # when identifying center
    
#--- x vars
xc_ND = np.linspace(0,geom_ND['w'],len(hr_nm[0])*10)    
x_lD = x_ND*(geom_ND['l0_m']/(lD_nm*1e-9))
xc_lD = xc_ND*(geom_ND['l0_m']/(lD_nm*1e-9))

def h2z_lD(h_nm):
    return (h_nm+geom_nm['a'])/lD_nm

def h2y_lD(h_nm,L_nm):
    return (L_nm*0.5-h_nm-geom_nm['a'])/lD_nm


sigma_ND_inf,phi_ND_inf,phi_inf_mV = ph.solve_isol_pl_Qreg(pH,surf_ND,lD_nm,Tdeg)    

#--- fit 
h_min = np.zeros(len(Lr_nm))
h_hiT = h_min.copy()
h_loT = h_min.copy()
hci_r,Fci_r = list(),list()

sim='s_reg'

for q in range(len(Lr_nm)):
    
    RES = ph.eval_x_thresh_x_vs_F(hr_nm[q],F[sim]['surf'][q,:]+F[sim]['vol'][q,:]\
         ,doPlots=False)       
    
    h_loT[q],h_hiT[q],h_min[q] = RES['x_lo'],RES['x_hi'],RES['x_min']
    hci_r.append(RES['xc'])
    Fci_r.append(RES['Fc'])   

    print('%d/%d'%(q+1,len(Lr_nm)))

Nplot = 4
fF,axs1 = fh.pfig()

#--- anl 
for p in range(len(Lr_nm)):

    fS,axs0 = fh.pfig()

    
    #--- sigma curves
    #--- calc Iplot
    Iplot=range(len(hr_nm[p]))    
    if len(hr_nm[p]) > Nplot:
        Iplot = np.asarray(np.linspace(0,len(hr_nm[p])-1,Nplot),dtype=int)
    
    for q in range(len(Iplot)):
        col = fh.chot(q,len(Iplot))
        sci = fh.spline_con_gradient(x_ND,sigma_NDrr[p][Iplot[q],:],[],[],k=3)                 
        
        axs0.plot(x_lD,sigma_NDrr[p][Iplot[q],:]/surf_ND['Gamma'],'+',color=col)
        axs0.plot(xc_lD,sci(xc_ND)/surf_ND['Gamma'],'-',\
            color=col,label='h: %.1f$\lambda_D$'%(hr_nm[p][Iplot[q]]/lD_nm))
    
    axs0.plot(x_lD,sigma_ND_inf/surf_ND['Gamma']*np.ones(np.shape(x_lD)),\
                '--b',label=r'$L \rightarrow \infty$')
    axs0.set_title('$a$: %.1f$\lambda_D$, $L$: %.1f$\lambda_D$'%(geom_nm['a']/lD_nm,Lr_nm[p]/lD_nm))
    axs0.legend(frameon=False)      
    fh.cxlbl(axs0,'x/\lambda_D','\sigma / (e \Gamma)')        
      
    #--- all f curves
    axs1.plot(h2y_lD(hr_nm[p],Lr_nm[p]),F['s_reg']['surf'][p,:]+F['s_reg']['vol'][p,:],'o',\
        color=fh.chot(p,len(Lr_nm)),label='$L: %.1f \lambda_D$'%(Lr_nm[p]/lD_nm))
    axs1.plot(h2y_lD(hci_r[p],Lr_nm[p]),Fci_r[p],'-',color=fh.chot(p,len(Lr_nm)))
    fh.cxlbl(axs1,'y/\lambda_D','(F-F_{min})/k_B T')      
    

#--- final BS fig
f2,axs2=fh.pfig()
axs2.plot(Lr_nm/lD_nm,(geom_nm['a']+h_min-Lr_nm*0.5)/lD_nm,'-k',label='min')
axs2.plot(Lr_nm/lD_nm,(geom_nm['a']+h_hiT-Lr_nm*0.5)/lD_nm,'--k',label='5\%-95\%')
axs2.plot(Lr_nm/lD_nm,(geom_nm['a']+h_loT-Lr_nm*0.5)/lD_nm,'--k',label='')
axs2.legend(loc=2,frameon=False)
fh.cxlbl(axs2,r'L/\lambda_D',r'x/\lambda_D')    

    
    
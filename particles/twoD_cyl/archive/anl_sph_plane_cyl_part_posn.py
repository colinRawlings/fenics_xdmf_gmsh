# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 18:10:44 2016

@author: cra

analysis of center position for particle as a f of gap distance

"""

execfile('/home/cra/dolfin/useful_f/init_cra.py')


#---- unpack data
#fName_save = r'/home/cra/dolfin/particles/twoD_cyl/vary_L_h_symQ_er1_symh.pickle'
#fName_save = r'/home/cra/dolfin/particles/twoD_cyl/vary_L_h_AsymQ_er1_symh.pickle'
fName_save = r'/home/cra/dolfin/particles/twoD_cyl/vary_L_h_adapt_02.pickle'


flds = ['geom_ND','params_ND','geom_nm','Lr_nm','hr_nm','Fr_int_gap',\
            'Fr_int_sph_pl','Fr_sph','hr_sph_pl_nm','Fr_int_pl_pl']
Rl = pickle.load(open(fName_save,'rb'))

for key,val in Rl.iteritems():
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
kT_lim = 3
y_min = np.zeros(len(Lr_nm))
y_plusKT = y_min.copy()
y_minusKT = y_min.copy()
yci_r = list()
Fci_r = list()


for q in range(len(Lr_nm)):
    
    yi = hr_nm[q]+geom_nm['a']    
    Fi = Fr_int_gap[q,:]
    
    #--- interpolate up to kT_lim
    ind_min = (Fi < np.min(Fi)+kT_lim)
    yci = np.linspace(np.min(yi[ind_min]),np.max(yi[ind_min]),Ni)
    yci_r.append(yci)    
            
    P = np.poly1d(np.polyfit(yi[ind_min],Fi[ind_min],2))
    Fci_r.append(P(yci))
    
    
    Fci = si.interp1d(yi,Fi,kind='quadratic')
    
#    ind_min = np.argmin(Fci(yci))
#    
#    y_min[q]=yci[ind_min]
#    
#    ind_subKT = (yci >y_min[q]) & (Fci(yci) < (np.min(Fci(yci))+kT_lim))    
#    y_plusKT[q]=np.max(yci[ind_subKT])
#    
#    ind_subKT = (yci <y_min[q]) & (Fci(yci) < (np.min(Fci(yci))+kT_lim))    
#    y_minusKT[q]=np.min(yci[ind_subKT])
    
    
#    y_plusKT[q]=np.max(yci[ind_subKT])
    
    
axs[1].plot(Lr_nm/lD_nm,y_min/lD_nm,'-b',label='argmin $F$')    
axs[1].plot([0,np.max(Lr_nm/lD_nm)],[0,np.max(Lr_nm/lD_nm)],'-k',label='$L$')    
axs[1].plot(Lr_nm/lD_nm,0.5*Lr_nm/lD_nm,':r',label='$L/2$')
axs[1].plot(Lr_nm/lD_nm,y_plusKT/lD_nm,'--b',label='$F- \min F < %.1f k_B T$'%(kT_lim))    
axs[1].plot(Lr_nm/lD_nm,y_minusKT/lD_nm,'--b')    

axs[1].set_xlabel('$L/\lambda_D$')
axs[1].set_ylabel('$y/\lambda_D$')
axs[1].set_ylim([0, np.max(Lr_nm/lD_nm)])
axs[1].set_xlim([0, np.max(Lr_nm/lD_nm)])
axs[1].legend(loc=2)

#========= plot (some of the) curves
N_plot = 8
if len(Lr_nm) > N_plot:
    I_plot = np.linspace(0,len(Lr_nm)-1,N_plot).astype(np.int64)
else:
    I_plot = range(len(Lr_nm))
    
    
for q in I_plot:

    geom_nm['L']=Lr_nm[q]
    col = fh.chot(q,len(Lr_nm))  
    
    #--- constant offset (all curves)
    axs[0].plot(xc(hr_nm[q],geom_nm)/lD_nm,Fr_int_gap[q,:],'o',color=col,\
            label='L-2a: %.1f$\lambda_D$'%((Lr_nm[q]-2*geom_nm['a'])/lD_nm))


    

axs[0].set_xlabel(r'$x/ \lambda_D$')
axs[0].set_ylabel(r'$F_{pl-sph-pl} / (k_B T) $')

axs[0].set_title(r'sphere-plane-sphere: '+t_str+'\n'+\
    'line: $F_{sph-pl}(h)+F_{sph-pl}(L-2a-h)-F_{sph}+F_{pl-pl}$')    
axs[0].legend(bbox_to_anchor=(1.0, 1))
#axs[3].set_visible(False)

fh.place_fig(f0)
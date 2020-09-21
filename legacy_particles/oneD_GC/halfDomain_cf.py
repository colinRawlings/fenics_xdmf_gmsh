# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 11:36:47 2015

@author: colin
"""

# utility functions for working with fenics and learning

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import dolfin as fn


######## MM
def th_to_M_2D(V,th,cell_no_mag,subdomains):
    # take a fenics function th and use it to compute
    # a new function M representing the magnetisation
    # V is the function space for the components so
    # that the function space of M is V x V
   
    Mxc = fn.project(fn.cos(th),V)
    Myc = fn.project(fn.sin(th),V)
    
    
    Mxdc = fn.interpolate(fn.Constant(0),V)
    Mydc = fn.interpolate(fn.Constant(0),V)
    
    sub_of_cell = subdomains.array()
    for cell_no in range(len(sub_of_cell)):
        if sub_of_cell[cell_no] != cell_no_mag:
            continue
        
        # else grab dofs
        dofs = V.dofmap().cell_dofs(cell_no)
        
        for p in dofs:
            Mxdc.vector()[p] = np.asscalar(Mxc.vector()[p])
            Mydc.vector()[p] = np.asscalar(Myc.vector()[p])
            
        
    # place in vector field
    mxr = Mxdc.vector().array()
    myr = Mydc.vector().array()
    
    mxym = np.array([mxr,myr])
    mxyr = np.squeeze(np.reshape(np.transpose(mxym),(2*np.alen(mxr),1)))
    
    Mxydc = fn.Function(V*V)
    Mxydc.vector()[:] = mxyr    
    
    return Mxydc
    
def XX(u_list,p0=[0,0,0],p1=[1,0,0],N=100):
    
    for pr in [p0,p1]:
        if np.alen(pr) < 3:
            for dim in range(np.alen(pr),3):
                pr.append(0)
        
    if not isinstance(u_list,list):
        u_list = [u_list]
        
    
    xr = np.linspace(p0[0],p1[0],N)
    yr = np.linspace(p0[1],p1[1],N)
    
    zr = np.linspace(p0[2],p1[2],N)
    
    ur = xr.copy()    
    sr = xr.copy()
    
    for u in u_list:
            
        for p in range(N):
            sr[p] = np.sqrt((xr[p]-xr[0])**2+(yr[p]-yr[0])**2+(zr[p]-zr[0])**2)
            if u.geometric_dimension()==1:
                ur[p] = u(xr[p])
            elif u.geometric_dimension()==2:
                ur[p] = u(xr[p],yr[p])
            elif u.geometric_dimension()==3:            
                ur[p] = u(xr[p],yr[p],zr[p])
                
        plt.plot(sr,ur,'-')
        
    
def is_sd_meshFunc_to_Form(sd_num,sd_meshFunc,mesh):
# convert a meshfunction used to mark subdomains into a form representation
# DG 0 elements are used for the representation, 

    V0 = fn.FunctionSpace(mesh,'DG',0) # for containing subdomain numbers
    sd = fn.Function(V0)
    
    sd_nums = np.unique(sd_meshFunc.array())
    
    sd_vals = np.zeros(np.shape(sd_nums))
    ind = plt.mlab.find(sd_nums==sd_num)    
    sd_vals[ind] = 1    
    
    dm = V0.dofmap()
    dof_r = np.zeros(np.shape(sd.vector().array()), dtype=np.int32) # init
    for cell in fn.cells(mesh):
       dof_r[dm.cell_dofs(cell.index())] = sd_meshFunc[cell] # copy cell's mf value to the elements of the dof array corresponding to that cell
    
    sd.vector()[:]=np.choose(dof_r,sd_vals)
    
    return sd
    
def sd_meshFunc_to_Form(sd_meshFunc,mesh):
# convert a meshfunction used to mark subdomains into a form representation
# DG 0 elements are used for the representation

    V0 = fn.FunctionSpace(mesh,'DG',0) # for containing subdomain numbers
    sd = fn.Function(V0)
    
    sd_nums = np.unique(sd_meshFunc.array())
    dm = V0.dofmap()
    dof_r = np.zeros(np.shape(sd.vector().array()), dtype=np.int32) # init
    for cell in fn.cells(mesh):
       dof_r[dm.cell_dofs(cell.index())] = sd_meshFunc[cell] # copy cell's mf value to the elements of the dof array corresponding to that cell
    
    sd.vector()[:]=np.choose(dof_r,sd_nums)
    
    return sd

######## figure helpers when using iPython with Qt backend
def place_fig(posn=(2818, 110, 640, 545)):
    #somewhere where it is not hidden by spyder
    mngr = plt.get_current_fig_manager()
    mngr.window.setGeometry(posn[0],posn[1],posn[2],posn[3])
    
    plt.draw()    
    
    mngr.full_screen_toggle() # ugly way to bring figure to front
    mngr.full_screen_toggle()
  
  
def raise_fig():
    mngr = plt.get_current_fig_manager()
    
    plt.draw()    
    
    mngr.full_screen_toggle() # ugly way to bring figure to front
    mngr.full_screen_toggle()
    
    
def cfigw(screen=2):
    plt.figure()
    
    mngr = plt.get_current_fig_manager()
    geom = mngr.window.geometry()
    x,y,dx,dy = geom.getRect()
    
    if screen==2:
        x = 1957
    else:
        x = 37
    
    
    place_fig(posn=(x,y,2*dx,dy))
    
  

########## element plotting
def plot_sub_els_2D(mesh,subdomains,ax=[],colors=[]):
    
    ## handle inputs
    if np.alen(ax)==0:
        plt.figure()
        ax=plt.subplot(111)
        
    if len(colors)==0:
        colors = std_colors()
        
    subs = np.unique(subdomains.array())
    N_subs = np.alen(subs)
    
    if N_subs>len(colors):
        raise ValueError('More subdomains than supplied colors')

    
    cellList = []
    for p in range(N_subs):
        I = matplotlib.mlab.find(subs[p] == subdomains.array())
        cellList.append(I)
        
    
    # display patches
    desc=[]
    for p in range(len(subs)):
        subList = cellList[p]
        subNum = subs[p] 
        
        plot_els_2D(mesh,subList,color=colors[subNum],ax=ax)
            
        desc_str = '%d=%s' %(subNum,colors[subNum])
        
        if len(desc) == 0:
            desc = 'subdomains: ' + desc_str
        else:
            desc = desc + ', ' + desc_str
            
    plt.title(desc)
    raise_fig()
    

def std_colors():
    return ['yellow','red','blue','green','cyan','magenta','yellow','black']          
            

def plot_els_2D(mesh,el_num_list,color='blue',ax=[]):

    ## handle inputs
    if np.alen(ax)==0:
        plt.figure()
        ax=plt.subplot(111)    
    
    coords = mesh.coordinates()
    cells = mesh.cells()    
    
    if np.alen(el_num_list)==0:
        return
        
    # initialise
    verts = cells[el_num_list[0],:]
    pts = coords[verts,:]
        
    # add together points
    for p in range(1,np.alen(el_num_list)):
        verts = cells[el_num_list[p],:]
        xy = coords[verts,:]   
        
        pts = np.append(pts,[[None, None]],axis=0)
        pts = np.append(pts,xy,axis=0)
        
    ax.add_patch(matplotlib.patches.Polygon(pts,facecolor=color))
    ax.plot(pts[0,0],pts[0,1],'-k')
    
    raise_fig()

        
def plot_el_2D(ax,mesh,el_num,color='blue'):
    coords = mesh.coordinates()
    cells = mesh.cells()    
    
    verts = cells[el_num,:]
    xy = coords[verts,:]
        
    ax.add_patch(matplotlib.patches.Polygon(xy,facecolor=color,edgecolor='Black'))
    
    ax.plot(xy[0,0],xy[0,1],'-k')
    plt.draw()  # draw line otherwise polygon doesn't seem to draw
   
   
   
# mesh inspection
def plot_cell_bnd(ax,cell_nums,mesh_obj,style,color=[]):
    coords = mesh_obj.coordinates()
    cells = mesh_obj.cells()
    x=np.zeros(4)
    y=np.zeros(4)   
    
    for p in np.arange(len(cell_nums)):
        verts = cells[cell_nums[p],:]
    
        for p in range(3):
            x[p]=coords[verts[p],0]
            y[p]=coords[verts[p],1]
        
        x[3]=coords[verts[0],0]
        y[3]=coords[verts[0],1]
        
        if len(color)==0:
            ax.plot(x,y,style)
        else:
            ax.plot(x,y,style,color=color)
        
    raise_fig()
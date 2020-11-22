# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 11:36:47 2015

@author: colin
"""

# utility functions for working with fenics and learning

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import dolfin as fn
import os
import pickle
import six
from six.moves import range


def dict_vals_to_Const(a_dict,use_dolfin_adjoint=False):
    ''' recursively convert the val of all keys in a dict which are floats
    or ints into fn.Constant (s) ready for use in a form
    
    formerly in particle_helpers module
    
    '''
    
    import dolfin as fn
    if use_dolfin_adjoint:
        import dolfin_adjoint as da
    
    new_dict = {}
    for key,val in six.iteritems(a_dict):
        if (isinstance(val,float) or isinstance(val,int)):
            if use_dolfin_adjoint:
                new_dict[key] = da.Constant(val)
            else:
                new_dict[key] = fn.Constant(val)
        elif isinstance(val,dict):
            new_dict[key] = dict_vals_to_Const(val)
        else:
            new_dict[key]=val
                       
    return new_dict   

##### time steps from modified gryphon code
def load_gryph_t_steps(folder):
    '''
    
    =======================================================================
    THE MODIFIED ESDIRK WAS LOST AND WILL HAVE TO BE REWRITTEN.  The dev scripts
    are save at /home/cra/dolfin/therm/heat_only/dev_HDF5_.*
    
    gryphon appears compatible with fenics v2016.1
    =======================================================================
    
    
    
    returns the time steps for the results stored in folder following a call to
    gryphon with 
    
    obj.parameters['output']['save_steps'] = True

    see also load_gryph_t_step which for a given time step (integer not the 
    value of t at that time step) returns the dolfin Function object for that 
    time step    
    
    following the call: folder = obj.savepath
    
    it may be preferrable to use the equivalent methods in the object:
        ESDIRK.load_t_steps
        ESDIRK.load_t_step
    
    '''

    fName_meta = folder+'/time_steps/meta_data.pickle'
    
    if not os.path.isfile(fName_meta):
        raise ValueError('the folder {} doesn"t contain any meta data'.format(fName_meta))
     
    meta = pickle.load(open(fName_meta,'rb'))
    
    return meta['t']


def load_gryph_t_step(folder,time_step,u=[]):
    
    '''
    =======================================================================
    THE MODIFIED ESDIRK WAS LOST AND WILL HAVE TO BE REWRITTEN.  The dev scripts
    are save at /home/cra/dolfin/therm/heat_only/dev_HDF5_.*
    
    gryphon appears compatible with fenics v2016.1
    =======================================================================
    
    
    returns the dolfin Function object for the requested time_step from the 
    results created following a call to gryphon with:
    
    obj.parameters['output']['save_steps'] = True

    see also load_gryph_t_steps which returns the list of time steps stored in 
    the folder
    
    following the call: folder = obj.savepath

    it may be preferrable to use the equivalent methods in the object:
        ESDIRK.load_t_steps
        ESDIRK.load_t_step    
    
    '''    
    
    fName_h5 = folder+'/time_steps/results.h5'
    fName_meta = folder+'/time_steps/meta_data.pickle'
    
    for f in [fName_h5,fName_meta]:
        if not os.path.isfile(f):
            raise ValueError('the folder {} doesn"t contain any meta data'.format(f))
    
    meta_o = pickle.load(open(fName_meta,'rb'))    
    
    fo = fn.HDF5File(fn.mpi_comm_world(),fName_h5,'r')
    mesho = fn.Mesh()
    fo.read(mesho,'mesh',False)
    
    # config the function space
    Vo = fn.FunctionSpace(mesho,meta_o['mesh']['type'],meta_o['mesh']['order'])
    uo = fn.Function(Vo)
    
    fo.read(uo,meta_o['fmt_str'].format(time_step))
    
    return uo

####### mesh scripting
def init():
    exec(compile(open('/home/cra/dolfin/useful_f/init_cra.py').read(), '/home/cra/dolfin/useful_f/init_cra.py', 'exec'))
    


def run_cmd(cmd,quiet=True,fName_out = '/tmp/py_proc_stdout',fName_err = '/tmp/py_proc_stderr'):
    # parts of cmd should be specified in the standard subprocess style 
    # i.e. with each part making up an element of a list
        
    import subprocess as sp
           
    fout = open(fName_out,'w')
    ferr = open(fName_err,'w')
    
    status = sp.call(cmd,stdout=fout,stderr=ferr)
    
    fout.close()
    ferr.close()
    
    # print output if there was a problem
    if (not status==0):
        print('Process returned with non zero exit status: '+'%d'%(status))
    
    if (not status==0) or (not quiet):
        print('stdout (saved in: '+fName_out+'):')
        sp.call(['less','-f',fName_out])
        print('-----')
        print('stderr (saved in: '+fName_err+'):')
        sp.call(['less','-f',fName_err])
    
    return status


def paraview_fld(u_list,meshName='tmp',Verbose=False,apply_stamp=False,pv_script=None,
                 refine_mesh=False,mesh=None, ref_num=None):
    '''
    display a solution

    :param u_list: list of solutions to export
    :param meshName: defaults to tmp
    :param Verbose: give op name (primarily for debug)
    :param apply_stamp: apply stamp, important if there is a chance of collision at load
    time in paraview
    :param ref_num: a number to use to tag the exported data.  Useful if you want to display
    several components of a solution.  This way they can be overwritten each time the script
    is run but will not collide and will be loaded into separate paraview sessions.
    :param pv_script: script to run after starting paraview
    :param refine_mesh: refine the mesh before projecting, this makes sense if the CG2
    elements are used since paraview uses CG1
    :param mesh: must be supplied if the mesh is to be refined ready for display
    :return: None
    '''

    # utility function for opening a solution in paraview
    
    import os    
    import random
    import dolfin as fn

    if not len(os.path.dirname(meshName)):
        meshName = os.getcwd()+ '/' + meshName
    
    if not isinstance(u_list,list):
        u_list = [u_list]        
    
    if apply_stamp:
        meshName = meshName+'{:.0f}'.format(random.random()*1000)

    if ref_num is not None:
        meshName = meshName + '{:.0f}'.format(ref_num)



    p=0
    for u in u_list:
        fName = meshName+'_{:.0f}.pvd'.format(p)
        if refine_mesh:
            if mesh is None:
                raise ValueError('mesh must be supplied for output of refined solution')

            meshr = fn.refine(mesh)
            El = fn.FiniteElement('CG', meshr.ufl_cell(), degree=1) # element used by paraview
            Vr = fn.FunctionSpace(meshr, El)

            u = fn.project(u,Vr)

        file = fn.File(fName)
        file << u
        if p == 0:
            fName_all=fName
        else:
            fName_all = fName_all+', '+fName
        p+=1    
    
    
    if Verbose:
        print(('Paraview file for fld inspection written to: '+fName_all))
    
    # --- build cmd
    if p==1:
        res=os.system('/home/cra/apps/ParaView-5.1.2-Qt4-OpenGL2-MPI-Linux-64bit/bin/paraview --data='+meshName+'_0.pvd'+' &')
    else:
        res=os.system('/home/cra/apps/ParaView-5.1.2-Qt4-OpenGL2-MPI-Linux-64bit/bin/paraview --data='+meshName+'_..pvd'+' &')
    
    if Verbose:
        print(res)
    
    

def paraview_mesh(mesh,sd_MF=[],meshName='tmp'):
    # utility function for viewing a mesh in paraview    

    V = fn.FunctionSpace(mesh,'Lagrange',1)
    u = fn.Function(V)
    
    if not sd_MF==[]:
        # project the subdomain numbers
        u = sd_meshFunc_to_Form(sd_MF,mesh) 
            
    paraview_fld(u,meshName=meshName)
    

def line_appender(filename,line):
    with open(filename, 'a') as f:
        f.write(line+'\n')
    

def line_prepender(filename, line):
    # from SO 5914627
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)


def params_to_mesh(params,fName_tail,fName_out,verbose=False,clear_files=True):
    
    '''
    create a dolfin mesh from a gmsh file (fName_tail) with the parameters 
    commented out and a dictionary containing those parameters and their
    desired values
    
    returns mesh,sd_MF,bnd_MF
    
    where sd_MF and bnd_MF are the subdomain and boundary mesh functions as
    defined for the physical surfaces in the gmsh description file
    '''    
    
    
    import shutil as sh
    import os.path      
    import os
    import random
            
    gmsh_bin = 'gmsh'

    def p2m_print(str):
        if verbose:
            print(str)
            
            
    # default to saving in same place as fName_tail
    #-- get full path for tail
    if not len(os.path.dirname(fName_tail)):
        fName_tail = os.getcwd()+ '/' + fName_tail
    
    if fName_out==[]:
        # add a random number to the timestamp to avoid collisions
        rand_tag = '{:.0f}'.format(random.random()*1000)
        fName_out='tmp'+tstamp()+'_'+rand_tag
        if not clear_files:
            print(('geometry files assigned tag: {}'.format(rand_tag)))
    
    
    if not len(os.path.dirname(fName_out)):
        dirSave = os.path.dirname(fName_tail)
        fName_out = dirSave+'/'+fName_out
    
    fName_parts = os.path.splitext(fName_out)
    fName_geo = fName_parts[0]+'.geo'
    
    # make space for output files
    extList = ['.geo','.msh','.xml']
    for ext in extList:
        fName = fName_parts[0]+ext
        if os.path.isfile(fName):
            os.remove(fName)
        
    # build geometry description file   
    sh.copyfile(fName_tail,fName_geo)
    
    for par,val in six.iteritems(params):
        line = par+' = DefineNumber[ '+ '%.12e'%(val)+', Name "Parameters/' + par + '" ];'
        line_prepender(fName_geo,line)
        p2m_print(line)
    p2m_print('Geometry file written to: '+fName_geo)
    
    # create mesh
    p2m_print('--> Generating mesh ...')
    run_cmd([gmsh_bin,fName_geo,'-3','-smooth','5'])
    p2m_print('Mesh written to: '+fName_parts[0]+'.msh')
    
    p2m_print('--> converting mesh to dolfin format ...')
    run_cmd(['dolfin-convert',fName_parts[0]+'.msh',fName_parts[0]+'.xml'])
    meshName = fName_parts[0]
    p2m_print('xml file(s) written to: ' + fName_parts[0]+'.xml')
    
    
    # load mesh
    mesh = fn.Mesh('%s.xml'%(meshName))
    fName_sd_MF = "%s_physical_region.xml"%(meshName)
    fName_sd_bnd = "%s_facet_region.xml"%(meshName)
    
    sd_MF = []
    if os.path.isfile(fName_sd_MF):
        sd_MF = fn.MeshFunction("size_t", mesh, fName_sd_MF)
        
    bnd_MF = []
    if os.path.isfile(fName_sd_bnd):
        bnd_MF = fn.MeshFunction("size_t", mesh, fName_sd_bnd)
    
    if clear_files:
        ext=['.geo','.msh','.xml','_facet_region.xml','_physical_region.xml']
        for ex in ext:
            os.remove(fName_parts[0]+ex)
    
    return mesh,sd_MF,bnd_MF
    
    

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
    
def calcXX(u_list,p0,p1,N=100):
    ''' 
    Calculate the values of the expressions in u_list at the points between
    p0 and p1.  Returned are the values of these expressions as well as the 
    a dictionary containing different options for the independent variables IVs:
    IVs['x0','x1',...,'xN','s'], where N is the geometry of the model, s is the
    Euclidean distance along the XX.
    
    return DVs,IVs
    
    '''
    
    import numpy as np    
    
    #--- check inputs
    if not isinstance(u_list,list):
        u_list = [u_list]    
    
    if not len(p0)==len(p1):
        raise ValueError('Inconsistent point dimensions')
        
    if not len(p0)==u_list[0].geometric_dimension():
        raise ValueError('point dimensions does not match form dimension')
        
    #--- calc independent variables        
    IVs = dict() 
    for p in range(len(p0)):
        IVs['x%d'%(p)] = np.linspace(p0[p],p1[p],N)
        
    s = np.zeros((N))
    for p in range(N):
        dd = 0
        for q in range(len(p0)):
            dd = dd+(IVs['x%d'%(q)][p]-IVs['x%d'%(q)][0])**2
        s[p]=np.sqrt(dd)
    IVs['s']=s
        
    #--- calc list of dependent variables
    DVs = list()
    
    for l in range(len(u_list)):
        #--- loop over all points
        arr = np.zeros(N)                    
        for p in range(N):
            #-- calc point
    
            pEval = tuple()                    
            for q in range(len(p0)):
                pEval = pEval + (IVs['x%d'%(q)][p],)
            
            #- eval
            arr[p] = u_list[l](pEval)
            
        DVs.append(arr)
        
    return DVs,IVs
        
        
    
def plotXX(u_list,p0,p1,N=100,axs=[],style='-',IV='s',label='',color=(1,0,0,1)):
    ''' plot cross sections for the forms in u_list between the specified
    start and end locations p0=(x,y,...) and p1 
    
    IV (independent variable) may be: 'x0', 'x1', ... 'xN' or 
    the default 's': the euclidean distance along the XX    
    '''    

    import numpy as np    
    
    
    #--- check inputs
    if not isinstance(u_list,list):
        u_list = [u_list]    
        
    if not isinstance(axs,list):
        axs = [axs] 
    
    if not len(p0)==len(p1):
        raise ValueError('Inconsistent point dimensions')
        
    if not len(p0)==u_list[0].geometric_dimension():
        raise ValueError('point dimensions does not match form dimension')
        
    if (not len(axs)==0) and (not len(axs)==len(u_list)):
        raise ValueError('number of supplied axes doesn''t match number of forms')
    
    DVs,IVs = calcXX(u_list,p0,p1,N=N) 
    
    #--- plottinge
    if len(axs)==0:
        f0,axs=plt.subplots(nrows=1,ncols=len(u_list))
        axs=np.reshape(axs,np.prod(np.shape(axs))) # convert to column vector of handles     
    
        if len(axs) > 1:
            make_cfigw()
        place_fig(f0)
    
    
    for p in range(len(DVs)):
        axs[p].plot(IVs[IV],DVs[p],style,label=label,color=color)
        axs[p].set_xlabel(IV)
        
    # raise_fig(axs[0])
    return DVs,IVs
    
def is_sd_meshFunc_to_Form(sd_num,sd_meshFunc,mesh):
    '''
    convert a meshfunction used to mark subdomains into a form representation
    DG 0 elements are used for the representation, 
    
    is_sd = is_sd_meshFunc_to_Form(sd_num,sd_meshFunc,mesh)
    '''


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
    '''
    convert a meshfunction used to mark subdomains into a form representation
    DG 0 elements are used for the representation
    
    sd = sd_meshFunc_to_Form(sd_meshFunc,mesh)
    
    n.b. the sd_MF takes one value per cell  while the Function takes one value 
    per dof, the DOF map relates the cells to the dofs, thus we copy the 
    value of sd_MF[cell] to each DOF in that cell
    '''



    V0 = fn.FunctionSpace(mesh,'DG',0) # for containing subdomain numbers
    sd = fn.Function(V0)
    
    #sd_nums = np.unique(sd_meshFunc.array())
    dm = V0.dofmap()
    dof_r = np.zeros(np.shape(sd.vector().array()), dtype=np.int32) # init
    for cell in fn.cells(mesh):
       dof_r[dm.cell_dofs(cell.index())] = sd_meshFunc[cell] # copy cell's mf value to the elements of the dof array corresponding to that cell
    
    sd.vector()[:]=dof_r#np.choose(dof_r,sd_nums)
    
    return sd

######## figure helpers when using iPython with Qt backend
def make_cfigw():
    '''
    turn the current figure into a wide figure:
    make_cfigw()
    
    '''    
    
    try:
        mngr = plt.get_current_fig_manager()
        g = mngr.window.geometry

        mngr.window.setGeometry(g().left(),g().bottom(),2.0*g().width(),g().height())
    except AttributeError:
        print("Could not resize figure")

######## figure helpers when using iPython with Qt backend
def make_cfigt():
    '''
    turn the current figure into a tall figure:
        
    '''    
    
    mngr = plt.get_current_fig_manager()
    g = mngr.window.geometry

    mngr.window.setGeometry(g().left(),g().bottom(),g().width(),2.0*g().height())

def make_cfigb():
    '''
    turn the current figure into a big figure:
        
    '''    
    
    mngr = plt.get_current_fig_manager()
    g = mngr.window.geometry

    mngr.window.setGeometry(g().left(),g().bottom(),2.0*g().width(),2.0*g().height())

def chot(p,N):
    ''' return a RGBa tuple for use as color=col argument, 
    
    - p is the line index
    - N is the number of lines to use'''
    
    import matplotlib 
    
    useRatio = 0.7 # don't use full colormap range as yellow-white colors hard to see
    
    if p > N:
        raise ValueError ('p (%.1f) should not exceed N (%.1f)')
    
    cm = matplotlib.cm.hot
    
    R = (1.*p)/(1.*N)

    Ract = R*useRatio
    
    return cm(Ract)
    


def place_current_fig(posn=2):
    '''  
    posn = [x,y,Lx,Ly] or 2 for screen 2 
    '''
    
    import six.moves.tkinter as tk
    
    #--- screen size
    root = tk.Tk()
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    #--- current position
    mngr = plt.get_current_fig_manager()    
    LHS = mngr.window.geometry().left()
    bottom = mngr.window.geometry().bottom()
    height = mngr.window.geometry().height()
    width = mngr.window.geometry().width()
    curr_posn = np.asarray([LHS,bottom-height,width,height])    
    
    
    #move to screen two
    set_posn = curr_posn[:]    
    if np.alen(posn)==4:
        set_posn=posn[:]
    
    if np.alen(posn)==1 and posn==2:
        if LHS < 1920 and screen_width > 1920:
            set_posn[0] = set_posn[0]+1920
    
    
    mngr.window.setGeometry(set_posn[0],set_posn[1],set_posn[2],set_posn[3])
    
    plt.draw()    
    
    mngr.full_screen_toggle() # ugly way to bring figure to front
    mngr.full_screen_toggle()
  
  
def place_fig(id,posn=(2)):
    return
    # raise_fig(id)
    # place_current_fig(posn=posn)
  
  
def raise_fig(ip):
    if type(ip) is int:
        plt.figure(ip)
    elif type(ip) is str:
        plt.figure(ip)
    elif type(ip) is matplotlib.figure.Figure:
        plt.figure(ip.number)
    elif issubclass(type(ip), matplotlib.axes.SubplotBase):
        fig = ip.get_figure()
        plt.figure(fig.number)
    else:
        raise ValueError('Argument should be an int, a str (specifying a fig) or an axes or a figure object')
            

    raise_current_fig()  
  
  
def raise_current_fig():

    try:
        mngr = plt.get_current_fig_manager()
        
        plt.draw()    
        plt.show()    
        
        mngr.full_screen_toggle() # ugly way to bring figure to front
        mngr.full_screen_toggle()
    except AttributeError:
        print("could not raise figure")
    
def cfig():
    
    ''' 
    create a normal sized figure
    
    return f0,axs
        
    '''
    
    import matplotlib.pyplot as plt
    import numpy as np 
    
    f0,axs = plt.subplots(nrows=1,ncols=1)
    f0.set_size_inches(6, 6, forward=True)

    return f0,axs
    
def pfig():
    
    ''' 
    create a paper sized figure 
    
    return f0,axs
    
    '''
    
    import matplotlib.pyplot as plt
    import numpy as np 
    
    f0,axs = plt.subplots(nrows=1,ncols=1)
    f0.set_size_inches(4, 4, forward=True)

    return f0,axs
        
def pfigw(ncols=2):
    '''
    create a wide figure for use in a paper
    :param ncols: number of columns in single row of figure to display
    :return: figure, axes
    '''

    import matplotlib.pyplot as plt
    import numpy as np

    f0, axs = plt.subplots(nrows=1, ncols=ncols)

    f0.set_size_inches(16, 6, forward=True)

    # --- ensures axs is a 1D list of axes
    if ncols == 1:
        axs = [axs]
    axs = np.reshape(axs, np.prod(np.shape(axs)))

    return f0, axs


def calc_meas(V, dl):
    """
    Calculate the measure of the dl over V

    :param V: FunctionSpace
    :param dl: Measure
    :return:
    """

    u_const = fn.project(fn.Constant(1.), V)
    return fn.assemble(u_const * dl)


def avg_val(u, V, dl):
    """
    Calculate the average value of the form u \in V on dl

    :param u: Form/Function
    :param V: FunctionSpace
    :param dl: Measure
    :return:
    """

    return fn.assemble(u * dl) / calc_meas(V, dl)

def cfigw(ncols=2):
    '''
    create a wide figure
    '''    
        
    
    import matplotlib.pyplot as plt
    import numpy as np    
    
    f0,axs = plt.subplots(nrows=1,ncols=ncols)

    f0.set_size_inches(16, 6, forward=True)

    #--- ensures axs is a 1D list of axes
    if ncols==1:
        axs=[axs]
    axs=np.reshape(axs,np.prod(np.shape(axs)))

    return f0,axs    
  
    
def cfigb(ncols=2,nrows=2):
    '''
    create a wide figure
    '''    
        
    
    import matplotlib.pyplot as plt
    import numpy as np    
    
    f0,axs = plt.subplots(nrows=nrows,ncols=ncols)

    f0.set_size_inches(18, 12, forward=True)

    # mngr = plt.get_current_fig_manager()
    # geom = mngr.window.geometry()
    # x,y,dx,dy = geom.getRect()
    #
    # place_current_fig(posn=(x+1920,y,2*dx,2*dy))
    
    plt.tight_layout(pad=6.)

    #--- ensures axs is a 1D list of axes
    #    if ncols==1:
    #        axs=[axs]
    #    axs=np.reshape(axs,np.prod(np.shape(axs)))

    return f0,axs     

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
    raise_current_fig()
    

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
    
    raise_current_fig()

        
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
        
    raise_current_fig()
    
    
def op_xy2File(fName,x,y,xfmt='%.5e',yfmt='%.5e',sep=', '):
    f = open(fName,'w')
    
    pat = xfmt+sep+yfmt+'\n'
    for p in range(np.alen(x)):
        str  = pat%(x[p],y[p])
        f.write(str)
    
    f.close()
        
def op_mat2File(fName,A,fmt='%.5e',sep=', '):
    file = open(fName,'w')

    (rows,cols) = np.shape(A)
    
    for p in range(rows):
        for q in range(cols):
            if q==0:
                str = fmt%(A[p,q])  
            else:
                str = str+sep+fmt%(A[p,q])
        
        file.write(str+'\n')
    
    
    file.close()
        
        
    
def solve_2D_osc_plate(params,Re,meshName_tail):
    # calculate the real and imaginary parts of the velocity flow field for the
    # non dimensional 2D plate flow problem
    #
    # R.['ur','ui','pr','pi'] = solve_2D_osc_plate()   
    
    #--- vars
    mesh,sd_MF,bnd_MF=params_to_mesh(params,meshName_tail,'out.geo')

    # Define function spaces
    V = fn.VectorFunctionSpace(mesh, "CG", 2, dim=4) # Trial: pair of velocity fields, Test: Cauchy Tensor
    Q = fn.VectorFunctionSpace(mesh, "CG", 1, dim=2) # Trial: pair of pressure fields, Test: Divergence condition
    W = V * Q
    S = fn.FunctionSpace(mesh,'CG',1) # for post only
    VV = fn.VectorFunctionSpace(mesh,'CG',1) # for post only
    
    # No-slip boundary condition for velocity
    noslip = fn.Constant((0., 0., 0., 0.))
    osc = fn.Expression(("0", "1","0.", "0.0"))
    freeS_p = fn.Constant((0,0))
    
    
    def is_origin(x,on_boundary):
        return x[0] < 10*fn.DOLFIN_EPS and x[1] < 10*fn.DOLFIN_EPS
    
    bc1 = fn.DirichletBC(W.sub(0), noslip, bnd_MF, 3)
    bc2 = fn.DirichletBC(W.sub(0), noslip, bnd_MF, 2)
    bc3 = fn.DirichletBC(W.sub(0), osc, bnd_MF, 1)
    bcp = fn.DirichletBC(W.sub(1),freeS_p,is_origin,'pointwise')
    
    # Collect boundary conditions
    bcs = [bc1,bc2,bc3,bcp]
    
    
    # Define variational problem
    (un, pn)  = fn.TrialFunctions(W)
    (vn, qn) = fn.TestFunctions(W)
    
    ur = fn.as_vector([un[0],un[1]])
    ui = fn.as_vector([un[2],un[3]])
    pr = pn[0]
    pi = pn[1]
    
    vr = fn.as_vector([vn[0],vn[1]])
    vi = fn.as_vector([vn[2],vn[3]])
    qr = qn[0]
    qi = qn[1]
    
    fr = fn.Constant((0,0))
    fi = fn.Constant((0,0))
    
    (j,l) = fn.indices(2)
    
    a = fn.Dx(ur[l],j)*fn.Dx(vr[l],j)*fn.dx-fn.Dx(vr[j],j)*pr*fn.dx + qr*fn.Dx(ur[j],j)*fn.dx+\
        fn.Dx(ui[l],j)*fn.Dx(vi[l],j)*fn.dx-fn.Dx(vi[j],j)*pi*fn.dx + qi*fn.Dx(ui[j],j)*fn.dx-\
        fn.Constant(Re)*ui[j]*vr[j]*fn.dx+fn.Constant(Re)*ur[j]*vi[j]*fn.dx
    
    L = fr[j]*vr[j]*fn.dx+fi[j]*vi[j]*fn.dx
    
    # Compute solution
    w = fn.Function(W)
    problem = fn.LinearVariationalProblem(a,L,w,bcs=bcs)
    solver = fn.LinearVariationalSolver(problem)
    
    solver.parameters['preconditioner']='ilu'
    solver.parameters['lu_solver']['report'] = True
    solver.parameters['linear_solver']='mumps'
    solver.parameters['krylov_solver']['monitor_convergence'] = True
    solver.parameters['krylov_solver']['maximum_iterations'] = 100
    solver.solve()
    
    # results: w = [uxy,uyr,uxi,uyi,pr,pi]
    return w
    

def outer_legend(ax):
    ''' 
    move an axis'' legend outside the axis (from SE4700614)

    assumes labels are already defined
    
    '''
        
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                     
                     
                 
                 
                 
    
def anl_2D_osc_plate(w,params,Re,fName_out=[],axs={'ReP':[],'ImP':[],'ReUy':[]},Npts=200):
    # display and optionally write to file the results of a 2D_osc_plate calculation:
    # w = solve_2D_osc_plate()  where w = [uxy,uyr,uxi,uyi,pr,pi]
    #     
    # 
    #   
    
    import numpy as np     
    import fenics_helpers as fh    
    
    #--- initialise figures
    for par,val in six.iteritems(axs):
        if np.alen(val)==0:
            f,axs[par] = plt.subplots()
            fh.place_fig(axs[par])
        
    
        
    #---- calc arrays    
    xp = np.linspace(0,1.1,Npts)
    Dpi_num = np.zeros(np.shape(xp))
    Dpr_num = np.zeros(np.shape(xp))
    
    xu = np.linspace(1+1e-5,1.4,Npts)
    uyr_num = np.zeros(np.shape(xu))
    
    (un,pn) = w.split()  #w = [uxy,uyr,uxi,uyi,pr,pi]
    
    for q in range(Npts):
        Dpi_num[q] = (pn(params['L']*0.5+xp[q],params['L']*0.5+params['t']*0.5+1e-4)[1]-\
                      pn(params['L']*0.5+xp[q],params['L']*0.5-params['t']*0.5-1e-4)[1])
        Dpr_num[q] = (pn(params['L']*0.5+xp[q],params['L']*0.5+params['t']*0.5+1e-4)[0]-\
                      pn(params['L']*0.5+xp[q],params['L']*0.5-params['t']*0.5-1e-4)[0])
        uyr_num[q] = un(params['L']*0.5+xu[q]+1e-4,params['L']*0.5)[1]
    
    #--- pressure
    axs['ReP'].plot(xp,Dpr_num,'-')    
    axs['ReP'].set_xlim(0,1.1)
    axs['ReP'].set_xlabel('$x/a$')
    axs['ReP'].set_ylabel('$\mathcal{R}(P)$')
    
    place_fig(axs['ReP'])
        
        
    axs['ImP'].plot(xp,Dpi_num,'-')
    axs['ImP'].set_xlim(0,1.1)
    axs['ImP'].set_xlabel('$x/a$')
    axs['ImP'].set_ylabel('$\mathcal{I}(p)$')
    
    place_fig(axs['ImP'])
    
    #-- velocity
    axs['ReUy'].plot(xu,uyr_num,'-')
    axs['ReUy'].set_xlabel('$x/a$')
    axs['ReUy'].set_ylabel('$\mathcal{R}(u_y)$')
    
    place_fig(axs['ReUy'])
        
    
    if np.alen(fName_out)==0:
        fName_out = 'Re=%.3e_'%(Re)
        
    op_xy2File(fName_out+'Im_P.txt',xp,Dpi_num)           
    op_xy2File(fName_out+'Real_P.txt',xp,Dpr_num)           
    op_xy2File(fName_out+'Real_Uy.txt',xu,uyr_num)    


def spline_con_gradient(x, y,x_con,Dy_con, k=3, s=0, w=None):
    '''
    fit a spline to the points in x,y subject to the constraints that 
    dy/dx|_{x = x_con[p]} = Dy_con[p] 
    
    n.b. k is the spline order which the documentation recommends to set to 3

    '''    
    
    import numpy as np
    from scipy.interpolate import UnivariateSpline, splev, splrep
    from scipy.optimize import minimize    
    
        #--- helpers 
    def guess(x, y, k, s, w=None):
        """Do an ordinary spline fit to provide knots"""
        return splrep(x, y, w, k=k, s=s)
    
    def err(c, x, y, t, k, w=None):
        """The error function to minimize"""
        diff = y - splev(x, (t, c, k))
        if w is None:
            diff = np.einsum('...i,...i', diff, diff)
        else:
            diff = np.dot(diff*diff, w)
        return np.abs(diff)    
    
    #--- build constraints
    def Dspline_err(c,t,k,x,Dy):
        ''' evaluate the derivative of the spline at x'''
        return splev(x, (t, c, k), der=1)-Dy
    
    #--- init
    t, c0, k = guess(x, y, k, s, w=w) 
    
    #--- constraints
    cons = list()
    for p in range(len(x_con)):
        cons.append({'type': 'eq','fun':Dspline_err,'args': (t,k,x_con[p],Dy_con[p])})
    cons = tuple(cons)

    #--- solve    
    opt = minimize(err, c0, (x, y, t, k, w), constraints=cons,method='SLSQP')
    copt = opt.x
    return UnivariateSpline._from_tck((t, copt, k))
    
#--- pretty plot production    
def cxlbl(ax,x_str,y_str):
    '''
    typeset x and y labels in latex mathmode
    '''    
    
    import matplotlib.pyplot as plt    
    
    ax.set_xlabel('$'+x_str+'$')    
    ax.set_ylabel('$'+y_str+'$')
    
    tight_ax(ax)
    
def cxulbl(ax,xv_str,xu_str,yv_str,yu_str):
    '''
    typeset x and y labels in latex mathmode with units
    '''    
    
    import matplotlib.pyplot as plt
    
    if len(xu_str)>0:
        ax.set_xlabel('$'+xv_str+'\quad (\mathrm{'+xu_str+'})$' )    
    else:
        ax.set_xlabel('$'+xv_str+'$' )    

    if len(yu_str)>0:        
        ax.set_ylabel('$'+yv_str+'\quad (\mathrm{'+yu_str+'})$' ) 
    else:
        ax.set_ylabel('$'+yv_str+'$' ) 
    
    tight_ax(ax)
    
def less_xticks(a):
    """
    :param a: the axis to modify
    """
    
    
    a.set_xticks(a.get_xticks()[::2])
    tight_ax(a)
        
def less_yticks(a):
    """
    :param a: the axis to modify
    """
    
    
    a.set_yticks(a.get_yticks()[::2])
    tight_ax(a)
    
def set_yticks(a,spacing):
    
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    
    a.yaxis.set_major_locator(ticker.MultipleLocator(spacing))    
    tight_ax(a)
    
def set_xticks(a,spacing):
    
    import matplotlib.pyplot as plt
    
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    
    a.xaxis.set_major_locator(ticker.MultipleLocator(spacing))    
    tight_ax(a)
    
    
    

def tight_ax(ax):
    import matplotlib.pyplot as plt    
    
    ax.figure.tight_layout() # ensure visible 
    plt.pause(0.001)
    ax.figure.canvas.draw()
    
def ctitle(ax,str):

    ax.set_title(str)    
    tight_ax(ax)
    
def tstamp():
    ''' 
    return a time stamp for use when constructing the name of a save file
    '''

    import time

    return  time.strftime("%Y%m%d-%H%M%S")   
    

    
#---- fitting of material parameters
def f_poly_gauss_np(u,v,Npoly):
    '''
    generic material fitting function
    
    a gaussian plus an Npoly order polynomial
    
    evaluates the current guess for the fitting parameters in the 
    dictionary v.  Uses numpy (i.e. for numerical eval)
        
    return f_np
    
    '''
    
    import numpy as np
    
    f_poly = v['N0']*np.ones(np.shape(u))
    for pm1 in range(Npoly):
        f_poly += u**(pm1+1)*v['N{:.0f}'.format(pm1+1)]

    return f_poly + v['Df_gauss']*np.exp(-(u-v['u_gauss'])**2/(2*v['sigma']**2))    

    
def f_poly_gauss_form(T,v,Npoly):
    '''
    generic material fitting function
    
    a gaussian plus an Npoly order polynomial for use when specifying a form
    
    evaluates the current guess for the fitting parameters in the 
    dictionary v.  Uses dolfin function
    
    return poly_gauss_form
    
    '''
    
    import dolfin as fn
    
    rho_poly = v['N0']*fn.Constant(1)
    for pm1 in range(Npoly):
        rho_poly += T**(pm1+1)*v['N{:.0f}'.format(pm1+1)]

    return (rho_poly + v['Df_gauss']*fn.exp(-(T-v['u_gauss'])**2/(2*v['sigma']**2)))
    
    
    
def fit_poly_gauss(u0,u1,f,Npoly,doPlots=False,Npts = 100,ax_plot=None):
    '''
    fit the function fu on the interval u: [u0,u1] using the poly_gauss 
    function.  fu should take a single argument which should be a numpy array.
    It can be defined via lambda if required.
    
    
    returns the dictionary of required parameters.
    
    return v
    '''
    
    import lmfit
    import numpy as np
        
    ur = np.linspace(u0,u1,Npts)
    fr = f(ur)
    
    #--- config fit
    p = lmfit.Parameters()
    u_gauss = ur[np.argmax(fr)] 
    Df_gauss = np.max(fr)-np.min(fr)
    sigma = 0.1*(np.max(ur)-ur[0])
    for q in range(Npoly+1):
        p.add('N{:.0f}'.format(q),value=0)
    p.add_many(('u_gauss',u_gauss),('Df_gauss',Df_gauss),('sigma',sigma))
    
    def residual(p):
        v = p.valuesdict()
        return f_poly_gauss_np(ur,v,Npoly)-fr
    
    mi = lmfit.minimize(residual,p,method='leastsq')
    v = mi.params.valuesdict()
    
    v['Npoly']=Npoly
    
    if doPlots:
        if ax_plot==None:
            f,ax_plot=cfig()
        fig,ax_plot=cfig()
        ax_plot.plot(ur,fr,'-k',label='exact')
        ax_plot.plot(ur,f_poly_gauss_np(ur,v,Npoly),'--g',label='fitted')
        ax_plot.legend()
        cxlbl(ax_plot,'u_r','f_r')

    return v
    
def save_var(fName,var):
    '''

    uses pickle to save var at fName


    :param fName: save filename
    :param var: variable to save
    :return: None
    '''

    
    import pickle
    import os
    
    
    folder = os.path.dirname(fName)
    if not os.path.isdir(folder):
        os.mkdir(folder)
    
    with open(fName, 'w') as f:
        pickle.dump(var, f)
        
    print(('wrote file: {}'.format(fName)))
    

def load_var(fName):
    ''' 
    helper/wrapper to load the variable pickled at fName
    '''
        
    import pickle

    return pickle.load(open(fName,'rb'))        
        
    
#--- dictionary access methods
def getFromDict(dataDict, mapList):  
    '''
    extract a value from a nested dictionary 
    
    val = dataDict['b']['v']['y']
    
    by defining mapList=['b','v','y']
    
    and calling val = getFromDict(dataDict,mapList)
    
    from SO14692690
    
    '''
    
    
    for k in mapList: dataDict = dataDict[k]
    return dataDict

    
def setInDict(dataDict, mapList, value):
    
    '''
    set a value in a nested dictionary 
    
    dataDict['b']['v']['y'] = value
    
    by defining mapList=['b','v','y']
    
    and calling setInDict(dataDict,mapList,value)
    
    from SO14692690
    
    '''
    
    getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value    
    
    

    
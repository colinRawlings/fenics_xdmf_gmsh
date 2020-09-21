# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 13:54:09 2015

@author: cra
"""


def calc_sph_surf_pot(params, R_balloon=4, R_mesh=0.05, doPlot=False):
    #    params={'a':a,'sigma',sigma,'n0',n0,'er',er}
    # solve for surface potential for CRA ND groups in units of KT/e

    import dolfin as fn
    import numpy as np

    k0 = np.sqrt(params['n0'] / params['er'])

    # -- build mesh and functions
    mesh = fn.IntervalMesh(int(np.ceil(R_balloon / R_mesh)), params['a'],
                           params['a'] + R_balloon / k0)
    V = fn.FunctionSpace(mesh, 'CG', 2)

    # -- setup for BC
    bnd_mf = fn.MeshFunction('size_t', mesh, 0)

    class initMF(fn.SubDomain):
        def inside(self, x, on_boundary):
            return True

    class LeftB(fn.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and np.abs(x[0] - params['a']) < 1e2 * fn.DOLFIN_EPS

    init_mf = initMF()
    lb = LeftB()
    init_mf.mark(bnd_mf, 0)
    lb.mark(bnd_mf, 1)
    ds = fn.Measure('ds')[bnd_mf]

    # -- specify problem
    u_ = fn.Function(V)
    u = fn.TrialFunction(V)
    v = fn.TestFunction(V)

    # -- solve problem
    u_ = fn.Function(V)
    u = fn.TrialFunction(V)
    v = fn.TestFunction(V)

    r = fn.SpatialCoordinate(mesh)

    F = r ** 2 * fn.Constant(params['sigma'] / params['er']) * v * ds(1) - \
        r ** 2 * fn.Dx(u, 0) * fn.Dx(v, 0) * fn.dx - \
        2 * fn.Constant(params['n0'] / params['er']) * r ** 2 * fn.sinh(u) * v * fn.dx

    F = fn.action(F, u_)
    J = fn.derivative(F, u_, u)

    problem = fn.NonlinearVariationalProblem(F, u_, [], J)
    solver = fn.NonlinearVariationalSolver(problem)
    solver.solve()

    if doPlot:
        fn.plot(u_)

    return u_(params['a'])


def calc_sph_charge(psi_s, params, R_balloon=4, R_mesh=0.05, doPlot=False, log=False):
    #    params={'a':a,'n0',n0,'er',er}
    # n.b. these are the ND params 
    #
    # find a surface charge to achieve the specified psi_s

    import dolfin as fn
    from scipy.optimize import minimize_scalar

    # --- search opt
    def f(s):
        params['sigma'] = s

        psi_surf = calc_sph_surf_pot(params, R_balloon=R_balloon, R_mesh=R_mesh, \
                                     doPlot=False)

        return (psi_surf - psi_s) ** 2

    fn.set_log_active(log)
    res = minimize_scalar(f)
    fn.set_log_active(True)

    return res.x


def calc_sphere_curves(ratchet, geom, Npts=200):
    #  calculate the curves directly above and below the sphere for both the sphere 
    # itself and its bounding curves

    import numpy as np

    th = np.linspace(0, np.pi, Npts)
    ym0 = np.interp(geom['x0'], ratchet['x'], ratchet['middle'])

    # --- calc for sphere
    x = geom['x0'] + np.cos(th) * geom['a']
    y = {}
    y['x'] = x
    y['top'] = ym0 + np.sin(th) * geom['a']
    y['bottom'] = ym0 - np.sin(th) * geom['a']

    y['top_rat'] = np.interp(x, ratchet['x'], ratchet['top'])
    y['bottom_rat'] = np.interp(x, ratchet['x'], ratchet['bottom'])

    # --- repeat for mesh
    yMesh = {}
    x = geom['x0'] + np.cos(th) * (geom['a'] + geom['De'])
    yMesh['x'] = x
    yMesh['top'] = ym0 + np.sin(th) * (geom['a'] + geom['De'])
    yMesh['bottom'] = ym0 - np.sin(th) * (geom['a'] + geom['De'])

    yMesh['top_rat'] = np.interp(x, ratchet['x'], ratchet['top_mesh'])
    yMesh['bottom_rat'] = np.interp(x, ratchet['x'], ratchet['bottom_mesh'])

    return y, yMesh


def load_const():
    ''' returns a dict of useful physical constants'''

    const = {'e': 1.602e-19, 'e0': 8.85e-12, 'kB': 1.38e-23, \
             'Na': 6.022e23}

    return const


def calc_ratchet_profile(geom_nm, Npts=200, doPlots=False, ax=[]):
    """
     calculate the shape of the ratchet and the densely meshed regions 
     surrounding it    
    
     return ratchet_nm,ySphMesh,ySph      
    """

    import numpy as np
    import matplotlib.pyplot as plt
    import fenics_helpers as fh

    # -- calc ratchet
    xc = np.linspace(0, geom_nm['w'] * 2, Npts)
    ytc = geom_nm['L'] * np.ones(np.shape(xc))

    xbi = np.cumsum([0, geom_nm['w'] / 2, geom_nm['p'], geom_nm['w'] - geom_nm['p'], \
                     geom_nm['p'], geom_nm['w'] / 2 - geom_nm['p']])
    ybi = geom_nm['D'] * np.asarray([0, 0.5, -0.5, 0.5, -0.5, 0])

    ybc = np.interp(xc, xbi, ybi)
    ymc = np.interp(xc, xbi, 0.5 * ybi + 0.5 * geom_nm['L'])

    # --- calc mesh surfaces surrounding charge layers
    ytcM = ytc - geom_nm['De']

    xbiM = xbi.copy()
    xbiM[1:5] = xbiM[1:5] + geom_nm['De']
    ybiM = ybi.copy() + geom_nm['De']
    ybcM = np.interp(xc, xbiM, ybiM)

    ratchet_nm = {'x': xc, 'top': ytc, 'bottom': ybc, 'middle': ymc, \
                  'top_mesh': ytcM, 'bottom_mesh': ybcM}

    ySph, ySphMesh = calc_sphere_curves(ratchet_nm, geom_nm, Npts=Npts)

    # --- plot
    if doPlots:
        if np.alen(ax) == 0:
            plt.figure()
            ax = plt.subplot(111)

        ax.plot(xbi, ybi, 'ok')
        ax.plot(xc, ybc, '-k')
        ax.plot(xc, ytc, '-k')
        ax.plot(xc, ymc, ':b')
        ax.plot(xc, ytcM, ':b')
        ax.plot(xc, ybcM, ':b')

        ax.plot(ySph['x'], ySph['top'], '-k')
        ax.plot(ySph['x'], ySph['bottom'], '-k')
        ax.plot(ySph['x'], ySph['top_rat'], '-k')
        ax.plot(ySph['x'], ySph['bottom_rat'], '-k')

        ax.plot(ySphMesh['x'], ySphMesh['top'], ':b')
        ax.plot(ySphMesh['x'], ySphMesh['bottom'], ':b')
        ax.plot(ySphMesh['x'], ySphMesh['top_rat'], ':b')
        ax.plot(ySphMesh['x'], ySphMesh['bottom_rat'], ':b')

        ax.set_xlabel('$x$ (nm)')
        ax.set_ylabel('$z$ (nm)')
        ax.axis('equal')
        plt.draw()
        fh.place_fig(ax)

    return ratchet_nm, ySphMesh, ySph


def calc_surf_ND_pars(surf, Tdeg):
    ''' 
        calc the ND surface pars given the dimensional surface pars and the
        relevant scaled par dicts
    '''

    const = load_const()
    Tkelv = Tdeg + 273
    l0_m = calc_l0_m(Tdeg)

    surf_ND = dict()
    surf_ND['pK'] = surf['pK']
    surf_ND['Gamma'] = surf['Gamma'] * l0_m ** 2
    surf_ND['C_stern'] = surf['C_stern'] * l0_m ** 2 * const['kB'] * Tkelv / const[
                                                                                 'e'] ** 2

    return surf_ND


def calc_phi_d_Grahame(sigma_ND, lD_nm, Tdeg, er_fluid):
    '''
        Invert the Grahame equation (ND) to obtain the solution for flat 
        isolated surfaces
        
        supplied sigma should be in ND units (scaled by l0**2/e)        
        
        return value is in ND units (scaled by e/kT)
    '''

    import numpy as np

    n0 = calc_n0_ND(lD_nm, Tdeg)

    return 2 * np.arcsinh(sigma_ND / np.sqrt(8 * n0 * er_fluid))


def calc_phi_d_Grahame_SI(sigma, n0_M, Tdeg, er=None):
    '''
    calculate the diffuse layer surface potential in terms of quantities in SI 
    units
    
    if er is not specified the value for water is used
    
    n0 is the concentration of a single ion species in moles
    
    returns the potential in volts
    '''

    import numpy as np

    const = load_const()

    if er is None:
        er = calc_er_fluid(Tdeg)

    lD_nm = n0_M_to_lD_nm(n0_M, Tdeg)
    Tkelv = Tdeg + 293

    return 2 * const['kB'] * (Tkelv) / const['e'] \
           * np.arcsinh(
        const['e'] * sigma * lD_nm * 1e-9 / (2 * er * const['e0'] * const['kB'] * Tkelv))


def calc_phi_d_SiOH_Beh(sigma_ND, pH, surf_ND):
    ''' 
    calc the relationship between the Stern potential and the surface charge dictated
    by the reaction coefficients for the dissociation of the Silanol groups
    '''

    import numpy as np

    # --- check no positive values for sigma_ND
    if np.sum(np.asarray(sigma_ND > 0, dtype=int)) > 0:
        raise ValueError('sigma_ND cannot exceed 0')

    return -sigma_ND / surf_ND['C_stern'] + (surf_ND['pK'] - pH) * np.log(10) \
           + np.log(-sigma_ND / (surf_ND['Gamma'] + sigma_ND))


def calc_init_opt(geom_nm, doPlots=True):
    # return the center position and acceptable bounds on h (in nm) for a given 
    # geom_nm    

    import numpy as np

    # ---- center position
    ratchet_nm, ySphMesh, ySph = calc_ratchet_profile(geom_nm, doPlots=doPlots)
    yc0 = np.interp(geom_nm['x0'], ratchet_nm['x'], ratchet_nm['middle'])

    hc = geom_nm['L'] - geom_nm['a'] - yc0

    # --- bounds on h
    dH_top = np.min(ySphMesh['top_rat'] - ySphMesh['top'])
    dH_bottom = np.min(ySphMesh['bottom'] - ySphMesh['bottom_rat'])

    h_bounds = np.zeros(2)
    h_bounds[0] = 1.05 * 2 * geom_nm['De']  # flat surface
    h_bounds[1] = 0.95 * (2 * geom_nm['De'] + dH_top + dH_bottom)

    # ensure mid-plane cannot cut the ratchet
    ylim = geom_nm['D'] / 2 + geom_nm['De'] * 1.02
    ymin = geom_nm['L'] - geom_nm['a'] - h_bounds[1]
    if ymin < ylim:
        ymin = ylim
        h_bounds[1] = geom_nm['L'] - geom_nm['a'] - ymin

    return hc, h_bounds


def lD_nm_to_n0_M(lD_nm, Tdeg):
    '''
    convert between a Debye length in nm and a concentration in Moles/l
    '''

    import numpy as np

    const = load_const()

    er_fluid = 87.74 - Tdeg * 0.4008  # dielectric constant for water at 20C: Malmberg, NIST, 1956
    Tkelv = Tdeg + 273

    return er_fluid * const['e0'] * const['kB'] * Tkelv / (
        2 * (lD_nm * 1e-9) ** 2 * const['e'] ** 2) \
           / (1e3 * const['Na'])


def n0_M_to_lD_nm(n0_M, Tdeg):
    '''
    convert between a Debye length in nm and a concentration in Moles/l
    '''

    import numpy as np

    const = load_const()

    er_fluid = 87.74 - Tdeg * 0.4008  # dielectric constant for water at 20C: Malmberg, NIST, 1956
    Tkelv = Tdeg + 273

    n0_SI = n0_M * const['Na'] * 1e3

    return 1e9 * np.sqrt(
        er_fluid * const['e0'] * const['kB'] * Tkelv / (2 * n0_SI * const['e'] ** 2))


def calc_er_fluid(Tdeg):
    '''
    dielectric constant for water at 20C: Malmberg, NIST, 1956
    '''

    return 87.74 - Tdeg * 0.4008 + 9.398e-4 * Tdeg ** 2 - 1.410e-6 * Tdeg ** 3


def calc_l0_m(Tdeg):
    '''
    calc the intrinsic length scale for thermal-electrostatic problems
    '''

    const = load_const()

    return const['e'] ** 2 / (const['e0'] * const['kB'] * (Tdeg + 273))


def calc_n0_ND(lD_nm, Tdeg):
    '''
    calc the non-dimensional concentration corresponding to a given Debye length
    for a 1:1 electrolyte
    '''

    er_fluid = calc_er_fluid(
        Tdeg)  # dielectric constant for water at 20C: Malmberg, NIST, 1956
    l0_m = calc_l0_m(Tdeg)

    n0 = er_fluid / (2 * (lD_nm * 1e-9) ** 2 * l0_m)

    return n0 * l0_m ** 3


def conv_to_ratchet_ND_pars(geom_nm, lD_nm, er_dry, Tdeg, psi_mV):
    '''     
     convert parameters to CRA ND groups for use in ratchet simulation
     
     returns a geom_ND dictionary for mesh generation and a params_ND used to
     define the equation
    
     er_dry contains the dielectric constant for the bounding volumes 
     (Si: upper surface, poly: lower surface, sph: sphere) 
      
     return geom_ND,params_ND
     '''

    import numpy as np

    const = load_const()

    er_fluid = calc_er_fluid(
        Tdeg)  # dielectric constant for water at 20C: Malmberg, NIST, 1956
    Tkelv = Tdeg + 273

    # -- calc scaled lengths
    l0_m = calc_l0_m(Tdeg)

    geom_ND = {}
    for par, val_nm in geom_nm.iteritems():
        geom_ND[par] = val_nm * 1e-9 / l0_m

        # CRA ND
    n0_ND = calc_n0_ND(lD_nm, Tdeg)
    n0 = n0_ND / l0_m ** 3

    k0 = np.sqrt(n0 * const['e'] ** 2 / (er_fluid * const['e0'] * const['kB'] * Tkelv))
    k0_ND = k0 * l0_m

    # calc surface charges
    psi_sph_ND = 1e-3 * psi_mV['sph'] * const['e'] / (const['kB'] * Tkelv)
    sigma_ND = {'sph': calc_sph_charge(psi_sph_ND,
                                       {'a': geom_ND['a'], 'n0': n0_ND, 'er': er_fluid},
                                       log=False)}
    for fld in ['Si', 'poly']:
        sigma_ND[fld] = np.sqrt(8 * n0_ND * er_fluid) * np.sinh(
            1e-3 * psi_mV[fld] / 2 * const['e'] / (
                const['kB'] * Tkelv))  # the Grahame equation

    # - build final params struct
    params_ND = {}
    for fld in ['Si', 'poly', 'sph']:
        params_ND['sigma_' + fld] = sigma_ND[fld]
        params_ND['psi_' + fld] = 1e-3 * psi_mV[fld] * const['e'] / (const['kB'] * Tkelv)

    params_ND['V_therm_mV'] = 1e3 * const['kB'] * Tkelv / const['e']

    if 'fluid' in er_dry:
        raise ValueError('Dielectric constant should not be specified in er_dry')

    params_ND['er_fluid'] = er_fluid
    params_ND['er_dry'] = er_dry
    params_ND['n0'] = n0_ND
    params_ND['k0'] = k0_ND
    params_ND['Tdeg'] = Tdeg
    geom_ND['l0_m'] = l0_m  # ugly!
    geom_ND['lD_nm'] = lD_nm  # ugly but convenient
    geom_ND['lD'] = lD_nm * 1e-9 / l0_m  # this is the Debye length in units of l0!

    return geom_ND, params_ND


def dict_vals_to_Const(a_dict, use_dolfin_adjoint=False):
    ''' recursively convert the val of all keys in a dict which are floats
    or ints into fn.Constant (s) ready for use in a form
    '''

    import fenics_helpers

    print('dict_vals_to_Const has been moved to the fenics_helpers module')

    return fenics_helpers.dict_vals_to_Const(a_dict,
                                             use_dolfin_adjoint=use_dolfin_adjoint)


def solve_isol_pl_Qreg(pH, surf_ND, lD_nm, Tdeg, s0_ND=-1):
    ''' 
    solve for the surface charge for an isolated plane surface, also returns the 
    corresponding surface potential
    '''

    import scipy.optimize

    const = load_const()

    def calc_err(sigma_ND, pH, surf_ND):
        err = calc_phi_d_Grahame(sigma_ND, lD_nm, Tdeg, calc_er_fluid(Tdeg)) \
              - calc_phi_d_SiOH_Beh(sigma_ND, pH, surf_ND)

        return err

    x, infodict, ier, mesg = scipy.optimize.fsolve(lambda s: calc_err(s, pH, surf_ND),
                                                   s0_ND, \
                                                   full_output=True, factor=1e-2)

    sigma_ND = float('nan')
    phi_d = float('nan')
    phi_d_mV = float('nan')
    if ier == 1:
        sigma_ND = x[0]
        phi_d = calc_phi_d_Grahame(sigma_ND, lD_nm, Tdeg, calc_er_fluid(Tdeg))
        phi_d_mV = 1e3 * const['kB'] * (Tdeg + 273) / const['e'] * phi_d
    return sigma_ND, phi_d, phi_d_mV


def calc_cylSym_ratchet_problem(geom_ND, params_ND, SDs, fName_tail, doPlots=True):
    ''' calc the ratchet problem for the case of cylindrical symmetry'''

    import dolfin as fn
    import fenics_helpers as fh

    # -- convert parameters to Constants to avoid unnecessary recompilation
    params_ND = dict_vals_to_Const(params_ND)

    # -- construct basis
    mesh, sd_MF, bnd_MF = fh.params_to_mesh(geom_ND, fName_tail,
                                            'out' + fh.tstamp() + '.geo' \
                                            , verbose=False)
    if doPlots:
        fh.paraview_mesh(mesh, sd_MF=sd_MF)

    V = fn.FunctionSpace(mesh, 'CG', 2)
    u_ = fn.Function(V)
    u = fn.TrialFunction(V)
    v = fn.TestFunction(V)

    # -- define equation
    ds = fn.Measure('dS')[bnd_MF]
    dx = fn.Measure('dx')[sd_MF]
    dx_wet = dx(SDs['vol']['fluid'])
    vol_dry = ['sph', 'poly', 'Si']

    r, z = fn.SpatialCoordinate(mesh)

    # -- fluid terms
    F = v / (r + 1e-6) * fn.Dx(u, 0) * params_ND['er_fluid'] * dx_wet \
        - 2 * params_ND['n0'] * fn.sinh(u) * v * dx_wet \
        - fn.dot(params_ND['er_fluid'] * fn.grad(u), fn.grad(v)) * dx_wet

    # -- add dry terms:derivative+charge
    for p in range(len(vol_dry)):
        F += v / (r + 1e-6) * fn.Dx(u, 0) * params_ND['er_dry'][vol_dry[p]] * dx(
            SDs['vol'][vol_dry[p]]) \
             - fn.dot(params_ND['er_dry'][vol_dry[p]] * fn.grad(u), fn.grad(v)) * dx(
            SDs['vol'][vol_dry[p]]) \
             + params_ND['sigma_' + vol_dry[p]] * v('+') * ds(SDs['surf'][vol_dry[p]])

    F = fn.action(F, u_)
    J = fn.derivative(F, u_, u)

    problem = fn.NonlinearVariationalProblem(F, u_, [], J)

    return problem, u_, mesh, sd_MF, bnd_MF


def calc_cylSim_F_int(params_ND, SDs, bnd_MF, sd_MF, u_, mesh):
    '''
    define the equation to be assembled for the evaluation of the free energy
    in the case of cylindrical symmetry
    
    this should be done after the solution so that any adaptive refinement is 
    taken into account
    
    return f => F = assemble(f)
    
    '''

    import dolfin as fn
    import fenics_helpers as fh
    import numpy as np

    # -- convert parameters to Constants to avoid unnecessary recompilation
    params_ND = dict_vals_to_Const(params_ND)

    # -- define measures
    ds = fn.Measure('dS')[bnd_MF.leaf_node()]
    dx = fn.Measure('dx')[sd_MF.leaf_node()]
    dx_wet = dx(SDs['vol']['fluid'])
    vol_dry = ['sph', 'poly', 'Si']

    r, z = fn.SpatialCoordinate(mesh.leaf_node())

    # - n.b. with |J| added for cylindrical coords
    u_el = params_ND['er_fluid'] * 0.5 * fn.dot(fn.grad(u_.leaf_node()), fn.grad(
        u_.leaf_node())) * 2 * fn.pi * r * dx_wet
    for p in range(len(vol_dry)):
        u_el += params_ND['er_dry'][vol_dry[p]] * 0.5 \
                * fn.dot(fn.grad(u_.leaf_node()),
                         fn.grad(u_.leaf_node())) * 2 * fn.pi * r * dx(
            SDs['vol'][vol_dry[p]])

    Tds = 2 * params_ND['n0'] * (fn.cosh(u_.leaf_node()) - 1.0 - u_.leaf_node() * fn.sinh(
        u_.leaf_node())) * 2 * fn.pi * r * dx_wet

    # --- calc surface charge contribution
    if not ('sigma_ND_inf' in params_ND):
        params_ND['sigma_ND_inf'] = params_ND['sigma_Si']
        params_ND['phi_ND_inf'] = u_.leaf_node()

    Fch = -(params_ND['sigma_Si'] - params_ND['sigma_ND_inf']) * 0.5 * (
        u_.leaf_node() + params_ND['phi_ND_inf']) \
          * fn.Constant(2 * np.pi) * r * ds(SDs['surf']['Si'])

    # --- outputs
    return u_el - Tds, Fch


def calc_ratchet_problem(geom, pars, SDs, meshName_tail, doPlots=True):
    r""" calc the ratchet problem for a given set of ND pars
    
      calc_ratchet_problem(geom,pars,SDs,meshName_tail,doPlots=True)"""

    import dolfin as fn
    import fenics_helpers as fh

    # --- mesh
    mesh, sd_MF, bnd_MF = fh.params_to_mesh(geom, meshName_tail, \
                                            'out' + fh.tstamp() + '.geo', verbose=True,
                                            clear_files=True)

    if doPlots:
        fh.paraview_mesh(mesh, sd_MF=sd_MF)

    # ----- specify problem
    # -- periodic boundary conditions:
    class PeriodicBoundary(fn.SubDomain):

        # Left boundary is "target domain" G
        def inside(self, x, on_boundary):
            return bool(x[0] < fn.DOLFIN_EPS and x[0] > -fn.DOLFIN_EPS and on_boundary)

        # Map right boundary (H) to left boundary (G)
        def map(self, x, y):
            y[0] = x[0] - geom['w'] * 2
            y[1] = x[1]
            y[2] = x[2]

    # Create periodic boundary condition
    pbc = PeriodicBoundary()

    V = fn.FunctionSpace(mesh, 'CG', 2, constrained_domain=pbc)
    u_ = fn.Function(V)
    u = fn.TrialFunction(V)
    v = fn.TestFunction(V)

    ds = fn.Measure('dS')[bnd_MF]
    dx = fn.Measure('dx')[sd_MF]
    dx_wet = dx(SDs['vol']['fluid'])
    vol_dry = ['sph', 'poly', 'Si']

    # -- fluid terms
    F = -2 * pars['n0'] * fn.sinh(u) * v * dx_wet \
        - fn.dot(pars['er_fluid'] * fn.grad(u), fn.grad(v)) * dx_wet

    # -- dry terms
    for p in range(len(vol_dry)):
        F = F - fn.dot(pars['er_dry'][vol_dry[p]] * fn.grad(u), fn.grad(v)) * dx(
            SDs['vol'][vol_dry[p]]) \
            + pars['sigma_' + vol_dry[p]] * v('+') * ds(SDs['surf'][vol_dry[p]])

    F = fn.action(F, u_)
    J = fn.derivative(F, u_, u)

    problem = fn.NonlinearVariationalProblem(F, u_, [], J)

    return problem, u_, mesh, sd_MF, bnd_MF


def calc_F_int(params_ND, SDs, bnd_MF, sd_MF, u_):
    '''
    define the equation to be assembled for the evaluation of the free energy
    
    this should be done after the solution so that any adaptive refinement is 
    taken into account
    
    return f => F = assemble(f)
    
    n.b. no consideration of charge regulation possible here so no Fch returned
    
    '''

    import dolfin as fn

    # -- convert parameters to Constants to avoid unnecessary recompilation
    params_ND = dict_vals_to_Const(params_ND)

    ds = fn.Measure('dS')[bnd_MF.leaf_node()]
    dx = fn.Measure('dx')[sd_MF.leaf_node()]
    dx_wet = dx(SDs['vol']['fluid'])
    vol_dry = ['sph', 'poly', 'Si']

    u_el = params_ND['er_fluid'] * 0.5 * fn.dot(fn.grad(u_.leaf_node()),
                                                fn.grad(u_.leaf_node())) * dx_wet
    for p in range(len(vol_dry)):
        u_el += params_ND['er_dry'][vol_dry[p]] * 0.5 * fn.dot(fn.grad(u_.leaf_node()),
                                                               fn.grad(
                                                                   u_.leaf_node())) * dx(
            SDs['vol'][vol_dry[p]])

    Tds = 2 * params_ND['n0'] * (
        fn.cosh(u_.leaf_node()) - 1.0 - u_.leaf_node() * fn.sinh(u_.leaf_node())) * dx_wet

    return u_el - Tds


def calc_phi_d_PB_fsig(x_ND, sigma_ND, geom_ND, params_ND, SDs, fName_tail, \
                       set_sigma_poly=False, cylSim=True, full_output=False,
                       doPlots=False, \
                       N_adapt=0, adapt_TOL=1e-2, disp_log=False):
    ''' 
    calc the relationship between the Stern potential and the surface charge dictated
    by the PB relationship
    
    the surface charge on the Si and optionally the poly boundary is defined by the 
    relationship sigma_ND(x_ND)
    '''

    import fenics_helpers as fh
    import numpy as np
    import dolfin as fn

    class pyExpression(fn.Expression):
        def __init__(self, sci):
            self.sci = sci

        def eval(self, value, x):
            value[0] = self.sci(x[0])


            # --- fit the spline

    sci = fh.spline_con_gradient(x_ND, sigma_ND, [], [], k=3)
    params_ND['sigma_Si'] = fn.avg(pyExpression(sci))

    if set_sigma_poly:
        params_ND['sigma_poly'] = fn.avg(pyExpression(sci))

    # --- get phi_d
    Fr_int, Fch_int, Fc_int, u_, mesh = solve_ratchet_problem(geom_ND, params_ND, SDs,
                                                              fName_tail, \
                                                              doPlots=doPlots,
                                                              cylSim=cylSim,
                                                              disp_log=disp_log,
                                                              N_adapt=N_adapt,
                                                              adapt_TOL=adapt_TOL)

    phi_d = np.zeros(len(x_ND))
    for p in range(len(x_ND)):
        phi_d[p] = u_(x_ND[p], geom_ND['L'])

    if full_output:
        return phi_d, Fr_int, Fch_int, u_, mesh
    else:
        return phi_d


def polyfit_minima(x, y, Dy_lim_init, Dy_contour, Ni=1e3, poly_order=4):
    '''
    fits a polynomial to the minima of y(x)
    
    uses only points for which y < Dy_lim+min(y)
    
    returns the argmin for the contour as well as the points:
    x = min(y)+Dy_contour
    
    '''

    import numpy as np
    import scipy.optimize

    # --- interpolate up to Dy_lim
    finished = False
    Dy_lim = Dy_lim_init
    while not finished:
        ind_min = (y < np.min(y) + Dy_lim)
        n_pts = sum(np.asarray(ind_min, dtype=int))

        if n_pts > poly_order:
            finished = True
        else:
            Dy_lim = Dy_lim + 0.1 * Dy_lim_init

    xc = np.linspace(np.min(x[ind_min]), np.max(x[ind_min]), Ni)

    P = np.polyfit(x[ind_min], y[ind_min], poly_order)
    fP = np.poly1d(P)

    # --- find 0 and \pm kT
    x_minDy = scipy.optimize.fsolve(fP - (np.min(y) + Dy_contour), xc[0])
    x_plusDy = scipy.optimize.fsolve(fP - (np.min(y) + Dy_contour), xc[len(xc) - 1])
    x_min = .5 * (x_minDy + x_plusDy)

    fit = dict()
    for var in ['xc', 'x_minDy', 'x_plusDy', 'x_min', 'fP']:
        fit[var] = eval(var)
    fit['yc'] = fP(xc)

    return fit


def solve_ratchet_problem(geom, pars, SDs, meshName_tail, doPlots=False, cylSim=False, \
                          adapt=False, N_adapt=0, adapt_TOL=1e-12, disp_log=True,
                          log_level=16):
    '''
    formulate and solve the ratchet problem based on the ND pars and the geometry
    specified by the subdomain listing (SDs) and the geometry file (meshName_tail)
    
    N_adapt>0 activates adaptive mesh refinement 
    cylSim=True turns on cylindrical symmetry (so no really a ratchet anymore)
    
    '''

    import fenics_helpers as fh
    import dolfin as fn
    import numpy as np

    if not fn.__version__ == '1.6.0':
        raise ValueError('wrong dolfin version: {}'.format(fn.__version__))

    fn.set_log_active(disp_log)
    fn.set_log_level(log_level)

    # --- calc problem
    if cylSim:
        problem, u_, mesh, sd_MF, bnd_MF = calc_cylSym_ratchet_problem(geom, pars, SDs,
                                                                       meshName_tail,
                                                                       doPlots=doPlots)
        f, fch = calc_cylSim_F_int(pars, SDs, bnd_MF, sd_MF, u_, mesh)
    else:
        problem, u_, mesh, sd_MF, bnd_MF = calc_ratchet_problem(geom, pars, SDs,
                                                                meshName_tail,
                                                                doPlots=doPlots)
        f = calc_F_int(pars, SDs, bnd_MF, sd_MF, u_)

        # --- solve problem
    if N_adapt == 0:
        solver = fn.NonlinearVariationalSolver(problem)
        solver.parameters['newton_solver']['linear_solver'] = 'gmres'
        solver.parameters['newton_solver']['error_on_nonconvergence'] = False
        solver.solve()
    else:
        solver = fn.AdaptiveNonlinearVariationalSolver(problem, f)
        solver.parameters['max_iterations'] = N_adapt
        solver.parameters['error_control']['dual_variational_solver'][
            'linear_solver'] = 'gmres'
        solver.parameters['nonlinear_variational_solver']['newton_solver'][
            'linear_solver'] = 'gmres'
        solver.solve(adapt_TOL)
    fn.set_log_active(True)

    if doPlots:
        fh.paraview_fld(u_.leaf_node())

        # --- calc free energy
    if cylSim:
        f, fch = calc_cylSim_F_int(pars, SDs, bnd_MF, sd_MF, u_, mesh)
        F_int = fn.assemble(f)
        Fch_int = fn.assemble(fch)
        Fc_int = F_int / pars['psi_sph'] ** 2 * np.sqrt(2) / pars[
                                                                 'er_fluid'] ** 1.5 * np.sqrt(
            pars['n0'])  # convert to Carnie ND variables
    else:
        f = calc_F_int(pars, SDs, bnd_MF, sd_MF, u_)
        F_int = 2 * fn.assemble(f)  # n.b. multiply by two as solve for a half space
        Fc_int = F_int / pars['psi_sph'] ** 2 * np.sqrt(2) / pars[
                                                                 'er_fluid'] ** 1.5 * np.sqrt(
            pars['n0'])  # convert to Carnie ND variables
        Fch_int = float('nan')

    return F_int, Fch_int, Fc_int, u_, mesh


def calc_fsig_grid_pts(geom_ND, Ra, Rx_xLTa, Rx_xGTa, smoothing=1e-4, doPlots=True):
    ''' 
    calculate a non-uniform grid for support of the spline use to convert the 
    fsig into a continuous function
    
    Ra fraction of a to continue with fine spacing (Rx_xLTa)
    Rx_xLTa - point spacing as a fraction of the Debye length for x<a (sphere radius)
    Rx_xGTa - point spacing as a fraction of the Debye length for x>a (sphere radius)
    
    '''

    import numpy as np
    import fenics_helpers as fh
    import matplotlib.pyplot as plt

    x_ND0 = np.concatenate((np.arange(0, geom_ND['a'] * Ra, geom_ND['lD'] * Rx_xLTa), \
                            np.arange(geom_ND['a'] * Ra, geom_ND['w'],
                                      geom_ND['lD'] * Rx_xGTa),
                            np.asarray([geom_ND['w']])))
    Nx = len(x_ND0)
    fP = fh.spline_con_gradient(np.arange(Nx), x_ND0, [0], [geom_ND['lD'] * Rx_xLTa],
                                s=smoothing, k=3)
    x_NDs = fP(np.arange(Nx))

    # --- rescale back to correct interval
    x_ND = (x_NDs - x_NDs[0]) / (x_NDs[Nx - 1] - x_NDs[0]) * geom_ND['w']

    # --- plot
    if doPlots:
        f, axs = fh.cfigw(ncols=2)
        axs[0].plot(np.arange(Nx), x_ND0, 'ok', label='init')
        axs[0].plot(np.arange(Nx), x_ND, 'o-r', label='smoothed')
        axs[0].legend(loc=2)
        fh.cxlbl(axs[0], '\mathrm{point}', 'x/l_0')
        axs[0].set_title('point distribution')
        axs[1].plot(x_ND, np.zeros(Nx), 'or')
        fh.cxlbl(axs[1], 'x/l_0', '1')
        plt.show()
        plt.pause(0.001)

    return x_ND, Nx


def eval_x_thresh_x_vs_F(x, F, threshold=0.05, doPlots=False, dF_lim_init=0.5):
    '''
    Evaluate the minimium value of F using the polyfit routine as well as the 
    threshold -> 1-threshold interval values by first evaluating the cumulative 
    density distribution based on the assumption that P \propto exp(-F)
    
    '''

    import scipy.interpolate as si
    import scipy.optimize as so
    import fenics_helpers as fh

    import numpy as np

    # --- ensure sorted (else interp fails)
    sort_index = np.argsort(x)
    x = x[sort_index]
    F = F[sort_index]

    x_eval = np.linspace(np.min(x), np.max(x), 1e3)
    fit = polyfit_minima(x, F, dF_lim_init, dF_lim_init)

    # --- extend using polyfit
    Fout = np.zeros(np.shape(x_eval))
    Imin = (x_eval >= np.min(fit['xc'])) & (x_eval < np.max(fit['xc']))
    Fout[Imin] = fit['fP'](x_eval[Imin])
    Fout[~Imin] = si.interp1d(x, F, kind='cubic', bounds_error=False)(x_eval[~Imin])

    # --- calc CDF
    P_unnorm = np.exp(-Fout + np.min(Fout))
    P_norm = P_unnorm / np.sum(P_unnorm)
    CDF = np.cumsum(P_norm)
    CDFf = si.interp1d(x_eval, CDF, kind='cubic')

    # --- find threshold values
    x_neg = so.fsolve(lambda x: CDFf(x) - threshold, fit['x_min'])
    x_pos = so.fsolve(lambda x: CDFf(x) - (1 - threshold), fit['x_min'])

    if doPlots:
        f, a = fh.cfigw(ncols=3)
        a[0].plot(x, F, 'ok')
        a[0].plot(fit['xc'], fit['yc'])
        a[0].plot(x_eval, Fout, '-r')

        a[1].plot(x_eval, P_norm, '-r')
        a[2].plot(x_eval, CDFf(x_eval), '-r')
        a[2].plot([x_neg, x_neg], [0, 1], '--b')
        a[2].plot([x_pos, x_pos], [0, 1], '--b')

    # --- build results
    RES = {}
    RES['x_lo'] = x_neg
    RES['x_hi'] = x_pos
    RES['x_min'] = fit['x_min']
    RES['xc'] = x_eval
    RES['Fc'] = Fout
    RES['P'] = P_norm
    RES['CDF'] = CDF

    return RES


def desc_part_sim(Rl, dgap_nm, lD_nm):
    '''
    returns a string detailing the various 
    parameters used in the simulation potential sweep parameters must be 
    specified explicitly to avoid confusion
    '''

    # --- simulation type
    desc = dict()
    desc['sim_type'] = 'unknown'
    if 'dgapr_nm' in Rl['dim_vars'].keys():
        desc['sim_type'] = 'gap size'
    elif 'lDr_nm' in Rl['dim_vars'].keys():
        desc['sim_type'] = 'Debye length'

    # --- describe ratchet geometry
    SFR_geom = Rl['dim_vars']['SFR_geom_nm']
    geom = Rl['dim_vars']['geom_nm']
    desc['geom'] = 'Period: %.0f nm, Steep slope: %.0f nm, Amplitude: %.0f nm' \
                   % (geom['w'], geom['p'], geom['D']) + '\n' \
                                                         'Recess: %.0f nm, Unpat PPA-Coverslip: %.0fnm\nparticle diameter: %.0f nm' \
                                                         % (SFR_geom['Recess'], dgap_nm,
                                                            geom['a'] * 2)

    # --- Debye length
    desc['lD'] = '%.1f nm' % (lD_nm)

    # --- surface charges
    desc['charge'] = ''

    for key, val in Rl['dim_vars']['psi_mV'].iteritems():
        desc['charge'] = desc['charge'] + key + ': %.0f mV' % (val) + ', '
    desc['charge'] = desc['charge'][0:len(desc['charge']) - 2]

    # -- build long string
    fld_desc = {'sim_type': 'Sweep parameter', 'geom': 'Geometry', \
                'charge': 'Charge (isol surfaces)', 'lD': r'$\lambda_D$'}
    str = '';
    for key, val in fld_desc.iteritems():
        str = str + r'\textbf{' + val + '}: ' + desc[key] + '\n'

    str = str[0:len(str) - 1]

    return str


def grid_gap_sim(xr_nm, L_nm, a_nm, hr_nm, Fr, xg_nm, yg_nm, method='cubic',
                 do3Dplot=False,
                 do2Dplot=True, axs_2D=None, title=''):
    """
    take the results of a ratchet simulation and upsample and grid the data so that it
    may be more easily visualised and converted into an effective potential

    :param xr_nm: The particle positions in the ratchet at which the free energy is
    evaluated
    :param a_nm: the particle radius
    :param L_nm: The absolute height of the silica (flat) surface
    :param hr_nm: A list of h_nm with one entry per particle position.  These are the set
                sph-flat surface spacings chosen at each x position
    :param Fr: A matrix of free energies corresponding to the particle position and
    :param xg_nm: Array of grid points (position of the particle center) for the
    upsampled data.
    :param yg_nm: Array of grid points (position of the particle center) for the
    upsampled data.
    :param title: title to display above the 2D plot axis
    :return: Fg - the re--gridded free energy
    """

    import fenics_helpers as fh
    import numpy as np
    import scipy.interpolate as si
    import matplotlib.pyplot as plt

    # --- Variables
    if do2Dplot and axs_2D is None:
        f, axs_2D = fh.cfig()

    if len(xr_nm) != len(hr_nm):
        raise ValueError('Inconsistent size of xr_nm and hr_nm')

    [xx, yy] = np.meshgrid(xg_nm, yg_nm)

    # --- Local Functions
    def yp(h, L):
        return L - a_nm - h

    # --- Construct vectors of data
    pts = np.zeros((0, 2))
    dF_all = np.zeros((0))
    for q in range(len(xr_nm)):
        hr_x_nm = hr_nm[q]
        yp_nm = yp(hr_x_nm, L_nm)
        xp_nm = xr_nm[q] * np.ones(np.shape(yp_nm))

        Fr_x = Fr[:, q]
        ind = ~np.isnan(Fr_x)

        Npts = np.sum(np.asarray(ind, dtype=int))
        ptsL = np.zeros((Npts, 2))

        ptsL[:, 0] = xp_nm[ind]
        ptsL[:, 1] = yp_nm[ind]

        pts = np.concatenate((pts, ptsL))
        dF_all = np.concatenate((dF_all, Fr_x[ind]))

    Fg = si.griddata(pts, dF_all, (xx, yy), method=method)

    # --- plotting
    if do3Dplot:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax3D = fig.add_subplot(111, projection='3d')
        ax3D.plot(pts[:, 0], pts[:, 1], dF_all, '-ok')

    if do2Dplot:
        from matplotlib import cm

        him = axs_2D.imshow(Fg - np.nanmin(Fg), extent=[np.min(xr_nm), np.max(xr_nm),
                                                        np.min(yg_nm), np.max(yg_nm)],
                            cmap=cm.afmhot)
        cbar = axs_2D.figure.colorbar(him, ax=axs_2D, orientation='horizontal')
        CS = axs_2D.contour(xx, yy, Fg - np.nanmin(Fg), [1, 3, 9], colors='w')
        axs_2D.clabel(CS, fmt=r'%.0fkT', fontsize=12)
        fh.cxulbl(axs_2D, 'x_p', 'nm', 'y_p', 'nm')
        axs_2D.set_title('$ \Delta F / k_B T $' + '\n' + title)

    return Fg


def calc_1D_V_eff(xg_nm, yg_nm, Fg, doPlots=True, axs=None, line_color=[0, 0, 0],
                  line_label=''):
    """

    Calculate the effective one dimensional potential corresponding to a free energy
    distribution F(xp,yp) (must be gridded).


    :param xg_nm: the gridded particle coordinates corresponding to the energies in Fg
    :param yg_nm: the gridded particle coordinates corresponding to the energies in Fg
    :param Fg: The free energies for the different particle positions
    :param doPlots: plots of p(xp,yp), p(xp), Veff(xp) [shifted to zero]
    :param axs: axes for plotting (list of length 3)
    :return: returns the effective potential shifted so that the minmium is at zero (not
    so that P=1)
    """

    import numpy as np
    import fenics_helpers as fh
    from matplotlib import cm
    import matplotlib.pyplot as plt

    # --- Variables
    if doPlots and axs is None:
        f, axs = fh.cfigw(ncols=3)

    dyg_nm = yg_nm[1] - yg_nm[0]
    dxg_nm = xg_nm[1] - xg_nm[0]

    # --- calc the effective potential: shifted to 0
    Pg_unnorm = np.exp(-(Fg - np.nanmin(Fg)))
    Pg_unnorm[np.isnan(Pg_unnorm)] = 0
    Pg_norm = Pg_unnorm / (np.sum(Pg_unnorm) * dyg_nm * dxg_nm)
    Pg_x = np.sum(Pg_norm, axis=0) * dyg_nm
    Pg_x_int = np.cumsum(Pg_x) * dxg_nm
    Feff_x = -np.log(Pg_x)
    Veff = Feff_x - np.min(Feff_x)

    if doPlots:
        him = axs[0].imshow(Pg_norm, extent=[np.min(xg_nm), np.max(xg_nm),
                                             np.min(yg_nm), np.max(yg_nm)],
                            cmap=cm.afmhot)
        cbar = axs[0].figure.colorbar(him, ax=axs[0], orientation='horizontal')

        fh.cxulbl(axs[0], 'x_p', 'nm', 'y_p', 'nm')
        axs[0].set_title('$P \propto \exp(-F/k_BT)$ (nm$^-2$)')

        axs[1].plot(xg_nm, Pg_x, '-', color=line_color, label=line_label)
        axs[2].plot(xg_nm, Veff, '-', color=line_color, label=line_label)

        fh.cxulbl(axs[1], 'x', 'nm', 'P(x)', 'nm^{-1}')
        fh.cxlbl(axs[2], 'x \ (\mathrm{nm})', '-\ln P(x)')

        axs[1].legend(loc=0, frameon=False)

        plt.draw()
        for p in range(len(axs)):
            axs[p].figure.tight_layout(pad=10)

    return Veff


def anl_Veff_dgap_scaling(x, Veff, d, lD_nm, doPlots=True, axs=None):
    """

    calculate the effective debye length using a simple LSA (no y integration) theory
    for a set of effective potentials measured
    at different gap heights d.

    :param x: an array of x values corresponding to the effective potentials
    :param Veff: one row per gap height
    :param d: array of gap heights
    :param lD_nm: the actual Debye length used in the simulation
    :param doPlots:
    :param axs: list of two axes for the scaled and scale factor (peak value)
    :return: dictionary of results containing lD_eff_nm, d, Veff_pk
    """

    import lmfit
    import fenics_helpers as fh
    import numpy as np
    import matplotlib.pyplot as plt

    # --- Variables
    if doPlots and axs is None:
        f, axs = fh.cfigw(ncols=3)

    Veff_pk = np.max(Veff, axis=1)
    dgapc_nm = np.linspace(np.min(d), np.max(d))

    # --- do fit
    k_guess = 1. / lD_nm
    A_guess = np.min(Veff_pk) * np.exp(k_guess * np.max(d) * 0.5)

    p = lmfit.Parameters()
    p.add_many(('A', A_guess), ('k', k_guess))

    def Vc(v, d):
        return v['A'] * np.exp(-v['k'] * d * 0.5)

    def residual(p):
        return Vc(p.valuesdict(), d) - Veff_pk

    mi = lmfit.minimize(residual, p, method='leastsq')
    R = dict(lD_nm_fit=1. / mi.params.valuesdict()['k'],
             A_fit = mi.params.valuesdict()['A'],
             lD_nm_actual=lD_nm,
             Veff_pk=Veff_pk,
             d_gap_nm=d)

    # --- plotting
    if doPlots:

        # --- display rescaled curves behaviour
        for p in range(len(d)):
            axs[0].plot(x, Veff[p, :] / Veff_pk[p], '-', color=fh.chot(p, len(d)))
        fh.cxulbl(axs[0], 'x', 'nm', '{V}_{eff} / \hat{V}_{eff}', '')

        # --- display scaling
        axs[1].plot(d, Veff_pk / np.min(Veff_pk), '.k', label='expt')
        axs[1].plot(dgapc_nm, Vc(mi.params.valuesdict(), dgapc_nm) / np.min(Veff_pk),
                    '-k', label='fit')
        fh.cxulbl(axs[1], 'd', 'nm', '\hat{V}_{eff} / \hat{V}_{eff,0}', '')
        axs[1].set_title('$\lambda_D$: est from slope: {:.1f} nm, actual: {:.1f} nm'
                         .format(1. / mi.params.valuesdict()['k'],
                                 lD_nm))
        axs[0].figure.tight_layout(pad=2)
        plt.pause(0.1)

    return R


def anl_ratchet_sim_vary_gap(Rl, doPlots=True, Nhi=300, Nxi=300):
    """

    :param Rl: results file generated by anl_ratchet_vary_gap
    :param doPlots:
    :param Nhi:  number of points used when re-gridding: y direction (h is particle
    surface spacing)
    :param Nxi:  number of points used when re-gridding: x direction
    :return: a dictionary containing the results of the analysis including the effective
     potentials and gridded/upsampled free energy
    """

    import fenics_helpers as fh
    import numpy as np


    # --- key vars
    Fr = Rl['Results']['Fr']
    hr_nm = Rl['Results']['hr_nm']
    xr_nm = Rl['Results']['xr_nm']

    Nh, Nx, Ngap = np.shape(Fr)
    dgapr_nm = Rl['dim_vars']['dgapr_nm']
    geom_nm = Rl['dim_vars']['geom_nm']
    SFR_geom_nm = Rl['dim_vars']['SFR_geom_nm']

    def yp(h, L):
        return L - geom_nm['a'] - h


    # --- Calc Veff(x)
    xpi_nm = np.linspace(np.min(xr_nm), np.max(xr_nm), Nxi)
    Veff = np.zeros((Ngap, Nxi))
    f = []
    for r in range(Ngap):
        SFR_geom_nm['d'] = dgapr_nm[r]
        L_nm = dgapr_nm[r] + SFR_geom_nm['Recess'] + SFR_geom_nm['Amplitude'] * 0.5

        title_str = desc_part_sim(Rl, dgapr_nm[r], Rl['dim_vars']['lD_nm'])

        # --- griddata (xp,yp)
        yp_min = yp(np.max(hr_nm[r * Nx + 0:(r + 1) * Nx - 1]), L_nm)
        yp_max = yp(np.min(hr_nm[r * Nx + 0:(r + 1) * Nx - 1]), L_nm)

        ypi_nm = np.linspace(yp_min, yp_max, Nhi)

        Fg = grid_gap_sim(xr_nm, L_nm, geom_nm['a'], hr_nm[r * Nx:(r + 1) * Nx],
                             Fr[:, :, r],
                             xpi_nm, ypi_nm,
                             do2Dplot =doPlots,
                             title=desc_part_sim(Rl, dgapr_nm[r],
                                                    Rl['dim_vars']['lD_nm']))

        Veff[r, :] = calc_1D_V_eff(xpi_nm, ypi_nm, Fg,
                                      line_color=fh.chot(r, Ngap),
                                      line_label='{} nm'.format(dgapr_nm[r]),
                                      doPlots=doPlots)
    R = anl_Veff_dgap_scaling(xpi_nm, Veff, dgapr_nm, Rl['dim_vars']['lD_nm'],
                                 doPlots=doPlots)


    Results = dict(FreeEnergy=dict(x=xpi_nm, y=ypi_nm, F=Fg),
                   EffectivePotential=dict(x=xpi_nm, gap=dgapr_nm, Veff=Veff, fit=R))

    return Results

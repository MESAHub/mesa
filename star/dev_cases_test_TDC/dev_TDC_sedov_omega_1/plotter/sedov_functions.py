import numpy as np
from globalvars import comvars as gv


def sedov_funcs(v):
    
##given the similarity variable v, return funtions, lambda, f, g, h, and th derivative of lambda with  v dlamdv

##although the ordinary differential equations are analytic, the sedov expressions
##appear to become singular for various combinations of parameters and at the lower
## limits of the integration range.  all these singularaties are rtemovable and done by this routine

##frequent combintaion and their derivative with v
##kamm equations 29 - 32, x4 a bit different to save a divide
##x1 is book's F

    x1 = gv.a_val * v
    dx1dv = gv.a_val
    
    cbag = max(gv.eps, gv.c_val * v - 1.0E0)
    x2 = gv.b_val * cbag
    dx2dv = gv.b_val * gv.c_val
    
    ebag = 1.0E0 - gv.e_val * v
    x3 = gv.d_val * ebag
    dx3dv = -gv.d_val * gv.e_val
    
    x4 = gv.b_val * (1.0E0 - 0.5E0 * gv.xg2 * v)
    dx4dv = - gv.b_val* 0.5E0 * gv.xg2

##transition region between standard and vacuum cases
##kamm page 15 or equations 88-92
##lambds = l_fun is book's zeta
##f_fun is book's V, g_fun is book's D, h_fun is book's P

    if gv.lsingular == True:
        l_fun = gv.rwant/gv.r2
        dlamdv = 0.0E0
        f_fun = l_fun
        g_fun = l_fun**(gv.xgeom - 2.0E0)
        h_fun = l_fun**gv.xgeom
        
##for the vacuum case in the hole
    elif gv.lvacuum and gv.rwant < gv.rvv:
        l_fun = 0.0E0
        dlamdv = 0.0E0
        f_fun = 0.0E0
        g_fun = 0.0E0
        h_fun = 0.0E0
        
##omega = omega2 = (2*(gamma - 1)+xgeom)/gamma case, denom2 = 0
##book expressions 20-22
    elif gv.lomega2 == True:
        beta0 = 1.0E0 / (2.0E0 * gv.e_val)
        pp1 = gv.gamm1 * beta0 
        c6 = 0.5E0 * gv.gamp1
        c2 = c6 / gv.gamma
        y = 1.0E0 / (x1 - c2)
        z = (1.0E0 - x1)*y
        pp2 = gv.gamp1 * beta0 * z
        dpp2dv = -gv.gamp1 * beta0 * dx1dv * y * (1.0E0 + z)
        pp3 = (4.0E0 - gv.xgeom - 2.0E0*gv.gamma) * beta0
        pp4 = -gv.xgeom * gv.gamma * beta0
        
        l_fun = x1**(-gv.a0) * x2**(pp1) * np.exp(pp2)
        dlamdv = (-gv.a0*dx1dv/x1 +pp1*dx2dv/x2 + dpp2dv) * l_fun
        f_fun = x1 * l_fun
        g_fun = x1**(gv.a0*gv.omega) * x2**pp3 * x4**gv.a5 * np.exp(-2.0E0*pp2)
        h_fun = x1**(gv.a0*gv.xgeom) * x2**pp4 * x4**(1.0E0 +gv.a5)
        
##omega = omega3 = xgeom*(2-gamma) case, denom3 = 0
##book expressions 23-25
    elif gv.lomega3 == True:
        beta0 = 1.0E0/ (2.0E0 * gv.e_val)
        pp1 = gv.a3 + gv.omega * gv.a2
        pp2 = 1.0E0 - 4.0E0 * beta0
        c6 = 0.5E0 * gv.gamp1
        pp3 = -gv.xgeom * gv.gamma * gv.gamp1 * beta0 * (1.0E0 - x1)/(c6 - x1)
        pp4 = 2.0E0 * (gv.xgeom * gv.gamm1 - gv.gamma) * beta0
        
        l_fun = x1**(-gv.a0) * x2**(-gv.a2) * x4**(-gv.a1)
        dlamdv = -(gv.a0*dx1dv/x1 + gv.a2*dx2dv/x2 + gv.a1*dx4dv/x4) * l_fun
        f_fun = x1 * l_fun
        g_fun = x1**(gv.a0*gv.omega) * x2**pp1 * x4**pp2 * np.exp(pp3)
        h_fun = x1**(gv.a0*gv.xgeom) * x4**pp4 * np.exp(pp3)
        
##for the standard or vacuum case not in the hole
##kamm equations 38 - 41
    else:
        l_fun = x1**(-gv.a0) * x2**(-gv.a2) * x3**(-gv.a1)
        dlamdv = -(gv.a0*dx1dv/x1 + gv.a2*dx2dv/x2 + gv.a1*dx3dv/x3) * l_fun
        f_fun = x1 * l_fun
        g_fun = x1**(gv.a0*gv.omega) * x2**(gv.a3 + gv.a2 * gv.omega) * x3**(gv.a4 +gv.a1 * gv.omega) * x4**gv.a5
        h_fun = x1**(gv.a0*gv.xgeom) * x3**(gv.a4 + gv.a1 *(gv.omega - 2.0E0)) * x4**(1.0E0 + gv.a5)
        
    return l_fun, dlamdv, f_fun, g_fun, h_fun
        
        
        
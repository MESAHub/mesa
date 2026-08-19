from globalvars import comvars as gv
from sedov_functions import sedov_funcs


def efun01(v):   
##evaluates teh first energy integrand, kamm equations 67 and 10.
##the (c_val*v - 1) term might be singular at v=vmin in the standard case
##the (1- c_val/gamma *v) term might be singular at v=vmin in the vacuum case
##due care should be taken for these removable singularities by the integrator

    l_fun, dlamdv, f_fun, g_fun, h_fun = sedov_funcs(v)
    efun1 = dlamdv * l_fun**(gv.xgeom + 1.0E0) * gv.gpogm * g_fun * v**2
    
    return efun1
    
 
   
def efun02(v): 
##evaluates teh first energy integrand, kamm equations 68 and 11.
##the (c_val*v - 1) term might be singular at v=vmin in the standard case
##the (1- c_val/gamma *v) term might be singular at v=vmin in the vacuum case
##due care should be taken for these removable singularities by the integrator
    l_fun, dlamdv, f_fun, g_fun, h_fun = sedov_funcs(v)
    z = 8.0E0/((gv.xgeom + 2.0E0 - gv.omega)**2 * gv.gamp1)
    efun2 = dlamdv * l_fun**(gv.xgeom - 1.0E0) * h_fun * z
    
    return efun2
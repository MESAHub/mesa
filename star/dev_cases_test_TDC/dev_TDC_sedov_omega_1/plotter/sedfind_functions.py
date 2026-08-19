from globalvars import comvars as gv
from sedov_functions import sedov_funcs


def sed_v_find(v):
##given corresponding physical distances, find the similarity variable v
##kamm equation 38 as a root find

    l_fun, dlamdv, f_fun, g_fun, h_fun = sedov_funcs(v)
    v_find = gv.r2*l_fun - gv.rwant
    
    return v_find
    
    
    
def sed_r_find(r):
##given the similarity variably v, find the sorrespoding physical distance
##kamm equation 38 as a root find
    
    l_fun, dlamdv, f_fun, g_fun, h_fun = sedov_funcs(gv.vwant)
    r_find = gv.r2*l_fun - r
    
    return r_find
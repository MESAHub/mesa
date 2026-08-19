import numpy as np
from globalvars import comvars as gv
from energy_functions import efun01, efun02
from sedfind_functions import sed_v_find, sed_r_find
from sedov_functions import sedov_funcs
from scipy import integrate
from scipy import optimize

def sed_1d(time, nstep, xpos, eblast, 
           rho0, vel0, ener0, pres0, cs0, gv):
#==============================================================================          
#..this routine produces 1d solutions for a sedov blast wave propagating
#..through a density gradient rho = rho**(-omega)
#..in planar, cylindrical or spherical geometry
#..for the standard, singular and vaccum cases.
#
#..standard case: a nonzero solution extends from the shock to the origin,
#.. where the pressure is finite.
#..singular case: a nonzero solution extends from the shock to the origin,
#.. where the pressure vanishes.
#..vacuum case : a nonzero solution extends from the shock to a boundary point,
#.. where the density vanishes making the pressure meaningless.
#
#..input:
#..time = temporal point where solution is desired seconds
#..xpos(i) = spatial points where solution is desired cm
#..eblast = energy of blast erg
#..rho0 = ambient density g/cm**3 rho = rho0 * r**(-omega_in)
#..omegain = density power law exponent rho = rho0 * r**(-omega_in)
#..vel0 = ambient material speed cm/s
#..pres0 = ambient pressure erg/cm**3
#..cs0 = ambient sound speed cm/s
#..gam0 = gamma law equation of state
#..xgeom_in = geometry factor, 3=spherical, 2=cylindircal, 1=planar
#
#..for efficiency reasons (doing the energy integrals only once),
#..this routine returns the solution for an array of spatial points
#..at the desired time point.
#
#..output:
#..den(i) = density 
#..enertot(i) = specific internal energy 
#..pres(i) = presssure 
#..vel(i) = velocity 
#..mach(i) = mach number 
#
#..although the ordinary differential equations are analytic,
#..the sedov expressions appear to become singular for various
#..combinations of parameters and at the lower limits of the integration
#..range. all these singularies are removable and done so by this routine.
#..these routines are written in real*16 precision because the
#..real*8 implementations simply run out of precision "near" the origin
#..in the standard case or the transition region in the vacuum case.
#
#..eps controls the integration accuracy, don’t get too greedy or the number
#..of function evaluations required kills.
#..eps2 controls the root find accuracy
#..osmall controls the size of transition regions
#==============================================================================
                
                            
                iprint = 1
                gv.eps = 1.0E-10
                eps2 = 1.0E-30
                osmall = 1.0E-4
                
#==============================================================================
##original code common block         
#==============================================================================
#                
#               global gamma, gamm1, gamp1, gpogm, xgeom, xg2, rwant
#                global r2, a0, a1, a2, a3, a4, a5, a_val, b_val, cval_
#                global d_val, e_val, omega, vv, xlam_want, vwant, rvv
#                global lsingular, lstandard, lvacuum, lomega2, lomega3
#
#create initial arrays of  density, velocity, pressure, energy, and sound speed.
#initally create arrays of zeros.                
                den = np.zeros(nstep)
                vel = np.zeros(nstep)
                pres = np.zeros(nstep)
                enertot = np.zeros(nstep)
                enertherm = np.zeros(nstep)
                enerkin = np.zeros(nstep)
                mach = np.zeros(nstep)
    
##return on uynphysical cases
##infinite mass
                if gv.omega >= gv.xgeom:
                    return
            
##transfer through comon block and create some frequent combonations.
#common block from original code is now emulated by global varibles.
                
                gv.gamm1 = gv.gamma - 1.0E0
                gv.gamp1 = gv.gamma + 1.0E0
                gv.gpogm = gv.gamp1 / gv.gamm1
                gv.xg2 = gv.xgeom + 2.0E0 - gv.omega
                denom2 = 2.0E0 * gv.gamm1 + gv.xgeom - gv.gamma * gv.omega
                denom3 = gv.xgeom * (2.0E0 - gv.gamma) - gv.omega
                
## post shock location v2 and location of singular point vstar
## kamm equations 18 and 19
                
                v2 = 4.0E0 / (gv.xg2 * gv.gamp1)
                vstar = 2.0E0 / (gv.gamm1*gv.xgeom + 2.0E0)
                
##set two logicals to determine type of solution
                gv.lstandard = False
                gv.lsingular = False
                gv.lavcuum = False
                
                if abs(v2 - vstar) <= osmall:
                    gv.lsingular = True
                elif v2 < vstar - osmall:
                    gv.lstandard = True
                elif v2 > vstar + osmall:
                    gv.lvacuum = True
                    
##two apparent singularities, books notation for omega2 and omega3
                gv.lomega2 = False
                gv.lomega3 = False
                
                if abs(denom2) <= osmall:
                    gv.lomega = True
                    denom2 = 1.0E-8
                    if iprint == 1:
                        print('omega2 case')
                elif abs(denom3) <= osmall:
                    gv.lomega3 = True
                    denom3 = 1.0E-8
                    if iprint == 1:
                        print('omega3 case')
                
##variuos exponents, kamm equations 42-47
                gv.a0 = 2.0E0/gv.xg2
                gv.a2 = -gv.gamm1/denom2
                gv.a1 = (gv.xg2*gv.gamma/(2.0E0 + gv.xgeom*gv.gamm1)) * (((2.0E0*(gv.xgeom*(2.0E0-gv.gamma)-gv.omega))/(gv.gamma*gv.xg2**2))-gv.a2)
                gv.a3 = (gv.xgeom - gv.omega) / denom2
                gv.a4 = gv.xg2 * ( gv.xgeom - gv.omega) * gv.a1/denom3
                gv.a5 = ((gv.omega * gv.gamp1) - (2.0E0 * gv.xgeom)) / denom3
                
##frequent combinations, kamm equations 33-37
                gv.a_val = 0.25E0 * gv.xg2 * gv.gamp1
                gv.b_val = gv.gpogm
                gv.c_val = 0.5E0 * gv.xg2 * gv.gamma
                gv.d_val = (gv.xg2 * gv.gamp1)/((gv.xg2 * gv.gamp1) - 2.0E0*(2.0E0 + (gv.xgeom*gv.gamm1)))
                gv.e_val = 0.5E0 * (2.0E0 + (gv.xgeom * gv.gamm1))
                
##evaluate energy intergrals
##the singular case can be done by hand; save some cpu cycles
##kamm equations 80, 81, and 85
                if gv.lsingular == True:
                    eval2 = gv.gamp1/(gv.xgeom*((gv.gamm1*gv.xgeom)+2.0E0)**2)
                    eval1 = 2.0E0 / gv.gamm1 * eval2
                    alpha = gv.gpogm* 2**(gv.xgeom) / (gv.xgeom*((gv.gamm1*gv.xgeom)+2.0E0)**2)
                    if int(gv.xgeom) != 1:
                        alpha = np.pi * alpha
                        
## for standard or vacuum cases
## v0 = post-shock orgin v0 and vv = vacuum boundry vv
## set the radius(rvv) coresponding to vv to 0 for now
## kamm equations 18, and 20
                else:
                    v0 = 2.0E0 / (gv.xg2 * gv.gamma)
                    gv.vv = 2.0E0 / gv.xg2
                    gv.rvv = 0.0E0
                    if gv.lstandard == True:
                        vmin = v0
                    if gv.lvacuum == True:
                        vmin = gv.vv
                    
##the first energy intergral
##in standard case the term (c_val*v - 1) mighht be singular at v = vmin
##kamm equations 18 and 28
#qromo subroutine from original code is replaced by romberg routine from scipy.integrate
                    eval1 = integrate.romberg(efun01, vmin, v2, tol = gv.eps, divmax = gv.its)                   
##int the vacuum case the term (1-c_val/gamma*v) might be singular at v=vmin
##in the standard case the term ( c_val * v-1) might be singular at v=vmin
                    eval2 = integrate.romberg(efun02, vmin, v2, tol = gv.eps, divmax = gv.its)
##in the vacuum case the term (1-c_val/gamma*v) might be singular at v=vmin
## kamm equations 57n and 58 for alpha in a slightly different form
                if gv.xgeom == 1.0:
                    alpha = (0.5E0 * eval1) + (eval2 / gv.gamm1)
                else:
                    alpha = (gv.xgeom - 1.0E0) * np.pi * (eval1 + 2.0E0 * eval2/gv.gamm1)
                    
##write what we have for the energy intergrals
                if iprint == 1:
                    print('xgeom=', gv.xgeom)
                    print('eblast=', eblast)
                    print('omega=', gv.omega)
                    print('alpha=', alpha)
                    print('j1=', eval1)
                    print('j2=', eval2)
                    
##immediate post-shock values
##kamm page 14 or equations 14, 16, 5, 13
##r2 = shock position, us = shock spees, rho1 = pre-shock density
##u2 = post-shock material speed, rho2 = post-shock density
##p2 = post-shock pressure, e2 = post-shock specific internal energy
##adn cs2 = post-shock sound speed
                    
                gv.r2 = (eblast/(alpha*rho0))**(1.0E0/gv.xg2) * time**(2.0E0/gv.xg2)
                us = (2.0E0/gv.xg2) * gv.r2/time
                rho1 = rho0 * gv.r2**(-gv.omega)
                u2 = 2.0E0 * us / gv.gamp1
                rho2 = gv.gpogm * rho1
                p2 = 2.0E0 * rho1 * us**2 / gv.gamp1
                e2 = p2/(gv.gamm1*rho2)
                cs2 = np.sqrt(gv.gamma*p2/rho2)
                
##find the radius corresponding to vv
#sub routine zeroin from original code is replaced by the brenth function from scipy.integrate
                if gv.lvacuum == True:
                    gv.vwant = gv.vv
                    gv.rvv = optimize.brenth(sed_r_find, 0.0E0, gv.r2, xtol = gv.eps)
                    
#print various values
                if gv.lstandard == True and iprint == 1:
                    print('r2=', gv.r2)
                    print('rho2=', rho2)
                    print('u2=', u2)
                    print('e2=', e2)
                    print('p2=', p2)
                    print('cs2=', cs2)
                    
                if gv.lvacuum == True and iprint ==1:
                    print('rv=', gv.rvv)
                    print('r2=', gv.r2)
                    print('rho2=', rho2)
                    print('u2=', u2)
                    print('e2=', e2)
                    print('p2=', p2)
                    print('cs2=', cs2)
                    
##now start loop over spatial positions
                for i in range(0,nstep):
                    gv.rwant = xpos[i]
                    
##if we are upstream from the shock
# area not yet disturbed by shock front
                    if gv.rwant >= gv.r2:
                        den[i] = rho0 * gv.rwant**(-gv.omega)
                        vel[i] = vel0
                        pres[i] = pres0
                        enertot[i] = ener0
                        enertherm[i] = pres[i]/(gv.gamm1)
                        enerkin[i] =  0.5E0 * den[i] * (vel[i])**2
                        mach[i] = vel[i] / cs0
                        
#check to see if solutins are valid in region,  if not delete the point from the array
#this will prevent python from stoping the code, and allow for solutions to be produced in valid regions.
#if statemant checks if sed_v_find are the same sign at two points( 0.9*v0 and v2)
                    elif sed_v_find(0.9E0*v0)/sed_v_find(v2) > 0:
                        xpos[i] = 0.0
                        den[i] = 0.0
                        vel[i] = 0.0
                        pres[i] = 0.0
                        enertot[i] = 0.0
                        enertherm[i] = 0.0
                        enerkin[i] = 0.0
                        mach[i] = 0.0  
                        
##if we are between the orgin and the shock front
                    else:
                        if gv.lstandard == True:
                            vat = optimize.brenth(sed_v_find, 0.90E0*v0, v2, xtol=eps2)
                        elif gv.lvacuum == True:
                            vat = optimize.brenth(sed_v_find, v2, 1.2E0*gv.vv, xtol=eps2)
                            
##the pysical solution
                        l_fun, dlamdv, f_fun, g_fun, h_fun = sedov_funcs(vat)
                        den[i] = rho2 * g_fun
                        vel[i] = u2 * f_fun
                        pres[i] = p2 * h_fun
                        enertot[i] = 0.0E0
                        enertherm[i] = 0.0E0
                        enerkin[i] = 0.0E0
                        mach[i] = vel[i] / cs0
                        if den[i] != 0.0:
                            enertot[i] = pres[i]/(gv.gamm1) + 0.5E0 * den[i] * (vel[i])**2
                            enertherm[i] = pres[i]/(gv.gamm1)
                            enerkin[i] = 0.5E0 * den[i] * (vel[i])**2

#remove void results from the arrays 
#void variable is to correct for indexing as values are deleted from the array
                void = 0
                for i in range(0, nstep):
                    if xpos[i] == 0:
                        den = np.delete(den, i - void)
                        vel = np.delete(vel, i - void)
                        pres = np.delete(pres, i - void)
                        enertot = np.delete(enertot, i - void)
                        enertherm = np.delete(enertherm, i - void)
                        enerkin = np.delete(enerkin, i - void)
                        mach = np.delete(mach, i - void)
                        void = void + 1
                        
                        
                xpos = np.trim_zeros(xpos, 'f')
                                
                
## end loop over spatial positions
                return den, vel, pres, enertot, enertherm, enerkin, mach, xpos
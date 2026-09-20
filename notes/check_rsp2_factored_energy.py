"""Check the factored RSP2 energy row without compiling or running MESA."""
from itertools import product
from pathlib import Path
import json

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(1976)


def residual(x, old, weights, mode, u_flag, seed, alfat, alfam,
             theta_p, theta_l, centered_velocity, factored):
    wm, w, wp, pi0, pi1, rho, lam, vel0, vel1 = x
    w_start, rho_start, old_lt0, old_lt1, old_vel0, old_vel1 = old
    a, b = weights
    face_w = np.array([a*w+(1-a)*wm, b*wp+(1-b)*w])
    face_w_div = np.array([a+(1-a)*wm/w, b*wp/w+(1-b)])
    face_energy = np.array([a*w*w+(1-a)*wm*wm, b*wp*wp+(1-b)*w*w])
    if mode:
        source_div = .5*(pi0/np.sqrt(face_energy[0])+pi1/np.sqrt(face_energy[1]))
        source = .5*(w/np.sqrt(face_energy[0])*pi0+w/np.sqrt(face_energy[1])*pi1)
        radiation_div = 0*w
    else:
        source_div = pi0
        source = (w+seed)*source_div
        source_div = source_div*(1+seed/w)
        radiation_div = .3*w/(rho*rho*lam*lam)
    damping_div = .7*w*w/lam
    dv = 1/rho-1/rho_start
    pressure = theta_p*.2*rho*w*w+(1-theta_p)*.2*rho_start*w_start*w_start
    velocity_work = np.array([vel0,vel1])
    if centered_velocity:
        velocity_work = .5*(velocity_work+np.array([old_vel0,old_vel1]))
    # Spatial/velocity factors are common to each row form. Signed work is
    # allowed, as with unequal current and time-centered velocity strains.
    heating_coefficient = alfam*lam*rho*rho*np.array([vel0,vel1])*velocity_work
    if u_flag:
        heating = .5*np.dot(heating_coefficient,face_w)
        heating_div = .5*np.dot(heating_coefficient,face_w_div)
    else:
        heating_div = heating_coefficient[0]
        heating = w*heating_div
    jumps=np.array([wm*wm-w*w,w*w-wp*wp])
    lt = -alfat*lam*rho*rho*face_w*jumps
    lt_div = -alfat*lam*rho*rho*face_w_div*jumps
    transport=theta_l*(lt[0]-lt[1])+(1-theta_l)*(old_lt0-old_lt1)
    transport_div=theta_l*(lt_div[0]-lt_div[1])+(1-theta_l)*(old_lt0-old_lt1)/w
    dt=.017
    if factored:
        return (w+theta_p*.2*rho*w*dv
                +((1-theta_p)*.2*rho_start*w_start*w_start*dv-w_start*w_start)/w
                +dt*(transport_div-source_div+damping_div+radiation_div-heating_div))
    return (w*w-w_start*w_start+pressure*dv
            +dt*(transport-source+w*damping_div+w*radiation_div-heating))/w


def check_equivalence():
    count=0; worst=0.
    for mode,u_flag,alfat,alfam,tp,tl,centered in product(
            (False,True),(False,True),(0.,.2),(0.,.25),(0.,.5,1.),(.5,1.),(False,True)):
        for seed in ((0.,) if mode else (0.,.1)):
            for _ in range(3):
                x=np.array([*10**RNG.uniform(-3,1,3),*RNG.uniform(-1,1,2),
                            *RNG.uniform(.3,2,2),*RNG.uniform(-1,1,2)])
                old=np.array([*RNG.uniform(.1,1,2),*RNG.uniform(-.1,.1,4)])
                weights=RNG.uniform(.01,.99,2)
                args=(old,weights,mode,u_flag,seed,alfat,alfam,tp,tl,centered)
                raw=residual(x,*args,False); new=residual(x,*args,True)
                assert abs(raw-new)<2e-12*max(1,abs(raw),abs(new))
                for j in range(len(x)):
                    z=x.astype(complex); z[j]+=1e-28j
                    d0=residual(z,*args,False).imag/1e-28
                    d1=residual(z,*args,True).imag/1e-28
                    error=abs(d0-d1)/max(1,abs(d0),abs(d1))
                    assert error<2e-10,(args,j,d0,d1)
                    worst=max(worst,error)
                count+=1
    return dict(cases=count,partials_per_case=9,max_scaled_derivative_difference=worst)


def check_failure_state():
    # Model 1975, cell 75 from the authorized 2026-09-19 crash reproduction.
    # Weights here multiply the target cell on each of its two bounding faces.
    faces=[(.46577287606805656,217022.45835707098,1.029286790271002e-6,-119784128863.06917),
           (.5342271239319434,1527.8244300009853,1.0991249375533113e-6,-21938482.096357785)]
    w_start=1.528268638841817e-23
    dt=.00022114478779988552
    cs=3321939.57231143; energy=27746698574758.57
    eq_div=9.993875541146074e-7; damping=1.1181261214750849e-10
    def factored(w):
        source_div=sum(.5*buoyancy*pi/np.sqrt(a*w*w+(1-a)*neighbor*neighbor)
                       for a,neighbor,buoyancy,pi in faces)
        return (w-w_start*w_start/w-dt*(source_div-damping*w*w+eq_div))*cs/energy
    h=w_start*1e-20
    derivative=factored(w_start+1j*h).imag/h
    reference=2.3944755541716443e-7
    assert abs(derivative/reference-1)<2e-15
    return dict(previous_native_AD=-7.739600595055e-5,
                corrected_derivative=derivative,high_precision_reference=reference)


def check_source():
    text=(ROOT/'star/private/hydro_rsp2.f90').read_text()
    assert 'resid_ad*max(1d0,s% csound(k)/w_00)' not in text
    assert 'div_by_w(s% RSP2_source_seed*source_div_w_ad,w_00)' in text
    assert '(1d0 - P_theta)*Ptrb_start*dV_ad - get_etrb_start(s,k), w_00)' in text
    assert 'compute_Eq_face(s, k+1, ierr, k)' in text
    assert 'compute_Lt(s, k+1, ierr, k)' in text
    assert 'if (.not. present(k_div_w)) s% Lt(k) = Lt%val' in text
    assert 'w_face = alfa + div_by_w(beta*wrap_w_m1(s,k),wrap_w_00(s,k))' in text
    assert 'w_face = div_by_w(alfa*wrap_w_00(s,k),wrap_w_m1(s,k)) + beta' in text
    assert 'get_etrb_start(s,k) == 0d0 .and. (.not. s% u_flag' in text
    return 'local factors, neighbor terms, history and flux-cache guards present'


if __name__=='__main__':
    result=dict(equivalence=check_equivalence(),failure_state=check_failure_state(),
                source=check_source(),limits='Standalone algebra and source checks; no native corrected executable run.')
    path=ROOT/'output/review/rsp3_small_w_fix_20260919/factored_energy_checks.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

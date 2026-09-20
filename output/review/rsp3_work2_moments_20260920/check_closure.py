"""Standalone covariance algebra and backward Euler checks, not MESA tests."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OUT=Path(__file__).resolve().parent
CD=(8/3)*np.sqrt(2/3)
CPHI=4*np.sqrt(2/3)
CPI_OLD=6*np.sqrt(2/3)
CPI_NEW=(CD+CPHI)/2


def rhs(t,y,corrected):
    energy,Pi,Phi=y
    w=np.sqrt(max(energy,0.))
    return [Pi-CD*w*energy,
            -(2/3)*energy+(Phi/3 if corrected else Phi)
            -(CPI_NEW if corrected else CPI_OLD)*w*Pi,
            -2*Pi-CPHI*w*Phi]


def phi_zero(t,y):return y[2]
phi_zero.terminal=True
phi_zero.direction=-1


def be_matrix(buoyancy,entropy_gradient,energy_decay,variance_decay,extra_decay,dt):
    # Unknowns are <u_r'^2> = (2/3)e_t, Pi, Phi.
    return np.array([[1+dt*energy_decay,-(2/3)*dt*buoyancy,0],
                     [-dt*entropy_gradient,1+dt*((energy_decay+variance_decay)/2+extra_decay),-dt*buoyancy/3],
                     [0,-2*dt*entropy_gradient,1+dt*variance_decay]])


def relative_min_eigenvalue(y):
    C=np.array([[y[0],y[1]],[y[1],y[2]]])
    return float(np.linalg.eigvalsh(C)[0]/max(np.linalg.norm(C,2),1e-300))


report={}
for corrected in [False,True]:
    sol=solve_ivp(lambda t,y:rhs(t,y,corrected),(0,10),[1e-12,0,1],
                  rtol=1e-11,atol=1e-14,events=phi_zero,max_step=.002)
    det=(2/3)*sol.y[0]*sol.y[2]-sol.y[1]**2
    report['corrected_ode' if corrected else 'old_ode']={
        'end_time':float(sol.t[-1]),'end_state':sol.y[:,-1].tolist(),
        'minimum_covariance_determinant':float(det.min())}
assert report['old_ode']['end_time']<2
assert report['corrected_ode']['end_time']==10
assert report['corrected_ode']['minimum_covariance_determinant']>-1e-12

rng=np.random.default_rng(20260920)
minimum=1.;max_residual=0.
for trial in range(3000):
    U,Phi=10**rng.uniform(-2,2,2)
    Pi=rng.uniform(-1,1)*np.sqrt(U*Phi)
    initial=np.array([U,Pi,Phi])
    buoyancy=10**rng.uniform(-1,1)
    entropy_gradient=-10**rng.uniform(-1,1)
    energy_decay,variance_decay,extra_decay=10**rng.uniform(-4,2,3)
    dt=10**rng.uniform(-6,6)
    A=be_matrix(buoyancy,entropy_gradient,energy_decay,variance_decay,extra_decay,dt)
    forcing=initial+dt*np.array([10**rng.uniform(-8,0),0.,0.])
    new=np.linalg.solve(A,forcing)
    minimum=min(minimum,relative_min_eigenvalue(new))
    error=np.max(abs(A@new-forcing))/max(np.max(abs(forcing)),np.max(abs(A)@abs(new)))
    max_residual=max(max_residual,float(error))
assert minimum>-2e-12
report['stable_backward_euler']={'cases':3000,'minimum_relative_covariance_eigenvalue':minimum,
                                'maximum_scaled_equation_error':max_residual}

# With the nonlinear turnover rates, freeze the trial w, solve the three
# linear moment equations, then find w^2 = 3U/2. This proves existence in
# these sampled stable cases; it is not a proposal to replace MESA's solver.
minimum=1.;max_residual=0.
for trial in range(300):
    e0,Phi0=10**rng.uniform(-2,2,2)
    Pi0=rng.uniform(-1,1)*np.sqrt((2/3)*e0*Phi0)
    initial=np.array([(2/3)*e0,Pi0,Phi0])
    b=10**rng.uniform(-1,1);h=-10**rng.uniform(-1,1);dt=10**rng.uniform(-5,5)
    def solve(w):
        return np.linalg.solve(be_matrix(b,h,CD*w,CPHI*w,0,dt),initial)
    def f(w):return w*w-1.5*solve(w)[0]
    hi=max(1.,np.sqrt(e0))
    while f(hi)<0:hi*=2
    w=brentq(f,0,hi,xtol=1e-13,rtol=1e-13)
    new=solve(w);minimum=min(minimum,relative_min_eigenvalue(new))
    max_residual=max(max_residual,abs(f(w))/max(w*w,1e-100))
assert minimum>-2e-12
assert max_residual<1e-9
report['nonlinear_stable_backward_euler']={'cases':300,'minimum_relative_covariance_eigenvalue':minimum,
                                          'maximum_relative_energy_consistency_error':max_residual}

# Continuous realizability does not remove a physical unstable growth
# timestep restriction for fixed coefficients without nonlinear saturation.
A=be_matrix(3,1,0,0,0,1)
new=np.linalg.solve(A,np.array([1.,0.,1.]))
report['unsaturated_unstable_BE_counterexample']=new.tolist()
assert new[0]<0 and new[2]<0

# Energy transport alone can leave the covariance cone; common implicit
# positive diffusion preserves it and conserves all three mass integrals.
initial=np.array([[1.,1.,1.],[0.,0.,0.]])
mix=np.array([[2.,-1.],[-1.,2.]])
energy_only=initial.copy();energy_only[:,0]=np.linalg.solve(mix,initial[:,0])
common=np.linalg.solve(mix,initial)
report['energy_only_transport_determinant']=float(np.linalg.det([[energy_only[0,0],energy_only[0,1]],[energy_only[0,1],energy_only[0,2]]]))
report['common_transport_minimum_eigenvalue']=min(relative_min_eigenvalue(x) for x in common)
assert report['energy_only_transport_determinant']<0
assert report['common_transport_minimum_eigenvalue']>=-1e-14
assert np.allclose(common.sum(axis=0),initial.sum(axis=0))

# Unit-control, homogeneous convective equilibrium remains unchanged.
energy=3/16;Pi=CD*energy**1.5;Phi=2*Pi/(CPHI*np.sqrt(energy))
old=(2/3)*energy+Phi-CPI_OLD*np.sqrt(energy)*Pi
new=(2/3)*energy+Phi/3-CPI_NEW*np.sqrt(energy)*Pi
assert max(abs(old),abs(new))<1e-14
report['shared_convective_equilibrium']={'energy':energy,'Pi':Pi,'Phi':Phi}
(OUT/'closure_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))

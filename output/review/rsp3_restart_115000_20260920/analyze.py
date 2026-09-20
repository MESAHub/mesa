"""Analyze the two authorized restarts and verify the local closure derivation."""
from pathlib import Path
from collections import Counter
import json, re
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent

def table(path):
    lines = path.read_text().splitlines()
    names = lines[5].split()
    data = np.array([list(map(float, s.split())) for s in lines[6:] if len(s.split()) == len(names)])
    return dict(zip(names, data.T))

results = {}
fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True, constrained_layout=True)
for name, label, color in [('baseline', 'alfat = 0', '#1f77b4'), ('alfat01', 'alfat = 0.1', '#d55e00')]:
    h = table(ROOT/name/'LOGS/history.data')
    text = (ROOT/name/'output.txt').read_text()
    retries, retry_rows, iterations = [], Counter(), 0
    last = ''
    for line in text.splitlines():
        if 'coeff' in line and 'max resid' in line:
            last = line
            iterations += 1
        if line.strip().startswith('retry:'):
            retries.append(int(line.split()[-1]))
            match = re.search(r'max resid\s+(\w+)\s+(\d+)', last)
            retry_rows[' '.join(match.groups()) if match else 'unknown'] += 1
    steps = h['model_number']-115000
    dt = 10**h['log_dt_sec']
    i = int(np.argmin(dt))
    results[name] = dict(models=len(steps), first_model=int(h['model_number'][0]),
        final_model=int(h['model_number'][-1]), retries=len(retries), retry_rows=dict(retry_rows),
        minimum_dt_s=float(dt[i]), minimum_dt_model=int(h['model_number'][i]),
        final_dt_s=float(dt[-1]), elapsed_days=float(h['day'][-1]-h['day'][0]+dt[0]/86400),
        median_iterations=float(np.median(h['num_iters'])), mean_iterations=float(np.mean(h['num_iters'])),
        accepted_iterations=dict(Counter(map(int,h['num_iters']))), total_printed_iterations=iterations,
        zones=sorted(set(map(int,h['num_zones']))), retry_models=retries,
        log_rel_cumulative_energy_error=float(h['log_rel_cumulative_energy_error'][-1]))
    axs[0].plot(steps,dt,label=label,color=color,lw=1.2)
    axs[1].plot(steps,h['num_iters'],label=label,color=color,lw=.85,alpha=.8)
axs[0].set_ylabel('Accepted timestep (s)'); axs[0].legend()
axs[1].set_ylabel('Accepted solver iterations'); axs[1].set_xlabel('Steps after model 115000')
for ax in axs: ax.grid(alpha=.2)
fig.suptitle('RSP3 restart: fixed 149-zone mesh, unchanged solver tolerances')
fig.savefig(ROOT/'restart_summary.png',dpi=160)
plt.close(fig)

# Exact Phi residual from profiles, uniform composition, dynamical gradL.
# gradL = grada_face * resolved pressure coefficient / temperature coefficient.
# Thus this uses the mapped entropy gradient of the dPrad/dm row, not gradT-grada.
checks = []
for n,k in [(4,74),(550,76),(633,74),(650,74),(740,74),(1333,74)]:
    p=table(ROOT/'baseline'/f'LOGS/profile{n}.data')
    old=table(ROOT/'baseline'/f'LOGS/profile{n-1}.data')
    h=table(ROOT/'baseline'/'LOGS/history.data')
    j=k-1
    a=p['dm'][j-1]/(p['dm'][j-1]+p['dm'][j]); b=1-a
    face=lambda key:a*p[key][j]+b*p[key][j-1]
    pressure_gradient=4*np.pi*p['radius_cm'][j]**2*face('rho')/(.5*(p['dm'][j-1]+p['dm'][j]))*(p['pressure'][j]-p['pressure'][j-1])/face('pressure')
    entropy_gradient=face('cp')*face('grada')*pressure_gradient*p['Y_face'][j]/p['gradL'][j]
    ef=a*p['w'][j]**2+b*p['w'][j-1]**2
    decay=4*np.sqrt(2/3)*np.sqrt(ef)/p['mlt_mixing_length'][j]
    rhs=2*entropy_gradient*p['Pi'][j]-decay*p['Phi'][j]
    ao=old['dm'][j-1]/(old['dm'][j-1]+old['dm'][j]); bo=1-ao
    area_rho=4*np.pi*old['radius_cm'][j]**2*(ao*old['rho'][j]+bo*old['rho'][j-1])
    temp=ao*old['temperature'][j]+bo*old['temperature'][j-1]
    Lscale=max(1,abs(p['L_start'][j]),1e-3*np.max(abs(p['L_start'])))
    flux_ref=Lscale/(area_rho*temp); v_ref=(Lscale/area_rho)**(1/3)
    scale=max(old['Phi'][j],1.5*(flux_ref/v_ref)**2)
    dt=float(10**h['log_dt_sec'][n-1])
    residual=(p['Phi'][j]-old['Phi'][j]-dt*rhs)/scale
    checks.append(dict(model=115000+n,face=k,Pi=float(p['Pi'][j]),Phi=float(p['Phi'][j]),
        Phi_start=float(old['Phi'][j]),entropy_gradient=float(entropy_gradient),Phi_rhs=float(rhs),
        Phi_scale=float(scale),dt=dt,reconstructed_residual=float(residual),
        unconstrained_variance_at_current_flux=float((old['Phi'][j]+2*dt*entropy_gradient*p['Pi'][j])/(1+dt*decay))))
results['profile_residual_checks']=checks

# Local uniform model: buoyancy=1, -ds/dr=-1, Lambda=1, no mean strain,
# radiation, viscosity or transport. The initial covariance is admissible.
cd=(8/3)*np.sqrt(2/3); cphi=4*np.sqrt(2/3)
def rhs(corrected):
    def evaluate(t,state):
        energy,Pi,Phi=state
        w=np.sqrt(max(0,energy))
        buoyancy_factor=1/3 if corrected else 1
        cpi=.5*(cd+cphi) if corrected else 6*np.sqrt(2/3)
        return [Pi-cd*w*energy,-(2/3)*energy+buoyancy_factor*Phi-cpi*w*Pi,
                -2*Pi-cphi*w*Phi]
    return evaluate

def boundary(t,state):return state[2]
boundary.terminal=True;boundary.direction=-1
results['continuous_local_checks']=[]
for corrected in (False,True):
    for tol in (1e-9,1e-11):
        sol=solve_ivp(rhs(corrected),(0,10),(1e-12,0,1),method='DOP853',rtol=tol,
                      atol=tol*1e-3,events=boundary,max_step=.002)
        energy,Pi,Phi=sol.y
        cov=(2/3)*energy*Phi-Pi**2
        results['continuous_local_checks'].append(dict(corrected=corrected,rtol=tol,
            boundary_times=sol.t_events[0].tolist(),final_time=float(sol.t[-1]),
            final_state=sol.y[:,-1].tolist(),final_derivative=rhs(corrected)(sol.t[-1],sol.y[:,-1]),
            min_covariance_determinant=float(min(cov))))

# Stationary unstable limit with buoyancy=entropy driving=Lambda=1.
eq=[]
for corrected in (False,True):
    factor=1/3 if corrected else 1
    cpi=.5*(cd+cphi) if corrected else 6*np.sqrt(2/3)
    energy=((2/3)+factor*2*cd/cphi)/(cd*cpi)
    Pi=cd*energy**1.5;Phi=2*Pi/(cphi*np.sqrt(energy))
    eq.append([energy,Pi,Phi])
assert np.allclose(eq[0],eq[1],rtol=1e-14,atol=0)
results['local_equilibrium_old_new']=eq

# Analytic determinant derivative with the corrected closure and optional
# isotropic compression, radiative cooling, positive shear heating.
rng=np.random.default_rng(149)
errors=[]
for _ in range(1000):
    energy,Phi=10**rng.uniform(-2,2,2)
    Pi=rng.uniform(-1,1)*np.sqrt((2/3)*energy*Phi)
    entropy_gradient,buoyancy,divu=rng.normal(size=3)
    pressure_factor=rng.uniform(0,1);cooling=rng.uniform(0,2);heating=rng.uniform(0,2)
    alpha_pi=rng.uniform(1,2)
    decay_e=cd*np.sqrt(energy);decay_phi=cphi*np.sqrt(energy)+2*cooling
    compression=(2/3)*pressure_factor*divu
    extra_decay=(alpha_pi-1)*.5*(cd+cphi)*np.sqrt(energy)
    de=buoyancy*Pi-(decay_e+compression)*energy+heating
    dpi=(2/3)*energy*entropy_gradient+(buoyancy/3)*Phi-(.5*(decay_e+decay_phi+compression)+extra_decay)*Pi
    dphi=2*entropy_gradient*Pi-decay_phi*Phi
    lhs=(2/3)*(de*Phi+energy*dphi)-2*Pi*dpi
    determinant=(2/3)*energy*Phi-Pi**2
    expected=-(decay_e+decay_phi+compression)*determinant+(2/3)*heating*Phi+2*extra_decay*Pi**2
    errors.append(abs(lhs-expected)/max(1,abs(lhs),abs(expected)))
assert max(errors)<1e-12
results['determinant_identity_max_relative_error']=max(errors)
(ROOT/'analysis.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps({k:v for k,v in results.items() if k not in ('baseline','alfat01')},indent=2))

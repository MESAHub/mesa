"""Check face-volume conservation, LNA inertia and matrix storage without MESA."""
import json
from decimal import Decimal, getcontext
getcontext().prec = 36
from pathlib import Path
import numpy as np
from scipy.linalg import solve_banded

rng = np.random.default_rng(20920)
checks = {}

def error(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.max(abs(a-b))/max(1., np.max(abs(a)), np.max(abs(b))))

def project(f):
    return .5*(f[:-1]+f[1:])

def masses(dm):
    return np.r_[dm[0]/2, .5*(dm[:-1]+dm[1:]), dm[-1]/2]

def edges(dm):
    mass = Decimal(0)
    edge = [mass]
    for d in dm:
        d = Decimal.from_float(float(d))
        edge.append(mass+d/2)
        mass += d
    return np.array(edge+[mass],dtype=object)

def remap(dm, old, new_dm, first=1, last=None):
    a, b = edges(dm), edges(new_dm)
    b *= a[-1]/b[-1]
    result = np.zeros((len(b)-1,3), dtype=float)
    j = 0
    for k in range(len(result)):
        while j < len(old):
            width = min(b[k+1],a[j+1])-max(b[k],a[j])
            if width > 0: result[k] += float(width)*old[j]
            if a[j+1] >= b[k+1]: break
            j += 1
    if last is None: last = len(new_dm)-1
    for k in range(len(result)):
        if first <= k <= last: continue
        result[min(last,max(first,k))] += result[k]
        result[k] = 0.
    return np.asarray(result/masses(new_dm)[:,None],float)

for n in [5, 6, 37, 150, 350]:
    dm = 10.**rng.uniform(-9,0,n)
    dm /= dm.sum()
    mu = masses(dm)
    w = rng.uniform(.1,5.,n+1); w[[0,-1]] = 0.
    phi = rng.uniform(.1,3.,n+1); phi[[0,-1]] = 0.
    pi = rng.uniform(-.9,.9,n+1)*np.sqrt(2/3*w*w*phi)
    old = np.column_stack((w*w,pi,phi))
    same = remap(dm,old,dm)
    checks[f'remap_identity_{n}'] = error(same,old)
    for new_dm in [np.r_[dm[:2],dm[2]/3,2*dm[2]/3,dm[3:]],
                   np.r_[dm[:2],dm[2]+dm[3],dm[4:]],
                   np.full(n+7,1/(n+7))]:
        new = remap(dm,old,new_dm)
        checks[f'remap_integrals_{n}_{len(new_dm)}'] = error(mu@old,masses(new_dm)@new)
        determinant = (2/3)*new[:,0]*new[:,2]-new[:,1]**2
        assert min(determinant) >= -2e-12, determinant
        assert min(new[:,0]) >= 0 and min(new[:,2]) >= 0
    old_cell = rng.uniform(.1,5.,n)
    # Native RSP initialization deposits half of each cell energy on each face.
    face_energy = np.zeros(n+1)
    for k in range(n):
        for j in (k,k+1):
            face_energy[min(n-1,max(1,j))] += .5*dm[k]*old_cell[k]
    rsp_w = np.sqrt(face_energy/mu)
    checks[f'rsp_initial_energy_{n}'] = error(dm@old_cell,mu@(rsp_w*rsp_w))

# v_flag: stress in cells, force and energy on faces. Unequal masses and moving boundaries.
n = 31
dm = rng.uniform(.1,3.,n); mu = masses(dm)
r = np.linspace(3.,.5,n+1); v = rng.normal(size=n+1); v0 = rng.normal(size=n+1)
w = rng.uniform(.1,4.,n+1); w[[0,-1]] = 0.
vc = .5*(v+v0); mass_correction = rng.uniform(.8,1.2,n)
coefficient = rng.uniform(.1,2.,n)
strain = np.diff(-v/r); strain_mid = np.diff(-vc/r)
chi = coefficient*project(w)*strain
uq = 4*np.pi*(np.r_[0.,chi]-np.r_[chi,0.])/(r*masses(dm*mass_correction))
eq_div_w_cell = 4*np.pi*coefficient*strain*strain_mid/dm
eq_face = w*(np.r_[dm*eq_div_w_cell,0.]+np.r_[0.,dm*eq_div_w_cell])/(2*mu)
eq_cell = project(w)*eq_div_w_cell
checks['v_flag_global_viscous_work'] = error(dm@eq_cell,-np.dot(masses(dm*mass_correction)*vc,uq))
checks['v_flag_cell_face_heating_budget'] = error(dm@eq_cell,mu@eq_face)
# Direct cell stress heating does not depend on the second inner neighbor.
h = 1e-30; k = 7
wv = w.astype(complex); wv[k+2] += 1j*h
perturbed = project(wv)*eq_div_w_cell
checks['v_flag_cell_heating_no_p2'] = abs(perturbed[k].imag/h)

# u_flag: shared face viscosity; a moving inner wall supplies explicit boundary work.
u, u0 = rng.normal(size=(2,n)); um = .5*(u+u0); rmid = .5*(r[:-1]+r[1:])
u_wall, u_wall0 = .13,-.07
strain = np.r_[0.,np.diff(-u/rmid),u[-1]/rmid[-1]-u_wall/r[-1]]
strain_mid = np.r_[0.,np.diff(-um/rmid),um[-1]/rmid[-1]-.5*(u_wall+u_wall0)/r[-1]]
wall_w = w.copy(); wall_w[-1] = w[-2]
chi = rng.uniform(.1,2.,n+1)*wall_w*strain
force = 4*np.pi*(chi[:-1]-chi[1:])/rmid
power = 4*np.pi*chi*strain_mid
power[-2] += power[-1]; power[-1] = 0.
eq_face = power/mu
boundary_work = -4*np.pi*chi[-1]*.5*(u_wall+u_wall0)/r[-1]
checks['u_flag_global_viscous_work'] = error(dm@project(eq_face)+um@force,boundary_work)

# Linearize the actual time-centered pressure increment at a static background.
rho = rng.uniform(.2,2.,n); alfap = .73
dlnrho = rng.normal(size=n); dw = rng.normal(size=n+1); dw[[0,-1]] = 0.
for theta in [0.,.37,1.]:
    h = 1e-30
    rr = rho*np.exp(1j*h*dlnrho); ww = w+1j*h*dw
    dV = 1/rr-1/rho
    outer = alfap/3*(theta*rr*ww[:-1]**2+(1-theta)*rho*w[:-1]**2)*dV
    inner = alfap/3*(theta*rr*ww[1:]**2+(1-theta)*rho*w[1:]**2)*dV
    face_increment = (np.r_[dm*outer,0.]+np.r_[0.,dm*inner])/mu
    derivative = np.imag(face_increment)/h
    expected = -alfap/3*w*w*(np.r_[dm*dlnrho,0.]+np.r_[0.,dm*dlnrho])/mu
    checks[f'pressure_LNA_inertia_{theta}'] = error(derivative,expected)
    native = -2*alfap/3*project(w*w)*dlnrho
    checks[f'pressure_LNA_zero_integral_{theta}'] = error(dm@(project(expected)-native),0.)

# Riemann pressure and contact velocity retain the inner-cell face-w derivative.
n = 15; k = 6; h = 1e-30
dm = rng.uniform(.5,2.,n); area = np.linspace(9.,1.,n+1)
rho = rng.uniform(.8,1.2,n); peos = rng.uniform(8.,12.,n)
u = rng.uniform(-.1,.1,n); gamma = 5/3; alfap = .4
w = rng.uniform(.5,1.5,n+1); w[[0,-1]] = 0.
def pressure_and_velocity(w):
    pc = peos + 2/3*alfap*rho*project(w*w)
    pf = np.zeros(n+1,dtype=w.dtype); vf = pf.copy()
    for j in range(1,n):
        ul,ur = u[j],u[j-1]; pl,pr = pc[j],pc[j-1]
        rl,rr = rho[j],rho[j-1]
        cl,cr = np.sqrt(gamma*pl/rl),np.sqrt(gamma*pr/rr)
        left = [ul-cl,ur-cr]; right = [ur+cr,ul+cl]
        sl = left[int(np.argmin(np.real(left)))]; sr = right[int(np.argmax(np.real(right)))]
        vf[j] = (ur*rr*(sr-ur)+ul*rl*(ul-sl)+pl-pr)/(rr*(sr-ur)+rl*(ul-sl))
        pf[j] = .5*(rl*(ul-sl)*(ul-vf[j])+pl+rr*(ur-sr)*(ur-vf[j])+pr)
    return pc,pf,vf
wc = w.astype(complex); wc[k+2] += 1j*h
pc,pf,vf = pressure_and_velocity(wc)
force = (area[k+1]*pf[k+1]-area[k]*pf[k]+pc[k]*(area[k]-area[k+1]))/dm[k]
checks['u_pressure_momentum_p2'] = error(force.imag/h,area[k+1]*pf[k+1].imag/(h*dm[k]))
for theta in [0.,.5,1.]:
    p0,pf0,vf0 = pressure_and_velocity(w)
    pv = theta*pf+(1-theta)*pf0; vv = theta*vf+(1-theta)*vf0
    local = pc[k]*(area[k]*vv[k]-area[k+1]*vv[k+1])/dm[k]
    flux = (area[k]*pv[k]*vv[k]-area[k+1]*pv[k+1]*vv[k+1])/dm[k]
    checks[f'u_local_work_p2_{theta}'] = error(local.imag/h,-p0[k]*area[k+1]*vv[k+1].imag/(h*dm[k]))
    checks[f'u_flux_work_p2_{theta}'] = error(flux.imag/h,-area[k+1]*(pv[k+1]*vv[k+1]).imag/(h*dm[k]))
assert abs(pf[k+1].imag/h) > 1e-3 and abs(vf[k+1].imag/h) > 1e-3

# Exact extra-band mapping for odd/even zone counts, including the line-search transpose.
for n in [1,2,3,8,9,17]:
    nv = 9; nh = 9; ie = 3; size = nv*n
    blocks = rng.normal(size=(3,n,nv,nv))
    blocks[1] += 100*np.eye(nv)
    extra = rng.normal(size=(n,nh,nh))
    dense = np.zeros((size,size)); groups = (n+1)//2; ng = 2*nv
    grouped = np.zeros((groups*ng,groups*ng))
    pl = np.zeros((groups,ng,ng)); pd = pl.copy(); pu = pl.copy()
    for k in range(n):
        r0 = k*nv
        g, p0 = k//2,k%2*nv
        pd[g,p0:p0+nv,p0:p0+nv] = blocks[1,k]
        dense[r0:r0+nv,r0:r0+nv] = blocks[1,k]
        if k:
            dense[r0:r0+nv,r0-nv:r0] = blocks[0,k]
            if p0 == 0: pl[g,:nv,nv:] = blocks[0,k]
            else: pd[g,nv:,:nv] = blocks[0,k]
        if k < n-1:
            dense[r0:r0+nv,r0+nv:r0+2*nv] = blocks[2,k]
            if p0 == 0: pd[g,:nv,nv:] = blocks[2,k]
            else: pu[g,nv:,:nv] = blocks[2,k]
        if k < n-2:
            dense[r0:r0+nh,r0+2*nv:r0+2*nv+nh] = extra[k]
            pu[g,p0:p0+nh,p0:p0+nh] = extra[k]
    if n%2: pd[-1,nv:,nv:] = np.eye(nv)
    for g in range(groups):
        r0 = g*ng
        grouped[r0:r0+ng,r0:r0+ng] = pd[g]
        if g: grouped[r0:r0+ng,r0-ng:r0] = pl[g]
        if g < groups-1: grouped[r0:r0+ng,r0+ng:r0+2*ng] = pu[g]
    checks[f'paired_matrix_{n}'] = error(grouped[:size,:size],dense)
    rhs = rng.normal(size=size)
    x = np.linalg.solve(grouped,np.r_[rhs,np.zeros(groups*ng-size)])[:size]
    checks[f'paired_solution_{n}'] = error(dense@x,rhs)
    for have_extra in [False,True]:
        matrix = dense.copy()
        if not have_extra:
            for k in range(n-2): matrix[k*nv:k*nv+nh,(k+2)*nv:(k+2)*nv+nh] = 0
        lower,upper=2*nv,(3 if have_extra else 2)*nv
        band = np.zeros((lower+upper+1,size))
        for i in range(size):
            for j in range(max(0,i-lower),min(size,i+upper+1)):
                band[upper+i-j,j] = matrix[i,j]
        xb = solve_banded((lower,upper),band,rhs)
        checks[f'banded_solution_{n}_{have_extra}'] = error(matrix@xb,rhs)
    truncated = dense.copy()
    for k in range(n-2): truncated[k*nv:k*nv+nh,(k+2)*nv:(k+2)*nv+nh] = 0
    transpose = truncated.T@rhs
    for k in range(n-2): transpose[(k+2)*nv:(k+2)*nv+nh] += extra[k].T@rhs[k*nv:k*nv+nh]
    checks[f'line_search_transpose_{n}'] = error(transpose,dense.T@rhs)

for name,value in checks.items():
    assert value < 2e-11,(name,value)
root = Path(__file__).resolve().parents[1]
out = root/'output/review/rsp3_face_w_plan_20260920'
out.mkdir(parents=True,exist_ok=True)
result = dict(checks=checks,scope='Algebra, overlap and matrix storage only. No MESA compilation or run.')
(out/'implementation_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(f'{len(checks)} checks passed; largest scaled error {max(checks.values()):.3e}')

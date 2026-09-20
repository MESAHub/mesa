"""Independent algebra and source checks. Does not compile or run MESA."""
import json
import re
from pathlib import Path

import numpy as np
from scipy.optimize import root
from check_rsp2_transport_startup import face_flux
from check_rsp2_shared_energy_boundary import solve_divided_transport

ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(20260919)


def read(path):
    return (ROOT / path).read_text()


def moments(x, a, driving, buoyancy, rate_pi, rate_phi, cooling, strain):
    w0, w1, pi_face, phi_face = x
    eface = a*w0*w0 + (1-a)*w1*w1
    wf = np.sqrt(eface)
    return np.array([
        (2/3)*eface*driving + buoyancy*phi_face
        - (rate_pi*wf + cooling + strain)*pi_face,
        2*driving*pi_face - (rate_phi*wf + 2*cooling)*phi_face,
        .5*w0/wf*buoyancy*pi_face,
        .5*w1/wf*buoyancy*pi_face,
    ])


def check_derivatives():
    worst = 0.
    for _ in range(240):
        a = RNG.uniform(.001, .999)
        x = np.array([10**RNG.uniform(-4, 2), 10**RNG.uniform(-4, 2),
                      RNG.uniform(-1, 1), RNG.uniform(.001, 1)])
        if _ % 3 == 0:
            x[0] = 0.  # One empty cell beside a finite energy face.
        g, b, c1, c2, rad, strain = RNG.uniform(-1, 1, 6)
        c1, c2, rad = abs(c1), abs(c2), abs(rad)
        args = (a, g, b, c1, c2, rad, strain)
        w0, w1, p, phi = x
        e = a*w0*w0+(1-a)*w1*w1
        wf = np.sqrt(e)
        dwf = np.array([a*w0, (1-a)*w1])/wf
        analytic = np.zeros((4, 4))
        analytic[0, :2] = (4/3)*np.array([a*w0, (1-a)*w1])*g - c1*p*dwf
        analytic[0, 2:] = [-(c1*wf+rad+strain), b]
        analytic[1, :2] = -c2*phi*dwf
        analytic[1, 2:] = [2*g, -(c2*wf+2*rad)]
        for cell in (0, 1):
            analytic[cell+2, :2] = -.5*b*p*x[cell]*dwf/wf**2
            analytic[cell+2, cell] += .5*b*p/wf
            analytic[cell+2, 2] = .5*b*x[cell]/wf
        numerical = np.empty_like(analytic)
        for j in range(4):
            z = x.astype(complex)
            z[j] += 1e-30j
            numerical[:, j] = moments(z, *args).imag/1e-30
        error = np.max(abs(analytic-numerical)/np.maximum(1, abs(analytic)))
        assert error < 2e-12, error
        worst = max(worst, error)
        if x[0] == 0:
            assert moments(x, *args)[2] == 0
            assert np.isfinite(analytic).all()
    return {"cases": 240, "max_scaled_error": worst}


def check_linear_moments():
    # The continuous LNA block and backward Euler block must share the rhs.
    # Check against exact constant-coefficient matrix exponent eigenvalues.
    for _ in range(100):
        e, g, b, rate_pi, rate_phi, cooling = RNG.uniform(.01, 1, 6)
        strain = RNG.uniform(-1, 1)
        a = np.array([[-rate_pi*np.sqrt(e)-cooling-strain, b],
                      [2*g, -rate_phi*np.sqrt(e)-2*cooling]])
        forcing = np.array([(2/3)*e*g, 0.])
        equilibrium = np.linalg.solve(a, -forcing)
        dt = 1e-4
        x = RNG.uniform(-1, 1, 2)
        step = np.linalg.solve(np.eye(2)-dt*a, x+dt*forcing)
        np.testing.assert_allclose(step-x-dt*(a@step+forcing), 0, atol=1e-14)
        modes = np.linalg.eigvals(a)
        step_modes = np.linalg.eigvals(np.linalg.inv(np.eye(2)-dt*a))
        np.testing.assert_allclose(np.sort_complex(step_modes),
                                   np.sort_complex(1/(1-dt*modes)), rtol=1e-12)
        np.testing.assert_allclose(a@equilibrium+forcing, 0, atol=1e-13)
    return 100


def check_collective_start():
    # Include the earlier quiet u block missed by a diagonal-only test.
    worst = 0.
    for weights in ([1/3, 2/3], [.01, .09, .9], [.4, .6]):
        n = len(weights)
        matrix = .125*np.tile(weights, (n, 1))-.1*np.eye(n)
        dt = .1
        guess = dt*np.maximum(0, matrix.sum(axis=1))
        assert np.all(guess > 0)
        for order in (range(n), range(n-1, -1, -1)):
            for k in order:
                linear = dt*matrix[k, k]
                available = dt*(matrix[k]@guess-matrix[k,k]*guess[k])
                disc = np.sqrt(linear*linear+4*available)
                value = (linear+disc)/2 if linear >= 0 else 2*available/(disc-linear)
                guess[k] = max(guess[k], value)
        sol = root(lambda w: w-dt*(matrix@w)/w, guess)
        assert sol.success and np.all(sol.x > 0)
        error = np.max(abs(sol.x**2-dt*matrix@sol.x))
        assert error < 1e-14
        worst = max(worst, error)
    return {"cases": 3, "max_original_energy_residual": worst}


def check_mesh_and_strain():
    for _ in range(500):
        dr0, dr1 = 10**RNG.uniform(-4, 2, 2)
        radii = np.array([3+dr0+dr1, 3+dr1, 3.])
        slope, offset = RNG.uniform(-5, 5, 2)
        velocity = slope*radii+offset
        weight = RNG.uniform(0, 1)
        interp = weight*(velocity[1]-velocity[2])/(radii[1]-radii[2])
        interp += (1-weight)*(velocity[0]-velocity[1])/(radii[0]-radii[1])
        np.testing.assert_allclose(interp, slope, rtol=1e-9, atol=1e-9)
        # AMR face insertion is point interpolation; cell energy is conserved.
        mass = 10**RNG.uniform(-4, 2)
        fraction = RNG.uniform(.01, .99)
        pi_outer, pi_inner = RNG.uniform(-4, 4, 2)
        phi_outer, phi_inner = RNG.uniform(0, 4, 2)
        assert min(pi_outer, pi_inner) <= (1-fraction)*pi_outer+fraction*pi_inner <= max(pi_outer, pi_inner)
        assert (1-fraction)*phi_outer+fraction*phi_inner >= 0
        e0, e1 = RNG.uniform(0, 4, 2)
        merged = fraction*e0+(1-fraction)*e1
        np.testing.assert_allclose(mass*merged, mass*fraction*e0+mass*(1-fraction)*e1)
    return 500


def check_transport_guess():
    count, worst = 0, 0
    for ratio in (1e-4, .01, 1., 100., 1e4):
        dm = np.array([ratio, 1.])
        for start in (np.array([0., 1.]), np.array([1., 0.])):
            for weighted in (False, True):
                for theta in (.5, 1.):
                    for step in (1e-6, .01, .1, .5):
                        dt = step*min(dm)
                        saved = start.copy()
                        guess = start.copy()
                        old, _ = face_flux(start, dm, weighted)
                        for order in ((0, 1), (1, 0)):
                            for k in order:
                                flux, df = face_flux(guess, dm, weighted)
                                sign = 1 if k == 0 else -1
                                rhs = dt*sign*(theta*flux+(1-theta)*old)/dm[k]
                                linear = dt*sign*theta*df[k]/dm[k]
                                available = start[k]**2+rhs-linear*guess[k]
                                if available < 0:
                                    continue
                                disc = np.sqrt(linear*linear+4*available)
                                sol = (linear+disc)/2 if linear >= 0 else 2*available/(disc-linear)
                                guess[k] = max(guess[k], sol)
                        assert np.array_equal(saved, start)
                        assert np.all(guess > 0)
                        w, iterations = solve_divided_transport(start, dm, dt, theta, weighted, guess)
                        assert abs(dm@(w*w-start*start))/sum(dm) < 1e-12
                        worst = max(worst, iterations)
                        count += 1
    return {"cases": count, "max_Newton_updates": worst}


def check_sources():
    hydro = read('star/private/hydro_rsp2.f90')
    solver = read('star/private/solver_support.f90')
    unpack = read('star/private/hydro_vars.f90')
    lna = read('star/private/star_LNA_support.f90')
    assert 'RSP2_w_fix_if_neg' not in solver+unpack
    assert 'w_00%val > 0d0 .and. w_00%val < s% csound(k)' in hydro
    assert 'max(1d0,s% csound(k)/w_00)' not in hydro
    assert 'call setup_dt_dLt_dm_ad(ierr, .true.)' in hydro
    assert 'Eq_cell%d1Array(i_w_m1)' in hydro and 'Eq_cell%d1Array(i_w_p1)' in hydro
    assert 'call rsp2_moment_rhs(s,k,Pi_rhs,Phi_rhs,ierr)' in hydro
    assert 'call rsp2_moment_rhs(s,k,Pi_rhs,Phi_rhs,ierr,use_time_centering=.false.)' in lna
    assert 'mtx%B(row_Pi,row_Pi) = 1d0' in lna
    assert 'mtx%B(row_Phi,row_Phi) = 1d0' in lna
    assert 'rsp2_zero_w_for_star_LNA(s,k,ierr)' in lna
    assert 'rsp2_entropy_flux' not in hydro+solver+unpack+lna
    defaults = read('star/defaults/controls_dev.defaults')
    for name, constant, factor in (('pi', 'x_ALFAPI', 6), ('phi', 'x_ALFAPHI', 4)):
        assert re.search(r'^\s*RSP2_alfa_'+name+r'\s*=\s*1d0\s*$', defaults, re.M)
        assert constant+' = '+str(factor)+'d0*sqrt_2_div_3' in hydro
        assert 's% RSP2_alfa_'+name+'*'+constant+'*w_face/Lambda_face' in hydro
    photo_in = read('star/private/photo_in.f90')
    photo_out = read('star/private/photo_out.f90')
    assert 'version /= 20' in photo_in and 'if (version >= 21)' in photo_in
    assert photo_in.index('s% RSP2_3equation_flag') < photo_in.index('call set_var_info') < photo_in.index('s% xh(:,1:nz)')
    assert 'write(iounit) s% RSP2_3equation_flag' in photo_out
    model_in = read('star/private/read_model.f90')
    model_out = read('star/private/write_model.f90')
    assert 'bit_for_RSP2_3equation = 17' in model_in
    assert 'xh(s% i_Pi,k) = vec(j)' in model_in and 'xh(s% i_Phi,k) = vec(j)' in model_in
    assert 'write1(s% Pi(k),ierr)' in model_out and 'write1(s% Phi(k),ierr)' in model_out
    # Check the actual two-column insertion point with optional trailing fields.
    for tail in ([], ['w_div_wc'], ['w_div_wc', 'j_rot']):
        names = ['lnd', 'lnT', 'lnR', 'L', 'u', 'w', 'Y']+tail
        old = {name: RNG.random(5) for name in names}
        pos = names.index('Y')+1
        new_names = names[:pos]+['Pi', 'Phi']+names[pos:]
        restored = [name for name in new_names if name not in ('Pi', 'Phi')]
        assert restored == names
        assert all(np.array_equal(old[name], old[restored[i]]) for i, name in enumerate(names))
    return 'layout and dispatch source checks passed; no binary restart executed'


def check_timestep_varcontrol():
    # Dimensional entropy moments must not change the structural timestep norm.
    source = read('star/private/timestep.f90')
    body = source.split('real(dp) function eval_varcontrol', 1)[1]
    skip = body.split('if (j ==', 1)[1].split('cycle', 1)[0]
    excluded = set(re.findall(r's%\s*i_(\w+)', skip))

    def measure(names, current, old, excluded):
        nz = current.shape[1]
        total = scales = 0.
        nterms = 0
        for j, name in enumerate(names):
            if name in excluded:
                continue
            row = sum(abs(sum(current[j, k-2:k+3])-sum(old[j, k-2:k+3]))/5
                      for k in range(2, nz-2))
            row += abs(2*current[j, 0]+current[j, 1]-2*old[j, 0]-old[j, 1])/3
            row += abs(2*current[j, -1]+current[j, -2]-2*old[j, -1]-old[j, -2])/3
            row += abs(sum(current[j, :3])-sum(old[j, :3]))/3
            # Preserve the existing endpoint convention in this comparison.
            if name == 'lnd':
                row /= 3
            total += row
            scales += max(1., abs(old[j, 0]))
            nterms += nz
        return total/scales/nterms

    count = 0
    worst_legacy_ratio = 0.
    for velocity in ('v', 'u'):
        for tail in ([], ['w_div_wc'], ['w_div_wc', 'j_rot']):
            names = ['lnd', 'lnT', 'lnR', 'lum', velocity, 'w', 'Y']+tail
            nz = 149
            old = np.zeros((len(names), nz))
            old[:3] = np.array([-8., 10., 24.])[:, None]
            current = old.copy()
            current[1, 20:60] += 1e-4
            baseline = measure(names, current, old, excluded)
            assert baseline > 0
            pos = names.index('Y')+1
            names3 = names[:pos]+['Pi', 'Phi']+names[pos:]
            for units in (1e-12, 1., 1e12):
                moments_old = np.zeros((2, nz))
                moments_old[:, 10:80] = np.array([1e15, 1e18])[:, None]*units
                moments_new = moments_old.copy()
                moments_new[:, 20:60] *= 1.+1e-8
                old3 = np.concatenate((old[:pos], moments_old, old[pos:]))
                current3 = np.concatenate((current[:pos], moments_new, current[pos:]))
                actual = measure(names3, current3, old3, excluded)
                assert actual == baseline, (velocity, tail, units, actual, baseline)
                moments_only = np.concatenate((old[:pos], moments_new, old[pos:]))
                assert measure(names3, moments_only, old3, excluded) == 0.
                assert measure(names3, old3, old3, excluded) == 0.
                legacy = measure(names3, current3, old3, excluded-{'Pi', 'Phi'})
                worst_legacy_ratio = max(worst_legacy_ratio, legacy/baseline)
                count += 1
    assert worst_legacy_ratio > 1e12
    return {'cases': count, 'structural_norm_unchanged': True,
            'max_legacy_inflation': worst_legacy_ratio}


def check_moment_initialization():
    flags = read('star/private/set_flags.f90')
    hydro = read('star/private/hydro_rsp2.f90')
    activate = flags.split('subroutine set_RSP2_flag(', 1)[1].split('end subroutine set_RSP2_flag', 1)[0]
    convert = flags.split('subroutine set_RSP2_3equation_flag(', 1)[1].split('end subroutine set_RSP2_3equation_flag', 1)[0]
    init = hydro.split('subroutine init_rsp2_moments(', 1)[1].split('end subroutine init_rsp2_moments', 1)[0]
    saved = activate.index('gradT_old = s% gradT')
    refreshed = activate.index('call set_vars(s, s% dt, ierr)')
    assigned = activate.index('= gradT_old - s% gradL')
    assert saved < activate.index('s% RSP2_flag = RSP2_flag') < refreshed < assigned
    assert activate.index('call set_vars(s, s% dt, ierr)', refreshed+1) > assigned
    assert convert.index('if (s% RSP2_3equation_flag .eqv. enabled) return') < convert.index('call set_vars')
    assert convert.index('= gradT_old - s% gradL') < convert.index('call unpack_xh(s,ierr)') < convert.index('call init_rsp2_moments')
    assert init.index('call get_rsp2_thermal_gradient') < init.index('entropy_gradient%val > 0d0 .and. Lc_old(k) > 0d0')
    assert 'RSP2_w_fix_if_neg' not in init

    # Imported RSP gradL may be unset. Both conversions must preserve gradT
    # using their own neutral reference; copying the old Y does not do this.
    count = 0
    for old_gradT in (.22, .39, .55):
        for neutral_one in (.15, .39, .6):
            for neutral_three in (.2, .4, .7):
                y_one = old_gradT-neutral_one
                y_three = (neutral_one+y_one)-neutral_three
                assert abs(neutral_three+y_three-old_gradT) < 2e-16
                count += 1

    def seed(energy, luminosity, area_rho_T, entropy_driving, forced=False):
        if forced or entropy_driving <= 0 or luminosity <= 0:
            return 0., 0.
        if energy == 0:
            raise ValueError('heat flux without turbulent energy')
        flux = luminosity/area_rho_T
        return flux, 1.5*flux*flux/energy

    cases = 0
    # Vanishing imported velocities no longer leave a finite variance in
    # stable layers. Active unstable faces retain their flux and bound.
    for weight in (.01, .5, .99):
        for speed in (0., 1e-20, 1e-8, 1., 1e6):
            energy = weight*speed**2+(1-weight)*(2*speed)**2
            old_luminosity = 9.*speed
            for driving in (-1., 0., 1.):
                Pi, Phi = seed(energy,old_luminosity,3.,driving)
                if driving <= 0 or speed == 0:
                    assert Pi == Phi == 0.
                else:
                    np.testing.assert_allclose(3*Pi,old_luminosity,rtol=2e-15)
                    np.testing.assert_allclose(Pi**2,(2/3)*energy*Phi,rtol=2e-15)
                cases += 1
            assert seed(energy,old_luminosity,3.,1.,forced=True) == (0.,0.)
            assert seed(energy,-old_luminosity,3.,1.) == (0.,0.)
    try:
        seed(0.,1.,3.,1.)
    except ValueError:
        pass
    else:
        raise AssertionError('inconsistent active zero energy face accepted')
    return {'gradient_conversions': count, 'moment_seeds': cases,
            'same_mode_restart_guard': True, 'initialization_uses_updated_Y': True}


def check_dependencies():
    modules = {}
    for folder in ('star/private', 'star/public', 'star/job', 'star_data'):
        for path in (ROOT/folder).rglob('*.f90'):
            text = path.read_text()
            match = re.search(r'^\s*module\s+(\w+)', text, re.M | re.I)
            if match:
                modules[match[1].lower()] = set(x.lower() for x in re.findall(r'^\s*use\s+(\w+)', text, re.M | re.I))
    visited, stack = set(), []
    def visit(name):
        assert name not in stack, stack+[name]
        if name in visited or name not in modules:
            return
        stack.append(name)
        for dep in modules[name]:
            visit(dep)
        stack.pop()
        visited.add(name)
    for name in modules:
        visit(name)
    return {"modules": len(modules), "cycles": 0}


if __name__ == '__main__':
    result = {
        'moment_and_source_derivatives': check_derivatives(),
        'continuous_and_discrete_moment_blocks': check_linear_moments(),
        'quiet_coupled_heating': check_collective_start(),
        'transport_with_updated_guess': check_transport_guess(),
        'AMR_and_affine_velocity_cases': check_mesh_and_strain(),
        'source_wiring': check_sources(),
        'timestep_varcontrol': check_timestep_varcontrol(),
        'moment_initialization': check_moment_initialization(),
        'module_dependencies': check_dependencies(),
        'limits': 'Standalone algebra and source checks only. No Fortran AD, MESA interpolation, binary restart, coupled hydro or eigenmode execution.'
    }
    path = ROOT/'output/review/rsp2_three_equation_20260919/implementation_checks.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))

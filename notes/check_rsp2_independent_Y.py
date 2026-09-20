"""Algebra and source checks only: does not compile or execute MESA."""
from pathlib import Path
import math
import random
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
RNG = random.Random(622075)


def close(actual, expected, scale=1.0, tol=2e-12):
    assert abs(actual - expected) <= tol * max(scale, abs(actual), abs(expected)), (actual, expected)


def check_flux():
    # Eliminate the independent flux row at fixed structure.
    for _ in range(1000):
        Lrad_coeff, Lc_gradT_coeff = [10**RNG.uniform(-3, 3) for _ in range(2)]
        gradL, Lt, L = RNG.uniform(0, 1), RNG.uniform(-20, 20), RNG.uniform(-20, 20)
        Y_face = (L - Lt - Lrad_coeff*gradL)/(Lrad_coeff + Lc_gradT_coeff)
        close(Lrad_coeff*(gradL + Y_face) + Lc_gradT_coeff*Y_face + Lt, L, scale=1000)
        assert Lrad_coeff + Lc_gradT_coeff > 0
    # The start-of-step luminosity scale is finite for zero and signed profiles.
    for luminosities in ([0, 0, 0], [1e30, 0, 1e20], [-1e30, -1e20, 0]):
        for luminosity in luminosities:
            scale = max(1, abs(luminosity), 1e-3*max(map(abs, luminosities)))
            assert math.isfinite(scale) and scale >= 1
            if max(map(abs, luminosities)):
                assert scale >= 1e-3*max(map(abs, luminosities))


def check_viscosity():
    for _ in range(500):
        nz = 9
        dm = [10**RNG.uniform(-2, 2) for _ in range(nz)]
        mass_correction = [RNG.uniform(0.6, 1.4) for _ in range(nz)]
        radius = list(reversed(sorted(RNG.uniform(1, 20) for _ in range(nz+1))))
        velocity = [RNG.uniform(-3, 3) for _ in radius]
        velocity_start = [RNG.uniform(-3, 3) for _ in radius]
        velocity_start[-1] = velocity[-1]  # fixed inner boundary
        velocity_work = [(a+b)/2 for a, b in zip(velocity, velocity_start)]
        # An arbitrary cell stress tests cancellation independently of its closure.
        Chi_cell = [RNG.uniform(-4, 4) for _ in dm]
        Chi_cell[0] = Chi_cell[1] = 0  # forced nonturbulent exterior
        Eq = [4*math.pi*Chi_cell[k]*(velocity_work[k]/radius[k] -
              velocity_work[k+1]/radius[k+1])/dm[k] for k in range(nz)]
        dual_mass = [0.5*dm[k]*mass_correction[k] +
                     (0.5*dm[k-1]*mass_correction[k-1] if k else 0) for k in range(nz)]
        Uq = [4*math.pi*((Chi_cell[k-1] if k else 0)-Chi_cell[k]) /
              (radius[k]*dual_mass[k]) for k in range(nz)]
        work = math.fsum(dual_mass[k]*velocity_work[k]*Uq[k] for k in range(nz))
        heating = math.fsum(dm[k]*Eq[k] for k in range(nz))
        boundary = -4*math.pi*Chi_cell[-1]*velocity[-1]/radius[-1]
        close(work + heating, boundary, scale=abs(work)+abs(heating))
        # Cell velocities use face stress, half-face heating, and the same lagged radius.
        cell_radius = [(radius[k]+radius[k+1])/2 for k in range(nz)]
        u = velocity_work[:-1]
        Chi_face = [0.0, 0.0] + [RNG.uniform(-4, 4) for _ in range(nz-1)]
        face_mass = [dm[0]/2] + [(dm[k-1]+dm[k])/2 for k in range(1, nz)] + [dm[-1]/2]
        shear = [0.0] + [u[k-1]/cell_radius[k-1] - u[k]/cell_radius[k]
                         for k in range(1, nz)] + [u[-1]/cell_radius[-1]-velocity[-1]/radius[-1]]
        Eq_face = [4*math.pi*c*d/m for c, d, m in zip(Chi_face, shear, face_mass)]
        Eq_cell = [(Eq_face[k]+Eq_face[k+1])/2 for k in range(nz)]
        force = [4*math.pi*(Chi_face[k]-Chi_face[k+1])/cell_radius[k] for k in range(nz)]
        work = math.fsum(u[k]*force[k] for k in range(nz))
        heating = math.fsum(dm[k]*Eq_cell[k] for k in range(nz))
        boundary = -4*math.pi*Chi_face[-1]*velocity[-1]/radius[-1]
        close(work + heating, boundary, scale=abs(work)+abs(heating))


def check_pii_staggering():
    # Check signed PII and its complete product rule on unequal face weights.
    def face_pii(values, alfa):
        Y, Lambda, Hp, Cp0, Cp1 = values
        return 0.5*math.sqrt(2/3)*Lambda/Hp*(alfa*Cp0+(1-alfa)*Cp1)*Y

    for alfa in (0.001, 0.25, 0.5, 0.999):
        Cp_face = alfa*2+(1-alfa)*4
        coefficient = 0.5*math.sqrt(2/3)*3/2*Cp_face
        for Y in (-1e8, -2, -1e-12, 0, 1e-12, 0.2, 1, 2, 1e8):
            values = [Y, 3.0, 2.0, 2.0, 4.0]
            pii = face_pii(values, alfa)
            close(pii, coefficient*Y)
            derivative = [coefficient, pii/3, -pii/2,
                          pii*alfa/Cp_face, pii*(1-alfa)/Cp_face]
            for j in range(len(values)):
                perturbed = [complex(value) for value in values]
                perturbed[j] += 1e-30j
                numerical = face_pii(perturbed, alfa).imag/1e-30
                close(numerical, derivative[j])

        # The flux uses interpolated cell w; the source uses the cell's own w.
        pii_outer = face_pii([2, 3, 2, 2, 4], alfa)
        pii_inner = face_pii([0.2, 3, 2, 2, 4], alfa)
        Hp_outer, Hp_inner = 2.0, 5.0
        area, T_rho_face = 7.0, alfa*3*2+(1-alfa)*11*7
        source_div_w = 0.5*(pii_outer/Hp_outer+pii_inner/Hp_inner)*11/(2*2)
        for w_cell, w_neighbor in ((0, 0), (0, 4), (1e-20, 4), (3, 4)):
            Lc = area*T_rho_face*pii_outer*(alfa*w_cell+(1-alfa)*w_neighbor)
            Source = w_cell*source_div_w
            assert math.isfinite(Lc) and math.isfinite(Source)
            if w_cell == 0:
                assert Source == 0
                if w_neighbor:
                    assert Lc > 0
            if w_cell != w_neighbor:
                wrong_w_face = math.sqrt(alfa*w_cell**2+(1-alfa)*w_neighbor**2)
                assert wrong_w_face > alfa*w_cell+(1-alfa)*w_neighbor


def check_remap():
    for _ in range(500):
        old_faces = [0.0] + sorted(RNG.random() for _ in range(19)) + [1.0]
        new_faces = [0.0] + sorted(RNG.random() for _ in range(11)) + [1.0]
        old_w = [10**RNG.uniform(-2, 2) for _ in range(len(old_faces)-1)]
        new_etrb = []
        for lo, hi in zip(new_faces, new_faces[1:]):
            integral = math.fsum(max(0, min(hi, b)-max(lo, a))*w*w
                                for a, b, w in zip(old_faces, old_faces[1:], old_w))
            new_etrb.append(integral/(hi-lo))
        old_total = math.fsum((b-a)*w*w for a, b, w in zip(old_faces, old_faces[1:], old_w))
        new_total = math.fsum((b-a)*e for a, b, e in zip(new_faces, new_faces[1:], new_etrb))
        close(new_total, old_total)
        # Split/merge reconstruction conserves dm*w^2, including its positivity fallback.
        dm, fraction, etrb, grad = 10**RNG.uniform(-3, 3), RNG.uniform(.01, .99), RNG.random()*10, RNG.uniform(-50, 50)
        dmR, dmL = dm*fraction, dm*(1-fraction)
        etrb_R = etrb + grad/4
        etrb_L = (dm*etrb-dmR*etrb_R)/dmL
        if min(etrb_R, etrb_L) < 0:
            etrb_R = etrb_L = etrb
        close(dmR*etrb_R+dmL*etrb_L, dm*etrb)
        close(math.sqrt((dmR*etrb_R+dmL*etrb_L)/dm)**2, etrb)

        # The common composition limiter must retain every species mass and sum X.
        abundance = [RNG.random() for _ in range(8)]
        abundance = [x/sum(abundance) for x in abundance]
        increment = [RNG.uniform(-2, 2) for _ in abundance]
        increment[-1] = -sum(increment[:-1])
        inner_increment = [-dmR*x/dmL for x in increment]
        factor = 1.0
        for x, dxR, dxL in zip(abundance, increment, inner_increment):
            for dx in (dxR, dxL):
                if dx > 0:
                    factor = min(factor, (1-x)/dx)
                elif dx < 0:
                    factor = min(factor, -x/dx)
        outer = [x+factor*dx for x, dx in zip(abundance, increment)]
        inner = [x+factor*dx for x, dx in zip(abundance, inner_increment)]
        close(sum(outer), 1)
        close(sum(inner), 1)
        for x, xR, xL in zip(abundance, outer, inner):
            assert -1e-14 <= min(xR, xL) and max(xR, xL) <= 1+1e-14
            close(dmR*xR+dmL*xL, dm*x)

        # Native v-grid half-face kinetic quadrature, compensated in thermal energy.
        vR, vL = RNG.uniform(-10, 10), RNG.uniform(-10, 10)
        v_new = vR + fraction*(vL-vR)
        old_ke = dm*(vR*vR+vL*vL)/4
        new_ke = (dmR*(vR*vR+v_new*v_new)+dmL*(v_new*v_new+vL*vL))/4
        thermal = 1000.0
        thermal_new = thermal + (old_ke-new_ke)/dm
        close(dm*thermal_new+new_ke, dm*thermal+old_ke)


def check_work_forms():
    # Test grid placement and power for every PdV/time-centering combination.
    # These are independent discrete identities, not execution of the Fortran.
    for grid in ('v', 'u'):
        for simple_pdv in (False, True):
            for time_centered in (False, True):
                for _ in range(100):
                    nz = 7
                    dm = [10**RNG.uniform(-1, 1) for _ in range(nz)]
                    radius = list(reversed(sorted(RNG.uniform(1, 20) for _ in range(nz+1))))
                    velocity = [RNG.uniform(-2, 2) for _ in range(nz+1)]
                    old_velocity = [RNG.uniform(-2, 2) for _ in range(nz+1)]
                    old_velocity[-1] = velocity[-1] = 0
                    work_velocity = [(a+b)/2 if time_centered or not simple_pdv else a
                                     for a, b in zip(velocity, old_velocity)]
                    rsp2_w = [RNG.uniform(0, 3) for _ in dm]
                    tdc_w = [RNG.uniform(0, 3) for _ in dm] + [0]
                    tdc_old_w = [RNG.uniform(0, 3) for _ in dm] + [0]
                    for convection in ('RSP2', 'TDC', 'TDC_explicit_momentum'):
                        if grid == 'v':
                            mass = [dm[0]/2] + [(dm[k-1]+dm[k])/2 for k in range(1, nz)]
                            strain = [velocity[k]/radius[k]-velocity[k+1]/radius[k+1]
                                      for k in range(nz)]
                            strain_work = [work_velocity[k]/radius[k]-work_velocity[k+1]/radius[k+1]
                                           for k in range(nz)]
                            w = (rsp2_w if convection == 'RSP2' else
                                 [(tdc_w[k]+tdc_w[k+1])/2 for k in range(nz)])
                            w_momentum = ([(tdc_old_w[k]+tdc_old_w[k+1])/2 for k in range(nz)]
                                          if convection == 'TDC_explicit_momentum' else w)
                            # Positive geometry/mixing-length coefficients multiply the strain.
                            chi_div_w = [RNG.random()*d for d in strain]
                            chi = [a*b for a, b in zip(chi_div_w, w_momentum)]
                            heating = math.fsum(4*math.pi*a*b*d for a, b, d in
                                                zip(chi_div_w, w, strain_work))
                            force = [4*math.pi*((chi[k-1] if k else 0)-chi[k])/radius[k]
                                     for k in range(nz)]
                            power = math.fsum(v*f for v, f in zip(work_velocity, force))
                        else:
                            cell_radius = [(radius[k]+radius[k+1])/2 for k in range(nz)]
                            face_mass = [dm[0]/2] + [(dm[k-1]+dm[k])/2 for k in range(1, nz)] + [dm[-1]/2]
                            strain = [0] + [velocity[k-1]/cell_radius[k-1]-velocity[k]/cell_radius[k]
                                            for k in range(1, nz)] + [velocity[nz-1]/cell_radius[-1]]
                            strain_work = [0] + [work_velocity[k-1]/cell_radius[k-1]-work_velocity[k]/cell_radius[k]
                                                 for k in range(1, nz)] + [work_velocity[nz-1]/cell_radius[-1]]
                            if convection == 'RSP2':
                                # Mass interpolation of cell w, not of w^2.
                                w = [rsp2_w[0]] + [(dm[k-1]*rsp2_w[k]+dm[k]*rsp2_w[k-1]) /
                                     (dm[k-1]+dm[k]) for k in range(1, nz)] + [rsp2_w[-1]]
                            else:
                                w = tdc_w[:-1] + [tdc_w[-2]]
                            w_momentum = (tdc_old_w[:-1] + [tdc_old_w[-2]]
                                          if convection == 'TDC_explicit_momentum' else w)
                            chi_div_w = [RNG.random()*d for d in strain]
                            chi = [a*b for a, b in zip(chi_div_w, w_momentum)]
                            eq_face = [4*math.pi*a*b*d/m for a, b, d, m in
                                       zip(chi_div_w, w, strain_work, face_mass)]
                            eq_cell = [(eq_face[k]+eq_face[k+1])/2 for k in range(nz)]
                            heating = math.fsum(m*e for m, e in zip(dm, eq_cell))
                            force = [4*math.pi*(chi[k]-chi[k+1])/cell_radius[k] for k in range(nz)]
                            power = math.fsum(v*f for v, f in zip(work_velocity, force))
                        # The TDC explicit option deliberately uses old w only in momentum.
                        # Retain its existing mismatch; do not assert exact cancellation there.
                        lag_difference = math.fsum(4*math.pi*a*(b-c)*d for a, b, c, d in
                                                  zip(chi_div_w, w, w_momentum, strain_work))
                        close(heating+power, lag_difference, scale=abs(heating)+abs(power))
                        if not simple_pdv:
                            # Conservative work always uses the exact finite-step KE velocity.
                            kinetic_power = math.fsum((v+v0)*f/2 for v, v0, f in
                                                      zip(velocity, old_velocity, force))
                            close(kinetic_power, power, scale=abs(power))


def check_tdc_preservation():
    # The shared cell length was deliberately changed in the preceding revision.
    # All other TDC routine bodies retain the base implementation.
    path = 'star/private/tdc_hydro.f90'
    original = subprocess.check_output(['git', 'show', '622075fbf:'+path], cwd=ROOT, text=True)
    current = (ROOT/path).read_text()
    cell_length = r'   function get_TDC_mixing_length_cell\(.*?end function get_TDC_mixing_length_cell'
    assert re.sub(cell_length, '', original.split('\ncontains\n', 1)[1], flags=re.S) == \
        re.sub(cell_length, '', current.split('\ncontains\n', 1)[1], flags=re.S)
    path = 'star/private/hydro_energy.f90'
    original = subprocess.check_output(['git', 'show', '622075fbf:'+path], cwd=ROOT, text=True)
    current = (ROOT/path).read_text()
    marker = "            else if (s% TDC_alpha_M >0d0 .and. s% MLT_option == 'TDC'"
    assert original.split(marker, 1)[1].split('            if (have_v_viscous_work)', 1)[0] == \
        current.split(marker, 1)[1].split('            if (have_v_viscous_work)', 1)[0]
    marker = '         subroutine setup_d_turbulent_energy_dt('
    assert original.split(marker, 1)[1].split('         end subroutine', 1)[0] == \
        current.split(marker, 1)[1].split('         end subroutine', 1)[0]


def check_energy():
    # For an ideal gas the static eps_grav entropy inertia is de + P d(1/rho).
    for _ in range(500):
        gas_constant, temperature, gamma = RNG.uniform(.1, 10), RNG.uniform(1, 100), RNG.uniform(1.1, 1.8)
        Cp = gamma*gas_constant/(gamma-1)
        Cv = gas_constant/(gamma-1)
        grada = (gamma-1)/gamma
        close(Cp*temperature*(1-grada), Cv*temperature)
        close(-Cp*temperature*grada, -gas_constant*temperature)


def check_sources():
    rsp2 = (ROOT/'star/private/hydro_rsp2.f90').read_text()
    riemann = (ROOT/'star/private/hydro_riemann.f90').read_text()
    assert 'wrap_Y_00' in rsp2 and 'rsp2_flux_residual' in rsp2
    assert 'get_TDC_mixing_length_face' in rsp2 and 'get_TDC_mixing_length_cell' in rsp2
    grad = rsp2.split('      function compute_RSP2_gradT', 1)[1].split('      end function', 1)[0]
    assert 'gradT = gradL + wrap_Y_00(s, k)' in grad
    assert not re.search(r'\biter\b|flux_resid|Lc_gradT_coeff', grad)
    assert 'u_face_ad(k) = s% u_face_ad(k) + Uq' not in riemann
    assert 'get_TDC_Lambda_face' not in rsp2
    assert 'compute_tdc_' not in rsp2
    assert 'Chi_cell = Chi_cell*wrap_w_00(s, k)' in rsp2
    assert 'w_face = alfa*wrap_w_00(s,k) + beta*wrap_w_m1(s,k)' in rsp2
    pii = rsp2.split('      function compute_PII_from_Y', 1)[1].split('      end function', 1)[0]
    assert 'PII_face = x_ALFAS*(Lambda_face/Hp_face)*Cp_face*Y_face' in pii
    assert not re.search(r'wrap_w_|wrap_e_|flux_limiter_function|PII_max', pii)
    assert 'use_TDC_enthalpy_flux_limiter' not in rsp2
    assert 'flux_limiter_function' not in rsp2
    source = rsp2.split('      function compute_Source_div_w', 1)[1].split('      end function', 1)[0]
    assert 'PII_face_00 = s% PII_ad(k)' in source
    assert 'PII_face_p1 = shift_p1(s% PII_ad(k+1))' in source
    assert '0.5d0*(PII_face_00/Hp_face_00 + PII_face_p1/Hp_face_p1)' in source
    assert 'Uq_cell = compute_Uq_dm_cell(s, k, ierr)' in riemann
    hydro_vars = (ROOT/'star/private/hydro_vars.f90').read_text()
    assert 's% Y_face(1:nz) = 0d0' not in hydro_vars
    flux = rsp2.split('      function rsp2_flux_residual', 1)[1].split('      end function', 1)[0]
    assert '1d-3*maxval(abs(s% L_start(1:s% nz)))' in flux
    assert 'abs(L_expected%val)' not in flux
    solver = (ROOT/'star/private/solver_support.f90').read_text()
    assert 'skip4 = s% i_Y' not in solver
    assert 's% x_scale(i,k) = max(xscale_min, abs(s% xh_start(i,k)))' in solver
    for Y_start in (-1e6, -1e-12, 0, 1e-12, 1e6):
        assert math.isfinite(1/max(1, abs(Y_start)))
    closures = (ROOT/'star/private/star_LNA_turbulence_closures.f90').read_text()
    lc_lna = closures.split('      subroutine rsp2_convective_luminosity_for_star_LNA', 1)[1].split(
        '      end subroutine', 1)[0]
    assert 'PII_face_ad = s% PII_ad(k)' in lc_lna
    assert 'source_ad = compute_Source(s, k, ierr)' in closures
    assert 'integer :: num_outer, num_inner' not in closures
    assert 'real(dp) :: alpha_M, alfa, beta, num_outer, num_inner' in closures
    # Preserve real TDC mask thresholds: integer truncation changes the inner edge.
    assert 9 > 10-1.2 and not 9 > 10-int(1.2)
    amr = (ROOT/'star/private/adjust_mesh_split_merge.f90').read_text()
    assert 'if (i < nz_old) v_L = s% v(ip)' in amr
    # The penultimate cell's inner face is the last stored face, not the center.
    face_velocities, center_velocity = [3, -4, 5], -7
    dm_outer, dm_inner = 0.4, 0.6
    kinetic_before = (face_velocities[-2]**2+face_velocities[-1]**2)/4
    split_velocity = face_velocities[-2]+dm_outer*(face_velocities[-1]-face_velocities[-2])
    kinetic_after = (dm_outer*(face_velocities[-2]**2+split_velocity**2)+
                     dm_inner*(split_velocity**2+face_velocities[-1]**2))/4
    close(kinetic_after+(kinetic_before-kinetic_after), kinetic_before)
    assert kinetic_before != (face_velocities[-2]**2+center_velocity**2)/4
    for directory in ('star/private', 'star_data/public', 'star_data/private'):
        for path in (ROOT/directory).iterdir():
            if path.suffix not in ('.f90', '.inc'):
                continue
            text = path.read_text()
            assert not re.search(r'\bi_Hp\b|\bi_Hp_(m1|00|p1)\b|i_equ_Hp|wrap_Hp_|RSP2_assume_HSE', text, re.I), path


if __name__ == '__main__':
    for check in (check_flux, check_pii_staggering, check_viscosity, check_work_forms, check_tdc_preservation,
                  check_remap, check_energy, check_sources):
        check()
        print(check.__name__ + ': PASS')
    print('This script does not compile or run MESA. Runtime convergence and MESA AD checks remain user-run validation.')

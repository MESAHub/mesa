"""Independent cell-length/stencil and viscous-power algebra; does not run MESA."""

from pathlib import Path
import re

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
N = 7
DM = np.array([0.6, 0.9, 0.5, 1.1, 0.7, 0.8, 0.4])
G = np.linspace(0.9, 1.1, N)
ALPHA = 1.8


def length(pressure, density, outer, inner, gout_mass, gin_mass, grav_const, alt, beta):
    """HSE/alternate cell height followed by the stable harmonic length."""
    gout = gout_mass / outer**2
    gin = gin_mass / inner**2 if inner.real > 0 else 0.0
    height = pressure / (density * 0.5 * (gout + gin))
    if alt and beta <= 0:
        alternate = np.sqrt(pressure / grav_const) / density
        if alternate.real < height.real:
            height = alternate
    base = ALPHA * height
    if beta <= 0:
        return base
    radius = (0.5 * (outer**3 + inner**3)) ** (1 / 3)
    radial = beta * radius
    if base.real <= radial.real:
        return base / (1 + base / radial)
    return radial / (1 + radial / base)


def check_length_derivatives():
    count = 0
    worst = 0.0
    # Includes the full centre, a fixed excised boundary, and interior cells.
    for inner in (0.0, 0.2, 1.0):
        for pressure in (0.01, 0.8, 100.0):
            for alt in (False, True):
                for beta in (0.0, 1e-8, 0.7, 1e8):
                    state = np.array([pressure, 1.3, 1.7, inner])
                    gout = 2.0 / state[2]**2
                    gin = 0.1 / inner**2 if inner > 0 else 0.0
                    height = pressure / (state[1] * 0.5 * (gout + gin))
                    dlogh = np.array([1., -1., 2*gout/(gout+gin), 2*gin/(gout+gin)])
                    alternate = np.sqrt(pressure) / state[1]
                    if alt and beta <= 0 and alternate < height:
                        height = alternate
                        dlogh = np.array([0.5, -1., 0., 0.])
                    expected = dlogh
                    if beta > 0:
                        base = ALPHA * height
                        radial = beta * (0.5*(state[2]**3 + inner**3))**(1/3)
                        dlogr = np.array([0., 0., state[2]**3, inner**3])
                        dlogr /= state[2]**3 + inner**3
                        expected = radial/(base+radial)*dlogh + base/(base+radial)*dlogr
                    value = length(*state, 2., 0.1, 1., alt, beta)
                    for j in range(4):
                        shifted = state.astype(complex)
                        shifted[j] *= np.exp(1e-30j)
                        actual = length(*shifted, 2., 0.1, 1., alt, beta).imag / 1e-30 / value
                        worst = max(worst, abs(actual - expected[j]))
                        np.testing.assert_allclose(actual, expected[j], atol=2e-14, rtol=2e-14)
                    count += 1
    return count, worst


def stresses(state, centre, alt, beta, tdc):
    radius, density, temperature = np.exp(state[:3])
    velocity, w = state[3:]
    mass = np.cumsum(DM[::-1])[::-1] + (0.3 if centre > 0 else 0.)
    chi = np.zeros(N, dtype=state.dtype)
    coefficient = np.zeros_like(chi)
    for k in range(N):
        rp1 = radius[k+1] if k+1 < N else centre
        vp1 = velocity[k+1] if k+1 < N else 0.
        gm1 = G[k+1]*mass[k+1] if k+1 < N else G[k]*0.3
        cell_length = length(density[k]*temperature[k], density[k], radius[k], rp1,
                             G[k]*mass[k], gm1, G[k], alt, beta)
        w_cell = w[k]
        if tdc:
            w_cell = 0.5*(w[k] + (w[k+1] if k+1 < N else 0.))
        coefficient[k] = 16*np.pi/3*0.25*density[k]**2 * 0.5*(radius[k]**6 + rp1**6) \
            * cell_length*w_cell/DM[k]
        chi[k] = coefficient[k]*(velocity[k]/radius[k] - vp1/(rp1 or 1.))
    return chi, coefficient


def check_stencil_and_lna():
    radius = np.linspace(2.0, 0.4, N)
    state = np.array([np.log(radius), np.linspace(-0.2, 0.3, N),
                      np.linspace(0.3, -0.2, N), np.sin(np.arange(N)),
                      np.linspace(0.1, 0.6, N)])
    count = 0
    for centre in (0., 0.2):
        for alt in (False, True):
            for beta in (0., 0.7):
                for tdc in (False, True):
                    for var in range(5):
                        for j in range(N):
                            shifted = state.astype(complex)
                            shifted[var, j] += 1e-30j
                            dchi = stresses(shifted, centre, alt, beta, tdc)[0].imag/1e-30
                            # Each stress has only k and k+1 state dependencies.
                            for k in range(N):
                                if j not in (k, k+1):
                                    assert dchi[k] == 0.
                            # Shift_m1 loses nothing from these stresses in momentum.
                            duq = np.r_[0., dchi[:-1]] - dchi
                            for k in range(N):
                                if abs(j-k) > 1:
                                    assert duq[k] == 0.
                    static = state.copy()
                    static[3] = 0.
                    coefficient = stresses(static, centre, alt, beta, tdc)[1]
                    perturbed = static.astype(complex)
                    # Perturb geometry, thermodynamics, w and velocity together.
                    perturbed += 1e-30j * state
                    actual = stresses(perturbed, centre, alt, beta, tdc)[0].imag/1e-30
                    expected = coefficient*np.diff(np.r_[state[3]/radius, 0.])*-1
                    np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=2e-14)
                    count += 1
    return count


def check_viscous_power():
    rng = np.random.default_rng(1926)
    worst = 0.
    for _ in range(100):
        # Arbitrary cell stresses include masked turbulent regions. The identity
        # is independent of which mixing-length closure produced the stress.
        chi = rng.normal(size=N)
        chi[rng.random(N) < 0.2] = 0.
        mass_factor = rng.uniform(0.8, 1.2, N)
        mass_face = 0.5*(DM*mass_factor + np.r_[0., (DM*mass_factor)[:-1]])
        radius_work = rng.uniform(0.2, 2., N)
        velocity_work = rng.normal(size=N)
        inner_ratio = rng.normal()  # excised-boundary work; zero at a full centre
        shear = -np.diff(np.r_[velocity_work/radius_work, inner_ratio])
        uq = 4*np.pi*(np.r_[0., chi[:-1]] - chi)/(radius_work*mass_face)
        heating = 4*np.pi*chi*shear/DM
        kinetic = np.sum(mass_face*velocity_work*uq)
        thermal = np.sum(DM*heating)
        boundary = -4*np.pi*chi[-1]*inner_ratio
        error = abs(kinetic + thermal - boundary) / max(abs(kinetic), abs(thermal), 1.)
        worst = max(worst, error)
        assert error < 5e-14
    return worst


def check_source_paths():
    source = (ROOT/'star/private/tdc_hydro.f90').read_text()
    original = (ROOT/'output/review/cell_mixing_length_20260919/tdc_hydro.before.f90').read_text()
    def routine(text, name):
        return re.search(r'   function '+name+r'\(.*?end function '+name, text, re.S).group()
    for name in ('get_TDC_Hp_face', 'get_TDC_mixing_length_face', 'compute_tdc_Uq_face',
                 'compute_tdc_Eq_cell', 'compute_tdc_Uq_dm_cell'):
        assert routine(source, name) == routine(original, name)
    lna = (ROOT/'star/private/star_LNA_turbulence_closures.f90').read_text()
    assert 'Lambda_cell_for_tdc_chi_for_star_LNA' not in lna
    for name in ('rsp2_chi_coefficient_for_star_LNA', 'tdc_chi_coefficient_for_star_LNA'):
        body = lna.split('function '+name+'(', 1)[1].split('end function '+name, 1)[0]
        assert 'get_TDC_mixing_length_cell(s, k, ierr)' in body
    assert 'damping_ad = compute_D(s, k, ierr)' in lna
    assert 'rad_damping_ad = compute_Dr(s, k, ierr)' in lna


if __name__ == '__main__':
    cases, error = check_length_derivatives()
    print(f'Length derivatives: {cases} cases; max absolute log-derivative error {error:.3e}')
    print(f'Stress stencil and static LNA: {check_stencil_and_lna()} cases passed')
    print(f'Viscous power: 100 cases; max normalized error {check_viscous_power():.3e}')
    check_source_paths()
    print('Shared LNA calls and unchanged face/conservation routines: passed')
    print('Independent algebra/source checks only; no MESA compilation or execution.')

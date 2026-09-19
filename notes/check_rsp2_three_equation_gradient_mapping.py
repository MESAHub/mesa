"""Standalone checks of the proposed RSP2 thermal-gradient mapping.

No MESA build or execution. Compare the mapped source against the differential
implied by each original temperature residual, using 80-digit arithmetic.
These checks do not execute Fortran or validate its automatic differentiation.
"""

import itertools
import json
from pathlib import Path

import mpmath as mp
import numpy as np


mp.mp.dps = 80


def coefficients(x, row, dynamic, options, math):
    """Dimensionless manufactured face states; all dimensional factors retained."""
    ti, to, pi, po, density, radius, opacity, hp = x
    mass_weight, face_weight, radiation_limit, opacity_floor, mass_correction = options
    area = 4 * math.pi * radius**2
    dm = (1 + mass_weight / (1 - mass_weight)) / 2
    tf = face_weight * ti + (1 - face_weight) * to
    pf = face_weight * pi + (1 - face_weight) * po
    # Distinct current and time-centered geometry/pressure, as allowed by QHSE.
    area_qhse = area * 1.03
    ppoint = (mass_weight * pi + (1 - mass_weight) * po) * 1.07
    gravity = 0.31 / radius**2
    conversion = area * density / dm
    pressure = conversion * (pi - po) / pf
    pressure_reference = pressure if dynamic else conversion * gravity * dm * mass_correction / area_qhse / pf
    gradad = math.mpf('0.4') if math is mp else 0.4
    cp = 2.7 + 0.1 * tf
    # Set c=a_rad=1: the constants cancel without changing the identities.
    lrad_coeff = area * (4 / 3) * tf**4 / (opacity * density * hp)
    kap_row = opacity * 1.2
    if np.real(complex(kap_row)) < opacity_floor:
        kap_row = opacity_floor + 0 * kap_row
    radiation_factor = 1 + 0 * ti
    if radiation_limit:
        t4_jump = ti**4 - to**4
        sign = 1 if np.real(complex(t4_jump)) >= 0 else -1
        flux_ratio = area * sign * t4_jump / dm / (kap_row * (ti**4 + to**4) / 2)
        radiation_factor = (6 + 3 * flux_ratio) / (6 + (3 + flux_ratio) * flux_ratio)
    if row == 'standard':
        tpoint = mass_weight * ti + (1 - mass_weight) * to
        coefficient = area * density * gravity / (area_qhse * ppoint) * tpoint / tf
        grad_from_row = (ti - to) / tpoint / (dm * gravity / (area_qhse * ppoint))
        direct_temperature = (ti - to) / tf
    elif row == 'radiation':
        prad = tf**4 / 3
        coefficient = density * kap_row * lrad_coeff / (4 * area * radiation_factor * prad)
        delta_prad = (ti**4 - to**4) / 3
        grad_from_row = delta_prad * area**2 * radiation_factor / (dm * kap_row * lrad_coeff)
        direct_temperature = delta_prad / (4 * prad)
    else:
        pressure = conversion * math.log1p((pi - po) / po)
        pressure_reference = pressure
        coefficient = pressure
        grad_from_row = math.log(ti / to) / math.log(pi / po) if pi != po else 0
        direct_temperature = math.log(ti / to)
    return cp, gradad, coefficient, pressure, pressure_reference, grad_from_row, conversion * direct_temperature


def mapped(x, y, composition, row, dynamic, options, math):
    cp, ad, coefficient, pressure, reference, _, _ = coefficients(x, row, dynamic, options, math)
    if row == 'log':
        return ad + composition, cp * coefficient * (y + composition)
    gradl = (ad + composition) * reference / coefficient
    source = cp * (coefficient * y + reference * composition)
    if not dynamic:
        source += cp * ad * (reference - pressure)
    return gradl, source


def direct_eliminated(x, y, composition, row, dynamic, options):
    cp, ad, coefficient, pressure, reference, _, _ = coefficients(x, row, dynamic, options, mp)
    gradl = ad + composition if row == 'log' else (ad + composition) * reference / coefficient
    return cp * (coefficient * (gradl + y) - ad * pressure)


def relative_error(value, reference):
    return float(abs(value - reference) / max(abs(reference), mp.mpf('1e-50')))


def main():
    cases = neutral_cases = small_y_cases = coordinate_cases = derivative_cases = 0
    max_row_error = max_small_y_error = max_coordinate_error = max_derivative_error = 0.0
    for row, dynamic, weight, face_weight, limiter, floor, correction, hp in itertools.product(
        ('standard', 'log', 'radiation'), (False, True), (0.1, 0.5, 0.9),
        (0.3, 0.5, 0.8), (False, True), (0.1, 2.0), (0.7, 1.0, 1.3), (0.2, 1.0, 5.0)
    ):
        options = (weight, face_weight, limiter, floor, correction)
        x = list(map(mp.mpf, ('1.011', '1', '1.03', '1', '0.8', '1.2', '0.7', str(hp))))
        cp, ad, coefficient, pressure, reference, grad_from_row, measured_temperature = coefficients(
            x, row, dynamic, options, mp)
        composition = mp.mpf('0.017') if cases % 2 else mp.mpf('-0.009')
        gradl, _ = mapped(x, 0, composition, row, dynamic, options, mp)
        _, actual = mapped(x, grad_from_row - gradl, composition, row, dynamic, options, mp)
        expected = cp * (measured_temperature - ad * pressure)
        max_row_error = max(max_row_error, relative_error(actual, expected))
        cases += 1
        # Changing the reference is a Y-coordinate change at fixed physical gradT.
        if row != 'log':
            other_gradl, _ = mapped(x, 0, composition, row, not dynamic, options, mp)
            _, other_source = mapped(x, grad_from_row - other_gradl, composition, row, not dynamic, options, mp)
            max_coordinate_error = max(max_coordinate_error, relative_error(actual, other_source))
            coordinate_cases += 1
        if dynamic or row == 'log':
            xf = list(map(float, x))
            _, zero_source = mapped(xf, 0., 0., row, dynamic, options, np)
            assert zero_source == 0
            neutral_cases += 1
            for y in (-1e-20, 1e-20, 1e-12):
                _, small_source = mapped(xf, y, 0., row, dynamic, options, np)
                max_small_y_error = max(max_small_y_error, relative_error(small_source, cp * coefficient * mp.mpf(y)))
                small_y_cases += 1
    # Complex-step differentiation includes varying face state, QHSE, radiative
    # coefficient, opacity floor, and radiation factor. Compare with mp.diff of
    # the uncancelled thermodynamic differential, not the rearranged source.
    for row, dynamic, limiter, floor in itertools.product(
        ('standard', 'log', 'radiation'), (False, True), (False, True), (0.1, 2.0)
    ):
        options = (0.3, 0.6, limiter, floor, 1.1)
        x = [1.011, 1., 1.03, 1., .8, 1.2, .7, .2]
        y, composition = .002, -.001
        for j in range(len(x) + 1):
            xx = np.array(x, dtype=complex)
            yy = complex(y)
            if j < len(x):
                xx[j] += 1e-30j
            else:
                yy += 1e-30j
            _, result = mapped(xx, yy, composition, row, dynamic, options, np)
            actual = result.imag / 1e-30
            mx = list(map(mp.mpf, x))
            if j < len(x):
                def reference(v):
                    args = mx.copy()
                    args[j] = v
                    return direct_eliminated(args, mp.mpf(y), mp.mpf(composition), row, dynamic, options)
                expected = mp.diff(reference, mx[j])
            else:
                expected = mp.diff(lambda v: direct_eliminated(mx, v, mp.mpf(composition), row, dynamic, options), mp.mpf(y))
            error = float(abs(actual - expected) / max(1, abs(expected)))
            max_derivative_error = max(max_derivative_error, error)
            derivative_cases += 1
    # Actual-log row remains defined at zero and reversed pressure gradients.
    for pressure in (0.97, 1., 1.03):
        for y in (-1e-20, 0., 1e-20):
            _, source = mapped([1.011, 1., pressure, 1., .8, 1.2, .7, 1.], y, 0.,
                               'log', True, (.3, .6, False, .1, 1.), np)
            assert np.isfinite(source)
            assert np.sign(source) == np.sign((pressure - 1) * y)
    # Match the discrete HSE force exactly at high precision. Unlike the
    # dynamical-on cancellation, floating-point force error remains physical
    # input to this branch; it cannot be removed by an entropy-row rewrite.
    max_static_neutral_error = 0.0
    max_static_tiny_y_error = 0.0
    for row, weight, face_weight in itertools.product(('standard', 'radiation'), (.1, .5, .9), (.3, .5, .8)):
        x = list(map(mp.mpf, ('1.011', '1', '1.03', '1', '.8', '1.2', '.7', '1')))
        options = (weight, face_weight, True, 2., 1.1)
        dm = (1 + weight / (1 - weight)) / 2
        area_qhse = 4 * mp.pi * x[5]**2 * 1.03
        # gravity in coefficients uses the exactly represented float .31.
        x[2] = x[3] + (mp.mpf(.31) / x[5]**2) * dm * 1.1 / area_qhse
        cp, ad, coefficient, pressure, reference, _, _ = coefficients(x, row, False, options, mp)
        _, source = mapped(x, mp.mpf('1e-20'), 0, row, False, options, mp)
        max_static_tiny_y_error = max(max_static_tiny_y_error, relative_error(source, cp * coefficient * mp.mpf('1e-20')))
        xf = list(map(float, x))
        _, source = mapped(xf, 0., 0., row, False, options, np)
        max_static_neutral_error = max(max_static_neutral_error, float(abs(source) / (cp * ad * pressure)))
    assert max_row_error < 1e-60
    assert max_coordinate_error < 1e-60
    assert max_small_y_error < 5e-14
    assert max_derivative_error < 1e-12
    assert max_static_tiny_y_error < 1e-55
    assert max_static_neutral_error < 1e-12
    result = dict(row_equivalence_cases=cases, max_row_relative_error=max_row_error,
                  coordinate_invariance_cases=coordinate_cases, max_coordinate_relative_error=max_coordinate_error,
                  exact_neutral_cases=neutral_cases, tiny_y_cases=small_y_cases,
                  max_tiny_y_relative_error=max_small_y_error, derivative_cases=derivative_cases,
                  max_derivative_scaled_error=max_derivative_error,
                  static_neutral_cases=18, max_static_neutral_force_scaled_error=max_static_neutral_error,
                  max_static_tiny_y_high_precision_error=max_static_tiny_y_error,
                  limitations='Standalone mathematics only; no Fortran AD, EOS, restart, mesh, LNA or MESA execution.')
    path = Path(__file__).resolve().parents[1] / 'output/review/rsp2_three_equation_20260919/gradient_mapping_checks.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

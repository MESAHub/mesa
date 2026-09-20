"""Standalone spatial and time-discretization checks, not a MESA run.

All variables are dimensionless. The heat subsystem freezes the moment
coefficients and examines ds/dt = -dPi/dx, dPi/dt = -q*ds/dx - gamma*Pi.
It does not represent the complete stellar Jacobian or a calibrated closure.
"""

import json
from pathlib import Path

import numpy as np


def diffusion_operator(mass, conductance):
    matrix = np.zeros((len(mass), len(mass)))
    for i, value in enumerate(conductance):
        matrix[i, i] -= value / mass[i]
        matrix[i, i + 1] += value / mass[i]
        matrix[i + 1, i] += value / mass[i + 1]
        matrix[i + 1, i + 1] -= value / mass[i + 1]
    return matrix


def periodic_heat_checks():
    n = 32
    identity = np.eye(n)
    # Face i is between cells i and i+1. Cell widths are one.
    gradient = np.roll(identity, 1, axis=1) - identity
    divergence = -gradient.T
    to_face = 0.5 * (identity + np.roll(identity, 1, axis=1))
    to_cell = to_face.T
    face_heat = divergence @ gradient
    cell_heat = divergence @ to_face @ to_cell @ gradient
    alternating = (-1.0) ** np.arange(n)
    q = 1.0
    gamma = 1.0
    face_dynamic = np.block([
        [np.zeros((n, n)), -divergence],
        [-q * gradient, -gamma * identity],
    ])
    cell_dynamic = np.block([
        [np.zeros((n, n)), -divergence @ to_face],
        [-q * to_cell @ gradient, -gamma * identity],
    ])
    wave = 2 * np.pi / n
    result = {
        "cells": n,
        "face_Pi_diffusion_nullity": n - int(np.linalg.matrix_rank(face_heat)),
        "averaged_cell_Pi_diffusion_nullity": n - int(np.linalg.matrix_rank(cell_heat)),
        "alternating_entropy_response_face_Pi": float(np.linalg.norm(face_heat @ alternating, np.inf)),
        "alternating_entropy_response_cell_Pi": float(np.linalg.norm(cell_heat @ alternating, np.inf)),
        "face_Pi_dynamic_zero_modes": int(np.sum(np.abs(np.linalg.eigvals(face_dynamic)) < 1e-10)),
        "cell_Pi_dynamic_zero_modes": int(np.sum(np.abs(np.linalg.eigvals(cell_dynamic)) < 1e-10)),
        "face_Pi_wave_number_ratio_longest_wave": float(2 * np.sin(wave / 2) / wave),
        "cell_Pi_wave_number_ratio_longest_wave": float(np.sin(wave) / wave),
        "with_radiative_diffusivity_0.001_alternating_decay_face_Pi": 4 * 1.001,
        "with_radiative_diffusivity_0.001_alternating_decay_cell_Pi": 4 * 0.001,
        "nonlocal_Pi_diffusion_cannot_drive_checkerboard_Pi_if_cell_entropy_gradient_is_zero":
            bool(np.allclose(to_cell @ gradient @ alternating, 0)),
    }
    assert result["face_Pi_diffusion_nullity"] == 1
    assert result["averaged_cell_Pi_diffusion_nullity"] == 2
    return result


def covariance_and_time_checks(rng):
    mass = np.exp(rng.uniform(-2, 2, 17))
    conductance = np.exp(rng.uniform(-2, 2, 16))
    operator = diffusion_operator(mass, conductance)
    weights = np.linalg.solve(np.eye(len(mass)) - 100 * operator, np.eye(len(mass)))
    vectors = rng.normal(size=(len(mass), 2))
    # Moment matrix entries are (2e/3, Pi; Pi, Phi).
    moments = np.einsum("ni,nj->nij", vectors, vectors)
    mixed = np.einsum("ij,jab->iab", weights, moments)
    old_total = np.einsum("i,iab->ab", mass, moments)
    new_total = np.einsum("i,iab->ab", mass, mixed)
    three_cell_operator = diffusion_operator(np.ones(3), np.ones(2))
    crank_nicolson = np.linalg.solve(
        np.eye(3) - 5 * three_cell_operator, np.eye(3) + 5 * three_cell_operator)
    backward_euler = np.linalg.solve(np.eye(3) - 10 * three_cell_operator, np.eye(3))
    result = {
        "common_diffusion_min_BE_weight": float(weights.min()),
        "common_diffusion_row_sum_error": float(np.max(np.abs(weights.sum(axis=1) - 1))),
        "common_diffusion_min_covariance_eigenvalue": float(np.linalg.eigvalsh(mixed).min()),
        "common_diffusion_relative_integral_error": float(np.max(np.abs(new_total-old_total))/np.max(np.abs(old_total))),
        "three_cell_CN_dt10_min_weight": float(crank_nicolson.min()),
        "three_cell_BE_dt10_min_weight": float(backward_euler.min()),
        "CN_middle_variance_from_initial_0_1_0": float(crank_nicolson[1, 1]),
    }
    assert weights.min() >= 0
    assert np.linalg.eigvalsh(mixed).min() >= -1e-12
    assert crank_nicolson[1, 1] < 0
    return result


def dual_energy_checks(rng):
    # Faces at cell boundaries. End faces carry half a cell mass.
    mass = np.exp(rng.uniform(-3, 3, 20))
    face_mass = 0.5 * np.r_[mass[0], mass[:-1] + mass[1:], mass[-1]]
    face_energy = np.exp(rng.uniform(-2, 2, len(mass) + 1))
    cell_energy = 0.5 * (face_energy[:-1] + face_energy[1:])
    # Dual transport fluxes are at cell centers, with closed outer boundaries.
    center_flux = rng.normal(size=len(mass))
    face_rate = (np.r_[0, center_flux] - np.r_[center_flux, 0]) / face_mass
    mapped_rate = mass * 0.5 * (face_rate[:-1] + face_rate[1:])
    primal_flux = np.zeros(len(mass) + 1)
    primal_flux[1:-1] = (
        mass[1:] * center_flux[:-1] + mass[:-1] * center_flux[1:]
    ) / (mass[:-1] + mass[1:])
    flux_rate = primal_flux[:-1] - primal_flux[1:]
    result = {
        "cell_vs_dual_total_energy_error": float(abs(mass @ cell_energy - face_mass @ face_energy)),
        "local_projected_transport_flux_error": float(np.max(np.abs(mapped_rate - flux_rate))),
        "global_closed_transport_energy_error": float(abs(mapped_rate.sum())),
    }
    assert max(result.values()) < 1e-11
    return result


def remap_checks(rng):
    old_edges = np.r_[0, np.sort(rng.uniform(size=23)), 1]
    new_edges = np.r_[0, np.sort(rng.uniform(size=38)), 1]
    old_mass, new_mass = np.diff(old_edges), np.diff(new_edges)
    overlap = np.maximum(0, np.minimum(new_edges[1:, None], old_edges[None, 1:])
                         - np.maximum(new_edges[:-1, None], old_edges[None, :-1]))
    weights = overlap / new_mass[:, None]
    vectors = rng.normal(size=(len(old_mass), 2))
    moments = np.einsum("ni,nj->nij", vectors, vectors)
    remapped = np.einsum("ij,jab->iab", weights, moments)
    old_total = np.einsum("i,iab->ab", old_mass, moments)
    new_total = np.einsum("i,iab->ab", new_mass, remapped)
    assert np.linalg.eigvalsh(remapped).min() >= -1e-12
    assert np.max(np.abs(new_total-old_total)) < 1e-12
    return {
        "common_overlap_min_covariance_eigenvalue": float(np.linalg.eigvalsh(remapped).min()),
        "common_overlap_integral_error": float(np.max(np.abs(new_total-old_total))),
        "applies_to_primal_or_dual_control_volumes": True,
    }


if __name__ == "__main__":
    rng = np.random.default_rng(20260920)
    result = {
        "scope": "Standalone constant-coefficient and finite-volume identities; no MESA execution.",
        "heat_coupling": periodic_heat_checks(),
        "common_transport_and_time": covariance_and_time_checks(rng),
        "face_moment_energy_projection": dual_energy_checks(rng),
        "remapping": remap_checks(rng),
    }
    output = Path(__file__).with_name("layout_checks.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))

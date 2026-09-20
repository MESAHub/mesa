"""Source and numerical remesh regressions; does not compile or execute MESA."""
from bisect import bisect_right
from pathlib import Path
import math
import random
import re

ROOT = Path(__file__).resolve().parents[1]
RNG = random.Random(20260919)


def close(a, b):
    assert abs(a-b) <= 5e-12*max(1, abs(a), abs(b)), (a, b)


def overlap_average(old_faces, new_faces, values):
    # Independent finite-volume reference, with no donor-search assumptions.
    return [math.fsum(max(0, min(hi, b)-max(lo, a))*v
                     for a, b, v in zip(old_faces, old_faces[1:], values))/(hi-lo)
            for lo, hi in zip(new_faces, new_faces[1:])]


def remap_squared(old_faces, new_faces, old_values, dual_coordinates=False):
    # Numerical transcription of adjust1_u/adjust1_etrb's interval traversal.
    # Exercise the branch logic against overlap_average, not against itself.
    old_dq = [b-a for a, b in zip(old_faces, old_faces[1:])]
    new_dq = [b-a for a, b in zip(new_faces, new_faces[1:])]
    comes_from = [min(len(old_dq)-1, bisect_right(old_faces, x)-1)
                  for x in new_faces[:-1]]
    unchanged = [lo == old_faces[i] and hi == old_faces[i+1]
                 for lo, hi, i in zip(new_faces, new_faces[1:], comes_from)]
    old_xq, new_xq = old_faces, new_faces
    if dual_coordinates:
        # The former bug: face dual boundaries coupled to physical cell masses.
        old_xq = [0] + [(a+b)/2 for a, b in zip(old_faces[:-2], old_faces[1:-1])] + [1]
        new_xq = [0] + [(a+b)/2 for a, b in zip(new_faces[:-2], new_faces[1:-1])] + [1]
    result = []
    for k, width in enumerate(new_dq):
        if unchanged[k] and (k == 0 or unchanged[k-1]):
            result.append(old_values[comes_from[k]])
            continue
        outer = new_xq[k]
        inner = outer + width if k+1 < len(new_dq) else 1.0
        start = (len(old_dq)-1 if outer >= old_xq[-2] else
                 0 if k == 0 else comes_from[k-1])
        mass, integral = 0.0, 0.0
        for j in range(start, len(old_dq)):
            lo, hi = old_xq[j:j+2]
            if hi <= outer:
                continue
            if lo >= outer and hi <= inner:
                overlap = old_dq[j]
                mass += overlap
                if mass > width:
                    overlap -= mass-width
                    mass = width
            elif lo <= outer and hi >= inner:
                overlap = width
                mass += overlap
            elif inner <= hi:
                overlap = width-mass
                mass = width
            else:
                overlap = min(max(0, hi-outer), width-mass)
                mass += overlap
            assert overlap >= 0
            integral += old_values[j]**2*overlap
            if mass >= width:
                break
        if not dual_coordinates:
            close(mass, width)
        value = math.sqrt(max(0, integral/width))
        result.append(math.copysign(value, old_values[comes_from[k]]))
    return result


def check_remap():
    # Unequal cells; an unchanged inner cell follows two newly split cells.
    old_faces, new_faces, w = [0, .2, .6, 1], [0, .1, .2, .6, 1], [1, 3, 2]
    expected = [1, 1, 3, 2]
    for a, b in zip(remap_squared(old_faces, new_faces, w), expected):
        close(a, b)
    wrong = remap_squared(old_faces, new_faces, w, dual_coordinates=True)
    old_total = math.fsum((b-a)*v*v for a, b, v in zip(old_faces, old_faces[1:], w))
    wrong_total = math.fsum((b-a)*v*v for a, b, v in zip(new_faces, new_faces[1:], wrong))
    assert abs(wrong_total/old_total-1) > .01
    print(f"Former coordinate bug: {100*(wrong_total/old_total-1):.6g}% turbulent-energy error")

    meshes = [(old_faces, old_faces), (old_faces, new_faces), (new_faces, old_faces),
              ([0, 1], [0, .2, .9, 1]), ([0, .2, .9, 1], [0, 1]),
              ([0, .1, .2, .8, 1], [0, .1, .5, .8, 1])]
    for _ in range(500):
        old = [0] + sorted(RNG.random() for _ in range(19)) + [1]
        # Mix unchanged cells, split cells, and multi-cell merges.
        new = [0] + sorted([x for x in old[1:-1] if RNG.random() < .5] +
                           [RNG.random() for _ in range(8)]) + [1]
        meshes.append((old, new))
    for old, new in meshes:
        for values in ([0.0]*(len(old)-1), [3.0]*(len(old)-1),
                       [RNG.choice([-1, 1])*10**RNG.uniform(-3, 3) for _ in old[:-1]]):
            remapped = remap_squared(old, new, values)
            reference = overlap_average(old, new, [v*v for v in values])
            for v, e in zip(remapped, reference):
                close(v*v, e)
            close(math.fsum((b-a)*v*v for a, b, v in zip(old, old[1:], values)),
                  math.fsum((b-a)*v*v for a, b, v in zip(new, new[1:], remapped)))


def check_energy():
    # Ordinary remesh: w^2 already equals the cell overlap average, so the
    # thermal correction must not add turbulent energy a second time.
    for grid in ('u', 'v'):
        for _ in range(200):
            old = [0] + sorted(RNG.random() for _ in range(7)) + [1]
            new = [0] + sorted(RNG.random() for _ in range(11)) + [1]
            energy = [RNG.uniform(100, 200) for _ in old[:-1]]
            w = [RNG.random()*2 for _ in old[:-1]]
            velocity = [RNG.uniform(-2, 2) for _ in old]
            new_velocity = [RNG.uniform(-2, 2) for _ in new]
            new_velocity[-1] = velocity[-1]  # fixed inner boundary
            if grid == 'u':
                ke = [v*v/2 for v in velocity[:-1]]
                new_ke = [v*v/2 for v in new_velocity[:-1]]
            else:
                ke = [(a*a+b*b)/4 for a, b in zip(velocity, velocity[1:])]
                new_ke = [(a*a+b*b)/4 for a, b in zip(new_velocity, new_velocity[1:])]
            pe = [-RNG.random()*10 for _ in old[:-1]]
            new_pe = [-RNG.random()*10 for _ in new[:-1]]
            avg_e, avg_ke, avg_pe = [overlap_average(old, new, x) for x in (energy, ke, pe)]
            new_w = remap_squared(old, new, w)
            new_e = [e+k+p-kn-pn for e, k, p, kn, pn in zip(avg_e, avg_ke, avg_pe, new_ke, new_pe)]
            old_total = math.fsum((b-a)*(e+k+p+v*v) for a, b, e, k, p, v in
                                 zip(old, old[1:], energy, ke, pe, w))
            new_total = math.fsum((b-a)*(e+k+p+v*v) for a, b, e, k, p, v in
                                 zip(new, new[1:], new_e, new_ke, new_pe, new_w))
            close(new_total, old_total)


def check_sources():
    source = (ROOT/'star/private/mesh_adjust.f90').read_text()
    for kind in ('u', 'etrb'):
        # Connect the numerical check to the actual call's staggered coordinates.
        routine = source.split('subroutine do_'+kind+'(')[1].split('end subroutine')[0]
        assert 'comes_from, old_xq, new_xq, &' in routine
        assert 'xout_old' not in routine and 'xout_new' not in routine
        routine = source.split('subroutine adjust1_'+kind+'(')[1].split('end subroutine')[0]
        assert 'xq_outer = new_xq(k)' in routine and 'xq1 = old_xq(kk+1)' in routine
        assert 'dq = dq - (dq_sum - new_cell_dq)' in routine
    face = source.split('subroutine do_RSP2_face_var(')[1].split('end subroutine')[0]
    assert 'call interpolate_rsp2_face(' in face
    hydro = (ROOT/'star/private/hydro_rsp2.f90').read_text()
    interp = hydro.split('subroutine interpolate_rsp2_face(')[1].split('end subroutine')[0]
    assert 'call interpolate_vector_pm(' in interp
    assert 'face_new(k) = face_old(j)' in interp
    assert 'face_new(k) = face_old(j+1)' in interp
    assert 'min(face_old(j),face_old(j+1))' in interp
    assert 'max(face_old(j),face_old(j+1))' in interp
    assert "'invalid RSP2 Phi before remesh'" in interp
    assert 'face_old_plus1(nz_old+1) = face_old_plus1(nz_old)' in face
    assert 'xh(i_var,1) = 0d0' in face
    assert 's, s% i_Y, nz, nz_old' in source
    assert 'max(0' not in face and 'abs(' not in face
    normal = source.split('subroutine do_mesh_adjust(')[1].split('end subroutine')[0]
    assert normal.index('call do_etrb(') < normal.index('call remesh_rsp2_moments(s, nz, dq, xh, ierr)')
    envelope = (ROOT/'star/private/tdc_hydro_support.f90').read_text()
    assert envelope.index('call remap1_cell_average2') < envelope.index('call remesh_rsp2_moments(s, nz, s%dq, s%xh, ierr)')
    assert 'call interpolate_rsp2_face(' in envelope
    cleanup = source.split('subroutine dealloc\n')[1].split('end subroutine')[0]
    assert 'call do_work_arrays(.false.,ierr_dealloc)' in cleanup
    assert 'if (ierr == 0) ierr = ierr_dealloc' in cleanup
    amr = (ROOT/'star/private/adjust_mesh_split_merge.f90').read_text()
    assert amr.index('call amr(s,ierr)') < amr.index('call remesh_rsp2_moments(s, s% nz, s% dq, s% xh, ierr)')
    assert 's% Pi(1:s% nz) = s% xh(s% i_Pi,1:s% nz)' in amr
    assert 's% Phi(1:s% nz) = s% xh(s% i_Phi,1:s% nz)' in amr
    defaults = (ROOT/'star/defaults/controls_dev.defaults').read_text()
    assert re.search(r'^\s*RSP2_remesh_when_load\s*=\s*\.false\.', defaults, re.M)
    energy = (ROOT/'star/private/star_utils.f90').read_text()
    partial = energy.split('subroutine eval_deltaM_total_energy_integrals(')[1].split('end subroutine')[0]
    assert partial.index('sum_dm = sum_dm + dm') > partial.index('dm = deltaM - sum_dm')


def check_partial_mass_energy():
    # Compare the loop to independent geometric integration through a partial
    # final cell, as used by mass-change energy accounting.
    dm = [.2, .5, .3]
    faces = [0, .2, .7, 1]
    thermal, kinetic, potential, turbulent = [11, 13, 17], [2, 3, 5], [-9, -8, -7], [1, 4, 9]
    for limit in (0, .1, .2, .4, .7, .9, 1, 2):
        for values in (thermal, kinetic, potential, turbulent):
            total, used = 0.0, 0.0
            for mass, value in zip(dm, values):
                if used >= limit:
                    break
                if used+mass > limit:
                    mass = limit-used
                used += mass
                total += mass*value
            expected = math.fsum(max(0, min(b, limit)-a)*v
                                 for a, b, v in zip(faces, faces[1:], values))
            close(total, expected)


if __name__ == '__main__':
    for check in (check_remap, check_energy, check_partial_mass_energy, check_sources):
        check()
        print(check.__name__, 'PASS')

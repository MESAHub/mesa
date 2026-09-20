"""Isolated RSP2 transport algebra checks; does not compile or run MESA."""
import numpy as np


def face_flux(w, dm, mass_interp=True):
    """Outward Lt and its derivatives, with positive geometric factor = 1."""
    outer, inner = w
    alfa = dm[0]/sum(dm) if mass_interp else 0.5
    beta = 1-alfa
    face = alfa*inner + beta*outer
    difference = inner**2-outer**2
    return face*difference, np.array([
        beta*difference-2*outer*face,
        alfa*difference+2*inner*face,
    ])


def residual_jacobian(w, start, dm, dt, theta, mass_interp=True):
    """Original nonlinear, centered energy rows; no predictor approximation."""
    current, derivative = face_flux(w, dm, mass_interp)
    old, _ = face_flux(start, dm, mass_interp)
    signs = np.array([-1., 1.])/dm
    residual = w*w-start*start + dt*signs*(theta*current + (1-theta)*old)
    jacobian = np.diag(2*w) + dt*theta*np.outer(signs, derivative)
    return residual, jacobian


def predict(start, dm, dt, theta, mass_interp=True, current=None):
    """Retain quadratic storage and linearize the current transport term."""
    if current is None:
        current = start.copy()
    old_flux, _ = face_flux(start, dm, mass_interp)
    flux, derivative = face_flux(current, dm, mass_interp)
    old_incoming = dt*np.array([old_flux, -old_flux])/dm
    incoming = dt*theta*np.array([flux, -flux])/dm
    detrb_dw = dt*theta*np.array([derivative[0], -derivative[1]])/dm
    guess = current.copy()
    for k in range(2):
        if old_incoming[k] <= start[k]**2:
            continue
        etrb = (start[k]**2 + (1-theta)*old_incoming[k] + incoming[k]
                - detrb_dw[k]*current[k])
        if etrb <= 0:
            continue
        discr = np.sqrt(detrb_dw[k]**2 + 4*etrb)
        soln = ((detrb_dw[k]+discr)/2 if detrb_dw[k] > 0
                else 2*etrb/(discr-detrb_dw[k]))
        guess[k] = max(guess[k], min(max(start), soln))
    return guess


def solve(start, dm, dt, theta, mass_interp=True, guess=None):
    w = start.copy() if guess is None else guess.copy()
    for iteration in range(25):
        residual, jacobian = residual_jacobian(w, start, dm, dt, theta, mass_interp)
        if max(abs(residual)) < 1e-13:
            return w, iteration
        change = np.linalg.solve(jacobian, -residual)
        next_w = np.maximum(w+change, 0.)
        if np.array_equal(next_w, w):
            return None, iteration+1
        w = next_w
    return None, 25


def check_onset():
    start, dm = np.array([0., 1.]), np.ones(2)
    assert solve(start, dm, .1, .5)[0] is None
    root, iterations = solve(start, dm, .1, .5, guess=predict(start, dm, .1, .5))
    assert np.allclose(root, [0.22785876226600094, 0.9736941945285522],
                       rtol=0, atol=1e-13)
    print(f'Zero-w counterexample: stalls without predictor, converges in {iterations} updates with it.')

    # A flux-only square-root guess fails for a narrow receiving cell.
    dm = np.array([1e-4, 1.])
    dt = .01*min(dm)
    old_flux, _ = face_flux(start, dm)
    naive = np.array([np.sqrt(dt*old_flux/dm[0]), 1.])
    assert solve(start, dm, dt, .5, guess=naive)[0] is None
    assert solve(start, dm, dt, .5, guess=predict(start, dm, dt, .5))[0] is not None

    count, worst = 0, 0
    for ratio in (1e-4, .01, .1, 1., 10., 100., 1e4):
        dm = np.array([ratio, 1.])
        for reverse in (False, True):
            for quiet in (0., 1e-14, 1e-8):
                start = np.array([quiet, 1.])
                if reverse:
                    start = start[::-1].copy()
                for mass_interp in (False, True):
                    for theta in (.5, 1.):
                        for step in (1e-6, .01, .1, .5):
                            dt = step*min(dm)
                            guess = predict(start, dm, dt, theta, mass_interp)
                            root, iterations = solve(start, dm, dt, theta, mass_interp, guess)
                            assert root is not None, (ratio, reverse, quiet, theta, step)
                            residual, _ = residual_jacobian(root, start, dm, dt, theta, mass_interp)
                            assert max(abs(residual)) < 1e-13
                            # Conservation is checked at the converged root, not at the initial guess.
                            assert abs(np.dot(dm, root*root-start*start)) < 1e-12*sum(dm)
                            count += 1
                            worst = max(worst, iterations)
    print(f'{count} onset cases passed; maximum Newton updates: {worst}.')


def check_derivatives_and_bounds():
    rng = np.random.default_rng(20260919)
    for _ in range(500):
        dm = 10**rng.uniform(-2, 2, 2)
        w = rng.uniform(.05, 2, 2)
        dt, theta = rng.uniform(.01, .5)*min(dm), rng.choice([.5, 1.])
        residual, jacobian = residual_jacobian(w, np.ones(2), dm, dt, theta)
        for j in range(2):
            upper, lower = w.copy(), w.copy()
            upper[j] += 1e-6
            lower[j] -= 1e-6
            finite = (residual_jacobian(upper, np.ones(2), dm, dt, theta)[0]
                      - residual_jacobian(lower, np.ones(2), dm, dt, theta)[0])/2e-6
            assert np.allclose(finite, jacobian[:, j], rtol=2e-7, atol=2e-9)
        start, saved = w.copy(), w.copy()
        guess = predict(start, dm, dt, theta)
        assert np.array_equal(start, saved)
        assert np.all(guess >= start)
        assert np.all(guess <= max(start))
        assert np.array_equal(predict(start, dm, 0., theta), start)
        assert np.array_equal(predict(np.ones(2), dm, dt, theta), np.ones(2))
    # Long-step bounding applies to the guess, not the physical flux or final root.
    start = np.array([0., 1.])
    assert np.array_equal(predict(start, np.ones(2), 1e8, .5), np.ones(2))
    print('Finite-difference Jacobians, unchanged starting state, and predictor bounds passed.')


if __name__ == '__main__':
    check_onset()
    check_derivatives_and_bounds()
    print('Algebra only: no full stellar coupling, source terms, MESA AD execution, or evolution tested.')

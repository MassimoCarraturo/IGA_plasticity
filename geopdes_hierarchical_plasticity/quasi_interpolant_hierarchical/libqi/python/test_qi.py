"""End-to-end sanity check of the Python wrapper for libqi."""

import os, sys
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from qi import HSpace, compute_coeffs, eval_tp, version


def halton(i, b):
    f, r = 1.0, 0.0
    while i > 0:
        f /= b
        r += f * (i % b)
        i //= b
    return r


def f_smooth(x, y):
    return ((np.tanh(9*y - 9*x) + 1)/9
            + np.exp(-((10*x - 6)**2 + (10*y + 7)**2)) / 1.5)


def test_polynomial_reproduction():
    n = 2000
    xs = np.array([halton(i+1, 2) for i in range(n)])
    ys = np.array([halton(i+1, 3) for i in range(n)])
    f = 0.7 - 1.3*xs + 0.5*ys + 2.1*xs**2 - 1.7*xs*ys + 0.9*ys**2
    hs = HSpace(degree=(2, 2), nel=(8, 8))
    c = compute_coeffs(hs, xs, ys, f)

    G = 50
    gx, gy = np.meshgrid(np.linspace(0,1,G), np.linspace(0,1,G), indexing='ij')
    fhat = eval_tp(hs, 0, c, gx.ravel(), gy.ravel()).reshape(gx.shape)
    exact = 0.7 - 1.3*gx + 0.5*gy + 2.1*gx**2 - 1.7*gx*gy + 0.9*gy**2
    err = np.max(np.abs(fhat - exact))
    print(f"  polynomial reproduction max error = {err:.3e}")
    assert err < 1e-9, "polynomial reproduction FAILED"


def test_smooth():
    n = 4000
    xs = np.array([halton(i+1, 2) for i in range(n)])
    ys = np.array([halton(i+1, 3) for i in range(n)])
    f = f_smooth(xs, ys)
    hs = HSpace(degree=(3, 3), nel=(16, 16))
    c, diag = compute_coeffs(hs, xs, ys, f, return_diag=True)
    G = 100
    gx, gy = np.meshgrid(np.linspace(0,1,G), np.linspace(0,1,G), indexing='ij')
    fhat = eval_tp(hs, 0, c, gx.ravel(), gy.ravel()).reshape(gx.shape)
    err = np.max(np.abs(fhat - f_smooth(gx, gy)))
    print(f"  smooth max error          = {err:.3e}")
    print(f"  degrees used (unique)     = {sorted(set(diag['deg_used']))}")
    assert err < 5e-2, "smooth FAILED"


def test_vector():
    n = 2000
    xs = np.array([halton(i+1, 2) for i in range(n)])
    ys = np.array([halton(i+1, 3) for i in range(n)])
    f = np.column_stack([
        0.7 - 1.3*xs + 0.5*ys + 2.1*xs**2 - 1.7*xs*ys + 0.9*ys**2,
        1 + 2*xs - 3*ys + xs*ys,
        -0.5 + 1.1*xs**2 - 0.4*ys**2,
    ])
    hs = HSpace(degree=(2, 2), nel=(8, 8))
    c = compute_coeffs(hs, xs, ys, f)
    assert c.shape == (hs.ndof, 3), f"shape mismatch: {c.shape}"
    print(f"  vector coeffs shape       = {c.shape}")


if __name__ == "__main__":
    print(f"libqi Python wrapper -- {version()}")
    test_polynomial_reproduction()
    test_smooth()
    test_vector()
    print("ALL PYTHON TESTS PASSED")

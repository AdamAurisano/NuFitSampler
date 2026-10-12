"""Python mirror of the C++ surrogate posterior (NuFitModel), plus exact profiles and
quadrature marginals computed by variable elimination over its factors.

Per ordering o the surrogate log-likelihood (-Delta chi^2 / 2) is a sum of factors on
the NuFIT grids:
    A(s23, dma, dcp)            = -chi2_3D / 2
    B(s12, l21)                 = -(chi2_2D - m_o) / 2
    C(s13)                      = -(chi2_1D - m_o) / 2                 (factorized)
    D(s13, s23)                 = -(chi2_{T13/T23} - chi2_{T23}) / 2   (coupled)
Every factor is N-linear on its grid and the factors share at most one axis, so:
  * profiles (max over the other parameters) are attained at grid nodes and are computed
    exactly by max-sum elimination on the node grids;
  * marginals are computed by sum-product elimination on a refined midpoint grid
    (sub-cell midpoint rule, O(h^2) accurate).
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.special import logsumexp

from .tables import VARS, load_variant

LN10 = np.log(10.0)
ORDERINGS = ("NO", "IO")


# ----------------------------------------------------------------------------- factors
class Factor:
    def __init__(self, vars_, arr):
        self.vars, self.arr = tuple(vars_), np.asarray(arr, dtype=float)
        assert self.arr.ndim == len(self.vars)


def _expand(f, order):
    """Broadcastable view of factor f in the variable order `order`."""
    perm = [f.vars.index(v) for v in order if v in f.vars]
    a = np.transpose(f.arr, perm)
    shape = [f.arr.shape[f.vars.index(v)] if v in f.vars else 1 for v in order]
    return a.reshape(shape)


def _combine(factors):
    order = []
    for f in factors:
        order += [v for v in f.vars if v not in order]
    total = 0.0
    for f in factors:
        total = total + _expand(f, order)
    return Factor(order, np.broadcast_to(total, [_dim(factors, v) for v in order]).copy())


def _dim(factors, v):
    for f in factors:
        if v in f.vars:
            return f.arr.shape[f.vars.index(v)]
    raise KeyError(v)


def eliminate(factors, keep, mode, logw=None):
    """Max-sum (mode='max') or log-sum-exp (mode='sum', weights exp(logw[v])) elimination
    of all variables not in `keep`. Returns Factor over `keep` (in that order)."""
    factors = list(factors)
    allv = {v for f in factors for v in f.vars}
    todo = [v for v in allv if v not in keep]
    while todo:
        def cost(v):
            inv = [f for f in factors if v in f.vars]
            u = {w for f in inv for w in f.vars}
            return np.prod([_dim(inv, w) for w in u], dtype=float)
        v = min(todo, key=cost)
        todo.remove(v)
        inv = [f for f in factors if v in f.vars]
        factors = [f for f in factors if v not in f.vars]
        F = _combine(inv)
        ax = F.vars.index(v)
        if mode == "max":
            red = F.arr.max(axis=ax)
        else:
            w = logw[v].reshape([-1 if i == ax else 1 for i in range(F.arr.ndim)])
            red = logsumexp(F.arr + w, axis=ax)
        factors.append(Factor([u for u in F.vars if u != v], red))
    scal = sum(float(f.arr) for f in factors if not f.vars)
    factors = [f for f in factors if f.vars]
    if not keep:
        return Factor([], np.array(scal))
    F = _combine(factors)
    return Factor(keep, np.transpose(F.arr, [F.vars.index(v) for v in keep]) + scal)


def _refine(arr, ax, r):
    """Values at the r sub-cell midpoints of every cell along axis ax (linear interpolation)."""
    n = arr.shape[ax]
    lo = np.take(arr, np.arange(n - 1), axis=ax)
    hi = np.take(arr, np.arange(1, n), axis=ax)
    parts = [lo * (1 - (k + 0.5) / r) + hi * ((k + 0.5) / r) for k in range(r)]
    st = np.stack(parts, axis=ax + 1)
    shape = list(arr.shape)
    shape[ax] = (n - 1) * r
    return st.reshape(shape)


def _midpoints(axis, r):
    lo, hi = axis[:-1], axis[1:]
    x = np.concatenate([[lo + (hi - lo) * (k + 0.5) / r] for k in range(r)]).T.ravel()
    w = np.repeat((hi - lo) / r, r)
    return x, w


# --------------------------------------------------------------------------- surrogate
class Surrogate:
    def __init__(self, datadir, variant="TBoff", coupled=False, angle_prior="sin2", dm21_prior="linear"):
        self.variant, self.coupled = variant, coupled
        self.angle_prior, self.dm21_prior = angle_prior, dm21_prior
        self.tables = load_variant(datadir, variant)
        self.m = {o: self.tables[o]["T13"].values.min() for o in ORDERINGS}
        self._interp = {}
        for o in ORDERINGS:
            t = self.tables[o]
            for tag in ("T23/DMA/DCP", "T12/DMS", "T13", "T13/T23", "T23"):
                self._interp[o, tag] = RegularGridInterpolator(
                    tuple(t[tag].axes), t[tag].values, bounds_error=False, fill_value=np.nan)

    # ---- pointwise evaluation (mirrors NuFitModel::deltaChi2 / logPrior) -------------
    def delta_chi2(self, P):
        """P: (N, 6) array of (s12, s13, s23, dcp, dm3l[1e-3 eV^2], dm21[1e-5 eV^2])."""
        P = np.atleast_2d(P)
        out = np.full(len(P), np.nan)
        dcp = (P[:, 3] + 180.0) % 360.0 - 180.0
        l21 = np.log10(P[:, 5] * 1e-5)
        for o, sel in (("NO", P[:, 4] > 0), ("IO", P[:, 4] < 0)):
            if not sel.any():
                continue
            q = P[sel]
            c = self._interp[o, "T23/DMA/DCP"](np.c_[q[:, 2], q[:, 4], dcp[sel]])
            c += self._interp[o, "T12/DMS"](np.c_[q[:, 0], l21[sel]]) - self.m[o]
            if self.coupled:
                c += self._interp[o, "T13/T23"](np.c_[q[:, 1], q[:, 2]]) - self._interp[o, "T23"](q[:, 2:3])
            else:
                c += self._interp[o, "T13"](q[:, 1:2]) - self.m[o]
            out[sel] = c
        return out

    def log_prior(self, P):
        P = np.atleast_2d(P)
        lp = np.zeros(len(P))
        if self.angle_prior == "angle":
            for k in range(3):
                lp += -0.5 * np.log(P[:, k] * (1 - P[:, k])) - np.log(2.0)
        if self.dm21_prior == "log":
            lp -= np.log(P[:, 5])
        return lp

    # ---- factors ----------------------------------------------------------------------
    def factors(self, o):
        t = self.tables[o]
        f = [Factor(("s23", "dma", "dcp"), -0.5 * t["T23/DMA/DCP"].values),
             Factor(("s12", "l21"), -0.5 * (t["T12/DMS"].values - self.m[o]))]
        if self.coupled:
            f.append(Factor(("s13", "s23"), -0.5 * (t["T13/T23"].values - t["T23"].values[None, :])))
        else:
            f.append(Factor(("s13",), -0.5 * (t["T13"].values - self.m[o])))
        return f

    def axes(self, o):
        t = self.tables[o]
        ax = dict(zip(("s23", "dma", "dcp"), t["T23/DMA/DCP"].axes))
        ax.update(zip(("s12", "l21"), t["T12/DMS"].axes))
        ax["s13"] = t["T13"].axes[0]
        return ax

    def log_prior_axis(self, v, x):
        """Log prior density of variable v in table coordinates."""
        if v in ("s12", "s13", "s23") and self.angle_prior == "angle":
            with np.errstate(divide="ignore"):
                return -0.5 * np.log(x * (1 - x)) - np.log(2.0)
        if v == "l21" and self.dm21_prior == "linear":
            return x * LN10
        return np.zeros_like(x)

    # ---- profiles ----------------------------------------------------------------------
    def profile(self, keep, o=None):
        """Surrogate Delta chi^2 profile on the node grid of the kept variables.
        o=None minimizes over the orderings (only for keep without 'dma')."""
        if o is None:
            assert "dma" not in keep
            return np.minimum(self.profile(keep, "NO"), self.profile(keep, "IO"))
        return -2.0 * eliminate(self.factors(o), list(keep), "max").arr

    # ---- marginals by quadrature --------------------------------------------------------
    def _trimmed(self, o, r, dchi2_max=60.0):
        """Refined factors restricted to the region where the 1D profiles are below
        dchi2_max; returns (factors, logw, midpoint coords, kept node axes)."""
        ax = self.axes(o)
        fac = self.factors(o)
        keep_idx = {}
        for v in VARS:
            v_ = v if v != "dma" else "dma"
            prof = -2.0 * eliminate(fac, [v_], "max").arr
            ok = np.nonzero(prof < dchi2_max + self.m[o])[0]
            keep_idx[v] = np.arange(max(ok[0] - 1, 0), min(ok[-1] + 2, len(ax[v])))
        out, logw, mids, nodes = [], {}, {}, {}
        for f in fac:
            a = f.arr
            for d, v in enumerate(f.vars):
                a = np.take(a, keep_idx[v], axis=d)
            for d in range(a.ndim):
                a = _refine(a, d, r)
            out.append(Factor(f.vars, a))
        for v in VARS:
            nodes[v] = ax[v][keep_idx[v]]
            x, w = _midpoints(nodes[v], r)
            mids[v] = x
            with np.errstate(divide="ignore"):
                logw[v] = np.log(w) + self.log_prior_axis(v, x)
        return out, logw, mids, nodes

    def marginals(self, r=4):
        """Per ordering: log evidence and, for each variable, the probability mass in each
        table cell (node-to-node bin). Returns dict with keys 'logZ', 'pIO', 'cells'."""
        res = {"logZ": {}, "cells": {}}
        for o in ORDERINGS:
            fac, logw, mids, nodes = self._trimmed(o, r)
            logZ = float(eliminate(fac, [], "sum", logw).arr)
            res["logZ"][o] = logZ
            res["cells"][o] = {}
            for v in VARS:
                lm = eliminate(fac, [v], "sum", logw).arr + logw[v]
                mass = np.exp(lm - logZ).reshape(-1, r).sum(axis=1)
                res["cells"][o][v] = (nodes[v], mass)
        zmax = max(res["logZ"].values())
        zn, zi = (np.exp(res["logZ"][o] - zmax) for o in ORDERINGS)
        res["pIO"] = zi / (zn + zi)
        return res


def samples_to_table_coords(P):
    """(N, 6) sample parameters -> dict of table-coordinate arrays."""
    return {"s12": P[:, 0], "s13": P[:, 1], "s23": P[:, 2],
            "dcp": (P[:, 3] + 180.0) % 360.0 - 180.0, "dma": P[:, 4],
            "l21": np.log10(P[:, 5] * 1e-5)}

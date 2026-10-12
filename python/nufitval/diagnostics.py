"""MCMC convergence diagnostics (Vehtari, Gelman, Simpson, Carpenter & Buerkner 2021).

chains: array of shape (n_chains, n_draws) for one scalar quantity.
"""

import numpy as np
from scipy.stats import norm, rankdata


def _split(chains):
    n = chains.shape[1] // 2
    return np.concatenate([chains[:, :n], chains[:, -n:]], axis=0)


def _rank_normalize(chains):
    r = rankdata(chains, method="average").reshape(chains.shape)
    return norm.ppf((r - 3.0 / 8) / (chains.size + 1.0 / 4))


def _rhat_basic(chains):
    m, n = chains.shape
    if n < 2:
        return np.nan
    w = chains.var(axis=1, ddof=1).mean()
    b = n * chains.mean(axis=1).var(ddof=1)
    if w == 0:
        return 1.0 if b == 0 else np.inf
    return np.sqrt(((n - 1) / n * w + b / n) / w)


def rhat(chains):
    """Rank-normalized split R-hat (max of bulk and folded/tail versions)."""
    s = _split(np.asarray(chains, float))
    bulk = _rhat_basic(_rank_normalize(s))
    tail = _rhat_basic(_rank_normalize(np.abs(s - np.median(s))))
    return max(bulk, tail)


def autocorr(x, max_lag=None):
    """Normalized autocorrelation of a 1D series (FFT)."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n]
    ac = ac / ac[0] if ac[0] > 0 else np.zeros(n)
    return ac[:max_lag] if max_lag else ac


def _ess_basic(chains):
    """ESS with Geyer's initial monotone sequence, combining chains as in Stan."""
    m, n = chains.shape
    if np.all(chains == chains.flat[0]):
        return np.nan
    acov = np.array([autocorr(c) * c.var() for c in chains])
    w = chains.var(axis=1, ddof=1).mean()
    var_plus = w * (n - 1) / n + (chains.mean(axis=1).var(ddof=1) if m > 1 else 0)
    rho = 1 - (w - acov.mean(axis=0)) / var_plus
    rho[0] = 1
    # Geyer: sum pairs while positive, enforce monotone decrease
    t, prev, total = 0, np.inf, 0.0
    while t + 1 < n:
        p = rho[t] + rho[t + 1]
        if p < 0:
            break
        p = min(p, prev)
        total += p
        prev = p
        t += 2
    tau = -1 + 2 * total
    return m * n / max(tau, 1.0 / np.log10(m * n))


def ess_bulk(chains):
    return _ess_basic(_rank_normalize(_split(np.asarray(chains, float))))


def ess_tail(chains):
    s = _split(np.asarray(chains, float))
    out = []
    for q in (0.05, 0.95):
        iq = (s <= np.quantile(s, q)).astype(float)
        out.append(_ess_basic(iq))
    return np.nanmin(out)


def ess_mean(chains):
    return _ess_basic(np.asarray(chains, float))


def tau_int(chains):
    """Integrated autocorrelation time in units of stored draws."""
    c = np.asarray(chains, float)
    return c.size / ess_mean(c)


def mcse_mean(chains):
    c = np.asarray(chains, float)
    return c.std(ddof=1) / np.sqrt(ess_mean(c))


def round_trip_rate(labels, K):
    """labels: (n_draws, K) rung index of each label; returns round trips per stored draw."""
    trips, direction = 0, np.zeros(K, int)
    for row in labels:
        lb, lt = np.argmin(row), np.argmax(row)
        if direction[lb] == -1:
            trips += 1
        direction[lb] = 1
        if direction[lt] == 1:
            direction[lt] = -1
    return trips / len(labels)

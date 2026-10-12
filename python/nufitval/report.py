"""Validation report for the NuFIT surrogate samplers.

    cd python
    python -m nufitval.report --variant TBoff --surrogate factorized \
        --results ../results --data ../data/nufit --out ../results/report_TBoff_factorized

Inputs (written by build/sample_direct and build/sample_pt, see README):
    direct_<variant>_<surrogate>_samples.npy          required
    pt_<variant>_<surrogate>_{cold,rungs,labels,ladder}.npy, _meta.json   optional
    pt_<variant>_<surrogate>_nojumps_cold.npy          optional (reducibility demo)
    legacy_samples.npy                                 optional (student's chain)
    direct_<variant>_factorized_angleprior_samples.npy optional (prior sensitivity)
Outputs: report.pdf, one PNG per page, summary.json.
"""

import argparse
import json
import os

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter
from scipy.stats import chi2 as chi2dist

from . import diagnostics as D
from .style import (CMAP, DIRECT, DIVERGING, EXTRA, GRID, INK, INK2, LEGACY, LEVELS_2D, MASSES,
                    MUTED, PT, note, page, plt)
from .surrogate import ORDERINGS, Surrogate, samples_to_table_coords
from .tables import LABELS, VARS

REPLICA_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
PAR_NAMES = ["s12", "s13", "s23", "dcp", "dm3l", "dm21"]
PAR_LABELS = {
    "s12": LABELS["s12"], "s13": LABELS["s13"], "s23": LABELS["s23"], "dcp": LABELS["dcp"],
    "dm3l": LABELS["dma"], "dm21": r"$\Delta m^2_{21}$ [$10^{-5}$ eV$^2$]", "logL": r"$-\Delta\chi^2/2$",
}


# ------------------------------------------------------------------------------ helpers
def load_optional(path):
    return np.load(path) if os.path.exists(path) else None


def hpd_levels(H, masses):
    """Density thresholds of a 2D histogram enclosing the given probability masses."""
    h = np.sort(H.ravel())[::-1]
    c = np.cumsum(h) / h.sum()
    return [h[min(np.searchsorted(c, m), len(h) - 1)] for m in masses]


def region_box(table_vals, axes, level, pad=0.15):
    """Bounding box of {Delta chi^2 < level} on a 2D table, padded."""
    ok = np.argwhere(table_vals < level)
    if len(ok) == 0:
        return [(a[0], a[-1]) for a in axes]
    box = []
    for d, a in enumerate(axes):
        lo, hi = a[ok[:, d].min()], a[ok[:, d].max()]
        w = hi - lo if hi > lo else (a[-1] - a[0]) * 0.05
        box.append((max(a[0], lo - pad * w), min(a[-1], hi + pad * w)))
    return box


def quantities(cold):
    """Scalar chains used for convergence diagnostics, from PT cold rows."""
    P = cold[:, 2:8]
    rad = np.radians(P[:, 3])
    return {
        "s12": P[:, 0], "s13": P[:, 1], "s23": P[:, 2], "cos dcp": np.cos(rad), "sin dcp": np.sin(rad),
        "|dm3l|": np.abs(P[:, 4]), "dm21": P[:, 5], "IO": (P[:, 4] < 0).astype(float),
        "upper octant": (P[:, 2] > 0.5).astype(float), "logL": cold[:, 8],
    }


def by_replica(cold, x):
    reps = np.unique(cold[:, 0]).astype(int)
    n = min((cold[:, 0] == r).sum() for r in reps)
    return np.stack([x[cold[:, 0] == r][:n] for r in reps])


def merged_edges(nodes, mass, target_bins=50, floor=1e-7):
    """Coarsen node-to-node cells to about target_bins bins over the populated range."""
    ok = np.nonzero(mass > floor)[0]
    lo, hi = ok[0], ok[-1] + 1
    k = max(1, int(np.ceil((hi - lo) / target_bins)))
    idx = list(range(lo, hi + 1, k))
    if idx[-1] != hi:
        idx.append(hi)
    edges = nodes[idx]
    cm = np.concatenate([[0], np.cumsum(mass)])
    return edges, np.diff(cm[idx])


def binned_chi2(counts, expected_prob, n, n_eff=None, min_exp=20):
    """Pearson chi^2 of histogram counts against expected probabilities. With n_eff
    (autocorrelated chains) the statistic is scaled by n_eff/n."""
    e = expected_prob * n
    sel = e > min_exp
    stat = ((counts[sel] - e[sel]) ** 2 / e[sel]).sum()
    if n_eff is not None:
        stat *= n_eff / n
    ndf = max(sel.sum() - 1, 1)
    return stat, ndf, chi2dist.sf(stat, ndf)


# ------------------------------------------------------------------------------- pages
def page_summary(s, pages):
    """Headline numbers; the full set is in summary.json."""
    fig, ax = page(1, 1, f"Validation summary: {s['model']}")
    ax = ax[0, 0]
    ax.axis("off")
    L = []
    sec = lambda t: L.append(("h", t))  # noqa: E731
    row = lambda k, v, ok=None: L.append(("r", k, v, ok))  # noqa: E731

    sec("Implementation")
    row("max |logL(Python mirror) - logL(C++)|", f"{s['max_abs_logL_py_minus_cpp']:.1e}",
        s["max_abs_logL_py_minus_cpp"] < 1e-9)
    sec("Exact (direct) sampler vs quadrature of the surrogate")
    pull = (s["P(IO) direct"] - s["P(IO) quadrature"]) / s["P(IO) direct binomial error"]
    row("P(IO): quadrature / direct (pull)", f"{s['P(IO) quadrature']:.4f} / {s['P(IO) direct']:.4f} ({pull:+.1f} sd)",
        abs(pull) < 3)
    pd = [v["direct"]["p"] for v in s["marginal_chi2"].values()]
    row(f"1D marginal chi^2 p-values ({len(pd)} tests)", " ".join(f"{p:.2g}" for p in pd), min(pd) > 0.001 / len(pd))
    row("P(upper octant)", f"{s['P(upper octant) direct']:.4f}")
    if "pt_diagnostics" in s:
        sec("Parallel tempering (cold chain, all replicas)")
        dg = s["pt_diagnostics"]
        rh = max(v["rhat"] for v in dg.values())
        eb = min(v["ess_bulk"] for v in dg.values())
        et = min(v["ess_tail"] for v in dg.values())
        row("max R-hat (rank-normalized split)", f"{rh:.4f}", rh < 1.01)
        row("min ESS bulk / tail", f"{eb:.0f} / {et:.0f}", min(eb, et) > 400)
        row("P(IO) PT +- MCSE (pull vs exact)",
            f"{s['P(IO) PT']:.4f} +- {s['P(IO) PT MCSE']:.4f} ({s['P(IO) PT pull vs exact']:+.1f} sd)",
            abs(s["P(IO) PT pull vs exact"]) < 3)
        pp = [v["PT"]["p"] for v in s["marginal_chi2"].values() if "PT" in v]
        row(f"PT 1D marginal chi^2 p-values (ESS-scaled)", " ".join(f"{p:.2g}" for p in pp), min(pp) > 0.001 / len(pp))
        if "P(IO) PT without jumps, per replica" in s:
            row("P(IO) per replica WITHOUT mode jumps", " ".join(f"{x:.2f}" for x in s["P(IO) PT without jumps, per replica"])
                + f"   (R-hat {s['R-hat(IO) PT without jumps']:.3f})", False)
    sec("Interpolation (decimated-grid test, linear, nodes inside 3 sigma)")
    for t in ("T23/DMA/DCP", "T12/DMS", "T13"):
        a, b = s["interp_loo"][f"{t} NO linear"], s["interp_loo"][f"{t} NO pchip"]
        row(f"{t} (NO): rms / max error, linear vs pchip", f"{a['rms']:.3f} / {a['max']:.2f}  vs  {b['rms']:.3f} / {b['max']:.2f}")
    sec("Surrogate vs NuFIT 2D tables not used to build it (max |Delta chi^2 residual| inside 3 sigma)")
    for o in ORDERINGS:
        r = s[f"residuals_{o}"]
        held = sorted(((v["max_abs_in_3sigma"], k) for k, v in r.items() if not v["by_construction"]), reverse=True)
        row(f"{o}: largest held-out residuals", ",  ".join(f"{k} {v:.2f}" for v, k in held[:4]),
            held[0][0] < 1.0)
        ex = max(v["max_abs_in_3sigma"] for v in r.values() if v["by_construction"])
        row(f"{o}: projections exact by construction (max residual)", f"{ex:.1e}", ex < 1e-6)

    y = 0.98
    for item in L:
        if item[0] == "h":
            y -= 0.018
            ax.text(0, y, item[1], fontsize=9.5, weight="bold", color=INK, transform=ax.transAxes)
            y -= 0.04
        else:
            _, k, v, ok = item
            mark = {True: "PASS", False: "FLAG", None: ""}[ok]
            ax.text(0.01, y, k, fontsize=8.5, color=INK2, transform=ax.transAxes)
            ax.text(0.44, y, v, fontsize=8.5, color=INK, family="monospace", transform=ax.transAxes)
            ax.text(0.985, y, mark, fontsize=8.5, ha="right", weight="bold", transform=ax.transAxes,
                    color={"PASS": "#008300", "FLAG": "#e34948", "": INK}[mark])
            y -= 0.034
    note(fig, "Pages: " + ", ".join(pages) + ".\nPASS threshold for the marginal tests is Bonferroni-corrected "
              "(0.001/n); the exact sampler was also checked with 2e7 samples. FLAG on held-out residuals\nis a property "
              "of the surrogate (lost theta13 correlations), not of the samplers.")
    return fig


def page_interp_loo(S, summary):
    """Leave-one-out: drop every other node along each axis, predict the dropped nodes."""
    fig, axs = page(3, 3, "Interpolation accuracy: decimated-grid (leave-one-out) test")
    methods = [("linear", DIRECT), ("pchip", PT), ("cubic", EXTRA)]
    res = {}
    for row, tag in enumerate(["T23/DMA/DCP", "T12/DMS", "T13"]):
        for col, o in enumerate(ORDERINGS):
            t = S.tables[o][tag]
            m = t.values.min()
            sub_axes = [a[::2] for a in t.axes]
            sub = t.values[tuple(slice(None, None, 2) for _ in t.axes)]
            grids = np.meshgrid(*t.axes, indexing="ij")
            idx = np.meshgrid(*[np.arange(len(a)) for a in t.axes], indexing="ij")
            held = np.zeros(t.values.shape, bool)
            for g in idx:
                held |= g % 2 == 1
            mask = held & (t.values - m < LEVELS_2D[-1])
            pts = np.stack([g[mask] for g in grids], axis=-1)
            truth = t.values[mask]
            ax = axs[row, col]
            bins = np.logspace(-5, 1.5, 60)
            for meth, colr in methods:
                f = RegularGridInterpolator(sub_axes, sub, method=meth, bounds_error=False, fill_value=np.nan)
                err = f(pts) - truth
                ok = np.isfinite(err)
                a = np.abs(err[ok])
                res[f"{tag} {o} {meth}"] = {"rms": float(np.sqrt(np.mean(a ** 2))),
                                            "p95": float(np.quantile(a, 0.95)), "max": float(a.max())}
                ax.hist(np.clip(a, 1e-5, None), bins=bins, histtype="step", color=colr, lw=1.6,
                        label=f"{meth}: rms {np.sqrt(np.mean(a**2)):.3f}, max {a.max():.2f}")
                if meth == "linear" and col == 0 and row == 0:
                    lin_err, lin_pts = err, pts
            ax.set_xscale("log")
            ax.set_xlabel(r"|interpolated $-$ true| $\Delta\chi^2$ at held-out nodes")
            ax.set_title(f"{tag} ({o}), nodes with $\\Delta\\chi^2<11.83$")
            ax.legend(loc="upper left")
        # third column: where the linear errors sit (3D table, NO; 2D/1D: residual vs coordinate)
        ax = axs[row, 2]
        if row == 0:
            sc = ax.scatter(lin_pts[:, 0], lin_pts[:, 2], c=np.abs(lin_err), s=3, cmap=CMAP, vmin=0,
                            vmax=np.nanquantile(np.abs(lin_err), 0.99))
            fig.colorbar(sc, ax=ax, label=r"|error| (linear)")
            ax.set_xlabel(LABELS["s23"])
            ax.set_ylabel(LABELS["dcp"])
            ax.set_title("linear-interpolation error, T23/DMA/DCP (NO), all DMA")
        else:
            t = S.tables["NO"][tag]
            if row == 1:
                prof = t.values.min(axis=0)
                ax.plot(t.axes[1], prof - t.values.min(), color=INK, lw=1.2)
                ax.plot(t.axes[1], prof - t.values.min(), "o", ms=2.5, color=DIRECT)
                ax.set_xlabel(LABELS["l21"])
                ax.set_xlim(-4.4, -3.9)
            else:
                ax.plot(t.axes[0], t.values - t.values.min(), color=INK, lw=1.2)
                ax.plot(t.axes[0], t.values - t.values.min(), "o", ms=2.5, color=DIRECT)
                ax.set_xlabel(LABELS["s13"])
                ax.set_xlim(0.019, 0.025)
            ax.set_ylim(0, 15)
            ax.set_ylabel(r"$\Delta\chi^2$ (NO) with grid nodes")
            ax.set_title("node spacing near the minimum")
    summary["interp_loo"] = res
    note(fig, "Each held-out node is predicted from a grid with twice the spacing, so these errors overestimate the "
              "real ones (by about 4x for linear, more for cubic). Linear interpolation is exact at the nodes and cannot "
              "overshoot; pchip/cubic are shown for comparison.")
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    return fig


def page_slices(S, best):
    fig, axs = page(2, 3, "Surrogate $\\Delta\\chi^2$ slices through the global best fit (piecewise-linear by construction)")
    ranges = {"s12": (0.26, 0.36), "s13": (0.019, 0.025), "s23": (0.38, 0.65), "dcp": (-180, 180),
              "dm3l": None, "dm21": (6.4, 8.6)}
    for k, (name, ax) in enumerate(zip(PAR_NAMES, axs.ravel())):
        lo, hi = ranges[name] or (abs(best[4]) - 0.15, abs(best[4]) + 0.15)
        x = np.linspace(lo, hi, 2001)
        P = np.repeat(best[None, :], len(x), axis=0)
        P[:, k] = x if name != "dm3l" else np.sign(best[4]) * x
        ax.plot(x, S.delta_chi2(P), color=DIRECT)
        ax.set_xlabel(PAR_LABELS[name] if name != "dm3l" else r"$|\Delta m^2_{3\ell}|$ [$10^{-3}$ eV$^2$]")
        ax.set_ylabel(r"$\Delta\chi^2$")
        ax.set_ylim(0, 12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def page_traces(cold, title):
    fig, axs = page(4, 2, title)
    names = PAR_NAMES + ["logL"]
    reps = np.unique(cold[:, 0]).astype(int)
    for k, (name, ax) in enumerate(zip(names, axs.ravel())):
        col = 2 + k
        for r in reps:
            sel = cold[:, 0] == r
            it, y = cold[sel, 1], cold[sel, col]
            if name == "dm3l":
                y = np.abs(y)  # the sign (ordering) flips constantly; see the mode-mixing page
            step = max(1, len(it) // 4000)
            ax.plot(it[::step], y[::step], lw=0.5, color=REPLICA_COLORS[r % 8], alpha=0.8, label=f"replica {r}")
        ax.set_ylabel(PAR_LABELS[name] if name != "dm3l" else r"$|\Delta m^2_{3\ell}|$", fontsize=8)
        ax.set_xlabel("iteration")
    axs.ravel()[-1].axis("off")
    axs[0, 0].legend(ncol=len(reps), loc="lower left", bbox_to_anchor=(0, 1.02))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def page_diag_table(cold, summary, thin):
    Q = quantities(cold)
    rows = []
    out = {}
    for name, x in Q.items():
        ch = by_replica(cold, x)
        r = D.rhat(ch)
        eb, et = D.ess_bulk(ch), D.ess_tail(ch)
        tau = D.tau_int(ch) * thin
        mc = D.mcse_mean(ch)
        out[name] = {"rhat": r, "ess_bulk": eb, "ess_tail": et, "tau_int_iter": tau,
                     "mean": float(ch.mean()), "mcse": mc}
        flag = "  OK" if (r < 1.01 and eb > 400 and et > 400) else "  CHECK"
        rows.append([name, f"{r:.4f}", f"{eb:.0f}", f"{et:.0f}", f"{tau:.1f}", f"{ch.mean():.5g}", f"{mc:.2g}", flag])
    summary["pt_diagnostics"] = out
    fig, ax = page(1, 1, "PT cold-chain convergence diagnostics (rank-normalized split R-hat, bulk/tail ESS)", h=5)
    ax = ax[0, 0]
    ax.axis("off")
    tb = ax.table(cellText=rows, colLabels=["quantity", "R-hat", "ESS bulk", "ESS tail", "tau_int [iter]",
                                             "mean", "MCSE", ""], loc="center", cellLoc="right")
    tb.auto_set_font_size(False)
    tb.set_fontsize(8.5)
    tb.scale(1, 1.4)
    for (i, j), c in tb.get_celld().items():
        c.set_edgecolor(GRID)
        if i == 0:
            c.set_text_props(color=INK, weight="bold")
    note(fig, f"Criteria (Vehtari et al. 2021): R-hat < 1.01 and ESS > 400. tau_int is in sampler iterations "
              f"(stored draws are thinned by {thin}). IO and upper octant are 0/1 indicators: they test mode mixing.")
    return fig


def page_acf(cold, thin):
    Q = quantities(cold)
    fig, axs = page(2, 5, f"Cold-chain autocorrelation (lag in stored draws; 1 draw = {thin} iterations)", h=5.5)
    for (name, x), ax in zip(Q.items(), axs.ravel()):
        ch = by_replica(cold, x)
        for r, c in enumerate(ch):
            ax.plot(D.autocorr(c, 200), lw=1.0, color=REPLICA_COLORS[r % 8])
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.set_title(name)
        ax.set_ylim(-0.1, 1)
        ax.set_xlabel("lag")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def page_pt_internals(meta, ladder, labels, rungs):
    reps = meta["replicas"]
    K = int(meta["nRungs"])
    fig, axs = page(2, 4, "Parallel-tempering internals")
    ax = axs[0, 0]
    for r, R in enumerate(reps):
        ax.plot(np.arange(K - 1) + 0.5, R["swapAccept"], "o-", ms=4, color=REPLICA_COLORS[r % 8], lw=1, label=f"replica {r}")
    ax.set_ylim(0, 1)
    ax.set_xlabel("rung pair (k, k+1)")
    ax.set_title("swap acceptance (post burn-in)")
    ax.legend()
    ax = axs[0, 1]
    for r, R in enumerate(reps):
        ax.semilogy(range(K), R["betas"], "o", ms=4, color=REPLICA_COLORS[r % 8])
    ax.set_xlabel("rung")
    ax.set_title(r"adapted ladder $\beta_k$ (frozen after burn-in)")
    ax = axs[0, 2]
    if ladder is not None:
        sel = ladder[:, 0] == 0
        for k in range(K):
            ax.semilogy(ladder[sel, 1], ladder[sel, 2 + k], color=DIRECT, lw=1)
    ax.set_xlabel("burn-in iteration")
    ax.set_title(r"ladder adaptation, replica 0")
    ax = axs[0, 3]
    for r, R in enumerate(reps):
        ax.plot(range(K), R["rwAccept"], "o-", ms=4, lw=1, color=REPLICA_COLORS[r % 8])
    ax.axhline(0.234, color=MUTED, ls="--", lw=1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("rung")
    ax.set_title("random-walk acceptance (target 0.234)")
    ax = axs[1, 0]
    jumps = meta["jumps"]
    ja = np.mean([[j["accept"] for j in R["jumpAccept"]] for R in reps], axis=0)
    if ja.size and ja.ndim == 2 and ja.shape[1] == len(jumps):
        for j, name in enumerate(jumps):
            ax.plot(range(K), ja[:, j], "o-", ms=4, lw=1, color=REPLICA_COLORS[j], label=name)
        ax.legend()
    ax.set_ylim(0, 1)
    ax.set_xlabel("rung")
    ax.set_title("mode-jump acceptance (replica mean)")
    ax = axs[1, 1]
    sel = labels[:, 0] == 0
    lab = labels[sel, 2:]
    n = min(300, len(lab))
    for l in range(min(3, K)):
        ax.plot(labels[sel, 1][:n], lab[:n, l], lw=0.8, color=REPLICA_COLORS[l])
    ax.set_xlabel("iteration")
    ax.set_ylabel("rung")
    ax.set_title("rung of three state copies, replica 0")
    ax = axs[1, 2]
    rl = rungs[rungs[:, 0] == 0, 2:]
    bins = np.linspace(0, np.quantile(-2 * rl[:, -1], 0.99), 80)
    cm = plt.get_cmap(CMAP)
    for k in range(K):
        ax.hist(-2 * rl[:, k], bins=bins, histtype="step", density=True, color=cm(0.35 + 0.65 * (K - 1 - k) / max(K - 1, 1)), lw=1.2)
    ax.set_xlabel(r"$\Delta\chi^2$ of the state on each rung")
    ax.set_title("per-rung $\\Delta\\chi^2$ (dark = cold); must overlap")
    ax = axs[1, 3]
    ax.bar(range(len(reps)), [R["roundTrips"] for R in reps], color=DIRECT, width=0.6)
    ax.set_xlabel("replica")
    ax.set_title("round trips cold -> hot -> cold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def page_running_fractions(pt_runs, exact):
    fig, axs = page(len(pt_runs), 2, "Mode mixing: running P(IO) and P(upper octant) per replica vs the exact value",
                    h=3.4 * len(pt_runs) + 1)
    for row, (label, cold) in enumerate(pt_runs):
        for col, (qname, f) in enumerate([("P(IO)", lambda c: c[:, 6] < 0), ("P(upper octant)", lambda c: c[:, 4] > 0.5)]):
            ax = axs[row, col]
            reps = np.unique(cold[:, 0]).astype(int)
            for r in reps:
                c = cold[cold[:, 0] == r]
                x = f(c).astype(float)
                ax.plot(c[:, 1], np.cumsum(x) / np.arange(1, len(x) + 1), lw=1.2, color=REPLICA_COLORS[r % 8],
                        label=f"replica {r}: {x.mean():.3f}")
            ax.axhline(exact[qname], color=INK, lw=1.2, ls="--", label=f"exact {exact[qname]:.4f}")
            ax.set_ylim(0, 1)
            ax.set_xlabel("iteration")
            ax.set_title(f"{label}: {qname}")
            ax.legend(fontsize=7, ncol=2)
    note(fig, "Without the mode-jump moves the support is disconnected in Delta m^2_3l: swaps only permute the K states, "
              "so the number of states in each ordering is fixed at initialization and the chain is reducible.")
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    return fig


def marginal_specs(quad):
    """(label, table var, ordering filter, nodes, probability per cell) for the 1D comparisons."""
    pIO = quad["pIO"]
    specs = []
    for v in ["s12", "s13", "s23", "dcp", "l21"]:
        nN, mN = quad["cells"]["NO"][v]
        nI, mI = quad["cells"]["IO"][v]
        # the populated (trimmed) node ranges can differ between orderings: combine on the union
        if np.array_equal(nN, nI):
            specs.append((v, v, None, nN, (1 - pIO) * mN + pIO * mI))
        else:
            allx = np.union1d(nN, nI)
            m = np.zeros(len(allx) - 1)
            for nodes, mass, w in ((nN, mN, 1 - pIO), (nI, mI, pIO)):
                i0 = np.searchsorted(allx, nodes[0])
                m[i0:i0 + len(mass)] += w * mass
            specs.append((v, v, None, allx, m))
    for o, w in (("NO", 1 - pIO), ("IO", pIO)):
        n, m = quad["cells"][o]["dma"]
        specs.append((f"dma ({o})", "dma", o, n, w * m))
    return specs


def page_marginals(quad, direct, pt_cold, pt_ess, summary):
    specs = marginal_specs(quad)
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("1D marginals: quadrature of the surrogate vs exact sampler vs PT", fontsize=12, x=0.01, ha="left")
    gs = fig.add_gridspec(4, 4, height_ratios=[3, 1, 3, 1], hspace=0.08, wspace=0.3)
    Xd = samples_to_table_coords(direct[:, :6])
    Xp = samples_to_table_coords(pt_cold[:, 2:8]) if pt_cold is not None else None
    res = {}
    for k, (label, v, o, nodes, mass) in enumerate(specs):
        r, c = 2 * (k // 4), k % 4
        ax, axp = fig.add_subplot(gs[r, c]), fig.add_subplot(gs[r + 1, c])
        edges, prob = merged_edges(nodes, mass)
        w = np.diff(edges)
        x = 0.5 * (edges[1:] + edges[:-1])
        ax.stairs(prob / w, edges, color=INK, lw=1.2, label="quadrature")
        out = {}
        for name, X, colr, sel, neff in (
                ("direct", Xd, DIRECT, None if o is None else (direct[:, 4] > 0) == (o == "NO"), None),
                ("PT", Xp, PT, None if (o is None or Xp is None) else (pt_cold[:, 6] > 0) == (o == "NO"),
                 pt_ess.get(v))):
            if X is None:
                continue
            vals = X[v] if sel is None else X[v][sel]
            N = len(X[v])
            h = np.histogram(vals, bins=edges)[0]
            ax.stairs(h / N / w, edges, color=colr, lw=1.0, label=name)
            stat, ndf, p = binned_chi2(h, prob, N, n_eff=neff)
            out[name] = {"chi2": stat, "ndf": ndf, "p": p}
            sig = np.sqrt(prob * N * (N / neff if neff else 1.0))
            ok = prob * N > 20
            axp.plot(x[ok], ((h - prob * N) / np.where(sig > 0, sig, 1))[ok], "o", ms=2.5, color=colr)
        res[label] = out
        ax.set_title(label + "  " + "  ".join(f"{n}: p={d['p']:.2f}" for n, d in out.items()), fontsize=8)
        ax.tick_params(labelbottom=False)
        axp.axhline(0, color=MUTED, lw=0.8)
        axp.set_ylim(-4.5, 4.5)
        axp.set_xlabel(LABELS[v], fontsize=8)
        axp.set_ylabel("pull", fontsize=7)
        if k == 0:
            ax.legend()
    summary["marginal_chi2"] = res
    note(fig, "Pulls use Poisson errors; for PT they are inflated by sqrt(N/ESS_bulk) of that parameter. p-values: "
              "Pearson chi^2 on bins with > 20 expected entries (PT statistic scaled by ESS/N).")
    return fig


def page_contours(S, direct, title):
    pairs = [("s23", "dcp", None), ("s13", "dcp", None), ("s13", "s23", None), ("s12", "l21", None),
             ("s23", "dma", "NO"), ("s23", "dma", "IO"), ("s13", "dma", "NO"), ("dma", "dcp", "NO")]
    fig, axs = page(2, 4, title)
    X = samples_to_table_coords(direct[:, :6])
    for (a, b, o), ax in zip(pairs, axs.ravel()):
        tabs = S.tables["NO"]
        tag = next(t for t in tabs if t.count("/") == 1 and set(tabs[t].vars) == {a, b})
        def nufit(oo):
            t = S.tables[oo][tag]
            v = t.values if t.vars == (a, b) else t.values.T
            return t.axes[t.vars.index(a)], t.axes[t.vars.index(b)], v
        if o is None:
            xa, xb, vN = nufit("NO")
            _, _, vI = nufit("IO")
            V = np.minimum(vN, vI)
            prof = S.profile([a, b])
            sel = np.ones(len(direct), bool)
        else:
            xa, xb, V = nufit(o)
            prof = S.profile([a, b], o)
            sel = (direct[:, 4] > 0) == (o == "NO")
        box = region_box(V, [xa, xb], 16.0)
        ax.contour(xa, xb, V.T, levels=LEVELS_2D, colors=INK, linewidths=[1.8, 1.2, 0.8])
        ax.contour(xa, xb, prof.T, levels=LEVELS_2D, colors=DIRECT, linewidths=[1.8, 1.2, 0.8], linestyles="--")
        H, ea, eb = np.histogram2d(X[a][sel], X[b][sel], bins=60, range=box)
        H = gaussian_filter(H, 1.0)
        lv = sorted(set(hpd_levels(H, MASSES)))
        if len(lv) == 3:
            ax.contour(0.5 * (ea[1:] + ea[:-1]), 0.5 * (eb[1:] + eb[:-1]), H.T, levels=lv, colors=PT,
                       linewidths=[0.8, 1.2, 1.8], linestyles=":")
        ax.set_xlim(box[0])
        ax.set_ylim(box[1])
        ax.set_xlabel(LABELS[a])
        ax.set_ylabel(LABELS[b])
        ax.set_title(f"{tag}" + (f" ({o})" if o else " (min over orderings)"), fontsize=9)
    h = [plt.Line2D([], [], color=INK, lw=1.5, label=r"NuFIT 6.0 table, $\Delta\chi^2$ = 2.30/6.18/11.83"),
         plt.Line2D([], [], color=DIRECT, lw=1.5, ls="--", label="surrogate profile, same levels"),
         plt.Line2D([], [], color=PT, lw=1.5, ls=":", label="sampled posterior HPD 68.3/95.4/99.7%")]
    fig.legend(handles=h, loc="upper right", ncol=3, bbox_to_anchor=(0.99, 0.965))
    note(fig, "Profile contours (frequentist, NuFIT) and posterior credible regions (Bayesian, flat priors, marginalized) "
              "are different objects and need not coincide; the surrogate profile must match NuFIT exactly for "
              "projections of a single factor. Per-ordering panels: posterior restricted to that ordering.")
    fig.tight_layout(rect=(0, 0.03, 1, 0.93))
    return fig


def page_residuals(S, o, summary):
    exact_tags = {"T23/DMA/DCP", "T12/DMS", "T23/DMA", "T23/DCP", "DMA/DCP"}
    if S.coupled:
        exact_tags.add("T13/T23")
    tabs = S.tables[o]
    tags = [t for t in tabs if t.count("/") == 1]
    fig, axs = page(4, 4, f"Surrogate profile minus NuFIT 6.0 table, {o} (all 15 two-dimensional projections)")
    res = {}
    m = S.m[o]
    for tag, ax in zip(tags, axs.ravel()):
        t = tabs[tag]
        prof = S.profile(list(t.vars), o)
        d = prof - t.values
        inside = t.values < m + LEVELS_2D[-1]
        mx = float(np.abs(d[inside]).max())
        res[tag] = {"max_abs_in_3sigma": mx, "rms_in_3sigma": float(np.sqrt(np.mean(d[inside] ** 2))),
                    "by_construction": tag in exact_tags}
        box = region_box(t.values - m, t.axes, 25.0, pad=0.1)
        dm = np.ma.masked_where(t.values > m + 25, d)
        im = ax.pcolormesh(t.axes[0], t.axes[1], dm.T, cmap=DIVERGING, vmin=-3, vmax=3, shading="nearest")
        ax.contour(t.axes[0], t.axes[1], (t.values - m).T, levels=LEVELS_2D, colors=INK, linewidths=0.7)
        ax.set_xlim(box[0])
        ax.set_ylim(box[1])
        ax.set_xlabel(LABELS[t.vars[0]], fontsize=7)
        ax.set_ylabel(LABELS[t.vars[1]], fontsize=7)
        kind = "exact by construction" if tag in exact_tags else "held out"
        ax.set_title(f"{tag} [{kind}]\nmax|res| in 3$\\sigma$: {mx:.2f}", fontsize=8)
    axs.ravel()[-1].axis("off")
    fig.colorbar(im, ax=axs.ravel()[-1], fraction=0.8, label=r"$\Delta\chi^2_{surrogate} - \Delta\chi^2_{NuFIT}$")
    note(fig, "Contours: NuFIT Delta chi^2 (relative to this ordering's minimum) = 2.30/6.18/11.83. Negative residual = "
              "surrogate too permissive there. Residuals on held-out projections measure what the 3D x 2D x 1D "
              "factorization loses (chiefly theta13 correlations).")
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    summary[f"residuals_{o}"] = res
    return fig


def page_1d(S, direct):
    fig, axs = page(2, 4, r"1D: NuFIT $\Delta\chi^2$ vs surrogate profile vs $-2\ln$ (marginal posterior / max)")
    X = samples_to_table_coords(direct[:, :6])
    tagof = {"s12": "T12", "s13": "T13", "s23": "T23", "dcp": "DCP", "l21": "DMS", "dma": "DMA"}
    panels = [("s12", None), ("s13", None), ("s23", None), ("dcp", None), ("l21", None), ("dma", "NO"), ("dma", "IO")]
    for (v, o), ax in zip(panels, axs.ravel()):
        if o is None:
            t = S.tables["NO"][tagof[v]]
            V = np.minimum(t.values, S.tables["IO"][tagof[v]].values)
            prof = S.profile([v])
            vals = X[v]
            x = t.axes[0]
        else:
            t = S.tables[o][tagof[v]]
            V = t.values
            prof = S.profile([v], o)
            vals = X[v][(direct[:, 4] > 0) == (o == "NO")]
            x = t.axes[0]
        sel = V < 16
        lo, hi = x[sel].min(), x[sel].max()
        ax.plot(x, V, color=INK, lw=2.2, label="NuFIT table")
        ax.plot(x, prof, color=DIRECT, lw=1.2, ls="--", label="surrogate profile")
        h, e = np.histogram(vals, bins=80, range=(lo, hi), density=True)
        ok = h > 0
        y = -2 * np.log(h[ok] / h.max()) + (S.m[o] if o else 0)
        ax.plot(0.5 * (e[1:] + e[:-1])[ok], y, color=PT, lw=1.2, label=r"$-2\ln p/p_{max}$ (samples)")
        ax.set_xlim(lo, hi)
        ax.set_ylim(0, 16)
        ax.set_xlabel(LABELS[v] + (f" ({o})" if o else ""))
        ax.set_ylabel(r"$\Delta\chi^2$")
    axs[0, 0].legend()
    axs.ravel()[-1].axis("off")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def page_overlay(title, datasets, notes=None):
    """1D density overlays of several sample sets: datasets = [(label, P(N,6), colour)]."""
    fig, axs = page(2, 4, title)
    panels = [("s12", None), ("s13", None), ("s23", None), ("dcp", None), ("l21", None), ("dma", "NO"), ("dma", "IO")]
    rows = []
    for (v, o), ax in zip(panels, axs.ravel()):
        allv = []
        for _, P, _ in datasets:
            X = samples_to_table_coords(P)
            allv.append(X[v] if o is None else X[v][(P[:, 4] > 0) == (o == "NO")])
        ref = allv[0]
        lo, hi = np.quantile(ref, [0.0005, 0.9995]) if len(ref) else (0, 1)
        for (label, P, colr), vals in zip(datasets, allv):
            w = np.full(len(vals), 1.0 / len(P))
            ax.hist(vals, bins=70, range=(lo, hi), weights=w * 70 / (hi - lo), histtype="step", color=colr, lw=1.3,
                    label=label)
        ax.set_xlabel(LABELS[v] + (f" ({o})" if o else ""))
        ax.set_yticks([])
    for label, P, _ in datasets:
        rows.append(f"{label}: P(IO)={np.mean(P[:, 4] < 0):.3f}  P(UO)={np.mean(P[:, 2] > 0.5):.3f}  "
                    f"sd(s13)={P[:, 1].std():.5f}  sd(s23)={P[:, 2].std():.4f}  sd(dm21)={P[:, 5].std():.3f}")
    ax = axs.ravel()[-1]
    ax.axis("off")
    axs[0, 0].legend(fontsize=7)
    ax.text(0, 1, "\n\n".join(rows), va="top", fontsize=7, color=INK2, wrap=True)
    if notes:
        note(fig, notes)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    return fig


# --------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", default="TBoff", choices=["TBoff", "TByes"])
    ap.add_argument("--surrogate", default="factorized", choices=["factorized", "coupled"])
    ap.add_argument("--results", default="../results")
    ap.add_argument("--data", default="../data/nufit")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    tag = f"{args.variant}_{args.surrogate}"
    out = args.out or os.path.join(args.results, f"report_{tag}")
    os.makedirs(out, exist_ok=True)
    R = lambda name: os.path.join(args.results, name)  # noqa: E731

    S = Surrogate(args.data, args.variant, coupled=args.surrogate == "coupled")
    direct = np.load(R(f"direct_{tag}_samples.npy"))
    pt_cold = load_optional(R(f"pt_{tag}_cold.npy"))
    meta = json.load(open(R(f"pt_{tag}_meta.json"))) if os.path.exists(R(f"pt_{tag}_meta.json")) else None
    nojumps = load_optional(R(f"pt_{tag}_nojumps_cold.npy"))

    summary = {"model": f"NuFIT 6.0 {args.variant}, {args.surrogate} surrogate, flat priors (sin^2, dCP, Dm3l, Dm21)",
               "n_direct": len(direct)}
    # implementation cross-check: Python mirror vs C++ log-likelihood
    sub = direct[:: max(1, len(direct) // 200000)]
    summary["max_abs_logL_py_minus_cpp"] = float(np.nanmax(np.abs(-0.5 * S.delta_chi2(sub[:, :6]) - sub[:, 6])))
    quad = S.marginals(r=4)
    exact = {"P(IO)": quad["pIO"], "P(upper octant)": float(np.mean(direct[:, 2] > 0.5))}
    summary["P(IO) quadrature"] = quad["pIO"]
    summary["P(IO) direct"] = float(np.mean(direct[:, 4] < 0))
    summary["P(IO) direct binomial error"] = float(np.sqrt(quad["pIO"] * (1 - quad["pIO"]) / len(direct)))
    summary["P(upper octant) direct"] = exact["P(upper octant)"]
    best = direct[np.argmax(direct[:, 6])]
    summary["best sampled point (s12,s13,s23,dcp,dm3l,dm21)"] = [float(x) for x in best[:6]]

    pt_ess = {}
    thin = int(meta["options"]["thin"]) if meta else 1
    if pt_cold is not None:
        X = samples_to_table_coords(pt_cold[:, 2:8])
        for v in VARS:
            pt_ess[v] = D.ess_bulk(by_replica(pt_cold, X[v] if v != "dcp" else np.cos(np.radians(X[v]))))
        io = by_replica(pt_cold, (pt_cold[:, 6] < 0).astype(float))
        summary["P(IO) PT"] = float(io.mean())
        summary["P(IO) PT MCSE"] = float(D.mcse_mean(io))
        summary["P(IO) PT pull vs exact"] = float((io.mean() - quad["pIO"]) / D.mcse_mean(io))
        if nojumps is not None:
            ioj = by_replica(nojumps, (nojumps[:, 6] < 0).astype(float))
            summary["P(IO) PT without jumps, per replica"] = [float(x) for x in ioj.mean(axis=1)]
            summary["R-hat(IO) PT without jumps"] = float(D.rhat(ioj))

    figs = []
    pdf_path = os.path.join(out, "report.pdf")
    def emit(name, fig):
        figs.append((name, fig))
        print("page", name)

    emit("interp_loo", page_interp_loo(S, summary))
    emit("slices", page_slices(S, best[:6]))
    if pt_cold is not None:
        emit("traces", page_traces(pt_cold, "PT cold-chain traces (all replicas, post burn-in)"))
        emit("diagnostics", page_diag_table(pt_cold, summary, thin))
        emit("acf", page_acf(pt_cold, thin))
        emit("pt_internals", page_pt_internals(meta, load_optional(R(f"pt_{tag}_ladder.npy")),
                                               np.load(R(f"pt_{tag}_labels.npy")), np.load(R(f"pt_{tag}_rungs.npy"))))
        runs = [("PT with mode jumps", pt_cold)] + ([("PT without mode jumps", nojumps)] if nojumps is not None else [])
        emit("mode_mixing", page_running_fractions(runs, exact))
    emit("marginals", page_marginals(quad, direct, pt_cold, pt_ess, summary))
    emit("contours", page_contours(S, direct, "Key 2D projections: NuFIT 6.0 vs surrogate vs sampled posterior"))
    emit("residuals_NO", page_residuals(S, "NO", summary))
    emit("residuals_IO", page_residuals(S, "IO", summary))
    emit("one_d", page_1d(S, direct))

    legacy = load_optional(R("legacy_samples.npy"))
    if legacy is not None and args.variant == "TBoff":
        emit("legacy", page_overlay("Student's sampler (legacy) vs the exact posterior (TBoff)",
                                    [("exact (direct)", direct[:, :6], DIRECT), ("legacy chain", legacy, LEGACY)],
                                    "The legacy chain used exp(-chi^2) instead of exp(-chi^2/2) and a biased chain-selection "
                                    "swap: its marginals are too narrow and its ordering weights are wrong."))
    sens = [("default priors", direct[:, :6], DIRECT)]
    angle = load_optional(R(f"direct_{args.variant}_factorized_angleprior_samples.npy"))
    if angle is not None:
        sens.append(("flat in angles, log Dm21", angle[:, :6], PT))
    other = "coupled" if args.surrogate == "factorized" else "factorized"
    oth = load_optional(R(f"direct_{args.variant}_{other}_samples.npy"))
    if oth is not None:
        sens.append((f"{other} surrogate", oth[:, :6], LEGACY))
    ov = "TByes" if args.variant == "TBoff" else "TBoff"
    ovs = load_optional(R(f"direct_{ov}_{args.surrogate}_samples.npy"))
    if ovs is not None:
        sens.append((f"{ov} (other NuFIT variant)", ovs[:, :6], EXTRA))
    if len(sens) > 1:
        emit("sensitivity", page_overlay("Sensitivity of the posterior to priors, surrogate and NuFIT variant", sens,
                                         "These are the inputs whose choice should be propagated into the "
                                         "consistency p-value as a systematic (report the spread)."))
    figs.insert(0, ("summary", page_summary(summary, [n for n, _ in figs])))
    for old in os.listdir(out):
        if old.endswith(".png"):
            os.remove(os.path.join(out, old))
    with PdfPages(pdf_path) as pdf:
        for i, (name, fig) in enumerate(figs):
            fig.savefig(os.path.join(out, f"{i:02d}_{name}.png"), dpi=110)
            pdf.savefig(fig)
            plt.close(fig)

    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=float)
    print(f"report -> {pdf_path}")


if __name__ == "__main__":
    main()

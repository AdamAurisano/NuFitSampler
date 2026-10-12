"""Reader for the NuFIT 6.0 Delta chi^2 release files.

Each section ("# T23/DMA/DCP projection: ...") is turned into a dense array on its
rectilinear grid. Parsed files are cached as .npz next to the originals.
"""

import os
import re

import numpy as np

# Short variable names used throughout the validation code (table coordinates).
#   s12, s13, s23: sin^2(theta_ij); dcp [deg]; dma: Delta m^2_3l [1e-3 eV^2];
#   l21: log10(Delta m^2_21 / eV^2)
TAG_VARS = {"T12": "s12", "T13": "s13", "T23": "s23", "DCP": "dcp", "DMA": "dma", "DMS": "l21"}
VARS = ["s12", "s13", "s23", "dcp", "dma", "l21"]
LABELS = {
    "s12": r"$\sin^2\theta_{12}$",
    "s13": r"$\sin^2\theta_{13}$",
    "s23": r"$\sin^2\theta_{23}$",
    "dcp": r"$\delta_{CP}$ [deg]",
    "dma": r"$\Delta m^2_{3\ell}$ [$10^{-3}$ eV$^2$]",
    "l21": r"$\log_{10}(\Delta m^2_{21}/\mathrm{eV}^2)$",
}


class Table:
    """Delta chi^2 on a complete rectilinear grid; values[i, j, ...] at (axes[0][i], axes[1][j], ...)."""

    def __init__(self, tag, vars_, axes, values):
        self.tag, self.vars, self.axes, self.values = tag, tuple(vars_), list(axes), values

    def __repr__(self):
        return f"Table({self.tag}, {self.vars}, shape={self.values.shape})"


def _to_table(tag, rows):
    vars_ = tuple(TAG_VARS[t] for t in tag.split("/"))
    nd = len(vars_)
    rows = rows.reshape(-1, nd + 1)
    axes = [np.unique(rows[:, d]) for d in range(nd)]
    shape = tuple(len(a) for a in axes)
    if np.prod(shape) != len(rows):
        raise ValueError(f"{tag}: not a complete grid")
    idx = tuple(np.searchsorted(axes[d], rows[:, d]) for d in range(nd))
    values = np.full(shape, np.nan)
    values[idx] = rows[:, nd]
    if np.isnan(values).any():
        raise ValueError(f"{tag}: duplicate nodes")
    return Table(tag, vars_, axes, values)


def read_release(path):
    """Return {tag: Table} for one NuFIT release file (uses/creates a .npz cache)."""
    cache = path + ".cache.npz"
    if os.path.exists(cache) and os.path.getmtime(cache) >= os.path.getmtime(path):
        z = np.load(cache, allow_pickle=False)
        out = {}
        for tag in z["tags"]:
            nd = len(tag.split("/"))
            axes = [z[f"{tag}|ax{d}"] for d in range(nd)]
            out[str(tag)] = Table(str(tag), [TAG_VARS[t] for t in tag.split("/")], axes, z[f"{tag}|val"])
        return out
    text = open(path).read()
    heads = list(re.finditer(r"^#\s*(\S+) projection:.*$", text, flags=re.M))
    out = {}
    for i, m in enumerate(heads):
        end = heads[i + 1].start() if i + 1 < len(heads) else len(text)
        body = text[m.end():end]
        rows = np.array(body.split(), dtype=float)
        out[m.group(1)] = _to_table(m.group(1), rows)
    arrays = {"tags": np.array(list(out))}
    for tag, t in out.items():
        for d, a in enumerate(t.axes):
            arrays[f"{tag}|ax{d}"] = a
        arrays[f"{tag}|val"] = t.values
    np.savez(cache, **arrays)
    return out


def release_path(datadir, variant, ordering):
    return os.path.join(datadir, f"v60.release-{variant}-{ordering}.txt")


def load_variant(datadir, variant):
    """{'NO': {tag: Table}, 'IO': {tag: Table}} for variant 'TBoff' or 'TByes'."""
    return {o: read_release(release_path(datadir, variant, o)) for o in ("NO", "IO")}

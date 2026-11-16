"""Run metadata, wall extraction (where the configurations are) and loading of the stored walls.

Stored files (tracked by git, the configurations are not):
  flat walls  <run>/walls.npz      one period of the wall of every critical configuration
                                   (open curve, x from 0 to L) + the KFR column profile u(x)
  bubbles     <run>/contours.npz   every phi = 0 contour of every critical configuration
                                   (the bubble wall is the largest positive-area one)
Curves are stored as a start vertex + int8 steps at 1/8 site.
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
SCALE = 8
CFG_RE = re.compile(r"L(\d+)_a([0-9.]+)_s(\d+)\.npz$")


# ---------------------------------------------------------------------------
# metadata
# ---------------------------------------------------------------------------
def load_results(d):
    """{(L, a, seed): record} from all results*.jsonl shards of a run."""
    meta = {}
    for fn in sorted(glob.glob(os.path.join(d, "results*.jsonl"))):
        for line in open(fn):
            r = json.loads(line)
            meta[(int(r["L"]), float(r["aniso"]), int(r["seed"]))] = r
    return meta


def run_label(d, meta):
    """Short description of a run from its metadata."""
    r = next(iter(meta.values()))
    b = float(r.get("burgers_angle", 0.0))
    if r.get("geometry") == "flat":
        kind = {0.0: "screw walls (line ∥ b)", 90.0: "edge walls (line ⊥ b)"}.get(b, f"walls, b at {b:g}°")
        return f"{kind}, Ly/L = {r.get('Ly', 2 * r['L']) / r['L']:g}"
    R0 = {round(m["R0"] / m["L"], 4) for m in meta.values() if m.get("R0")}
    r0 = {m["R0"] for m in meta.values() if m.get("R0")}
    size = f"R0 = {r0.pop():g}" if len(r0) == 1 else (f"R0 = L/{1 / R0.pop():.0f}" if len(R0) == 1 else "R0 varies")
    return f"bubbles, {size}" + (f", b at {b:g}°" if b else "")


def default_out(dirs):
    """runs/<common name of the runs>_analysis, e.g. flat_D0.04_b0_A2 + flat_D0.04_b90_A2 ->
    flat_D0.04_A2_analysis, fss_D0.04 -> fss_D0.04_analysis."""
    names = [os.path.basename(os.path.normpath(d)) for d in dirs]
    pre = os.path.commonprefix(names)
    suf = os.path.commonprefix([n[::-1] for n in names])[::-1] if len(names) > 1 else ""
    if len(names) > 1:          # cut both at a "_" unless they end on a whole token already
        if not all(n == pre or n[len(pre):].startswith("_") for n in names):
            pre = pre[:pre.rfind("_") + 1] if "_" in pre else ""
        if not all(n == suf or n[:len(n) - len(suf)].endswith("_") for n in names):
            suf = suf[suf.find("_"):] if "_" in suf else ""
    name = (pre + suf.lstrip("_") if pre.endswith("_") else pre + suf).strip("_") or "analysis"
    return os.path.join(os.path.dirname(os.path.normpath(dirs[0])), name + "_analysis")


def _files(d):
    files = sorted(glob.glob(os.path.join(d, "configs", "*.npz")))
    if not files:
        sys.exit(f"no configurations in {d}/configs: run the extraction where the scan ran")
    return files


def _encode(curves):
    """list of (n, 2) float arrays -> start (m, 2) int32, steps (N, 2) int8, offsets (m + 1)."""
    q = [np.rint(np.asarray(c) * SCALE).astype(np.int64) for c in curves]
    steps = [np.diff(c, axis=0) for c in q]
    st = np.concatenate(steps) if steps else np.zeros((0, 2), np.int64)
    assert np.abs(st).max(initial=0) < 128, "contour step too long for int8"
    return (np.array([c[0] for c in q], np.int32).reshape(-1, 2), st.astype(np.int8),
            np.r_[0, np.cumsum([len(s) for s in steps])].astype(np.int64))


def _decoder(z):
    """Curve i of a walls / contours file (arrays read once: npz members are re-read per access)."""
    start, st, off, sc = z["start"].astype(np.int64), z["steps"].astype(np.int64), z["offsets"], float(z["scale"])

    def curve(i):
        return np.vstack([start[i], start[i] + np.cumsum(st[off[i]:off[i + 1]], axis=0)]) / sc
    return curve


# ---------------------------------------------------------------------------
# flat walls
# ---------------------------------------------------------------------------
def flat_wall_curve(phi):
    """One period of the wall (top boundary of the phi > 0 domain attached to the bottom,
    islands ignored, holes left out) as an (n, 2) array (x, y) from x = 0 to x = L:
    the phi = 0 contour of the field tiled 3 times in x."""
    import contourpy
    import bubble_depinning as bd
    L = phi.shape[1]
    m = bd.wall_mask(phi)
    f = np.where(m, phi, -np.abs(phi))
    lines = contourpy.contour_generator(z=np.hstack([f, f, f])).lines(0.0)
    if not lines:
        return None
    c = max(lines, key=lambda a: np.ptp(a[:, 0]))
    if np.ptp(c[:, 0]) < 2 * L:
        return None
    if c[0, 0] > c[-1, 0]:
        c = c[::-1]
    i1 = int(np.argmax(c[:, 0] >= L))
    target = c[i1] + np.array([L, 0.0])
    rest = np.arange(i1 + 1, len(c))
    i2 = rest[np.argmin(np.hypot(*(c[rest] - target).T))]
    return c[i1:i2] - np.array([L, 0.0])


def _flat_one(fn):
    import bubble_depinning as bd
    L, a, seed = CFG_RE.search(fn).groups()
    z = np.load(fn)
    phi = z["phi_c"].astype(np.float64)
    seg = flat_wall_curve(phi)
    an = bd.analyse_flat(phi)
    if seg is None or an is None:
        return None
    return dict(L=int(L), aniso=float(a), seed=int(seed), h_d=float(z["h_d"]), curve=seg,
                u=an["u"].astype(np.float32), mult=np.minimum(an["mult"], 127).astype(np.int8))


def extract_flat(d, workers):
    files = _files(d)
    with ProcessPoolExecutor(workers) as ex:
        rows = [r for r in ex.map(_flat_one, files, chunksize=4) if r is not None]
    rows.sort(key=lambda r: (r["aniso"], r["L"], r["seed"]))
    out = {k: np.array([r[k] for r in rows]) for k in ("L", "aniso", "seed", "h_d")}
    out["start"], out["steps"], out["offsets"] = _encode([r["curve"] for r in rows])
    for L in sorted({r["L"] for r in rows}):
        out[f"u_L{L}"] = np.stack([r["u"] for r in rows if r["L"] == L])
        out[f"mult_L{L}"] = np.stack([r["mult"] for r in rows if r["L"] == L])
    out["n_configs"] = len(files)
    out["scale"] = SCALE
    fn = os.path.join(d, "walls.npz")
    np.savez_compressed(fn, **out)
    print(f"saved {fn} ({os.path.getsize(fn) / 1e6:.1f} MB, {len(rows)} of {len(files)} walls)")


def load_flat(d):
    """Walls of a flat run: dict of per-wall arrays + lists 'curve', 'u', 'mult'."""
    z = np.load(os.path.join(d, "walls.npz"))
    W = {k: z[k] for k in ("L", "aniso", "seed", "h_d")}
    n = len(W["L"])
    dec = _decoder(z)
    W["curve"] = [dec(i) for i in range(n)]
    W["u"], W["mult"] = [None] * n, [None] * n
    for L in np.unique(W["L"]):
        for i, u, m in zip(np.nonzero(W["L"] == L)[0], z[f"u_L{L}"], z[f"mult_L{L}"]):
            W["u"][i], W["mult"][i] = u, m
    meta = load_results(d)
    rec = [meta.get((int(L), float(a), int(s)), {}) for L, a, s in zip(W["L"], W["aniso"], W["seed"])]
    W["travel"] = np.array([r.get("travel", np.nan) for r in rec], float)
    W["Ly"] = np.array([r.get("Ly", 2 * L) for r, L in zip(rec, W["L"])], float)
    W["burgers_angle"] = float(next(iter(meta.values())).get("burgers_angle", 0.0)) if meta else 0.0
    W["label"] = run_label(d, meta) if meta else os.path.basename(d)
    W["geometry"] = "flat"
    return W


# ---------------------------------------------------------------------------
# bubbles
# ---------------------------------------------------------------------------
def _bubble_one(fn):
    """All phi = 0 contours of one configuration, with the signed area (> 0: encloses phi > 0)."""
    import contourpy
    L, a, seed = CFG_RE.search(fn).groups()
    z = np.load(fn)
    phi = z["phi_c"].astype(np.float64)
    out = []
    for c in contourpy.contour_generator(z=phi).lines(0.0):
        if len(c) < 4:
            continue
        x, y = c[:, 0], c[:, 1]
        A = 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
        # which side is phi > 0? sample phi half a site to the left of a few segments
        k = np.linspace(0, len(c) - 2, min(9, len(c) - 1)).astype(int)
        t = c[k + 1] - c[k]
        nrm = np.c_[-t[:, 1], t[:, 0]] / np.maximum(np.hypot(*t.T), 1e-12)[:, None]
        q = 0.5 * (c[k] + c[k + 1]) + 0.5 * nrm
        qi = np.clip(np.rint(q).astype(int), 0, np.array(phi.shape[::-1]) - 1)
        left_pos = np.sign(np.sum(np.sign(phi[qi[:, 1], qi[:, 0]]))) or 1.0
        out.append((c, abs(A) * np.sign(A) * left_pos))
    return dict(L=int(L), aniso=float(a), seed=int(seed), h_d=float(z["h_d"]), contours=out)


def extract_bubbles(d, workers):
    files = _files(d)
    with ProcessPoolExecutor(workers) as ex:
        rows = list(ex.map(_bubble_one, files, chunksize=4))
    rows.sort(key=lambda r: (r["aniso"], r["L"], r["seed"]))
    out = {k: np.array([r[k] for r in rows]) for k in ("L", "aniso", "seed", "h_d")}
    cs = [c for r in rows for c in r["contours"]]
    out["start"], out["steps"], out["offsets"] = _encode([c[0] for c in cs])
    out["start"] = out["start"].astype(np.uint16)
    out["area"] = np.array([c[1] for c in cs], np.float32)
    out["sample"] = np.repeat(np.arange(len(rows)), [len(r["contours"]) for r in rows]).astype(np.int32)
    out["scale"] = SCALE
    out["n_configs"] = len(files)
    fn = os.path.join(d, "contours.npz")
    np.savez_compressed(fn, **out)
    print(f"saved {fn} ({os.path.getsize(fn) / 1e6:.1f} MB, {len(cs)} contours)")


def load_bubbles(d):
    """Bubbles of a run: dict of per-bubble arrays + list 'curve' (the wall, closed) and the
    number of holes (phi < 0 inside) and islands (phi > 0 domains apart from the bubble)."""
    z = np.load(os.path.join(d, "contours.npz"))
    W = {k: z[k] for k in ("L", "aniso", "seed", "h_d")}
    n = len(W["L"])
    smp, area = z["sample"], z["area"].astype(float)
    main = np.full(n, -1)
    best = np.full(n, -np.inf)
    for j, (i, A) in enumerate(zip(smp, area)):
        if A > best[i]:
            best[i], main[i] = A, j
    keep = main >= 0
    W = {k: v[keep] for k, v in W.items()}
    dec = _decoder(z)
    W["curve"] = [dec(j) for j in main[keep]]
    W["n_holes"] = np.bincount(smp[area < 0], minlength=n)[keep]
    W["n_islands"] = (np.bincount(smp[area > 0], minlength=n) - 1)[keep]
    meta = load_results(d)
    rec = [meta.get((int(L), float(a), int(s)), {}) for L, a, s in zip(W["L"], W["aniso"], W["seed"])]
    W["R0"] = np.array([r.get("R0") or np.nan for r in rec], float)
    W["burgers_angle"] = float(next(iter(meta.values())).get("burgers_angle", 0.0)) if meta else 0.0
    W["label"] = run_label(d, meta) if meta else os.path.basename(d)
    W["geometry"] = "bubble"
    return W

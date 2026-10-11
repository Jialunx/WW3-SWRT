"""Fold WRT's working set in its quadrature database quad99999.bqf into the row table that SWRT selects from, for one grid.

Xnl makes the database for the grid (xnl_ref.init) and gives its run-time constants (xnl_ref.consts).
WRT's run-time window, from Xnl's filter settings on this grid: k3/k1 <= qf_krat = 2.5 (di = jk3 - jk1 = 0..di_max) and
|theta3 - theta1| <= qf_dmax = 75 deg (amem - 1 = 0..dj_max); on the paper's grid (ratio 1.1, 15 deg) di_max = 4, dj_max = 5.
The window must be the same at every frequency (it is unless k3/k1 or the direction difference falls on a threshold).
Deep-water loci on a grid with a constant frequency ratio are self-similar, so one configuration (di, amem) has the same node
offsets at every jk1. Row = (config, node): geometry (u2, v2, u4, v4) from the reference bin, weights zz(jk1) for all nk bins,
matched by node index when the node count equals the reference count, otherwise by geometry (tolerance 0.25 bin). The
reference bin is nk // 2 (bin 20 on the paper's grid), or higher if the loci reach below the lowest frequency from there
(ref_bin). The fold stops if a node cannot be placed, which happens on grids too short to hold WRT's loci. A node whose k2 or k4 lies just beyond the highest bin at the reference bin is stored there in WRT's out-of-grid
encoding (ik = nk, fraction (k - k_nk)/50), which is not its true offset; its geometry is taken from the nearest lower jk1
where it is inside the grid (at the reference bin itself both corners clamp to bin nk either way).
Direction senses exactly as WRT visits them (q_xnl4v4 loops ik3 = ik1..nk, ia3 = 1..na):
  amem = 2..na/2: +dj and -dj (the -dj row is the mirror image, v -> -v);  amem = 1 (di > 0) and amem = na/2 + 1: once;
  di = 0 is included with both senses, i.e. every equal-frequency pair is visited from both ends, as in WRT.
Also stored: WRT's run-time constants, so the SWRT kernel never needs the WRT library.
usage: python3 build_table.py --f1 0.03453 --ratio 1.1 --nf 40 --ndir 24 --work DIR [--constants paper|xnl]
       ->  DIR/table.npz, DIR/table_build.json"""
import os, json, hashlib, argparse
import numpy as np
from bqf_decode import decode, offsets, raw_nodes, NAMES

TOL = 0.25
UNMATCHED_MAX = 1e-6        # largest share of WRT's weight that may be left unmatched (3.5e-18 on the paper's grid)
TOO_SHORT = "the grid is too short for WRT's loci (see 'Grid length' in tools/README.md): use more frequencies"


def valid_hash(Q):
    """sha256 of the slots of the database that WRT reads (those past nloc are uninitialised)"""
    m = np.arange(Q["nl"])[None, None, :] < Q["nloc"][:, :, None]; h = hashlib.sha256(np.ascontiguousarray(Q["nloc"], "<i4").tobytes())
    for nm in NAMES: h.update(np.ascontiguousarray(Q[nm][m], "<f4" if Q[nm].dtype.kind == "f" else "<i4").tobytes())
    return h.hexdigest()


def window(Q, C):
    """di_max and dj_max of WRT's run-time filter on this grid (Xnl, q_xnl4v4: skip if k3/k1 > qf_krat or the direction
    difference > qf_dmax, in single precision), checked to be the same at every frequency and direction"""
    nk, na = Q["nk"], Q["na"]; k = np.asarray(C["k"], np.float32); krat = np.float32(C["qf_krat"])
    inc = [bool(np.all(k[d:] / k[:nk - d] <= krat)) for d in range(nk)]
    exc = [bool(np.all(k[d:] / k[:nk - d] > krat)) for d in range(nk)]
    di_max = inc.index(False) - 1 if False in inc else nk - 1
    if not all(inc[:di_max + 1]) or not all(exc[di_max + 1:]):
        raise SystemExit("the k3/k1 window of WRT is not the same at every frequency on this grid")
    ad = Q["angle"].astype(np.float32); dmax = np.float32(C["qf_dmax"]); f180 = np.float32(180.0); dj_max = -1
    for dj in range(na // 2 + 1):
        a_dif = f180 - np.abs(f180 - np.abs(ad - np.roll(ad, -dj)))
        if np.all(a_dif <= dmax): dj_max = dj
        elif np.any(a_dif <= dmax): raise SystemExit("the direction window of WRT is not the same at every direction on this grid")
        else: break
    if 2 * dj_max >= na: raise SystemExit("the direction window of WRT reaches the opposite direction; this is not supported")
    return di_max, dj_max


def ref_bin(Q, di_max, dj_max):
    """reference bin of the geometry: nk // 2, or the lowest bin from which the lowest locus point (as stored at the highest
    anchor bin) stays inside the grid, if that is higher"""
    nk = Q["nk"]; jh = nk - di_max
    low = min(min(o[:, 0].min(), o[:, 2].min()) for di in range(di_max + 1) for dj in range(dj_max + 1)
              if not (di == 0 and dj == 0) for o in [offsets(Q, jh, jh + di, dj + 1)] if len(o))
    j = max(nk // 2, int(np.ceil(1 - low)))
    if j > jh: raise SystemExit(TOO_SHORT + f" (the loci reach {-low:.1f} bins below k1, and no bin up to {jh} holds them)")
    return j


def fold(Q, C, di_max, dj_max):
    NK, NA = Q["nk"], Q["na"]; JREF = ref_bin(Q, di_max, dj_max)
    rows, st = [], dict(fallback_bins=0, unmatched_nodes=0, unmatched_zz=0.0, total_zz=0.0, edge_fixed=[])
    for di in range(0, di_max + 1):
        for dj in range(0, dj_max + 1):
            if di == 0 and dj == 0: continue                       # k1 = k3 is skipped by WRT
            jr = JREF if JREF + di <= NK else NK - di
            ref = offsets(Q, jr, jr + di, dj + 1); rr = raw_nodes(Q, jr, jr + di, dj + 1)
            for m in np.flatnonzero((rr[:, 0] >= NK) | (rr[:, 2] >= NK)) if len(rr) else []:
                for j in range(jr - 1, 0, -1):                       # nearest lower bin with the node inside the grid
                    r2 = raw_nodes(Q, j, j + di, dj + 1)
                    if len(r2) == len(rr) and r2[m, 0] < NK and r2[m, 2] < NK:
                        st["edge_fixed"].append((di, dj, int(m) + 1, float(ref[m, 0]), float(offsets(Q, j, j + di, dj + 1)[m, 0]), j))
                        ref[m, :4] = offsets(Q, j, j + di, dj + 1)[m, :4]; break
                else:
                    raise SystemExit(TOO_SHORT + f" (configuration di = {di}, dj = {dj}: node {m + 1} is above the grid at every bin)")
            W = np.zeros((len(ref), NK))
            for j1 in range(1, NK + 1 - di):
                nd = offsets(Q, j1, j1 + di, dj + 1); st["total_zz"] += nd[:, 4].sum() if len(nd) else 0.0
                if len(nd) == len(ref): W[:, j1 - 1] = nd[:, 4]; continue
                st["fallback_bins"] += 1
                for q in nd:                                        # match each stored node to its reference node
                    d = np.abs(ref[:, :4] - q[:4]).max(1) if len(ref) else np.array([9.0])
                    if d.min() < TOL: W[d.argmin(), j1 - 1] = q[4]
                    else: st["unmatched_nodes"] += 1; st["unmatched_zz"] += q[4]
            for m in range(len(ref)):
                for s in ((1,) if dj in (0, NA // 2) else (1, -1)):
                    rows.append((di, s * dj, dj + 1, m, ref[m, 0], s * ref[m, 1], ref[m, 2], s * ref[m, 3], W[m]))
    t = {k: np.array([r[i] for r in rows], int) for i, k in enumerate(("di", "dj", "amem", "node"))}
    t.update({k: np.array([r[i] for r in rows], float) for i, k in enumerate(("u2", "v2", "u4", "v4"), 4)})
    t["zz"] = np.array([r[8] for r in rows])
    t.update(freq=C["freq"], sigma=C["sigma"], k=C["k"], cg=C["cg"], dk=C["dk"], dtheta=C["dtheta"],
             a=C["k"] * C["dtheta"] * C["dk"],                        # area factor k dtheta dk
             b=C["k"] / C["cg"] * 2 * np.pi * C["sigma"],             # N(kx,ky) -> E(f,theta) factor
             scale=2 * C["q_scale"], qf_frac=C["qf_frac"], kpow=C["kpow"], depth_ref=C["depth"])
    return t, st


def mirror_units(t):
    """mirror partner (theta -> -theta about k1): for dj != 0 the same node with the opposite sense; for dj = 0 the node
    of the same configuration at (u2, -v2, u4, -v4) (dj = 0 loci are symmetric, nodes match to 1e-5 bin).
    A unit = {row, mirror}; SWRT selects units, so the operator keeps WRT's left/right symmetry."""
    di, dj, node = t["di"], t["dj"], t["node"]; n = len(di)
    key = {(a, b, c): i for i, (a, b, c) in enumerate(zip(di, dj, node))}; mirror = np.full(n, -1)
    for i in range(n):
        if dj[i] != 0: mirror[i] = key[(di[i], -dj[i], node[i])]; continue
        same = np.flatnonzero((di == di[i]) & (dj == 0)); G = np.column_stack([t[k][same] for k in ("u2", "v2", "u4", "v4")])
        d = np.abs(G - np.array([t["u2"][i], -t["v2"][i], t["u4"][i], -t["v4"][i]])).max(1); mirror[i] = same[d.argmin()]
        assert d.min() < 1e-4
    assert (mirror[mirror] == np.arange(n)).all() and (mirror != np.arange(n)).all()
    return mirror, np.unique(np.minimum(np.arange(n), mirror), return_inverse=True)[1]


def build(bqf, C):
    """the table (dict of arrays) and a summary, from the database file and WRT's constants of the same grid"""
    Q = decode(bqf)
    if Q["nk"] != len(C["freq"]) or Q["na"] != C["na"] or not np.allclose(Q["sigma"], C["sigma"], rtol=1e-5):   # Xnl: single precision
        raise SystemExit(f"{bqf} does not belong to the grid of the constants")
    di_max, dj_max = window(Q, C)
    t, st = fold(Q, C, di_max, dj_max)
    if st["unmatched_zz"] > UNMATCHED_MAX * st["total_zz"]:
        raise SystemExit(TOO_SHORT + f" ({st['unmatched_zz'] / st['total_zz']:.1e} of WRT's weight has no place in the table)")
    t["mirror"], t["unit"] = mirror_units(t)
    t.update(nk=Q["nk"], na=Q["na"], di_max=di_max, dj_max=dj_max, nlocus=Q["nl"],
             **{k: C[k] for k in ("qf_krat", "qf_dmax", "iq_filt", "iq_lump", "iq_gauleg", "nlocus0", "g", "constants")})
    dm = [float(max(abs(t["u2"][i] - t["u2"][t["mirror"][i]]), abs(t["v2"][i] + t["v2"][t["mirror"][i]]),
                    abs(t["u4"][i] - t["u4"][t["mirror"][i]]), abs(t["v4"][i] + t["v4"][t["mirror"][i]]))) for i in np.flatnonzero(t["dj"] == 0)]
    summary = dict(nk=Q["nk"], na=Q["na"], header=Q["header"], di_max=di_max, dj_max=dj_max, ref_bin=ref_bin(Q, di_max, dj_max), rows=len(t["di"]),
                   configs=len(set(zip(t["di"], t["dj"]))), rows_dj0=int((t["dj"] == 0).sum()), units=int(t["unit"].max() + 1),
                   fallback_bins=st["fallback_bins"], unmatched_nodes=st["unmatched_nodes"], unmatched_zz_share=st["unmatched_zz"] / st["total_zz"],
                   edge_fixed=st["edge_fixed"], dj0_mirror_position_mismatch_max=max(dm) if dm else 0.0, constants=C["constants"],
                   bqf_valid_sha256=valid_hash(Q))
    return t, summary


def save(t, summary, work):
    np.savez(os.path.join(work, "table.npz"), **t,
             note="rows of WRT's own quadrature (quad99999.bqf), WRT's working window; zz raw (N(kx,ky) units)")
    json.dump(summary, open(os.path.join(work, "table_build.json"), "w"), indent=1)
    print(f"table: {summary['rows']} rows, {summary['units']} configurations (mirror pairs of rows) at {summary['configs']} offsets (di, dj) "
          f"(di <= {summary['di_max']}, |dj| <= {summary['dj_max']}), reference bin {summary['ref_bin']}, constants '{summary['constants']}'; "
          f"nodes matched by geometry in {summary['fallback_bins']} bins, unmatched {summary['unmatched_nodes']} (weight share {summary['unmatched_zz_share']:.1e})")


if __name__ == "__main__":
    import xnl_ref
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--f1", type=float, required=True); ap.add_argument("--ratio", type=float, required=True)
    ap.add_argument("--nf", type=int, required=True); ap.add_argument("--ndir", type=int, required=True)
    ap.add_argument("--work", required=True); ap.add_argument("--constants", default="paper", choices=("paper", "xnl"))
    a = ap.parse_args(); os.makedirs(a.work, exist_ok=True)
    freq, _ = xnl_ref.grid(a.f1, a.ratio, a.nf, a.ndir); xnl_ref.init(freq, a.ndir, os.path.join(a.work, "xnl"))
    t, s = build(os.path.join(a.work, "xnl", "quad99999.bqf"), xnl_ref.consts(a.constants)); save(t, s, a.work)

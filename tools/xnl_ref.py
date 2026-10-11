"""WRT for one spectral grid: Van Vledder's Xnl (WAVEWATCH III NL2), through the Python module xnl_wrt of build_xnl.sh.

Settings of the paper: iquad = 2 (deep-water loci, Webb's coupling coefficient), full circle (iq_grid = 3), f^-5 tail,
g = 9.81, and Xnl's defaults otherwise (30 points per locus, no lumping, no Gauss-Legendre, filtering on), which init()
checks. For iquad = 2, Xnl computes its loci and run-time arrays at its deep-water depth of 2000 m (q_maxdepth) whatever
depth it is given. The transfer of wrt() is called with a depth of 10 km, so that the WAM depth factor of iquad = 2 is 1 and
the result is the deep-water transfer that the table reproduces.
Xnl keeps its state in the module, so one process serves one grid and one work directory. It reads and writes its files
there: quad99999.bqf, the quadrature database that build_table.py folds, and xnl4v5.log/.prt/.tst. A database already
there is deleted first, because Xnl reuses a file whose header, directions and sigma (to 0.01 rad/s) agree with the grid,
so that a file left by another grid with as many frequencies and directions could be taken for this one. A configuration
file xnl4v5.cfg there would change Xnl's settings, and init() stops if it finds one.

    import xnl_ref
    xnl_ref.init(freq, ndir, work)    # freq (Hz, constant ratio), ndir directions over the full circle
    S = xnl_ref.wrt(E)                # E (nf, ndir) or (n, nf, ndir) energy density E(f, theta) -> S_nl, same layout
    C = xnl_ref.consts()              # WRT's run-time constants, for build_table.py"""
import os, sys, glob, contextlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TPI, G, TAIL, IQUAD, IQGRID = 2 * np.pi, 9.81, -5.0, 2, 3
DEPTH_INIT = 1000.0          # depth given to xnl_init, as in the paper (Xnl uses 2000 m for iquad = 2)
DEPTH_CONST = 1000.0         # depth of the paper's constants k, cg, dk (consts mode 'paper')
DEPTH_WRT = 1.0e4            # depth given to xnl_main: the WAM depth factor is then 1
DEFAULTS = dict(iq_lump=0, iq_gauleg=0, iq_filt=1, nlocus0=30)
_S = {}


@contextlib.contextmanager
def _cd(d):
    cwd = os.getcwd(); os.chdir(d)
    try: yield
    finally: os.chdir(cwd)


def _lib():
    sys.path.insert(0, os.path.join(HERE, "build"))
    try:
        import xnl_wrt
    except ImportError as e:
        raise SystemExit(f"xnl_wrt cannot be imported ({e}); build it first: bash tools/build_xnl.sh")
    return xnl_wrt


def grid(f1, ratio, nf, ndir):
    """frequencies (Hz) and direction centres (deg) of a grid, as the paper's scripts define them"""
    return f1 * ratio ** np.arange(nf), (360.0 / ndir) * (np.arange(ndir) + 0.5)


def init(freq, ndir, work):
    """set Xnl up for the grid (makes quad99999.bqf in work); later calls must give the same grid and work directory"""
    freq = np.asarray(freq, float); work = os.path.abspath(work)
    if _S:
        if not (np.array_equal(_S["freq"], freq) and _S["ndir"] == ndir and _S["work"] == work):
            raise SystemExit("xnl_ref serves one grid and one work directory per process")
        return
    lib = _lib(); os.makedirs(work, exist_ok=True)
    if os.path.exists(os.path.join(work, "xnl4v5.cfg")):
        raise SystemExit(f"{work}/xnl4v5.cfg would change Xnl's settings; remove it")
    for f in glob.glob(os.path.join(work, "quad*.bqf")): os.remove(f)
    sig = (TPI * freq).astype(np.float32)
    dirs = (360.0 / ndir) * (np.arange(ndir) + 0.5)            # deg; only the spacing matters to WRT
    with _cd(work):
        lib.m_constants.init_constants()
        ierr = lib.m_xnldata.xnl_init(sig, dirs.astype(np.float32), np.float32(TAIL), np.float32(G), np.array([DEPTH_INIT], np.float32),
                                      IQUAD, IQGRID, 0)
    if ierr != 0: raise SystemExit(f"xnl_init failed (ierr = {ierr}), see {work}/xnl4v5.log")
    if not os.path.exists(os.path.join(work, "quad99999.bqf")): raise SystemExit(f"Xnl wrote no quad99999.bqf in {work}")
    xn = lib.m_xnldata; got = {k: int(getattr(xn, k)) for k in DEFAULTS}
    if got != DEFAULTS: raise SystemExit(f"Xnl does not run with its default settings {DEFAULTS}: {got}")
    _S.update(lib=lib, freq=freq, ndir=ndir, sig=sig, dirs=dirs, ang=(np.arange(ndir) * TPI / ndir).astype(np.float32), work=work)


def wrt(E):
    """deep-water WRT transfer dE/dt(f, theta) of energy spectra E(f, theta)"""
    xn, SIG = _S["lib"].m_xnldata, _S["sig"]
    E = np.asarray(E, float); one = E.ndim == 2; E = E[None] if one else E
    out = np.empty_like(E)
    with _cd(_S["work"]):
        for n, e in enumerate(E):
            A = (e / (TPI * SIG[:, None])).astype(np.float32)          # A(sigma, theta) = E(f, theta) / (2 pi sigma)
            x, _, ierr = xn.xnl_main(A, SIG, _S["ang"], np.float32(DEPTH_WRT), IQUAD, 0)
            if ierr != 0: raise SystemExit(f"xnl_main failed (ierr = {ierr}), see {_S['work']}/xnl4v5.log")
            out[n] = x * (TPI * SIG[:, None])
    return out[0] if one else out


def consts(mode="paper"):
    """WRT's run-time constants and settings for build_table.py. mode 'paper': k, cg and dk from Xnl's functions at 1000 m,
    as in the paper (they differ from Xnl's run-time arrays by 2e-4 at 0.0345 Hz and more at lower frequencies); 'xnl':
    Xnl's own run-time arrays q_k, q_cg, q_dk and q_kpow (at 2000 m). q_scale is only set inside xnl_main; it is 1 here."""
    lib, freq, dirs = _S["lib"], _S["freq"], _S["dirs"]; xn, sv = lib.m_xnldata, lib.serv_xnl4v5
    fp = min(max(0.1, freq[2]), freq[-3])
    wrt(np.outer(np.exp(-((freq - fp) / (0.2 * fp)) ** 2), np.cos(np.deg2rad(dirs - 180) / 2) ** 2))
    SIG = _S["sig"]
    if mode == "paper":
        k = np.array([sv.z_wnumb(s, np.float32(DEPTH_CONST), np.float32(G)) for s in SIG], float)
        cg = np.array([sv.z_cmpcg(s, np.float32(DEPTH_CONST), np.float32(G)) for s in SIG], float)
        dk = np.asarray(sv.z_steps(k.astype(np.float32)), float); kpow = (k[0] / k) ** 7.5; depth = DEPTH_CONST
    elif mode == "xnl":
        k, cg, dk, kpow = (np.asarray(getattr(xn, v), float).copy() for v in ("q_k", "q_cg", "q_dk", "q_kpow")); depth = float(xn.q_depth)
    else:
        raise ValueError(mode)
    return dict(freq=freq, sigma=SIG.astype(float), k=k, cg=cg, dk=dk, dtheta=float(xn.q_delta), q_scale=float(xn.q_scale),
                qf_frac=float(xn.qf_frac), qf_krat=float(xn.qf_krat), qf_dmax=float(xn.qf_dmax), kpow=kpow,
                iq_dscale=int(xn.iq_dscale), iq_filt=int(xn.iq_filt), iq_lump=int(xn.iq_lump), iq_gauleg=int(xn.iq_gauleg),
                nlocus0=int(xn.nlocus0), depth=depth, g=G, na=_S["ndir"], constants=mode)

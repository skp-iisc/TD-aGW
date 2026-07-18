import numpy as np
import matplotlib.pyplot as plt

from qtm.constants import RYDBERG
from qtm.gspace import GSpace, GkSpace
from qtm.kpts import gen_monkhorst_pack_grid
from qtm.mpi import QTMComm
from qtm.dft import DFTCommMod, scf, KSHam, KSWfn
from qtm.containers import get_WavefunG, get_FieldG, FieldGType, FieldRType
from qtm.pseudo import NonlocGenerator
from qtm.pseudo.loc import loc_generate_pot_rhocore
from qtm.pot import hartree, xc
from qtm.mpi.utils import scatter_slice

from qtm.io_utils.dft_printers import print_scf_status

from cryst import sc_crystal

from qtm.config import MPI4PY_INSTALLED
if MPI4PY_INSTALLED:
    from mpi4py.MPI import COMM_WORLD
else:
    COMM_WORLD = None

comm_world = QTMComm(COMM_WORLD)

# Only k-pt parallelization:
dftcomm = DFTCommMod(comm_world, comm_world.size, 1)

crystal = sc_crystal('Si')

mpgrid_shape = (4,4,4)
kpts = gen_monkhorst_pack_grid(
    crystal, mpgrid_shape, (True, True, True),
    use_symm=False, is_time_reversal=False
)

ecut_wfn = 25 * RYDBERG
ecut_rho = 4 * ecut_wfn
grho = GSpace(crystal.recilat, ecut_rho)
gwfn = grho

numbnd = crystal.numel//2 + 16    # 16 vb + __ cb
conv_thr = 1e-8 * RYDBERG
diago_thr_init = 1e-2 * RYDBERG

E0 = 1e-3     # initial pulse in Ha/Bohr

# ─────────────────────────────────────────────────────────────────────────────
# Ground-state SCF (zero external field)
# ─────────────────────────────────────────────────────────────────────────────
out0 = scf(
    dftcomm,
    crystal,
    kpts,
    grho,
    gwfn,
    numbnd,
    is_spin=False,
    is_noncolin=False,
    symm_rho=True,
    rho_start=None,
    occ_typ="fixed",
    conv_thr=conv_thr,
    diago_thr_init=diago_thr_init,
    iter_printer=print_scf_status,
)
_, rho0, l_wfn0, en0 = out0
# l_wfn0 : list[list[KSWfn]]
#   outer list: k-points handled by this MPI rank
#   inner list: spin channels

i_kpts_kgrp = list(
    range(kpts.numkpts)[scatter_slice(kpts.numkpts, dftcomm.n_kgrp, dftcomm.i_kgrp)]
)

full_wfn_gs = [None] * kpts.numkpts
for ik_local, kswfn_k in enumerate(l_wfn0):
    full_wfn_gs[i_kpts_kgrp[ik_local]] = kswfn_k[0]

# ─────────────────────────────────────────────────────────────────────────────
# Reconstructing KSHam at each step
# ─────────────────────────────────────────────────────────────────────────────
FieldG_rho = get_FieldG(grho)
v_ion_g = FieldG_rho.zeros(())
rho_core = FieldG_rho.zeros(1)
l_nloc = []
for sp in crystal.l_atoms:
    v_ion_sp, rho_core_sp = loc_generate_pot_rhocore(sp, grho)
    v_ion_g += v_ion_sp
    rho_core += rho_core_sp
    l_nloc.append(NonlocGenerator(sp, gwfn))
v_ion = v_ion_g.to_r()

# XC functional identifier
libxc_func = xc.get_libxc_func(crystal)

def compute_vloc(rho_in: FieldGType) -> FieldRType:
    v_hart, _ = hartree.compute(rho_in)
    v_xc_r, _ = xc.compute(rho_in, rho_core, *libxc_func)
    vloc = v_ion + v_hart + v_xc_r
    vloc /= np.prod(gwfn.grid_shape)
    return vloc

def build_ksham(
    kswfn: KSWfn,
    vloc: FieldRType,
    efield_cart=None,
    ik_global: int = None,
) -> KSHam:
    """
    Construct a KSHam for a single k-point (spin-unpolarized).

    Parameters
    ----------
    kswfn : KSWfn
        Provides the GkSpace (gkspc) for this k-point.
    vloc : FieldRType
        Full local potential, shape () in real space (scalar, spin-up channel).
        Must already be divided by Nfft.
    efield_cart : array-like of length 3, or None
    ik_global : int or None
    """
    # vloc[0] selects the spin-up (only) channel
    if efield_cart is not None:
        return KSHam(
            kswfn.gkspc, False, vloc[0], l_nloc,
            efield_cart=list(efield_cart),
            full_wfn=full_wfn_gs,
            kpts=kpts,
            kgrid_shape=mpgrid_shape,
            ik=ik_global,
            crystal=crystal,
            current_kswfn=kswfn,
        )
    return KSHam(kswfn.gkspc, False, vloc[0], l_nloc)


def hamiltonian_matrix(ksham: KSHam, psi0: np.ndarray) -> np.ndarray:
    """
    Compute H_{nm,k} at a k-point.

    Parameters
    ----------
    ksham : KSHam
        Kohn-Sham Hamiltonian for this k-point.
    psi0 : ndarray, shape (numbnd, nG), complex128
        Normalised ground-state KS orbitals in the plane-wave basis.

    Returns
    -------
    H_k : ndarray, shape (numbnd, numbnd), complex128
    """
    gkspc = ksham.gkspc
    WavefunG = get_WavefunG(gkspc, 1)

    psi_wfn  = WavefunG(psi0.copy())
    hpsi_wfn = WavefunG(np.zeros_like(psi0))

    ksham.h_psi(psi_wfn, hpsi_wfn)
    H_k = psi0.conj() @ hpsi_wfn.data.T   # shape (numbnd, numbnd)
    return H_k

def drho_dt(rho_k: np.ndarray, H_k: np.ndarray) -> np.ndarray:
    return -1j * (H_k @ rho_k - rho_k @ H_k)

def rk4_step(rho_k: np.ndarray, H_k: np.ndarray, dt: float) -> np.ndarray:
    k1 = drho_dt(rho_k, H_k)
    k2 = drho_dt(rho_k + 0.5 * dt * k1, H_k)
    k3 = drho_dt(rho_k + 0.5 * dt * k2, H_k)
    k4 = drho_dt(rho_k + dt * k3,  H_k)
    return rho_k + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

# ─────────────────────────────────────────────────────────────────────────────
# Density matrix in the ground-state KS basis
# ─────────────────────────────────────────────────────────────────────────────
nk_local = len(l_wfn0)   # number of k-points on this MPI rank

# Ground-state wavefunctions
psi0_list = []      # psi0_list[ik_local] : (numbnd, nG) complex128
for kswfn_k in l_wfn0:
    kswfn = kswfn_k[0]
    psi = kswfn.evc_gk.data.copy()   # (numbnd, nG)
    
    norms = np.sqrt(np.einsum('ng,ng->n', psi.conj(), psi).real)
    psi /= norms[:, np.newaxis]     # Ensure normalisation
    psi0_list.append(psi)

rho_now = []
for kswfn_k in l_wfn0:
    occ = kswfn_k[0].occ.copy()
    rho_now.append(np.diag(occ.astype(complex)))   # (numbnd, numbnd)

# ─────────────────────────────────────────────────────────────────────────────
# Linear Response
# ─────────────────────────────────────────────────────────────────────────────

vel_x_list = []   # v_x[ik], shape (numbnd, numbnd)
evl_list   = []   # evl[ik], shape (numbnd,)
occ_list   = []   # occ[ik], shape (numbnd,)

print(f"\n[Computing momentum matrices for linear response...]\n")

for ik_local, kswfn_k in enumerate(l_wfn0):
    kswfn = kswfn_k[0]
    psi   = psi0_list[ik_local]
    
    evl = np.array(kswfn.evl, dtype=float)
    occ = np.array(kswfn.occ, dtype=float)
    evl_list.append(evl)
    occ_list.append(occ)

    gk_x = np.array(kswfn.gkspc.gk_cart[0], dtype=float)
    v_x  = (psi.conj() * gk_x) @ psi.T
    v_x  = 0.5 * (v_x + v_x.conj().T)  # Ensure exact Hermiticity
    vel_x_list.append(v_x)

# ─────────────────────────────────────────────────────────────────────────────
# Density Evolution and calculation of P(t)
# ─────────────────────────────────────────────────────────────────────────────
dt_prop = 0.02       # propagation dt [Ha]
n_t     = 10000      # T = 200 Ha
T_total = n_t * dt_prop

t_prop = np.arange(n_t + 1, dtype=float) * dt_prop
V_cell = grho.reallat_cellvol

print(f"[Density evolution: n_t={n_t}, dt={dt_prop}, T={T_total:.1f} Ha]")

P_t = np.zeros(n_t + 1)

for ik_local, kswfn_k in enumerate(l_wfn0):
    w_k = kswfn_k[0].k_weight
    evl = evl_list[ik_local]
    occ = occ_list[ik_local]
    v_x = vel_x_list[ik_local]
    
    deps = evl[:, None] - evl[None, :]
    df   = occ[None, :] - occ[:, None]
    
    with np.errstate(divide='ignore', invalid='ignore'):
        x_nm = -1j * v_x / deps
    x_nm[deps == 0] = 0.0
    
    drho_0 = -1j * E0 * x_nm * df
    
    obs_mat = -x_nm.T
    
    P_k = np.zeros(n_t + 1)
    for i, t in enumerate(t_prop):
        drho_t = drho_0 * np.exp(-1j * deps * t)
        P_k[i] = np.real(np.sum(drho_t * obs_mat))
        
    P_t += (2.0 / V_cell) * w_k * P_k

print(f"  P_x range: [{P_t.min():.4e}, {P_t.max():.4e}]")

# ─────────────────────────────────────────────────────────────────────────────
# Fourier transform for getting P_tilde(omega)
# ─────────────────────────────────────────────────────────────────────────────

sigma_smear  = 60.0
gaussian_win = np.exp(-0.5 * (t_prop / sigma_smear)**2)

omega     = np.linspace(1e-4, 1.0, 3000)   # frequency grid [Ha] (up to ~27 eV)
phase_ft  = np.exp(np.outer(1j * omega, t_prop))
P_tilde   = dt_prop * (phase_ft @ (P_t * gaussian_win))

eps   = 1.0 + 4.0 * np.pi * P_tilde / E0
eps1  = np.real(eps)
eps2  = np.imag(eps)

# ─────────────────────────────────────────────────────────────────────────────
# Outputs
# ─────────────────────────────────────────────────────────────────────────────
# import os
# os.makedirs("outputs", exist_ok=True)
# os.makedirs("save", exist_ok=True)

HA2EV    = 27.2114
omega_eV = omega * HA2EV

plt.rcParams.update({
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})

fig1, ax1 = plt.subplots(figsize=(6, 4), dpi=300)
ax1.plot(omega_eV, eps1, lw=1.5, color="C0")
ax1.set_xlabel(r"$\omega$ (eV)", fontsize=14)
ax1.set_ylabel(r"$\epsilon_1(\omega)$", fontsize=14)
fig1.savefig("outputs/si_dyn_eps1.png", bbox_inches="tight")
print("\n[Plot saved to outputs/si_dyn_eps1.png]")
plt.close(fig1)

fig2, ax2 = plt.subplots(figsize=(6, 4), dpi=300)
ax2.plot(omega_eV, eps2, lw=1.5, color="C0")
ax2.set_xlabel(r"$\omega$ (eV)", fontsize=14)
ax2.set_ylabel(r"$\epsilon_2(\omega)$", fontsize=14)
fig2.savefig("outputs/si_dyn_eps2.png", bbox_inches="tight")
print("[Plot saved to outputs/si_dyn_eps2.png]")
plt.close(fig2)

# try:
#     import pickle
#     with open("save/si_dyn_eps.pkl", "wb") as f:
#         pickle.dump({
#             "t_prop": t_prop, "P_t": P_t,
#             "omega": omega, "eps1": eps1, "eps2": eps2,
#             "E0": E0, "sigma_smear": sigma_smear,
#         }, f)
#     print("[Saved to save/si_dyn_eps.pkl]")
# except Exception as e:
#     print(f"[ERROR: {e}]")

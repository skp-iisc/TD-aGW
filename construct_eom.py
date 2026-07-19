import numpy as np
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

# ─────────────────────────────────────────────────────────────────────────────
# System setup
# ─────────────────────────────────────────────────────────────────────────────
crystal = sc_crystal('Si', alat=10.2612)

mpgrid_shape = (10,10,10)
kpts = gen_monkhorst_pack_grid(
    crystal, mpgrid_shape, (True, True, True),
    use_symm=False, is_time_reversal=False
)

ecut_wfn = 25 * RYDBERG
ecut_rho = 4 * ecut_wfn
grho = GSpace(crystal.recilat, ecut_rho)
gwfn = grho

numbnd = crystal.numel // 2 + 16
conv_thr = 1e-8 * RYDBERG
diago_thr_init = 1e-2 * RYDBERG

E0 = 1e-4     # kick amplitude [Ha/Bohr]

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Ground-state SCF (zero external field)
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
#   inner list: spin channels (length 1 for spin-unpolarized)

i_kpts_kgrp = list(
    range(kpts.numkpts)[scatter_slice(kpts.numkpts, dftcomm.n_kgrp, dftcomm.i_kgrp)]
)

full_wfn_gs = [None] * kpts.numkpts
for ik_local, kswfn_k in enumerate(l_wfn0):
    full_wfn_gs[i_kpts_kgrp[ik_local]] = kswfn_k[0]

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Precompute static objects needed to reconstruct KSHam at each step
# ─────────────────────────────────────────────────────────────────────────────
# Ionic local potential + core density
FieldG_rho = get_FieldG(grho)
v_ion_g = FieldG_rho.zeros(())
rho_core = FieldG_rho.zeros(1)
l_nloc = []
for sp in crystal.l_atoms:
    v_ion_sp, rho_core_sp = loc_generate_pot_rhocore(sp, grho)
    v_ion_g += v_ion_sp
    rho_core += rho_core_sp
    l_nloc.append(NonlocGenerator(sp, gwfn))
v_ion = v_ion_g.to_r()   # real-space ionic potential (unnormalised by Nfft)

# XC functional identifier
libxc_func = xc.get_libxc_func(crystal)


def compute_vloc(rho_in: FieldGType) -> FieldRType:
    v_hart, _ = hartree.compute(rho_in)
    v_xc_r, _ = xc.compute(rho_in, rho_core, *libxc_func)
    vloc = v_ion + v_hart + v_xc_r
    vloc /= np.prod(gwfn.grid_shape)   # normalise to match KSHam convention
    return vloc

def build_ksham(
    kswfn: KSWfn,
    vloc: FieldRType,
    efield_cart=None,
    ik_global: int = None,
) -> KSHam:
    # vloc[0] selects the spin-up (only) channel; shape () after indexing
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
    Compute H_{nm,k}.

    Parameters
    ----------
    ksham : KSHam
    psi0 : ndarray, shape (numbnd, nG), complex128

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
    k1 = drho_dt(rho_k,                  H_k)
    k2 = drho_dt(rho_k + 0.5 * dt * k1, H_k)
    k3 = drho_dt(rho_k + 0.5 * dt * k2, H_k)
    k4 = drho_dt(rho_k +       dt * k3,  H_k)
    return rho_k + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

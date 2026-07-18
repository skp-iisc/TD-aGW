import numpy as np
import os, sys

from qtm.constants import RYDBERG
from qtm.lattice import RealLattice
from qtm.crystal import BasisAtoms, Crystal
from qtm.pseudo import UPFv2Data

from materials import MATERIALS

def sc_crystal(material_name, alat=None):
    """
    Set up crystal for III-V semiconductors.

    Inputs
    ----------
    material_name : str - 'Si', 'AlAs', 'GaAs', 'GaP', or 'AlP'
    alat : float or None - lattice constant in Bohr

    Returns
    ----------
    Crystal(reallat, [cation_atoms, anion_atoms])
    """
    mat = MATERIALS[material_name]
    if alat is None:
        alat = mat['alat']

    # Simple Cubic (SC) primitive vectors (in units of alat)
    a1 = np.array([1.0, 0.0, 0.0])
    a2 = np.array([0.0, 1.0, 0.0])
    a3 = np.array([0.0, 0.0, 1.0])

    reallat = RealLattice.from_alat(alat=alat, a1=a1, a2=a2, a3=a3)

    rc1 = np.array([-0.125, -0.125, -0.125])
    rc2 = np.array([0.375,  0.375, -0.125])
    rc3 = np.array([0.375, -0.125,  0.375])
    rc4 = np.array([-0.125,  0.375,  0.375])
    ra1 = np.array([0.125,  0.125,  0.125])
    ra2 = np.array([0.625,  0.625,  0.125])
    ra3 = np.array([0.625,  0.125,  0.625])
    ra4 = np.array([0.125,  0.625,  0.625])

    rc = np.column_stack([rc1, rc2, rc3, rc4])
    ra = np.column_stack([ra1, ra2, ra3, ra4])

    pp_dir = os.path.dirname(os.path.abspath(__file__))
    cation_pp = UPFv2Data.from_file(os.path.join(pp_dir, mat['cation_pp']))
    anion_pp = UPFv2Data.from_file(os.path.join(pp_dir, mat['anion_pp']))

    cations = BasisAtoms(mat['cation'], cation_pp, mat['cation_mass'],
                        reallat, rc)
    anions = BasisAtoms(mat['anion'], anion_pp, mat['anion_mass'],
                        reallat, ra)

    return Crystal(reallat, [cations, anions])

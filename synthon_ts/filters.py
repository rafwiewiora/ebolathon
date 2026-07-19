"""Local drug-like pre-dock filter.

onepot's *analog* search does not constrain physicochemical properties, so the
molecules it returns can drift out of a sensible drug-like window. Docking is the
expensive step, so we filter analogs locally (RDKit) BEFORE they reach the oracle
— this saves dock cost and keeps hits drug-like.

Default window (hard cuts):
    MW              <= 550     (soft-ideal <= 500, not enforced)
    cLogP           1 .. 5     (soft-ideal < 4.5, not enforced)
    HBD             <= 5
    HBA             <= 10
    TPSA            <= 140
    rotatable bonds <= 8

CAVEAT: the medicinal-chemistry target is cLogD(pH 7.4) 1-5, but neither RDKit nor
onepot expose cLogD. We use RDKit ``Crippen.MolLogP`` (a cLogP estimate) as a
proxy for the cLogD window here, mirroring onepot ``sample_space``'s ``clogp``
filter. cLogP ignores ionisation, so this is an approximation.
"""
from __future__ import annotations

from dataclasses import dataclass

try:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
    _HAVE_RDKIT = True
except Exception:  # pragma: no cover
    _HAVE_RDKIT = False


@dataclass
class DrugLikeLimits:
    """Configurable drug-like window. `logp_*` is a cLogP proxy for cLogD 1-5."""
    mw_max: float = 550.0
    logp_min: float = 1.0
    logp_max: float = 5.0
    hbd_max: int = 5
    hba_max: int = 10
    tpsa_max: float = 140.0
    rotatable_max: int = 8


DEFAULT_LIMITS = DrugLikeLimits()


def passes_druglike(mol_or_smiles, limits: DrugLikeLimits = DEFAULT_LIMITS) -> bool:
    """True if the molecule sits inside the drug-like window.

    Accepts an RDKit ``Mol`` or a SMILES string. If RDKit is unavailable, or the
    SMILES cannot be parsed, returns True (fail-open — never silently drop a
    ligand just because we couldn't measure it)."""
    if not _HAVE_RDKIT:
        return True
    mol = mol_or_smiles
    if isinstance(mol_or_smiles, str):
        mol = Chem.MolFromSmiles(mol_or_smiles)
    if mol is None:
        return True
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    return (mw <= limits.mw_max
            and limits.logp_min <= logp <= limits.logp_max
            and hbd <= limits.hbd_max
            and hba <= limits.hba_max
            and tpsa <= limits.tpsa_max
            and rot <= limits.rotatable_max)


def make_filter(limits: DrugLikeLimits = DEFAULT_LIMITS):
    """Return a ``smiles -> bool`` callable bound to `limits`, suitable for
    passing as ``run_loop(..., mol_filter=...)``."""
    return lambda smi: passes_druglike(smi, limits)

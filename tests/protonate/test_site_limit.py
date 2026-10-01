"""Checks the per-pattern limit on how many sites are protonated."""

from rdkit import Chem

from dimorphite_dl import protonate_smiles

VERY_BASIC_PH = 10000000.0


def test_site_limit_is_inclusive() -> None:
    """Checks that a molecule with exactly max_sites_per_molecule (50)
    carboxyls has all of them deprotonated. The limit was checked with >=,
    so the 50th group was left neutral."""

    smiles = "CC(C(=O)O)" * 50

    output = protonate_smiles(
        smiles, ph_min=VERY_BASIC_PH, ph_max=VERY_BASIC_PH, precision=0.5
    )

    assert len(output) == 1, output
    mol = Chem.MolFromSmiles(output[0])
    assert mol is not None
    assert Chem.GetFormalCharge(mol) == -50, output[0]

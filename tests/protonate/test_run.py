import pytest
from conftest import compare_smiles  # type: ignore

from dimorphite_dl.protonate import protonate_smiles


# Every molecule should be protonated
@pytest.mark.parametrize(
    ("smiles_input", "smiles_correct"),
    [
        ("C#CCO", "C#CCO"),  # alcohol
        ("C(=O)N", "NC=O"),  # Amide,
        ("CC(=O)NOC(C)=O", "CC(=O)NOC(C)=O"),  # Amide_electronegative,
        ("COC(=N)N", "COC(N)=[NH2+]"),  # AmidineGuanidine2,
        (
            "Brc1ccc(C2NCCS2)cc1",
            "Brc1ccc(C2[NH2+]CCS2)cc1",
        ),  # Amines_primary_secondary_tertiary,
        ("CC(=O)[n+]1ccc(N)cc1", "CC(=O)[n+]1ccc([NH3+])cc1"),  # Anilines_primary,
        ("CCNc1ccccc1", "CC[NH2+]c1ccccc1"),  # Anilines_secondary,
        ("Cc1ccccc1N(C)C", "Cc1ccccc1[NH+](C)C"),  # Anilines_tertiary,
        ("BrC1=CC2=C(C=C1)NC=C2", "Brc1ccc2[nH]ccc2c1"),  # Indole_pyrrole,
        ("O=c1cc[nH]cc1", "O=c1cc[nH]cc1"),  # Aromatic_nitrogen_protonated,
        ("C-N=[N+]=[N@H]", "CN=[N+]=N"),  # Azide,
        ("BrC(C(O)=O)CBr", "O=C(O)C(Br)CBr"),  # Carboxyl,
        ("NC(NN=O)=N", "NC(=[NH2+])NN=O"),  # AmidineGuanidine1,
        ("C(F)(F)(F)C(=O)NC(=O)C", "CC(=O)NC(=O)C(F)(F)F"),  # Imide,
        ("O=C(C)NC(C)=O", "CC(=O)NC(C)=O"),  # Imide2,
        ("CC(C)(C)C(N(C)O)=O", "CN(O)C(=O)C(C)(C)C"),  # N-hydroxyamide,
        ("C[N+](O)=O", "C[N+](=O)O"),  # Nitro,
        ("O=C1C=C(O)CC1", "O=C1C=C(O)CC1"),  # O=C-C=C-OH,
        ("C1CC1OO", "OOC1CC1"),  # Peroxide2,
        ("C(=O)OO", "O=COO"),  # Peroxide1,
        ("Brc1cc(O)cc(Br)c1", "Oc1cc(Br)cc(Br)c1"),  # Phenol,
        ("CC(=O)c1ccc(S)cc1", "CC(=O)c1ccc(S)cc1"),  # Phenyl_Thiol,
        ("C=CCOc1ccc(C(=O)O)cc1", "C=CCOc1ccc(C(=O)O)cc1"),  # Phenyl_carboxyl,
        ("COP(=O)(O)OC", "COP(=O)(O)OC"),  # Phosphate_diester,
        ("CP(C)(=O)O", "CP(C)(=O)O"),  # Phosphinic_acid,
        ("CC(C)OP(C)(=O)O", "CC(C)OP(C)(=O)O"),  # Phosphonate_ester,
        ("CC1(C)OC(=O)NC1=O", "CC1(C)OC(=O)NC1=O"),  # Ringed_imide1,
        ("O=C(N1)C=CC1=O", "O=C1C=CC(=O)N1"),  # Ringed_imide2,
        ("O=S(OC)(O)=O", "COS(=O)(=O)O"),  # Sulfate,
        ("COc1ccc(S(=O)O)cc1", "COc1ccc(S(=O)O)cc1"),  # Sulfinic_acid,
        ("CS(N)(=O)=O", "CS(N)(=O)=O"),  # Sulfonamide,
        ("CC(=O)CSCCS(O)(=O)=O", "CC(=O)CSCCS(=O)(=O)O"),  # Sulfonate,
        ("CC(=O)S", "CC(=O)S"),  # Thioic_acid,
        ("C(C)(C)(C)(S)", "CC(C)(C)S"),  # Thiol,
        ("Brc1cc[nH+]cc1", "Brc1cc[nH+]cc1"),  # Aromatic_nitrogen_unprotonated,
        ("C=C(O)c1c(C)cc(C)cc1C", "C=C(O)c1c(C)cc(C)cc1C"),  # Vinyl_alcohol,
        ("CC(=O)ON", "CC(=O)O[NH3+]"),  # Primary_hydroxyl_amine,
        ("O=P(O)(O)OCCCC", "CCCCOP(=O)(O)O"),  # Phosphate
        ("CC(P(O)(O)=O)C", "CC(C)P(=O)(O)O"),  # Phosphonate
    ],
)
def test_very_acidic_single(smiles_input, smiles_correct):
    ph_min = -10000000
    ph_max = -10000000

    output = protonate_smiles(smiles_input, ph_min=ph_min, ph_max=ph_max, precision=0.5)
    assert len(output) == 1
    smiles_output = output[0]

    compare_smiles(smiles_output, smiles_correct)


# Every molecule should be deprotonated
@pytest.mark.parametrize(
    ("smiles_input", "smiles_correct"),
    [
        ("C#CCO", "C#CC[O-]"),  # Alcohol
        ("C(=O)N", "[NH-]C=O"),  # Amide
        ("CC(=O)NOC(C)=O", "CC(=O)[N-]OC(C)=O"),  # Amide_electronegative
        ("COC(=N)N", "COC(=N)N"),  # AmidineGuanidine2
        (
            "Brc1ccc(C2NCCS2)cc1",
            "Brc1ccc(C2NCCS2)cc1",
        ),  # Amines_primary_secondary_tertiary
        ("CC(=O)[n+]1ccc(N)cc1", "CC(=O)[n+]1ccc(N)cc1"),  # Anilines_primary
        ("CCNc1ccccc1", "CCNc1ccccc1"),  # Anilines_secondary
        ("Cc1ccccc1N(C)C", "Cc1ccccc1N(C)C"),  # Anilines_tertiary
        ("BrC1=CC2=C(C=C1)NC=C2", "Brc1ccc2[n-]ccc2c1"),  # Indole_pyrrole
        ("O=c1cc[nH]cc1", "O=c1cc[n-]cc1"),  # Aromatic_nitrogen_protonated
        ("C-N=[N+]=[N@H]", "CN=[N+]=[N-]"),  # Azide
        ("BrC(C(O)=O)CBr", "O=C([O-])C(Br)CBr"),  # Carboxyl
        ("NC(NN=O)=N", "N=C(N)NN=O"),  # AmidineGuanidine1
        ("C(F)(F)(F)C(=O)NC(=O)C", "CC(=O)[N-]C(=O)C(F)(F)F"),  # Imide
        ("O=C(C)NC(C)=O", "CC(=O)[N-]C(C)=O"),  # Imide2
        ("CC(C)(C)C(N(C)O)=O", "CN([O-])C(=O)C(C)(C)C"),  # N-hydroxyamide
        ("C[N+](O)=O", "C[N+](=O)[O-]"),  # Nitro
        ("O=C1C=C(O)CC1", "O=C1C=C([O-])CC1"),  # O=C-C=C-OH
        ("C1CC1OO", "[O-]OC1CC1"),  # Peroxide2
        ("C(=O)OO", "O=CO[O-]"),  # Peroxide1
        ("Brc1cc(O)cc(Br)c1", "[O-]c1cc(Br)cc(Br)c1"),  # Phenol
        ("CC(=O)c1ccc(S)cc1", "CC(=O)c1ccc([S-])cc1"),  # Phenyl_Thiol
        ("C=CCOc1ccc(C(=O)O)cc1", "C=CCOc1ccc(C(=O)[O-])cc1"),  # Phenyl_carboxyl
        ("COP(=O)(O)OC", "COP(=O)([O-])OC"),  # Phosphate_diester
        ("CP(C)(=O)O", "CP(C)(=O)[O-]"),  # Phosphinic_acid
        ("CC(C)OP(C)(=O)O", "CC(C)OP(C)(=O)[O-]"),  # Phosphonate_ester
        ("CC1(C)OC(=O)NC1=O", "CC1(C)OC(=O)[N-]C1=O"),  # Ringed_imide1
        ("O=C(N1)C=CC1=O", "O=C1C=CC(=O)[N-]1"),  # Ringed_imide2
        ("O=S(OC)(O)=O", "COS(=O)(=O)[O-]"),  # Sulfate
        ("COc1ccc(S(=O)O)cc1", "COc1ccc(S(=O)[O-])cc1"),  # Sulfinic_acid
        ("CS(N)(=O)=O", "CS([NH-])(=O)=O"),  # Sulfonamide
        ("CC(=O)CSCCS(O)(=O)=O", "CC(=O)CSCCS(=O)(=O)[O-]"),  # Sulfonate
        ("CC(=O)S", "CC(=O)[S-]"),  # Thioic_acid
        ("C(C)(C)(C)(S)", "CC(C)(C)[S-]"),  # Thiol
        ("Brc1cc[nH+]cc1", "Brc1ccncc1"),  # Aromatic_nitrogen_unprotonated
        ("C=C(O)c1c(C)cc(C)cc1C", "C=C([O-])c1c(C)cc(C)cc1C"),  # Vinyl_alcohol
        ("CC(=O)ON", "CC(=O)ON"),  # Primary_hydroxyl_amine
        ("O=P(O)(O)OCCCC", "CCCCOP(=O)([O-])[O-]"),  # Phosphate
        ("CC(P(O)(O)=O)C", "CC(C)P(=O)([O-])[O-]"),  # Phosphonate
    ],
)
def test_very_basic(smiles_input, smiles_correct):
    ph_min = 10000000
    ph_max = 10000000

    output = list(
        protonate_smiles(smiles_input, ph_min=ph_min, ph_max=ph_max, precision=0.5)
    )
    assert len(output) == 1
    smiles_output = output[0]

    compare_smiles(smiles_output, smiles_correct)


@pytest.mark.parametrize(
    ("smiles_input", "smiles_protonated", "smiles_deprotonated", "pka_avg"),
    [
        ["C#CCO", "C#CCO", "C#CC[O-]", 14.780384615384616],  # alcohol
        ["C(=O)N", "NC=O", "[NH-]C=O", 12.00611111111111],  # amide
        [
            "CC(=O)NOC(C)=O",
            "CC(=O)NOC(C)=O",
            "CC(=O)[N-]OC(C)=O",
            3.4896,
        ],  # Amide_electronegative
        [
            "COC(=N)N",
            "COC(N)=[NH2+]",
            "COC(=N)N",
            10.035538461538462,
        ],  # AmidineGuanidine2"
        [
            "Brc1ccc(C2NCCS2)cc1",
            "Brc1ccc(C2[NH2+]CCS2)cc1",
            "Brc1ccc(C2NCCS2)cc1",
            8.159107682388349,
        ],  # Amines_primary_secondary_tertiary
        [
            "CC(=O)[n+]1ccc(N)cc1",
            "CC(=O)[n+]1ccc([NH3+])cc1",
            "CC(=O)[n+]1ccc(N)cc1",
            3.899298673194805,
        ],  # Anilines_primary
        [
            "CCNc1ccccc1",
            "CC[NH2+]c1ccccc1",
            "CCNc1ccccc1",
            4.335408163265306,
        ],  # Anilines_secondary
        [
            "Cc1ccccc1N(C)C",
            "Cc1ccccc1[NH+](C)C",
            "Cc1ccccc1N(C)C",
            4.16690685045614,
        ],  # Anilines_tertiary
        [
            "BrC1=CC2=C(C=C1)NC=C2",
            "Brc1ccc2[nH]ccc2c1",
            "Brc1ccc2[n-]ccc2c1",
            14.52875,
        ],  # Indole_pyrrole
        [
            "O=c1cc[nH]cc1",
            "O=c1cc[nH]cc1",
            "O=c1cc[n-]cc1",
            7.17,
        ],  # Aromatic_nitrogen_protonated
        ["C-N=[N+]=[N@H]", "CN=[N+]=N", "CN=[N+]=[N-]", 4.65],  # Azide
        [
            "BrC(C(O)=O)CBr",
            "O=C(O)C(Br)CBr",
            "O=C([O-])C(Br)CBr",
            3.456652971502591,
        ],  # Carboxyl
        [
            "NC(NN=O)=N",
            "NC(=[NH2+])NN=O",
            "N=C(N)NN=O",
            12.025333333333334,
        ],  # AmidineGuanidine1
        [
            "C(F)(F)(F)C(=O)NC(=O)C",
            "CC(=O)NC(=O)C(F)(F)F",
            "CC(=O)[N-]C(=O)C(F)(F)F",
            2.466666666666667,
        ],  # Imide
        ["O=C(C)NC(C)=O", "CC(=O)NC(C)=O", "CC(=O)[N-]C(C)=O", 10.23],  # Imide2
        [
            "CC(C)(C)C(N(C)O)=O",
            "CN(O)C(=O)C(C)(C)C",
            "CN([O-])C(=O)C(C)(C)C",
            9.301904761904762,
        ],  # N-hydroxyamide
        ["C[N+](O)=O", "C[N+](=O)O", "C[N+](=O)[O-]", -1000.0],  # Nitro
        ["C1CC1OO", "OOC1CC1", "[O-]OC1CC1", 11.978235294117647],  # Peroxide2
        ["C(=O)OO", "O=COO", "O=CO[O-]", 8.738888888888889],  # Peroxide1
        [
            "Brc1cc(O)cc(Br)c1",
            "Oc1cc(Br)cc(Br)c1",
            "[O-]c1cc(Br)cc(Br)c1",
            7.065359866910526,
        ],  # Phenol
        [
            "CC(=O)c1ccc(S)cc1",
            "CC(=O)c1ccc(S)cc1",
            "CC(=O)c1ccc([S-])cc1",
            4.978235294117647,
        ],  # Phenyl_Thiol
        [
            "C=CCOc1ccc(C(=O)O)cc1",
            "C=CCOc1ccc(C(=O)O)cc1",
            "C=CCOc1ccc(C(=O)[O-])cc1",
            3.463441968255319,
        ],  # Phenyl_carboxyl
        [
            "COP(=O)(O)OC",
            "COP(=O)(O)OC",
            "COP(=O)([O-])OC",
            2.7280434782608696,
        ],  # Phosphate_diester
        ["CP(C)(=O)O", "CP(C)(=O)O", "CP(C)(=O)[O-]", 2.9745],  # Phosphinic_acid
        [
            "CC(C)OP(C)(=O)O",
            "CC(C)OP(C)(=O)O",
            "CC(C)OP(C)(=O)[O-]",
            2.0868,
        ],  # Phosphonate_ester
        [
            "CC1(C)OC(=O)NC1=O",
            "CC1(C)OC(=O)NC1=O",
            "CC1(C)OC(=O)[N-]C1=O",
            6.4525,
        ],  # Ringed_imide1
        [
            "O=C(N1)C=CC1=O",
            "O=C1C=CC(=O)N1",
            "O=C1C=CC(=O)[N-]1",
            8.681666666666667,
        ],  # Ringed_imide2
        ["O=S(OC)(O)=O", "COS(=O)(=O)O", "COS(=O)(=O)[O-]", -2.36],  # Sulfate
        [
            "COc1ccc(S(=O)O)cc1",
            "COc1ccc(S(=O)O)cc1",
            "COc1ccc(S(=O)[O-])cc1",
            1.7933333333333332,
        ],  # Sulfinic_acid
        [
            "CS(N)(=O)=O",
            "CS(N)(=O)=O",
            "CS([NH-])(=O)=O",
            7.9160326086956525,
        ],  # Sulfonamide
        [
            "CC(=O)CSCCS(O)(=O)=O",
            "CC(=O)CSCCS(=O)(=O)O",
            "CC(=O)CSCCS(=O)(=O)[O-]",
            -1.8184615384615386,
        ],  # Sulfonate
        ["CC(=O)S", "CC(=O)S", "CC(=O)[S-]", 0.678267],  # Thioic_acid
        ["C(C)(C)(C)(S)", "CC(C)(C)S", "CC(C)(C)[S-]", 9.12448275862069],  # Thiol
        [
            "Brc1cc[nH+]cc1",
            "Brc1cc[nH+]cc1",
            "Brc1ccncc1",
            4.3535441240733945,
        ],  # Aromatic_nitrogen_unprotonated
        [
            "C=C(O)c1c(C)cc(C)cc1C",
            "C=C(O)c1c(C)cc(C)cc1C",
            "C=C([O-])c1c(C)cc(C)cc1C",
            8.871850714285713,
        ],  # Vinyl_alcohol
        [
            "CC(=O)ON",
            "CC(=O)O[NH3+]",
            "CC(=O)ON",
            4.035714285714286,
        ],  # Primary_hydroxyl_amine
    ],
)
def test_pka_average(smiles_input, smiles_protonated, smiles_deprotonated, pka_avg):
    """Test that when the pH is equal to the average pKa, the protonation
    state is always both"""

    output = list(
        protonate_smiles(smiles_input, ph_min=pka_avg, ph_max=pka_avg, precision=0.5)
    )
    assert len(output) == 2
    smiles_output_sorted = tuple(sorted(output))
    smiles_correct_sorted = tuple(sorted((smiles_protonated, smiles_deprotonated)))
    for smiles_output, smiles_correct in zip(
        smiles_output_sorted, smiles_correct_sorted
    ):
        compare_smiles(smiles_output, smiles_correct)


def test_no_carbanion():
    smi = (
        "Cc1nc2cc(-c3[nH]c4cc5ccccc5c5c4c3CCN(C(=O)O)[C@@H]5O)cc3c(=O)[nH][nH]c(n1)c23"
    )
    output = list(protonate_smiles(smi))

    if "[C-]" in "".join(output).upper():
        msg = "Processing " + smi + " produced a molecule with a carbanion!"
        raise RuntimeError(msg)
    else:
        print("(CORRECT) No carbanion: " + smi)


def test_max_variants():
    # Make sure max number of variants is limited (old bug).
    smi = "CCCC[C@@H](C(=O)N)NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](C(C)C)NC(=O)[C@@H](NC(=O)[C@H](Cc1c[nH]c2c1cccc2)NC(=O)[C@@H](NC(=O)[C@@H](Cc1ccc(cc1)O)N)CCC(=O)N)C)C)Cc1nc[nH]c1)Cc1ccccc1"
    output = list(protonate_smiles(smi))

    assert len(output) == 128, f"Should produce 128 mol, but produced {len(output)}"


@pytest.mark.parametrize(
    ("smiles_input", "ph", "smiles_correct"),
    [
        (
            "O=P(O)(OP(O)(OP(O)(OCC1OC(C(C1O)O)N2C=NC3=C2N=CN=C3N)=O)=O)O",
            0.5,
            "[NH3+]c1[nH+]c[nH+]c2c1[nH+]cn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O",
        ),
        (
            "O=P(O)(OP(O)(OP(O)(OCC1OC(C(C1O)O)N2C=NC3=C2N=CN=C3N)=O)=O)O",
            1.0,
            "[NH3+]c1[nH+]c[nH+]c2c1[nH+]cn2C1OC(COP(=O)(O)OP(=O)([O-])OP(=O)(O)O)C(O)C1O",
        ),
        (
            "O=P(O)(OP(O)(OP(O)(OCC1OC(C(C1O)O)N2C=NC3=C2N=CN=C3N)=O)=O)O",
            2.6,
            "[NH3+]c1[nH+]c[nH+]c2c1[nH+]cn2C1OC(COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])O)C(O)C1O",
        ),
        (
            "O=P(O)(OP(O)(OP(O)(OCC1OC(C(C1O)O)N2C=NC3=C2N=CN=C3N)=O)=O)O",
            7.0,
            "Nc1ncnc2c1ncn2C1OC(COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])C(O)C1O",
        ),
        (
            "O=P(O)(OP(O)(OCC1C(O)C(O)C(N2C=NC3=C(N)N=CN=C32)O1)=O)OCC(O4)C(O)C(O)C4[N+]5=CC=CC(C(N)=O)=C5",
            0.5,
            "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c([NH3+])ncnc54)C(O)C3O)C(O)C2O)c1",
        ),
        (
            "O=P(O)(OP(O)(OCC1C(O)C(O)C(N2C=NC3=C(N)N=CN=C32)O1)=O)OCC(O4)C(O)C(O)C4[N+]5=CC=CC(C(N)=O)=C5",
            2.5,
            "NC(=O)c1ccc[n+](C2OC(COP(=O)([O-])OP(=O)([O-])OCC3OC(n4cnc5c([NH3+])ncnc54)C(O)C3O)C(O)C2O)c1",
        ),
        (
            "O=P(O)(OP(O)(OCC1C(O)C(O)C(N2C=NC3=C(N)N=CN=C32)O1)=O)OCC(O4)C(O)C(O)C4[N+]5=CC=CC(C(N)=O)=C5",
            7.4,
            "NC(=O)c1ccc[n+](C2OC(COP(=O)([O-])OP(=O)([O-])OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1",
        ),
    ],
)
def test_multiple_ph(smiles_input, ph, smiles_correct):
    output = protonate_smiles(
        smiles_input, ph_min=ph, ph_max=ph, precision=0.0, validate_output=True
    )
    print(output)
    assert len(output) == 1

    compare_smiles(output[0], smiles_correct)

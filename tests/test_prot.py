from dimorphite_dl.protonate import protonate_smiles


def check_tests(smiles, kwargs, expected_output, labels):
    """Tests most ionizable groups. The ones that can only loose or gain a single proton.

    Args:
        smiles: SMILES to protonate.
        kwargs: The keyword arguments to pass to `protonate_smiles`.
        expected_output: A list of the expected SMILES-strings output.
        labels: The labels. A list containing combo of BOTH, PROTONATED,
                DEPROTONATED.
    """

    output = list(protonate_smiles(smiles, **kwargs))
    output = [o.split() for o in output]

    num_states = len(expected_output)

    if len(output) != num_states:
        msg = (
            "smiles"
            + " should have "
            + str(num_states)
            + " states at at pH "
            + str(kwargs["min_ph"])
            + ": "
            + str(output)
        )
        raise Exception(msg)

    if len(set([l[0] for l in output]) - set(expected_output)) != 0:
        msg = (
            smiles
            + " is not "
            + " AND ".join(expected_output)
            + " at pH "
            + str(kwargs["min_ph"])
            + " to "
            + str(kwargs["max_ph"])
            + "; it is "
            + " AND ".join([l[0] for l in output])
        )
        raise Exception(msg)

    if len(set([l[1] for l in output]) - set(labels)) != 0:
        msg = (
            smiles
            + " not labeled as "
            + " AND ".join(labels)
            + "; it is "
            + " AND ".join([l[1] for l in output])
        )
        raise Exception(msg)

    ph_range = sorted(list(set([kwargs["min_ph"], kwargs["max_ph"]])))
    ph_range_str = "(" + " - ".join("{0:.2f}".format(n) for n in ph_range) + ")"
    print(
        "(CORRECT) "
        + ph_range_str.ljust(10)
        + " "
        + smiles
        + " => "
        + " AND ".join([l[0] for l in output])
    )


def test_very_acidic(kwargs_default, smiles_groups, smiles_phosphates):
    kwargs = kwargs_default
    kwargs["min_ph"] = -10000000
    kwargs["max_ph"] = -10000000

    for smi, protonated, deprotonated, category in smiles_groups:
        check_tests(smi, kwargs, [protonated], ["PROTONATED"])

    # Test phosphates separately
    for smi, protonated, mix, deprotonated, category in smiles_phosphates:
        check_tests(smi, kwargs, [protonated], ["PROTONATED"])


def test_very_basic(kwargs_default, smiles_groups, smiles_phosphates):
    kwargs = kwargs_default
    kwargs["min_ph"] = 10000000
    kwargs["max_ph"] = 10000000

    for smi, protonated, deprotonated, category in smiles_groups:
        check_tests(smi, kwargs, [deprotonated], ["DEPROTONATED"])

    for smi, protonated, mix, deprotonated, category in smiles_phosphates:
        check_tests(smi, kwargs, [deprotonated], ["DEPROTONATED"])


def test_avg_pka(
    kwargs_default,
    smiles_groups,
    smiles_phosphates,
    average_pkas_groups,
    average_pkas_phosphates,
):
    kwargs = kwargs_default

    for smi, protonated, deprotonated, category in smiles_groups:
        avg_pka = average_pkas_groups[category]

        kwargs["min_ph"] = avg_pka
        kwargs["max_ph"] = avg_pka

        check_tests(smi, kwargs, [protonated, deprotonated], ["BOTH"])

    for smi, protonated, mix, deprotonated, category in smiles_phosphates:
        kwargs["smiles"] = smi

        avg_pka = average_pkas_phosphates[category][0]
        kwargs["min_ph"] = avg_pka
        kwargs["max_ph"] = avg_pka

        check_tests(smi, kwargs, [mix, protonated], ["BOTH"])

        avg_pka = average_pkas_phosphates[category][1]
        kwargs["min_ph"] = avg_pka
        kwargs["max_ph"] = avg_pka

        check_tests(smi, kwargs, [mix, deprotonated], ["DEPROTONATED", "DEPROTONATED"])

        avg_pka = 0.5 * (
            average_pkas_phosphates[category][0] + average_pkas_phosphates[category][1]
        )
        kwargs["min_ph"] = avg_pka
        kwargs["max_ph"] = avg_pka
        kwargs["pka_precision"] = 5  # Should give all three

        check_tests(smi, kwargs, [mix, deprotonated, protonated], ["BOTH", "BOTH"])


def test_no_carbanion():
    smi = (
        "Cc1nc2cc(-c3[nH]c4cc5ccccc5c5c4c3CCN(C(=O)O)[C@@H]5O)cc3c(=O)[nH][nH]c(n1)c23"
    )
    output = list(protonate_smiles(smi))

    if "[C-]" in "".join(output).upper():
        msg = "Processing " + smi + " produced a molecule with a carbanion!"
        raise Exception(msg)
    else:
        print("(CORRECT) No carbanion: " + smi)


def test_max_variants():
    # Make sure max number of variants is limited (old bug).
    smi = "CCCC[C@@H](C(=O)N)NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](C(C)C)NC(=O)[C@@H](NC(=O)[C@H](Cc1c[nH]c2c1cccc2)NC(=O)[C@@H](NC(=O)[C@@H](Cc1ccc(cc1)O)N)CCC(=O)N)C)C)Cc1nc[nH]c1)Cc1ccccc1"
    output = list(protonate_smiles(smi))
    if len(output) != 128:
        msg = "Processing " + smi + " produced more than 128 variants!"
        raise Exception(msg)
    else:
        print("(CORRECT) Produced 128 variants: " + smi)


def test_atp_nad():
    # Make sure ATP and NAD work at different pHs (because can't test
    # Internal_phosphate_polyphos_chain and
    # Initial_phosphate_like_in_ATP_ADP with monoprotic examples.
    specific_examples = [
        [
            "O=P(O)(OP(O)(OP(O)(OCC1OC(C(C1O)O)N2C=NC3=C2N=CN=C3N)=O)=O)O",  # input, ATP
            (
                0.5,
                "[NH3+]c1[nH+]c[nH+]c2c1[nH+]cn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O",
            ),
            (
                1.0,
                "[NH3+]c1[nH+]c[nH+]c2c1[nH+]cn2C1OC(COP(=O)(O)OP(=O)([O-])OP(=O)(O)O)C(O)C1O",
            ),
            (
                2.6,
                "[NH3+]c1[nH+]c[nH+]c2c1[nH+]cn2C1OC(COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])O)C(O)C1O",
            ),
            (
                7.0,
                "Nc1ncnc2c1ncn2C1OC(COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])C(O)C1O",
            ),
        ],
        [
            "O=P(O)(OP(O)(OCC1C(O)C(O)C(N2C=NC3=C(N)N=CN=C32)O1)=O)OCC(O4)C(O)C(O)C4[N+]5=CC=CC(C(N)=O)=C5",  # input, NAD
            (
                0.5,
                "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c([NH3+])ncnc54)C(O)C3O)C(O)C2O)c1",
            ),
            (
                2.5,
                "NC(=O)c1ccc[n+](C2OC(COP(=O)([O-])OP(=O)([O-])OCC3OC(n4cnc5c([NH3+])ncnc54)C(O)C3O)C(O)C2O)c1",
            ),
            (
                7.4,
                "NC(=O)c1ccc[n+](C2OC(COP(=O)([O-])OP(=O)([O-])OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1",
            ),
        ],
    ]
    for example in specific_examples:
        smi = str(example[0])
        for ph, expected_output in example[1:]:
            ph = float(ph)
            output = list(
                protonate_smiles(smi, min_ph=ph, max_ph=ph, pka_precision=0.0)
            )
            if output[0].strip() == expected_output:
                print(
                    "(CORRECT) "
                    + smi
                    + " at pH "
                    + str(ph)
                    + " is "
                    + output[0].strip()
                )
            else:
                msg = (
                    smi
                    + " at pH "
                    + str(ph)
                    + " should be "
                    + expected_output
                    + ", but it is "
                    + output[0].strip()
                )
                raise Exception(msg)

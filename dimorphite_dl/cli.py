import argparse

from loguru import logger

from dimorphite_dl import __version__
from dimorphite_dl.protonate import protonate_smiles


def run_cli() -> None:
    """The main definition run when you call the script from the commandline."""
    parser = argparse.ArgumentParser(description=f"dimorphite_dl v{__version__}")
    parser.add_argument(
        "--ph_min",
        metavar="MIN",
        type=float,
        default=6.4,
        help="Minimum pH to consider (default: 6.4)",
    )
    parser.add_argument(
        "--ph_max",
        metavar="MAX",
        type=float,
        default=8.4,
        help="Maximum pH to consider",
    )
    parser.add_argument(
        "--precision",
        metavar="PRE",
        type=float,
        default=1.0,
        help="pKa precision factor (i.e., number of standard devations)",
    )
    parser.add_argument(
        "--output_file",
        metavar="FILE",
        type=str,
        help="Output file to write protonated SMILES (optional)",
    )
    parser.add_argument(
        "--max_variants",
        metavar="MXV",
        type=int,
        default=128,
        help="Limit number of variants per input compound (default: 128)",
    )
    parser.add_argument(
        "--label_states",
        action="store_true",
        help="label protonated SMILES with target state "
        + '(i.e., "DEPROTONATED", "PROTONATED", or "BOTH").',
    )
    parser.add_argument(
        "--silent", action="store_true", help="do not print any messages to the screen"
    )
    parser.add_argument(
        "smiles", metavar="SMI", type=str, help="SMILES or path to SMILES to protonate"
    )

    args = parser.parse_args()

    if args.output_file is not None:
        logger.info("Writing smiles to {}", args.output_file)
        f = open(args.output_file, "w", encoding="utf-8")

    for smiles_protonated in protonate_smiles(
        smiles_input=args.smiles,
        ph_min=args.ph_min,
        ph_max=args.ph_max,
        precision=args.precision,
        label_states=args.label_states,
        max_variants=args.max_variants,
    ):
        if args.output_file is not None:
            f.write(smiles_protonated + "\n")
        else:
            print(smiles_protonated)

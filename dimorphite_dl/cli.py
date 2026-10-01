import argparse
import os

from loguru import logger

from dimorphite_dl import __version__, enable_logging, protonate_smiles

LOG_LEVEL_TO_INT = {"debug": 10, "info": 20, "warning": 30, "error": 40, "critical": 50}


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
        help="Maximum pH to consider (default: 8.4)",
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
        "--log_level",
        choices=["none", "debug", "info", "warning", "error", "critical"],
        default="none",
        help="Enable and set logging level. Defaults to none (i.e., no logging)",
    )
    parser.add_argument(
        "smiles", metavar="SMI", type=str, help="SMILES or path to SMILES to protonate"
    )

    args = parser.parse_args()
    if args.log_level != "none":
        enable_logging(LOG_LEVEL_TO_INT[args.log_level])

    # Writing would replace the input with its own protonated forms.
    # samefile also catches symlinks and hard links.
    if (
        args.output_file is not None
        and os.path.exists(args.smiles)
        and os.path.exists(args.output_file)
        and os.path.samefile(args.smiles, args.output_file)
    ):
        parser.error(f"--output_file ({args.output_file}) is the input file")

    # Protonated before the output file is opened, so that invalid arguments
    # fail without truncating an existing file.
    smiles_protonated_all = protonate_smiles(
        smiles_input=args.smiles,
        ph_min=args.ph_min,
        ph_max=args.ph_max,
        precision=args.precision,
        label_states=args.label_states,
        max_variants=args.max_variants,
    )

    if args.output_file is not None:
        logger.info("Writing smiles to {}", args.output_file)
        with open(args.output_file, "w", encoding="utf-8") as f:
            for smiles_protonated in smiles_protonated_all:
                f.write(smiles_protonated + "\n")
    else:
        for smiles_protonated in smiles_protonated_all:
            print(smiles_protonated)

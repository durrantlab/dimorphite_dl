from typing import Any

import argparse
import sys
from io import StringIO

from rdkit import Chem

from dimorphite_dl.io import LoadSMIFile
from dimorphite_dl.mol import Protonate


def print_header():
    """Prints out header information."""
    # Always let the user know a help file is available.
    print("\nFor help, use: dimorphite_dl --help")

    # And always report citation information.
    print("\nIf you use Dimorphite-DL in your research, please cite:")
    print("Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An")
    print(
        "open-source program for enumerating the ionization states of drug-like small"
    )
    print("molecules. J Cheminform 11:14. doi:10.1186/s13321-019-0336-9.\n")


def main(params: None | dict[str, Any] = None) -> list[str] | None:
    """The main definition run when you call the script from the commandline.

    Args:
        params: The parameters to use. Entirely optional. If absent,
            defaults to None, in which case arguments will be taken from
            those given at the command line.
    Returns:
         Returns a list of the SMILES strings return_as_list parameter is
             True. Otherwise, returns None.
    """

    parser = ArgParseFuncs.get_args()
    args = vars(parser.parse_args())
    args = ArgParseFuncs.clean_args(args)

    if not args["silent"]:
        print_header()

    # Add in any parameters in params.
    if params is not None:
        for k, v in params.items():
            args[k] = v

    # If being run from the command line, print out all parameters.
    if __name__ == "__main__":
        if not args["silent"]:
            print("\nPARAMETERS:\n")
            for k in sorted(args.keys()):
                print(k.rjust(13) + ": " + str(args[k]))
            print("")

    # Run protonation
    if "output_file" in args and args["output_file"] is not None:
        # An output file was specified, so write to that.
        with open(args["output_file"], "w") as file:
            for protonated_smi in Protonate(args):
                file.write(protonated_smi + "\n")
    elif "return_as_list" in args and args["return_as_list"] == True:
        return list(Protonate(args))
    else:
        # No output file specified. Just print it to the screen.
        for protonated_smi in Protonate(args):
            print(protonated_smi)


class MyParser(argparse.ArgumentParser):
    """Overwrite default parse so it displays help file on error. See
    https://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu"""

    def error(self, message: str) -> None:
        """Overwrites the default error message.

        Args:
            message: The default error message.
        """

        self.print_help()
        msg = "ERROR: %s\n\n" % message
        print(msg)
        raise Exception(msg)

    def print_help(self, file=None) -> None:
        """Overwrite the default print_help function

        Args:
            file: Output file, defaults to None
        """

        print("")

        if file is None:
            file = sys.stdout
        self._print_message(self.format_help(), file)
        print(
            """
examples:
  dimorphite_dl --smiles_file sample_molecules.smi
  dimorphite_dl --smiles "CCC(=O)O" --min_ph -3.0 --max_ph -2.0
  dimorphite_dl --smiles "CCCN" --min_ph -3.0 --max_ph -2.0 --output_file output.smi
  dimorphite_dl --smiles_file sample_molecules.smi --pka_precision 2.0 --label_states"""
        )
        print("")


class ArgParseFuncs:
    """A namespace for storing functions that are useful for processing
    command-line arguments. To keep things organized."""

    @staticmethod
    def get_args():
        """Gets the arguments from the command line.

        Returns:
            A parser object.
        """

        parser = MyParser(
            description="Dimorphite 1.2.4: Creates models of "
            + "appropriately protonated small moleucles. "
            + "Apache 2.0 License. Copyright 2020 Jacob D. "
            + "Durrant."
        )
        parser.add_argument(
            "--min_ph",
            metavar="MIN",
            type=float,
            default=6.4,
            help="minimum pH to consider (default: 6.4)",
        )
        parser.add_argument(
            "--max_ph",
            metavar="MAX",
            type=float,
            default=8.4,
            help="maximum pH to consider (default: 8.4)",
        )
        parser.add_argument(
            "--pka_precision",
            metavar="PRE",
            type=float,
            default=1.0,
            help="pKa precision factor (number of standard devations, default: 1.0)",
        )
        parser.add_argument(
            "--smiles", metavar="SMI", type=str, help="SMILES string to protonate"
        )
        parser.add_argument(
            "--smiles_file",
            metavar="FILE",
            type=str,
            help="file that contains SMILES strings to protonate",
        )
        parser.add_argument(
            "--output_file",
            metavar="FILE",
            type=str,
            help="output file to write protonated SMILES (optional)",
        )
        parser.add_argument(
            "--max_variants",
            metavar="MXV",
            type=int,
            default=128,
            help="limit number of variants per input compound (default: 128)",
        )
        parser.add_argument(
            "--label_states",
            action="store_true",
            help="label protonated SMILES with target state "
            + '(i.e., "DEPROTONATED", "PROTONATED", or "BOTH").',
        )
        parser.add_argument(
            "--silent",
            action="store_true",
            help="do not print any messages to the screen",
        )

        return parser

    @staticmethod
    def clean_args(args):
        """Cleans and normalizes input parameters

        Args:
            args: A dictionary containing the arguments.
        """

        defaults = {
            "min_ph": 6.4,
            "max_ph": 8.4,
            "pka_precision": 1.0,
            "label_states": False,
            "max_variants": 128,
        }

        for key in defaults:
            if key not in args:
                args[key] = defaults[key]

        keys = list(args.keys())
        for key in keys:
            if args[key] is None:
                del args[key]

        if "smiles" not in args and "smiles_file" not in args:
            msg = "Error: No SMILES in params. Use the -h parameter for help."
            print(msg)
            raise Exception(msg)

        # If the user provides a smiles string, turn it into a file-like StringIO
        # object.
        if "smiles" in args:
            if isinstance(args["smiles"], str):
                args["smiles_file"] = StringIO(args["smiles"])

        args["smiles_and_data"] = LoadSMIFile(args["smiles_file"], args)

        return args


def run(**kwargs):
    """A helpful, importable function for those who want to call Dimorphite-DL
    from another Python script rather than the command line. Note that this
    function accepts keyword arguments that match the command-line parameters
    exactly. If you want to pass and return a list of RDKit Mol objects, import
    run_with_mol_list() instead.

    Args:
        **kwargs: For a complete description, run dimorphite_dl from the
            command line with the -h option.
    """
    main(kwargs)


def run_with_mol_list(mol_lst, **kwargs):
    """A helpful, importable function for those who want to call Dimorphite-DL
    from another Python script rather than the command line. Note that this
    function is for passing Dimorphite-DL a list of RDKit Mol objects, together
    with command-line parameters. If you want to use only the same parameters
    that you would use from the command line, import run() instead.

    Args:
        mol_lst: A list of rdkit.Chem.rdchem.Mol objects.

    Returns:
        A list of properly protonated rdkit.Chem.rdchem.Mol objects.
    """

    # Do a quick check to make sure the user input makes sense.
    for bad_arg in ["smiles", "smiles_file", "output_file", "test"]:
        if bad_arg in kwargs:
            msg = (
                "You're using Dimorphite-DL's run_with_mol_list(mol_lst, "
                + '**kwargs) function, but you also passed the "'
                + bad_arg
                + '" argument. Did you mean to use the '
                + "run(**kwargs) function instead?"
            )
            print(msg)
            raise Exception(msg)

    # Set the return_as_list flag so main() will return the protonated smiles
    # as a list.
    kwargs["return_as_list"] = True

    # Having reviewed the code, it will be very difficult to rewrite it so
    # that a list of Mol objects can be used directly. Intead, convert this
    # list of mols to smiles and pass that. Not efficient, but it will work.
    protonated_smiles_and_props = []
    for m in mol_lst:
        props = m.GetPropsAsDict()
        kwargs["smiles"] = Chem.MolToSmiles(m, isomericSmiles=True)
        protonated_smiles_and_props.extend(
            [(s.split("\t")[0], props) for s in main(kwargs)]
        )

    # Now convert the list of protonated smiles strings back to RDKit Mol
    # objects. Also, add back in the properties from the original mol objects.
    mols = []
    for s, props in protonated_smiles_and_props:
        m = Chem.MolFromSmiles(s)
        if m:
            for prop, val in props.items():
                if type(val) is int:
                    m.SetIntProp(prop, val)
                elif type(val) is float:
                    m.SetDoubleProp(prop, val)
                elif type(val) is bool:
                    m.SetBoolProp(prop, val)
                else:
                    m.SetProp(prop, str(val))
            mols.append(m)
        else:
            print(
                "WARNING: Could not process molecule with SMILES string "
                + s
                + " and properties "
                + str(props)
            )

    return mols

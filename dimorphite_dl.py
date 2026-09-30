# Copyright 2018 Jacob D. Durrant
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
This script identifies and enumerates the possible protonation sites of SMILES
strings.
"""

from __future__ import print_function
import copy
import os
import argparse
import sys

try:
    # Only needed by type checkers; the annotations below are type comments so
    # the module still parses under Python 2.
    from typing import List, Tuple
except ImportError:
    pass

try:
    # Python2
    from StringIO import StringIO
except ImportError:
    # Python3
    from io import StringIO

try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    raise ImportError("Dimorphite-DL requires RDKit. See https://www.rdkit.org/")

# What RDKit raises for a molecule it cannot handle: sanitization failures
# derive from ValueError, and C++ invariant violations surface as
# RuntimeError. Those mean "skip this input". Anything else is a bug in this
# module and should propagate rather than be reported as a bad SMILES.
RDKIT_ERRORS = (ValueError, RuntimeError)

# Caps the protonation states kept per input molecule. Each BOTH site doubles
# the count, so without a cap a molecule with 20 such sites would need about a
# million states in memory.
DEFAULT_MAX_VARIANTS = 128


def main(params=None):
    """The main definition run when you call the script from the commandline.

    :param params: The parameters to use. Entirely optional. If absent,
                   defaults to None, in which case argments will be taken from
                   those given at the command line.
    :param params: dict, optional
    :return: Returns a list of the SMILES strings return_as_list parameter is
             True. Otherwise, returns None.
    """

    parser = ArgParseFuncs.get_args()
    # Library callers pass params, so the host process's sys.argv (e.g., a
    # Jupyter kernel's "-f kernel.json") must not be parsed.
    args = vars(parser.parse_args([] if params is not None else None))

    # Add in any parameters in params.
    if params is not None:
        for k, v in params.items():
            args[k] = v

    # If being run from the command line, print out all parameters.
    if __name__ == "__main__":
        print("\nPARAMETERS:\n", file=sys.stderr)
        for k in sorted(args.keys()):
            print(k.rjust(13) + ": " + str(args[k]), file=sys.stderr)
        print("", file=sys.stderr)

    # Built before any output file is opened, so that invalid arguments raise
    # without truncating an existing file.
    protonator = Protonate(args)

    # Run protonation
    if "output_file" in args and args["output_file"] is not None:
        # An output file was specified, so write to that.
        with open(args["output_file"], "w") as file:
            for protonated_smi in protonator:
                file.write(protonated_smi + "\n")
    elif "return_as_list" in args and args["return_as_list"] == True:
        return list(protonator)
    else:
        # No output file specified. Just print it to the screen.
        for protonated_smi in protonator:
            print(protonated_smi)


class MyParser(argparse.ArgumentParser):
    """Overwrite default parse so it displays help file on error. See
    https://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    """

    def error(self, message):
        """Overwrites the default error message.

        :param message: The default error message.
        """

        self.print_help()
        msg = "ERROR: %s\n\n" % message
        print(msg)
        raise Exception(msg)

    def print_help(self, file=None):
        """Overwrite the default print_help function

        :param file: Output file, defaults to None
        """

        print("")

        if file is None:
            file = sys.stdout
        self._print_message(self.format_help(), file)
        print(
            """
examples:
  python dimorphite_dl.py --smiles_file sample_molecules.smi
  python dimorphite_dl.py --smiles "CCC(=O)O" --min_ph -3.0 --max_ph -2.0
  python dimorphite_dl.py --smiles "CCCN" --min_ph -3.0 --max_ph -2.0 --output_file output.smi
  python dimorphite_dl.py --smiles_file sample_molecules.smi --pka_precision 2.0 --label_states"""
        )
        print("")


class ArgParseFuncs:
    """A namespace for storing functions that are useful for processing
    command-line arguments. To keep things organized."""

    @staticmethod
    def get_args():
        """Gets the arguments from the command line.

        :return: A parser object.
        """

        parser = MyParser(
            description="Dimorphite 1.1: Creates models of "
            + "appropriately protonated small moleucles. "
            + "Apache 2.0 License. Copyright 2018 Jacob D. "
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
            "--max_variants",
            metavar="MXV",
            type=int,
            default=DEFAULT_MAX_VARIANTS,
            help="limit number of variants per input compound (default: %d)"
            % DEFAULT_MAX_VARIANTS,
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
            "--label_states",
            action="store_true",
            help="label protonated SMILES with target state "
            + '(i.e., "DEPROTONATED", "PROTONATED", or "BOTH").',
        )
        return parser

    @staticmethod
    def clean_args(args):
        """Cleans and normalizes input parameters

        :param args: A dictionary containing the arguments.
        :type args: dict
        :raises Exception: No SMILES in params.
        """

        defaults = {
            "min_ph": 6.4,
            "max_ph": 8.4,
            "pka_precision": 1.0,
            "max_variants": DEFAULT_MAX_VARIANTS,
            "label_states": False,
        }

        for key in defaults:
            if key not in args:
                args[key] = defaults[key]

        keys = list(args.keys())
        for key in keys:
            if args[key] is None:
                del args[key]

        if not "smiles" in args and not "smiles_file" in args:
            msg = "Error: No SMILES in params. Use the -h parameter for help."
            print(msg)
            raise Exception(msg)

        if "smiles" in args and "smiles_file" in args:
            msg = (
                "Error: Both smiles and smiles_file were given. Only one input "
                + "can be used."
            )
            UtilFuncs.eprint(msg)
            raise ValueError(msg)

        if args["min_ph"] > args["max_ph"]:
            msg = "Error: min_ph (%s) is greater than max_ph (%s)." % (
                args["min_ph"],
                args["max_ph"],
            )
            UtilFuncs.eprint(msg)
            raise ValueError(msg)

        if args["pka_precision"] < 0:
            msg = "Error: pka_precision (%s) must not be negative." % (
                args["pka_precision"]
            )
            UtilFuncs.eprint(msg)
            raise ValueError(msg)

        if args["max_variants"] < 1:
            msg = "Error: max_variants (%s) must be at least 1." % (
                args["max_variants"]
            )
            UtilFuncs.eprint(msg)
            raise ValueError(msg)

        # If the user provides a smiles string, turn it into a file-like StringIO
        # object.
        if "smiles" in args:
            if isinstance(args["smiles"], str):
                args["smiles_file"] = StringIO(args["smiles"])

        args["smiles_and_data"] = LoadSMIFile(args["smiles_file"])

        return args


class UtilFuncs:
    """A namespace to store functions for manipulating mol objects. To keep
    things organized."""

    @staticmethod
    def neutralize_mol(mol):
        """All molecules should be neuralized to the extent possible. The user
        should not be allowed to specify the valence of the atoms in most cases.

        :param rdkit.Chem.rdchem.Mol mol: The rdkit Mol objet to be neutralized.
        :return: The neutralized Mol object.
        """

        # Get the reaction data
        rxn_data = [
            [
                "[Ov1-1:1]",
                "[Ov2+0:1]-[H]",
            ],  # To handle O- bonded to only one atom (add hydrogen).
            [
                "[#7v4+1:1]-[H]",
                "[#7v3+0:1]",
            ],  # To handle N+ bonded to a hydrogen (remove hydrogen).
            [
                "[Ov2-:1]",
                "[Ov2+0:1]",
            ],  # To handle O- bonded to two atoms. Should not be Negative.
            [
                "[#7v3+1:1]",
                "[#7v3+0:1]",
            ],  # To handle N+ bonded to three atoms. Should not be positive.
            [
                "[#7v2-1:1]",
                "[#7+0:1]-[H]",
            ],  # To handle N- Bonded to two atoms. Add hydrogen.
            # ['[N:1]=[N+0:2]=[N:3]-[H]', '[N:1]=[N+1:2]=[N+0:3]-[H]'],  # To
            # handle bad azide. Must be protonated. (Now handled elsewhere, before
            # SMILES converted to Mol object.)
            [
                "[H]-[N:1]-[N:2]#[N:3]",
                "[N:1]=[N+1:2]=[N:3]-[H]",
            ],  # To handle bad azide. R-N-N#N should be R-N=[N+]=N
        ]

        # Add substructures and reactions (initially none)
        for i, rxn_datum in enumerate(rxn_data):
            rxn_data[i].append(Chem.MolFromSmarts(rxn_datum[0]))
            rxn_data[i].append(None)

        # Add hydrogens (respects valence, so incomplete).
        # Chem.calcImplicitValence(mol)
        mol.UpdatePropertyCache(strict=False)
        mol = Chem.AddHs(mol)

        while True:  # Keep going until all these issues have been resolved.
            current_rxn = None  # The reaction to perform.
            current_rxn_str = None

            for i, rxn_datum in enumerate(rxn_data):
                (
                    reactant_smarts,
                    product_smarts,
                    substruct_match_mol,
                    rxn_placeholder,
                ) = rxn_datum
                if mol.HasSubstructMatch(substruct_match_mol):
                    if rxn_placeholder is None:
                        current_rxn_str = reactant_smarts + ">>" + product_smarts
                        current_rxn = AllChem.ReactionFromSmarts(current_rxn_str)
                        rxn_data[i][3] = current_rxn  # Update the placeholder.
                    else:
                        current_rxn = rxn_data[i][3]
                    break

            # Perform the reaction if necessary
            if current_rxn is None:  # No reaction left, so break out of while loop.
                break
            else:
                mol = current_rxn.RunReactants((mol,))[0][0]
                mol.UpdatePropertyCache(strict=False)  # Update valences

        # The mols have been altered from the reactions described above, we need
        # to resanitize them. Make sure aromatic rings are shown as such This
        # catches all RDKit Errors. without the catchError and sanitizeOps the
        # Chem.SanitizeMol can crash the program.
        sanitize_string = Chem.SanitizeMol(
            mol,
            sanitizeOps=rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
            catchErrors=True,
        )

        return mol if sanitize_string.name == "SANITIZE_NONE" else None

    @staticmethod
    def convert_smiles_str_to_mol(smiles_str):
        """Given a SMILES string, check that it is actually a string and not a
        None. Then try to convert it to an RDKit Mol Object.

        :param string smiles_str: The SMILES string.
        :return: A rdkit.Chem.rdchem.Mol object, or None if it is the wrong type or
            if it fails to convert to a Mol Obj
        """

        if smiles_str is None or type(smiles_str) is not str:
            return None

        # Check that there are no type errors, ie Nones or non-string
        # A non-string type will cause RDKit to hard crash
        try:
            # Try to fix azides here. They are just tricky to deal with.
            smiles_str = smiles_str.replace("N=N=N", "N=[N+]=N")
            smiles_str = smiles_str.replace("NN#N", "N=[N+]=N")
            mol = Chem.MolFromSmiles(smiles_str)
        except RDKIT_ERRORS:
            return None

        # Check that there are None type errors Chem.MolFromSmiles has sanitize on
        # which means if there is even a small error in the SMILES (kekulize,
        # nitrogen charge...) then mol=None. ie.
        # Chem.MolFromSmiles("C[N]=[N]=[N]") = None this is an example of an
        # nitrogen charge error. It is cased in a try statement to be overly
        # cautious.

        return None if mol is None else mol

    @staticmethod
    def eprint(*args, **kwargs):
        """Error messages should be printed to STDERR. See
        https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python"""

        print(*args, file=sys.stderr, **kwargs)


class LoadSMIFile(object):
    """A generator class for loading in the SMILES strings from a file, one at
    a time."""

    def __init__(self, filename):
        """Initializes this class.

        :param filename: The filename or file object (i.e., StringIO).
        :type filename: str or StringIO
        """

        if type(filename) is str:
            # It's a filename
            self.f = open(filename, "r")
        else:
            # It's a file object (i.e., StringIO)
            self.f = filename

    def __iter__(self):
        """Returns this generator object.

        :return: This generator object.
        :rtype: LoadSMIFile
        """

        return self

    def __next__(self):
        """Ensure Python3 compatibility.

        :return: A dict, where the "smiles" key contains the canonical SMILES
                 string and the "data" key contains the remaining information
                 (e.g., the molecule name).
        :rtype: dict
        """

        return self.next()

    def next(self):
        """Get the data associated with the next line.

        :raises StopIteration: If there are no more lines left iin the file.
        :return: A dict, where the "smiles" key contains the canonical SMILES
                 string and the "data" key contains the remaining information
                 (e.g., the molecule name).
        :rtype: dict
        """

        # A loop rather than recursion, so a long run of blank or bad lines
        # cannot exhaust the stack.
        while True:
            line = self.f.readline()

            if line == "":
                # EOF
                self.f.close()
                raise StopIteration()

            # Divide line into smi and data
            splits = line.split()
            if len(splits) == 0:
                # Blank line? Go to next one.
                continue

            # Generate mol object
            smiles_str = splits[0]

            # Convert from SMILES string to RDKIT Mol. This series of tests is
            # to make sure the SMILES string is properly formed and to get it
            # into a canonical form. Filter if failed.
            mol = UtilFuncs.convert_smiles_str_to_mol(smiles_str)
            if mol is None:
                UtilFuncs.eprint(
                    "WARNING: Skipping poorly formed SMILES string: " + line
                )
                continue

            # Handle nuetralizing the molecules. Filter if failed.
            mol = UtilFuncs.neutralize_mol(mol)
            if mol is None:
                UtilFuncs.eprint(
                    "WARNING: Skipping poorly formed SMILES string: " + line
                )
                continue

            # Remove the hydrogens.
            try:
                mol = Chem.RemoveHs(mol)
            except RDKIT_ERRORS:
                UtilFuncs.eprint(
                    "WARNING: Skipping poorly formed SMILES string: " + line
                )
                continue

            if mol is None:
                UtilFuncs.eprint(
                    "WARNING: Skipping poorly formed SMILES string: " + line
                )
                continue

            # Regenerate the smiles string (to standardize).
            new_mol_string = Chem.MolToSmiles(mol, isomericSmiles=True)

            return {"smiles": new_mol_string, "data": splits[1:]}


class Protonate(object):
    """A generator class for protonating SMILES strings, one at a time."""

    def __init__(self, args):
        """Initialize the generator.

        :param args: A dictionary containing the arguments.
        :type args: dict
        """

        # Make the args an object variable variable.
        self.args = args

        # A list to store the protonated SMILES strings associated with a
        # single input model.
        self.cur_prot_SMI = []

        # Clean and normalize the args
        self.args = ArgParseFuncs.clean_args(args)

        # Load the substructures that can be protonated.
        self.subs = ProtSubstructFuncs.load_protonation_substructs_calc_state_for_ph(
            self.args["min_ph"], self.args["max_ph"], self.args["pka_precision"]
        )

    def __iter__(self):
        """Returns this generator object.

        :return: This generator object.
        :rtype: Protonate
        """

        return self

    def __next__(self):
        """Ensure Python3 compatibility.

        :return: A dict, where the "smiles" key contains the canonical SMILES
                 string and the "data" key contains the remaining information
                 (e.g., the molecule name).
        :rtype: dict
        """

        return self.next()

    def next(self):
        """Get the next protonated SMILES string.

        :raises StopIteration: If there are no more lines left iin the file.
        :return: TODO A dict, where the "smiles" key contains the canonical SMILES
                 string and the "data" key contains the remaining information
                 (e.g., the molecule name).
        :rtype: dict
        """

        # If there are any SMILES strings in self.cur_prot_SMI, just return
        # the first one and update the list to include only the remaining.
        if len(self.cur_prot_SMI) > 0:
            first, self.cur_prot_SMI = self.cur_prot_SMI[0], self.cur_prot_SMI[1:]
            return first

        # self.cur_prot_SMI is empty, so try to add more to it.

        # Get the next SMILES string from the input file.
        try:
            smile_and_datum = self.args["smiles_and_data"].next()
        except StopIteration:
            # There are no more input smiles strings...
            raise StopIteration()

        smi = smile_and_datum["smiles"]
        data = smile_and_datum["data"]  # Everything on SMILES line but the
        # SMILES string itself (e.g., the
        # molecule name).

        # Collect the data associated with this smiles (e.g., the molecule
        # name).
        tag = " ".join(data)

        # sites is a list of (atom index, "PROTONATED|DEPROTONATED|BOTH").
        # Note that the second entry indicates what state the site SHOULD be
        # in (not the one it IS in per the SMILES string). It's calculated
        # based on the probablistic distributions obtained during training.
        sites = ProtSubstructFuncs.get_prot_sites_and_target_states(smi, self.subs)

        # This is the same parse the site indices came from, so they line up
        # with its atom ordering.
        mol = UtilFuncs.convert_smiles_str_to_mol(smi)
        if mol is None:
            new_smis = [smi]
        else:
            # Not re-parsed on purpose: some states, such as a protonated
            # pyridone's O=c...[nH+], cannot be kekulized but are still
            # reported.
            new_smis = [
                smi_and_mol[0]
                for smi_and_mol in ProtSubstructFuncs.protonate_sites(
                    mol, sites, self.args["max_variants"]
                )
            ]

        # If the user wants to see the target states, add those
        # to the ends of each line.
        if self.args["label_states"]:
            states = "\t".join([x[1] for x in sites])
            new_lines = [x + "\t" + tag + "\t" + states for x in new_smis]
        else:
            new_lines = [x + "\t" + tag for x in new_smis]

        self.cur_prot_SMI = new_lines
        return self.next()


class ProtSubstructFuncs:
    """A namespace to store functions for loading the substructures that can
    be protonated. To keep things organized."""

    @staticmethod
    def load_protonation_substructs_calc_state_for_ph(
        min_ph=6.4, max_ph=8.4, pka_std_range=1
    ):
        """A pre-calculated list of R-groups with protonation sites, with their
        likely pKa bins.

        :param float min_ph:  The lower bound on the pH range, defaults to 6.4.
        :param float max_ph:  The upper bound on the pH range, defaults to 8.4.
        :param pka_std_range: Basically the precision (stdev from predicted pKa to
                            consider), defaults to 1.
        :return: A dict of the protonation substructions for the specified pH
                range.
        """

        subs = []
        pwd = os.path.dirname(os.path.realpath(__file__))

        site_structures_file = "{}/{}".format(pwd, "site_substructures.smarts")
        with open(site_structures_file, "r") as substruct:
            for line in substruct:
                line = line.strip()
                sub = {}
                if line != "":
                    splits = line.split()
                    sub["name"] = splits[0]
                    sub["smart"] = splits[1]
                    sub["mol"] = Chem.MolFromSmarts(sub["smart"])

                    # NEED TO DIVIDE THIS BY 3s
                    pka_ranges = [
                        splits[i : i + 3] for i in range(2, len(splits) - 1, 3)
                    ]

                    prot = []
                    for pka_range in pka_ranges:
                        site = pka_range[0]
                        std = float(pka_range[2]) * pka_std_range
                        mean = float(pka_range[1])
                        protonation_state = ProtSubstructFuncs.define_protonation_state(
                            mean, std, min_ph, max_ph
                        )

                        prot.append([site, protonation_state])

                    sub["prot_states_for_pH"] = prot
                    subs.append(sub)
        return subs

    @staticmethod
    def define_protonation_state(mean, std, min_ph, max_ph):
        """Updates the substructure definitions to include the protonation state
        based on the user-given pH range. The size of the pKa range is also based
        on the number of standard deviations to be considered by the user param.

        :param float mean:   The mean pKa.
        :param float std:    The precision (stdev).
        :param float min_ph: The min pH of the range.
        :param float max_ph: The max pH of the range.
        :return: A string describing the protonation state.
        :raises ValueError: If min_ph > max_ph or std < 0. Either inverts an
            interval, and the overlap test below then silently gives the wrong
            state.
        """

        if min_ph > max_ph:
            raise ValueError(
                "min_ph (%s) is greater than max_ph (%s)." % (min_ph, max_ph)
            )
        if std < 0:
            raise ValueError("The pKa standard deviation (%s) is negative." % std)

        min_pka = mean - std
        max_pka = mean + std

        # This needs to be reassigned, and 'ERROR' should never make it past the
        # next set of checks.
        if min_pka <= max_ph and min_ph <= max_pka:
            protonation_state = "BOTH"
        elif mean > max_ph:
            protonation_state = "PROTONATED"
        else:
            protonation_state = "DEPROTONATED"

        return protonation_state

    @staticmethod
    def get_prot_sites_and_target_states(smi, subs):
        """For a single molecule, find all possible matches in the protonation
        R-group list, subs. Items that are higher on the list will be matched
        first, to the exclusion of later items.

        :param string smi: A SMILES string.
        :param list subs: Substructure information.
        :return: A list of protonation sites and their pKa bin. ('PROTONATED',
            'BOTH', or  'DEPROTONATED')
        """

        # Convert the Smiles string (smi) to an RDKit Mol Obj
        mol = UtilFuncs.convert_smiles_str_to_mol(smi)

        # Check Conversion worked
        if mol is None:
            UtilFuncs.eprint("ERROR:   ", smi)
            return []

        # Try to Add hydrogens. if failed return []
        try:
            mol = Chem.AddHs(mol)
        except RDKIT_ERRORS:
            UtilFuncs.eprint("ERROR:   ", smi)
            return []

        # Check adding Hs worked
        if mol is None:
            UtilFuncs.eprint("ERROR:   ", smi)
            return []

        ProtectUnprotectFuncs.unprotect_molecule(mol)
        protonation_sites = []

        for item in subs:
            smart = item["mol"]
            if mol.HasSubstructMatch(smart):
                matches = ProtectUnprotectFuncs.get_unprotected_matches(mol, smart)
                prot = item["prot_states_for_pH"]
                context_idxs = []
                for match in matches:
                    # Re-checked per match so that two matches of one pattern
                    # cannot claim the same site, as happens with the three
                    # interchangeable OH groups of phosphoric acid.
                    if not ProtectUnprotectFuncs.is_match_unprotected(mol, match):
                        continue

                    # We want to move the site from being relative to the
                    # substructure, to the index on the main molecule.
                    site_idxs = []
                    for site in prot:
                        proton = int(site[0])
                        category = site[1]
                        site_idxs.append(match[proton])
                        new_site = (match[proton], category, item["name"])

                        if not new_site in protonation_sites:
                            # Because sites must be unique.
                            protonation_sites.append(new_site)

                    h_idxs = [
                        neighbor.GetIdx()
                        for idx in site_idxs
                        for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors()
                        if neighbor.GetAtomicNum() == 1
                    ]
                    ProtectUnprotectFuncs.protect_molecule(mol, site_idxs + h_idxs)

                    context_idxs.extend(
                        idx
                        for idx in match
                        if idx not in site_idxs
                        and mol.GetAtomWithIdx(idx).GetAtomicNum() != 6
                    )

                # Heteroatom context (an amidine's second N, a carbonyl O, a
                # phosphate's ester O) is part of this functional group, so
                # later patterns must not claim it. It is locked only after
                # the whole pattern is done because one bridging O serves both
                # phosphates of a pyrophosphate. Carbon context is just the
                # attachment point and is never a site, so it stays free: the
                # amine of NCP(=O)(O)O needs the carbon the phosphonate used.
                ProtectUnprotectFuncs.protect_molecule(mol, context_idxs)

        return protonation_sites

    @staticmethod
    def protonate_sites(mol, sites, max_variants):
        # type: (Chem.Mol, List[Tuple[int, str, str]], int) -> List[Tuple[str, Chem.Mol]]
        """Enumerates the protonation states of every site on one molecule.

        Works on Mol objects end to end because the site indices refer to
        mol's atom ordering. Re-serializing to canonical SMILES between sites
        can reorder the atoms once a charge has been set. Duplicates are
        dropped and the list is capped after every site, not just at the end,
        so memory stays bounded when many sites are BOTH.

        Args:
            mol: The molecule the site indices were computed from.
            sites: (idx, target_prot_state, prot_site_name) tuples.
            max_variants: The most states to keep. Once reached, later sites
                are still assigned, but only on the states already kept, and
                a warning is printed.

        Returns:
            Unique (canonical SMILES, Mol) pairs, in enumeration order.
        """

        input_smi = Chem.MolToSmiles(mol, isomericSmiles=True)
        new_mols = [(input_smi, mol)]
        truncated = False
        for site in sites:
            # Note that new_mols is a growing list. This is how multiple
            # protonation sites are handled.
            protonated = ProtSubstructFuncs.protonate_site(
                [m for _, m in new_mols], site
            )

            # In some cases, the script might generate redundant molecules.
            # Phosphonates, when the pH is between the two pKa values and the
            # stdev value is big enough, for example, will generate two
            # identical BOTH states. Removing them before capping keeps them
            # from taking the place of distinct states.
            new_mols = []
            seen = set()
            for protonated_mol in protonated:
                new_smi = Chem.MolToSmiles(protonated_mol, isomericSmiles=True)
                if new_smi not in seen:
                    seen.add(new_smi)
                    new_mols.append((new_smi, protonated_mol))

            if len(new_mols) > max_variants:
                new_mols = new_mols[:max_variants]
                truncated = True

        if truncated:
            UtilFuncs.eprint(
                "WARNING: Limited number of variants to %d (see --max_variants): %s"
                % (max_variants, input_smi)
            )

        return new_mols

    @staticmethod
    def protonate_site(mols, site):
        # type: (List[Chem.Mol], Tuple[int, str, str]) -> List[Chem.Mol]
        """Given a list of Mol objects, we protonate the site.

        :param list mols:  The list of Mol objects.
        :param tuple site: Information about the protonation site.
                        (idx, target_prot_state, prot_site_name)
        :return: A list of the appropriately protonated Mol objects.
        """

        # Decouple the atom index and its target protonation state from the site
        # tuple
        idx, target_prot_state, prot_site_name = site

        # Initialize the output list
        output_mols = []

        state_to_charge = {"DEPROTONATED": [-1], "PROTONATED": [0], "BOTH": [-1, 0]}

        charges = state_to_charge[target_prot_state]

        # Now make the actual molecules match the target protonation state.
        output_mols = ProtSubstructFuncs.set_protonation_charge(
            mols, idx, charges, prot_site_name
        )

        return output_mols

    @staticmethod
    def set_protonation_charge(mols, idx, charges, prot_site_name):
        # type: (List[Chem.Mol], int, List[int], str) -> List[Chem.Mol]
        """Sets the atomic charge on a particular site for a set of Mol objects.

        :param list mols:             A list of the Mol objects. They are not
                                    modified.
        :param int idx:               The index of the atom to consider.
        :param list charges:          A list of the charges (ints) to assign at
                                    this site.
        :param string prot_site_name: The name of the protonation site.
        :return: A list of the processed Mol objects.
        """

        # Sets up the output list and the Nitrogen charge
        output = []

        for charge in charges:
            # The charge for Nitrogens is 1 higher than others (i.e., protonated
            # state is positively charged).
            nitro_charge = charge + 1

            # But there are a few nitrogen moieties where the acidic group is the
            # neutral one. Amides are a good example. I gave some thought re. how
            # to best flag these. I decided that those nitrogen-containing
            # moieties where the acidic group is neutral (rather than positively
            # charged) will have "*" in the name.
            if "*" in prot_site_name:
                nitro_charge = nitro_charge - 1  # Undo what was done previously.

            for parent_mol in mols:

                # Copy, because the parent is reused for the other charges.
                mol = Chem.Mol(parent_mol)

                atom = mol.GetAtomWithIdx(idx)

                # Assign the protonation charge, with special care for Nitrogens
                element = atom.GetAtomicNum()
                if element == 7:
                    target_charge = nitro_charge
                else:
                    target_charge = charge

                delta = target_charge - atom.GetFormalCharge()
                num_hs = atom.GetTotalNumHs()
                atom.SetFormalCharge(target_charge)

                # On deprotonation, pin the hydrogen count rather than letting
                # RDKit recompute it. Aromatic [nH] keeps its H as an explicit
                # count even though GetNoImplicit() is False, so neither
                # SetFormalCharge nor UpdatePropertyCache removes it, which
                # gives an invalid [nH-].
                if delta < 0:
                    atom.SetNumExplicitHs(max(0, num_hs + delta))
                    atom.SetNoImplicit(True)

                mol.UpdatePropertyCache(strict=False)

                output.append(mol)

        return output


class ProtectUnprotectFuncs:
    """A namespace for storing functions that are useful for protecting and
    unprotecting molecules. To keep things organized. We need to identify and
    mark groups that have been matched with a substructure."""

    @staticmethod
    def unprotect_molecule(mol):
        """Sets the protected property on all atoms to 0. This also creates the
        property for new molecules.

        :param rdkit.Chem.rdchem.Mol mol: The rdkit Mol object.
        :type mol: The rdkit Mol object with atoms unprotected.
        """

        for atom in mol.GetAtoms():
            atom.SetProp("_protected", "0")

    @staticmethod
    def protect_molecule(mol, match):
        """Given a 'match', a list of molecules idx's, we set the protected status
        of each atom to 1. This will prevent any matches using that atom in the
        future.

        :param rdkit.Chem.rdchem.Mol mol: The rdkit Mol object to protect.
        :param list match: A list of molecule idx's.
        """

        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            atom.SetProp("_protected", "1")

    @staticmethod
    def get_unprotected_matches(mol, substruct):
        """Finds substructure matches with atoms that have not been protected.
        Returns list of matches, each match a list of atom idxs.

        :param rdkit.Chem.rdchem.Mol mol: The Mol object to consider.
        :param string substruct: The SMARTS string of the substructure ot match.
        :return: A list of the matches. Each match is itself a list of atom idxs.
        """

        matches = mol.GetSubstructMatches(substruct)
        unprotected_matches = []
        for match in matches:
            if ProtectUnprotectFuncs.is_match_unprotected(mol, match):
                unprotected_matches.append(match)
        return unprotected_matches

    @staticmethod
    def is_match_unprotected(mol, match):
        """Checks a molecule to see if the substructure match contains any
        protected atoms.

        :param rdkit.Chem.rdchem.Mol mol: The Mol object to check.
        :param list match: The match to check.
        :return: A boolean, whether the match is present or not.
        """

        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            protected = atom.GetProp("_protected")
            if protected == "1":
                return False
        return True


def run(**kwargs):
    """A helpful, importable function for those who want to call Dimorphite-DL
    from another Python script rather than the command line. Note that this
    function accepts keyword arguments that match the command-line parameters
    exactly. If you want to pass and return a list of RDKit Mol objects, import
    run_with_mol_list() instead.

    :param **kwargs: For a complete description, run dimorphite_dl.py from the
        command line with the -h option.
    :type kwargs: dict
    :return: The list of protonated SMILES lines if return_as_list is True.
        Otherwise, None.
    :rtype: list or None
    """

    # Run the main function with the specified arguments.
    return main(kwargs)


def run_with_mol_list(mol_lst, **kwargs):
    """A helpful, importable function for those who want to call Dimorphite-DL
    from another Python script rather than the command line. Note that this
    function is for passing Dimorphite-DL a list of RDKit Mol objects, together
    with command-line parameters. If you want to use only the same parameters
    that you would use from the command line, import run() instead.

    :param mol_lst: A list of rdkit.Chem.rdchem.Mol objects.
    :type mol_lst: list
    :raises Exception: If the **kwargs includes "smiles", "smiles_file", or
                       "output_file" parameters.
    :return: A list of properly protonated rdkit.Chem.rdchem.Mol objects.
    :rtype: list
    """

    # Do a quick check to make sure the user input makes sense.
    for bad_arg in ["smiles", "smiles_file", "output_file"]:
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
    # list of mols to smiles and pass that. Not memory efficient, but it will
    # work.
    protonated_smiles = []
    for m in mol_lst:
        smiles = Chem.MolToSmiles(m, isomericSmiles=True)
        kwargs["smiles"] = smiles
        protonated_smiles.extend([s.split("\t")[0] for s in main(kwargs)])

    # Now convert the list of protonated smiles strings back to RDKit Mol
    # objects.
    mols = [Chem.MolFromSmiles(s) for s in protonated_smiles]

    return mols


if __name__ == "__main__":
    # Printed only when run as a script, so that importing the module has no
    # side effects. These go to stderr because stdout carries the protonated
    # SMILES.

    # Always let the user know a help file is available.
    UtilFuncs.eprint("\nFor help, use: python dimorphite_dl.py --help")

    # And always report citation information.
    UtilFuncs.eprint("\nIf you use Dimorphite-DL in your research, please cite:")
    UtilFuncs.eprint(
        "Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An"
    )
    UtilFuncs.eprint(
        "open-source program for enumerating the ionization states of drug-like small"
    )
    UtilFuncs.eprint("molecules. J Cheminform 11:14. doi:10.1186/s13321-019-0336-9.\n")

    main()

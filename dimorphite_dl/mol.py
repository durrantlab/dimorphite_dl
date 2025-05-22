import copy
import importlib.resources as pkg_resources
from io import StringIO

from rdkit import Chem

from dimorphite_dl import data
from dimorphite_dl.io import LoadSMIFile, convert_smiles_str_to_mol


class Protonate(object):
    """A generator class for protonating SMILES strings, one at a time."""

    def __init__(self, args):
        """Initialize the generator.

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

        # Make the args an object variable variable.
        self.args = args
        # If the user provides a smiles string, turn it into a file-like StringIO
        # object.
        if "smiles" in args:
            if isinstance(args["smiles"], str):
                args["smiles_file"] = StringIO(args["smiles"])

        args["smiles_and_data"] = LoadSMIFile(args["smiles_file"], args)

        # A list to store the protonated SMILES strings associated with a
        # single input model.
        self.cur_prot_SMI = []

        # Make sure functions in ProtSubstructFuncs have access to the args.
        ProtSubstructFuncs.args = args

        # Load the substructures that can be protonated.
        self.subs = ProtSubstructFuncs.load_protonation_substructs_calc_state_for_ph(
            self.args["min_ph"], self.args["max_ph"], self.args["pka_precision"]
        )

    def __iter__(self):
        """Returns this generator object.

        Returns:
            This generator object.
        """

        return self

    def __next__(self):
        """Ensure Python3 compatibility.

        Returns:
            A dict, where the "smiles" key contains the canonical SMILES
                 string and the "data" key contains the remaining information
                 (e.g., the molecule name).
        """

        return self.next()

    def next(self):
        """Return the next protonated SMILES string.

        Returns:
            A dict, where the "smiles" key contains the canonical SMILES
                 string and the "data" key contains the remaining information
                 (e.g., the molecule name).
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

        # Keep track of the original smiles string for reporting, starting the
        # protonation process, etc.
        orig_smi = smile_and_datum["smiles"]

        # Dimorphite-DL may protonate some sites in ways that produce invalid
        # SMILES. We need to keep track of all smiles so we can "rewind" to
        # the last valid one, should things go south.
        properly_formed_smi_found = [orig_smi]

        # Everything on SMILES line but the SMILES string itself (e.g., the
        # molecule name).
        data = smile_and_datum["data"]

        # Collect the data associated with this smiles (e.g., the molecule
        # name).
        tag = " ".join(data)

        # sites is a list of (atom index, "PROTONATED|DEPROTONATED|BOTH",
        # reaction name, mol). Note that the second entry indicates what state
        # the site SHOULD be in (not the one it IS in per the SMILES string).
        # It's calculated based on the probablistic distributions obtained
        # during training.
        (sites, mol_used_to_idx_sites) = (
            ProtSubstructFuncs.get_prot_sites_and_target_states(orig_smi, self.subs)
        )

        new_mols = [mol_used_to_idx_sites]
        if len(sites) > 0:
            for site in sites:
                # Make a new smiles with the correct protonation state. Note that
                # new_smis is a growing list. This is how multiple protonation
                # sites are handled.
                new_mols = ProtSubstructFuncs.protonate_site(new_mols, site)
                if len(new_mols) > self.args["max_variants"]:
                    new_mols = new_mols[: self.args["max_variants"]]
                    if "silent" in self.args and not self.args["silent"]:
                        print(
                            "WARNING: Limited number of variants to "
                            + str(self.args["max_variants"])
                            + ": "
                            + orig_smi
                        )

                # Go through each of these new molecules and add them to the
                # properly_formed_smi_found, in case you generate a poorly
                # formed SMILES in the future and have to "rewind."
                properly_formed_smi_found += [Chem.MolToSmiles(m) for m in new_mols]
        else:
            # Deprotonate the mols (because protonate_site never called to do
            # it).
            mol_used_to_idx_sites = Chem.RemoveHs(mol_used_to_idx_sites)
            new_mols = [mol_used_to_idx_sites]

            # Go through each of these new molecules and add them to the
            # properly_formed_smi_found, in case you generate a poorly formed
            # SMILES in the future and have to "rewind."
            properly_formed_smi_found.append(Chem.MolToSmiles(mol_used_to_idx_sites))

        # In some cases, the script might generate redundant molecules.
        # Phosphonates, when the pH is between the two pKa values and the
        # stdev value is big enough, for example, will generate two identical
        # BOTH states. Let's remove this redundancy.
        new_smis = list(
            set(
                [
                    Chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
                    for m in new_mols
                ]
            )
        )

        # Sometimes Dimorphite-DL generates molecules that aren't actually
        # possible. Simply convert these to mol objects to eliminate the bad
        # ones (that are None).
        new_smis = [s for s in new_smis if convert_smiles_str_to_mol(s) is not None]

        # If there are no smi left, return the input one at the very least.
        # All generated forms have apparently been judged
        # inappropriate/malformed.
        if len(new_smis) == 0:
            properly_formed_smi_found.reverse()
            for smi in properly_formed_smi_found:
                if convert_smiles_str_to_mol(smi) is not None:
                    new_smis = [smi]
                    break

        # If the user wants to see the target states, add those to the ends of
        # each line.
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

    args = {}

    @staticmethod
    def load_substructre_smarts_file():
        """Loads the substructure smarts file. Similar to just using readlines,
        except it filters out comments (lines that start with "#").

        Returns:
            A list of the lines in the site_substructures.smarts file,
                 except blank lines and lines that start with "#"
        """

        # Read the file using pkg_resources
        with pkg_resources.open_text(data, "site_substructures.smarts") as f:
            # Filter out empty lines and comments
            lines = [l for l in f if l.strip() != "" and not l.startswith("#")]
        return lines

    @staticmethod
    def load_protonation_substructs_calc_state_for_ph(
        min_ph=6.4, max_ph=8.4, pka_std_range=1
    ):
        """A pre-calculated list of R-groups with protonation sites, with their
        likely pKa bins.

        Returns
            min_ph:  The lower bound on the pH range, defaults to 6.4.
            max_ph:  The upper bound on the pH range, defaults to 8.4.
            pka_std_range: Basically the precision (stdev from predicted pKa to
                              consider), defaults to 1.
        Returns:
            A dict of the protonation substructions for the specified pH
                 range.
        """

        subs = []

        for line in ProtSubstructFuncs.load_substructre_smarts_file():
            line = line.strip()
            sub = {}
            if line != "":
                splits = line.split()
                sub["name"] = splits[0]
                sub["smart"] = splits[1]
                sub["mol"] = Chem.MolFromSmarts(sub["smart"])

                pka_ranges = [splits[i : i + 3] for i in range(2, len(splits) - 1, 3)]

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

        Args:
            mean:   The mean pKa.
            std:    The precision (stdev).
            min_ph: The min pH of the range.
            max_ph: The max pH of the range.

        Returns:
            A string describing the protonation state.
        """

        min_pka = mean - std
        max_pka = mean + std

        # This needs to be reassigned, and 'ERROR' should never make it past
        # the next set of checks.
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

        Returns:
            smi: A SMILES string.
            subs: Substructure information.

        Returns:
            A list of protonation sites (atom index), pKa bin.
                ('PROTONATED', 'BOTH', or  'DEPROTONATED'), and reaction name.
                Also, the mol object that was used to generate the atom index.
        """

        # Convert the Smiles string (smi) to an RDKit Mol Obj
        mol_used_to_idx_sites = convert_smiles_str_to_mol(smi)

        # Check Conversion worked
        if mol_used_to_idx_sites is None:
            print("ERROR:   ", smi)
            return []

        # Try to Add hydrogens. if failed return []
        try:
            mol_used_to_idx_sites = Chem.AddHs(mol_used_to_idx_sites)
        except:
            print("ERROR:   ", smi)
            return []

        # Check adding Hs worked
        if mol_used_to_idx_sites is None:
            print("ERROR:   ", smi)
            return []

        ProtectUnprotectFuncs.unprotect_molecule(mol_used_to_idx_sites)
        protonation_sites = []

        for item in subs:
            smart = item["mol"]
            if mol_used_to_idx_sites.HasSubstructMatch(smart):
                matches = ProtectUnprotectFuncs.get_unprotected_matches(
                    mol_used_to_idx_sites, smart
                )
                prot = item["prot_states_for_pH"]
                for match in matches:
                    # We want to move the site from being relative to the
                    # substructure, to the index on the main molecule.
                    for site in prot:
                        proton = int(site[0])
                        category = site[1]
                        new_site = (match[proton], category, item["name"])

                        if not new_site in protonation_sites:
                            # Because sites must be unique.
                            protonation_sites.append(new_site)

                    ProtectUnprotectFuncs.protect_molecule(mol_used_to_idx_sites, match)

        return protonation_sites, mol_used_to_idx_sites

    @staticmethod
    def protonate_site(mols, site):
        """Given a list of molecule objects, we protonate the site.

        Returns:
            mols:  The list of molecule objects.
            site: Information about the protonation site.
                (idx, target_prot_state, prot_site_name)

        Returns:
            A list of the appropriately protonated molecule objects.
        """

        # Decouple the atom index and its target protonation state from the
        # site tuple
        idx, target_prot_state, prot_site_name = site

        state_to_charge = {"DEPROTONATED": [-1], "PROTONATED": [0], "BOTH": [-1, 0]}

        charges = state_to_charge[target_prot_state]

        # Now make the actual smiles match the target protonation state.
        output_mols = ProtSubstructFuncs.set_protonation_charge(
            mols, idx, charges, prot_site_name
        )

        return output_mols

    @staticmethod
    def set_protonation_charge(mols, idx, charges, prot_site_name):
        """Sets the atomic charge on a particular site for a set of SMILES.

        Returns:
            mols: A list of the input molecule objects.
            idx: The index of the atom to consider.
            charges: A list of the charges (ints) to assign at this site.
            prot_site_name: The name of the protonation site.

        Returns:
            A list of the processed (protonated/deprotonated) molecule
                 objects.
        """

        # Sets up the output list and the Nitrogen charge
        output = []

        for charge in charges:
            # The charge for Nitrogens is 1 higher than others (i.e.,
            # protonated state is positively charged).
            nitrogen_charge = charge + 1

            # But there are a few nitrogen moieties where the acidic group is
            # the neutral one. Amides are a good example. I gave some thought
            # re. how to best flag these. I decided that those
            # nitrogen-containing moieties where the acidic group is neutral
            # (rather than positively charged) will have "*" in the name.
            if "*" in prot_site_name:
                nitrogen_charge = nitrogen_charge - 1  # Undo what was done previously.

            for mol in mols:
                # Make a copy of the molecule.
                mol_copy = copy.deepcopy(mol)

                # Remove hydrogen atoms.
                try:
                    mol_copy = Chem.RemoveHs(mol_copy)
                except:
                    if (
                        "silent" in ProtSubstructFuncs.args
                        and not ProtSubstructFuncs.args["silent"]
                    ):
                        print(
                            "WARNING: Skipping poorly formed SMILES string: "
                            + Chem.MolToSmiles(mol_copy)
                        )
                    continue

                atom = mol_copy.GetAtomWithIdx(idx)

                explicit_bond_order_total = sum(
                    [b.GetBondTypeAsDouble() for b in atom.GetBonds()]
                )

                # Assign the protonation charge, with special care for
                # nitrogens
                element = atom.GetAtomicNum()
                if element == 7:
                    atom.SetFormalCharge(nitrogen_charge)

                    # Need to figure out how many hydrogens to add.
                    if nitrogen_charge == 1 and explicit_bond_order_total == 1:
                        atom.SetNumExplicitHs(3)
                    elif nitrogen_charge == 1 and explicit_bond_order_total == 2:
                        atom.SetNumExplicitHs(2)
                    elif nitrogen_charge == 1 and explicit_bond_order_total == 3:
                        atom.SetNumExplicitHs(1)
                    elif nitrogen_charge == 0 and explicit_bond_order_total == 1:
                        atom.SetNumExplicitHs(2)
                    elif nitrogen_charge == 0 and explicit_bond_order_total == 2:
                        atom.SetNumExplicitHs(1)
                    elif nitrogen_charge == -1 and explicit_bond_order_total == 2:
                        atom.SetNumExplicitHs(0)
                    elif nitrogen_charge == -1 and explicit_bond_order_total == 1:
                        atom.SetNumExplicitHs(1)
                    #### JDD
                else:
                    atom.SetFormalCharge(charge)
                    if element == 8 or element == 16:  # O and S
                        if charge == 0 and explicit_bond_order_total == 1:
                            atom.SetNumExplicitHs(1)
                        elif charge == -1 and explicit_bond_order_total == 1:
                            atom.SetNumExplicitHs(0)

                # Deprotonating protonated aromatic nitrogen gives [nH-]. Change this
                # to [n-].
                if "[nH-]" in Chem.MolToSmiles(mol_copy):
                    atom.SetNumExplicitHs(0)

                mol_copy.UpdatePropertyCache(strict=False)
                # prod.UpdatePropertyCache(strict=False)

                output.append(mol_copy)

        return output


class ProtectUnprotectFuncs:
    """A namespace for storing functions that are useful for protecting and
    unprotecting molecules. To keep things organized. We need to identify and
    mark groups that have been matched with a substructure."""

    @staticmethod
    def unprotect_molecule(mol):
        """Sets the protected property on all atoms to 0. This also creates the
        property for new molecules.

        Args:
            mol: The rdkit Mol object.
        """

        for atom in mol.GetAtoms():
            atom.SetProp("_protected", "0")

    @staticmethod
    def protect_molecule(mol, match):
        """Given a 'match', a list of molecules idx's, we set the protected status
        of each atom to 1. This will prevent any matches using that atom in the
        future.

        Returns
            mol: The rdkit Mol object to protect.
            match: A list of molecule idx's.
        """

        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            atom.SetProp("_protected", "1")

    @staticmethod
    def get_unprotected_matches(mol, substruct):
        """Finds substructure matches with atoms that have not been protected.
        Returns list of matches, each match a list of atom idxs.

        Args:
            mol: The Mol object to consider.
            substruct: The SMARTS string of the substructure ot match.

        Returns:
            A list of the matches. Each match is itself a list of atom idxs.
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

        Args:
            mol: The Mol object to check.
            match: The match to check.

        Returns:
            A boolean, whether the match is present or not.
        """

        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            protected = atom.GetProp("_protected")
            if protected == "1":
                return False
        return True

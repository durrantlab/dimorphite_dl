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

import copy
import os

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def protonate(args):
    """Protonates a set of molecules as given by the user inputs.
    
    :param dict args: A dictionary containing the arguments.
    :return: A list of the protonated smiles strings.
    """

    args = clean_args(args)
    subs = load_protonation_substructs(args["min_ph"], args["max_ph"], args["st_dev"])
    smiles = args["smiles"]

    # Note that args["data"] includes everything on SMILES line but the SMILES
    # string itself (e.g., the molecule name). It is set in clean_args().
    data = args["data"]  

    output = []
    for i, smi in enumerate(smiles):
        if smi.startswith("NONE|"):
            print "ERROR: Skipping poorly formed SMILES string: " + smi[5:] + "\t" + " ".join(data[i])
            continue

        # Collect the data associated with this smiels (e.g., the molecule
        # name).
        tag = " ".join(data[i])
        
        # sites is a list of (atom index, "PROTONATED|DEPROTONATED|BOTH").
        # Note that the second entry indicates what state the site SHOULD be
        # in (not the one it IS in per the SMILES string). It's calculated
        # based on the probablistic distributions obtained during training.
        sites = get_prot_sites_and_target_states(smi, subs)
        
        new_smis = [smi]
        for site in sites:
            # Make a new smiles with the correct protonation state. Note that
            # new_smis is a growing list. This is how multiple protonation
            # sites are handled.
            new_smis = protonate_site(new_smis, site)
        new_lines = [x + '\t' + tag for x in new_smis]
        output.extend(new_lines)

    return output

def clean_args(args):
    """Cleans and normalizes input parameters
    
    :param dict args: A dictionary containing the arguments.
    :raises Exception: No SMILES in params.
    :return: A dict of the arguments, now fixed.
    """

    defaults = {'min_ph' : 6.4,
                'max_ph' : 8.4,
                'st_dev' : 1.5}

    for key in defaults:
        if key not in args:
            args[key] = defaults[key]

    keys = list(args.keys())
    for key in keys:
        if args[key] is None:
            del args[key]

    if "smiles" in args:
        if isinstance(args["smiles"], str):
            splits = args["smiles"].strip().split()
            args["smiles"] = [splits[0]]
            args["data"] = [splits[1:]]
    elif "smiles_file" in args:
        args["smiles"], args["data"] = load_files(args["smiles_file"])
    else:
        raise Exception("Error: No SMILES in params.")

    mol_str_list = []
    for i, smiles_str in enumerate(args["smiles"]):
        
        # Convert from SMILES string to RDKIT Mol
        # Filter if failed.

        mol = convert_smiles_str_to_mol(smiles_str)
        if mol is None:
            mol_str_list.append("NONE|" + smiles_str)
            continue

        # Handle nuetralizing the molecules
        # Filter if failed.
        mol = neutralize_mol(mol)
        if mol is None:
            mol_str_list.append("NONE|" + smiles_str)
            continue

        try:
            mol = Chem.RemoveHs(mol)
        except:
            mol_str_list.append("NONE|" + smiles_str)
            continue

        if mol is None:
            mol_str_list.append("NONE|" + smiles_str)
            continue


        new_mol_string = Chem.MolToSmiles(mol)
        mol_str_list.append(new_mol_string)

    args["smiles"] = [x for x in mol_str_list]

    return args

def neutralize_mol(mol):
    """All molecules need to be neuralized to the extent possible. The user
    should not be allowed to specify the valence of the atoms in most cases.

    :param rdkit.Chem.rdchem.Mol mol: The rdkit Mol objet to be neutralized.
    :return: The neutralized Mol object.
    """

    # Get the reaction data
    rxn_data = [
        ['[Ov1-1:1]', '[Ov2+0:1]-[H]'],  # To handle O- bonded to only one atom (add hydrogen).
        ['[#7v4+1:1]-[H]', '[#7v3+0:1]'],  # To handle N+ bonded to a hydrogen (remove hydrogen).
        ['[Ov2-:1]', '[Ov2+0:1]'],  # To handle O- bonded to two atoms. Should not be Negative.
        ['[#7v3+1:1]', '[#7v3+0:1]'],  # To handle N+ bonded to three atoms. Should not be positive.
        ['[#7v2-1:1]', '[#7+0:1]-[H]'],  # To handle N- Bonded to two atoms. Add hydrogen.
        # ['[N:1]=[N+0:2]=[N:3]-[H]', '[N:1]=[N+1:2]=[N+0:3]-[H]'],  # To handle bad azide. Must be protonated. (Now handled elsewhere, before SMILES converted to Mol object.)
        ['[H]-[N:1]-[N:2]#[N:3]', '[N:1]=[N+1:2]=[N:3]-[H]']  # To handle bad azide. R-N-N#N should be R-N=[N+]=N
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
        rxn = None  # The reaction to perform.
        rxn_str = None

        for i, rxn_datum in enumerate(rxn_data):
            reactant_smarts, product_smarts, substruct_match_mol, rxn2 = rxn_datum
            if mol.HasSubstructMatch(substruct_match_mol):
                if rxn2 is None:
                    rxn_str = reactant_smarts + '>>' + product_smarts
                    rxn = AllChem.ReactionFromSmarts(rxn_str)
                    rxn_data[i][3] = rxn
                break

        # Perform the reaction if necessary
        if rxn is None:  # No reaction left, so break out of while loop.
            break
        else:
            mol = rxn.RunReactants((mol,))[0][0]
            mol.UpdatePropertyCache(strict=False)  # Update valences

    # The mols have been altered from the reactions described above, we need to resanitize them.
    # Make sure aromatic rings are shown as such
    # This catches all RDKit Errors. without the catchError and sanitizeOps
    # the Chem.SanitizeMol can crash the program.
    sanitize_string =  Chem.SanitizeMol(
        mol, 
        sanitizeOps=rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, 
        catchErrors = True
    )

    return mol if sanitize_string.name == "SANITIZE_NONE" else None

def load_files(smile_file):
    """Loads smiles from file.
    
    :param string smile_file: The smiles file name.
    :return: (list, list), the smiles strings and associated data, 
             respectively.
    """

    smiles = []
    data = []
    with open(smile_file, 'r') as smis:
        for line in smis:
            splits = line.split()
            if len(splits) != 0:
                smiles.append(splits[0])
                data.append(splits[1:])
    return smiles, data

def load_protonation_substructs(min_ph=6.4, max_ph=8.4, pka_std_range=1):
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
    pwd = os.path.dirname(__file__)
    site_structures_file = "{}/{}".format(pwd,"site_substructures.smarts")
    with open(site_structures_file, 'r') as substruct:
        for line in substruct:
            line = line.strip()
            sub = {}
            if line is not "":
                splits = line.split()
                sub["name"] = splits[0]
                sub["smart"] = splits[1]
                sub["mol"] = Chem.MolFromSmarts(sub["smart"])

                #NEED TO DIVIDE THIS BY 3s
                pka_ranges = [splits[i:i+3] for i in range(2, len(splits)-1, 3)]

                prot = []
                for pka_range in pka_ranges:
                    site = pka_range[0]
                    std = float(pka_range[2]) * pka_std_range
                    mean = float(pka_range[1])
                    protonation_state = define_protonation_state(mean, std, min_ph, \
                        max_ph)

                    prot.append([site, protonation_state])

                sub["prot"] = prot
                subs.append(sub)
    return subs

def define_protonation_state(mean, std, min_ph, max_ph):
    """Updates the substructure definitions to include the protonation state
    based on the user-given pH range. The size of the pKa range is also based
    on the number of standard deviations to be considered by the user param.
    
    :param float mean:   The mean pKa.
    :param float std:    The precision (stdev).
    :param float min_ph: The min pH of the range.
    :param float max_ph: The max pH of the range.
    :raises Exception:   HORRIBLE NONSENSE HAS OCCURED.
    :return: A string describing the protonation state.
    """

    min_pka = mean - std
    max_pka = mean + std

    # This needs to be reassigned, and 'ERROR' should never make it past the
    # next set of checks.
    protonation_state = 'ERROR'

    if min_pka <= max_ph and min_ph <= max_pka:
        protonation_state = 'BOTH'
    elif mean > max_ph:
        protonation_state = 'PROTONATED'
    elif mean < min_ph:
        protonation_state = 'DEPROTONATED'

    # We are error handling here
    if protonation_state == 'ERROR':
        raise Exception("HORRIBLE NONSENSE HAS OCCURED.")

    return protonation_state


###
# We need to identify and mark groups that have been matched with a substructure.
###

def unprotect_molecule(mol):
    """Sets the protected property on all atoms to 0. This also creates the
    property for new molecules.
    
    :param rdkit.Chem.rdchem.Mol mol: The rdkit Mol object.
    :type mol: The rdkit Mol object with atoms unprotected.
    """

    for atom in mol.GetAtoms():
        atom.SetProp('_protected', '0')

def protect_molecule(mol, match):
    """Given a 'match', a list of molecules idx's, we set the protected status
    of each atom to 1. This will prevent any matches using that atom in the
    future.
    
    :param rdkit.Chem.rdchem.Mol mol: The rdkit Mol object to protect.
    :param list match: A list of molecule idx's.
    """

    for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        atom.SetProp('_protected', '1')

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
        if is_match_unprotected(mol, match):
            unprotected_matches.append(match)
    return unprotected_matches

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

def neutralize_molecule(mol):
    """Neutralize things. Maybe?
    
    :param rdkit.Chem.rdchem.Mol mol: The Mol object to conside.r
    """

    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)

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
    mol = convert_smiles_str_to_mol(smi)

    # Check Conversion worked
    if mol is None:
        print("ERROR:   ", smi)
        return []

    # Try to Add hydrogens. if failed return []
    try:
        mol =  Chem.AddHs(mol)
    except:
        print("ERROR:   ", smi)
        return []

    # Check adding Hs worked
    if mol is None:
        print("ERROR:   ", smi)
        return []

    unprotect_molecule(mol)
    protonation_sites = []

    for item in subs:
        smart = item['mol']
        if mol.HasSubstructMatch(smart):
            matches = get_unprotected_matches(mol, smart)
            prot = item['prot']
            for match in matches:
                # We want to move the site from being relative to the
                # substructure, to the index on the main molecule.
                for site in prot:
                    proton = int(site[0])
                    category = site[1]
                    new_site = (match[proton], category, item["name"])
                    protonation_sites.append(new_site)
                protect_molecule(mol, match)
    return protonation_sites

def protonate_site(smis, site):
    """Given a list of SMILES strings, we protonate the site.
    
    :param list smis:  The list of SMILES strings.
    :param tuple site: Information about the protonation site. 
                       (idx, target_prot_state, prot_site_name)
    :return: A list of the appropriately protonated SMILES.
    """

    # Decouple the atom index and its target protonation state from the site
    # tuple
    idx, target_prot_state, prot_site_name = site

    # Initialize the output list
    output_smis = []

    state_to_charge = {"DEPROTONATED": [-1],
                       "PROTONATED": [0],
                       "BOTH": [-1, 0]}

    charges = state_to_charge[target_prot_state]

    # Now make the actual smiles match the target protonation state.
    output_smis = set_protonation_charge(smis, idx, charges, prot_site_name)

    return output_smis

def set_protonation_charge(smis, idx, charges, prot_site_name):
    """Sets the atomic charge on a particular site for a set of SMILES.
    
    :param list smis:             A list of the SMILES strings.
    :param int idx:               The index of the atom to consider.
    :param list charges:          A list of the charges (ints) to assign at
                                  this site.
    :param string prot_site_name: The name of the protonation site.
    :return: A list of the processed SMILES strings.
    """

    # Sets up the output list and the Nitrogen charge
    output = []

    for charge in charges:
        # The charge for Nitrogens is 1 higher than others (i.e., protonated state
        # is positively charged).
        nitro_charge = charge + 1

        # But there are a few nitrogen moieties where the acidic group is the
        # neural one. Amides are a good example. I gave some thought re. how
        # to best flag these. I decided that those nitrogen-containing
        # moieties where the acidic group is neutral (rather than positively
        # charged) will have "*" in the name.
        if "*" in prot_site_name:
            nitro_charge = nitro_charge - 1  # Undo what was done previously.

        for smi in smis:
            
            # Convert smilesstring (smi) into a RDKit Mol
            mol = convert_smiles_str_to_mol(smi)
             
            # Check that the conversion worked, skip if it fails
            if mol is None:
                continue    
                            
            atom = mol.GetAtomWithIdx(idx)

            # Assign the protonation charge, with special care for Nitrogens
            element = atom.GetAtomicNum()
            if element == 7:
                atom.SetFormalCharge(nitro_charge)
            else:
                atom.SetFormalCharge(charge)

            # Convert back to SMILE and add to output
            out_smile = Chem.MolToSmiles(mol, isomericSmiles=True,canonical=True)
            output.append(out_smile)

    return output

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
    except:
        return None

    # Check that there are None type errors Chem.MolFromSmiles has sanitize on
    # which means if there is even a small error in the SMILES (kekulize,
    # nitrogen charge...) then mol=None. ie.
    # Chem.MolFromSmiles("C[N]=[N]=[N]") = None this is an example of an
    # nitrogen charge error. It is cased in a try statement to be overly
    # cautious.
     
    return None if mol is None else mol


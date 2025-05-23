from loguru import logger


def unprotect_molecule(mol):
    """Sets the protected property on all atoms to 0."""
    logger.trace("Unprotecting each atom")
    for atom in mol.GetAtoms():
        atom.SetProp("_protected", "0")


def protect_molecule(mol, match):
    """Given a 'match', a list of molecules idx's, we set the protected status
    of each atom to 1."""
    logger.trace("Protecting atom(s)")
    for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        atom.SetProp("_protected", "1")


def get_unprotected_matches(mol, substruct):
    """Finds substructure matches with atoms that have not been protected."""
    matches = mol.GetSubstructMatches(substruct)
    logger.debug("Found {} substructure match(es)", len(matches))
    unprotected_matches = []
    for match in matches:
        if is_match_unprotected(mol, match):
            unprotected_matches.append(match)
    logger.debug(
        "{}/{} matches were unprotected", len(unprotected_matches), len(matches)
    )
    return unprotected_matches


def is_match_unprotected(mol, match):
    """Checks a molecule to see if the substructure match contains any
    protected atoms."""
    for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        protected = atom.GetProp("_protected")
        if protected == "1":
            return False
    return True

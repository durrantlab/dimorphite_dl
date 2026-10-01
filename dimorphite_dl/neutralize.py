from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem

RXN_DATA = (
    # Charge-separated sulfones and sulfonates (e.g., C[S+2](C)([O-])[O-]).
    # Must run before the O- rule, which would otherwise protonate these
    # oxygens. The resulting S+-O- is collapsed by the next rule.
    ("[#16+2:1]-[Ov1-1:2]", "[#16+1:1]=[O+0:2]"),
    # Charge-separated sulfoxides (e.g., C[S+](C)[O-]).
    ("[#16+1:1]-[Ov1-1:2]", "[#16+0:1]=[O+0:2]"),
    # Charge-separated phosphine oxides and phosphates (e.g.,
    # C[P+](C)(C)[O-]). The phosphate sites expect P=O.
    ("[#15+1:1]-[Ov1-1:2]", "[#15+0:1]=[O+0:2]"),
    # O- bonded to only one atom (add hydrogen). The O- of an N-oxide or
    # nitrone is left alone: no site pattern restores it from N+-OH, whereas a
    # nitro O-H is deprotonated again by Nitro. Nitrate gets only one H,
    # giving nitric acid (O[N+](=O)[O-]); a second OH would match Nitro as a
    # second site and survive as neutral HNO3.
    (
        "[Ov1-1;!$([O-]-[#7+;!$([#7+]=O)]);!$([O-]-[#7+](=O)-[#8]-[H]):1]",
        "[Ov2+0:1]-[H]",
    ),
    # S- bonded to only one atom (add hydrogen). The thiol site patterns all
    # require S-H.
    ("[Sv1-1:1]", "[Sv2+0:1]-[H]"),
    # To handle N+ bonded to a hydrogen (remove hydrogen).
    ("[#7v4+1:1]-[H]", "[#7v3+0:1]"),
    # To handle O- bonded to two atoms. Should not be Negative.
    ("[Ov2-:1]", "[Ov2+0:1]"),
    # To handle N+ bonded to three atoms. Should not be positive.
    ("[#7v3+1:1]", "[#7v3+0:1]"),
    # N- bonded to two atoms (add hydrogen). The terminal N- of a diazo group
    # (C=[N+]=[N-]) is excluded: no site pattern removes that proton, so it
    # would stay a +1 cation at every pH. Azides still match, since their far
    # atom is nitrogen.
    ("[#7v2-1;!$([#7-]=[#7+]=[#6]):1]", "[#7+0:1]-[H]"),
    # To handle bad azide. R-N-N#N should be R-N=[N+]=N.
    ("[H]-[N:1]-[N:2]#[N:3]", "[N:1]=[N+1:2]=[N:3]-[H]"),
    # Amidinium/guanidinium drawn with the charge on a substituted N. Moving
    # the charge onto the N-H lets the N+-H rule remove the proton on the next
    # pass.
    (
        "[#7+1;H0;!$([#7+]~[O-]):1]=[#6:2]-[#7+0;!H0:3]",
        "[#7+0:1]-[#6:2]=[#7+1:3]",
    ),
    # Same, for imidazolium-type rings (e.g., C[n+]1cc[nH]c1).
    ("[n+1;H0;!$([n+]~[O-]):1]:[c:2]:[n+0;!H0:3]", "[n+0:1]:[c:2]:[n+1:3]"),
    # Same, for pyrazolium-type rings (e.g., C[n+]1ccc[nH]1).
    ("[n+1;H0;!$([n+]~[O-]):1]:[n+0;!H0:2]", "[n+0:1]:[n+1:2]"),
)

# Each pass fixes one charged atom, so a sane molecule needs at most a few
# passes per atom. More than this means a rule keeps recreating its own
# reactant, which would otherwise loop forever.
NEUTRALIZE_PASSES_PER_ATOM = 10


class NeutralizationReaction:
    """
    Represents a single neutralization reaction defined by a pair of SMARTS strings
    """

    def __init__(self, smarts_reactant: str, smarts_product: str):
        """
        Args:
            smarts_reactant: SMARTS for detecting the reactants of a defined
                neutralization reaction.
            smarts_product: SMARTS for what the detected `smarts_reactant` should
                be transformed to.
        """
        self.smarts_reactant = smarts_reactant
        self.smarts_product = smarts_product
        self._pattern = Chem.MolFromSmarts(smarts_reactant)
        self._rxn = AllChem.ReactionFromSmarts(f"{smarts_reactant}>>{smarts_product}")

    def __str__(self) -> str:
        return f"{self.smarts_reactant} >> {self.smarts_product}"

    def __repr__(self) -> str:
        return self.__str__()

    def matches(self, mol: Chem.Mol) -> bool:
        """Check if this reaction can be applied to the given molecule."""
        return mol.HasSubstructMatch(self._pattern)

    def apply(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Apply the neutralization reaction to the molecule. Returns the first product.
        If multiple products are generated, only the first is returned.
        """
        # Only the first product is used, so enumerating RDKit's default of up
        # to 1000 is wasted work on symmetric molecules.
        products = self._rxn.RunReactants((mol,), 1)
        if products:
            # products is a tuple of tuples; take the first product set, first product
            return products[0][0]
        return mol


class ReactionRegistry:
    """
    Holds a collection of NeutralizationReaction objects and applies them repeatedly
    until no further matches are found.
    """

    def __init__(self, rxn_data: tuple[tuple[str, str]]):
        self.reactions = []
        for reactant, product in rxn_data:
            self.reactions.append(NeutralizationReaction(reactant, product))

    def neutralize(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Apply all registered neutralization reactions to the molecule in a loop
        until no further transformations are possible. Assumes explicit H atoms
        have already been added.
        """
        mol.UpdatePropertyCache(strict=False)
        input_smiles = Chem.MolToSmiles(mol)
        max_passes = NEUTRALIZE_PASSES_PER_ATOM * mol.GetNumAtoms()
        passes = 0
        changed = True
        while changed:
            changed = False
            for reaction in self.reactions:
                if reaction.matches(mol):
                    logger.debug("Found reaction match: {}", str(reaction))
                    passes += 1
                    if passes > max_passes:
                        raise RuntimeError(
                            f"Neutralization did not converge after {max_passes} "
                            f"passes on {input_smiles}; rule {reaction} keeps "
                            "matching its own product."
                        )
                    mol = reaction.apply(mol)
                    mol.UpdatePropertyCache(strict=False)
                    changed = True
                    break  # restart scanning from first reaction
                else:
                    logger.trace("No match to reaction: {}", str(reaction))
        # Final sanitization
        sanitized = Chem.SanitizeMol(
            mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors=True
        )
        if sanitized.name == "SANITIZE_NONE":
            logger.debug("After neutralizing: {}", Chem.MolToSmiles(mol))
            return mol
        raise RuntimeError("Ran into issue sanitizing mol")


class MoleculeNeutralizer:
    """
    High-level class to take SMILES, handle preprocessing, add Hs,
    run neutralization, and return a clean SMILES.
    """

    def __init__(self, rxn_data: tuple[tuple[str, str]] | None = None):
        if rxn_data is None:
            rxn_data = RXN_DATA
        self.registry = ReactionRegistry(rxn_data)

    def neutralize_smiles(self, smiles: str) -> str | None:
        logger.debug("Neutralizing {}", smiles)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        # Add explicit Hs
        mol = Chem.AddHs(mol)
        logger.debug("After adding hydrogens: {}", Chem.MolToSmiles(mol))
        # Run neutralization
        mol = self.registry.neutralize(mol)
        # Remove explicit Hs
        mol = Chem.RemoveHs(mol)
        logger.debug("After removing hydrogens: {}", Chem.MolToSmiles(mol))
        # Generate final SMILES
        return Chem.MolToSmiles(mol)

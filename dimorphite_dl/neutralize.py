from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem

RXN_DATA = (
    # To handle O- bonded to only one atom (add hydrogen).
    ("[Ov1-1:1]", "[Ov2+0:1]-[H]"),
    # To handle N+ bonded to a hydrogen (remove hydrogen).
    ("[#7v4+1:1]-[H]", "[#7v3+0:1]"),
    # To handle O- bonded to two atoms. Should not be Negative.
    ("[Ov2-:1]", "[Ov2+0:1]"),
    # To handle N+ bonded to three atoms. Should not be positive.
    ("[#7v3+1:1]", "[#7v3+0:1]"),
    # To handle N- Bonded to two atoms. Add hydrogen.
    ("[#7v2-1:1]", "[#7+0:1]-[H]"),
    # To handle bad azide. R-N-N#N should be R-N=[N+]=N.
    ("[H]-[N:1]-[N:2]#[N:3]", "[N:1]=[N+1:2]=[N:3]-[H]"),
)


class NeutralizationReaction:
    """
    Represents a single neutralization reaction defined by a pair of SMARTS strings
    """

    def __init__(self, reactant_smarts: str, product_smarts: str):
        self.reactant_smarts = reactant_smarts
        self.product_smarts = product_smarts
        self._pattern = Chem.MolFromSmarts(reactant_smarts)
        self._rxn = AllChem.ReactionFromSmarts(f"{reactant_smarts}>>{product_smarts}")

    def __str__(self) -> str:
        return f"{self.reactant_smarts} >> {self.product_smarts}"

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
        products = self._rxn.RunReactants((mol,))
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
        changed = True
        while changed:
            changed = False
            for reaction in self.reactions:
                if reaction.matches(mol):
                    logger.debug("Found reaction match: {}", str(reaction))
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
    High-level class to take SMILES, handle preprocessing (like azides), add Hs,
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

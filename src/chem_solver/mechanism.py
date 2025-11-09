"""
Chemical Mechanism Parser

Parses chemistry mechanisms from JSON files and prepares them for solving.
"""

import json
import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class Reaction:
    """Represents a single chemical reaction"""
    reactants: Dict[str, float]  # species name -> stoichiometric coefficient
    products: Dict[str, float]
    rate_constant: Optional[float] = None
    arrhenius_A: Optional[float] = None
    arrhenius_Ea: Optional[float] = None  # Activation energy in K
    arrhenius_n: Optional[float] = 0.0
    photolysis: Optional[str] = None  # Photolysis rate key if applicable

    def compute_rate_constant(self, temperature: float, photolysis_rates: Dict[str, float] = None) -> float:
        """
        Compute rate constant for this reaction.

        Args:
            temperature: Temperature in K
            photolysis_rates: Dictionary of photolysis rates (if photolysis reaction)

        Returns:
            Rate constant
        """
        if self.photolysis:
            if photolysis_rates and self.photolysis in photolysis_rates:
                return photolysis_rates[self.photolysis]
            return 0.0

        if self.rate_constant is not None:
            return self.rate_constant

        if self.arrhenius_A is not None:
            # Arrhenius equation: k = A * T^n * exp(-Ea/T)
            return self.arrhenius_A * (temperature ** self.arrhenius_n) * np.exp(-self.arrhenius_Ea / temperature)

        return 0.0


class ChemicalMechanism:
    """
    Chemical mechanism handler that parses and manages reaction mechanisms.
    """

    def __init__(self, mechanism_file: str):
        """
        Initialize mechanism from JSON file.

        Args:
            mechanism_file: Path to JSON mechanism file
        """
        self.mechanism_file = mechanism_file
        self.species: List[str] = []
        self.reactions: List[Reaction] = []
        self.species_idx: Dict[str, int] = {}

        self._load_mechanism()

    def _load_mechanism(self):
        """Load and parse the mechanism file"""
        with open(self.mechanism_file, 'r') as f:
            data = json.load(f)

        # Load species
        self.species = data.get('species', [])
        self.species_idx = {sp: idx for idx, sp in enumerate(self.species)}

        # Load reactions
        reactions_data = data.get('reactions', [])
        for rxn_data in reactions_data:
            reactants = rxn_data.get('reactants', {})
            products = rxn_data.get('products', {})

            # Parse rate information
            rate_info = rxn_data.get('rate', {})

            reaction = Reaction(
                reactants=reactants,
                products=products,
                rate_constant=rate_info.get('constant'),
                arrhenius_A=rate_info.get('A'),
                arrhenius_Ea=rate_info.get('Ea'),
                arrhenius_n=rate_info.get('n', 0.0),
                photolysis=rate_info.get('photolysis')
            )

            self.reactions.append(reaction)

    def get_num_species(self) -> int:
        """Return number of species"""
        return len(self.species)

    def get_species_index(self, species_name: str) -> int:
        """Get index of a species"""
        return self.species_idx.get(species_name, -1)

    def compute_reaction_rates(self, concentrations: np.ndarray, temperature: float,
                              pressure: float, photolysis_rates: Dict[str, float]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute production and loss rates for all species.

        Args:
            concentrations: Species concentrations (molecules/cm^3)
            temperature: Temperature in K
            pressure: Pressure in Pa
            photolysis_rates: Dictionary of photolysis rates

        Returns:
            Tuple of (production_rates, loss_rates) for each species
        """
        n_species = len(self.species)
        production = np.zeros(n_species)
        loss = np.zeros(n_species)

        for reaction in self.reactions:
            # Compute rate constant
            k = reaction.compute_rate_constant(temperature, photolysis_rates)

            # Compute reaction rate
            rate = k
            for reactant, stoich in reaction.reactants.items():
                idx = self.species_idx.get(reactant, -1)
                if idx >= 0:
                    rate *= concentrations[idx] ** stoich

            # Add to production and loss
            for reactant, stoich in reaction.reactants.items():
                idx = self.species_idx.get(reactant, -1)
                if idx >= 0:
                    loss[idx] += stoich * rate

            for product, stoich in reaction.products.items():
                idx = self.species_idx.get(product, -1)
                if idx >= 0:
                    production[idx] += stoich * rate

        return production, loss

    def compute_derivatives(self, concentrations: np.ndarray, temperature: float,
                          pressure: float, photolysis_rates: Dict[str, float]) -> np.ndarray:
        """
        Compute time derivatives (dC/dt) for all species.

        Args:
            concentrations: Species concentrations
            temperature: Temperature in K
            pressure: Pressure in Pa
            photolysis_rates: Dictionary of photolysis rates

        Returns:
            Time derivatives for each species
        """
        production, loss = self.compute_reaction_rates(concentrations, temperature, pressure, photolysis_rates)
        return production - loss

    def compute_jacobian(self, concentrations: np.ndarray, temperature: float,
                        pressure: float, photolysis_rates: Dict[str, float],
                        epsilon: float = 1e-6) -> np.ndarray:
        """
        Compute Jacobian matrix for the system (numerical approximation).

        Args:
            concentrations: Species concentrations
            temperature: Temperature in K
            pressure: Pressure in Pa
            photolysis_rates: Dictionary of photolysis rates
            epsilon: Perturbation for numerical derivative

        Returns:
            Jacobian matrix (n_species x n_species)
        """
        n_species = len(self.species)
        jacobian = np.zeros((n_species, n_species))

        f0 = self.compute_derivatives(concentrations, temperature, pressure, photolysis_rates)

        for i in range(n_species):
            conc_perturbed = concentrations.copy()
            delta = epsilon * max(concentrations[i], 1e-10)
            conc_perturbed[i] += delta

            f_perturbed = self.compute_derivatives(conc_perturbed, temperature, pressure, photolysis_rates)
            jacobian[:, i] = (f_perturbed - f0) / delta

        return jacobian

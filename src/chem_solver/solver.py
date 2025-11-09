"""
SMVGEAR-like Chemical Solver

Implementation of a chemical kinetics solver using implicit methods
for stiff atmospheric chemistry problems.
"""

import numpy as np
from typing import Dict, List, Optional, Union
from multiprocessing import Pool, cpu_count
import warnings

from .mechanism import ChemicalMechanism


class ChemicalSolver:
    """
    SMVGEAR-like chemical kinetics solver for atmospheric chemistry.

    This solver uses implicit integration methods (Backward Euler with iterations)
    to solve stiff chemical kinetics problems efficiently.
    """

    def __init__(self, mechanism_file: str, time_step: float,
                 rtol: float = 1e-3, atol: float = 1e-10,
                 max_iterations: int = 100):
        """
        Initialize the chemical solver.

        Args:
            mechanism_file: Path to JSON mechanism file
            time_step: Default time step in seconds
            rtol: Relative tolerance for convergence
            atol: Absolute tolerance for convergence
            max_iterations: Maximum number of Newton iterations
        """
        self.mechanism = ChemicalMechanism(mechanism_file)
        self.time_step = time_step
        self.rtol = rtol
        self.atol = atol
        self.max_iterations = max_iterations

    def _backward_euler_step(self, concentrations: np.ndarray, dt: float,
                            temperature: float, pressure: float,
                            photolysis_rates: Dict[str, float]) -> np.ndarray:
        """
        Perform one Backward Euler step with Newton iterations.

        Solves: C_new = C_old + dt * f(C_new, T, P, J)

        Args:
            concentrations: Initial concentrations
            dt: Time step
            temperature: Temperature in K
            pressure: Pressure in Pa
            photolysis_rates: Photolysis rates dictionary

        Returns:
            Updated concentrations
        """
        C_old = concentrations.copy()
        C_new = concentrations.copy()

        for iteration in range(self.max_iterations):
            # Compute derivatives and Jacobian at current guess
            f = self.mechanism.compute_derivatives(C_new, temperature, pressure, photolysis_rates)
            J = self.mechanism.compute_jacobian(C_new, temperature, pressure, photolysis_rates)

            # Residual: R = C_new - C_old - dt * f(C_new)
            residual = C_new - C_old - dt * f

            # Check convergence
            error = np.linalg.norm(residual) / (np.linalg.norm(C_new) + 1e-30)
            if error < self.rtol:
                break

            # Newton step: (I - dt*J) * delta_C = -R
            # delta_C = C_new^(k+1) - C_new^(k)
            try:
                A = np.eye(len(C_new)) - dt * J
                delta_C = np.linalg.solve(A, -residual)
                C_new = C_new + delta_C

                # Ensure non-negative concentrations
                C_new = np.maximum(C_new, 0.0)

            except np.linalg.LinAlgError:
                warnings.warn(f"Singular matrix in Newton iteration {iteration}, using explicit step")
                # Fall back to simple explicit step if matrix is singular
                C_new = C_old + dt * self.mechanism.compute_derivatives(C_old, temperature, pressure, photolysis_rates)
                C_new = np.maximum(C_new, 0.0)
                break

        return C_new

    def solve_step(self, concentrations: Dict[str, float],
                   temperature: float, pressure: float,
                   photolysis_rates: Dict[str, float],
                   time_step: Optional[float] = None) -> Dict[str, float]:
        """
        Perform a single integration step.

        Args:
            concentrations: Dictionary of species concentrations (molecules/cm^3)
            temperature: Temperature in K
            pressure: Pressure in Pa
            photolysis_rates: Dictionary of photolysis rates
            time_step: Time step in seconds (uses default if None)

        Returns:
            Dictionary of updated concentrations
        """
        dt = time_step if time_step is not None else self.time_step

        # Convert concentration dict to array
        conc_array = np.zeros(self.mechanism.get_num_species())
        for species, value in concentrations.items():
            idx = self.mechanism.get_species_index(species)
            if idx >= 0:
                conc_array[idx] = value

        # Perform Backward Euler step
        new_conc_array = self._backward_euler_step(conc_array, dt, temperature, pressure, photolysis_rates)

        # Convert back to dictionary
        result = {}
        for idx, species in enumerate(self.mechanism.species):
            result[species] = new_conc_array[idx]

        return result

    def solve(self, concentrations: Dict[str, float],
              temperature: float, pressure: float,
              photolysis_rates: Dict[str, float],
              total_time: float,
              time_step: Optional[float] = None,
              output_steps: Optional[int] = None) -> Dict[str, float]:
        """
        Solve chemistry over multiple time steps.

        Args:
            concentrations: Initial species concentrations (molecules/cm^3)
            temperature: Temperature in K
            pressure: Pressure in Pa
            photolysis_rates: Dictionary of photolysis rates
            total_time: Total integration time in seconds
            time_step: Time step in seconds (uses default if None)
            output_steps: Number of intermediate outputs (None for final only)

        Returns:
            Dictionary of final concentrations (or list if output_steps specified)
        """
        dt = time_step if time_step is not None else self.time_step
        n_steps = int(np.ceil(total_time / dt))

        current_conc = concentrations.copy()
        outputs = []

        for step in range(n_steps):
            # Adjust last time step if needed
            current_dt = dt if step < n_steps - 1 else total_time - step * dt

            current_conc = self.solve_step(current_conc, temperature, pressure,
                                          photolysis_rates, current_dt)

            # Store output if requested
            if output_steps and (step + 1) % (n_steps // output_steps) == 0:
                outputs.append(current_conc.copy())

        if output_steps:
            return outputs
        return current_conc

    def solve_grid(self, concentrations_grid: List[Dict[str, float]],
                   temperatures: Union[float, List[float]],
                   pressures: Union[float, List[float]],
                   photolysis_rates_grid: Union[Dict[str, float], List[Dict[str, float]]],
                   total_time: float,
                   time_step: Optional[float] = None,
                   parallel: bool = True,
                   n_processes: Optional[int] = None) -> List[Dict[str, float]]:
        """
        Solve chemistry for multiple grid cells.

        Args:
            concentrations_grid: List of concentration dictionaries for each grid cell
            temperatures: Temperature(s) in K (scalar or list)
            pressures: Pressure(s) in Pa (scalar or list)
            photolysis_rates_grid: Photolysis rates (single dict or list of dicts)
            total_time: Total integration time in seconds
            time_step: Time step in seconds (uses default if None)
            parallel: Use parallel processing
            n_processes: Number of processes (None for auto)

        Returns:
            List of final concentration dictionaries for each grid cell
        """
        n_cells = len(concentrations_grid)

        # Broadcast scalar inputs to lists
        if isinstance(temperatures, (int, float)):
            temperatures = [temperatures] * n_cells
        if isinstance(pressures, (int, float)):
            pressures = [pressures] * n_cells
        if isinstance(photolysis_rates_grid, dict):
            photolysis_rates_grid = [photolysis_rates_grid] * n_cells

        if parallel:
            # Use multiprocessing for parallel execution
            if n_processes is None:
                n_processes = min(cpu_count(), n_cells)

            # Prepare arguments for parallel processing
            args_list = [
                (concentrations_grid[i], temperatures[i], pressures[i],
                 photolysis_rates_grid[i], total_time, time_step)
                for i in range(n_cells)
            ]

            with Pool(processes=n_processes) as pool:
                results = pool.starmap(self._solve_single_cell, args_list)

            return results
        else:
            # Sequential execution
            results = []
            for i in range(n_cells):
                result = self.solve(
                    concentrations_grid[i],
                    temperatures[i],
                    pressures[i],
                    photolysis_rates_grid[i],
                    total_time,
                    time_step
                )
                results.append(result)
            return results

    def _solve_single_cell(self, concentrations: Dict[str, float],
                          temperature: float, pressure: float,
                          photolysis_rates: Dict[str, float],
                          total_time: float,
                          time_step: Optional[float]) -> Dict[str, float]:
        """
        Helper method for parallel execution.

        Args:
            concentrations: Initial concentrations
            temperature: Temperature in K
            pressure: Pressure in Pa
            photolysis_rates: Photolysis rates
            total_time: Total time
            time_step: Time step

        Returns:
            Final concentrations
        """
        return self.solve(concentrations, temperature, pressure, photolysis_rates,
                         total_time, time_step)

    def get_species_list(self) -> List[str]:
        """Get list of all species in the mechanism"""
        return self.mechanism.species.copy()

    def get_num_species(self) -> int:
        """Get number of species"""
        return self.mechanism.get_num_species()

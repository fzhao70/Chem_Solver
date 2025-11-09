#!/usr/bin/env python3
"""
Example: Multiple Grid Cells with Parallel Computing

Demonstrates solving chemistry for multiple grid cells using parallel processing.
"""

import sys
import os
import time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from chem_solver import ChemicalSolver


def main():
    print("=" * 60)
    print("SMVGEAR-like Solver: Multi-Grid Parallel Example")
    print("=" * 60)

    # Initialize solver
    mechanism_file = os.path.join(os.path.dirname(__file__), '..', 'mechanisms', 'simple_chemistry.json')
    solver = ChemicalSolver(mechanism_file, time_step=60.0)

    # Create multiple grid cells with varying conditions
    n_cells = 100
    print(f"\nNumber of grid cells: {n_cells}")

    # Generate initial concentrations for each grid cell
    # Vary concentrations across grid
    concentrations_grid = []
    for i in range(n_cells):
        # Add some spatial variation
        factor = 0.5 + 1.0 * (i / n_cells)
        concentrations_grid.append({
            'O3': 1.0e12 * factor,
            'NO': 1.0e11 * factor,
            'NO2': 5.0e11 * factor,
            'O1D': 0.0,
            'O3P': 0.0,
            'M': 2.5e19
        })

    # Temperature varies with grid cell (e.g., altitude or latitude)
    temperatures = [298.0 - 10.0 * (i / n_cells) for i in range(n_cells)]

    # Pressure varies with grid cell
    pressures = [101325.0 * (1.0 - 0.2 * i / n_cells) for i in range(n_cells)]

    # Photolysis rates (same for all cells in this example)
    photolysis_rates = {
        'JNO2': 0.008,
        'JO3_O1D': 2.0e-5
    }

    total_time = 3600.0  # 1 hour

    print(f"Temperature range: {min(temperatures):.1f} - {max(temperatures):.1f} K")
    print(f"Pressure range: {min(pressures):.0f} - {max(pressures):.0f} Pa")
    print(f"Integration time: {total_time / 60.0} minutes")

    # Solve with parallel processing
    print("\n" + "-" * 60)
    print("Running with PARALLEL processing...")
    print("-" * 60)
    start_time = time.time()

    results_parallel = solver.solve_grid(
        concentrations_grid,
        temperatures,
        pressures,
        photolysis_rates,
        total_time=total_time,
        parallel=True
    )

    parallel_time = time.time() - start_time
    print(f"Parallel execution time: {parallel_time:.3f} seconds")

    # Solve without parallel processing for comparison
    print("\n" + "-" * 60)
    print("Running with SEQUENTIAL processing...")
    print("-" * 60)
    start_time = time.time()

    results_sequential = solver.solve_grid(
        concentrations_grid,
        temperatures,
        pressures,
        photolysis_rates,
        total_time=total_time,
        parallel=False
    )

    sequential_time = time.time() - start_time
    print(f"Sequential execution time: {sequential_time:.3f} seconds")

    # Compare results
    print("\n" + "=" * 60)
    print("Performance Comparison")
    print("=" * 60)
    print(f"Parallel time:   {parallel_time:.3f} s")
    print(f"Sequential time: {sequential_time:.3f} s")
    print(f"Speedup factor:  {sequential_time / parallel_time:.2f}x")

    # Show sample results
    print("\n" + "=" * 60)
    print("Sample Results (first 5 cells)")
    print("=" * 60)
    for i in range(min(5, n_cells)):
        print(f"\nGrid Cell {i}:")
        print(f"  Temperature: {temperatures[i]:.1f} K")
        print(f"  Initial O3: {concentrations_grid[i]['O3']:.2e} molecules/cm³")
        print(f"  Final O3:   {results_parallel[i]['O3']:.2e} molecules/cm³")
        change = (results_parallel[i]['O3'] - concentrations_grid[i]['O3']) / concentrations_grid[i]['O3'] * 100
        print(f"  Change:     {change:+.2f}%")

    print("\n" + "=" * 60)


if __name__ == '__main__':
    main()

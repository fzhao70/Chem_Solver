#!/usr/bin/env python3
"""
Basic tests for the chemical solver
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
from chem_solver import ChemicalSolver


def test_single_step_solver():
    """Test single step solver"""
    print("Testing single step solver...")

    mechanism_file = os.path.join(os.path.dirname(__file__), '..', 'mechanisms', 'simple_chemistry.json')
    solver = ChemicalSolver(mechanism_file, time_step=60.0)

    concentrations = {
        'O3': 1.0e12,
        'NO': 1.0e11,
        'NO2': 5.0e11,
        'O1D': 0.0,
        'O3P': 0.0,
        'M': 2.5e19
    }

    temperature = 298.0
    pressure = 101325.0
    photolysis_rates = {'JNO2': 0.008, 'JO3_O1D': 2.0e-5}

    result = solver.solve_step(concentrations, temperature, pressure, photolysis_rates)

    # Check that result is a dictionary
    assert isinstance(result, dict), "Result should be a dictionary"

    # Check that all species are present
    for species in concentrations.keys():
        assert species in result, f"Species {species} missing from result"

    # Check that concentrations are non-negative
    for species, conc in result.items():
        assert conc >= 0, f"Concentration of {species} is negative"

    print("  ✓ Single step solver test passed")
    return True


def test_multi_step_solver():
    """Test multi-step solver"""
    print("Testing multi-step solver...")

    mechanism_file = os.path.join(os.path.dirname(__file__), '..', 'mechanisms', 'simple_chemistry.json')
    solver = ChemicalSolver(mechanism_file, time_step=60.0)

    concentrations = {
        'O3': 1.0e12,
        'NO': 1.0e11,
        'NO2': 5.0e11,
        'O1D': 0.0,
        'O3P': 0.0,
        'M': 2.5e19
    }

    temperature = 298.0
    pressure = 101325.0
    photolysis_rates = {'JNO2': 0.008, 'JO3_O1D': 2.0e-5}

    result = solver.solve(
        concentrations,
        temperature,
        pressure,
        photolysis_rates,
        total_time=3600.0
    )

    # Check that result is a dictionary
    assert isinstance(result, dict), "Result should be a dictionary"

    # Check that all species are present
    for species in concentrations.keys():
        assert species in result, f"Species {species} missing from result"

    # Check that concentrations are non-negative
    for species, conc in result.items():
        assert conc >= 0, f"Concentration of {species} is negative"

    print("  ✓ Multi-step solver test passed")
    return True


def test_grid_solver():
    """Test grid solver"""
    print("Testing grid solver...")

    mechanism_file = os.path.join(os.path.dirname(__file__), '..', 'mechanisms', 'simple_chemistry.json')
    solver = ChemicalSolver(mechanism_file, time_step=60.0)

    concentrations = {
        'O3': 1.0e12,
        'NO': 1.0e11,
        'NO2': 5.0e11,
        'O1D': 0.0,
        'O3P': 0.0,
        'M': 2.5e19
    }

    n_cells = 10
    concentrations_grid = [concentrations.copy() for _ in range(n_cells)]
    temperatures = [298.0] * n_cells
    pressures = [101325.0] * n_cells
    photolysis_rates = {'JNO2': 0.008, 'JO3_O1D': 2.0e-5}

    # Test sequential
    results_seq = solver.solve_grid(
        concentrations_grid,
        temperatures,
        pressures,
        photolysis_rates,
        total_time=3600.0,
        parallel=False
    )

    assert len(results_seq) == n_cells, f"Expected {n_cells} results, got {len(results_seq)}"

    # Test parallel
    results_par = solver.solve_grid(
        concentrations_grid,
        temperatures,
        pressures,
        photolysis_rates,
        total_time=3600.0,
        parallel=True
    )

    assert len(results_par) == n_cells, f"Expected {n_cells} results, got {len(results_par)}"

    # Results should be similar (within tolerance due to numerical precision)
    for i in range(n_cells):
        for species in concentrations.keys():
            diff = abs(results_seq[i][species] - results_par[i][species])
            rel_diff = diff / max(results_seq[i][species], 1e-30)
            assert rel_diff < 1e-6, f"Parallel and sequential results differ for {species} in cell {i}"

    print("  ✓ Grid solver test passed")
    return True


def test_mechanism_loading():
    """Test mechanism loading"""
    print("Testing mechanism loading...")

    mechanism_file = os.path.join(os.path.dirname(__file__), '..', 'mechanisms', 'simple_chemistry.json')
    solver = ChemicalSolver(mechanism_file, time_step=60.0)

    species = solver.get_species_list()
    assert len(species) > 0, "No species loaded"

    num_species = solver.get_num_species()
    assert num_species == len(species), "Species count mismatch"

    print(f"  ✓ Loaded {num_species} species")
    return True


def main():
    """Run all tests"""
    print("=" * 60)
    print("Running Chemical Solver Tests")
    print("=" * 60)
    print()

    tests = [
        test_mechanism_loading,
        test_single_step_solver,
        test_multi_step_solver,
        test_grid_solver,
    ]

    passed = 0
    failed = 0

    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"  ✗ Test failed: {e}")
            failed += 1

    print()
    print("=" * 60)
    print(f"Tests passed: {passed}/{len(tests)}")
    if failed > 0:
        print(f"Tests failed: {failed}/{len(tests)}")
        return 1
    else:
        print("All tests passed!")
        return 0


if __name__ == '__main__':
    sys.exit(main())

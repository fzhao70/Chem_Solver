# Chem_Solver

A SMVGEAR-like chemical kinetics solver for atmospheric chemistry simulations.

## Overview

Chem_Solver is a Python implementation of a chemical kinetics solver inspired by SMVGEAR (Sparse-Matrix Vectorized Gear), designed for solving stiff atmospheric chemistry problems. It uses implicit integration methods (Backward Euler with Newton iterations) to efficiently handle the wide range of timescales present in atmospheric chemical mechanisms.

## Features

- **Two API Modes**:
  - **One-step solver**: `solve_step()` - Advances chemistry by a single time step
  - **Multi-step solver**: `solve()` - Integrates chemistry over multiple time steps

- **Flexible Grid Support**:
  - **Single grid cell**: Solve chemistry for a single location
  - **Multiple grid cells**: Solve chemistry for multiple locations with parallel computing support

- **Custom Chemistry Mechanisms**:
  - Define chemistry using JSON mechanism files
  - Support for Arrhenius rate constants
  - Support for photolysis reactions
  - Automatic Jacobian computation

- **Performance**:
  - Implicit Backward Euler method for stiff systems
  - Newton iteration for solving nonlinear systems
  - Parallel processing for multi-grid simulations using multiprocessing

## Installation

### From source

```bash
git clone https://github.com/yourusername/Chem_Solver.git
cd Chem_Solver
pip install -r requirements.txt
pip install -e .
```

## Quick Start

### Single Step Solver

```python
from chem_solver import ChemicalSolver

# Initialize solver
solver = ChemicalSolver('mechanisms/simple_chemistry.json', time_step=60.0)

# Set initial conditions
concentrations = {
    'O3': 1.0e12,
    'NO': 1.0e11,
    'NO2': 5.0e11,
    'M': 2.5e19
}

temperature = 298.0  # K
pressure = 101325.0  # Pa
photolysis_rates = {'JNO2': 0.008, 'JO3_O1D': 2.0e-5}

# Perform one step
result = solver.solve_step(concentrations, temperature, pressure, photolysis_rates)
```

### Multi-Step Solver

```python
# Solve over multiple steps (e.g., 1 hour)
final_concentrations = solver.solve(
    concentrations,
    temperature,
    pressure,
    photolysis_rates,
    total_time=3600.0  # seconds
)
```

### Multiple Grid Cells with Parallel Computing

```python
# Create concentration arrays for multiple grid cells
concentrations_grid = [concentrations.copy() for _ in range(100)]
temperatures = [298.0] * 100
pressures = [101325.0] * 100

# Solve in parallel
results = solver.solve_grid(
    concentrations_grid,
    temperatures,
    pressures,
    photolysis_rates,
    total_time=3600.0,
    parallel=True
)
```

## Chemistry Mechanism Format

Chemistry mechanisms are defined in JSON format:

```json
{
  "name": "Mechanism Name",
  "description": "Description",
  "species": ["O3", "NO", "NO2", "O1D"],
  "reactions": [
    {
      "name": "NO2 photolysis",
      "reactants": {"NO2": 1},
      "products": {"NO": 1, "O3P": 1},
      "rate": {
        "photolysis": "JNO2"
      }
    },
    {
      "name": "NO + O3 -> NO2",
      "reactants": {"NO": 1, "O3": 1},
      "products": {"NO2": 1},
      "rate": {
        "A": 3.0e-12,
        "n": 0.0,
        "Ea": 1500.0
      }
    }
  ]
}
```

### Rate Constant Formats

1. **Constant rate**: `{"constant": 1.0e-12}`
2. **Arrhenius**: `{"A": 3.0e-12, "n": 0.0, "Ea": 1500.0}`
   - Rate = A × T^n × exp(-Ea/T)
3. **Photolysis**: `{"photolysis": "JNO2"}`
   - Rate provided via photolysis_rates dictionary

## Examples

The `examples/` directory contains several complete examples:

1. **example_single_step.py** - Demonstrates one-step solver
2. **example_multi_step.py** - Demonstrates multi-step solver
3. **example_single_grid.py** - Single grid cell simulation
4. **example_multi_grid_parallel.py** - Multi-grid with parallel computing

Run examples:

```bash
python examples/example_single_step.py
python examples/example_multi_step.py
python examples/example_multi_grid_parallel.py
```

## API Reference

### ChemicalSolver Class

#### Initialization

```python
ChemicalSolver(mechanism_file, time_step, rtol=1e-3, atol=1e-10, max_iterations=100)
```

**Parameters:**
- `mechanism_file` (str): Path to JSON mechanism file
- `time_step` (float): Default time step in seconds
- `rtol` (float): Relative tolerance for Newton convergence
- `atol` (float): Absolute tolerance for Newton convergence
- `max_iterations` (int): Maximum Newton iterations

#### Methods

##### solve_step()

Perform a single integration step.

```python
solve_step(concentrations, temperature, pressure, photolysis_rates, time_step=None)
```

**Parameters:**
- `concentrations` (dict): Species concentrations (molecules/cm³)
- `temperature` (float): Temperature in K
- `pressure` (float): Pressure in Pa
- `photolysis_rates` (dict): Photolysis rates (1/s)
- `time_step` (float, optional): Time step (uses default if None)

**Returns:** Dictionary of updated concentrations

##### solve()

Solve chemistry over multiple time steps.

```python
solve(concentrations, temperature, pressure, photolysis_rates, total_time, time_step=None)
```

**Parameters:**
- `concentrations` (dict): Initial concentrations (molecules/cm³)
- `temperature` (float): Temperature in K
- `pressure` (float): Pressure in Pa
- `photolysis_rates` (dict): Photolysis rates (1/s)
- `total_time` (float): Total integration time in seconds
- `time_step` (float, optional): Time step (uses default if None)

**Returns:** Dictionary of final concentrations

##### solve_grid()

Solve chemistry for multiple grid cells.

```python
solve_grid(concentrations_grid, temperatures, pressures, photolysis_rates_grid,
          total_time, time_step=None, parallel=True, n_processes=None)
```

**Parameters:**
- `concentrations_grid` (list): List of concentration dictionaries
- `temperatures` (float or list): Temperature(s) in K
- `pressures` (float or list): Pressure(s) in Pa
- `photolysis_rates_grid` (dict or list): Photolysis rates
- `total_time` (float): Total integration time
- `time_step` (float, optional): Time step
- `parallel` (bool): Use parallel processing
- `n_processes` (int, optional): Number of processes (auto-detect if None)

**Returns:** List of final concentration dictionaries

## Algorithm Details

The solver uses a Backward Euler method with Newton iterations:

1. **Backward Euler**: C(t+Δt) = C(t) + Δt × f(C(t+Δt))
2. **Newton Iteration**: Solves the nonlinear system iteratively
3. **Jacobian**: Computed numerically using finite differences
4. **Convergence**: Checks relative error against tolerance

This implicit method is well-suited for stiff chemical systems where reactions occur on vastly different timescales.

## Performance Tips

1. **Time step selection**: Choose time steps appropriate for your chemistry
   - Faster chemistry requires smaller time steps
   - Typical range: 1-300 seconds

2. **Parallel processing**: Enable for multi-grid simulations
   - Automatically uses all available CPU cores
   - Most efficient for 10+ grid cells

3. **Tolerance settings**: Adjust for accuracy vs. speed
   - `rtol=1e-3`: Good balance for most applications
   - Decrease for higher accuracy (slower)
   - Increase for faster (less accurate) solutions

## Comparison to SMVGEAR

This implementation is inspired by SMVGEAR but uses simplified approaches:

**Similarities:**
- Implicit integration for stiff systems
- Support for custom chemistry mechanisms
- Multi-grid parallel processing

**Differences:**
- Uses Backward Euler instead of Gear's method (BDF)
- Dense matrix operations instead of sparse matrices
- Numerical Jacobian instead of analytical
- Python implementation instead of Fortran

For production atmospheric models requiring maximum performance, consider using the original SMVGEAR or modern alternatives like KPP (Kinetic PreProcessor).

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use this solver in your research, please cite:

```
Chem_Solver: A SMVGEAR-like Chemical Kinetics Solver
https://github.com/yourusername/Chem_Solver
```

## References

1. Jacobson, M. Z., & Turco, R. P. (1994). SMVGEAR: A sparse-matrix, vectorized Gear code for atmospheric models. Atmospheric Environment, 28(2), 273-284.

2. Sandu, A., et al. (2006). Technical note: Simulating chemical systems in Fortran90 and Matlab with the Kinetic PreProcessor KPP-2.1. Atmospheric Chemistry and Physics, 6(1), 187-195.

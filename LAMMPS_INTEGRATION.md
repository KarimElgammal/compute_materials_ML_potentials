# LAMMPS Integration Guide for ML Potentials

This guide provides detailed instructions for using Machine Learning Potentials (MLPs) with LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator).

## Prerequisites

- LAMMPS installation (latest stable version recommended)
- Python 3.11
- ML Potential packages (as needed: ORB, CHGNet, MatGL, MACE, or SevenNet)

## Installation Steps

### 1. ORB Models Integration

```bash
# Clone the ORB-LAMMPS patch repository
git clone https://github.com/stefanbringuier/ORB-LAMMPS-PATCH
cd ORB-LAMMPS-PATCH

# Install ORB models if not already installed
pip install orb-models

# Copy the pair_orb files to your LAMMPS src directory
cp pair_orb.* /path/to/lammps/src/

# Recompile LAMMPS with the new pair style
cd /path/to/lammps/src
make yes-python    # Enable Python support
make serial       # or make mpi for parallel support
```

### 2. Other ML Potentials Integration

Most ML potentials can be used through the LAMMPS Python interface or dedicated pair styles.

## Available ML Packages in LAMMPS

LAMMPS supports several ML potential packages that can be integrated:

1. **MLMOD-PYTORCH**
   - Machine learning methods for data-driven modeling with PyTorch
   - Provides time-step integrators for dynamics and interactions
   - Supports Neural Networks, Kernel Regression, and other ML models
   - Can import models trained in PyTorch or other ML frameworks
   ```bash
   # Installation via pip
   pip install mlmod-pytorch
   ```

2. **USER-MLIP**
   - Supports Moment Tensor Potential (MTP)
   - Handles linearized ML potentials
   - Efficient implementation for large-scale simulations

3. **USER-DEEPMD**
   - Integration with DeePMD-kit ML potentials
   - Can be compiled as a plugin for existing LAMMPS installations
   ```bash
   # Installation as plugin
   cmake -D PKG_USER-DEEPMD=yes ../cmake
   ```

4. **USER-AENET**
   - Atomic Energy Network potentials using artificial neural networks
   - Flexible conditions using conventional pair_style formats
   - Compatible with EAM-like syntax
   ```bash
   git clone https://github.com/HidekiMori-CIT/aenet-lammps
   ```

## Usage Examples

### ORB Models

```lammps
# LAMMPS input script for ORB models
units metal
atom_style atomic

# Read your structure
read_data your_structure.data

# Set up ORB potential
pair_style orb /path/to/orb_driver.py
pair_coeff * * orb-d3-v2 Si O  # Example for SiO2

# Run your simulation
minimize 1.0e-4 1.0e-6 100 1000
```

### CHGNet Example

```lammps
# LAMMPS input script for CHGNet
units metal
atom_style atomic

# Read your structure
read_data your_structure.data

# Set up CHGNet potential
pair_style mlip_model
pair_coeff * * /path/to/chgnet_model.pb

# Run your simulation
minimize 1.0e-4 1.0e-6 100 1000
```

## Performance Tips

1. **GPU Acceleration**
   ```lammps
   # Enable GPU acceleration when available
   package gpu 1
   pair_style orb /path/to/orb_driver.py gpu
   ```

2. **Parallel Processing**
   ```bash
   # Run with MPI
   mpirun -np 4 lmp -in input.lammps
   ```

3. **Model Selection**
   - Use smaller model variants for faster calculations
   - Consider `orb-d3-sm-v2` or `orb-d3-xs-v2` for large-scale simulations

4. **Plugin-based Integration**
   ```bash
   # Load ML potential as plugin at runtime
   plugin load libdeepmd_lmp.so
   ```

5. **Package Selection**
   - Choose the appropriate ML package based on your needs:
     - MLMOD-PYTORCH: For PyTorch-based models and custom implementations
     - USER-MLIP: For moment tensor potentials
     - USER-DEEPMD: For DeePMD models
     - USER-AENET: For ANN-based potentials

## Common Issues and Solutions

### Memory Issues
```lammps
# Adjust these settings if you encounter memory problems
package gpu 1 neigh no split 2.0 gpuID 0
```

### Performance Issues
- Ensure proper GPU drivers are installed
- Monitor memory usage with `nvidia-smi`
- Use appropriate cutoff distances

## Advanced Features

### Energy and Force Analysis
```lammps
# Output energy and forces
compute pe all pe
compute forces all property/atom fx fy fz
dump forces all custom 100 forces.dump id type x y z fx fy fz
```

### Temperature Control
```lammps
# NVT ensemble example
fix nvt all nvt temp 300.0 300.0 100.0
```

## Validation and Testing

1. **Energy Conservation**
```lammps
# Check energy conservation in NVE
fix nve all nve
thermo_style custom step temp pe ke etotal
thermo 100
```

2. **Force Validation**
```python
# Python script to compare forces with DFT
from ase.io import read
from ase.calculators.lammpsrun import LAMMPS
# ... validation code ...
```

## References

1. [LAMMPS external packages documentation](https://www.lammps.org/external.html)
2. [ORB-LAMMPS-PATCH](https://github.com/stefanbringuier/ORB-LAMMPS-PATCH)
3. [Advance Soft Corp LAMMPS](https://github.com/advancesoftcorp/lammps)
4. Yang, Y., Zhang, S., Ranasinghe, K. D., Isayev, O., & Roitberg, A. E. (2024). Machine Learning of Reactive Potentials. Annual Review of Physical Chemistry, 75, 371-395. [DOI](https://doi.org/10.1146/annurev-physchem-062123-024417)

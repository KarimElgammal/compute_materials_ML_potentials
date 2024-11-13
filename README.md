# Materials ML Potentials Guide

A comprehensive guide for setting up and using various Machine Learning Potentials for materials science calculations. This guide covers installation and setup across different platforms.

## Prerequisites
- Python 3.11 (Python 3.12 not supported)
- pip (package installer for Python)
- git

## ML Models discussed with examples here

- **CHGNet**: Universal deep learning potential for materials
- **MatGL**: Graph Learning for Materials
- **MACE**: Message Passing Neural Networks
- **SevenNet**: Deep learning potential for materials
- **ORB Models**: Pretrained models for atomic simulations

## Installation Guide

```bash
# 1. Set up Python environment
pyenv install 3.11
pyenv virtualenv 3.11 MLpotentials
pyenv activate MLpotentials

# 2. Clean existing installations
pip install --upgrade pip
pip uninstall -y torch torchvision torchaudio lightning pytorch-lightning chgnet matgl dgl

# 3. Install PyTorch ecosystem (M1-specific versions)
pip install torch==2.1.0 torchvision==0.16.0 torchaudio==2.1.0

# 4. Install DGL (CPU version for M1)
pip install dgl -f https://data.dgl.ai/wheels/cpu/repo.html

# 5. Install base dependencies
pip install numpy==1.26.4 scipy scikit-learn==1.3.1 pandas
pip install mp-api matplotlib tqdm joblib==1.4.2

# 6. Install MatCalc and its dependencies
git clone https://github.com/materialsvirtuallab/matcalc.git
cd matcalc
pip install -r requirements.txt
pip install -e .

# 7. Install ML Potential Models
pip install git+https://github.com/CederGroupHub/chgnet  # CHGNet
pip install git+https://github.com/materialsvirtuallab/matgl.git  # MatGL
pip install git+https://github.com/ACEsuit/mace.git  # MACE
pip install git+https://github.com/MDIL-SNU/SevenNet.git  # SevenNet

# 8. Install ORB Models (I advice if you install it on different )
pip install orb-models
pip install "pynanoflann@git+https://github.com/dwastberg/pynanoflann#egg=af434039ae14bedcbb838a7808924d6689274168",
```

## example using ORB models
```python
from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator
# initialize model
orbff = pretrained.orb_d3_v2(device="cpu") # or device="cuda"
calc = ORBCalculator(orbff, device="cpu")
#Set up structure and calculate
atoms.set_calculator(calc)
energy = atoms.get_potential_energy()
forces = atoms.get_forces()
```


## Available ORB Models
- **orb-v2**: Trained on MPTraj + Alexandria
- **orb-mptraj-only-v2**: Trained on MPTraj dataset only
- **orb-d3-v2**: Trained with integrated D3 corrections (recommended)
- **orb-d3-{sm,xs}-v2**: Smaller versions of orb-d3-v2

## Materials Project API Setup
1. Register at [Materials Project](https://materialsproject.org/)
2. Get your API key from your dashboard
3. Create `mp_api_key.txt` in your working directory
4. Add your API key to the file

## Common Issues and Solutions

### M1 Mac Issues
- If DGL fails: Try reinstalling with CPU-only version
- PyTorch issues: Stick to versions specified above
- Memory errors: Reduce batch sizes in calculations

### Linux/Windows Issues
- CUDA version conflicts: Match PyTorch and DGL CUDA versions
- DGL installation fails: Try CPU version first
- ORB Models: Windows support not guaranteed

### General Troubleshooting
1. Always use Python 3.11
2. Create fresh virtual environment if conflicts occur
3. Install packages in the order specified
4. Check GPU compatibility if using CUDA versions

## References
- [MatCalc Documentation](https://materialsvirtuallab.github.io/matcalc)
- [CHGNet Paper](https://www.nature.com/articles/s43588-022-00349-3)
- [Materials Project](https://materialsproject.org/)
- [ORB Models Repository](https://github.com/ur-whitelab/orb-models)


### For Linux Users
```bash
# 1. Set up Python environment
python -m venv MLpotentials
source MLpotentials/bin/activate

# 2. Install PyTorch with CUDA support (if GPU available)
pip install torch torchvision torchaudio

# 3. Install DGL with CUDA support
pip install dgl -f https://data.dgl.ai/wheels/cu117/repo.html

# 4. Follow steps 5-7 from M1 Mac installation
```

### For Windows Users
```bash
# 1. Set up Python environment
python -m venv MLpotentials
MLpotentials\Scripts\activate

# 2. Install PyTorch
pip install torch torchvision torchaudio

# 3. Install DGL
pip install dgl -f https://data.dgl.ai/wheels/cu117/repo.html

# 4. Follow steps 5-7 from M1 Mac installation
```

## Materials Project API Setup
1. Register at [Materials Project](https://materialsproject.org/)
2. Get your API key from your dashboard
3. Create `mp_api_key.txt` in your working directory
4. Add your API key to the file

## Package Overview

### Core Packages
- **MatCalc**: Main framework for materials calculations
- **PyTorch**: Deep learning framework
- **DGL**: Deep Graph Library for graph neural networks

### ML Potential Models
- **CHGNet**: Universal deep learning potential for materials
- **MatGL**: Graph Learning for Materials
- **MACE**: Message Passing Neural Networks
- **SevenNet**: Deep learning potential for materials

### Supporting Libraries
- **NumPy**: Numerical computing
- **SciPy**: Scientific computing
- **Pandas**: Data manipulation
- **mp-api**: Materials Project API client

## Verification
Test your installation:
```python
# Test MatCalc and ML models
from matcalc.utils import get_universal_calculator
models = [(name, get_universal_calculator(name)) 
          for name in ("M3GNet", "CHGNet", "MACE", "SevenNet")]

# Test Materials Project API
with open('mp_api_key.txt', 'r') as f:
    api_key = f.read().strip()
from mp_api.client import MPRester
mpr = MPRester(api_key)
```

## Common Issues and Solutions

### M1 Mac Issues
- If DGL fails: Try reinstalling with CPU-only version
- PyTorch issues: Stick to versions specified above
- Memory errors: Reduce batch sizes in calculations

### Linux/Windows Issues
- CUDA version conflicts: Match PyTorch and DGL CUDA versions
- DGL installation fails: Try CPU version first

### General Troubleshooting
1. Always use Python 3.11
2. Create fresh virtual environment if conflicts occur
3. Install packages in the order specified
4. Check GPU compatibility if using CUDA versions

## Usage Examples
```python
# Basic structure optimization example
from matcalc.utils import get_universal_calculator
from pymatgen.core import Structure

# Load calculator
calc = get_universal_calculator("CHGNet")

# Load structure
structure = Structure.from_file("your_structure.cif")

# Calculate properties
energy = calc.get_potential_energy(structure)
forces = calc.get_forces(structure)
```

## References
- [MatCalc Documentation](https://materialsvirtuallab.github.io/matcalc)
- [CHGNet Paper](https://www.nature.com/articles/s43588-022-00349-3)
- [Materials Project](https://materialsproject.org/)

## License
This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.

I made sure that '6_use_matCalc/matcalc/examples/Calculating MLIP properties.ipynb' is working by setting the MP_API_KEY in the environment. matcalc is a great package for computing materials properties using machine learning potentials.
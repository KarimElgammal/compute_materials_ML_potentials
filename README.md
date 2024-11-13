# Compute Materials at atomistic scale using ML Potentials

This repository contains tools and scripts for computing materials properties using machine learning potentials. I gave this seminar talk at the American University in Cairo in November 2024, and this repository contains the code for the talk. The talk is about using machine learning potentials to compute materials properties at the atomistic scale. It was invited talk as part of course given by [Dr. Mostafa Youssef](https://sites.google.com/aucegypt.edu/materialstheorygroup/home).

## Repository Structure

### 0_system_setup/
- Setup scripts and initial configuration files for system preparation
- Contains base structure files and input parameters

### 1_geometry_optimisation/
- Scripts for performing geometry optimization calculations
- Optimizes atomic positions to find minimum energy configurations
- Outputs energy and structure files at each optimization step

### 3_substrate_energy/
- Calculations for isolated substrate systems
- Computes reference energies for substrate surfaces

### 4_adsorbate_energy/
- Calculations for isolated adsorbate molecules
- Determines reference energies for adsorbate species

### 5_calc_binding_energy/
- Scripts for computing binding energies
- Combines results from substrate and adsorbate calculations
- Outputs final binding energy analysis

### 7_example_with_Cu/
- Complete example calculation for Cu(111) surface using orbModel
- Includes:
  - System setup
  - Geometry optimization
  - Energy calculations
  - Results analysis

### 99_seminar_slides/
- Presentation slides for the talk

## Key Features

- Automated binding energy calculations
- Geometry optimization with force convergence
- Support for periodic boundary conditions
- Analysis tools for energy and structure data
- Example calculations with Cu(111) surface

## Dependencies

### OrbModels
Install OrbModels and its dependencies:
```bash
pip install orb-models
pip install "pynanoflann@git+https://github.com/dwastberg/pynanoflann#egg=af434039ae14bedcbb838a7808924d6689274168"
```

### MatCalc
Install MatCalc using:
```bash
pip install git+https://github.com/materialsvirtuallab/matcalc.git
```
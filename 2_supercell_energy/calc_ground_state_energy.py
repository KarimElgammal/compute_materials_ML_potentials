import ase
from ase.io import read, write
import numpy as np
from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator
from ase import Atoms
import re
import sys
import argparse
from datetime import datetime
import os

def main():
    parser = argparse.ArgumentParser(description='Calculate ground state energy using ORB.')
    parser.add_argument('input_file', help='Path to the input XYZ file')
    parser.add_argument('--device', default='cpu', choices=['cpu', 'cuda'], 
                       help='Device to run calculations on (default: cpu)')
    
    args = parser.parse_args()

    # Create output directory with timestamp and input file name
    basename = os.path.splitext(os.path.basename(args.input_file))[0]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"orb_energy_{basename}_{timestamp}"
    os.makedirs(output_dir, exist_ok=True)

    # Set up logging
    energy_file = os.path.join(output_dir, "energy.txt")

    def parse_lattice(lattice_str):
        numbers = [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+", lattice_str)]
        return np.array(numbers).reshape(3, 3)

    # Read structure
    try:
        with open(args.input_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: Could not find file '{args.input_file}'")
        sys.exit(1)

    # Parse header
    try:
        header = lines[1]
        lattice_str = re.search(r'Lattice="([^"]+)"', header).group(1)
        pbc_str = re.search(r'pbc="([^"]+)"', header).group(1)
    except (IndexError, AttributeError):
        print("Error: File format incorrect. Expected extended XYZ format with Lattice and pbc information.")
        sys.exit(1)

    cell = parse_lattice(lattice_str)
    pbc = [x.strip().lower() == 't' for x in pbc_str.split()]

    # Parse atoms
    symbols = []
    positions = []
    tags = []

    for line in lines[2:]:
        if line.strip():
            try:
                parts = line.split()
                symbols.append(parts[0])
                positions.append([float(x) for x in parts[1:4]])
                tags.append(int(parts[4]) if len(parts) > 4 else 0)
            except (IndexError, ValueError) as e:
                print(f"Error parsing line: {line.strip()}")
                print(f"Error details: {e}")
                sys.exit(1)

    atoms = Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell,
        pbc=pbc,
        tags=tags
    )

    print(f"Successfully loaded structure with {len(atoms)} atoms")
    print(f"Output directory created: {output_dir}")

    # Initialize model and calculate energy
    print(f"Initializing ORB-D3 model on {args.device}...")
    orbff = pretrained.orb_d3_v2(device=args.device)
    calc = ORBCalculator(orbff, device=args.device)
    atoms.calc = calc

    # Calculate ground state energy
    energy = atoms.get_potential_energy()
    
    # Save energy and system information to file
    with open(energy_file, 'w') as f:
        f.write(f"System: {basename}\n")
        f.write(f"Number of atoms: {len(atoms)}\n")
        f.write(f"Composition: {atoms.symbols}\n")
        f.write(f"Ground state energy: {energy:.6f} eV\n")

    print(f"\nGround state energy: {energy:.6f} eV")
    print(f"Results saved to {energy_file}")

if __name__ == "__main__":
    main()
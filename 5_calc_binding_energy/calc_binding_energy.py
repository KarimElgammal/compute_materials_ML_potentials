import os
import glob
import re

def find_latest_energy_file(directory):
    """Find the most recent energy.txt file in the given directory"""
    pattern = os.path.join(directory, "orb_energy_*", "energy.txt")
    energy_files = glob.glob(pattern)
    
    if not energy_files:
        raise FileNotFoundError(f"No energy.txt files found in {directory}")
    
    # Sort by modification time and get the most recent
    latest_file = max(energy_files, key=os.path.getmtime)
    return latest_file

def extract_energy(filename):
    """Extract the ground state energy from energy.txt file"""
    with open(filename, 'r') as f:
        for line in f:
            if "Ground state energy:" in line:
                energy = float(re.search(r"[-+]?\d*\.\d+", line).group())
                return energy
    raise ValueError(f"Could not find energy value in {filename}")

def main():
    # Define directories
    supercell_dir = "../2_supercell_energy"
    substrate_dir = "../3_substrate_energy"
    adsorbate_dir = "../4_adsorbate_energy"
    
    try:
        # Get the most recent energy files
        supercell_file = find_latest_energy_file(supercell_dir)
        substrate_file = find_latest_energy_file(substrate_dir)
        adsorbate_file = find_latest_energy_file(adsorbate_dir)
        
        # Extract energies
        E_total = extract_energy(supercell_file)
        E_substrate = extract_energy(substrate_file)
        E_adsorbate = extract_energy(adsorbate_file)
        
        E_binding = (E_total - (E_substrate + E_adsorbate))
        
        # Create output directory with results
        os.makedirs("binding_energy_results", exist_ok=True)
        
        # Save results
        with open("binding_energy_results/binding_energy.txt", 'w') as f:
            f.write("Binding Energy Analysis\n")
            f.write("======================\n\n")
            f.write(f"Total system energy:     {E_total:.6f} eV\n")
            f.write(f"Substrate energy:        {E_substrate:.6f} eV\n")
            f.write(f"Adsorbate energy:        {E_adsorbate:.6f} eV\n")
            f.write(f"Binding energy:          {E_binding:.6f} eV\n\n")
            f.write("Files used:\n")
            f.write(f"Total system:    {supercell_file}\n")
            f.write(f"Substrate:       {substrate_file}\n")
            f.write(f"Adsorbate:       {adsorbate_file}\n")
        
        # Print results to console
        print("\nBinding Energy Analysis")
        print("======================")
        print(f"\nTotal system energy:     {E_total:.6f} eV")
        print(f"Substrate energy:        {E_substrate:.6f} eV")
        print(f"Adsorbate energy:        {E_adsorbate:.6f} eV")
        print(f"Binding energy:          {E_binding:.6f} eV\n")
        print(f"Results saved to: binding_energy_results/binding_energy.txt")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Make sure you have run the energy calculations in all three directories.")
    except ValueError as e:
        print(f"Error: {e}")
        print("Make sure the energy.txt files contain valid energy values.")

if __name__ == "__main__":
    main() 
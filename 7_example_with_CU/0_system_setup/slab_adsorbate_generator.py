from ase.build import fcc111
from ase.io import write
from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
import numpy as np

# Dictionary of common metal lattice constants (in Å)
METAL_LATTICE_CONSTANTS = {
    'Al': 4.05,
    'Cu': 3.61,
    'Ag': 4.09,
    'Au': 4.08,
    'Pt': 3.92,
    'Pd': 3.89,
    'Ni': 3.52,
    'Fe': 2.87,
}

def create_slab_with_adsorbate(metal='Cu', surface_index=(1,1,1), sizes=[(4,4)], 
                              adsorbate_smiles='O=C=O', layers=4, vacuum=10.0, 
                              a=None, molecule_height=3.0):
    """
    Create metal surface slabs with an adsorbate molecule on top.
    
    Parameters:
    metal (str): Chemical symbol of the metal (default is 'Cu')
    surface_index (tuple): Miller indices of the surface (default is (1,1,1))
    sizes (list of tuples): List of (x, y) sizes for the slab, e.g. [(2,2), (3,3)]
    adsorbate_smiles (str): SMILES string of the adsorbate molecule (default is CO2)
    layers (int): Number of atomic layers in the slab (default is 4)
    vacuum (float): Vacuum size in Angstroms (default is 10.0)
    a (float): Optional - Lattice constant for the metal in Angstroms. If None, uses value from METAL_LATTICE_CONSTANTS
    molecule_height (float): Height of adsorbate above surface (default is 3.0)
    
    Returns:
    None, but saves .xyz files and .png visualizations
    """
    # Get lattice constant if not specified
    if a is None:
        if metal not in METAL_LATTICE_CONSTANTS:
            raise ValueError(f"Lattice constant for {metal} not found. Please provide 'a' parameter or add to METAL_LATTICE_CONSTANTS")
        a = METAL_LATTICE_CONSTANTS[metal]

    # Generate adsorbate molecule
    mol = Chem.MolFromSmiles(adsorbate_smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    
    # Convert RDKit molecule to ASE Atoms object
    positions = mol.GetConformer().GetPositions()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    adsorbate = Atoms(symbols=symbols, positions=positions)
    
    # Create surface based on Miller indices
    if surface_index == (1,1,1):
        surface_function = fcc111
    # Add more surface types as needed
    # elif surface_index == (1,0,0):
    #     surface_function = fcc100
    
    for size in sizes:
        # Create the slab
        slab = surface_function(metal, size=(size[0], size[1], layers), a=a, vacuum=vacuum)
        
        # Calculate the center of the slab in xy plane
        center_x = np.mean(slab.positions[:, 0])
        center_y = np.mean(slab.positions[:, 1])
        top_z = np.max(slab.positions[:, 2])
        
        # Position the adsorbate molecule
        adsorbate.translate([center_x, center_y, top_z + molecule_height])
        
        # Combine slab and molecule
        combined = slab + adsorbate
        
        # Generate filenames
        base_name = f"{metal}_{surface_index[0]}{surface_index[1]}{surface_index[2]}_{size[0]}x{size[1]}x{layers}"
        xyz_filename = f"{base_name}.xyz"
        png_filename = f"{base_name}.png"
        
        # Save the combined structure as an XYZ file
        write(xyz_filename, combined)
        
        # Create visualizations of the combined structure (top and side views)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        
        # Top view
        plot_atoms(combined, ax1, radii=0.5, rotation=('0x,0y,0z'))
        ax1.set_title('Top View')
        ax1.axis('off')
        
        # Side view
        plot_atoms(combined, ax2, radii=0.5, rotation=('90x,0y,0z'))
        ax2.set_title('Side View')
        ax2.axis('off')
        
        plt.tight_layout()
        plt.savefig(png_filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Created {xyz_filename} with {len(combined)} atoms")
        print(f"Saved visualization as {png_filename}")

# Example usage for CO2 on Cu(111)
if __name__ == "__main__":
    sizes_to_create = [(4,4)]
    create_slab_with_adsorbate(
        metal='Cu', # or 'Al'
        surface_index=(1,1,1),
        sizes=sizes_to_create,
        adsorbate_smiles='O=C=O',  # SMILES for CO2 # CO: 'C#O' H2O: 'O' CH4: 'C' NH3: 'N' CH3OH: 'CO'
        layers=6,
        vacuum=10.0,
        molecule_height=3.0
    )
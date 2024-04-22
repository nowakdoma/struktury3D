# Processing 3D Structures Using the Biopython Package

## Description:
This repository contains two Python projects for analyzing and visualizing protein and RNA structures using the Biopython package. Each project consists of scripts dedicated to specific structural analyses and visualizations.

## Dependencies:
- Python 3.x
- Biopython
- numpy
- matplotlib
- DSSP (for running Ramachandran plot on Linux)

## Protein Contact Map Calculator and Visualizer

### Description:
This Python script calculates and visualizes the contact map of a protein structure based on a provided PDB (Protein Data Bank) file. A contact map is a graphical representation illustrating the proximity of amino acid residues within a protein structure. Residues are considered in contact if the distance between their alpha carbon atoms falls below a specified threshold.

### Features:
- **Input:** Accepts a PDB file containing the 3D coordinates of a protein structure.
- **Chain Selection:** Allows the user to specify a particular chain within the protein structure for analysis.
- **Threshold Setting:** Enables customization of the threshold distance for defining contacts (default: 8.0 angstroms).
- **Contact Map Calculation:** Utilizes numpy to compute the contact map matrix based on pairwise distances between residues.
- **Visualization:** Utilizes matplotlib to create a heatmap visualization of the contact map, highlighting residue contacts.
- **Error Handling:** Provides warnings and error messages for invalid input or missing data.

### Usage:
1. Ensure Python and the required libraries are installed.
2. Run the script from the command line.
3. Provide the path to the input PDB file using the `--f` argument.
4. Optionally, specify the chain ID using the `--c` argument (default: first chain).
5. Optionally, adjust the threshold distance using the `--t` argument (default: 8.0 angstroms).

#### Example:
```
python protein_contact_map.py --f structure.pdb --c A --t 6.0
```

This command calculates and visualizes the contact map for chain 'A' in the 'structure.pdb' file, considering contacts within 6.0 angstroms.

### Output:
A graphical representation of the protein contact map will be displayed, showing residue contacts based on the specified threshold.

## Protein and RNA Structure Analysis

### Description:
This project consists of Python scripts for analyzing and visualizing protein and RNA structures.

### Features:

#### Protein Structure Analysis:
- **Ramachandran Plotting:** The `plotRamachandran` script generates a Ramachandran plot, illustrating the distribution of phi and psi angles for protein residues.
- **Phi-Psi Calculation:** Calculates phi and psi angles for each residue in the protein structure.
- **Secondary Structure Prediction:** Predicts the secondary structure of the protein using DSSP.

#### RNA Structure Analysis:
- **RNA Torsion Angle Calculation:** The `getTorsionAngles` script calculates torsion angles (alpha, beta, gamma, delta, epsilon, zeta, chi) for RNA residues.
- **HETATM Removal:** Removes non-standard residues (HETATM) from the RNA structure before analysis.

### Usage:

#### Protein Structure Analysis:
1. Run the `katytorsyjne` script.
2. Provide the path to the input PDB file using the `--f` argument.
3. Use the `--p` flag to specify protein structure analysis.

#### RNA Structure Analysis:
1. Run the `katytorsyjne` script.
2. Provide the path to the input PDB file using the `--f` argument.
3. Use the `--r` flag to specify RNA structure analysis.

### Example:

#### Protein Structure Analysis:
```
python katytorsyjne.py --f protein_structure.pdb --p
```

#### RNA Structure Analysis:
```
python katytorsyjne.py --f rna_structure.pdb --r
```

### Output:
- For protein structure analysis, the `katytorsyjne` script generates a Ramachandran plot showing the distribution of phi and psi angles.
- For RNA structure analysis, the `katytorsyjne` script calculates torsion angles (alpha, beta, gamma, delta, epsilon, zeta, chi) and saves them to a CSV file.


### Note:
- The script assumes that the input PDB file contains only protein atoms and follows the standard PDB format.
- Ensure proper installation of Biopython, numpy, and matplotlib libraries and DSSP before running the script.

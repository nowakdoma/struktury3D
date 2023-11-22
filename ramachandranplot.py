from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB.Polypeptide import Polypeptide
from Bio.PDB.DSSP import DSSP


def calPhiPsi(structure):
    phiPsiAngles = []

    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            for pp in polypeptides:
                angles = []

                phi_psi_list = Polypeptide(pp).get_phi_psi_list()

                for phi, psi in phi_psi_list:
                    if phi is not None and psi is not None:
                        phi_degrees = np.degrees(phi)
                        psi_degrees = np.degrees(psi)
                        angles.append((phi_degrees, psi_degrees))

                phiPsiAngles.extend(angles)

    return phiPsiAngles


def getSS(pdb_file_path, structure):
    model = structure[0]
    dssp = DSSP(model, pdb_file_path, dssp='/usr/bin/dssp')

    secondaryStructure = []
    for key, value in dssp.property_dict.items():
        #resId = key[1]
        ss = value[2]
        secondaryStructure.append(ss)


    return secondaryStructure


def plotRamachandran(pdb_file_path):
    structure = PDB.PDBParser(QUIET=True).get_structure("protein", pdb_file_path)
    phiPsiAngles = calPhiPsi(structure)
    secondary_structure = getSS(pdb_file_path, structure)

    plt.figure(figsize=(10, 8))

    color_dict = {'H': 'pink', 'G': 'orange', 'I': 'yellow', 'E': 'green', 'B': 'blue', 'T': 'purple', 'S': 'brown', 'c': 'gray'}

    for phi_psi, ss in zip(phiPsiAngles, secondary_structure):
        color = color_dict.get(ss, 'black')  # black for unknown ss
        plt.scatter(phi_psi[0], phi_psi[1], marker='o', alpha=0.6, color=color, label=ss)

    plt.title('Ramachandran Plot')
    plt.xlabel('Phi Angle (degrees)')
    plt.ylabel('Psi Angle (degrees)')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.grid(True)

    legend_labels = {'H': 'α-helix', 'G': '310-helix', 'I': 'π-helix', 'E': 'Extended Strand', 'B': 'Isolated β-bridge', 'T': 'Turn', 'S': 'Bend', 'c': 'Coil'}
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_dict[ss], markersize=10, label=legend_labels[ss]) for ss in color_dict.keys()]
    plt.legend(handles=legend_elements, loc='upper right')

    plt.show()

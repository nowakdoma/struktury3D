from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB.Polypeptide import Polypeptide


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


def plotRamachandran(structure):
    phiPsiAngles = calPhiPsi(structure)

    plt.figure(figsize=(10, 8))

    plt.scatter(*zip(*phiPsiAngles), marker='o', alpha=0.6)

    plt.title('Ramachandran Plot')
    plt.xlabel('Phi Angle (degrees)')
    plt.ylabel('Psi Angle (degrees)')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.grid(True)
    plt.show()


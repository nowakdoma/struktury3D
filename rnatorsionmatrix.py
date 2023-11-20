from Bio import PDB
from Bio.PDB.vectors import calc_dihedral
import numpy as np
import csv


def modifiedResidues(filename):
    modifRes = set()

    with open(filename, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("SEQADV"):
                modifRes.add(line[12:15])
            elif modifRes:
                break

    return list(modifRes)


def removeHetatm(structure, filename):
    modiRes = modifiedResidues(filename)
    for model in structure:
        for chain in model:
            # Remove HETATM residues
            hetatmResidues = [residue for residue in chain if residue.id[0] != ' ' and residue.resname not in modiRes]
            for residue in hetatmResidues:
                chain.detach_child(residue.id)


def calTorsionAngles(structure):
    torsionAngles = []

    for model in structure:
        for chain in model:
            residues = list(chain.get_residues())
            for i in range(0, len(residues)):
                if i >= 1:
                    alfaAt = [atom.get_vector() for atom in
                              [residues[i - 1]["O3'"], residues[i]["P"], residues[i]["O5'"], residues[i]["C5'"]]]
                    alfa = np.degrees(calc_dihedral(*alfaAt))
                else:
                    alfa = np.nan

                betaAt = [atom.get_vector() for atom in
                          [residues[i]["P"], residues[i]["O5'"], residues[i]["C5'"], residues[i]["C4'"]]]
                beta = np.degrees(calc_dihedral(*betaAt))

                gammaAt = [atom.get_vector() for atom in
                           [residues[i]["O5'"], residues[i]["C5'"], residues[i]["C4'"], residues[i]["C3'"]]]
                gamma = np.degrees(calc_dihedral(*gammaAt))

                deltaAt = [atom.get_vector() for atom in
                           [residues[i]["C5'"], residues[i]["C4'"], residues[i]["C3'"], residues[i]["O3'"]]]
                delta = np.degrees(calc_dihedral(*deltaAt))

                if i < (len(residues) - 1):
                    epsilonAt = [atom.get_vector() for atom in
                                 [residues[i]["C4'"], residues[i]["C3'"], residues[i]["O3'"], residues[i + 1]["P"]]]
                    epsilon = np.degrees(calc_dihedral(*epsilonAt))

                    zetaAt = [atom.get_vector() for atom in
                              [residues[i]["C3'"], residues[i]["O3'"], residues[i + 1]["P"], residues[i + 1]["O5'"]]]
                    zeta = np.degrees(calc_dihedral(*zetaAt))
                else:
                    epsilon, zeta = np.nan, np.nan

                if ("A" or "G") in residues[i].resname:  # pur
                    chiAt = [atom.get_vector() for atom in
                             [residues[i]["O4'"], residues[i]["C1'"], residues[i]["N9"],
                              residues[i]["C4"]]]
                else:
                    chiAt = [atom.get_vector() for atom in
                             [residues[i]["O4'"], residues[i]["C1'"], residues[i]["N1"],
                              residues[i]["C2"]]]
                chi = np.degrees(calc_dihedral(*chiAt))

                torsionAngles.append((alfa, beta, gamma, delta, epsilon, zeta, chi))

    return torsionAngles


def saveTorsionAnglesToCSV(torsion_angles, output_file="torsion_angles.csv"):
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ["Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Chi"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for angles in torsion_angles:
            writer.writerow({
                "Alpha": angles[0],
                "Beta": angles[1],
                "Gamma": angles[2],
                "Delta": angles[3],
                "Epsilon": angles[4],
                "Zeta": angles[5],
                "Chi": angles[6]
            })


def getTorsionAngles(fileName):
    structure = PDB.PDBParser(QUIET=True).get_structure("RNA", fileName)

    removeHetatm(structure, fileName)

    torsion_angles = calTorsionAngles(structure)

    saveTorsionAnglesToCSV(torsion_angles)

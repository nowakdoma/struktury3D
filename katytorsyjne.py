from ramachandranplot import plotRamachandran, calPhiPsi
from Bio import PDB


def main():
    pdb_file_path = "1hhb.pdb"
    structure = PDB.PDBParser(QUIET=True).get_structure("protein", pdb_file_path)
    plotRamachandran(structure)


if __name__ == "__main__":
    main()

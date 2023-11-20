import argparse
from ramachandranplot import plotRamachandran
from Bio import PDB
from rnatorsionmatrix import getTorsionAngles


def parseArgs():
    parser = argparse.ArgumentParser(description="Plot Ramachandran or RNA torsion angles.")
    parser.add_argument("--f", help="Path to the PDB file", required=True)
    parser.add_argument("--p", help="Flag for protein structure", action="store_true")
    parser.add_argument("--r", help="Flag for RNA structure", action="store_true")

    return parser.parse_args()


def main():
    args = parseArgs()
    pdb_file_path = args.f

    if args.p:
        structure = PDB.PDBParser(QUIET=True).get_structure("protein", pdb_file_path)
        plotRamachandran(structure)
    elif args.r:
        getTorsionAngles(pdb_file_path)
    else:
        print("Please specify either --p or --r flag.")


if __name__ == "__main__":
    main()

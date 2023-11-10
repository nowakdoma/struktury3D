import sys
import argparse
from Bio.PDB import *
import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


def parseArgs():
    parser = argparse.ArgumentParser(description="Calculate and visualize protein contact map.")
    parser.add_argument("--f", metavar="PDB_FILE", type=str, help="Input PDB file", required=True)
    parser.add_argument("--c", metavar="CHAIN_ID", type=str, help="Chain ID (default: first chain in PDB file)")
    parser.add_argument("--t", metavar="THRESHOLD", type=float, default=8.0, help="Threshold for contact map ("
                                                                                  "default: 8.0)")

    return parser.parse_args()


def calContactMap(chain, threshold):
    distance = numpy.zeros((len(chain), len(chain)))

    for row, res1 in enumerate(chain):
        for col, res2 in enumerate(chain):
            vectordiff = res1["CA"].coord - res2["CA"].coord
            distance[row, col] = numpy.sqrt(numpy.sum(vectordiff * vectordiff))

    return distance < threshold


def visContactMap(conMap):
    plt.imshow(conMap, cmap=ListedColormap([(1, 0.5, 1), (0.71, 0, 0.71)]), interpolation='none')
    plt.title('Protein Contact Map')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.show()


def main():
    args = parseArgs()

    structure = PDBParser(QUIET=True).get_structure("structure", args.f)

    # check --c
    if args.c is not None:
        try:
            chain = structure[0][args.c]
        except KeyError:
            print(f"Warning: Chain '{args.c}' not found. Choosing the first chain instead.")
            chain = next(structure[0].get_chains(), None)
    else:
        # If --c 0
        chain = next(structure[0].get_chains(), None)

    if chain is None:
        sys.exit("Error: No chains found in the PDB file.")

    # remove hetatm
    for residue in chain:
        if residue.id[0] != ' ':
            chain.detach_child(residue.id)

    conMap = calContactMap(chain, args.t)

    visContactMap(conMap)


if __name__ == "__main__":
    main()

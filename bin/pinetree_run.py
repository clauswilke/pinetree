#! /usr/bin/env python3

import argparse

from pinetree.parser import Parser

def run(params_file, outfile, seed):
    """
    Run simulation.
    """
    with open(params_file, "r") as my_file:
        params_parser = Parser(my_file)

    if seed:
        params_parser.set_seed(seed)

    params_parser.simulation.run(outfile)

    with open(outfile + ".log", "w") as logfile:
        logfile.write("seed, " + str(seed))


def main():
    """
    Parse input file and run simulation.
    """
    parser = argparse.ArgumentParser(description="Simulate transcription \
        and translation.")

    # Add parameter file argument
    parser.add_argument("params", metavar="params", type=str,
                        help="parameter file")
    # Add optional output file
    parser.add_argument("-o", "--outfile", metavar="outfile", type=str, default="results",
                        help="prefix for output files")
    parser.add_argument("-s", "--seed", metavar="seed", type=int)

    args = parser.parse_args()
    my_params_file = args.params
    outfile = args.outfile

    run(my_params_file, outfile, args.seed)



if __name__ == "__main__":
    main()

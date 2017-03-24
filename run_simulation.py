#! /usr/bin/env python3

import argparse

from pysinthe.parser import Parser


def main(my_params_file=""):
    """
    Parse input file and run simulation.
    """
    if my_params_file == "":
        parser = argparse.ArgumentParser(description="Simulate transcription \
            and translation.")

        # Add parameter file argument
        parser.add_argument("params", metavar="params", type=str,
                            help="parameter file")
        # Add optional output file
        parser.add_argument("-o", "--outfile", metavar="outfile", type=str,
                            help="output file")

        args = parser.parse_args()
        my_params_file = args.params
        outfile = args.outfile
    else:
        outfile = False

    with open(my_params_file, "r") as my_file:
        params_parser = Parser(my_file)

    if outfile:
        handle = open(args.outfile, "w")

    for output in params_parser.simulation.run():
        if outfile:
            handle.write(output)
            handle.write("\n")
        else:
            print(output)

    if outfile:
        handle.close()


if __name__ == "__main__":
    main()

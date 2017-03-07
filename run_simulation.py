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
        parser.add_argument("params", metavar="p", type=str, nargs=1,
                            help="parameter file")

        args = parser.parse_args()
        my_params_file = args.params[0]

    with open(my_params_file, "r") as my_file:
        params_parser = Parser(my_file)

    params_parser.simulation.run()


if __name__ == "__main__":
    main()

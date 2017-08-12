#! /usr/bin/env python3

import argparse
import os
import random

from pysinthe_run import run as runner

def main():
    parser = argparse.ArgumentParser(description="Serially run multiple iterations of pysinthe from a single parameter file using different random seeds.")
    parser.add_argument("params", metavar="params", type=str,
                        help="parameter file")
    parser.add_argument("n", metavar="iterations", type=int, help="number of iteractions to run")
    # Add optional output file
    parser.add_argument("-o", "--outdir", metavar="outdir", type=str, default="results",
                        help="output directory")

    args = parser.parse_args()

    os.mkdir(args.outdir)

    for i in range(args.n):
        runner(args.params, args.outdir + "/" + "run" + str(i), int(random.random()*1000000))





if __name__ == "__main__":
    main()
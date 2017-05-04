from pysinthe import core as pysinthe

from pysinthe.parser import Parser

with open("old_tests/params/single_gene.yml", "r") as f:
    parse = Parser(f)

parse.simulation.run()

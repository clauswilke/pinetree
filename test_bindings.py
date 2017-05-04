from pysinthe import core as pysinthe

from pysinthe.parser import Parser

with open("examples/T7_033117.yml", "r") as f:
    parse = Parser(f)

parse.simulation.run()

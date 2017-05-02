import build.pysinthe as pysinthe

sim = pysinthe.Simulation()
test = pysinthe.SpeciesReaction(1.5, ["reactant"], ["product"])
sim.register_reaction(test)

prom = pysinthe.Promoter("test", 1, 10, ["ecolipol"])
print(prom.start, prom.stop)
prom.start = 3
print(prom.start)

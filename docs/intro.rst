Introduction
============

Installation
------------


Construct a simulation
----------------------

All pinetree simulations begin with the construction of a Model object. At a minimum, Simulation must define the volume in which the simulation will take place.

```
import pinetree as pt

model = pt.Model()
``` 

Next, we'll define a genome and register it with Simulation.

Defining a genome
-----------------

Pinetree supports linear genomes of any size, represented by Genome objects. A Genome object must be given a name and a length, in base pairs.

```
plasmid = pt.Genome(name="myplasmid", length=300)
```

After defining a Genome object, we can add promoters, terminators, and genes. These genetic elements can be defined in any order and may even overlap.

```
plasmid.add_promoter()
plasmid.add_terminator()
plasmid.add_gene()
plasmid.add_rnase_site()
```

When all genetic elements have been added, register the Genome object with our Simulation object.

```
model.register_genome(plasmid)
```

At this point we could run the simulation, but nothing would happen because we have not defined any polymerases or ribosomes that interact with the Genome object.

Defining polymerases and ribosomes
----------------------------------

To simulate both transcription and translation, we'll add polymerases and then add ribosomes. Since these enzymes may interact with more than one type of genome, we add them to the Simulation object.

```
model.add_polymerase(name="T7rnapol", speed="", footprint="", copies=100)
model.add_ribosome(speed="", footprint="", copies="")
```

Run the simulation
------------------

model.compile(cell_volume="", degrade_init="", degrade_speed="")
model.simulate(time_limit="", time_step="", output="")

Interpretting results
---------------------


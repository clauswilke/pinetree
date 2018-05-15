Introduction
============

Installation
------------

Pinetree has no requirements except Python and CMake. Python 3 is recommend. To install pinetree from pip, run the following:

.. code-block:: bash
   
   pip install pinetree 

The latest development build may be installed from GitHub:

.. code-block:: bash
   
   git clone https://github.com/benjaminjack/pinetree.git
   pinetree/setup.py install


Construct a simulation
----------------------

All pinetree simulations begin with the construction of a Model object. At a minimum, Model must define the volume in which the simulation will take place.

.. code-block:: python

   import pinetree as pt

   model = pt.Model(cell_volume=8e-16) 

Next, we'll define a genome and register it with Model.

Defining a genome
-----------------

Pinetree supports linear genomes of any size, represented by Genome objects. A Genome object must be given a name and a length, in base pairs.

.. code-block:: python

   plasmid = pt.Genome(name="myplasmid", length=300)

After defining a Genome object, we can add promoters, terminators, and genes. These genetic elements can be defined in any order.

.. code-block:: python

   plasmid.add_promoter(name="phi1", start=1, stop=10, 
                        interactions={"rnapol": 2e8})

   plasmid.add_terminator(name="t1", start=299, stop=300,
                          efficiency={"rnapol": 1.0})

   plasmid.add_gene(name="rnapol", start=26, stop=225,
                    rbs_start=11, rbs_stop=26, rbs_strength=1e7)

   plasmid.add_gene(name="proteinX", start=241, stop=280,
                    rbs_start=226, rbs_stop=241, rbs_strength=1e7)


Here we've defined a plasmid with two genes, 'rnapol', an RNA polymerase which binds to promoter 'phi1' and some other 'proteinX'. Each genetic element has a name, a start position, and a stop position. For more information on the other arguments of these methods, please see the full Python documentation. 

When all genetic elements have been added, register the Genome object with our Model object.

.. code-block:: python

   model.register_genome(plasmid)

At this point we could run the simulation, but nothing would happen because we have not defined any polymerases or ribosomes that interact with the Genome object.

Defining polymerases and ribosomes
----------------------------------

To simulate both transcription and translation, we'll add polymerases and then add ribosomes. Since these enzymes may interact with more than one type of genome, we add them to the Model object.

.. code-block:: python

   model.add_polymerase(name="rnapol", speed=40, footprint=10, copies=10)
   model.add_ribosome(speed=30, footprint=10, copies=100)

Polymerases and ribosomes may move at any speed. Their respective footprints, however, must be smaller than the sites to which they bind. For example, if 'rnapol' has a footprint of 10 bp, then the promoter it binds to must also be at least 10 bp in length. Likewise, if a ribosome has a footprint of 10 bp, the ribosome binding site must be at least 10 bp. 

Define species reactions
------------------------

Pinetree supports option species reactions between one or two molecular species. For example, we may define a reaction such that proteinX forms a complex with rnapol called rnapol-X.

.. code-block:: python
   
   model.add_reaction(reactants=['proteinX', 'rnapol'], 
                      products=['rnapol-X'],
                      rate=1e-7)

Run the simulation
------------------

To simulate gene expression, specify a time limit and a time step at which to output data. All protein and transcript counts will be output in tab seperated format.

.. code-block:: python

   model.simulate(time_limit=60, time_step=1, output="simulation.tsv")


Interpretting results
---------------------

A pinetree simulation produces an output file with 5 columns.


time
    Current time of simulation in seconds.

species
    Name of a molecular species, derived from a polymerase name, a gene name,
    or an explicitly defined molecular species. Any name with a '__' double 
    underscore prefix is used internally by pinetree. 

protein
    Quantity of *free* proteins corresponding to a species name. For example, the number in this column corresponding to 'rnapol' would represent free RNA polymerases that are not actively transcribing.

transcript
    Quantity of transcripts for corresponding to a species name. If a species only exists as a protein or otherwise has no transcript precursor, this value will be 0. 

ribo_density (experimental)
    Average quantity of ribosomes actively translating on a transcript.



Explanation of paramaters
=========================

All input parameters are defined in a YAML file. The parameter is divided into the following sections: ``simulation``, ``species``, ``reactions``, and ``genome``. The only required section is ``simulation``. Each section is described below, followed by a properly formated example. Sections can appear in any order in the parameter file. 

To see a complete parameter file example, see the **github repository**. For a minimal working example, please see the **introduction**.

simulation
----------

Defines basic global parameters for the simulation.

``seed``
    seed for the random number generator
``runtime``
    simulated runtime in seconds
``time_step``
    time step, in seconds, in which to report output
``cell_volume``
    volume of cell, in liters, in which reactions will take place

*Example* ::

    simulation:
        seed: 34
        runtime: 60  # Simulate 60 seconds
        time_step: 1  # Report output every second
        cell_volume: 8e-16

polymerases
-----------

Defines a list of different polymerases in the simulation. Each type of polymerase should define the following fields:

``name``
    name of the polymerase which can be referred to in ``reactions`` and ``elements``
``copy_number``
    the initial number of copies of the polymerase
``speed``
    the speed, in base pairs per second, at which the polymerase transcribes
``footprint`` 
    the footprint, in base pairs, of the polymerase on the genome

*Example* ::

    polymerases:
    - name: rna_pol
      copy_number: 10
      speed: 200
      footprint: 10
    - name: ecoli_pol
      copy_number: 20
      speed: 40
      footprint: 10


ribosomes
---------

Defines the ribosomes in the simulation. Each type of ribosome should define the following fields:

``name``
    name of ribosome which can be referred to in ``reactions``
``copy_number``
    the initial number of copies of the ribosome
``speed``
    the default speed, in base pairs per second, at which the ribosome translates transcripts into proteins; this value may be scaled if codon-specific translation rates are provided
``footprint``
    the footprint, in base pairs, of the ribosome on a transcript
``binding_constant``
    the binding constant, in **units**, of all ribosome-to-ribosome-binding-site interactions

*Example* ::

    ribosomes:
    - name: ribosome
      copy_number: 100
      speed: 30
      footprint: 10
      binding_constant: 1e7

species
-------

reactions
---------

genome
------

Defines basic parameters of the genome.

``name``
    name of genome
``copy_number``
    initial number of copies of the genome (*values other than 1 are not currently supported*)

*Example* ::

    genome:
        name: my_plasmid
        copy_number: 1

elements
--------

Defines all genomic elements in the simulation, including promoters, terminators, and transcripts. `elements` defines a list of these elements, and each element will include a `type`, `name`, `start`, and `stop` fields.

``name``
    name of the genomic element

``type``
    There are three supported types:
    
    - ``promoter``: a transcriptional promoter
    - ``terminator``: a transcriptional terminator
    - ``transcript``: a gene, which will be transcribed and translated into a protein

``start``
    start position of element, in genomic coordinates

``stop``
    stop position of element, in genomic coordinates

.. note::
   All elements can overlap with one another, so define your elements carefully. For example, if you define two overlapping elements in the same reading frame (determined from the start position), the downstream element may not be transcribed completely.


Each of the three element types several additional required fields:

promoter
^^^^^^^^

``interactions``
    defines polymerases and a ``binding_constant`` in which this promoter interacts

*Example* ::
    
    - type: promoter
      name: phi1
      start: 1
      stop: 10
      interactions:
          rnapol:
              binding_constant: 2e8

terminator
^^^^^^^^^^

``interactions``
    defines the polymerases in which this terminator interacts, and the termination ``efficiency``, which ranges from 0.0 (no termination) to 1.0 (all polymerases terminate)

*Example* ::

    - type: terminator
      name: t1
      start: 604
      stop: 605
      interactions:
          rnapol:
              efficiency: 1.0

transcript
^^^^^^^^^^

``rbs``
    location of the ribosome binding site, relative to the start position of the transcript

*Example* ::

    - type: transcript
      name: rnapol
      start: 26
      stop: 225
      rbs: -15



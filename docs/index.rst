Pinetree -- a flexible gene expression simulator with codon-specific translation rates
==============================================================================

.. image:: pinetree-logo.svg

Pinetree is a stochastic gene expression simulator that tracks polymerases
and ribosomes at the single-molecule level. It includes nucleotide-resolution transcription and translation rates. This granularity allows users to simulate the effects of codon-usage bias and dinucleotide bias on gene expression. 

Pinetree is a stochastic simulation that is implemented in Python with a C++ back-end and designed to be highly efficient. On a desktop-class CPU, simulating gene expression in a 40 kilobase viral genome for 30 minutes takes about 3 hours of computation time.

.. toctree:: 
   :caption: User Documentation

   intro
   python_ref
   design
   faq

.. toctree::
   :caption: Developer Documentation
   
   cpp_ref



Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

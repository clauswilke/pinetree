Frequently asked questions
==========================

Do genomic coordinates start from 1 or from 0? Are start and stop coordinates inclusive?
-----------------------------------------------------------------------------

Genomic coordinates *start from 1* and start and stop coordinates are *inclusive*. These conventions follow those established by the GenBank format. Internally, Pinetree converts these coordinates to a 0-based indices, but the end user need not worry about these conversions.

How do I specify transcript degradation rates?
----------------------------------------------

Transcript degradation rates are specified when constructing a new Genome object. Please see the documentation for Genome for further information.

Where did the name come from?
-----------------------------

The name Pinetree is a reference to the "Christmas tree-like" electron micrographs of transcription and translation. Strands of mRNA transcripts branch out from the DNA, and each transcript is studded with ribosomes. This pattern looks like a Christmas tree.
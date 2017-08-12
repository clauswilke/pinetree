Frequently asked questions
==========================

Do genomic coordinates start from 1 or from 0? Are start and stop coordinates inclusive?
-----------------------------------------------------------------------------

Genomic coordinates *start from 1* and start and stop coordinates are *inclusive*. These conventions follow those established by the GenBank format. Internally, Pinetree converts these coordinates to a 0-based indices, but the end user need not worry about these conversions.

How do I specify transcript degradation rates?
----------------------------------------------

Pinetree does not currently support transcript degradation. Any transcript produced during a simulation never degrades.

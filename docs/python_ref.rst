Python reference
===================

.. toctree::
   :maxdepth: 3

This is an reference for the public python pinetree interface.

.. currentmodule:: pinetree

Model
----------
.. autoclass:: Model
    :members:

Genome
------
.. autoclass:: Genome
    :members: add_gene, add_mask, add_promoter, add_terminator, add_weights, add_sequence
    
    .. method:: add_rnase_site(name, start, stop, rate)

        Define an Rnase cleavage site.
        
        :param name: A unique identifier for this RNase cleavage site.
        :type name: string
        :param start: Start position of RNase cleavage site.
        :type start: int
        :param stop: Stop position of RNase cleavage site.
        :type stop: int
        :param rate: Binding rate constant between RNase and this cleavage site.
        :type rate: float

        .. note::

            The internal RNase binding rate constant can alternatively be supplied as an 
            argument during Genome initialization (see ``Genome`` class description above).
            In this case, ``add_rnase_site()`` will only accept values for start and stop,
            and all rnase binding sites will be initialized with the same rate constant. 
            Warning: This method is deprecated and may be removed in the future.

Transcript
----------
.. autoclass:: Transcript
    :members:
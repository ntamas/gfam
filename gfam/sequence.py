"""This module contains a simple `Sequence` and `SeqRecord` (sequence record)
class. `Sequence` instances are just strings for the time being, while a
`SeqRecord` encapsulates a `Sequence` with some associated metadata.

.. _BioPython: http://www.biopython.org
"""

class Sequence(str):
    """String representing a sequence.
    
    This class tries to be compatible with `BioPython`_'s `Sequence`
    class; however, there is no guarantee that it will always stay
    so.

    For the time being, this class is essentially equivalent to a Python
    string.
    """
    pass

class SeqRecord(object):
    """A sequence with its associated metadata.

    This class is very close to `BioPython`_'s `SeqRecord` class;
    however, there is no guarantee about API compatibility.
    
    The class has the following fields:

    - ``seq``: the sequence itself, an instance of `Sequence`.
    - ``id``: a short, unique ID for the sequence
    - ``name``: the human-readable name of the sequence
    - ``description``: the whole description line of the sequence
      as parsed from the original data source (such as a FASTA
      file). The FASTA writer (`gfam.fasta.Writer`) uses this
      when writing the sequence record to a file (unless if there
      is no description, in which case the ID is used).
    """

    __slots__ = ("seq", "id", "name", "description")

    def __init__(self, seq, id="<unknown id>", name="<unknown name>", \
            description="<no description>"):
        self.seq = seq
        self.id = id
        self.name = name
        self.description = description

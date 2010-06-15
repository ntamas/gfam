"""Sequences, domain assignments and related routines."""

class Sequence(str):
    """String representing a sequence.
    
    This class tries to be compatible with BioPython's ``Sequence``
    class; however, there is no guarantee that it will always stay
    so.
    """
    pass

class SeqRecord(object):
    """A sequence with its associated metadata.

    This class is very close to BioPython's ``SeqRecord`` class;
    however, there is no guarantee about API compatibility. API
    compatibility is ensured only to an extent which is necessary
    to keep ``gfam`` working with or without ``BioPython``.
    """

    __slots__ = ("seq", "id", "name", "description")

    def __init__(self, seq, id="<unknown id>", name="<unknown name>", \
            description="<no description>"):
        self.seq = seq
        self.id = id
        self.name = name
        self.description = description

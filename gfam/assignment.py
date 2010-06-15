"""Classes corresponding to domain assignments (`Assignment`) and
sequences with corresponding domain assignments (`SequenceWithAssignment`).
"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["Assignment", "bidict", "EValueFilter", "open_anything", \
           "redirected", "SequenceWithAssignments", "search_file", \
           "temporary_dir", "UniqueIdGenerator", "UniversalSet"]

try:
    from collections import namedtuple
except ImportError:
    # For Python 2.5 and older
    from gfam.compat import namedtuple

from gfam.enum import Enum

# pylint: disable-msg=C0103,E1101
# C0103: invalid name
# E1101: instance has no 'foo' member. Pylint does not know namedtuple
# trickery.
class Assignment(namedtuple("Assignment", \
    "id length start end source domain evalue interpro_id comment")):
    """Class representing a record in an InterPro ``iprscan`` output.
    
    An InterPro domain assignment has the following fields:
        
        - ``id``: the ID of the sequence

        - ``length``: the length of the sequence

        - ``start``: the starting position of the region in the sequence
          that is assigned to some InterPro domain (inclusive)

        - ``end``: the ending position (inclusive)

        - ``source``: the assignment source as reported by ``iprscan``

        - ``domain``: the ID of the domain being assigned, according to the
          assignment source

        - ``evalue``: E-value of the assignment if that makes sense,
          ``None`` otherwise

        - ``interpro_id``: the InterPro ID corresponding to ``domain`` in
          ``source``.

        - ``comment``: an arbitrary comment
        """

    __slots__ = ()

    def resolve_interpro_ids(self, interpro):
        """If the assignment has an InterPro ID, this method makes sure
        that the domain is equal to the highest common ancestor of the
        InterPro ID in the InterPro tree. If the assignment does not have
        an InterPro ID yet, this method tries to look it up.
        
        Returns a new tuple which might or might not be equal to this one.
        """
        if self.interpro_id:
            anc = interpro.tree.get_most_remote_ancestor(self.interpro_id)
        else:
            anc = interpro.mapping.get(self.domain)
        if self.domain == anc:
            return self
        return self._replace(domain=anc)

    def short_repr(self):
        """Short representation of this assignment, used in error messages"""
        return "%s(%d-%d)" % (self.domain, self.start, self.end)


class OverlapType(Enum):
    """Enum describing the different overlap types that can be detected
    by `AssignmentOverlapChecker`. See the documentation of
    `AssignmentOverlapChecker._check_single` for more details.
    """
    NO_OVERLAP = "NO_OVERLAP"
    DUPLICATE = "DUPLICATE"
    INSERTION_DIFFERENT = "INSERTION_DIFFERENT"
    DIFFERENT = "DIFFERENT"
    OVERLAP = "OVERLAP"


class AssignmentOverlapChecker(object):
    """Static class that contains the central logic of determining
    whether an assignment can be added to a partially assigned
    `SequenceWithAssignments`_.
    """

    max_overlap = 20

    @classmethod
    def check(cls, sequence, assignment):
        """Checks whether an `assignment` can be added to a partially
        assigned `sequence`. `sequence` must be an instance of
        `SequenceWithAssignments`, `assignment` must be an instance
        of `Assignment`.

        The output is equivalent to the output of the first `_check_single`
        that returns anything different from `OverlapType.NO_OVERLAP`,
        or `OverlapType.NO_OVERLAP` otherwise.
        """
        for other_assignment in sequence.assignments:
            result = cls._check_single(assignment, other_assignment)
            if result != OverlapType.NO_OVERLAP:
                return result
        return OverlapType.NO_OVERLAP

    @classmethod
    def _check_single(cls, assignment, other_assignment):
        """Checks whether the given `assignment` overlaps with another
        assignment `other_assignment`. Returns one of the following:

            - `OverlapType.NO_OVERLAP`: there is no overlap between the
              two given assignments

            - `OverlapType.DUPLICATE`: `assignment` is a duplicate of
              `other_assignment` (same starting and ending positions)

            - `OverlapType.INSERTION_DIFFERENT`: `assignment` is inserted into
              `other_assignment` or vice versa, but they have different data
              sources.

            - `OverlapType.DIFFERENT`: `other_assignment` overlaps with `assignment`
              partially, but they have different data sources

            - `OverlapType.OVERLAP`: `other_assignment` overlaps with `assignment`
              partially, they have the same data source, but the size of
              the overlap is larger than the maximum allowed overlap specified
              in `OverlapChecker.max_overlap`.
        """
        start, end = assignment.start, assignment.end
        other_start, other_end = other_assignment.start, other_assignment.end

        if other_start == start and other_end == end:
            # This is a duplicate assignment, so we must skip it
            return OverlapType.DUPLICATE

        if other_start <= start and other_end >= end:
            if other_assignment.source == assignment.source:
                # This is a valid domain insertion, the new domain is the
                # one which was inserted into the old one
                return OverlapType.NO_OVERLAP
            return OverlapType.INSERTION_DIFFERENT

        if other_start >= start and other_end <= end:
            if other_assignment.source == assignment.source:
                # This is a valid domain insertion, the new domain is the
                # one which contains the old one completely
                return OverlapType.NO_OVERLAP
            return OverlapType.INSERTION_DIFFERENT

        if other_start <= start and other_end <= end and other_end >= start:
            if other_assignment.source == assignment.source:
                # This is a partial overlap
                overlap_size = other_end-start+1
                if overlap_size > cls.max_overlap:
                    return OverlapType.OVERLAP
            else:
                return OverlapType.DIFFERENT

        if other_start >= start and other_end >= end and other_start <= end:
            if other_assignment.source == assignment.source:
                # This is a partial overlap
                overlap_size = end-other_start+1
                if overlap_size > cls.max_overlap:
                    return OverlapType.OVERLAP
            else:
                return OverlapType.DIFFERENT

        return OverlapType.NO_OVERLAP


class SequenceWithAssignments(object):
    """Class representing a sequence for which some parts are assigned to
    InterPro domains.

    The class has the following fields:

        - ``name``: the name of the sequence

        - ``length``: the number of amino acids in the sequence

        - ``assignments``: a list of `Assignment`_ instances that describe
          the domain architecture of the sequence
    """

    __slots__ = ("name", "length", "assignments")
    overlap_checker = AssignmentOverlapChecker

    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.assignments = []

    def __len__(self):
        return self.length

    def assign_(self, start, end, domain, source="Novel", *args, **kwds):
        """Assigns a fragment of this sequence to the given domain.
        `start` and `end` are the starting and ending positions, inclusive.
        `domain_name` is the name of the domain, `source` is the assignment
        source (``Novel`` by default).
        """
        assignment = Assignment(id=self.name, start=start, end=end, \
                interpro_id=None, source=source, domain=domain, \
                evalue=None, length=self.length, comment=None)
        return self.assign(assignment, *args, **kwds)

    def assign(self, assignment, overlap_check=True):
        """Assigns a fragment of this sequence using the given assignment.
        If `overlap_check` is ``False``, we will not check for overlaps or
        conflicts with existing assignments.

        Returns ``True`` if the assignment was added, ``False`` if it
        wasn't due to an overlap conflict.
        """
        if ":SF" in assignment.domain:
            # Remove subfamily info from domain
            # pylint: disable-msg=W0212
            new_domain = assignment.domain[0:assignment.domain.index(":SF")]
            assignment = assignment._replace(domain=new_domain)

        if overlap_check:
            overlap_state = self.overlap_checker.check(self, assignment)
            if overlap_state != OverlapType.NO_OVERLAP:
                return False

        self.assignments.append(assignment)
        return True

    def coverage(self):
        """Returns the coverage of the sequence, i.e. the fraction of residues
        covered by at least one assignment."""
        ok = [0] * self.length
        for a in self.assignments:
            ok[a.start:(a.end+1)] = [1] * ((a.end+1)-a.start)
        return sum(ok) / float(self.length)

    def is_completely_unassigned(self, start, end):
        """Checks whether the given region is completely unassigned.
        start and end positions are both inclusive"""
        return all(a.end < start or a.start > end for a in self.assignments)

    def resolve_interpro_ids(self, interpro):
        """Calls resolve_interpro_ids() on each assignment of this sequence"""
        self.assignments = [assignment.resolve_interpro_ids(interpro) \
                            for assignment in self.assignments]

    def unassigned_regions(self):
        """Returns a generator that iterates over the unassigned regions
        of the sequence. Each entry yielded by the generator is a tuple
        containing the start and end positions"""
        ok = [True] * (self.length+1)
        for a in self.assignments:
            ok[a.start:(a.end+1)] = [False] * ((a.end+1)-a.start)
        i = 1
        while i <= self.length:
            while i <= self.length and not ok[i]:
                i += 1
            start = i
            if start == self.length:
                break
            while i <= self.length and ok[i]:
                i += 1
            yield start, i - 1




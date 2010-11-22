"""Classes corresponding to domain assignments (`Assignment`) and
sequences with corresponding domain assignments (`SequenceWithAssignments`).
"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["Assignment", "AssignmentOverlapChecker", "OverlapType",
           "SequenceWithAssignments", "EValueFilter"]

try:
    from collections import namedtuple
except ImportError:
    # For Python 2.5 and older
    from gfam.compat import namedtuple

from gfam.enum import Enum

import operator

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

    def get_assigned_length(self):
        """Returns the number of amino acids covered by the assignment
        within the sequence."""
        return self.end - self.start + 1

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
    `AssignmentOverlapChecker.check_single` for more details.
    """
    NO_OVERLAP = "NO_OVERLAP"
    DUPLICATE = "DUPLICATE"
    INSERTION = "INSERTION"
    INSERTION_DIFFERENT = "INSERTION_DIFFERENT"
    DIFFERENT = "DIFFERENT"
    OVERLAP = "OVERLAP"


class AssignmentOverlapChecker(object):
    """Static class that contains the central logic of determining
    whether an assignment can be added to a partially assigned
    `SequenceWithAssignments`.

    The class has a class variable named `max_overlap` which stores
    the maximum allowed overlap size. This is 20 by default.
    """

    #: The maximum allowed overlap size.
    max_overlap = 20

    @classmethod
    def check(cls, sequence, assignment):
        """Checks whether an `assignment` can be added to a partially
        assigned `sequence`. `sequence` must be an instance of
        `SequenceWithAssignments`, `assignment` must be an instance
        of `Assignment`.

        The output is equivalent to the output of the first `check_single`
        that returns anything different from `OverlapType.NO_OVERLAP`,
        or `OverlapType.NO_OVERLAP` otherwise.
        """
        for other_assignment in sequence.assignments:
            result = cls.check_single(assignment, other_assignment)
            if result != OverlapType.NO_OVERLAP:
                return result
        return OverlapType.NO_OVERLAP

    @classmethod
    def check_single(cls, assignment, other_assignment):
        """Checks whether the given `assignment` overlaps with another
        assignment `other_assignment`. Returns one of the following:

        - `OverlapType.NO_OVERLAP`: there is no overlap between the
          two given assignments

        - `OverlapType.DUPLICATE`: `assignment` is a duplicate of
          `other_assignment` (same starting and ending positions)

        - `OverlapType.INSERTION`: there is a complete domain insertion in
          either direction

        - `OverlapType.INSERTION_DIFFERENT`: `assignment` is inserted into
          `other_assignment` or vice versa, but they have different data
          sources.

        - `OverlapType.DIFFERENT`: `other_assignment` overlaps with `assignment`
          partially, but they have different data sources

        - `OverlapType.OVERLAP`: `other_assignment` overlaps with `assignment`
          partially, they have the same data source, but the size of
          the overlap is larger than the maximum allowed overlap specified
          in `AssignmentOverlapChecker.max_overlap`.
        """
        start, end = assignment.start, assignment.end
        other_start, other_end = other_assignment.start, other_assignment.end

        if other_start == start and other_end == end:
            # This is a duplicate assignment, so we must skip it
            return OverlapType.DUPLICATE

        if other_start <= start and other_end >= end:
            if other_assignment.source == assignment.source:
                # This is a valid domain insertion, assignment is inserted
                # into other_assignment
                return OverlapType.INSERTION
            return OverlapType.INSERTION_DIFFERENT

        if other_start >= start and other_end <= end:
            if other_assignment.source == assignment.source:
                # This is a valid domain insertion, other_assignment is
                # inserted into assignment
                return OverlapType.INSERTION
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

    @classmethod
    def get_overlap_size(cls, assignment, other_assignment):
        """Returns the length of the overlap between the given `assignment`
        and another (`other_assignment`). It is assumed (and not checked)
        that the two assignments refer to the same sequence.
        """
        start, end = assignment.start, assignment.end
        other_start, other_end = other_assignment.start, other_assignment.end

        if other_start <= start and other_end >= end:
            return end-start+1

        if other_start >= start and other_end <= end:
            return other_end-other_start+1

        if other_start <= start and other_end <= end and other_end >= start:
            return other_end-start+1

        if other_start >= start and other_end >= end and other_start <= end:
            return end-other_start+1

        return 0


class SequenceWithAssignments(object):
    """Class representing a sequence for which some parts are assigned to
    InterPro domains.

    The class has the following fields:

    - ``name``: the name of the sequence

    - ``length``: the number of amino acids in the sequence

    - ``assignments``: a list of `Assignment` instances that describe
      the domain architecture of the sequence
    """

    __slots__ = ("name", "length", "assignments")

    #: The overlap checker used by this instance. This points to
    #: `AssignmentOverlapChecker` by default.
    overlap_checker = AssignmentOverlapChecker

    #: Acceptable overlap types which we allow in assignments.
    #: By default, we consider complete domain insertions and the
    #: absence of overlaps acceptable.
    acceptable_overlaps = set([OverlapType.NO_OVERLAP, OverlapType.INSERTION])

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
            if overlap_state not in self.acceptable_overlaps:
                return False

        self.assignments.append(assignment)
        return True

    def coverage(self, sources=None):
        """Returns the coverage of the sequence, i.e. the fraction of residues
        covered by at least one assignment.
        
        `sources` specifies the data sources to be included in the coverage
        calculation. If `None`, all the data sources will be considered; otherwise
        it must be a set containing the accepted sources."""
        ok = [0] * self.length
        if sources is None:
            for a in self.assignments:
                ok[a.start:(a.end+1)] = [1] * ((a.end+1)-a.start)
        else:
            if isinstance(sources, basestring):
                sources = [sources]
            for a in self.assignments:
                if a.source in sources:
                    ok[a.start:(a.end+1)] = [1] * ((a.end+1)-a.start)
        return sum(ok) / float(self.length)

    def data_sources(self):
        """Returns the list of data sources that were used in this assignment."""
        return sorted(set(a.source for a in self.assignments))

    def domain_architecture(self, sources=None):
        """Returns the domain architecture of the assignment.

        The domain architecture is a list which contains the IDs of the assigned
        regions (domains) in ascending order of their starting positions. If
        `sources` is ``None``, all data sources will be considered; otherwise it
        must be a set or iterable which specifies the data sources to be
        included in the result.
        """
        sorted_assignments = sorted(self.assignments, key=operator.attrgetter("start"))
        if sources is None:
            return [a.domain for a in sorted_assignments]
        if isinstance(sources, basestring):
            sources = [sources]
        return [a.domain for a in sorted_assignments if a.source in sources]

    def is_completely_unassigned(self, start, end):
        """Checks whether the given region is completely unassigned.
        start and end positions are both inclusive"""
        return all(a.end < start or a.start > end for a in self.assignments)

    def resolve_interpro_ids(self, interpro):
        """Calls `Assignment.resolve_interpro_ids` on each assignment of this
        sequence"""
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


class EValueFilter(object):
    """Given an `Assignment`, this filter tells whether the assignment's
    E-value is satisfactory to accept it.

    The filter supports different E-values for different data sources. By
    default, the E-value threshold is infinity for all data sources.
    """

    def __init__(self):
        self.default_e_value = float('inf')
        self.thresholds = {}

    def set_threshold(self, source, evalue):
        """Sets the E-value threshold for the given data source."""
        self.thresholds[source] = evalue

    def is_acceptable(self, assignment):
        """Checks whether the given assignment is acceptable.
        
        This method looks up the E-value threshold corresponding to the
        ``source`` of the assignment and returns ``True`` if the E-value
        of the assignment is less than the threshold, ``False``
        otherwise."""
        threshold = self.thresholds.get(assignment.source, self.default_e_value)
        return assignment.evalue <= threshold

    @classmethod
    def FromString(cls, description):
        """Constructs an E-value filter from a string description that
        can be used in command line arguments and configuration files.
        
        The string description is a semicolon-separated list of
        source-threshold pairs. For instance, the following is a valid
        description giving an E-value of 0.001 for HMMPfam sources,
        0.005 for HMMSmart sources and 0.007 for everything else::

            HMMPfam=0.001;HMMSmart=0.005;0.007
        
        The last entry denotes the default E-value; in particular, if a
        number is not preceded by a source name, it is assumed to be a
        default E-value. If no default E-value is given, infinity will
        be used.
        """
        result = cls()
        for part in description.split(";"):
            part = part.strip()
            if "=" in part:
                source, evalue = part.split("=", 1)
                evalue = float(evalue)
                result.set_threshold(source, evalue)
            else:
                result.default_e_value = float(part)
        return result


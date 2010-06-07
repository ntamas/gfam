"""Common routines and utility classes for GFam"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["Assignment", "EValueFilter", "open_anything", \
           "Sequence", "UniqueIdGenerator", "UniversalSet"]

import bz2
import gzip
import sys
import urllib2

from collections import namedtuple

# pylint: disable-msg=C0103
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


class EValueFilter(object):
    """Given an assignment, this filter tells whether the assignment's
    E-value is satisfactory to accept it.

    Different E-values may be given for different data sources.
    """

    def __init__(self):
        self.default_e_value = float('inf')
        self.thresholds = {}

    def set_threshold(self, source, evalue):
        """Sets the E-value threshold for the given data source"""
        self.thresholds[source] = evalue

    def is_acceptable(self, assignment):
        """Checks whether the given assignment is acceptable"""
        threshold = self.thresholds.get(assignment.source, self.default_e_value)
        return assignment.evalue <= threshold

    @classmethod
    def FromString(cls, description):
        """Constructs an E-value filter from a string description that
        can be used in command line arguments.
        
        The string description is a semicolon-separated list of source-threshold
        pairs; for instance: HMMPfam=0.001;HMMSmart=0.005;0.007. The last
        entry denotes the default E-value; in particular, if a number is
        not preceded by a source name, it is assumed to be a default E-value.
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


class Sequence(object):
    """Class representing a sequence for which some parts are assigned to
    InterPro domains.

    The class has the following fields:

        - ``name``: the name of the sequence

        - ``length``: the number of amino acids in the sequence

        - ``assignments``: a list of `Assignment`_ instances that describe
          the domain architecture of the sequence
    """

    __slots__ = ("name", "length", "assignments")

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

    def assign(self, assignment, silent=False, overlap_check=True):
        """Assigns a fragment of this sequence using the given assignment.
        If `silent` is ``True``, no error messages are printed about overlapping
        assignments or conflicts with existing assignments. If `overlap_check`
        is ``False``, we will not check for overlaps or conflicts with existing
        assignments.
        """
        if ":SF" in assignment.domain:
            # Remove subfamily info from domain
            # pylint: disable-msg=W0212
            new_domain = assignment.domain[0:assignment.domain.index(":SF")]
            assignment = assignment._replace(domain=new_domain)

        skip, overlap_with = None, None
        if overlap_check:
            for a in self.assignments:
                skip = self._overlap_check(assignment, a)
                if skip:
                    overlap_with = a
                    break

        if not skip:
            self.assignments.append(assignment)
            return True

        if skip == "duplicate":
            return False

        if silent:
            return False

        if skip == "overlap":
            print >>sys.stderr, "WARNING: partially overlapping " +\
                     "assignment for %s: %s and %s" % (self.name, \
                     assignment.short_repr(), \
                     overlap_with.short_repr())
        elif skip == "overlap_different":
            print >>sys.stderr, "WARNING: partially overlapping " +\
                     "assignment for %s from different sources: "+\
                     "%s and %s" % (self.name, \
                     assignment.short_repr(), \
                     overlap_with.short_repr())
        elif skip == "insertion_different":
            print >>sys.stderr, "WARNING: fully overlapping " +\
                     "assignment for %s from different sources: "+\
                     "%s and %s" % (self.name, \
                     assignment.short_repr(), \
                     overlap_with.short_repr())
        else:
            print >>sys.stderr, "WARNING: skipping assignment %s, reason: %s"\
                    % (assignment.short_repr(), skip)
        return False

    def _overlap_check(self, assignment, a):
        """Checks whether the given `assignment` overlaps with another assignment
        `a`. Returns ``None`` if there is no overlap, or the type of overlap
        as a string. The type of overlap is one of the following:

            - ``duplicate``: `assignment` is a duplicate of `a` (same starting and
              ending positions)

            - ``insertion_different``: `a` is inserted into `assignment` or vice
              versa, but they have different data sources.

            - ``overlap_different``: `a` overlaps with `assignment` partially,
              but they have different data sources.

            - ``overlap``: `a` overlaps with `assignment` partially, and although
              they have the same data source, the overlap is too large, so `a`
              cannot be allowed.
        """
        domain, start, end = assignment.domain, assignment.start, assignment.end

        if a.start == start and a.end == end:
            # This is a duplicate assignment, so we must skip it
            return "duplicate"

        if a.start <= start and a.end >= end:
            if a.source == assignment.source:
                # This is a valid domain insertion, the new domain is the
                # one which was inserted into the old one
                return None
            return "insertion_different"

        if a.start >= start and a.end <= end:
            if a.source == assignment.source:
                # This is a valid domain insertion, the new domain is the
                # one which contains the old one completely
                return None
            return "insertion_different"

        if a.start <= start and a.end <= end and a.end >= start:
            if a.source == assignment.source:
                # This is a partial overlap
                overlap_size = a.end-start+1
                if overlap_size > 20:
                    return "overlap"
            else:
                return "overlap_different"

        if a.start >= start and a.end >= end and a.start <= end:
            if a.source == assignment.source:
                # This is a partial overlap
                overlap_size = end-a.start+1
                if overlap_size > 20:
                    return "overlap"
            else:
                return "overlap_different"

        return None

    def coverage(self):
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


def open_anything(fname, *args, **kwds):
    """Opens the given file. The file may be given as a file object
    or a filename. If the filename ends in .bz2 or .gz, it will
    automatically be decompressed on the fly. If the filename starts
    with ``http://``, ``https://`` or ``ftp://`` and there is no
    other argument given, the remote URL will be opened for reading.
    A single dash in place of the filename means the standard input.
    """
    if isinstance(fname, file):
        infile = fname
    elif fname == "-":
        infile = sys.stdin
    elif (fname.startswith("http://") or fname.startswith("ftp://") or \
         fname.startswith("https://")) and not kwds and not args:
        infile = urllib2.urlopen(fname)
    elif fname[-4:] == ".bz2":
        infile = bz2.BZ2File(fname, *args, **kwds)
    elif fname[-3:] == ".gz":
        infile = gzip.GzipFile(fname, *args, **kwds)
    else:
        infile = open(fname, *args, **kwds)
    return infile


class UniqueIdGenerator(object):
    """A dictionary-like class that can be used to assign unique integer IDs to
    names.

    Usage:
    
    >>> gen = UniqueIdGenerator()
    >>> gen["A"]
    0
    >>> gen["B"]
    1
    >>> gen["C"]
    2
    >>> gen["A"]      # Retrieving already existing ID
    0
    >>> len(gen)      # Number of already used IDs
    3
    """

    def __init__(self, id_generator=None):
        """Creates a new unique ID generator. `id_generator` specifies how do we
        assign new IDs to elements that do not have an ID yet. If it is `None`,
        elements will be assigned integer identifiers starting from 0. If it is
        an integer, elements will be assigned identifiers starting from the given
        integer. If it is an iterator or generator, its `next` method will be
        called every time a new ID is needed."""
        if id_generator is None: id_generator = 0
        if isinstance(id_generator, int):
            import itertools
            self._generator = itertools.count(id_generator)
        else:
            self._generator = id_generator
        self._ids = {}

    def __getitem__(self, item):
        """Retrieves the ID corresponding to `item`. Generates a new ID for `item`
        if it is the first time we request an ID for it."""
        try:
            return self._ids[item]
        except KeyError:
            self._ids[item] = self._generator.next()
            return self._ids[item]

    def __len__(self):
        """Retrieves the number of added elements in this UniqueIDGenerator"""
        return len(self._ids)

    def reverse_dict(self):
        """Returns the reversed mapping, i.e., the one that maps generated IDs to their
        corresponding items"""
        return dict((v, k) for k, v in self._ids.iteritems())

    def values(self):
        """Returns the list of items added so far. Items are ordered according to
        the standard sorting order of their keys, so the values will be exactly
        in the same order they were added if the ID generator generates IDs in
        ascending order. This hold, for instance, to numeric ID generators that
        assign integers starting from a given number."""
        return sorted(self._ids.keys(), key = self._ids.__getitem__)


class UniversalSet(object):
    """Does what it says on the tin: it is a fake set containing everything.
    
    This class can be used in places where a set is expected, and the only
    operation performed on the set will be membership checking. The class
    returns ``True`` when the user searches for a given member in the set,
    regardless of what that member is.
    """

    def __contains__(self, what):
        return True




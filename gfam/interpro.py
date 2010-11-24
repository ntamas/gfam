"""Classes related to handling InterPro-related files in HyFam"""

import re

from collections import defaultdict
from gfam.assignment import Assignment
from gfam.utils import bidict, open_anything

try:
    from collections import Mapping
except ImportError:
    from gfam.compat import Mapping

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["AssignmentReader", "InterPro", "InterProIDMapper", \
           "InterProNames", "InterProTree", "InterPro2GOMapping"]


class AssignmentReader(object):
    """Iterates over assignments in an InterPro domain assignment file.
    
    This reader parses the output of ``iprscan`` and yields appropriate
    `Assignment` instances for each line.
    """

    def __init__(self, filename):
        self._fp = open_anything(filename)

    def assignments(self):
        """A generator that yields the assignments in the InterPro domain
        assignment file one by one. Each object yielded by this generator
        will be an instance of `Assignment`."""
        for line in self._fp:
            assignment = self.parse_line(line)
            if assignment is not None:
                yield assignment

    def assignments_and_lines(self):
        """A generator that yields the assignments in the InterPro domain
        assignment file and the corresponding raw lines one by one. Each
        object yielded by this generator will be a tuple containing an
        instance of `Assignment` and the corresponding line."""
        for line in self._fp:
            assignment = self.parse_line(line)
            if assignment is not None:
                yield assignment, line

    def parse_line(self, line):
        """Parses a single line from an InterPro domain assignment file and
        returns a corresponding `Assignment` instance."""
        line = line.strip()
        if not line:
            return None

        parts = line.split("\t")

        if len(parts) < 7:
            raise ValueError(repr(line))

        if len(parts) < 15:
            parts.extend([None] * (15-len(parts)))

        try:
            evalue = float(parts[8])
        except (ValueError, TypeError):
            evalue = None

        if parts[11] == 'NULL' or not parts[11]:
            parts[11] = None

        assignment = Assignment( \
                id = parts[0],
                length = int(parts[2]),
                source = parts[3],
                domain = parts[4],
                start = int(parts[6]),
                end = int(parts[7]),
                evalue = evalue,
                interpro_id = parts[11],
                comment = parts[14]
        )

        return assignment


    def __iter__(self):
        return self.assignments()


class InterProTree(Mapping):
    """Dict-like object that tells the parent ID corresponding to every InterPro
    domain ID.
    
    This class is usually not constructed directly; an instance of this class
    is a member of every instance of `InterPro`.

    The class behaves much like a dictionary that returns the parent ID for
    every InterPro domain ID. For unknown domain IDs, the domain ID itself is
    returned. For sub-subfamilies, the returned value is the ID of the subfamily,
    not the family; use `get_most_remote_ancestor` if you always need the family
    ID no matter what.

    Usage example::

        >>> interpro = InterPro.FromFile("data/ParentChildTreeFile.txt")
        >>> tree = interpro.tree
        >>> tree["IPR000010"]      # this is a family ID
        'IPR000010'
        >>> tree["IPR001713"]      # this is a subfamily of IPR000010
        'IPR000010'
        >>> tree["IPR001321"]      # a sub-subfamily ID
        'IPR013655'
        >>> tree["IPR9999"]        # this is an unknown ID
        'IPR9999'
        >>> "IPR9999" in tree
        False
        >>> tree.get_most_remote_ancestor("IPR001321")
        'IPR000014'
    """

    def __init__(self):
        self._data = {}

    def __contains__(self, item):
        return self._data.__contains__(item)

    def __getitem__(self, item):
        return self._data.get(item, item)

    def get_most_remote_ancestor(self, item):
        """Returns the most remote ancestor of the given item.
        
        For family IDs, this returns the ID itself. For subfamily IDs and below,
        this returns the corresponding family ID."""
        parent, child = self._data.get(item, item), item
        while parent != child:
            parent, child = self._data.get(parent, parent), parent
        return parent

    def __iter__(self):
        return self._data.__iter__()

    def __len__(self):
        return len(self._data)

    def __setitem__(self, child, parent):
        self._data[child] = parent

    def __delitem__(self, child):
        del self._data[child]


class InterProIDMapper(object):
    """Dict-like object that maps domain IDs from various data sources
    to their corresponding InterPro IDs.
    
    This class is usually not constructed directly; an instance of this class
    is a member of every instance of `InterPro`.

    The class can generally be used like a dictionary::

        >>> interpro = InterPro.FromFile("data/ParentChildTreeFile.txt")
        >>> mapper = interpro.mapping
        >>> mapper["PF00031"]          # an alias for IPR000010
        'IPR000010'
        >>> mapper["IPR9999"]          # no such ID
        Traceback (most recent call last):
            ...
        KeyError: 'IPR9999'
        >>> mapper.get("IPR9999")      # no such ID
        'IPR9999'
        >>> mapper["IPR9999"] = "Fake ID for testing"
        >>> mapper["IPR9999"]
        'Fake ID for testing'
        >>> del mapper["IPR9999"]
        >>> "IPR9999" in mapper
        False
    """

    def __init__(self):
        self._data = {}

    def __getitem__(self, key):
        return self._data.__getitem__(key)

    def __setitem__(self, key, value):
        self._data.__setitem__(key, value)

    def __contains__(self, key):
        return self._data.__contains__(key)

    def __delitem__(self, key):
        self._data.__delitem__(key)

    def __len__(self):
        return len(self._data)

    def get(self, key, default=None):
        """Returns the InterPro ID corresponding to the given key or the
        given default value if no InterPro ID exists for that key. If
        the default value is None, returns the key itself if there is
        no InterPro ID for the key"""
        return self._data.get(key, default or key)


class InterProNames(object):
    """Dict-like object mapping IDs to human-readable names, the only
    difference being that unknown IDs are handled gracefully instead of
    raising a `KeyError`.

    The name of this class is a bit of a misnomer as it works for *any*
    type of IDs, not only for InterPro IDs.

    .. todo:: refactor and fix the class name

    Usage example::

        >>> names = InterProNames.FromFile("data/names.dat.gz")
        >>> "IPR015503" in names
        True
        >>> "no-such-name" in names
        False
        >>> names["IPR015503"]
        'Cortactin'
        >>> names["no-such-name"]
        'no-such-name'
    """

    def __init__(self):
        self.names = {}

    def load(self, filename):
        """Loads ID-name assignments from a simple tab-separated flat file.
        
        Lines not containins any tab characters are silently ignored."""
        for line in open_anything(filename):
            parts = line.strip().split("\t", 1)
            if len(parts) > 1:
                self.names[parts[0]] = parts[1]

    def __contains__(self, name):
        return name in self.names

    def __getitem__(self, name):
        return self.names.get(name, name)

    @classmethod
    def FromFile(cls, filename):
        """Shortcut method that does exactly what the following snippet does::
        
            >>> names = InterProNames()
            >>> names.load(filename)           #doctest: +SKIP
        """
        result = cls()
        result.load(filename)
        return result


class InterPro(object):
    """Class that encapsulates the InterPro parent-child tree (an instance
    of `InterProTree`) and the InterPro ID mapper (an instance of
    `InterProIDMapper`) under the same hood. This makes it easier to pass
    both of them around in the code.

    When parsing the parent-child tree, subfamiy IDs are removed from
    aliases starting with ``PTHR``; i.e. ``PTHR10829:SF4`` will be
    stored as ``PTHR10829`` (and mapped to ``IPR015503`` at the time of
    writing).

    For usage examples, see `InterProTree` and `InterProIDMapper`.
    """

    def __init__(self):
        self.tree = InterProTree()
        self.mapping = InterProIDMapper()

    @classmethod
    def FromFile(cls, filename):
        """Constructs this object from an InterPro parent-child mapping file,
        pointed to by the given filename. Both the tree and the ID-name mapping
        will be built from the same file.
        """
        result = cls()
        path_to_root = []

        for line in open_anything(filename):
            line = line.strip()
            dash_count = 0
            while line[dash_count] == "-":
                dash_count += 1
            if dash_count % 2 != 0:
                raise ValueError("dash count in InterPro file not even")

            line = line[dash_count:]
            parts = line.split("::")
            interpro_id, aliases = parts[0], parts[2:]

            level = dash_count // 2 + 1
            if level <= len(path_to_root):
                path_to_root = path_to_root[:level]
                path_to_root[-1] = interpro_id
            else:
                path_to_root.append(interpro_id)
                if level != len(path_to_root):
                    raise ValueError("tree depth increased by more than "
                                     "one between two lines")
            if len(path_to_root) > 1:
                result.tree[interpro_id] = path_to_root[-2]

            for alias in aliases:
                if alias[0:4] == "PTHR" and ":SF" in alias:
                    alias = alias[0:alias.index(":SF")]
                result.mapping[alias] = interpro_id

        return result


class InterPro2GOMapping(bidict):
    """Bidirectional dictionary (`bidict`) that tells the corresponding Gene
    Ontology terms of every InterPro ID and vice versa.

    This object encapsulates two dictionaries: `self.terms` gives the Gene
    Ontology terms of a given InterPro domain, and `self.domains` gives the
    InterPro domains annotated by a given GO term."""

    def __init__(self):
        bidict.__init__(self)
        self.terms = self.left
        self.domains = self.right

    def add_annotation(self, interpro_id, go_id):
        """Adds a single Gene Ontology ID to the list of Gene
        Ontology IDs for a given InterPro ID"""
        self.add_left(interpro_id, go_id)

    @classmethod
    def from_file(cls, filename, tree):
        """Constructs a mapping from a mapping file. The format of this
        file should be identical to the official ``interpro2go`` file
        provided by the Gene Ontology project. `tree` is a Gene Ontology
        tree object (see `gfam.go.Tree`) that will be used to
        look up terms from IDs."""

        regex = re.compile("InterPro:([A-Z0-9]+) .* > .* ; (GO:[0-9]+)")
        result = cls()
        for line in open_anything(filename):
            if line[0] == '!':
                continue
            match = regex.match(line)
            if not match:
                continue
            result.add_annotation(match.group(1), tree.lookup(match.group(2)))
        return result


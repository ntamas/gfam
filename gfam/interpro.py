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

__all__ = ["AssignmentParser", "AssignmentReader", \
           "InterPro", "InterProNames", \
           "InterPro2GOMapping"]

class AssignmentParser(object):
    """Parses lines from an InterPro domain assignment file"""

    def __init__(self):
        pass

    def parse(self, line):
        """Parses a line from an InterPro domain assignment file and
        returns the corresponding assignment"""
        parts = line.strip().split("\t")

        try:
            evalue = float(parts[8])
        except ValueError:
            evalue = None

        if len(parts) < 15:
            parts.extend([None] * (15-len(parts)))

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


class AssignmentReader(object):
    """Iterates over assignments in an InterPro domain assignment file"""

    def __init__(self, filename):
        self._fp = open_anything(filename)

    def assignments(self):
        """Yields the assignments in the InterPro domain assignment file
        one by one"""
        parser = AssignmentParser()
        for line in self._fp:
            yield parser.parse(line)

    def __iter__(self):
        return self.assignments()


class InterProTree(Mapping):
    """Dict-like object that tells the parent of every InterPro ID"""

    def __init__(self):
        self._data = {}

    def __contains__(self, item):
        return self._data.__contains__(item)

    def __getitem__(self, item):
        return self._data.get(item, item)

    def get_most_remote_ancestor(self, item):
        """Returns the most remote ancestor of the given item"""
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
    """Dict-like object that maps domain IDs to their corresponding InterPro
    IDs."""

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
    def __init__(self, name):
        self.names = {}
        if name is None:
            return
        self.load(name)

    def load(self, name):
        for line in open_anything(name):
            parts = line.strip().split("\t", 1)
            if len(parts) > 1:
                self.names[parts[0]] = parts[1]

    def __contains__(self, name):
        return name in self.names

    def __getitem__(self, name):
        return self.names.get(name, name)

    @classmethod
    def FromFile(cls, filename):
        return cls(filename)


class InterPro(object):
    def __init__(self):
        self.tree = InterProTree()
        self.mapping = InterProIDMapper()

    @classmethod
    def FromFile(cls, filename):
        result = cls()
        path_to_root = []

        for line in open_anything(filename):
            line = line.strip()
            dash_count = 0
            while line[dash_count] == "-":
                dash_count += 1
            if dash_count % 2 != 0:
                raise ValueError, "dash count in InterPro file not even"

            line = line[dash_count:]
            parts = line.split("::")
            interpro_id, aliases = parts[0], parts[2:]

            level = dash_count / 2 + 1
            if level <= len(path_to_root):
                path_to_root = path_to_root[:level]
                path_to_root[-1] = interpro_id
            else:
                path_to_root.append(interpro_id)
                if level != len(path_to_root):
                    raise ValueError, "tree depth increased by more than " +\
                                      "one between two lines"
            if len(path_to_root) > 1:
                result.tree[interpro_id] = path_to_root[-2]

            for alias in aliases:
                if alias[0:4] == "PTHR" and ":SF" in alias:
                    alias = alias[0:alias.index(":SF")]
                result.mapping[alias] = interpro_id

        return result


class InterPro2GOMapping(bidict):
    """Bidirectional dictionary that tells the corresponding Gene
    Ontology terms of every InterPro ID and vice versa.

    This object encapsulates two dictionaries: `self.terms`
    gives the Gene Ontology terms of a given InterPro domain,
    and `self.domains` gives the InterPro domains annotated by a
    given GO term."""

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
        tree object (see `gene_ontology.Tree`_) that will be used to
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


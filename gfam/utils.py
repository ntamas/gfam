"""Common routines and utility classes for GFam that fit nowhere else."""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["bidict", "complementerset", "open_anything", "redirected",
           "search_file", "temporary_dir", "UniqueIdGenerator"]

try:
    import bz2
except ImportError:
    pass

try:
    import gzip
except ImportError:
    pass

import os
import platform
import sys

from contextlib import contextmanager
from shutil import rmtree
from tempfile import mkdtemp

class bidict(object):
    """Bidirectional many-to-many mapping.
    
    This class models many-to-many mappings that are used in some places in
    GFam. For instance, the mapping of GO term identifiers to InterPro
    domain identifiers is typically a many-to-many mapping (domains are
    annotated by multiple GO terms, and GO terms may belong to multiple
    domains).
    
    Being a general many-to-many mapping class, instances contain two
    member dictionaries: `left` (which maps items from one side of the
    relationship to *sets* of items from the other side of the relationship)
    and `right` (which contains the exact opposite). You should only query
    these dictionaries directly, manipulation should be done by the methods
    provided by the `bidict` class to ensure that the two dicts are kept
    in sync.

    Example::

        >>> bd = bidict()
        >>> bd.add_left("foo", "bar")
        >>> bd.add_left("foo", "baz")
        >>> bd.get_left("foo") == set(['bar', 'baz'])
        True
        >>> bd.add_right("baz", "frob")
        >>> bd.get_right("bar")
        set(['foo'])
        >>> bd.get_right("baz") == set(['foo', 'frob'])
        True
        >>> bd.len_left()
        2
        >>> bd.len_right()
        2
    """

    def __init__(self, items=None):
        object.__init__(self)

        self.left = {}
        self.right = {}

        if items is None:
            return

        if isinstance(items, bidict):
            # Copying an existing bidict
            self.left = dict(items.left)
            self.right = dict(items.right)
        elif isinstance(items, dict):
            # items assigns an element from the left dict to a
            # set of elements from the right dict
            for key, values in items.iteritems():
                self.add_left_multi(key, values)
        else:
            raise TypeError("items must be dict or bidict, got %r" % \
                            type(items))

    def add_left(self, v1, v2):
        """Adds a pair of items `v1` and `v2` to the mapping s.t. `v1` as a
        left item is mapped to `v2` as a right item."""
        try:
            self.left[v1].add(v2)
        except KeyError:
            self.left[v1] = set([v2])
        try:
            self.right[v2].add(v1)
        except KeyError:
            self.right[v2] = set([v1])

    def add_right(self, v1, v2):
        """Adds a pair of items `v1` and `v2` to the mapping s.t. `v1` as a
        right item is mapped to `v2` as a left item."""
        return self.add_left(v2, v1)

    def add_left_multi(self, v1, v2s):
        """Associates multiple items in `v2s` to `v1` when `v1` is interpreted
        as a left item"""
        try:
            self.left[v1].update(v2s)
        except KeyError:
            self.left[v1] = set(v2s)
        for v2 in v2s: 
            try:
                self.right[v2].add(v1)
            except KeyError:
                self.right[v2] = set([v1])

    def add_right_multi(self, v2, v1s):
        """Associates multiple items in `v1s` to `v2` when `v2` is interpreted
        as a right item"""
        self.right[v2].update(v1s)
        for v1 in v1s:
            self.left[v1].add(v2)

    def get_left(self, v1, default = None):
        """Returns the items associated to `v1` when `v1` is looked up from the
        left dictionary. `default` will be returned if `v1` is not in the left
        dictionary."""
        return self.left.get(v1, default)

    def get_right(self, v1, default=None):
        """Returns the items associated to `v1` when `v1` is looked up from the
        right dictionary. `default` will be returned if `v2` is not in the right
        dictionary."""
        return self.right.get(v1, default)

    def len_left(self):
        """Returns the number of unique left items"""
        return len(self.left)

    def len_right(self):
        """Returns the number of unique right items"""
        return len(self.right)

    def iteritems_left(self):
        """Iterates over the left dictionary"""
        return self.left.iteritems()

    def iteritems_right(self):
        """Iterates over the right dictionary"""
        return self.right.iteritems()


def open_anything(fname, *args, **kwds):
    """Opens the given file. The file may be given as a file object
    or a filename. If the filename ends in ``.bz2`` or ``.gz``, it will
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
        import urllib2
        infile = urllib2.urlopen(fname)
    elif fname[-4:] == ".bz2":
        infile = bz2.BZ2File(fname, *args, **kwds)
    elif fname[-3:] == ".gz":
        infile = gzip.GzipFile(fname, *args, **kwds)
    else:
        infile = open(fname, *args, **kwds)
    return infile


# pylint:disable-msg=W0613
# W0613: unused argument.
@contextmanager
def redirected(stdin=None, stdout=None, stderr=None):
    """Temporarily redirects some of the input/output streams.

    `stdin` is the new standard input, `stdout` is the new standard output,
    `stderr` is the new standard error. ``None`` means to leave the
    corresponding stream unchanged.
    
    Example::
    
        with redirected(stdout=open("test.txt", "w")):
            print "Yay, redirected to a file!"
    """
    loc = locals()
    stream_names = ["stdin", "stdout", "stderr"]
    old_streams = {}
    try:
        for sname in stream_names:
            stream = loc.get(sname, None)
            if stream is not None and stream != getattr(sys, sname):
                old_streams[sname] = getattr(sys, sname)
                setattr(sys, sname, loc.get(sname, None))
        yield
    finally:
        for sname, stream in old_streams.iteritems():
            setattr(sys, sname, stream)

def search_file(filename, search_path=None, executable=True):
    """Finds the given `filename` in the given search path.
    If `executable` and we are on Windows, ``.exe`` will be
    appended to the filename. Returns the full path of the
    file or ``None`` if it is not found on the path."""
    if not search_path:
        search_path = os.environ["PATH"]

    if executable and platform.system() == "Windows":
        filename = filename+".exe"

    exists = os.path.exists
    join = os.path.join

    for path in search_path.split(os.pathsep):
        fullpath = join(path, filename)
        if exists(fullpath):
            return os.path.abspath(fullpath)

    return None


@contextmanager
def temporary_dir(*args, **kwds):
    """Context manager that creates a temporary directory when entering the
    context and removes it when exiting.
    
    Every argument is passed on to `tempfile.mkdtemp` except a keyword argument
    named `change` which tells whether we should change to the newly created
    temporary directory or not. The current directory will be restored when
    exiting the context manager."""
    change = "change" in kwds
    if change:
        del kwds["change"]

    name = mkdtemp(*args, **kwds)
    try:
        if change:
            old_dir = os.getcwd()
            os.chdir(name)
        yield name
    finally:
        if change:
            os.chdir(old_dir)
        rmtree(name)


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
        if id_generator is None:
            id_generator = 0
        if isinstance(id_generator, int):
            import itertools
            self._generator = itertools.count(id_generator)
        else:
            self._generator = id_generator
        self._ids = {}

    def __getitem__(self, item):
        """Retrieves the ID corresponding to `item`. Generates a new ID for
        `item` if it is the first time we request an ID for it."""
        try:
            return self._ids[item]
        except KeyError:
            self._ids[item] = self._generator.next()
            return self._ids[item]

    def __len__(self):
        """Retrieves the number of added elements in this `UniqueIDGenerator`"""
        return len(self._ids)

    def reverse_dict(self):
        """Returns the reversed mapping, i.e., the one that maps generated IDs
        to their corresponding items"""
        return dict((v, k) for k, v in self._ids.iteritems())

    def values(self):
        """Returns the list of items added so far. Items are ordered according
        to the standard sorting order of their keys, so the values will be
        exactly in the same order they were added if the ID generator generates
        IDs in ascending order. This hold, for instance, to numeric ID
        generators that assign integers starting from a given number."""
        return sorted(self._ids.keys(), key = self._ids.__getitem__)


class complementerset(object):
    """This object behaves more or less like a set, with one exception,
    the membership checking. For a `complementerset` object, you can
    define the elements which are *not* in the set, everything else is
    contained in it. The semantics of the operators are the same as for
    sets.

    Usage example::

        >>> s = complementerset()
        >>> "abc" in s
        True
        >>> s in s
        True
        >>> s -= set(["abc"])
        >>> s
        complementerset(['abc'])
        >>> "abc" in s
        False
    """

    __slots__ = ("_set", )

    def __init__(self, iterable=()):
        """Constructs a complementer set that contains everything except
        the members of the given iterable."""
        self._set = set(iterable)

    def __and__(self, other):
        """Example::

            >>> s = complementerset()
            >>> s2 = set([1,2,3])
            >>> (s & s2) == s2
            True
            >>> s = complementerset([3,4])
            >>> (s & s2) == set([1,2])
            True
            >>> (s & complementerset([1,2,3])) == complementerset([1,2,3,4])
            True
        """
        if isinstance(other, self.__class__):
            return complementerset(self._set | other._set)

        self._ensure_set(other)
        return set(other) - self._set

    def __contains__(self, what):
        """Example::

            >>> s = complementerset([1,2])
            >>> s in s
            True
            >>> "foo" in s
            True
            >>> 1 in s
            False
        """
        return what not in self._set

    def __eq__(self, other):
        """Example::

            >>> s = complementerset([1,2])
            >>> s == complementerset([1,2])
            True
            >>> s == complementerset([1,2,3])
            False
            >>> s == 1
            False
        """
        return isinstance(other, self.__class__) and \
               self._set == other._set

    def __ge__(self, other):
        """Example::

            >>> s = complementerset()
            >>> s2 = complementerset()
            >>> s2 >= s
            True
            >>> s >= s2
            True
            >>> s >= complementerset([1,2])
            True
            >>> s >= set([1,2,3])
            True
            >>> s >= 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
            >>> complementerset([1,2,3]) >= complementerset([1,2,3])
            True
            >>> complementerset([1,2]) >= complementerset([2,3,4])
            False
        """
        if isinstance(other, self.__class__):
            return self._set <= other._set

        self._ensure_set(other)
        return True

    def __gt__(self, other):
        """Example::

            >>> s = complementerset()
            >>> s2 = complementerset()
            >>> s2 > s
            False
            >>> s > s2
            False
            >>> s > complementerset([1,2])
            True
            >>> s > set([1,2,3])
            True
            >>> s > 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
            >>> complementerset([1,2,3]) > complementerset([1,2,3])
            False
            >>> complementerset([1,2]) > complementerset([1,2,3])
            True
            >>> complementerset([1,2,3]) > complementerset([1,2])
            False
            >>> complementerset([1,2]) > complementerset([2,3,4])
            False
        """
        return self != other and self >= other

    def __hash__(self):
        return hash(self._set)

    def __iand__(self, other):
        """Example::

            >>> s = complementerset([3,4,5])
            >>> s &= s
            >>> s == complementerset([3, 4, 5])
            True
            >>> s &= complementerset([5,6])
            >>> s
            complementerset([3, 4, 5, 6])
            >>> s &= set([1,2,3])
            >>> s == set([1, 2])
            True
            >>> s = complementerset()
            >>> s &= 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        if isinstance(other, self.__class__):
            self._set |= other._set
            return self

        self._ensure_set(other)
        return other - self._set

    def __ior__(self, other):
        """Example::

            >>> s = complementerset([3,4,5])
            >>> s |= s
            >>> s == complementerset([3, 4, 5])
            True
            >>> s |= set([1,2,3])
            >>> s == complementerset([4, 5])
            True
            >>> s |= complementerset([1, 4])
            >>> s
            complementerset([4])
            >>> s |= 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        if isinstance(other, self.__class__):
            self._set &= other._set
            return self

        self._ensure_set(other)
        self._set -= other
        return self

    def __isub__(self, other):
        """Example::

            >>> s = complementerset()
            >>> s -= set([1,2,3])
            >>> s
            complementerset([1, 2, 3])
            >>> s -= 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
            >>> s -= complementerset([1,2,4])
            >>> s
            set([4])
        """
        if isinstance(other, self.__class__):
            return other._set - self._set

        self._ensure_set(other)
        self._set |= other
        return self

    def __le__(self, other):
        """Example::

            >>> s = complementerset()
            >>> s2 = complementerset()
            >>> s2 <= s
            True
            >>> s <= s2
            True
            >>> s <= complementerset([1,2])
            False
            >>> s <= set([1,2,3])
            False
            >>> s <= 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
            >>> complementerset([1,2,3]) <= complementerset([1,2,3])
            True
            >>> complementerset([1,2]) <= complementerset([2,3,4])
            False
        """
        if isinstance(other, self.__class__):
            return self._set >= other._set

        self._ensure_set(other)
        return False

    def __lt__(self, other):
        """Example::

            >>> s = complementerset()
            >>> s2 = complementerset()
            >>> s2 < s
            False
            >>> s < s2
            False
            >>> s < complementerset([1,2])
            False
            >>> complementerset([1,2]) < s
            True
            >>> s < set([1,2,3])
            False
            >>> s < 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
            >>> complementerset([1,2,3]) < complementerset([1,2,3])
            False
            >>> complementerset([1,2]) < complementerset([1,2,3])
            False
            >>> complementerset([1,2,3]) < complementerset([1,2])
            True
            >>> complementerset([1,2]) < complementerset([2,3,4])
            False
        """
        return self != other and self <= other

    def __ne__(self, other):
        return not self == other
    
    def __or__(self, other):
        """Example::

            >>> s = complementerset([2,4,5])
            >>> s | set([1,2,3]) == complementerset([4, 5])
            True
            >>> s | complementerset([1,2])
            complementerset([2])
            >>> s | 2
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        if isinstance(other, self.__class__):
            return complementerset(self._set & other._set)

        self._ensure_set(other)
        return complementerset(self._set - other)

    def __rand__(self, other):
        """Example::

            >>> s = complementerset([2,4,5])
            >>> set([1,2,3]) & s == set([1,3])
            True
            >>> complementerset([1,2,3]) & s == complementerset([1,2,3,4,5])
            True
            >>> 2 & s
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        return self & other

    def __ror__(self, other):
        """Example::

            >>> s = complementerset([2,4,5])
            >>> set([1,2,3]) | s == complementerset([4, 5])
            True
            >>> complementerset([1, 2]) | s
            complementerset([2])
            >>> 2 | s
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        return self | other

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, list(self._set))

    def __rsub__(self, other):
        """Example::

            >>> set([1,2,3]) - complementerset()
            set([])
            >>> (set([1,2,3]) - complementerset([1,2,4])) == set([1,2])
            True
            >>> 2 - complementerset()
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        if isinstance(other, self.__class__):
            raise NotImplementedError

        self._ensure_set(other)
        return other.intersection(self._set)

    def __sub__(self, other):
        """Example::

            >>> complementerset() - set([1,2,3])
            complementerset([1, 2, 3])
            >>> complementerset([1,2,3]) - complementerset([2,4])
            set([4])
            >>> complementerset([3,4]) - set([1,2,3])
            complementerset([1, 2, 3, 4])
            >>> complementerset() - 2 
            Traceback (most recent call last):
              File "<stdin>", line 1, in ?
            NotImplementedError
        """
        if isinstance(other, self.__class__):
            return other._set - self._set

        self._ensure_set(other)
        return complementerset(self._set | other)

    @staticmethod
    def _ensure_set(obj):
        if not isinstance(obj, (set, frozenset, complementerset)):
            raise NotImplementedError

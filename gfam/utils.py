"""Common routines and utility classes for GFam that fit nowhere else."""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["bidict", "complementerset", "Histogram",
           "open_anything", "redirected", "RunningMean",
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
from math import ceil
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


class Histogram(object):
    """Generic histogram class for real numbers
    
    Example:
        
        >>> h = Histogram(5)     # Initializing, bin width = 5
        >>> h << [2,3,2,7,8,5,5,0,7,9]     # Adding more items
        >>> print h
        N = 10, mean +- sd: 4.8000 +- 2.9740
        [ 0,  5): **** (4)
        [ 5, 10): ****** (6)
    """

    def __init__(self, bin_width = 1, data = None):
        """Initializes the histogram with the given data set.

        :param bin_width: the bin width of the histogram.
        :param data: the data set to be used. Must contain real numbers.
        """
        self._bin_width = float(bin_width)
        self._bins = None
        self._min, self._max = None, None
        self._running_mean = RunningMean()
        self.clear()

        if data:
            self.add_many(data)

    def _get_bin(self, num, create = False):
        """Returns the bin index corresponding to the given number.

        :param num: the number for which the bin is being sought
        :param create: whether to create a new bin if no bin exists yet.

        :Returns:
          the index of the bin or ``None`` if no bin exists yet and
          `create` is ``False``."""
        if len(self._bins) == 0:
            if not create:
                result = None
            else: 
                self._min = int(num/self._bin_width)*self._bin_width
                self._max = self._min+self._bin_width
                self._bins = [0]
                result = 0
            return result

        if num >= self._min:
            binidx = int((num-self._min)/self._bin_width)
            if binidx < len(self._bins):
                return binidx
            if not create:
                return None
            extra_bins = binidx-len(self._bins)+1
            self._bins.extend([0]*extra_bins)
            self._max = self._min + len(self._bins)*self._bin_width
            return binidx

        if not create:
            return None

        extra_bins = int(ceil((self._min-num)/self._bin_width))
        self._bins[0:0] = [0]*extra_bins
        self._min -= extra_bins*self._bin_width
        self._max = self._min + len(self._bins)*self._bin_width
        return 0

    @property
    def bin_width(self):
        """Returns the bin width of the histogram"""
        return self._bin_width

    @property
    def n(self):
        """Returns the number of elements in the histogram"""
        return len(self._running_mean)

    @property
    def mean(self):
        """Returns the mean of the elements in the histogram"""
        return self._running_mean.mean

    # pylint: disable-msg=C0103
    @property
    def sd(self):
        """Returns the standard deviation of the elements in
        the histogram"""
        return self._running_mean.sd

    @property
    def var(self):
        """Returns the variance of the elements in the histogram"""
        return self._running_mean.var

    def add(self, num, repeat=1):
        """Adds a single number to the histogram.
        
        :param num: the number to be added
        :param repeat: number of repeated additions
        """
        num = float(num)
        binidx = self._get_bin(num, True)
        self._bins[binidx] += repeat 
        self._running_mean.add(num, repeat)

    def add_many(self, data):
        """Adds a single number or the elements of an iterable to the
        histogram.

        :param data: an iterable containing the data to be added
        """
        try:
            iterator = iter(data)
        except TypeError:
            iterator = iter([data])
        for x in iterator:
            self.add(x)
    __lshift__ = add_many

    def clear(self):
        """Clears the collected data"""
        self._bins = []
        self._min, self._max = None, None
        self._running_mean = RunningMean()

    def bins(self):
        """Generator returning the bins of the histogram in increasing
        order
        
        Returns a tuple with the following elements: left bound, right
        bound, number of elements in the bin"""
        x = self._min
        for elem in self._bins:
            yield (x, x+self._bin_width, elem)
            x += self._bin_width

    def to_string(self, max_width=78, show_bars=True, show_counts=True):
        """Returns the string representation of the histogram.

        `max_width` is the maximal width of each line of the string
        representation. `show_bars` specify whether the histogram bars
        should be shown, `show_counts` specify whether the histogram counts
        should be shown. If both are `False`, only a general descriptive
        statistics (number of elements, mean and standard deviation) is
        shown.
        
        `max_width` may not be obeyed if it is too small.
        """

        if self._min is None or self._max is None:
            return "N = 0"

        # Determine how many decimal digits should we use
        if int(self._min) == self._min and int(self._bin_width) == self._bin_width:
            number_format = "%d"
        else:
            number_format = "%.3f"
        num_length = max(len(number_format % self._min), \
                         len(number_format % self._max))
        number_format = "%" + str(num_length) + number_format[1:]
        format_string = "[%s, %s): %%s" % (number_format, number_format)

        # Calculate the scale of the bars on the histogram
        if show_bars:
            maxval = max(self._bins)
            if show_counts:
                maxval_length = len(str(maxval))
                scale = maxval // (max_width-2*num_length-maxval_length-9)
            else:
                scale = maxval // (max_width-2*num_length-6)
            scale = max(scale, 1)

        result = ["N = %d, mean +- sd: %.4f +- %.4f" % \
            (self.n, self.mean, self.sd)]

        if show_bars:
            # Print the bars
            if scale > 1:
                result.append("Each * represents %d items" % scale)
            if show_counts:
                format_string += " (%d)"
                for left, right, cnt in self.bins():
                    result.append(format_string % (left, right, '*'*(cnt//scale), cnt))
            else:
                for left, right, cnt in self.bins():
                    result.append(format_string % (left, right, '*'*(cnt//scale)))
        elif show_counts:
            # Print the counts only
            for left, right, cnt in self.bins():
                result.append(format_string % (left, right, cnt))

        return "\n".join(result)

    def __str__(self):
        return self.to_string()

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
    """Context manager that temporarily redirects some of the input/output
    streams.

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


class RunningMean(object):
    """Running mean calculator.
    
    This class can be used to calculate the mean of elements from a
    list, tuple, iterable or any other data source. The mean is
    calculated on the fly without explicitly summing the values,
    so it can be used for data sets with arbitrary item count. Also
    capable of returning the standard deviation (also calculated on
    the fly)
    """

    # pylint: disable-msg=C0103
    def __init__(self, n = 0.0, mean = 0.0, sd = 0.0):
        """RunningMean(n=0.0, mean=0.0, sd=0.0)
        
        Initializes the running mean calculator. Optionally the
        number of already processed elements and an initial mean
        can be supplied if we want to continue an interrupted
        calculation.

        @param n: the initial number of elements already processed
        @param mean: the initial mean
        @param sd: the initial standard deviation"""
        self._nitems = float(n)
        self._mean = float(mean)
        if n > 1:
            self._sqdiff = float(sd) ** 2 * float(n-1)
            self._sd = float(sd)
        else:
            self._sqdiff = 0.0
            self._sd = 0.0
        
    def add(self, value, repeat=1):
        """RunningMean.add(value, repeat=1)
        
        Adds the given value to the elements from which we calculate
        the mean and the standard deviation.

        @param value: the element to be added
        @param repeat: number of repeated additions
        @return: the new mean and standard deviation as a tuple"""
        repeat = int(repeat)
        self._nitems += repeat
        delta = value - self._mean
        self._mean += (repeat*delta / self._nitems)
        self._sqdiff += (repeat*delta) * (value - self._mean)
        if self._nitems > 1:
            self._sd = (self._sqdiff / (self._nitems-1)) ** 0.5
        return self._mean, self._sd

    def add_many(self, values):
        """RunningMean.add(values)
        
        Adds the values in the given iterable to the elements from
        which we calculate the mean. Can also accept a single number.
        The left shift (C{<<}) operator is aliased to this function,
        so you can use it to add elements as well:
            
          >>> rm=RunningMean()
          >>> rm << [1,2,3,4]           # doctest:+ELLIPSIS
          (2.5, 1.290994...)
        
        @param values: the element(s) to be added
        @type values: iterable
        @return: the new mean"""
        try:
            iterator = iter(values)
        except TypeError:
            iterator = iter([values])
        for value in iterator:
            self.add(value)
        return self._mean, self._sd

    @property
    def result(self):
        """Returns the current mean and standard deviation as a tuple"""
        return self._mean, self._sd

    @property
    def mean(self):
        """Returns the current mean"""
        return self._mean

    @property
    def sd(self):
        """Returns the current standard deviation"""
        return self._sd

    @property
    def var(self):
        """Returns the current variation"""
        return self._sd ** 2

    def __str__(self):
        return "Running mean (N=%d, %f +- %f)" % \
            (self._nitems, self._mean, self._sd)
    
    __lshift__ = add_many
    
    def __float__(self):
        return float(self._mean)

    def __int__(self):
        return int(self._mean)

    def __long__(self):
        return long(self._mean)

    def __complex__(self):
        return complex(self._mean)

    def __len__(self):
        return self._nitems


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

    def difference_update(self, *args):
        """Removes all elements of another set from this set.

        Example::

            >>> s = complementerset([1,2])
            >>> s.difference_update([4,5])
            >>> print s
            complementerset([1, 2, 4, 5])
            >>> s.difference_update([2], [1,6], [7,5,"spam"])
            >>> print any(item in s for item in [2,1,6,7,5,"spam",4])
            False
        """
        self._set.update(*args)

    def discard(self, member):
        """Removes an element from the complementer set if it is a member.

        Example::

            >>> s = complementerset()
            >>> s.discard(2)
            >>> print s
            complementerset([2])
            >>> s.discard(2)
            >>> print s
            complementerset([2])
        """
        self._set.add(member)

    def iterexcluded(self):
        """Iterates over the items excluded from the complementer set.
        
        Example::
            
            >>> s = complementerset([5, 7, 4])
            >>> print sorted(list(s.iterexcluded()))
            [4, 5, 7]
        """
        return iter(self._set)

    def remove(self, member):
        """Removes an element from the complementer set; it must be a member.
        If the element is not a member, raises a ``KeyError``.
        
        Example::
            
            >>> s = complementerset()
            >>> s.remove(2)
            >>> print s
            complementerset([2])
            >>> s.remove(2)
            Traceback (most recent call last):
              File "<stdin>", line 4, in ?
            KeyError: 2
        """
        if member in self._set:
            raise KeyError(member)
        self._set.add(member)

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
        if hasattr(what, "__hash__"):
            return what not in self._set
        else:
            return True

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

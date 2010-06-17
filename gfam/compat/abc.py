"""Implementation of some of the abstract base classes
introduced in Python 2.6.

This module is not meant to be a complete replacement for
all the abstract base classes in Python 2.6, only the parts
necessary in order to make ``gfam`` compatible with Python 2.5
are implemented.

The API of the classes are identical to the ones in Python 2.6,
but they have been reimplemented from scratch to avoid licensing
issues (``_abcoll.py`` in Python 2.6 is implemented by Google and
licensed to the Python Software Foundation under a special
agreement).
"""

__all__ = ["Mapping"]


class Sized(object):
    def __len__(self):
        raise NotImplementedError

class Iterable(object):
    def __iter__(self):
        raise NotImplementedError

class Container(object):
    def __contains__(self, x):
        raise NotImplementedError

class Mapping(Sized, Iterable, Container):
    """An abstract base class similar to `collections.Mapping` in Python 2.6.
    This will be used in place of `collections.Mapping` if necessary."""

    def __getitem__(self, key):
        raise KeyError

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def items(self):
        return [(key, self[key]) for key in self]

    def iteritems(self):
        for key in self:
            yield key, self[key]

    def iterkeys(self):
        return iter(self)

    def itervalues(self):
        for key in self:
            yield self[key]

    def keys(self):
        return list(self)

    def values(self):
        return [self[key] for key in self]

    def __contains__(self, key):
        try:
            self[key]
            return True
        except KeyError:
            return False

    def __eq__(self, other):
        return isinstance(other, Mapping) and \
               dict(self.items()) == dict(other.items())

    def __hash__(self):
        raise NotImplementedError

    def __ne__(self, other):
        return not isinstance(other, Mapping) or \
               dict(self.items()) != dict(other.items())



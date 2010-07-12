# -*- coding: UTF-8 -*-
"""
Simple enumeration class and metaclass.
"""

__all__ = ["Enum"]

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"


def first(iter):
    """Helper function that takes an iterable and returns the
    first element. No big deal, but it makes the code readable
    in some cases. It is typically useful when `iter` is a
    generator expression as you can use the indexing operator
    for ordinary lists and tuples."""
    for item in iter:
        return item
    raise ValueError("iterable is empty")

class EnumMeta(type):
    """Metaclass for enum classes.

    Do not use this class directly, there's no need to do that.
    Derive your enum class from `Enum` instead.
    """

    def __init__(cls, name, bases, attrs):
        # Call the super-metaclass constructor
        super(EnumMeta, cls).__init__(name, bases, attrs)

        # This dict will contain the enum values
        enum_values = {}

        # For all the base classes, fetch the inherited items
        for base in bases:
            if hasattr(base, "__enum__"):
                enum_values.update(base.__enum__)

        # Extend enum_values with the items directly declared here
        for key, value in attrs.iteritems():
            # Skip internal methods, properties and callables
            if key[:2] != "__" and not hasattr(value, "__call__") \
                    and not isinstance(value, property):
                inst = cls(key, value, override=True)
                enum_values[key] = inst
                super(EnumMeta, cls).__setattr__(key, inst)

        # Store enum_values in the class
        cls.__enum__ = enum_values


    def __setattr__(self, name, value):
        """Raises an `AttributeError` to prevent users from messing
        around with enum values"""
        if name == "__enum__" or name == "__doc__" \
                or not self._finalized:
            return super(EnumMeta, self).__setattr__(name, value)
        raise AttributeError("Enum attributes are read-only")

    def __delattr__(self, name, value):
        """Raises an `AttributeError` to prevent users from messing
        around with enum values"""
        raise AttributeError("Enum attributes cannot be deleted")

    def __repr__(self):
        return "<Enum '%s'>" % self.__name__

    def __len__(self):
        return len(self.__enum__)

    def __iter__(self):
        return iter(self.__enum__)

    def __contains__(self, key):
        return key in self.__enum__

    def from_name(self, name):
        """Constructs an instance of this enum from its name"""
        try:
            return self.__enum__[name]
        except KeyError:
            raise NameError("no enum item with the given name: %r" % name)

    def from_value(self, value):
        """Constructs an instance of this enum from its value"""
        try:
            return first(val for val in self.__enum__.itervalues() \
                         if val.value == value)
        except ValueError:
            raise ValueError("no enum item with the given value: %r" % value)

    def has_key(self, key):
        """Returns whether this enum has the given key or not"""
        return self.__enum__.has_key(key)

    def iteritems(self):
        """Returns an iterator over key-value pairs of this enum"""
        return self.__enum__.iteritems()

    def itervalues(self):
        """Returns an iterator over key-value pairs of this enum"""
        return self.__enum__.itervalues()

    def keys(self):
        """Returns the keys in this enum"""
        return self.__enum__.keys()

    def values(self):
        """Returns the values in this enum"""
        return self.__enum__.values()


class Enum(object):
    """An instance of an enumeration value and a class representing a
    whole enum.

    This is mainly used as a superclass for enumerations. There is
    a clear distinction between using the class itself or using one
    of its instances. Using the class means that you are referring to
    the enum as a whole (with all its possible keys and values).
    Using one of the instances means that you are using a single
    key-value pair from the enum. Instances should never be created
    directly, as all the valid instances are accessible as attributes
    of the class itself.

    Usage example::

        >>> class Spam(Enum):
        ...     SPAM  = "Spam spam spam"
        ...     EGGS  = "Eggs"
        ...     BACON = "Bacon"

    After you have defined an enum class like the one above, you
    can make use of it this way::

        >>> Spam.BACON
        Spam.BACON
        >>> Spam.BACON.value
        'Bacon'

    Think about enums as Python dictionaries that map symbolic names to
    values. Enums even provide methods similar to the non-mutating methods
    of Python dictionaries::

        >>> sorted(Spam.keys())
        ['BACON', 'EGGS', 'SPAM']
    
    You can also get an instance of the enum from its symbolic name or value:

        >>> Spam.from_name("BACON")
        Spam.BACON
        >>> Spam.from_value("Spam spam spam")
        Spam.SPAM
    """

    __metaclass__ = EnumMeta
    __slots__ = ("_key", "_value")

    def __init__(self, key, value, **kwds):
        if "override" not in kwds:
            raise TypeError("Enums should not be instantiated directly")
        self._key = key
        self._value = value

    def __repr__(self):
        """Returns an executable representation of this instance.

        Strictly speaking, this method should return a string that allows
        one to construct this instance, but we don't want anyone to
        start instantiating `Enum`s by hand, so we return an expression
        that refers to this very same instance instead.
        """
        return "%s.%s" % (self.__class__.__name__, self._key)

    def __str__(self):
        return "%s.%s" % (self.__class__.__name__, self._key)

    @property
    def value(self):
        return self._value

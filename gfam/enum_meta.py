###############################################################################
# two_way_dict python module
# Copyright (c) 2005, RADLogic Pty. Ltd. 
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#    * Neither the name of RADLogic nor the names of its contributors
#      may be used to endorse or promote products derived from this
#      software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
###############################################################################
"""Provide enumeration classes using an Enum metaclass.

See the Enum class documentation for details.

This requires Python >= 2.2.

After writing this module I discovered this:
http://www.python.org/doc/essays/metaclasses/Enum.py

It is anonymous, and takes a slightly different approach, making the enum
values instances of another class to enforce certain restrictions.

The thing I grabbed from it was the instantiation of the metaclass, so
that things could subclass from Enum, rather than use the ugly __metaclass__.

Also, as a result of doing that, extendable enumerations came almost for free!

Key differences in this version:
- enumeration classes cannot be modified after they have been created
- a dictionary-style interface is made available to access enumeration values
- no restriction on base classes (might not be a good idea though)
- no special enumeration value class

"""

__author__ = 'Tim Wegener <twegener@radlogic.com.au>'
__version__ = '$Revision: 0.3 $'
__date__ = '$Date: 2005/09/20 08:35:51 $'
__credits__ = '''http://www.python.org/doc/essays/metaclasses/Enum.py
(For metaclass instantiation trick, but overall idea arrived at independently.)
'''

# Here is the magic code snippet that explains what goes on with metaclasses:
# (from http://www.python.org/2.2.3/descrintro.html#metaclasses)
#
# cls = M.__new__(M, name, bases, dict)
# assert cls.__class__ is M
# M.__init__(cls, name, bases, dict)
#
#
# Just remember:
# - a class is an object
# - every object has a class (at least I think!)
# - the metaclass is the class of a class object


class EnumMeta(type):
    """Metaclass for enum classes that keeps track of class variables.

    A metaclass is used, so that the class does not need to be instantiated.
    (In fact, it should not be instantiated.)

    All you need to know to use this is summed up with an example:

    >>> class Symbol(Enum):
    ... 
    ...     A = 'A'
    ...     B = 'B'
    ...     A_B = 'A/B' 
    ...     PLUS = '+' 

    Then use the Symbol enumeration like this:

    >>> symbol = Symbol.A
    >>> symbol == Symbol.A
    1
    >>> symbol == Symbol.B
    0

    >>> Symbol.B
    'B'
    >>> Symbol.A_B
    'A/B'
    >>> Symbol.PLUS
    '+'

    >>> values = Symbol.values()
    >>> values.sort()
    >>> values
    ['+', 'A', 'A/B', 'B']

    >>> Symbol.D = 'D'
    Traceback (most recent call last):
    AttributeError: Cannot modify enumeration attributes.

    Not sure if this is good or not, but this checks the key rather than the
    value.
    >>> print 'A' in Symbol
    1
    >>> print '+' in Symbol
    0

    >>> len(Symbol)
    4

    >>> repr(Symbol)
    "<Enum 'Symbol'>"

    >>> str(Symbol)
    "<Enum 'Symbol'>"
    
    Note: You cannot have multiple enumeration attributes with the same value.
    >>> class Bad(Enum):
    ... 
    ...     A = 'A'
    ...     B = 'A'
    Traceback (most recent call last):
    ValueError: Enumeration value used more than once: 'A'

    Enumerations can be subclassed (inheritance):

    >>> class ExtendedSymbol(Symbol):
    ... 
    ...     C = 'C'

    >>> symbol_ext = ExtendedSymbol.A
    >>> symbol_ext == ExtendedSymbol.A
    1
    >>> symbol_ext == Symbol.A
    1
    >>> ExtendedSymbol.A == Symbol.A
    1
    >>> ExtendedSymbol.C
    'C'

    >>> len(ExtendedSymbol)
    5

    Enumerations can be merged (multiple-inheritance):

    >>> class Men(Enum):
    ... 
    ...     JACK = 'Jack'
    ...     PETE = 'Pete'
    ...     GREG = 'Greg'
    
    >>> class Women(Enum):
    ... 
    ...     JILL = 'Jill'
    ...     PETA = 'Peta'
    ...     JOAN = 'Joan'
    
    >>> class People(Men, Women):
    ...     pass
    
    >>> People.JACK == Men.JACK
    1
    >>> People.JILL == Women.JILL
    1

    >>> len(Men)
    3
    >>> len(People)
    6
    
    """

    def __init__(class_, name, bases, dict_):

        # Call the *bound* __init__ method of the super class (which in this
        # case is type)
        super(EnumMeta, class_).__init__(name, bases, dict_)

        # Store the non-special items in a class-wide __enum__ dict.
        enum_dict = {}

        # Grab the inherited items first.
        # Note: This relies on Enum subclasses not changing once they have
        #       have been created.
        source_dicts = []
        for base in bases:
            if issubclass(base, Enum) and base is not Enum:
                source_dicts.append(base.__enum__)

        source_dicts.append(dict_)

        for source_dict in source_dicts:
            for k, v in source_dict.items():
                if k[:2] == '__':
                    continue
                # Check for duplicate values.
                if v in enum_dict.values():
                    raise ValueError(
                        "Enumeration value used more than once: %r" % v)
                enum_dict[k] = v

        # Set the enum dict as a class variable.
        class_.__enum__ = enum_dict

    # Prevent new enumeration values from being added.
    def __setattr__(self, name, val):

##         if hasattr(self, '__enum__'):
##             raise AttributeError('Cannot modify enumeration attributes.')

        # Have to check this way to allow subclassing Enum subclasses.
        # This means that users are able to fiddle with the __enum__ dict,
        # but the double-underscores warn them not to. 
        if name != '__enum__':
            raise AttributeError('Cannot modify enumeration attributes.')

        type.__setattr__(self, name, val)

    # Define special methods here, to expose the __enum__ dictionary methods
    # to classes that have Enum as their metaclass.

    # Note: We don't want __getitem__, __setitem__, __delitem__,
    #       clear, update, get, setdefault, pop, popitem, fromkeys.

    def __repr__(self):

        # todo: Not sure what is best for this. Should it include the values?
        return "<Enum '%s'>" % (self.__name__,)
##                                ','.join(self.__enum__.keys()))

    def __len__(self):

        return len(self.__enum__)

    def __iter__(self):

        return iter(self.__enum__)

    def values(self):

        return self.__enum__.values()

    def keys(self):

        return self.__enum__.keys()

    def items(self):

        return self.__enum__.items()

    def has_key(self, key):

        return self.__enum__.has_key(key)

    def itervalues(self):

        return self.__enum__.itervalues()

    def iteritems(self):

        return self.__enum__.iteritems()


# Make a class out of the metaclass, so that it can be subclassed.
# Grabbed this trick from:
# http://www.python.org/doc/essays/metaclasses/Enum.py
Enum = EnumMeta('Enum', (), {})


if __name__ == '__main__':
    # Run doctests.
    import sys
    import doctest
    module = sys.modules.get('__main__')
    doctest.testmod(module)


    # Debug
##     class Symbol(Enum):
        
##         A = 'A'
##         B = 'B'
##         A_B = 'A/B'


##     class ExtSymbol(Symbol):

##         C = 'C'

##     print ExtSymbol.C
##     print ExtSymbol.A
##     print ExtSymbol.A_B
    

##     Symbol.D = 2

##     print Symbol.D


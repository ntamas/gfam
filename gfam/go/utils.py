"""Utility classes and routines that do not fit anywhere else"""

class ParseError(Exception):
    """Exception thrown when a parsing error occurred"""

    def __init__(self, msg, lineno = 1):
        Exception.__init__(self, "%s near line %d" % (msg, lineno))
        self.lineno = lineno

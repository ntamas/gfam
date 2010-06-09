"""
Hashing functions for Modula
"""

import sys

if sys.version_info >= (2, 5):
    from hashlib import sha1
else:
    from sha import new as sha1

__all__ = ["sha1"]


"""This is the main module of GFam.

On its own, this module contains nothing, all the functionality is implemented
in one of the following submodules:
    
:mod:`gfam.assignment`
  Contains routines related to assignments, i.e. the objects describing the fact
  that a given section of a sequence is assigned to a given domain.

:mod:`gfam.blast`
  Classes and functions related to handling BLAST output files and external BLAST
  utilities.

:mod:`gfam.compat`
  Classes and functions to maintain compatibility with Python 2.5.

:mod:`gfam.config`
  An extension of Python's built-in :mod:`optparse` module to allow supplying
  default values for command line options from a configuration file.

:mod:`gfam.enum`
  A simple enumeration class using some metaclass magic. Sadly enough, Python
  does not have a built-in and flexible enumeration class like Java does.

:mod:`gfam.fasta`
  A simple parser and emitter for the FASTA sequence format. Yes, we could
  have used BioPython_ instead, but this way we saved a dependency.

:mod:`gfam.go`
  Classes and functions related to the Gene Ontology.

:mod:`gfam.interpro`
  Classes and functions related to parsing InterPro files, especially the
  output of ``iprscan``.

:mod:`gfam.modula`
  A simple module that manages tasks and dependencies of the GFam pipeline.
  It will ensure that a GFam execution does not run long calculations
  unnecessarily if the same results are available from a previous interrupted
  GFam run.

:mod:`gfam.scripts`
  Implementations for the main steps of the GFam pipeline. Each submodule of
  this module can be executed on its own as a command-line utility.

:mod:`gfam.sequence`
  A simple replacement for the ``Sequence`` and ``SeqRecord`` classes of
  BioPython_.

:mod:`gfam.utils`
  Various utility routines that did not fit anywhere else.

.. _BioPython: http://www.biopython.org

General comments that apply for the whole GFam API:

- Routines that accept files usually accept either filenames or file-like
  objects. If the filename ends in ``.bz2`` or ``.gz``, it will be
  decompressed on-the-fly in memory. If the filename starts with ``http://``,
  ``https://`` or ``ftp://``, it is assumed to be remote object and will
  be downloaded accordingly. This is achieved by `gfam.utils.open_anything`,
  which is called whenever a filename or a file-like object is passed into
  a function.
"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

__version__ = "1.0"

import gfam.fasta as fasta

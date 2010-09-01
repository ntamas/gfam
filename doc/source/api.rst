API documentation
=================

.. default-role:: obj

:mod:`gfam` -- The main module
------------------------------

.. automodule:: gfam

:mod:`gfam.assignment` -- Routines related to sequence-domain annotations
-------------------------------------------------------------------------

.. automodule:: gfam.assignment
   :members:

:mod:`gfam.blast` -- Handling BLAST file formats and utilities
--------------------------------------------------------------

.. automodule:: gfam.blast
   :members:

:mod:`gfam.compat` -- Compatibility classes for Python 2.5
----------------------------------------------------------

.. automodule:: gfam.compat
   :members:

:mod:`gfam.config` -- Configuration file handling
-------------------------------------------------

.. automodule:: gfam.config
   :members:

:mod:`gfam.enum` -- A simple enumeration class
----------------------------------------------

.. automodule:: gfam.enum
   :members:

:mod:`gfam.fasta` -- FASTA parser and emitter
---------------------------------------------

.. automodule:: gfam.fasta
   :members:

:mod:`gfam.go` -- Handling the Gene Ontology
--------------------------------------------

.. automodule:: gfam.go
   :members:

:mod:`gfam.go.obo` -- Parsing OBO ontology files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: gfam.go.obo
   :members:

:mod:`gfam.go.overrepresentation` -- Overrepresentation analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: gfam.go.overrepresentation
   :members:

:mod:`gfam.interpro` -- Handling InterPro-related files
-------------------------------------------------------

.. automodule:: gfam.interpro
   :members:

:mod:`gfam.modula` -- Modular calculation framework
---------------------------------------------------

Modula is a modular calculation framework for Python that allows you to
define tasks that depend on input files and on each others. Modula will
figure out which tasks have to be executed in which order in order to
calculate the final results -- this is done by a simple depth first
search on the task dependency graph.

Modula started out as a separate project, and you don't have to know its
internals in order to use the GFam API. The only reason why it has been
placed as a submodule of GFam is to avoid foricng users to install Modula
separately. The Modula API is not documented here as it is not an internal
part of GFam. Modula is used only by the GFam master script (see
:mod:`gfam.scripts.master`) to execute the calculation steps in the
proper order.

:mod:`gfam.sequence` -- Simple sequence and sequence record classes
-------------------------------------------------------------------

.. automodule:: gfam.sequence
   :members:

:mod:`gfam.scripts` -- Command line scripts
-------------------------------------------

.. automodule:: gfam.scripts
   :members:

Each step in the `GFam pipeline <pipeline>`_ is implemented in a separate
submodule of :mod:`gfam.scripts`. These submodules contain only a single
class per submodule, derived from `gfam.scripts.CommandLineApp`. The
submodules are invoked automatically in the right order by a master
script in `gfam.scripts.master`. In general, you only have to invoke
the master script and it will do the rest for you, but some of the
steps might be useful on their own, so they can be invoked independently
from the command line as:

.. code-block:: sh

   $ python -m gfam.scripts.modulename

where *modulename* is the name of the submodule to be executed. You can
get usage information for each submodule by typing:

.. code-block:: sh

   $ python -m gfam.scripts.modulename --help

:mod:`gfam.utils` -- Utility classes and functions
--------------------------------------------------

.. automodule:: gfam.utils
   :members:



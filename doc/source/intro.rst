Introduction
============

What is GFam?
-------------

GFam (or ``gfam``) is a Python module to aid the automatic annotation of gene
families based on consensus domain architecture. ``gfam`` started out as a
collection of loosely coupled Python scripts that process the output of
``iprscan`` (a tool to obtain domain assignments of individual genes from
InterPro) and conduct some analyses using BLAST to detect novel, previously
uncharacterised domains. The original domains and the detected novel domain
candidates are then used to create a consensus domain assignment for each gene
sequence. Genes are then finally assigned to families based on their domain
architectures, and a Gene Ontology overrepresentation analysis is conducted on
the GO annotations of individual domains in the same sequence to come up with a
set of functional labels for each sequence.

Requirements
------------

You will need the following tools to run ``gfam``:

* `Python 2.5`_ or latest. Python 3 is not supported yet.

* `NCBI BLAST`_; in particular, the ``formatdb`` and ``blastall`` tools
  from the legacy C-based BLAST distribution.

.. _`Python 2.5`: http://www.python.org
.. _`NCBI BLAST`: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST

The latest release of `SciPy`_ is recommended, but not necessary.
``gfam`` uses `SciPy`_ for calculating the logarithm of the gamma
function in the overrepresentation analysis routines, but it falls
back to a (somewhat slower) Python implementation if `SciPy`_ is
not installed.

.. _`SciPy`: http://www.scipy.org





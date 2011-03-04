Introduction
============

What is GFam?
-------------

GFam (or ``gfam``) is a Python module to aid the automatic annotation of gene
families based on consensus domain architectures. ``gfam`` started out as a
collection of loosely coupled Python scripts that process the output of
``iprscan`` (a tool to obtain domain assignments of individual genes from
InterPro) and conduct some analyses using BLAST to detect novel, previously
uncharacterised domains. The original domains and the detected novel domain
candidates are then used to create a consensus domain assignment for each gene
sequence. Genes are then finally assigned to families based on their domain
architectures. Finally, the tool derives functional labels for families based
on the Gene Ontology and an assignment between InterPro domains and Gene
Ontology terms. Optionally, a Gene Ontology overrepresentation analysis can also
be conducted on the GO annotations of individual domains in the same sequence to
reinforce the functional labels.

Requirements
------------

You will need the following tools to run ``gfam``:

* Python_ 2.5 or later. Python 3 is not supported yet. ``gfam`` was also
  tested with Jython_ 2.5.1.

* `NCBI BLAST`_; in particular, the ``formatdb`` and ``blastall`` tools
  from the legacy C-based BLAST distribution. You can still use ``gfam``
  wit the newer, C++-based BLAST if you have the ``legacy_blast.pl``
  wrapper script in the BLAST folder.

.. _Python: http://www.python.org
.. _`NCBI BLAST`: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST
.. _Jython: http://www.jython.org

The latest release of `SciPy`_ is recommended, but not necessary.
``gfam`` uses `SciPy`_ for calculating the logarithm of the gamma
function in the overrepresentation analysis routines, but it falls
back to a (somewhat slower) Python implementation if `SciPy`_ is
not installed.

.. _`SciPy`: http://www.scipy.org

For the impatient
-----------------

.. highlight:: sh

``gfam`` is driven by a master configuration file named ``gfam.cfg``.
A sample configuration file is given in the distribution. The sample
file works fine for the gene sequences of *Arabidopsis thaliana*; for
other species, you might have to tweak some of the parameters, and you
will surely have to modify the paths to the data files. The configuration
file is documented and mostly self-explanatory.

If you do not have a configuration file for some reason, or you want to
generate a new one from scratch, you can ask ``gfam`` to do it::

    $ bin/gfam init

This will create a default configuration file named ``gfam.cfg`` (if it does
not exist already) and lists the configuration options you have to set in
the file before starting GFam.

If the configuration file is well in order, you can launch ``gfam`` by typing::

    $ bin/gfam

This will run the whole ``gfam`` analysis pipeline using the configuration
specified in ``gfam.cfg``. If your configuration file is named otherwise,
you can run it by typing::

    $ bin/gfam -c my_config.cfg

The results will be put into whatever work directory you specified in the
configuration file. By default, this is named ``work``. See :ref:`output-files`
for more details on what will be calculated and where you can find them.

Questions, comments
-------------------

If you have a question or a comment about ``gfam`` or you think you have
found a bug, feel free to contact me using the email address given in the
header of this document.



Running GFam
============

Input files
-----------

The input files can be grouped into three large groups: *data files*, *mapping
files* and the *configuration file*.  Data files contain the actual input data
that is specific to a given organism. Mapping files usually map between IDs of
different data sources (for instance, from InterPro domain IDs to Gene Ontology
terms) or IDs to human-readable descriptions. The :ref:`configuration file
<config-file>` tells GFam where to find the data files and the mapping files.
When one wants to process a new organism with GFam, it is therefore usually
enough to replace the paths of the data files in the configuration only, as the
mapping files can be re-used for multiple analyses.

GFam requires the following data files:

**Sequence file**
    This file must contain all the sequences that are being analysed and
    annotated by GFam.

**Domain assignment file**
    This file is produced by running ``iprscan``, the command-line variant of
    `InterProScan`_ and it assigns sections of each segment in the sequence
    file to known domains in `InterPro`_. The sequence IDs in this file must be
    identical to the ones in the sequence file; if not, one can specify a
    regular expression in the :ref:`configuration file <config-file>` to
    extract the sequence ID from the FASTA defline.

.. _InterProScan: http://www.ebi.ac.uk/Tools/InterProScan
.. _InterPro: http://www.ebi.ac.uk/interpro

Besides the data files, the following mapping files are also needed:

**InterPro -- GO mapping**
    This file maps InterPro IDs to their corresponding GO terms, and it
    can be obtained from <http://www.geneontology.org/external2go/interpro2go>. 

**Mapping of domain IDs to human-readable names**
    Fairly self-explanatory; a tab-separated flat file with two columns, the
    first being the domain ID and the second being the corresponding
    human-readable name. It is advisable to construct a file which contains at
    least the InterPro, Pfam, SMART and Superfamily IDs as these are the most
    common (and many Pfam, SMART and Superfamily IDs do not have corresponding
    InterPro IDs yet). If you want to create such a mapping file easily, please
    refer to :ref:`Updating the mapping of IDs to human-readable names
    <updating-mappings>`.

**Parent-child relationships of InterPro terms**
    This file contains the parent/child relationships between InterPro
    accession numbers to indicate family/subfamily relationships. This file
    is used to map each InterPro subfamily ID to the corresponding family
    ID, and it can be obtained from `EBI`_.

**The Gene Ontology**
    This file contains the `Gene Ontology`_ in OBO format, and it is
    required only for the `label assignment`_ and `overrepresentation analysis`_
    steps. The latest version of the file can be obtained from the homepage of
    the `Gene Ontology`_ project.

.. _EBI: ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt
.. _Gene Ontology: http://www.geneontology.org
.. _label assignment: pipeline-step-label-assignment
.. _overrepresentation analysis: pipeline-step-overrep

GFam accepts uncompressed files or files compressed with ``gzip`` or ``bzip2``
for both the data and the mapping files. Compressed files will be decompressed
on-the-fly in memory when needed.

.. _config-file:

The configuration file
----------------------

The default configuration file of GFam is called ``gfam.cfg``, but you can
specify an alternative configuration file name on the command line using the
``-c`` switch. A sample configuration file is included in the GFam distribution;
however, you can always generate a new one by running the following command::

    $ bin/gfam init

This will generate a file named ``gfam.cfg`` in the current directory and list
the configuration keys you have to modify before starting your analyses.

The configuration file consists of sections, led by a ``[section]`` header and
followed by ``name=value`` entries. Lines beginning with ``#`` or ``;`` are
ignored and used to provide comments. Lines containing whitespace characters
only are also ignored. For more details about the configuration file format,
please refer to the `ConfigParser module`_ in the documentation of Python.

.. _ConfigParser module: http://docs.python.org/library/configparser.html

The full list of supported configuration keys and their default values is as
follows:

.. include:: generated/config.inc

.. _output-files:

Output files
------------

GFam produces four output files in the output folder specified in the
`configuration file`. These files are as follows:

``domain_architectures.tab``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A simple tab-separated flat file that contains the inferred domain architecture
for each sequence in a simple, summarised format. The file is sorted in a way
such that more frequent domain architectures are placed at the top. Sequences
having the same domain architecture are sorted according to their IDs.

The file has six columns. The first column is the ID of the sequence (e.g.,
``AT1G09650.1``), the second is the sequence length (e.g, ``382``). The third
column contains a summary of the domain architecture of the sequence, where
domains are ordered according to the starting position, and consecutive domain
IDs are separated by semicolons (e.g., ``IPR022364;IPR017451``). The InterPro
domain ID is used whenever possible. Novel domains identified by GFam are
denoted by ``NOVELxxxxx``, where ``xxxxx`` is a five-digit identifier.  The
fourth column contains the frequency of this domain architecture (i.e. the
number of sequences that have the same domain architecture). The fifth column
is the same as the third, but the exact starting and ending positions of the
domain are also added in parentheses after the domain ID (e.g.,
``IPR022364(9-57);IPR017451(112-357)``). The sixth column contains the
concatenated human-readable descriptions of the domains (for instance, ``F-box
domain, Skp2-like;F-box associated interaction domain``).

``domain_architecture_details.txt``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. highlight:: none

This file is the human-readable variant of ``domain_architectures.tab`` (which
is more suitable for machine parsing). It contains blocks separated by two
newline characters; each block corresponds to a sequence and has the following
format::

    AT1G09650.1
        Primary assignment source: HMMTigr
        Number of data sources used: 2
        Data sources: superfamily, HMMTigr
        Coverage: 0.772
        Coverage w/o novel domains: 0.772
           9-  57: SSF81383 (superfamily, stage: 2) (InterPro ID: IPR022364)
                   F-box domain, Skp2-like
         112- 357: TIGR01640 (HMMTigr, stage: 1) (InterPro ID: IPR017451)
                   F-box associated interaction domain

The first line of each block is unindentend and contains the sequence ID. The
remaining lines are indented by at least four spaces. The second line contains
the name of the InterPro data source that was used to come up with the primary
assignment in :ref:`step 2 of the pipeline <pipeline-step-preliminary>` (see
more details later in :ref:`pipeline`), followed by the number of data sources
used to construct the final assignment, and of course the data sources
themselves. The fifth and sixth lines contain the fraction of positions in the
sequence that are covered by at least one domain; the fifth line takes into
account novel domains (``NOVELxxxxx``), while the sixth line does not. The
remaining lines list the domains themselves along with the data source they
came from and the stage in which they were selected.  For more details about
the stages, see :ref:`pipeline`.

``assigned_labels.txt``
^^^^^^^^^^^^^^^^^^^^^^^

TODO

``overrepresentation_analysis.txt``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains the results of the Gene Ontology overrepresentation analysis
for the domain architecture of each sequence. Note that since the results of
the overrepresentation analysis depend only on the domain architecture,
the results of sequences having the same domain architecture will be completely
identical.

The file consists of blocks separated by two newlines, and each block
corresponds to one sequence. Each block has the following format::

    AT1G61040.1
      0.0009: GO:0016570 (histone modification)
      0.0009: GO:0016569 (covalent chromatin modification)
      0.0024: GO:0016568 (chromatin modification)
      0.0036: GO:0006325 (chromatin organization)
      0.0049: GO:0051276 (chromosome organization)
      0.0055: GO:0006352 (transcription initiation)
      0.0095: GO:0006461 (protein complex assembly)
      0.0109: GO:0065003 (macromolecular complex assembly)
      0.0111: GO:0006996 (organelle organization)
      0.0126: GO:0043933 (macromolecular complex subunit organization)

In each block, the first number is the p-value obtained from the
overrepresentation analysis, the second column is the GO ID. The name
corresponding to the GO label is contained in parentheses.  Blocks containing a
sequence ID only represent sequences with no significant overrepresented GO
labels in their domain architecture.

Command line options
--------------------

.. highlight:: sh

GFam is started by the master script in ``bin/`` as follows::

    $ bin/gfam

The exact command line syntax is ``bin/gfam [options] [command]``, where
``command`` is one of the following:

``init``
    Generates a configuration file for GFam from scratch. The name of the
    configuration file will be ``gfam.cfg`` by default, but you can change
    it with the ``-c`` switch. GFam will refuse to overwrite existing
    configuration files. Example::

        $ bin/gfam -c a_lyrata.cfg init

``run``
    Runs the whole GFam pipeline. This is the default command.

``clean``
    Removes the temporary directory used to store the intermediate results.
    The name of the temporary directory is determined by the ``folder.work``
    configuration option in the `configuration file`.

    .. WARNING::
       If the output directory is the same as the temporary directory
       (``folder.work`` is equal to ``folder.output`` in the configuration),
       the ``clean`` command will also delete the final results from the output
       folder!

The default configuration file used is always ``gfam.cfg``, but it can be
overridden with the ``-c`` switch. For example, the following command will
clean the work directory specified in ``a_lyrata.cfg``::

    $ bin/gfam -c a_lyrata.cfg clean

The following extra command line switches are also available:

-h, --help                   shows a help message and then exits
-c FILE, --config-file=FILE  specifies the name of the configuration ``FILE``
-v, --verbose                enables verbose logging
-d, --debug                  shows debug messages as well
-f, --force                  forces the recalculation of the results of
                             intermediary steps in the GFam pipeline even
                             when GFam thinks everything is up-to-date.

Besides the master script, there are scripts for re-running individual steps of
the GFam pipeline. These scripts are separate Python modules in
``gfam/scripts`` and they correspond to the :ref:`steps of the GFam pipeline
<pipeline>`. It is unlikely that you will have to run them by hand, but if you
do, you have to supply the necessary input on the standard input stream of the
scripts.  For instance, if you want to do some custom filtering on a BLAST
tabular result file, you can use ``gfam/scripts/blast_filter.py`` as follows::

    $ python -m gfam.scripts.blast_filter -e 1e-5 <input.blast

This will filter ``input.blast`` and remove all entries with an E-value
larger than 10\ :sup:`-5`. The result will be written to the standard
output.

You can get a summary of the usage of each script in ``gfam/scripts``
as follows::

    $ python -m gfam.scripts.blast_filter --help

Of course replace ``blast_filter`` with the name of the script you are
interested in.  The default values of the command line switches of these
scripts come from the :ref:`configuration file <config-file>`, and they
also support ``-c`` to change the name of the configuration file.

In 99.9999% of the cases, you will only have to do ``bin/gfam init``
to create a new configuration file, ``bin/gfam`` to run the pipeline and
``bin/gfam clean`` to clean up the results.

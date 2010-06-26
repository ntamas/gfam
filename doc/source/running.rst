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
    human-readable name. It is advisable to construct a file which contains
    at least the InterPro, Pfam and Superfamily IDs as these are the most
    common (and many Pfam and Superfamily IDs do not have corresponding
    InterPro IDs yet). If you want to create such a mapping file easily,
    please refer to `Updating the mapping of IDs to human-readable names`_.

**Parent-child relationships of InterPro terms**
    This file contains the parent/child relationships between InterPro
    accession numbers to indicate family/subfamily relationships. This file
    is used to map each InterPro subfamily ID to the corresponding family
    ID, and it can be obtained from `EBI`_.

**The Gene Ontology**
    This file contains the `Gene Ontology`_ in OBO format, and it is
    required only for the `overrepresentation analysis`_ step. The latest
    version of the file can be obtained from the homepage of the
    `Gene Ontology`_ project.

.. _EBI: ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt
.. _Gene Ontology: http://www.geneontology.org
.. _overrepresentation analysis: pipeline-step-overrep

GFam accepts uncompressed files or files compressed with ``gzip`` or ``bzip2``
for both the data and the mapping files. Compressed files will be decompressed
on-the-fly in memory when needed.

Updating the mapping of IDs to human-readable names
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. highlight:: sh

As stated above, GFam relies on an external tab-separated flat file to map
domain IDs to human-readable descriptions when producing the final output.
Such a file should contain at least the InterPro, Pfam and Superfamily IDs.
The GFam distribution contains a script that can download the mappings
automatically from known sources on the Internet. The script can be
invoked as follows::

    $ bin/download_names.py >data/names.dat

This will download the InterPro, Pfam and Superfamily IDs from the Internet
and prepare the appropriate name mapping file in ``data/names.dat``. If you
wish to put it elsewhere, simply specify a different output file name. If
you omit the trailing ``>data/names.dat`` part, the mapping will be written
into the standard output. You can also compress the mapping file on-the-fly
using ``gzip`` or ``bzip2`` and use the compressed file directly in the
configuration file as GFam will uncompress it when needed. The following
command constructs a compressed name mapping file::

    $ bin/download_names.py | gzip -9 >data/names.dat.gz

Note that the script relies on the following locations to download data:

- <ftp://ftp.ebi.ac.uk/pub/databases/interpro/names.dat> for the InterPro
  name mapping

- <http://pfam.sanger.ac.uk/families?output=text> for the PFam name mapping

- <http://smart.embl-heidelberg.de/smart/descriptions.pl> for the Smart
  name mapping

- <http://scop.mrc-lmb.cam.ac.uk/scop/parse/> for the SCOP description files
  (named ``dir.des.scop.txt_X.XX``, where ``X.XX`` stands for the SCOP
  version number). It also relies on the most recent version of the SCOP
  description file being linked from the above page. The script will simply
  scan the links of the above page to determine what is the most recent
  version of SCOP. If the version number cannot be determined, the script
  will silently skip downloading the SCOP IDs.

.. _config-file:

The configuration file
----------------------

.. _output-files:

Output files
------------

Command line options
--------------------


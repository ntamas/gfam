Steps of the GFam pipeline
==========================

The GFam pipeline consists of multiple steps. In this section, we will
describe what input files does the GFam pipeline operate on, how the steps
are executed in order one by one and what output files are produced in the
end. First, a short overview of the whole process will be given, followed
by a more detailed description of each step.

Input files
-----------

The input files can be grouped into three large groups: *data files*, *mapping
files* and the *configuration file*.  Data files contain the actual input data
that is specific to a given organism. Mapping files usually map between IDs of
different data sources (for instance, from InterPro domain IDs to Gene Ontology
terms) or IDs to human-readable descriptions. The configuration file tells GFam
where to find the data files and the mapping files.  When one wants to process
a new organism with GFam, it is therefore usually enough to replace the paths
of the data files in the configuration only, as the mapping files can be
re-used for multiple analyses.

GFam requires the following data files:

**Sequence file**
    This file must contain all the sequences that are being analysed and
    annotated by GFam.

**Domain assignment file**
    This file is produced by running ``iprscan``, the command-line variant
    of `InterProScan`_ and it assigns sections of each segment in the
    sequence file to known domains in `InterPro`_. The sequence IDs in
    this file must be identical to the ones in the sequence file; if not,
    one can specify a regular expression in the :ref:`config-file` to
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
    at least the InterPro, HMMPfam and Superfamily IDs as these are the most
    common (and many HMMPfam and Superfamily IDs do not have corresponding
    InterPro IDs yet). An up-to-date assignment containing the mapping for
    these three data sources is available from the authors upon request.

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

.. _pipeline-step-extract:

Overview of a GFam analysis
---------------------------

GFam infers annotations for sequences by first finding a consensus domain
architecture for each step, then collecting Gene Ontology terms for each domain
in a given domain architecture, and performing a Gene Ontology
overrepresentation analysis on the terms to determine whether a domain
architecture is annotated by some GO terms that occur more frequently than
expected by random chance. Out of these three steps, the calculation of the
consensus domain architecture is the most complicated one, as GFam has to
account for not only the known domain assignments from InterPro, but also
for the possible existence of novel, previously uncharacterised domains.
The whole pipeline can be broken to ten steps as follows:

1. Extracting valid gene IDs from the sequence file.

2. Determining a preliminary domain architecture for each sequence by
   considering known domains from the domain assignment file only.

3. Finding the unassigned regions of each sequence; i.e. the regions
   that are not assigned to any domain in the preliminary domain
   architecture.

4. Slicing out the unassigned regions of each sequence into a separate
   FASTA file.

5. Running an all-against-all BLAST comparison of the unassigned sequence
   fragments.

6. Filtering BLAST results to determine which fragments may correspond to
   the same novel domain. Such filtering is based primarily on E-values
   and alignment lengths. At this point, we obtain a graph on the sequence
   fragments where two fragments are connected if they passed the BLAST
   filter.

7. Calculating the Jaccard similarity of the sequence fragments based on
   the connection patterns and removing those connections which have a
   low Jaccard similarity.

8. Finding the connected components of the remaining graph. Each connected
   component will correspond to a tentative novel domain.

9. Calculating the consensus domain architecture by merging the
   preliminary domain architecture with the newly detected novel domains.

10. Conducting a Gene Ontology overrepresentation analysis on each of
    the sequences and their domain architectures to derive the final
    annotations.

These steps will be described more in detail in the next few subsections.

Step 1 -- Extracting valid gene IDs
-----------------------------------


.. _pipeline-step-overrep:

Step X -- Overrepresentation analysis
-------------------------------------


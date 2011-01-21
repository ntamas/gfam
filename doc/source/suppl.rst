.. _suppl:

Supplementary scripts
=====================

.. highlight:: sh

The scripts described in this chapter are not parts of the main GFam pipeline,
but they provide useful extra functionality nevertheless. These scripts are
found in the ``bin`` subdirectory of GFam, and they can be run separately from
the command line, provided that GFam itself is on the Python path. Since the
current directory is always on the Python path, it is best to run these
scripts from the root directory of GFam.

Plotting descriptive statistics of the input file
-------------------------------------------------

The GFam pipeline has several parameters (e.g., E-value and domain length
thresholds) with sensible default values, but in order to achieve the best
results on a given dataset, these parameters can be adapted to the properties
of the input data if necessary. Such decisions are made by humans after
inspecting the distribution of E-values and domain lengths, and the distribution
of overlap sizes between different domains in the InterPro domain assignment
file. GFam can readily generate these plots using Matplotlib_, a plotting library
for Python. Matplotlib is available as a package in all major Linux distributions,
and the project provides an installer for Microsoft Windows and Mac OS X.

.. _Matplotlib: http://matplotlib.sourceforge.net

The script can be invoked as follows (assuming that the configuration file is
named ``gfam.cfg``)::

    $ bin/plot.py -c gfam.cfg figurename

where *figurename* is the name of the figure to be plotted. You may also save
figures to an output file::

    $ bin/plot.py -c gfam.cfg -o *outfile*.pdf figurename1 figurename2 ...

The supported output formats include PDF, PNG, JPG and SVG, provided that the
corresponding Matplotlib_ backends are installed. ASCII art representations
of the histograms may also be printed if the extension of the output file is
``.txt``.

To get a list of the supported figure names, specify ``list`` in place of the
figure name::

    $ bin/plot.py -c gfam.cfg list

The supported figures are as follows:

``evalue_distribution``
    Plots the count or relative frequency of domains with a given log E-value,
    sorted by different data sources in the InterPro input file, using a bin
    for each integer log E-value.

``length_distribution``
    Plots the count or relative frequency of domains with a given length,
    sorted by different data sources in the InterPro input file, using 25
    bins up to a length of at most 750. The rightmost column contains all
    domains with length greater than 750.

``overlap_distribution``
    Plots the count or relative frequency of overlap length between all pairs
    of domains that overlap by at least one residue and have the same data
    source. The plots are sorted by data sources and use a bin width of 5.

Command line options
^^^^^^^^^^^^^^^^^^^^

-a FILE, --assignment-file=FILE
                    If you don't have a GFam configuration file or you want to run
                    the script on a different InterPro assignment file (not the one
                    specified in the configuration file), you may specify the name
                    of the InterPro file directly using this switch. In this case,
                    ``-c`` is not needed.

--cumulative        Plot cumulative distributions (if that makes sense for the
                    selected plot).

-o FILE, --output=FILE
                    Specify the name of the file to save the plots to. The desired
                    format of the file is inferred from its extension. Supported
                    formats: PNG, JPG, SVG and PDF (assuming that the required
                    Matplotlib_ backends are installed). You may also use a
                    ``.txt`` extension here, which turns on `--text-mode`
                    automatically.

--relative          Plot relative frequencies instead of absolute counts on the Y
                    axis (if that makes sense for the selected plot), and use a
                    line chart instead of a bar chart.

--survival          Plot survival distributions (if that makes sense for the
                    selected plot), and use a line chart instead of a bar chart.

-t, --text-mode     Print an ASCII art representation of each histogram. This
                    option is useful if you are sitting at a non-graphical
                    terminal (e.g., an ssh shell) or if you want to dump the
                    histograms to a text file that you can analyze later. This
                    option is turned on automatically if the extension of the
                    output file is ``.txt``.

.. _updating-mappings:

Updating the mapping of IDs to human-readable names
---------------------------------------------------

GFam relies on an external tab-separated flat file to map domain IDs to
human-readable descriptions when producing the final output.  Such a file
should contain at least the InterPro, Pfam, SMART and Superfamily IDs.  The
GFam distribution contains a script that can download the mappings
automatically from known sources on the Internet. The script can be invoked as
follows::

    $ bin/download_names.py >data/names.dat

This will download the InterPro, Pfam, SMART and Superfamily IDs from the
Internet and prepare the appropriate name mapping file in ``data/names.dat``.
If you wish to put it elsewhere, simply specify a different output file name.
If you omit the trailing ``>data/names.dat`` part, the mapping will be written
into the standard output. You can also compress the mapping file on-the-fly
using ``gzip`` or ``bzip2`` and use the compressed file directly in the
configuration file as GFam will uncompress it when needed. The following
command constructs a compressed name mapping file::

    $ bin/download_names.py | gzip -9 >data/names.dat.gz

Note that the script relies on the following locations to download data:

- <ftp://ftp.ebi.ac.uk/pub/databases/interpro/names.dat> for the InterPro
  name mapping

- <http://pfam.sanger.ac.uk/families?output=text> for the Pfam name mapping

- <http://smart.embl-heidelberg.de/smart/descriptions.pl> for the SMART
  name mapping

- <http://scop.mrc-lmb.cam.ac.uk/scop/parse/> for the SCOP description files
  (named ``dir.des.scop.txt_X.XX``, where ``X.XX`` stands for the SCOP
  version number). It also relies on the most recent version of the SCOP
  description file being linked from the above page. The script will simply
  scan the links of the above page to determine what is the most recent
  version of SCOP. If the version number cannot be determined, the script
  will silently skip downloading the SCOP IDs.



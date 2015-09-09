#!/usr/bin/env python
"""Master script driving the whole of GFam

This script makes use of ``modula``, a simple Python framework that facilitates
the execution of (typically scientific) calculations where a result of a
calculation may depend on one or more other calculations. ``modula`` makes sure
that only those calculations are performed which are required to obtain the
final result; if a ``modula`` calculation is interrupted for whatever reason,
the pre-requisites that are already complete will not be re-calculated again
when one tries to resume the calculation.

``modula`` is integrated into ``gfam``, but it can be used as a standalone
Python package as well.
"""

from __future__ import with_statement

import gfam.modula as modula
import logging
import os
import shutil
import sys
import textwrap

from ConfigParser import ConfigParser
from cStringIO import StringIO
from functools import wraps
from gfam.modula.hash import sha1
from gfam.modula.module import CalculationModule
from gfam.modula.storage import DiskStorageEngine, NotFoundError
from gfam.scripts import CommandLineApp
from gfam.utils import redirected

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["GFamMasterScript"]

class GFamCalculation(CalculationModule):
    """Class representing a GFam calculation step. This is a subclass of
    `modula.CalcuationModule`_ and it assumes that the name of the module
    refers to a valid Python module in `gfam.scripts`_."""

    def run(self):
        """Runs the calculation"""
        self.logger.info("Starting module %s" % self.name)

        self.prepare()

        # Search for the CommandLineApp object in the module
        app = []
        for value in self.module.__dict__.itervalues():
            if isinstance(value, type) and value != CommandLineApp \
                    and issubclass(value, CommandLineApp):
                app.append(value)

        if len(app) != 1:
            raise ValueError("more than one CommandLineApp in %s" % self.name)

        # Create the application
        app = app[0](logger=self.logger)
        args = ["-c", self.config.get("@global.config_file")]

        for param, value in self.parameters.iteritems():
            if not param.startswith("switch."):
                continue
            switch, value = value.split(" ", 1)
            value = modula.storage_engine.get_filename(value.strip())
            args.extend([switch, value])

        if "infile" in self.parameters:
            infiles = self.parameters["infile"].split(",")
            for infile in infiles:
                infile = modula.storage_engine.get_filename(infile.strip())
                args.append(infile)

        if "stdin" in self.parameters:
            stdin = modula.storage_engine.get_source(self.parameters["stdin"])
        else:
            stdin = None

        out_fname = modula.storage_engine.get_filename(self.name)
        stdout = modula.storage_engine.get_result_stream(self, mode="wb")
        try:
            with redirected(stdin=stdin, stdout=stdout):
                retcode = app.run(args)
            stdout.close()
            if retcode:
                raise RuntimeError("non-zero return code from child module")
        except:
            # If an error happens, remove the output file and re-raise
            # the exception
            stdout.close()
            os.unlink(out_fname)
            raise

        self.logger.info("Finished module %s" % self.name)


class GFamDiskStorageEngine(DiskStorageEngine):
    """Disk storage engine for GFam that has an empty `store` method.  This is
    because `GFamCalculation` writes directly into the results file and always
    returns ``None``, so there is no need to store the results explicitly.
    """

    def get_filename(self, source_or_module):
        """Retrieves the filename corresponding to a data source."""
        try:
            module = self.module_manager.get(source_or_module)
            if hasattr(module, "filename"):
                return module.filename
            return self._get_module_result_filename(module)
        except Exception, ex:
            raise NotFoundError(source_or_module, source_or_module, ex)


    def store(self, module, result):
        """Empty, does nothing"""
        pass


def needs_config(func):
    """Decorator for methods in `GFamMasterScript` that require a valid
    configuration file.
    """
    @wraps(func)
    def wrapper(self, *args, **kwds):
        if self.config is None:
            self.error("Configuration file not found: %s. Try gfam init "
                       "to generate a new one." % self.options.config_file)
        return func(self, *args, **kwds)
    return wrapper


class GFamMasterScript(CommandLineApp):
    """\
    Usage: %prog [options] [command]

    Runs the whole GFam pipeline, driven by the given configuration
    file (specified with the `-c` option). The given command must be
    one of the following:

        - init: generates a configuration file from scratch.

        - run: runs the whole pipeline. This is the default.

        - clean: removes the temporary directory used to store
          intermediate results.
    """

    short_name = "gfam"

    def __init__(self, *args, **kwds):
        super(GFamMasterScript, self).__init__(*args, **kwds)
        self.modula = None
        self.config = None

    def create_parser(self):
        """Creates the command line parser for the GFam master script"""
        parser = super(GFamMasterScript, self).create_parser()
        parser.add_option("-f", "--force", dest="force", action="store_true",
                help="force recalculation of results even when gfam thinks "\
                     "everything is up-to-date")
        return parser

    def get_modula_config(self, config):
        """Based on the given `ConfigParser` instance in `config`, constructs
        another `ConfigParser` that tells Modula which tasks to execute and
        where each of the input files are to be found."""
        modula_config_str = textwrap.dedent("""\
        [@global]
        [@paths]
        [@inputs]

        [extract_gene_ids]
        depends=file.input.sequences
        infile=file.input.sequences

        [assignment_source_filter]
        depends: file.input.iprscan, file.input.sequences,
            file.mapping.interpro_parent_child, extract_gene_ids
        infile=file.input.iprscan
        switch.0=-g extract_gene_ids

        [find_unassigned]
        depends=assignment_source_filter, file.input.sequences
        infile=assignment_source_filter

        [seqslicer]
        depends=find_unassigned, file.input.sequences
        infile=find_unassigned, file.input.sequences

        [blast_all]
        depends=seqslicer
        infile=seqslicer
        switch.0=-o blast_all

        [blast_filter]
        depends=blast_all, seqslicer
        infile=blast_all
        switch.0=-S seqslicer

        [jaccard]
        depends=blast_filter
        infile=blast_filter

        [cca]
        depends=jaccard
        infile=jaccard

        [find_domain_arch]
        depends=assignment_source_filter, cca
        infile=assignment_source_filter, cca

        [label_assignment]
        depends=file.mapping.gene_ontology, file.mapping.interpro2go, find_domain_arch
        infile=file.mapping.gene_ontology, file.mapping.interpro2go, find_domain_arch

        [overrep]
        depends=file.mapping.gene_ontology, file.mapping.interpro2go, find_domain_arch
        infile=file.mapping.gene_ontology, file.mapping.interpro2go, find_domain_arch

        """)

        modula_config = ConfigParser()
        modula_config.readfp(StringIO(modula_config_str))

        # Store the name of the config file
        modula_config.set("@global", "config_file", self.options.config_file)

        # Store the hash of the configuration as a default parameter for
        # all the algorithms
        config_str = StringIO()
        config.write(config_str)
        modula_config.set("DEFAULT", "config_file_hash", \
                sha1(config_str.getvalue()).hexdigest())

        # Set up the module and storage path
        modula_config.set("@paths", "modules", \
            os.path.dirname(sys.modules[__name__].__file__))
        modula_config.set("@paths", "storage", \
            config.get("DEFAULT", "folder.work"))

        # Add the input files
        for name, value in config.items("generated"):
            if name.startswith("file."):
                modula_config.set("@inputs", name, value)

        return modula_config

    def read_config(self):
        """Reads the configuration from the given file and returns an
        appropriate `ConfigParser` instance."""
        self.options.config_file = self.options.config_file or "gfam.cfg"

        config_file = self.options.config_file
        if not os.path.exists(config_file):
            return None

        config = ConfigParser()
        config.read([config_file])
        return config

    def run_real(self):
        """Runs the application"""

        # Shut up the root logger, logging will be done by modula
        logging.getLogger('').handlers = []

        # Read the configuration file
        self.config = self.read_config()

        if self.config:
            # Initialize Modula
            modula_config = self.get_modula_config(self.config)
            modula.init(modula_config, debug = self.options.debug,
                        storage_engine_factory = GFamDiskStorageEngine)
            modula.module_manager.module_factory = GFamCalculation
            self.modula = modula

            # Use the logger from Modula
            self.log = self.modula.logger

        # Run the commands
        if not self.args:
            self.args = ["run"]
        for command in self.args:
            try:
                method = getattr(self, "do_%s" % command)
            except AttributeError:
                self.log.fatal("No such command: %s" % command)
                return 1
            error_code = method()
            if error_code:
                return error_code

    @needs_config
    def do_clean(self):
        """Clears the temporary directory used for intermediate results."""
        # Get the output folder name
        outfolder = self.config.get("DEFAULT", "folder.work")

        # Remove the folder
        self.log.info("Removing temporary folder: %s" % outfolder)
        shutil.rmtree(outfolder)
        self.log.info("Temporary folder removed successfully.")

    def do_init(self):
        """Initializes a configuration file in the current directory."""
        # Check whether we already have a config file in the current directory.
        cfgfile = self.options.config_file

        if os.path.exists(cfgfile):
            self.error("%s already exists, please remove it before "
                       "generating a new one." % cfgfile)
            return 1

        with open(cfgfile, "w") as fp:
            fp.write(CONFIGURATION_FILE_TEMPLATE)

        print "Configuration file generated successfully."
        print
        print "Please open the configuration file in the text editor and provide"
        print "the correct values for the following configuration keys:"
        print

        # Check which config keys are empty in the config file, these have
        # to be filled
        config = ConfigParser()
        config.read([cfgfile])
        seen_keys = set()
        for section in ["DEFAULT"] + config.sections():
            for name, value in config.items(section, raw=True):
                if not value and name not in seen_keys:
                    print "  - %s (in section %r)" % (name, section)
                    if section == "DEFAULT":
                        seen_keys.add(name)

    @needs_config
    def do_run(self):
        """Runs the whole GFam pipeline"""
        # Get the output folder name
        outfolder = os.path.abspath(self.config.get("DEFAULT", "folder.output"))

        # Create the output folder if needed
        if outfolder and not os.path.isdir(outfolder):
            self.modula.logger.warning("Creating output folder: %s" % outfolder)
            os.makedirs(outfolder)

        # Run and export the inferred domain architectures
        outfile = os.path.join(outfolder, "domain_architectures.tab")
        self.modula.run("find_domain_arch", force=self.options.force)
        shutil.copy(self.modula.storage_engine.get_filename("find_domain_arch"),
                outfile)
        self.log.info("Exported domain architectures to %s." % outfile)

        # Run and export the label assignment
        outfile = os.path.join(outfolder, "assigned_labels.txt")
        self.modula.run("label_assignment", force=self.options.force)
        shutil.copy(self.modula.storage_engine.get_filename("label_assignment"),
                outfile)
        self.log.info("Exported label assignment to %s." % outfile)

        # Run and export the overrepresentation analysis
        outfile = os.path.join(outfolder, "overrepresentation_analysis.txt")
        self.modula.run("overrep", force=self.options.force)
        shutil.copy(self.modula.storage_engine.get_filename("overrep"), outfile)
        self.log.info("Exported overrepresentation analysis to %s." % outfile)

###########################################################################

CONFIGURATION_FILE_TEMPLATE = """\
###########################################################################
## Configuration file for GFam                                           ##
###########################################################################

[DEFAULT]

###########################################################################
## Input files                                                           ##
###########################################################################

# Raw output file from IPRScan that contains domain assignments for all the
# sequences we are interested in
file.input.iprscan=

# A FASTA file containing all the sequences being analysed
file.input.sequences=

# Log file to store why given sequences have been rejected during the
# filtering of the input IPRScan file. Feel free to leave it empty.
file.log.iprscan_exclusions=

# A file containing the Gene Ontology in OBO format
file.mapping.gene_ontology=data/gene_ontology.obo

# File containing the mapping of GO terms to InterPro entries, as downloaded
# from geneontology.org
file.mapping.interpro2go=data/interpro2go

# File containing a tab-separated list of InterPro IDs and their corresponding
# human-readable descriptions. This can be constructed by the following
# command::
#
#     $ bin/download_names.py | gzip -9 >data/names.dat.gz
file.mapping.interpro2name=data/names.dat.gz

# File containing the parent-child relationships of InterPro terms
file.mapping.interpro_parent_child=data/ParentChildTreeFile.txt

###########################################################################
## Folders                                                               ##
###########################################################################

# The working folder in which to put intermediary files
folder.work=work

# The output folder in which to put the final results
folder.output=work

###########################################################################
## External utilities used by GFam                                       ##
###########################################################################

[utilities]

# The folder containing the BLAST executables. It does not matter whether you
# have the old C-based or the newer C++-based tools, GFam can use both if you
# also have the ``legacy_blast.pl`` script that adapts the new tools to the
# command line syntax used by the older ones.
folder.blast=

# The path to formatdb. You may use the name of the folder containing the BLAST
# executables or the full path (including the name of the tool). If you have the
# newer, C++-based BLAST tools (which do not have formatdb), pass the name of the
# folder containing the BLAST executables here, and if you have the
# ``legacy_blast.pl`` script in the same folder (plus a working Perl setup), GFam
# will detect the situation and run ``legacy_blast.pl`` accordingly.
util.formatdb=%(folder.blast)s

# The path to blastall. You may use the name of the folder containing the BLAST
# executables or the full path (including the name of the tool). If you have the
# newer, C++-based BLAST tools (which do not have blastall), pass the name of the
# folder containing the BLAST executables here, and if you have the
# ``legacy_blast.pl`` script in the same folder (plus a working Perl setup), GFam
# will detect the situation and run ``legacy_blast.pl`` accordingly.
util.blastall=%(folder.blast)s

###########################################################################
## Analysis parameters                                                   ##
###########################################################################

# Some of the parameters in the next few section refer to other configuration
# options in the form of %(item)s. In general, you should keep these
# references and modify the referred keys instead.

[DEFAULT]

# Python regular expression that matches gene IDs from the sequence file.
# This is necessary if the IPRScan output uses a sequence ID that is
# only a part of the sequence ID in the input FASTA file. If it is
# empty, no ID transformation will be done, the sequence IDs in the
# input FASTA file will be matched intact to the IPRScan output.
# If it is not empty, it must be a valid Python regular expression
# with a *named* group "id" that matches the gene ID that is used
# in the IPRScan output. If you don't know what named groups are,
# check the documentation of the Python ``re`` module.
sequence_id_regexp=

# Which assignment sources NOT to trust from InterPro? (space separated)
untrusted_sources=HAMAP PatternScan FPrintScan Seg Coil

# Maximum overlap allowed between assignments originating from the
# same data source
max_overlap=20

# Hint on the number of CPU cores to use during the analysis. Currently
# the BLAST invocation uses this hint to select the number of threads
# used by BLAST to speed up calculations.
#
# The default value is 1 since it is not possible to auto-detect the
# number of CPU cores in a platform independent way. Feel free to raise this
# value if your computer has multiple CPU cores.
num_cpu_cores=1

[analysis:iprscan_filter]

# E-value thresholds to use when processing the initial InterPro file.
# This entry is a semicolon-separated list of source=threshold pairs, the given
# threshold will be used for the given data source. If an entry contains
# only a threshold, this will be considered as a default threshold for
# sources not explicitly mentioned here.
e_value_thresholds=1e-3;superfamily=inf;HMMPanther=inf;Gene3D=inf;HMMPIR=inf

# File containing the parent-child relationships of InterPro terms
interpro_parent_child_mapping=%(file.mapping.interpro_parent_child)s

# These configuration keys specify which assignment sources are to be taken
# into account at each stage of the analysis. For more information about what
# these stages are, please refer to the documentation, especially the
# description of Step 2 in section "Steps of the GFam pipeline"
stages.1=ALL-HMMPanther-Gene3D
stages.2=ALL-HMMPanther-Gene3D
stages.3=ALL

[analysis:find_unassigned]

# Minimum number of amino acids in a sequence in order to consider it further
# (i.e. to calculate its unassigned fragments)
min_seq_length=30

# Minimum number of amino acids in a sequence fragment in order to
# consider that fragment as a novel domain candidate
min_fragment_length=75

# Input FASTA file containing all the sequences of the representative
# gene model being analysed
sequences_file=%(file.input.sequences)s

[analysis:blast_filter]

# Minimum sequence identity between two fragments to consider them
# as being in the same domain
min_seq_identity=45

# Minimum alignment length between two fragments to consider them
# as being in the same domain. If normalization_method is not off,
# this must be the normalized alignment length threshold according
# to the chosen normalization method.
min_alignment_length=0.7

# Maximum E-value between two fragments in order to consider them
# as being in the same domain
max_e_value=1e-3

# Normalization method to use for calculating normalized alignment length
# Must be one of: off, smaller, larger, query, hit
normalization_method=query

[analysis:jaccard]

# Minimum Jaccard similarity between the neighbour sets of two
# fragments in order to consider them as being in the same domain
min_similarity=0.66

# Whether to assume that a protein is connected to itself or not
assume_loops=1

# Whether to consider only those protein pairs which are linked in
# the input file. If this is 1, protein pairs not in the input file
# will not be returned even if their Jaccard similarity is larger
# than the given threshold
only_linked=1

[analysis:cca]

# Empty section, this script has no individual parameters to tune, but
# this might change in the future

[analysis:find_domain_arch]

# A novel domain occur in at least this number of sequences
min_novel_domain_size=4

[analysis:label_assignment]

# Empty section, this script has no individual parameters to tune, but
# this might change in the future

[analysis:overrep]

# The p-value threshold in the hypergeometric test used in the
# overrepresentation analysis process
confidence=0.05

# The method used to account for multiple hypothesis testing.
# Valid choices: bonferroni (Bonferroni correction), Sidak
# (Sidak correction), fdr (Benjamini-Hochberg method), none (off).
# The Bonferroni and Sidak correction methods control the
# family-wise error rate (FWER), while the Benjamini-Hochberg method
# control the false discovery rate.
correction=fdr

# The minimum number of annotated domains a GO term must have in order
# to be considered in the overrepresentation analysis
min_term_size=1

[analysis:coverage]

# Empty section, this script has no individual parameters to tune, but
# this might change in the future

## WARNING: you should not have to modify anything beyond this point.
## Make sure that you know what you are doing.

###########################################################################

###########################################################################
## Paths for the generated files                                         ##
###########################################################################

[generated]

# File in which the unassigned sequence fragments are stored
file.unassigned_fragments=%(folder.work)s/unassigned_fragments.ffa

# File containing a list of valid gene IDs (extracted from the input file)
file.valid_gene_ids=%(folder.work)s/gene_ids.txt

# File containing the detailed final domain architecture for each sequence
file.domain_architecture_details=%(folder.output)s/domain_architecture_details.txt

# File containing genome-level domain architecture statistics
file.domain_architecture_stats=%(folder.output)s/domain_architecture_stats.txt
"""

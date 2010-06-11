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
from gfam.modula.hash import sha1
from gfam.modula.module import CalculationModule
from gfam.modula.storage import DiskStorageEngine, NotFoundError
from gfam.scripts import CommandLineApp
from gfam.utils import redirected

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"


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
                app.run(args)
            stdout.close()
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


class GFamMasterScript(CommandLineApp):
    """\
    Usage: %prog [options]

    Runs the whole GFam pipeline, driven by the given configuration
    file (specified with the `-c` option).
    """

    short_name = "gfam"

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

        [overrep]
        depends=file.input.gene_ontology, file.mapping.interpro2go, find_domain_arch
        infile=file.input.gene_ontology, file.mapping.interpro2go, find_domain_arch

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
            self.error("Configuration file not found: %s" % config_file)

        config = ConfigParser()
        config.read([config_file])
        return config

    def run_real(self):
        """Runs the application"""

        # Shut up the root logger, logging will be done by modula
        logging.getLogger('').handlers = []

        # Read the configuration file
        config = self.read_config()

        # Initialize Modula
        modula_config = self.get_modula_config(config)
        modula.init(modula_config, debug = self.options.debug,
                    storage_engine_factory = GFamDiskStorageEngine)
        modula.module_manager.module_factory = GFamCalculation

        # Get the output folder name
        outfolder = config.get("DEFAULT", "folder.work")

        # Run and export the inferred domain architectures
        outfile = os.path.join(outfolder, "domain_architectures.txt")
        modula.run("find_domain_arch", force=self.options.force)
        shutil.copy(modula.storage_engine.get_filename("find_domain_arch"),
                outfile)
        self.log.info("Exported domain architectures to %s." % outfile)

        # Run and export the overrepresentation analysis
        outfile = os.path.join(outfolder, "overrepresentation_analysis.txt") 
        modula.run("overrep", force=self.options.force)
        shutil.copy(modula.storage_engine.get_filename("overrep"), outfile)
        self.log.info("Exported overrepresentation analysis to %s." % outfile)


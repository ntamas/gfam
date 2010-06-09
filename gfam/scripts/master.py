#!/usr/bin/env python
"""Master script driving the whole of GFam"""

import gfam.modula as modula
import logging
import os
import sys
import textwrap

from ConfigParser import ConfigParser
from cStringIO import StringIO
from gfam.scripts import CommandLineApp
from gfam.utils import redirected

from gfam.modula.configuration import Configuration
from gfam.modula.module import CalculationModule, DefaultModuleManager
from gfam.modula.storage import DiskStorageEngine

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"


class GFamCalculation(CalculationModule):
    """Class representing a GFam calculation step."""

    def __init__(self, *args, **kwds):
        super(GFamCalculation, self).__init__(*args, **kwds)
        self._verbose = False

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, verbose):
        self._verbose = bool(verbose)

    def run(self):
        self.logger.info("Starting module %s" % self.name)

        self.prepare()

        # Search for the CommandLineApp object in the module
        app = []
        for key, value in self.module.__dict__.iteritems():
            if isinstance(value, type) and value != CommandLineApp \
                    and issubclass(value, CommandLineApp):
                app.append(value)

        if len(app) != 1:
            raise ValueError("more than one CommandLineApp in %s" % self.name)

        # Create the application
        app = app[0](logger=self.logger)
        args = ["-c", self.config.get("@global.config_file")]

        if self.verbose:
            args.append("-v")

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
    """Disk storage engine for GFam that has an empty `store` method.
    This is because `GFamCalculation` writes directly into the results
    file.
    """

    def get_filename(self, source_or_module):
        """Retrieves the filename corresponding to a data source."""
        try:
            module = self.module_manager.get(source_or_module)
            if hasattr(module, "filename"):
                return module.filename
            return self._get_module_result_filename(module)
        except Exception, ex:
            raise NotFoundError(source_name, source_name, ex)


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
        parser = super(GFamMasterScript, self).create_parser()
        return parser

    def run_real(self):
        """Runs the application"""

        # Shut up the root logger, logging will be done by modula
        logging.getLogger('').handlers = []

        config_file = self.options.config_file or "gfam.cfg"

        if not os.path.exists(config_file):
            self.error("Configuration file not found: %s" % config_file)

        config = ConfigParser()
        config.read([config_file])

        modula_config_str = textwrap.dedent("""\
        [@global]
        config_file=%s
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

        """ % config_file)
        modula_configuration = ConfigParser()
        modula_configuration.readfp(StringIO(modula_config_str))

        modula_configuration.set("@paths", "modules", \
            os.path.dirname(sys.modules[__name__].__file__))
        modula_configuration.set("@paths", "storage", \
            config.get("DEFAULT", "folder.work"))

        for name, value in config.items("generated"):
            if name.startswith("file."):
                modula_configuration.set("@inputs", name, value)

        modula.init(modula_configuration, \
                    storage_engine_factory = GFamDiskStorageEngine,
                    debug = self.options.debug)
        modula.module_manager.module_factory = GFamCalculation

        modula.run("overrep")

if __name__ == "__main__":
    sys.exit(GFamMasterScript().run())


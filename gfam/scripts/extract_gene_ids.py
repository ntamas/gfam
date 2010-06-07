#!/usr/bin/env python

import sys

from gfam import fasta
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

class ExtractGeneIDsApp(CommandLineApp):
    """\
    Usage: %prog [options] [result_file]

    Extracts the gene IDs from a FASTA file.
    """

    def create_parser(self):
        """Creates the command line parser for this application"""
        parser = super(ExtractGeneIDsApp, self).create_parser()
        parser.add_option("-r", "--seq-id-regexp", metavar="REGEXP",
                help="remap sequence IDs using REGEXP",
                config_key="sequence_id_regexp",
                dest="sequence_id_regexp")
        return parser

    def run_real(self):
        """Runs the application"""
        if not self.args:
            infiles = ["-"]
        else:
            infiles = self.args

        for infile in infiles:
            self.process_file(infile)

    def process_file(self, filename):
        """Processes the given input file"""
        self.log.info("Processing %s..." % filename)

        parser = fasta.Parser(open_anything(filename))
        parser = fasta.regexp_remapper(parser,
                self.options.sequence_id_regexp)
        for seq in parser:
            print seq.id


if __name__ == "__main__":
    sys.exit(ExtractGeneIDsApp().run())


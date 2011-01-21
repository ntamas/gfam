#!/usr/bin/env python

import optparse
import sys

from collections import defaultdict
from gfam.interpro import AssignmentReader, InterPro
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

class InterProCoverageApp(CommandLineApp):
    """\
    Usage: %prog [options] [assignment_file]

    For each data source in the given assignment file, calculates the number
    of unique domain IDs for each data source and the fraction of those IDs
    that are also integrated into InterPro. For the latter, the script
    needs the ParentChildTreeFile.txt from InterPro, whose path may be
    passed either from a config file or from the command line.
    """

    short_name = "interpro_coverage"

    def create_parser(self):
        parser = super(InterProCoverageApp, self).create_parser()
        parser.add_option("-i", "--interpro-parent-child-file",
                dest="parent_child_file",
                metavar="FILE",
                help="use the given InterPro parent-child FILE",
                config_key="file.mapping.interpro_parent_child",
                default=None)
        return parser

    def run_real(self):
        """Runs the application"""
        if not self.options.parent_child_file:
            self.parser.error("must specify the InterPro parent-child "
                              "file using -i")

        self.interpro = InterPro.FromFile(self.options.parent_child_file)

        if not self.args:
            self.args = ["-"]

        for arg in self.args:
            self.process_infile(arg)

    def process_infile(self, fname):
        """Processes the given input file `fname`, which must be either a filename
        or a stream. If the filename is ``-``, it is assumed to be the standard
        input."""
        mapper = self.interpro.mapping
        self.log.info("Processing %s..." % fname)

        domain_ids = defaultdict(set)
        domain_ids_without_interpro_ids = defaultdict(set)

        for assignment in AssignmentReader(fname):
            if not assignment.source or not assignment.domain:
                continue
            if ":SF" in assignment.domain:
                domain = assignment.domain[0:assignment.domain.index(":SF")]
            else:
                domain = assignment.domain
            domain_ids[assignment.source].add(domain)

        print "Source\tInInterPro\tNotInInterPro\tTotal\tPercentage"
        for source in sorted(domain_ids):
            ids = domain_ids[source]
            total = len(ids)
            unknown = sum(1 for item in ids if item not in mapper)
            percentage = 100.0 * unknown / total
            print "%s\t%d\t%d\t%d\t%.2f%%" % (source, total-unknown, unknown, total, percentage)


if __name__ == "__main__":
    sys.exit(InterProCoverageApp().run())


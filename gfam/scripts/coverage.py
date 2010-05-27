#!/usr/bin/env python

import optparse
import sys

from gfam.interpro import AssignmentParser
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything, Sequence, UniversalSet

class CoverageApp(CommandLineApp):
    """\
    Usage: %prog [options] [file]

    Calculates the overall coverage of each sequence in the given
    InterPro assignment file"""

    short_name = "coverage"

    def create_parser(self):
        parser = super(CoverageApp, self).create_parser()
        parser.add_option("-x", "--exclude", dest="ignored",
                metavar="SOURCE",
                help="ignore the given assignment source",
                action="append", default=[])
        parser.add_option("-g", "--gene-ids", dest="gene_id_file",
                metavar="FILE", help="only consider those gene IDs which "+
                   "are present in the list in the given FILE",
                default=None)
        return parser

    def run_real(self):
        """Runs the application"""
        if self.options.gene_id_file:
            self.log.info("Loading sequence IDs from %s..." % \
                          self.options.gene_id_file)
            self.valid_sequence_ids = set()
            for line in open_anything(self.options.gene_id_file):
                self.valid_sequence_ids.add(line.strip())
        else:
            self.valid_sequence_ids = UniversalSet()

        self.ignored = set()
        for source in self.options.ignored:
            self.ignored.update(source.split())

        if self.ignored:
            self.log.info("Ignoring sources: %s" % ", ".join(self.ignored))

        if not self.args:
            self.args = ["-"]

        for arg in self.args:
            self.process_infile(arg)

    def process_infile(self, fname):
        """Processes the given input file `fname`, which must be either a filename
        or a stream. If the filename is ``-``, it is assumed to be the standard
        input."""
        self.log.info("Processing %s..." % fname)
        current_id, assignments = None, []
        infile = open_anything(fname)
        parser = AssignmentParser()
        valid_ids = self.valid_sequence_ids
        for line in infile:
            line = line.strip()
            assignment = parser.parse(line)
            if assignment.source in self.ignored:
                continue
            if assignment.id not in valid_ids:
                continue
            if assignment.id != current_id:
                self.process_sequence(current_id, assignments)
                current_id = assignment.id
                assignments = list()
            assignments.append(assignment)
        self.process_sequence(current_id, assignments)

    def process_sequence(self, name, assignments):
        """Processes the given sequence `name` with the given `assignments`."""
        if not assignments:
            return

        seq = Sequence(name, assignments[0].length)
        for assignment in assignments:
            seq.assign(assignment, overlap_check=False)

        print "%s %.4f" % (name, seq.coverage())


if __name__ == "__main__":
    sys.exit(CoverageApp().run())

#!/usr/bin/env python

import optparse
import sys

from collections import defaultdict
from gfam.assignment import SequenceWithAssignments
from gfam.interpro import AssignmentReader
from gfam.scripts import CommandLineApp
from gfam.utils import complementerset, open_anything

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"


class SequenceLevelOutputFormatter(object):
    """Output formatter class that prints sequence-level statistics for each
    sequence in the input data file.

    For each sequence and each data source used in the domain assignments of
    that sequence, the sequence-level output formatter prints the following
    statistics in a row, separated by tabs:

        - the ID of the sequence
        - the length of the sequence
        - the data source
        - the number of residues covered by at least one domain from the
          data source
        - the fraction of residues covered by at least one domain from the
          data source

    A special data source named ``ALL`` represents all the data sources
    combined (excluding those that were not read from the input file according
    to the command line options).
    """

    def process_assignments(self, seq):
        for source in seq.data_sources():
            cov = seq.coverage(sources=[source])
            print "%s\t%d\t%s\t%d\t%.4f" % (name, len(seq), source, round(len(seq)*cov), cov)
        cov = seq.coverage()
        print "%s\t%d\tALL\t%d\t%.4f" % (name, len(seq), round(len(seq)*cov), cov)

    def finish(self):
        pass


class GenomeLevelOutputFormatter(object):
    """Output formatter class that prints genome-level statistics for the
    input data file.
    
    The output contains one line for each data source in the input data file.
    The rows contain the following information, separated by tabs:
        
        - The name of the data source
        - The number of sequences annotated by the data source
        - The number of distinct domain architectures obtained from the
          data source
        - The fraction of residues covered by at least one domain from the
          data source
          
    ALL as a data source name denotes the union of all data sources."""

    def __init__(self):
        self.families_by_source = defaultdict(set)
        self.num_sequences_by_source = defaultdict(int)
        self.total_covered_by_source = defaultdict(int)
        self.total_residues = 0

    def process_assignments(self, seq):
        self.total_residues += len(seq)

        for source in seq.data_sources():
            family = tuple(seq.domain_architecture(sources=[source]))
            self.num_sequences_by_source[source] += 1
            self.total_covered_by_source[source] += round(len(seq) * seq.coverage(sources=[source]))
            if family:
                self.families_by_source[source].add(family)

        family = tuple(seq.domain_architecture())

        self.num_sequences_by_source["ALL"] += 1
        self.total_covered_by_source["ALL"] += round(len(seq) * seq.coverage())
        if family:
            self.families_by_source["ALL"].add(family)

    def finish(self):
        sources = set(self.num_sequences_by_source.iterkeys())
        for source in sorted(sources):
            num_seqs = self.num_sequences_by_source[source]
            num_families = len(self.families_by_source[source])
            print "%s\t%d\t%d\t%.4f" % (source, num_seqs, num_families,
                    self.total_covered_by_source[source] / self.total_residues)


class CoverageApp(CommandLineApp):
    """\
    Usage: %prog [options] [file]

    Calculates some coverage statistics in the given InterPro assignment file.
    The output includes one line per sequence with the following information
    in tab-separated columns: sequence ID, sequence length, coverage.
    """

    short_name = "coverage"

    def create_parser(self):
        parser = super(CoverageApp, self).create_parser()
        parser.add_option("-i", "--include", dest="include_sources",
                metavar="SOURCE", help="use the given assignment source",
                config_key="analysis:coverage/include_sources",
                action="append", default=[])
        parser.add_option("-x", "--exclude", dest="exclude_sources",
                metavar="SOURCE",
                help="ignore the given assignment source",
                config_key="analysis:coverage/untrusted_sources",
                action="append", default=[])
        parser.add_option("-g", "--gene-ids", dest="gene_id_file",
                metavar="FILE", help="only consider those gene IDs which "+
                   "are present in the list in the given FILE",
                config_key="generated/file.valid_gene_ids",
                default=None)
        parser.add_option("--totals", dest="print_totals", action="store_true",
                help="print genome-level statistics instead of sequence-level "
                     "statistics",
                default=False,
                config_key="analysis:coverage/print_totals")
        return parser

    def run_real(self):
        """Runs the application"""

        # Load valid sequence IDs (if necessary)
        if self.options.gene_id_file:
            self.log.info("Loading sequence IDs from %s..." % \
                          self.options.gene_id_file)
            self.valid_sequence_ids = set()
            for line in open_anything(self.options.gene_id_file):
                self.valid_sequence_ids.add(line.strip())
        else:
            self.valid_sequence_ids = complementerset()

        # Find which sources will be allowed
        if not self.options.include_sources:
            self.sources = complementerset()
        else:
            self.sources = set(self.options.include_sources)
        self.sources.difference_update(self.options.exclude_sources)
        if isinstance(self.sources, complementerset):
            self.log.info("Ignored sources: %s" % ", ".join(self.sources.iterexcluded()))
        else:
            self.log.info("Accepted sources: %s" % ", ".join(self.sources))

        if not self.args:
            self.args = ["-"]

        for arg in self.args:
            # Set up the output formatter
            if self.options.print_totals:
                self.output_formatter = GenomeLevelOutputFormatter()
            else:
                self.output_formatter = SequenceLevelOutputFormatter()
            # Process the file
            self.process_infile(arg)
            # Print the results
            self.output_formatter.finish()

    def process_infile(self, fname):
        """Processes the given input file `fname`, which must be either a filename
        or a stream. If the filename is ``-``, it is assumed to be the standard
        input."""
        self.log.info("Processing %s..." % fname)
        current_id, assignments = None, []
        valid_ids = self.valid_sequence_ids
        for assignment in AssignmentReader(fname):
            if assignment.source not in self.sources:
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

        seq = SequenceWithAssignments(name, assignments[0].length)
        for assignment in assignments:
            seq.assign(assignment, overlap_check=False)
        self.output_formatter.process_assignments(seq)


if __name__ == "__main__":
    sys.exit(CoverageApp().run())

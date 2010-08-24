#!/usr/bin/env python

import sys

from gfam.blast import BlastFilter
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

class BlastFilterApp(CommandLineApp):
    """\
    Usage: %prog [options] [result_file]

    Filters the given BLAST result file (which must be in tabular
    format) and drops matches not satisfying the criteria given
    in the options.
    """

    def create_parser(self):
        """Creates the parser that parses the command line options"""
        parser = super(BlastFilterApp, self).create_parser()
        parser.add_option("-s", "--sequence-identity",
                dest="sequence_identity", metavar="PERCENTAGE",
                type=float, default=0,
                config_key="analysis:blast_filter/min_seq_identity",
                help="drop matches with sequence identity smaller "
                     "than the given PERCENTAGE")
        parser.add_option("-l", "--alignment-length",
                dest="alignment_length", metavar="LENGTH",
                type=float, default=0,
                config_key="analysis:blast_filter/min_alignment_length",
                help="drop matches with alignment length less "
                     "than the given LENGTH")
        parser.add_option("-n", "--normalize",
                dest="normalize_alignment_length",
                metavar="METHOD", type="choice",
                choices=["off", "query", "hit", "smaller", "larger"], 
                default="off",
                config_key="analysis:blast_filter/normalization_method",
                help="normalize the alignment length to between 0 and 1 "
                     "by dividing it with the length of the sequence given "
                     "by METHOD (smaller, larger, query or hit)")
        parser.add_option("-e", "--max-e-value", dest="max_e_value",
                metavar="VALUE", type=float, default=100,
                config_key="analysis:blast_filter/max_e_value",
                help="drop matches with an E-value larger than the given VALUE")
        parser.add_option("-S", "--sequences",
                dest="sequences_file", metavar="FILE",
                help="FASTA file containing all the sequences used by BLAST. "
                     "This is necessary if -n is used.",
                config_key="generated/file.unassigned_fragments",
                default=None)
        return parser

    def run_real(self):
        """Runs the BLAST filtering application"""
        if not self.args:
            infiles = [sys.stdin]
        else:
            infiles = self.args

        filter = self.construct_blast_filter()
        for infile in infiles:
            self.process_file(infile, filter)


    def construct_blast_filter(self):
        """Constructs an appropriate `BlastFilter` from the command line
        options and returns it."""
        options = self.options

        filter = BlastFilter()
        filter.set_normalize_func(options.normalize_alignment_length)
        filter.max_e_value = options.max_e_value
        filter.min_sequence_identity = options.sequence_identity
        filter.min_alignment_length = options.alignment_length

        if options.normalize_alignment_length != "off":
            if not options.sequences_file:
                self.error("must specify sequences file using "\
                           "-S when -n is given")
            filter.load_sequences_from_file(options.sequences_file)

        return filter

    def process_file(self, filename, filter):
        """Processes the given file using the given `filter`."""
        self.log.info("Processing %s..." % filename)
        for line in self.process_lines(open_anything(filename), filter):
            sys.stdout.write(line)

    def process_lines(self, lines, filter):
        """Processes the lines yielded by the given generator.
        The input generator must yield lines in BLAST's tabular
        format. The result of this method is another generator that
        yields only those lines that pass the given `filter` (an
        instance of `BlastFilter`). Comments and empty lines are kept.
        """
        return (line for line in lines if filter.accepts(line))


if __name__ == "__main__":
    sys.exit(BlastFilterApp().run())

#!/usr/bin/env python
"""Sequence slicer application"""

import optparse
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything

class SeqSlicerApp(CommandLineApp):
    """\
    Usage: %prog [options] [slice_file] sequences_file

    Given a sequence database in FASTA format and a list of regions
    to be sliced from those sequences, generates another FASTA file
    that contains the sliced regions.

    slice_file must be a file defining which slices we need. The
    file must be in the following format:

      seqID1 start1 end1
      seqID2 start2 end2
      ...

    When the end position is omitted, the default is the length of
    the sequence. When the start position is also omitted, the whole
    sequence will be returned.

    sequences_file must be a sequence database in FASTA format.
    """

    short_name = "seqslicer"

    def __init__(self):
        super(SeqSlicerApp, self).__init__()
        self.seqs = None

    def create_parser(self):
        """Creates the command line parser"""
        parser = super(SeqSlicerApp, self).create_parser()

        parser.add_option("-a", "--try-alternative-splicing",
                dest="try_alternative_splicing", default=False,
                action="store_true",
                help="try sequenceID.1 if sequenceID is not found")
        parser.add_option("-i", "--ignore-unknown", dest="ignore_unknown",
                default=False, action="store_true",
                help="ignore unknown sequence IDs")

        return parser

    def load_sequences(self, seq_file):
        """Loads the sequences from the given sequence file in FASTA format"""
        if isinstance(seq_file, (str, unicode)) and hasattr(SeqIO, "index"):
            self.seqs = SeqIO.index(seq_file, "fasta")
        else:
            seq_file = open_anything(seq_file)
            self.seqs = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))

    def run_real(self):
        """Runs the application and returns the exit code"""

        if len(self.args) == 1:
            slice_file = "-"
            seq_file = self.args[0]
        elif len(self.args) == 2:
            slice_file = self.args[0]
            seq_file = self.args[1]
        else:
            self.parser.print_help()
            return 1

        self.load_sequences(seq_file)
        self.process_file(slice_file)

    def process_file(self, slice_file):
        """Processes the given slice file"""
        for line in open_anything(slice_file):
            parts = line.strip().split()
            if not parts:
                continue

            seq_id, record = parts[0], None
            try:
                record = self.seqs[seq_id]
            except KeyError:
                if self.options.try_alternative_splicing:
                    try:
                        record = self.seqs[seq_id+".1"]
                    except KeyError:
                        pass

            if record is None:
                if self.options.ignore_unknown:
                    self.log.info("Ignoring unknown sequence ID: %s" % seq_id)
                    continue
                self.log.fatal("Unknown sequence ID in input file: %s" % seq_id)
                return 1

            if len(parts) == 1:
                start, end = 1, len(record.seq)
                new_id = record.id
            else:
                start = int(parts[1])
                if len(parts) == 2:
                    end = len(record.seq)
                else:
                    end = int(parts[2])
                new_id = "%s:%d-%d" % (record.id, start, end)

            new_record = SeqRecord(record.seq[(start-1):end],
                    id=new_id, name=record.name, description="")
            SeqIO.write([new_record], sys.stdout, "fasta")


if __name__ == "__main__":
    sys.exit(SeqSlicerApp().run())

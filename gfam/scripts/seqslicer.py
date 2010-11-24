#!/usr/bin/env python
"""Sequence slicer application"""

import sys

from gfam import fasta
from gfam.sequence import SeqRecord
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

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
    sequence will be returned. Sequence positions are defined by
    integers starting from 1, that is, the first amino acid at the
    N-terminus is at position 1. Both the start and end positions are
    inclusive. If you specify a negative number, this will be interpreted
    as a position from the C-terminus, that is, position -60 is position
    60 counted from the C-terminus.

    sequences_file must be a sequence database in FASTA format.
    """

    short_name = "seqslicer"

    def __init__(self, *args, **kwds):
        super(SeqSlicerApp, self).__init__(*args, **kwds)
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
        parser.add_option("-r", "--seq-id-regexp", metavar="REGEXP",
                help="remap sequence IDs using REGEXP",
                config_key="sequence_id_regexp",
                dest="sequence_id_regexp")
        parser.add_option("-k", "--keep-ids", dest="keep_ids",
                action="store_true",
                help="keep original sequence IDs even if this will duplicate "
                     "existing IDs in the output file.")

        return parser

    def load_sequences(self, seq_file):
        """Loads the sequences from the given sequence file in FASTA format"""
        self.log.info("Loading sequences from %s..." % seq_file)

        parser = fasta.Parser(open_anything(seq_file))
        parser = fasta.regexp_remapper(parser,
                self.options.sequence_id_regexp)

        self.seqs = dict(((seq.id, seq) for seq in parser))

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
        self.log.info("Processing input file: %s..." % slice_file)

        writer = fasta.Writer(sys.stdout)

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
                    self.log.warning("Ignoring unknown sequence ID: %s" % seq_id)
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

            if start == 0:
                self.log.warning("Ignoring sequence ID: %s, "
                        "requested start position is zero" % seq_id)
            elif end == 0:
                self.log.warning("Ignoring sequence ID: %s, "
                        "requested end position is zero" % seq_id)

            if start < 0:
                start = len(record.seq) + start + 1
            if end < 0:
                end = len(record.seq) + end + 1

            if not self.options.keep_ids:
                new_id = "%s:%d-%d" % (record.id, start, end)
            else:
                new_id = seq_id

            new_record = SeqRecord(record.seq[(start-1):end],
                    id=new_id, name=record.name, description="")
            writer.write(new_record)


if __name__ == "__main__":
    sys.exit(SeqSlicerApp().run())

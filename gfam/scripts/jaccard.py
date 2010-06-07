#!/usr/bin/env python

import sys

from collections import defaultdict
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

class JaccardSimilarityApp(CommandLineApp):
    """\
    Usage: %prog [options] [input_file]

    Calculates the Jaccard similarity of gene pairs from an input
    file. It is assumed that each gene is linked to itself.
    The input file must be in the following format:

        id1 id2
        id3 id4
        ...

    IDs must not contain whitespace. Everything that's after the
    second ID in a line is ignored. The output does not contain
    redundant lines, i.e. id2-id1 is not reported if id1-id2 was
    already reported. Columns of the output are separated by tab
    characters.
    """

    short_name = "jaccard"

    def create_parser(self):
        parser = super(JaccardSimilarityApp, self).create_parser()
        parser.add_option("-n", "--no-loops", dest="add_loops",
                action="store_false", default=True,
                config_key="analysis:jaccard/assume_loops",
                help="don't assume that a protein is connected to itself")
        parser.add_option("-l", "--only-linked", dest="only_linked",
                action="store_true", default=False,
                config_key="analysis:jaccard/only_linked",
                help="report only those pairs that are linked in the input file")
        parser.add_option("-m", "--min-similarity", dest="min_similarity",
                default=0, type=float, metavar="VALUE",
                config_key="analysis:jaccard/min_similarity",
                help="report only pairs with similarity not less than VALUE")
        return parser

    def run_real(self):
        """Runs the application"""
        for infile in (self.args or ["-"]):
            self.process_file(infile)

    def process_file(self, filename):
        """Processes the input file with the given filename"""
        self.log.info("Processing %s..." % filename)
        infile = open_anything(filename)
        neis = defaultdict(set)
        for line_no, line in enumerate(infile):
            parts = line.strip().split()
            if not parts:
                continue
            if len(parts) < 2:
                raise ValueError("line %d contains only a single ID" % line_no)
            neis[parts[0]].add(parts[1])
            neis[parts[1]].add(parts[0])

        if self.options.add_loops:
            for k, v in neis.iteritems():
                v.add(k)

        all_ids = sorted(neis.keys())
        lens = dict((id, len(neis1)) for id, neis1 in enumerate(neis))
        for id1 in all_ids:
            neis1 = neis[id1]
            len1 = float(len(neis1))
            if self.options.only_linked:
                others = sorted(neis1)
            else:
                others = all_ids
            for id2 in others:
                if id2 < id1:
                    continue
                neis2 = neis[id2]
                isect = len(neis2.intersection(neis1))
                sim = isect / (len1+len(neis2)-isect)
                if sim < self.options.min_similarity:
                    continue
                print "%s\t%s\t%.8f" % (id1, id2, sim)


if __name__ == "__main__":
    sys.exit(JaccardSimilarityApp().run())

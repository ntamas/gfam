#!/usr/bin/env python
"""Gene Ontology overrepresentation analysis application"""

import optparse
import sys

from collections import defaultdict
from gfam.go import Tree as GOTree
from gfam.go.overrepresentation import OverrepresentationAnalyser
from gfam.interpro import InterPro2GOMapping
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

class OverrepresentationAnalysisApp(CommandLineApp):
    """\
    Usage: %prog [options] [go_tree_file] [go_mapping_file] [input_file]

    Conducts a Gene Ontology overrepresentation analysis on the
    domain architecture of the sequences given in the input file.

    go_tree_file must represent the Gene Ontology tree in OBO format.

    go_mapping_file establishes the mapping between InterPro domain IDs
    and Gene Ontology IDs.

    Each line in the input file must be in a tab-separated format
    where the first column contains the sequence ID and the
    *third* column contains a semicolon-separated list of InterPro
    domain IDs. The output of ``find_domain_architecture.py`` can be
    used directly.
    """

    short_name = "overrep"

    def __init__(self, *args, **kwds):
        super(OverrepresentationAnalysisApp, self).__init__(*args, **kwds)
        self.go_tree = None

    def create_parser(self):
        """Creates the command line parser for this application"""
        parser = super(OverrepresentationAnalysisApp, self).create_parser()
        parser.add_option("-p", "--confidence", dest="confidence",
                metavar="PVALUE",
                help="set the p-value cutoff to the given PVALUE. "
                     "Default: %default",
                config_key="analysis:overrep/confidence",
                default=0.05, type=float)
        parser.add_option("--correction", dest="correction",
                metavar="METHOD", default="fdr",
                choices=("none", "bonferroni", "sidak", "fdr"),
                config_key="analysis:overrep/correction",
                help="use the given correction METHOD for multiple "
                     "hypothesis testing. Default: %default ")
        parser.add_option("-s", "--min-size", dest="min_size",
                metavar="SIZE", default=5,
                config_key="analysis:overrep/min_term_size",
                help="don't test for overrepresentation of GO terms "
                     "with less than SIZE annotations. Default: %default")
        return parser

    def run_real(self):
        """Runs the overrepresentation analysis application"""
        if len(self.args) != 3:
            self.error("expected exactly three input file names")

        go_tree_file, go_mapping_file, input_file = self.args

        self.log.info("Loading GO tree from %s..." % go_tree_file)
        self.go_tree = GOTree.from_obo(go_tree_file)

        self.log.info("Loading InterPro --> GO mapping from %s..." % \
                go_mapping_file)
        self.go_mapping = InterPro2GOMapping.from_file(go_mapping_file, self.go_tree)

        self.log.info("Processing domain architectures from %s..." % \
                input_file)
        return self.process_file(input_file)

    def process_file(self, input_file):
        """Processes the given input file that contains the domain
        architectures."""
        self.log.info("Running overrepresentation analysis")
        self.log.info("p-value = %.4f, correction method = %s" % \
                      (self.options.confidence, self.options.correction))

        overrep = OverrepresentationAnalyser(self.go_tree, self.go_mapping,
                confidence = self.options.confidence,
                min_count = self.options.min_size,
                correction = self.options.correction)
        cache = {}

        num_no_annotations = 0
        num_no_domains = 0
        total_seqs = 0

        for line in open_anything(input_file):
            parts = line.strip().split("\t")
            gene_id, arch = parts[0], tuple(parts[2].split(";"))
            total_seqs += 1

            if arch == ("NO_ASSIGNMENT", ):
                num_no_domains += 1
                num_no_annotations += 1
                continue

            if arch not in cache:
                cache[arch] = overrep.test_group(arch)

            print gene_id
            for term, p_value in cache[arch]:
                print "  %.4f: %s (%s)" % (p_value, term.id, term.name)
            print
            if len(cache[arch]) == 0:
                num_no_annotations += 1

        self.log.info("Total number of sequences processed: %d" % total_seqs)
        if num_no_annotations:
            self.log.info("%d sequences have no overrepresented annotations :("
                    % num_no_annotations)
        if num_no_domains:
            self.log.info("%d sequences have no domains at all :(" % num_no_domains)


if __name__ == "__main__":
    sys.exit(OverrepresentationAnalysisApp().run())

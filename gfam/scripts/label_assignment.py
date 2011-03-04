#!/usr/bin/env python
"""Gene Ontology label assignment application"""

import optparse
import sys

from collections import defaultdict
from gfam.go import Tree as GOTree
from gfam.interpro import InterPro2GOMapping
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2011, Tamas Nepusz"
__license__ = "GPL"

class LabelAssignmentApp(CommandLineApp):
    """\
    Usage: %prog [options] [go_tree_file] [go_mapping_file] [input_file]

    Assigns a few Gene Ontology labels to each sequence given in the input file
    based on their domain architecture and a mapping file that associates
    functional labels to domains.

    go_tree_file must represent the Gene Ontology tree in OBO format.

    go_mapping_file establishes the mapping between InterPro domain IDs
    and GO labels.

    Each line in the input file must be in a tab-separated format
    where the first column contains the sequence ID and the
    *third* column contains a semicolon-separated list of InterPro
    domain IDs. The output of ``find_domain_architecture.py`` can be
    used directly.
    """

    short_name = "label_assignment"

    def __init__(self, *args, **kwds):
        super(LabelAssignmentApp, self).__init__(*args, **kwds)
        self.go_tree = None

    def run_real(self):
        """Runs the label assignment application"""
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
                all_terms = set()
                for domain in arch:
                    all_terms.update(self.go_mapping.get_left(domain, []))
                for path in self.go_tree.paths_to_root(*list(all_terms)):
                    all_terms.difference_update(path[1:])
                all_terms = sorted(all_terms, key =
                        lambda x: len(self.go_mapping.get_right(x, [])))
                cache[arch] = all_terms

            print gene_id
            for term in cache[arch]:
                print "  %s (%s)" % (term.id, term.name)
            print

            if len(cache[arch]) == 0:
                num_no_annotations += 1

        self.log.info("Total number of sequences processed: %d" % total_seqs)
        if num_no_annotations:
            self.log.info("Could not assign functional label to %d sequences :("
                    % num_no_annotations)
        if num_no_domains:
            self.log.info("%d sequences have no domains at all :(" % num_no_domains)


if __name__ == "__main__":
    sys.exit(LabelAssignmentApp().run())

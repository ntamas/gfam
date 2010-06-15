#!/usr/bin/env python

import sys

from collections import defaultdict
from gfam.interpro import AssignmentParser, InterPro
from gfam.scripts import CommandLineApp
from gfam.utils import AssignmentOverlapChecker, EValueFilter, \
                       SequenceWithAssignments, UniversalSet, \
                       open_anything

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

class AssignmentSourceFilterApp(CommandLineApp):
    """\
    Usage: %prog [options] [assignment_file]

    Tries to determine a consensus domain architecture for sequences based
    on the output of IPRScan. Domains that have a corresponding InterPro ID
    have priority over those who don't have one.

    This program expects incoming assignments from the standard input or from
    a given file and print the selected ones to the standard output.

    The assignment process has three stages:

        1. For a given sequence, find all the data sources that have InterPro
           IDs for all their domains assigned to that sequence. Take the source
           having the maximal coverage, and use it as a primary assignment.

        2. Loop over the unused domains and try to augment the primary assignment
           with them. Domain insertions and overlaps are allowed only if both
           domains (the one being inserted and the one which is already there in
           the assignment) have the same data source. In this step, HMMPanther
           and Gene3D domains are excluded.

        3. Try step 2 again with HMMPanther and Gene3D domains.
    """

    short_name = "assignment_source_filter"

    def create_parser(self):
        """Creates the command line parser used by this script"""
        parser = super(AssignmentSourceFilterApp, self).create_parser()
        parser.add_option("-x", "--exclude", dest="ignored",
                metavar="SOURCE",
                help="add SOURCE to the list of ignored sources",
                config_key="analysis:iprscan_filter/untrusted_sources",
                action="append", default=[])
        parser.add_option("-e", "--e-value", dest="max_e",
                metavar="THRESHOLD",
                help="E-value THRESHOLD to filter assignments",
                config_key="analysis:iprscan_filter/e_value_thresholds",
                default="inf")
        parser.add_option("-i", "--interpro-file", dest="interpro_file",
                metavar="FILE",
                help="use the InterPro parent-child FILE to remap IDs",
                config_key="analysis:iprscan_filter/interpro_parent_child_mapping",
                default=None)
        parser.add_option("-g", "--gene-ids", dest="gene_id_file",
                metavar="FILE", help="only consider those IDs which "+
                   "are present in the list in the given FILE",
                config_key="generated/file.valid_gene_ids",
                default=None)
        parser.add_option("--max-overlap", metavar="SIZE",
                help="sets the maximum overlap size allowed between "
                     "assignments of the same data source. Default: %default",
                config_key="max_overlap",
                dest="max_overlap", type=int, default=20)
        return parser

    def run_real(self):
        """Runs the application"""
        AssignmentOverlapChecker.max_overlap = self.options.max_overlap

        if self.options.interpro_file:
            self.log.info("Loading known InterPro IDs from %s..." % \
                    self.options.interpro_file)
            self.interpro = InterPro.FromFile(self.options.interpro_file)
        else:
            self.interpro = InterPro()

        if self.options.gene_id_file:
            self.log.info("Loading sequence IDs from %s..." % \
                          self.options.gene_id_file)
            self.valid_sequence_ids = set()
            for line in open_anything(self.options.gene_id_file):
                self.valid_sequence_ids.add(line.strip())
        else:
            self.valid_sequence_ids = UniversalSet()

        self.ignored = set()
        for ignored_source in self.options.ignored:
            parts = ignored_source.split()
            self.ignored.update(parts)

        if not self.args:
            self.args = ["-"]
        if len(self.args) > 1:
            self.error("Only one input file may be given")

        self.process_infile(self.args[0])

    def process_infile(self, fname):
        self.log.info("Processing %s..." % fname)

        current_id, assignments_by_source = None, defaultdict(list)
        infile = open_anything(fname)
        parser = AssignmentParser()
        valid_ids = self.valid_sequence_ids
        evalue_filter = EValueFilter.FromString(self.options.max_e)
        for line in infile:
            line = line.strip()
            assignment = parser.parse(line)
            if assignment.source in self.ignored:
                continue
            if assignment.id not in valid_ids:
                continue
            if assignment.evalue is not None and \
                    not evalue_filter.is_acceptable(assignment):
                continue
            if assignment.id != current_id:
                for row in self.filter_assignments(current_id, assignments_by_source):
                    print row
                current_id = assignment.id
                assignments_by_source = defaultdict(list)
            assignments_by_source[assignment.source].append((assignment, line))

        for row in self.filter_assignments(current_id, assignments_by_source):
            print row


    def filter_assignments(self, name, assignments_by_source):
        """Given a sequence name and its assignments ordered in a dict by
        their sources, selects a representative assignment set based on the
        rules outlined in the documentation of `FindUnassignedApp`.
        """

        if not assignments_by_source:
            return []

        result = []

        # First, find the data source which is fully covered by InterPro IDs and
        # covers the most of the sequence
        coverage = {}
        assignments_without_interpro_ids = []
        all_sources = set()
        for source, assignments in assignments_by_source.iteritems():
            all_sources.add(source)
            # Check if all the assignments have InterPro IDs
            ok = True
            for a, line in assignments:
                if a.interpro_id:
                    continue
                if source not in ["PatternScan", "FPrintScan"]:
                    assignments_without_interpro_ids.append((a, line))
                ok = False

            # This step temporarily disabled
            # if not ok:
            #     continue

            # At this stage, we don't consider HMMPanther or Gene3D
            if source == "HMMPanther" or source == "Gene3D":
                continue

            # Calculate the coverage
            seq = SequenceWithAssignments(name, a.length)
            for a, _ in assignments:
                seq.assign(a)
            coverage[source] = seq.coverage()

        a = assignments_by_source.keys()[0]

        seq = SequenceWithAssignments(name, assignments_by_source[a][0][0].length)
        unused_assignments = []
        if coverage:
            best_source = max(coverage.keys(), key = coverage.__getitem__)
            for a, line in assignments_by_source[best_source]:
                seq.assign(a)
                tab_count = list(line).count("\t")
                if tab_count < 13:
                    line = line + "\t" * (13-tab_count)
                result.append("%s\t%s" % (line, 1))
        else:
            best_source = None

        for source, assignments in assignments_by_source.iteritems():
            if source == best_source:
                continue
            unused_assignments.extend(assignments)

        # Try to fill the unassigned regions with the rest of the assignments
        # that were unused so far, starting from the longest assignment.
        # In the first stage, we exclude HMMPanther

        unused_assignments.sort(\
              key = lambda x: x[0].end-x[0].start, reverse = True
        )

        if "HMMPanther" in all_sources or "Gene3D" in all_sources:
            stages = [all_sources - set(["HMMPanther", "Gene3D"]), all_sources]
        else:
            stages = [all_sources]

        selected_idxs = set()
        idx_to_stage = {}
        for stage_no, sources in enumerate(stages):
            for idx, (a, _) in enumerate(unused_assignments):
                if a.source not in sources:
                    continue
                if seq.assign(a):
                    selected_idxs.add(idx)
                    idx_to_stage[idx] = stage_no+2
        for idx in sorted(selected_idxs):
            row = unused_assignments[idx][1]
            tab_count = list(row).count("\t")
            if tab_count < 13:
                row = row + "\t" * (13-tab_count)
            result.append("%s\t%s" % (row, idx_to_stage[idx]))

        return result


if __name__ == "__main__":
    sys.exit(AssignmentSourceFilterApp().run())

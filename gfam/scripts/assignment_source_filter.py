#!/usr/bin/env python

import re
import sys

from collections import defaultdict
from gfam.assignment import Assignment, AssignmentOverlapChecker, \
                            EValueFilter, SequenceWithAssignments
from gfam.interpro import AssignmentReader, InterPro
from gfam.scripts import CommandLineApp
from gfam.utils import complementerset, open_anything

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

    The assignment process has three stages by default; stages can be
    configured in the config file. The default setup is as follows:

        1. For a given sequence, take the source having the maximal coverage
           and use it as a primary assignment. In this step, HMMPanther and
           Gene3D domains are excluded.

        2. Loop over the unused domains and try to augment the primary assignment
           with them. Domain insertions and overlaps are allowed only if both
           domains (the one being inserted and the one which is already there in
           the assignment) have the same data source. In this step, HMMPanther
           and Gene3D domains are still excluded.

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
        parser.add_option("--log-exclusions", dest="exclusions_log_file",
                metavar="FILE", help="log excluded sequences to the given FILE "
                   "for debugging purposes", default=None,
                config_key="DEFAULT/file.log.iprscan_exclusions")
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
            self.valid_sequence_ids = complementerset()

        if self.options.exclusions_log_file:
            self.log.info("Logging excluded sequences to %s." %
                    self.options.exclusions_log_file)
            self.exclusion_log = open(self.options.exclusions_log_file, "a+")
        else:
            self.exclusion_log = None

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
        valid_ids = self.valid_sequence_ids
        evalue_filter = EValueFilter.FromString(self.options.max_e)

        reader = AssignmentReader(fname)
        for assignment, line in reader.assignments_and_lines():
            if assignment.id != current_id:
                self.filter_and_print_assignments(current_id, assignments_by_source)
                current_id = assignment.id
                assignments_by_source = defaultdict(list)

            if assignment.source in self.ignored:
                continue
            if assignment.evalue is not None and \
                    not evalue_filter.is_acceptable(assignment):
                continue
            assignments_by_source[assignment.source].append((assignment, line))

        # ...and the last batch
        self.filter_and_print_assignments(current_id, assignments_by_source)

    def filter_assignments(self, name, assignments_by_source):
        """Given a sequence name and its assignments ordered in a dict by
        their sources, selects a representative assignment set based on the
        rules outlined in the documentation of `FindUnassignedApp`.
        """

        if not assignments_by_source:
            self.log_exclusion(name, "no assignments in the input data file "
                    "passed the filters")
            return []

        # Determine the length of the sequence (and check that the length is
        # the same across all assignments; if not, then the input file is
        # inconsistent and the sequence will be skipped).
        source = assignments_by_source.keys()[0]
        seq_length = assignments_by_source[source][0][0].length
        for source, assignments in assignments_by_source.iteritems():
            if any(assignment.length != seq_length \
                   for assignment, _ in assignments):
                self.log.warning("Sequence %s has multiple assignments with "
                                 "different sequence lengths in the "
                                 "input file, skipping" % name)
                self.log_exclusion(name, "ambiguous sequence length in input file")
                return []

        # Initially, the result is empty
        result = []

        # Set up the stages
        stages = self.get_stages_from_config()
        """
        stages = [complementerset(["HMMPanther", "Gene3D"]),
                  complementerset(["HMMPanther", "Gene3D"]),
                  complementerset()]
        """

        # The first stage is treated specially as we have to select a single
        # source thas has the largest coverage. In the remaining stages, we
        # are allowed to cherrypick from different sources.

        # First, find the data source which covers the most of the sequence
        # and is allowed in stage 1
        first_stage = stages.pop(0)
        coverage = {}
        for source, assignments in assignments_by_source.iteritems():
            # Exclude those sources that we don't consider in the first stage
            if source not in first_stage:
                continue

            # Calculate the coverage
            seq = SequenceWithAssignments(name, seq_length)
            for a, _ in assignments:
                seq.assign(a)
            coverage[source] = seq.coverage()

        # Find the source giving the best coverage, add its domains into
        # the current assignment.
        seq = SequenceWithAssignments(name, seq_length)
        if coverage:
            best_source = max(coverage.keys(), key = coverage.__getitem__)
            for a, line in assignments_by_source[best_source]:
                line = line.strip()
                seq.assign(a)
                tab_count = list(line).count("\t")
                if tab_count < 13:
                    line = line + "\t" * (13-tab_count)
                result.append("%s\t%s" % (line, 1))
        else:
            best_source = None

        # Collect the unused assignments (not from the best source)
        # into unused_assignments
        unused_assignments = []
        for source, assignments in assignments_by_source.iteritems():
            if source == best_source:
                continue
            unused_assignments.extend(assignments)

        # Try to fill the unassigned regions with the rest of the assignments
        # that were unused so far, starting from the longest assignment.
        unused_assignments.sort(key = lambda x: -x[0].get_assigned_length())

        # Okay, we're done with the first stage, process the rest.
        # idx_to_stage will contain the indices of the selected
        # assignments as keys and the number of the corresponding
        # stage in which they were selected as values.
        idx_to_stage = {}
        for stage_no, sources in enumerate(stages):
            for idx, (a, _) in enumerate(unused_assignments):
                if a.source in sources and seq.assign(a):
                    idx_to_stage[idx] = stage_no+2
        for idx in sorted(idx_to_stage.keys()):
            row = unused_assignments[idx][1].strip()
            tab_count = list(row).count("\t")
            if tab_count < 13:
                row = row + "\t" * (13-tab_count)
            result.append("%s\t%s" % (row, idx_to_stage[idx]))

        if not result:
            self.log_exclusion(name, "no assignments were selected after executing "
                    "all the stages")

        return result

    def filter_and_print_assignments(self, name, assignments_by_source):
        """Filters and prints the list of assignments of the gene with the
        given `name`. `assignments_by_source` must contain the list of
        domain assignments, sorted by data source."""
        if name is None:
            return
        if name not in self.valid_sequence_ids:
            self.log_exclusion(name, "not in the list of valid gene IDs")
            return
        for row in self.filter_assignments(name, assignments_by_source):
            print row

    def get_stages_from_config(self):
        """Turns to the configuration file specified at startup to
        fetch the data sources to be used in each stage of the algorithm.
        If there is no configuration file specified or it does not
        contain the corresponding keys, it will simply use a default
        stage setup which ignores HMMPanther and Gene3D in the first
        and second steps, but uses all sources in the third step.
        
        The method will be looking for configuration keys named like
        ``stages.1``, ``stages.2`` and so on in the ``analysis:iprscan_filter``
        section of the config file. The value of each such config key must
        be an expression consisting of assignment source names and the
        operators ``+`` and ``-``, with their usual meaning of addition
        and exclusion. The special source name ``ALL`` means all possible
        data sources, enabling us to write expressions like ``ALL-HMMPanther``
        (meaning all the sources except HMMPanther). Some examples:
            
        - ``HMMPanther`` means HMMPanther only.
        - ``ALL`` means all possible data sources.
        - ``HMMPanther+HMMPfam`` means HMMPanther or HMMPfam.
        - ``ALL-HMMPanther-Gene3D`` means all possible data sources but
          HMMPanther or Gene3D.
        - ``ALL+HMMPanther`` does not really make sense as you are extending
          all data sources with HMMPanther, so it is equivalent to ``ALL``.
          GFam will figure out what you meant anyway.
        """
        cfg = self.parser.config
        if cfg is None:
            spec = ["ALL-HMMPanther-Gene3D", "ALL-HMMPanther-Gene3D",
                    "ALL"]
        else:
            spec, idx = [], 1
            section = "analysis:iprscan_filter"
            while cfg.has_option(section, "stages.%d" % idx):
                spec.append(cfg.get(section, "stages.%d" % idx))
                idx += 1

        regexp = re.compile("([-+])?\s*([^-+]+)")
        result = []
        for item in spec:
            sources = set()
            for match in regexp.finditer(item):
                sign, source = match.groups()
                if source == "ALL":
                    source = complementerset()
                else:
                    source = set([source.strip()])
                if sign == "-":
                    sources -= source
                else:
                    sources |= source
            result.append(sources)

        return result

    def log_exclusion(self, name, reason):
        """Adds an entry to the exclusions log file, noting that the
        sequence with the given `name` was excluded from further consideration
        because of the given `reason`.

        This method works only if an exclusion log file was specified when
        calling the application.
        """
        if self.exclusion_log is None:
            return
        self.exclusion_log.write("%s: %s\n" % (name, reason))


if __name__ == "__main__":
    sys.exit(AssignmentSourceFilterApp().run())

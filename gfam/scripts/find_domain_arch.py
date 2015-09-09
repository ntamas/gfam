#!/usr/bin/env python
"""Application that calculates the domain architecture of each gene and
outputs them in a simple text-based format.
"""

from collections import defaultdict

import operator
import optparse
import sys

from gfam.assignment import AssignmentOverlapChecker, SequenceWithAssignments
from gfam.interpro import InterPro, InterProNames
from gfam.scripts import CommandLineApp
from gfam.utils import complementerset, redirected

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

class FindDomainArchitectureApp(CommandLineApp):
    """\
    Usage: %prog [options] interpro_file clustering_file

    Application that calculates the domain architecture of each gene and
    oputs them in a simple text-based format, given the filtered InterPro
    assignments and a clustering of the unknown regions in a separate file.
    """

    short_name = "find_domain_arch"

    def __init__(self, *args, **kwds):
        super(FindDomainArchitectureApp, self).__init__(*args, **kwds)
        self.seqcat = {}

    def create_parser(self):
        """Creates the command line parser for this application"""
        parser = super(FindDomainArchitectureApp, self).create_parser()
        parser.add_option("-s", "--min-size", dest="min_size",
                metavar="N",
                help="consider only clusters with at least N elements as novel domains",
                config_key="analysis:find_domain_arch/min_novel_domain_size",
                default=2, type=int)
        parser.add_option("-S", "--sequences",
                dest="sequences_file", metavar="FILE",
                help="FASTA file containing all the sequences of the representative gene model",
                config_key="file.input.sequences", default=None)
        parser.add_option("-i", "--interpro-parent-child-file",
                dest="interpro_parent_child_file",
                metavar="FILE",
                help="use the InterPro parent-child FILE to remap IDs",
                config_key="file.mapping.interpro_parent_child",
                default=None)
        parser.add_option("--details",
                dest="details", metavar="FILE",
                help="print more details about the domain architecture into FILE",
                config_key="generated/file.domain_architecture_details",
                default=None)
        parser.add_option("--stats",
                dest="stats", metavar="FILE",
                help="print genome-level statistics about the domain architectures into FILE",
                config_key="generated/file.domain_architecture_stats",
                default=None)
        parser.add_option("-n", "--names",
                dest="interpro_names_file",
                metavar="FILE",
                help="use the given FILE to assign InterPro IDs to names",
                config_key="file.mapping.interpro2name",
                default=None)
        parser.add_option("-r", "--seq-id-regexp", metavar="REGEXP",
                help="remap sequence IDs using REGEXP",
                config_key="sequence_id_regexp",
                dest="sequence_id_regexp")
        parser.add_option("--max-overlap", metavar="SIZE",
                help="sets the maximum overlap size allowed between "
                     "assignments of the same data source. Default: %default",
                config_key="max_overlap",
                dest="max_overlap", type=int, default=20)
        return parser

    def run_real(self):
        """Runs the applications"""
        if len(self.args) != 2:
            self.error("exactly two input files are expected")

        AssignmentOverlapChecker.max_overlap = self.options.max_overlap

        if self.options.interpro_parent_child_file:
            self.log.info("Loading InterPro parent-child assignments from %s..." % \
                    self.options.interpro_parent_child_file)
            self.interpro = InterPro.FromFile(self.options.interpro_parent_child_file)
        else:
            self.interpro = InterPro()

        self.interpro_names = InterProNames.FromFile(self.options.interpro_names_file)

        if self.options.details:
            self.details_file = open(self.options.details, "w")
        else:
            self.details_file = None

        interpro_file, clustering_file = self.args
        self.process_interpro_file(interpro_file)
        self.process_clustering_file(clustering_file)
        self.sort_by_domain_architecture()

        for seqs in self.domain_archs.itervalues():
            seqs.sort()

        self.domain_archs = self.domain_archs.items()
        self.domain_archs.sort(key=lambda x: len(x[1]), reverse=True)

        for domain_arch, members in self.domain_archs:
            if domain_arch:
                arch_str = ";".join(domain_arch)
            else:
                arch_str = "NO_ASSIGNMENT"
                arch_str_pos = "NO_ASSIGNMENT"
                arch_desc = "NO_DESCRIPTION"

            family_length = len(members)
            for member in members:
                seq = self.seqcat[member]
                if domain_arch:
                    arch_str_pos = ";".join(assignment.short_repr() \
                            for assignment in seq.assignments)
                    arch_desc = ";".join( \
                            self.interpro_names[assignment.domain]
                            for assignment in seq.assignments
                    )
                print "%s\t%d\t%s\t%d\t%s\t%s" % (member, seq.length, arch_str, \
                                              family_length, arch_str_pos, \
                                              arch_desc)

        self.details_file.close()

        if self.options.stats:
            stats_file = open(self.options.stats, "w")

            total_residues, covered_residues, covered_residues_nonnovel = 0.0, 0, 0
            nonnovel_sources = complementerset(["Novel"])

            for seq in self.seqcat.itervalues():
                total_residues += seq.length
                covered_residues += round(seq.coverage() * seq.length)
                covered_residues_nonnovel += round(seq.coverage(sources=nonnovel_sources) * seq.length)

            all_archs = set(arch for arch, _ in self.domain_archs)
            num_archs = len(all_archs)
            if "" in self.domain_archs:
                num_archs -= 1

            def exclude_novel_domains(domain_architecture):
                """Excludes novel domains from a domain architecture and returns
                the filtered domain architecture as a tuple."""
                return tuple(a for a in domain_architecture if not a.startswith("NOVEL"))

            archs_without_novel = set(exclude_novel_domains(arch)
                    for arch in all_archs)
            archs_without_novel.discard(())
            num_archs_without_novel = len(archs_without_novel)

            num_seqs_with_nonempty_domain_arch = \
                    sum(len(value) for key, value in self.domain_archs if key)
            num_seqs_with_nonempty_domain_arch_ignore_novel = \
                    sum(len(value) for key, value in self.domain_archs
                        if exclude_novel_domains(key) in archs_without_novel)
            num_seqs_with_nonempty_nonnovel_domain_arch = \
                    sum(len(value) for key, value in self.domain_archs
                            if key and not any(a.startswith("NOVEL") for a in key))

            with redirected(stdout=stats_file):
                print "Domain architectures"
                print "===================="
                print ""
                print "Non-empty: %d" % num_archs
                print "Non-empty (when ignoring novel domains): %d" % num_archs_without_novel
                print ""
                print "Sequences"
                print "========="
                print ""
                print "Total: %d" % len(self.seqcat)
                print "With at least one domain: %d (%.4f%%)" %\
                        (num_seqs_with_nonempty_domain_arch,
                         100.0 * num_seqs_with_nonempty_domain_arch / len(self.seqcat))
                print "With at least one non-novel domain: %d (%.4f%%)" %\
                        (num_seqs_with_nonempty_domain_arch_ignore_novel,
                         100.0 * num_seqs_with_nonempty_domain_arch_ignore_novel / len(self.seqcat))
                print "With at least one domain and no novel domains: %d (%.4f%%)" %\
                        (num_seqs_with_nonempty_nonnovel_domain_arch,
                         100.0 * num_seqs_with_nonempty_nonnovel_domain_arch / len(self.seqcat))
                print ""
                print "Residues"
                print "========"
                print ""
                print "Total: %d" % total_residues
                print "Covered: %d (%.4f%%)" % (covered_residues, 100.0*covered_residues/total_residues)
                print "Covered by non-novel: %d (%.4f%%)" % (covered_residues_nonnovel, 100.0*covered_residues_nonnovel/total_residues)
            stats_file.close()


    def process_interpro_file(self, interpro_file):
        from gfam.scripts.find_unassigned import FindUnassignedApp
        unassigned_app = FindUnassignedApp()
        unassigned_app.set_sequence_id_regexp(self.options.sequence_id_regexp)
        unassigned_app.process_sequences_file(self.options.sequences_file)
        unassigned_app.process_infile(interpro_file)
        self.seqcat = unassigned_app.seqcat
        for seq_id in set(unassigned_app.seq_ids_to_length.keys()) - set(self.seqcat.keys()):
            self.seqcat[seq_id] = SequenceWithAssignments(seq_id, \
                                  unassigned_app.seq_ids_to_length[seq_id])

    def process_clustering_file(self, fname):
        f = open(fname)
        idx = 1
        for line in f:
            ids = line.strip().split()
            if len(ids) < self.options.min_size:
                continue
            domain_name = "NOVEL%05d" % idx
            idx += 1
            for id in ids:
                seq_id, _, limits = id.rpartition(":")
                start, end = map(int, limits.split("-"))
                self.seqcat[seq_id].assign_(start, end, domain_name)
        f.close()

    def sort_by_domain_architecture(self):
        self.domain_archs = defaultdict(list)
        for seq_id, seq in self.seqcat.iteritems():
            assignments = sorted(seq.assignments, key=operator.attrgetter("start"))
            domains = []
            if self.details_file:
                print >>self.details_file, seq_id

            primary_source = set()

            new_assignments = []
            for assignment in assignments:
                new_assignment = assignment.resolve_interpro_ids(self.interpro)
                if assignment.comment == "1":
                    primary_source.add(assignment.source)
                domains.append(new_assignment.domain)
                new_assignments.append(new_assignment)
            self.domain_archs[tuple(domains)].append(seq_id)

            if not primary_source:
                primary_source = None
            else:
                primary_source = ", ".join(primary_source)

            if self.details_file:
                seq2 = SequenceWithAssignments(seq.name, seq.length)
                seq2.assignments = [assignment for assignment in assignments \
                                    if assignment.source != "Novel"]
                sources = sorted(set(assignment.source \
                        for assignment in assignments \
                        if assignment.source != "Novel"))

                print >>self.details_file, "    Primary assignment source:", primary_source
                print >>self.details_file, "    Number of data sources used:", len(sources)
                print >>self.details_file, "    Data sources: %s" % ", ".join(sources)
                print >>self.details_file, "    Coverage: %.3f" % seq.coverage()
                print >>self.details_file, "    Coverage w/o novel domains: %.3f" % seq2.coverage()
                for assignment in assignments:
                    attrs = assignment._asdict()
                    if assignment.comment is None and \
                       assignment.domain.startswith("NOVEL"):
                        attrs["comment"] = "novel"
                    row = "    %(start)4d-%(end)4d: %(domain)s "\
                          "(%(source)s, stage: %(comment)s)" % attrs
                    print >>self.details_file, row,
                    interpro_id = assignment.interpro_id
                    if not interpro_id and assignment.domain in self.interpro.mapping:
                        interpro_id = self.interpro.mapping[assignment.domain]
                    if interpro_id:
                        anc = self.interpro.tree.get_most_remote_ancestor(interpro_id)
                        if interpro_id == anc:
                            print >>self.details_file, "(InterPro ID: %s)" % anc
                        else:
                            print >>self.details_file, "(InterPro ID: %s --> %s)" % (interpro_id, anc)
                        if anc in self.interpro_names:
                            print >>self.details_file, " "*(row.index(":")+1), self.interpro_names[anc]
                    else:
                        print >>self.details_file, ""
                        if assignment.domain in self.interpro_names:
                            print >>self.details_file, " "*(row.index(":")+1), self.interpro_names[assignment.domain]
                print >>self.details_file, ""

            seq.assignments = new_assignments

if __name__ == "__main__":
    sys.exit(FindDomainArchitectureApp().run())


#!/usr/bin/env python
"""Script that runs a breadth-first search on a weighted undirected graph"""

import sys

from collections import defaultdict, deque
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything, UniqueIdGenerator

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"


class BreadthFirstSearch(object):
    """Breadth-first search implementation using an adjacency list"""

    def __init__(self, adj_list):
        self.adj_list = adj_list

    def run(self, start_vertex):
        """Runs a breadth-first search from the given start vertex and
        yields the visited vertices one by one."""
        queue = deque([start_vertex])
        visited = set([start_vertex])
        adj_list = self.adj_list

        while queue:
            vertex = queue.popleft()
            yield vertex
            unseen_neis = adj_list[vertex]-visited
            visited.update(unseen_neis)
            queue.extend(unseen_neis)


class ConnectedComponentAnalysisApp(CommandLineApp):
    """\
    Usage: %prog [options] [input_file]

    Finds the connected components of a given weighted input
    graph, optionally applying a weight threshold before the
    analysis. The input file must be in the following format:

        id1 id2 weight
        id3 id4 weight
        ...

    IDs must not contain whitespace. Everything that's after the
    weight in a line is ignored. The output will contain the list
    of connected components, one per line. Items in the same
    component are separated by tabs.
    """

    short_name = "cca"

    def create_parser(self):
        """Creates the command line parser for the CCA analysis application"""
        parser = super(ConnectedComponentAnalysisApp, self).create_parser()
        parser.add_option("-t", "--threshold", dest="threshold",
                type=float, default=0, metavar="EPS",
                config_key="analysis:cca/threshold",
                help="ignores edges with weight less than EPS")
        return parser

    def run_real(self):
        """Runs the application"""
        for infile in (self.args or ["-"]):
            self.process_file(infile)

    def process_file(self, filename):
        """Processes the input file with the given filename"""
        threshold = self.options.threshold
        adj_list = defaultdict(set)
        idgen = UniqueIdGenerator()

        self.log.info("Processing %s..." % filename)
        for line_no, line in enumerate(open_anything(filename)):
            parts = line.strip().split()
            if not parts:
                continue
            if len(parts) < 2:
                raise ValueError("line %d contains only a single ID" % line_no)
            if len(parts) < 3:
                parts.append(1.0)
            else:
                parts[2] = float(parts[2])

            id1, id2, weight = parts[:3]
            if weight < threshold:
                continue

            id1, id2 = idgen[id1], idgen[id2]
            if id1 == id2:
                continue
            if id1 > id2:
                id1, id2 = id2, id1

            adj_list[id1].add(id2)
            adj_list[id2].add(id1)

        names = idgen.values()

        bfs = BreadthFirstSearch(adj_list)
        not_seen = set(range(len(idgen)))
        while not_seen:
            component = list(bfs.run(not_seen.pop()))
            print "\t".join(names[idx] for idx in component)
            not_seen.difference_update(component)


if __name__ == "__main__":
    sys.exit(ConnectedComponentAnalysisApp().run())

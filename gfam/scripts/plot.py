"""
Plot script that is used to generate figures showing general descriptive
statistics of E-value distributions, overlap sizes and domain lengths.
"""

from __future__ import division

import matplotlib
import os
import sys

from collections import defaultdict
from gfam.assignment import Assignment, AssignmentOverlapChecker
from gfam.interpro import AssignmentReader
from gfam.scripts import CommandLineApp
from gfam.utils import Histogram, open_anything
from math import ceil, log10
from textwrap import TextWrapper

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"


def figure_name(name):
    """Returns a decorator that assigns a human-readable figure name to a
    function that generates a plot."""
    def decorator(func):
        func.figure_name = name
        return func
    return decorator


def get_subplot_sizes(num):
    """Assuming that a figure should have `num` subplots, calculates
    how many subplots should be in one row and how many rows should
    there be in total in order to have a nice arrangement. Returns
    the number of rows and columns in a tuple. Supports at most 24
    subplots."""
    if num <= 0:
        return 1, 1
    if num <= 3:
        return 1, num
    if num <= 8:
        return 2, int(ceil(num / 2.))
    if num <= 15:
        return 3, int(ceil(num / 3.))
    if num <= 24:
        return 4, int(ceil(num / 4.))
    raise ValueError("cannot place %d subplots" % num)


def match(fragment, names):
    """Given a possible list of `names` and a name `fragment` (or maybe a
    complete name), returns the name that starts with the given fragment
    if there is a unique match; otherwise, raises `ValueError`.
    """
    matches = [name for name in names if name.startswith(fragment)]
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise ValueError("fragment does not match any of the "
                         "names: %s" % fragment)
    raise ValueError("multiple matches for fragment: %s" % fragment)


class PlotApp(CommandLineApp):
    """\
    Usage: %prog [options] figure_name [figure_name] ...

    Creates one or more figures showing some descriptive statistics of the
    input data. figure_name specifies the name of the figure to be plotted.
    For a list of supported figures, use "list" as the figure name.  """

    def create_parser(self):
        """Creates the parser that parses the command line options"""
        parser = super(PlotApp, self).create_parser()
        parser.add_option("-a", "--assignment-file", dest="assignment_file",
                metavar="FILE", config_key="file.input.iprscan",
                help="read InterPro domain assignments from FILE")
        parser.add_option("-o", "--output", dest="output", metavar="FILE",
                help="save the plot to the given FILE. Supported extensions: "
                     "png, pdf, jpg, txt")
        parser.add_option("--cumulative", dest="cumulative",
                action="store_true",
                help="plot cumulative distributions (if that makes sense "
                     "for a given plot)")
        parser.add_option("--relative", dest="relative", action="store_true",
                help="plot relative frequencies instead of absolute counts "
                     "on the Y axis (if that makes sense for a given plot)")
        return parser

    def calculate_histograms_from_assignments(self, funcs, bin_size=1):
        """Calculates some measures derived from individual assignments
        and returns a histogram for each of them. `funcs` is a dict of
        names and callables. Each callable is called for each InterPro
        assignment, and their return values will be stored and sorted
        in histograms. The return value will be a dict that maps
        the function names to sub-dicts, each sub-dict mapping assignment
        sources to `Histogram` instances.

        This may sound too abstract so far, so here's an example. Let
        us assume that `funcs` is a dict with two keys: ``evalue``
        maps to a function that returns the E-value of an `Assignment`,
        while ``length`` maps to a function that returns the length of
        a domain assignment. The result will be a dict with the same
        two keys (i.e. ``evalue`` and ``length``). ``evalue`` in the
        result will correspond to a dict that maps assignment sources
        to E-value histograms, while ``length`` in the result will
        correspond to another dict that maps assignment sources to
        domain length histograms.

        `bin_size` is the bin size of the histogram.
        """
        result = defaultdict(lambda: defaultdict(
            lambda: Histogram(bin_size)
        ))
        for assignment in self.get_assignment_reader():
            for name, func in funcs.iteritems():
                value = func(assignment)
                if value is None:
                    continue
                result[name][assignment.source].add(value)
        return result

    def get_assignment_reader(self):
        """Returns an `AssignmentReader` that reads from the assignment
        file given in the command line."""
        if self.options.assignment_file is None:
            self.parser.error("must specify assignment file, use -a")
        return AssignmentReader(self.options.assignment_file)

    def get_available_figures(self):
        """Returns the list of available figures"""
        return sorted((method[5:], func) \
                      for method, func in self.__class__.__dict__.iteritems() \
                      if method.startswith("plot_") and callable(func))

    def get_barplot_from_histograms(self, histograms, xlabel=None,
            ylabel="Count", xlim=None):
        """Given a dict mapping assignment source names to `Histogram`
        instances, creates a nice plot that shows each histogram in a
        separate panel. Returns an instance of a Matplotlib `Figure`.

        `xlabel` and `ylabel` gives the labels for the X and Y axes,
        respectively.
        """
        sources = sorted(k for k in histograms.keys())

        figure = self.get_empty_figure()
        figure.subplots_adjust(wspace=0.4, hspace=0.3)
        num_fig_rows, num_fig_cols = get_subplot_sizes(len(sources))

        for idx, source in enumerate(sources):
            histogram = histograms[source]

            if self.options.cumulative:
                cumsum = 0
                points = []
                for left, _, count in histogram.bins():
                    cumsum += count
                    points.append((left, cumsum))
            else:
                points = [(left, count) for left, _, count in histogram.bins()]

            if points:
                if self.options.relative:
                    total = points[-1][1]
                else:
                    total = sum(y for _, y in points)
                points = [(x, y/total) for x, y in points]

            axes = figure.add_subplot(num_fig_rows, num_fig_cols, idx+1)
            args = zip(*points)
            args.append(histogram.bin_width)
            axes.bar(*args)
            axes.set_title(source.capitalize())
            if idx % num_fig_cols == 0:
                axes.set_ylabel(ylabel)
            if idx / num_fig_cols == num_fig_rows-1:
                axes.set_xlabel(xlabel)
            if xlim:
                axes.set_xlim(xlim)

        return figure

    def get_empty_figure(self, figsize=(12, 8), dpi=96):
        """Returns an empty Matplotlib `Figure` with the given size
        and DPI."""
        from matplotlib import pyplot
        return pyplot.figure(figsize=figsize, dpi=dpi)

    def run_real(self):
        """Runs the plotting application"""
        if not self.args:
            self.parser.error("please specify at least one figure "
                              "to be plotted")

        if self.options.output:
            _, output_ext = os.path.splitext(self.options.output)
            if output_ext:
                output_ext = output_ext[1:].lower()
            if output_ext != "txt":
                matplotlib.use("agg")
        else:
            output_ext = None

        available_figures = self.get_available_figures()
        figure_names = [name for name, _ in available_figures]
        for arg in self.args:
            try:
                method_name = match(arg, figure_names)
            except ValueError as ex:
                self.log.warning(ex)
                continue

            method = getattr(self, "plot_%s" % method_name)
            figure = method()
            if figure is None:
                continue

            if self.options.output:
                if output_ext == "txt":
                    # TODO
                    pass
                else:
                    figure.savefig(self.options.output)
            else:
                from matplotlib import pyplot
                pyplot.show()

    def plot_list(self):
        """Lists the names of the available figures"""
        wrapper = TextWrapper(subsequent_indent = " " * 22,
                              width = 78)
        for method, func in self.get_available_figures():
            if method != "list":
                wrapper.initial_indent = ("%-20s " % method).ljust(22)
                print wrapper.fill(func.figure_name)

    @figure_name("Distribution of log E-values, sorted by data sources")
    def plot_evalue_distribution(self):
        """Plots the distribution of E-values for domains from each one
        of the data sources"""

        def trimmed_log_evalue(assignment):
            if assignment.evalue is None:
                return None
            if assignment.evalue <= 1e-20:
                return -20
            return log10(assignment.evalue)

        histograms = self.calculate_histograms_from_assignments(
            {"log_evalue": trimmed_log_evalue}
        )["log_evalue"]
        return self.get_barplot_from_histograms(histograms,
                xlabel="log(E-value)", xlim=(-20, 5))

    @figure_name("Distribution of domain lengths, sorted by data sources")
    def plot_length_distribution(self):
        """Plots the distribution of domain lengths from each one
        of the data sources"""

        max_length = 750

        def trimmed_length(assignment):
            length = assignment.get_assigned_length()
            return min(length, max_length-1)

        histograms = self.calculate_histograms_from_assignments(
            {"length": trimmed_length},
            bin_size=int(max_length / 25)
        )["length"]
        return self.get_barplot_from_histograms(histograms,
                xlabel="Domain length", xlim=(0, max_length))

    @figure_name("Distribution of overlap lengths, sorted by data sources")
    def plot_overlap_distribution(self):
        """Plots the distribution of overlaps between domain assignments
        of the same data source, for each of the data sources."""
        histograms = defaultdict(lambda: Histogram(5))
        get_overlap_size = AssignmentOverlapChecker.get_overlap_size

        same_seq_assignments, prev_id = [], None
        for assignment in self.get_assignment_reader():
            if assignment.id == prev_id:
                # Calculate overlap with all the previous
                # assignments
                for other in same_seq_assignments:
                    # We care only about assignments from the same source
                    if assignment.source != other.source:
                        continue
                    overlap = get_overlap_size(assignment, other)
                    # If the overlap size equals the domain size,
                    # we skip it -- it is likely to be a duplicate entry
                    if overlap >= assignment.get_assigned_length():
                        continue
                    if overlap > 99:
                        overlap = 99
                    if overlap:
                        histograms[assignment.source].add(overlap)
                same_seq_assignments.append(assignment)
            else:
                prev_id = assignment.id
                same_seq_assignments = [assignment]
        return self.get_barplot_from_histograms(histograms,
                xlabel="Overlap length", xlim=(0, 100))

if __name__ == "__main__":
    sys.exit(PlotApp().run())


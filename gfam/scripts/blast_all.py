#!/usr/bin/env python
"""All-against-all BLASTing of a given set of sequences.

This is a convenience wrapper around `formatdb` and `blastall` that
runs an all-against-all BLAST on a given set of sequences and removes
the temporary files afterwards. It also detects whether the legacy
C-based BLAST tools or the newer C++-based tools are installed, and
adjusts the command line accordingly.
"""

from __future__ import with_statement

import os
import subprocess
import sys

from gfam import fasta
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything, search_file, temporary_dir

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

class AllAgainstAllBLASTApp(CommandLineApp):
    """\
    Usage: %prog [options] sequences_file

    Given a sequence database in FASTA format, runs an all-against-all
    BLAST query and returns the output in tabular format per default.

    sequences_file must be a sequence database in FASTA format.
    """

    short_name = "allblast"

    def create_parser(self):
        """Creates the command line parser"""
        parser = super(AllAgainstAllBLASTApp, self).create_parser()

        parser.add_option("-m", dest="blast_output_format",
                default=8, type=int, metavar="FORMAT",
                help="use FORMAT as the output format. This option "
                "is passed on intact to blastall. Default: %default",
                config_key="analysis:allblast/output_format"
        )
        parser.add_option("-a", dest="num_threads", default=1,
                type=int, metavar="N",
                help="instruct BLAST to use N threads. This option "
                "is passed on intact to blastall. Default: %default",
                config_key="analysis:allblast/num_cpu_cores"
        )
        parser.add_option("-o", dest="output_file", metavar="FILE",
                help="send the BLAST output to FILE. This option "
                "is passed on intact to blastall if present.",
                config_key="analysis:allblast/output_file"
        )
        parser.add_option("-p", dest="blast_tool",
                default="blastp", metavar="TOOL",
                help="run the given BLAST tool after formatdb. "
                "This option is passed on intact to blastall. "
                "Default: %default",
                config_key="analysis:allblast/blast_tool"
        )
        parser.add_option("--formatdb-path", dest="formatdb_path",
                help="uses PATH as the path to the formatdb executable",
                config_key="utilities/util.formatdb", metavar="PATH"
        )
        parser.add_option("--blastall-path", dest="blastall_path",
                help="uses PATH as the path to the blastall executable",
                config_key="utilities/util.blastall", metavar="PATH"
        )

        return parser

    def get_blast_cmdline(self, tool_name, args = None):
        """Given the name of a BLAST tool in `tool_name`, looks up the
        value of the corresponding option from the command line, and tries
        to figure out whether it is a path to a folder containing the tools
        or the path to the tool itself. Returns a list of arguments that
        should be passed to `subprocess.Popen` in order to execute the tool.
        `args` is a set of extra arguments that should be passed to the tool.

        If the newer, C++-based BLAST tools are installed, it automatically
        prepends the list of arguments with ``legacy_blast.pl`` to adapt
        the new tools to the calling interface of the old ones.
        """
        if not args:
            args = []

        tool_path = getattr(self.options, tool_name + "_path")
        if not tool_path:
            tool_path = os.getcwd()

        if os.path.isfile(tool_path):
            # The path exists and points to a real file, so we simply return
            return [tool_path] + args

        if not os.path.exists(tool_path):
            # The path does not exist. Assume that this is a file referring to
            # the old BLAST tools, but the user has the new one. Try to extract
            # the folder and see if the folder exists.
            base, tool_name = os.path.split(tool_path)

        if os.path.isdir(tool_path):
            # The path exists and points to a folder
            base = os.path.normpath(tool_path)

        # Okay, first check for the actual tool in the folder
        full_path = os.path.join(tool_path, tool_name)
        if os.path.isfile(full_path):
            return [full_path] + args

        # Nope, try the legacy script
        full_path = os.path.join(tool_path, "legacy_blast.pl")
        if os.path.isfile(full_path):
            return [full_path, tool_name] + args + ["--path", tool_path]

        # Nothing succeeded
        return None


    def run_real(self):
        """Runs the application and returns the exit code"""
        if not self.args:
            self.args = ["-"]

        # Find formatdb and blastall in the current path if needed, ensure that
        # the paths are absolute.
        for util in ["formatdb", "blastall"]:
            optkey = "%s_path" % util
            path = getattr(self.options, optkey)
            if not path:
                path = search_file(util, executable=True)
                if not path:
                    self.error("cannot find %s in system path" % util)
            setattr(self.options, optkey, os.path.abspath(path))

        for arg in self.args:
            if not self.process_file(arg):
                return 1

    def process_file(self, sequence_file):
        """Processes the given sequence file"""
        self.log.info("Processing file: %s..." % sequence_file)
        sequence_file = os.path.abspath(sequence_file)

        with temporary_dir(change=True) as dirname:
            self.log.debug("Using temporary directory: %s" % dirname)
            ok = self.run_formatdb(sequence_file)
            ok = ok and self.run_blastall(sequence_file)
            if not ok:
                return False

        return True

    def run_formatdb(self, sequence_file):
        """Runs ``formatdb`` on the given sequence file.

        Returns ``True`` if the execution was successful, ``False`` otherwise.
        """
        self.log.info("Invoking formatdb...")

        args = ["-n", "database", "-i", sequence_file, "-p", "T", "-o", "F"]
        args = self.get_blast_cmdline("formatdb", args)
        if not args:
            self.log.fatal("cannot find formatdb in %s" % self.options.formatdb_path)
            return False

        formatdb = subprocess.Popen(args, stdin=open(os.devnull), \
                stdout=sys.stderr)
        retcode = formatdb.wait()
        if retcode != 0:
            self.log.fatal("formatdb exit code was %d, exiting..." % retcode)
            return False

        self.log.info("formatdb returned successfully.")
        return True

    def run_blastall(self, sequence_file):
        """Runs ``blastall`` on the given sequence file.

        Returns ``True`` if the execution was successful, ``False`` otherwise.
        """
        self.log.info("Invoking blastall, this might take a long time...")

        args = []
        args.extend(["-p", self.options.blast_tool])
        args.extend(["-d", "database", "-i", sequence_file])
        args.extend(["-m", str(self.options.blast_output_format)])
        args.extend(["-a", str(self.options.num_threads)])
        if self.options.output_file:
            args.extend(["-o", self.options.output_file])

        args = self.get_blast_cmdline("blastall", args)
        if not args:
            self.log.fatal("cannot find blastall in %s" % self.options.blastall_path)
            return False

        blastall = subprocess.Popen(args, stdin=open(os.devnull))
        retcode = blastall.wait()
        if retcode != 0:
            self.log.fatal("blastall exit code was %d, exiting..." % retcode)
            return False

        self.log.info("blastall returned successfully.")
        return True

if __name__ == "__main__":
    sys.exit(AllAgainstAllBLASTApp().run())

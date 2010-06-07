"""
Classes and utilities commonly used in GFam command line scripts
"""

from gfam.config import ConfigurableOptionParser
from textwrap import dedent

import logging

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["CommandLineApp"]

class CommandLineApp(object):
    """Generic command line application class"""

    def __init__(self):
        self.log = self.create_logger()
        self.options, self.args = None, None

    def create_parser(self):
        """Creates a command line parser for the application"""
        doc = self.__class__.__doc__
        parser = ConfigurableOptionParser(usage=dedent(doc).strip())
        parser.add_option("-v", "--verbose", dest="verbose",
                action="store_true", help="verbose logging")
        return parser

    def create_logger(self):
        """Creates a logger for the application"""
        if hasattr(self.__class__, "short_name"):
            log_name = self.__class__.short_name
        else:
            log_name = self.__class__.__module__
        log = logging.getLogger(log_name)
        log.setLevel(logging.WARNING)
        logging.basicConfig(level=logging.DEBUG, format="%(message)s")
        return log

    def run(self, args=None):
        """Runs the application. This method processes the command line using the
        command line parser and as such, it should not be overridden in child
        classes unless you know what you are doing. If you want to implement
        the actual logic of your application, override `run_real` instead."""
        self.parser = self.create_parser()
        self.options, self.args = self.parser.parse_args(args)

        if self.options.verbose:
            self.log.setLevel(logging.INFO)

        return self.run_real()

    def run_real(self):
        self.log.info("Nothing to do.")
        return 0


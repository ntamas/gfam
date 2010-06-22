"""
Classes and utilities commonly used in GFam command line scripts.

The scripts implementing the individual steps of the GFam pipeline
are designed in a way that they can either be invoked in standalone
mode (like any other command line script) or they can be instantiated
and used from other Python modules. Each command line script is
derived from the `CommandLineApp` class, which takes care of
implementing functionality common for all the scripts, such as:

- providing a command line parser (an instance of `ConfigurableOptionParser`)
- providing a logger instance
- defining methods for extending the default option parser and for
  signaling fatal errors to the caller
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
    """Generic command line application class that provides
    common functionality for all GFam command line scripts.
    """

    def __init__(self, logger=None):
        """Creates a command line script instance that will use
        the given logger to log messages. If `logger` is ``None``,
        it will create a logger exclusively for the script.

        The following fields will be set up:

        - `self.log`: the logger where messages should be logged to
        - `self.parser`: an instance of `ConfigurableOptionParser`
          to parse the command line options. This is set up in
          `create_parser()`, not in the constructor.
        - `self.options`: the parsed command line options. Parsing
          occurs when `run()` is called, so this field contains
          ``None`` when `run()` has not been called.
        - `self.args`: the positional arguments from the parsed
          command line.
        """
        if logger:
            self.log = logger
        else:
            self.log = self.create_logger()
        self.options, self.args = None, None

    def create_parser(self):
        """Creates a command line parser for the application.

        By default, the parser will use the class docstring as
        help string and add two options to the parser: ``-v``
        specifies verbose logging mode, while ``-d`` specifies
        debug mode. In verbose mode, log messages having a level
        above or equal to `logging.INFO` are printed. In debug
        mode, all log messages (including debug messages) are
        printed. The default is to print warnings and errors
        only.

        The created command line parser will be returned to the
        caller.
        """
        doc = self.__class__.__doc__
        parser = ConfigurableOptionParser(usage=dedent(doc).strip())
        parser.add_option("-v", "--verbose", dest="verbose",
                action="store_true", help="verbose logging")
        parser.add_option("-d", "--debug", dest="debug",
                action="store_true", help="show debug messages")
        return parser

    def create_logger(self):
        """Creates a logger for the application and returns it."""
        if hasattr(self.__class__, "short_name"):
            log_name = self.__class__.short_name
        else:
            log_name = self.__class__.__module__

        log = logging.getLogger(log_name)
        log.setLevel(logging.WARNING)
        logging.basicConfig(level=logging.DEBUG, format="%(message)s")

        return log

    def error(self, message):
        """Signals a fatal error and shuts down the application."""
        self.parser.error(message)

    def run(self, args=None):
        """Runs the application. This method processes the command line using
        the command line parser and as such, it should not be overridden in
        child classes unless you know what you are doing. If you want to
        implement the actual logic of your application, override `run_real`
        instead.
        
        `args` contains the command line arguments that should be parsed.
        If `args` is ``None``, the arguments will be obtained from
        ``sys.argv[1:]``."""
        self.parser = self.create_parser()
        self.options, self.args = self.parser.parse_args(args)

        if self.options.verbose:
            self.log.setLevel(logging.INFO)
        if self.options.debug:
            self.log.setLevel(logging.DEBUG)

        return self.run_real()

    def run_real(self):
        self.log.info("Nothing to do.")
        return 0


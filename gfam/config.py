"""Configuration-related classes for GFam.

This module provides `ConfigurableOptionParser`, an extension of Python's
built-in `optparse.OptionParser` that lets you provide default values for
command-line options from a given configuration file.
"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["ConfigurableOption", "ConfigurableOptionParser"]

from ConfigParser import ConfigParser
from optparse import OptionParser, Option

import sys


class ConfigurableOption(Option):
    """An extension of the `Option` class of Python that stores a
    configuration key which can be used to fetch the value of the
    option if it is not given in the command line.
    """

    def __init__(self, *args, **kwds):
        try:
            self.config_key = kwds["config_key"]
            del kwds["config_key"]
        except KeyError:
            self.config_key = None

        Option.__init__(self, *args, **kwds)
        self._seen = False

    def get_config_section_and_item(self):
        """Returns the section and item in which we have to look for
        the value of this option."""
        if self.config_key is None:
            return None, None
        parts = self.config_key.split("/", 1)
        if len(parts) == 1:
            return "DEFAULT", parts[0]
        if not parts[0]:
            parts[0] = "DEFAULT"
        return tuple(parts)

    def process(self, *args, **kwds):
        """Method invoked by the command line parser when it sees
        this option on the command line.

        This method simply makes note of the fact that the option
        was seen on the command line, and calls the superclass.
        """
        self._seen = True
        return Option.process(self, *args, **kwds)

    @property
    def seen(self):
        """Returns whether this option was seen on the command line."""
        return self._seen


class ConfigurableOptionParser(OptionParser):
    """An extension of the `OptionParser` class of Python that also uses
    a config file.

    Basically, the scripts of GFam can be driven both from command line
    arguments and configuration files. This class lets you specify the
    command line arguments using the same syntax you would use with
    `OptionParser`, but it also lets you attach a configuration key
    to each of the options. It also registers an option ``-c`` and a
    long option ``--config-file`` that can be used to specify the
    input configuration file. If an option is not present on the
    command line, this option parser will try to look up the corresponding
    configuration key in the given configuration file.

    Configuration keys must be specified in slashed format (i.e.
    ``section/item``).

    This class also exposes an instance attribute named `config`, which
    contains the parsed configuration values from the specified configuration
    file. `config` will be an instance of `ConfigParser` or ``None`` if
    `parse_args()` was not called so far.
    """

    def __init__(self, *args, **kwds):
        if "option_class" not in kwds:
            kwds["option_class"] = ConfigurableOption

        OptionParser.__init__(self, *args, **kwds)
        self.config = None

        self.add_option("-c", "--config-file", dest="config_file",
                help="name of the configuration FILE", metavar="FILE",
                default=None)

    def parse_args(self, *args, **kwds):
        """Parses the command line and returns a tuple ``(options, args)``,
        where ``options`` contains the values of the parsed command line
        options (after adding the values from the config file for missing
        options) and ``args`` contains the list of positional arguments.
        """
        options, args = OptionParser.parse_args(self, *args, **kwds)

        if not options.config_file:
            return options, args

        # Load the configuration file
        config = ConfigParser()
        config.read(options.config_file)

        # Extend the options with the ones given in the config file
        for option in self.option_list:
            if not hasattr(option, "seen") or option.seen:
                continue

            section, item = option.get_config_section_and_item()
            if section is None:
                continue

            if config.has_option(section, item):
                option.process(option.config_key, config.get(section, item), \
                        options, self)

        self.config = config

        return options, args


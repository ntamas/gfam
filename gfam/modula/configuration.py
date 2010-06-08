"""
Configuration file reader for modula
"""

import os

from ConfigParser import SafeConfigParser as Parser

class Configuration(object):
    """Class storing configuration details for Modula"""

    def __init__(self, rootdir=".", cfg=None):
        defaults = {
            "@paths.modules": "modules",
            "@paths.storage": "storage"
        }

        if cfg:
            self.cfg = cfg
        else:
            cfgfiles = ["modules.cfg", "modula.cfg"]
            cfgfiles = [os.path.join(rootdir, file) for file in cfgfiles]

            self.cfg = Parser()
            self.cfg.read(cfgfiles)

        # Add defaults
        for k, v in defaults.iteritems():
            if k not in self:
                self[k] = v

        # Make sure all paths are absolute
        for k, v in self.items("@paths"):
            self["@paths.%s" % k] = os.path.abspath(os.path.expanduser(v))

        # Make sure all input files are absolute
        for k, v in self.items("@inputs"):
            self["@inputs.%s" % k] = os.path.abspath(os.path.expanduser(v))


    def _parse_name(self, name):
        if "." in name:
            section, option = name.split(".", 1)
        else:
            section = "@global"
            option = name
        return section, option

    def get(self, name): return self.cfg.get(*self._parse_name(name))
    def getInt(self, name): return self.cfg.getint(*self._parse_name(name))
    def getFloat(self, name): return self.cfg.getfloat(*self._parse_name(name))
    def getBoolean(self, name): return self.cfg.getboolean(*self._parse_name(name))
    def items(self, section): return self.cfg.items(section)

    def __contains__(self, name):
        return self.cfg.has_option(*self._parse_name(name))

    def __getitem__(self, name):
        try:
            return self.cfg.get(*self._parse_name(name))
        except:
            raise KeyError(name)

    def __setitem__(self, name, value):
        section, option = self._parse_name(name)
        if not self.cfg.has_section(section):
            self.cfg.add_section(section)
        self.cfg.set(section, option, value)

    def __delitem__(self, name):
        section, option = self._parse_name(name)
        try:
            ok = self.cfg.remove_option(section, option)
            if not ok: raise NoOptionError
        except:
            raise KeyError, name

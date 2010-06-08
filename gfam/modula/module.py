"""
Classes for module handling in Modula
"""

from gfam.modula.log import get_logger
from gfam.modula.storage import NotFoundError

import inspect
import os
import sys

from ConfigParser import NoSectionError

class Module(object):
    """This class represents a Modula module (calculation, local file etc)"""

    def __init__(self, name, config):
        self.name = name
        self.config = config
        self.logger = get_logger("module.%s" % self.name)

    def get_dependencies(self):
        raise NotImplementedError

    def get_last_updated_at(self):
        raise NotImplementedError


class LocalFile(Module):
    """This class represents a local file used as a constant input in Modula"""

    def __init__(self, name, config):
        super(LocalFile, self).__init__(name, config)
        if config is None:
            self.filename = name
        else:
            self.filename = self.config["@inputs.%s" % name]

    def get_dependencies(self):
        """Returns an empty list - local files never depend on anything"""
        return []

    def get_last_updated_at(self):
        """Returns the modification time of the file as a UNIX timestamp"""
        return os.stat(self.filename).st_mtime

    def get_handle(self, mode="rb"):
        """Returns a handle to this file with the given mode"""
        return file(self.filename, mode)


class CalculationModule(Module):
    """This class represents a Modula module that performs some calculation"""

    def __init__(self, name, config, module_dir=None, parameters=None):
        Module.__init__(self, name, config)

        if module_dir:
            self.module_dir = module_dir
        else:
            self.module_dir = self.config["@paths.modules"]

        self.module = None

        if parameters is None:
            parameters = self.get_default_parameters()
        self.parameters = parameters

    def __del__(self):
        if inspect is None:
            return
        if inspect.ismodule(self.module):
            key = None
            for key, mod in sys.modules.iteritems():
                if mod == self.module:
                    del sys.modules[key]
                    break

    def get_dependencies(self):
        """Returns the names of the modules this module depends on"""
        try:
            depends = self.config.get("%s.depends" % self.name)
        except:
            depends = ""
        depends = [str(x).strip() for x in depends.split(",") if x]
        return depends

    def get_last_updated_at(self):
        """Returns when the module result file was last updated for the
        current parameter set"""
        # At this point, the Modula engine should be initialized
        from gfam.modula import storage_engine as storage
        
        try:
            return storage.get_result_modification_time(self)
        except NotFoundError:
            return -1

    def get_default_parameters(self):
        """Retrieves the default module parameter dict for this module"""
        try:
            params = dict(self.config.items(self.name))
            params["depends"] = self.get_dependencies()
        except NoSectionError:
            params = {}
        return params

    def prepare(self):
        if not self.module:
            sys.path.insert(0, self.module_dir)
            self.module = __import__(self.name, fromlist=[], level=0)
            del sys.path[0]

    def run(self):
        self.logger.info("Starting module %s" % self.name)
        self.prepare()
        result = self.module.run(self.parameters)
        self.logger.info("Finished module %s" % self.name)
        return result


class ModuleManager(object):
    """Abstract module manager object.

    Module managers are responsible for resolving module names to
    :class:`Module` objects, checking whether a given module exists,
    and listing the available modules in the system.
    """

    def __init__(self, config):
        """Initializes the module manager from the given configuration"""
        self.config = config
        self.logger = get_logger("module")

    def get(self, module_name):
        """Retrieves a module by name and returns a :class:`Module`
        instance"""
        raise NotImplementedError

    def has_module(self, module_name):
        """Checks if a given module exists"""
        raise NotImplementedError

    def list_modules(self):
        """Returns a list containing all the valid module names"""
        raise NotImplementedError


class DefaultModuleManager(ModuleManager):
    """Default module manager implementation."""

    def __init__(self, config, module_factory=CalculationModule):
        """Initializes the module manager from the given configuration.
        
        `config` is the `Configuration` object to be used.
        `module_factory` is a factory method that creates `Module`
        instances. This will be used by `DefaultModuleManager.get`."""
        super(DefaultModuleManager, self).__init__(config)
        self.module_dir = config["@paths.modules"]
        self.module_factory = module_factory
        self.logger.debug("Using module path: %s" % self.module_dir)

    def get(self, module_name):
        """Retrieves a module by name and returns a :class:`Module`
        instance"""
        if self.has_module(module_name):
            result = self.module_factory(module_name, self.config)
        elif self.config["@inputs.%s" % module_name]:
            result = LocalFile(module_name, self.config)
        else:
            raise KeyError, module_name
        return result

    def has_module(self, module_name):
        """Checks if a given module exists"""
        if os.path.exists(os.path.join(self.module_dir, "%s.py" % module_name)):
            return True
        return False

    def list_modules(self):
        """Returns a list containing all the valid module names"""
        import fnmatch
        result = []
        for f in os.listdir(self.module_dir):
            if fnmatch.fnmatch(f, '*.py'):
                result.append(f[0:-3])
        return result


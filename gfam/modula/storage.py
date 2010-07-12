"""
Storage classes for Modula
"""

from gfam.modula.hash import sha1
from gfam.modula.log import get_logger

import os

try:
    import cPickle as pickle
except ImportError:
    # Python 3
    import pickle

class NotFoundError(Exception):
    def __init__(self, name, resource=None, nestedExc = None):
        if resource:
            msg = "%s (module/source name: %s)" % (resource, name)
        else:
            msg = "module %s" % resource

        if nestedExc:
            if isinstance(nestedExc, Exception):
                nestedExc = nestedExc.__class__.__name__
            msg += ", original cause: %s" % nestedExc
        Exception.__init__(self, msg)


class AbstractStorageEngine(object):
    """Abstract storage engine.

    This class defines the interface of storage engines. Each
    storage must implement the methods in this class.
    """

    def __init__(self, config, module_manager):
        self.config = config
        self.module_manager = module_manager
        self.logger = get_logger("storage")

    def get_source(self, source_name, mode="rb"):
        """Retrieves a data source. This method is common for
        all storage engines and it shouldn't be overridden.
        Returns a handle (a file-like object) that can be used to
        read from the source (but not to write to it). `mode`
        can be used to override the mode used to open files."""
        try:
            module = self.module_manager.get(source_name)
            return module.get_handle(mode)
        except Exception, ex:
            raise NotFoundError(source_name, source_name, ex)

    def get_result(self, module, parameters=None):
        """Retrieves the result object of a given module with the
        given parameters. If the parameters are not given,
        the defaults are fetched from the current configuration."""
        raise NotImplementedError

    def get_result_modification_time(self, module, parameters=None):
        """Retrieves the last modification time of the result object
        of a given module with the given parameters. If the parameters
        are not given, the defaults are fetched from the current
        configuration. The result is returned as a UNIX timestamp."""
        raise NotImplementedError

    def get_result_stream(self, module, parameters=None, mode="rb"):
        """Retrieves the result object of a given module with the
        given parameters. If the parameters are not given,
        the defaults are fetched from the current configuration.
        """
        raise NotImplementedError

    def get(self, module_or_source, *args, **kwds):
        """Retrieves the result object of a given module or the
        contents of a given input as a file-like object. `parameters`
        are only used when querying a module. If the parameters are not
        given and the name specified refers to a module, the defaults
        are fetched from the current configuration."""
        module = self.module_manager.get(module_or_source)
        if hasattr(module, "parameters"):
            return self.get_result(module, *args, **kwds)
        else:
            return self.get_source(module_or_source, *args, **kwds)

    def store(self, module, parameters, result):
        """Stores the result of the given module with the given parameter
        set."""
        raise NotImplementedError


class DiskStorageEngine(AbstractStorageEngine):
    def __init__(self, *args, **kwds):
        """Creates a new disk storage engine based on the given
        configuration"""
        AbstractStorageEngine.__init__(self, *args, **kwds)
        self.storage_dir = self.config["@paths.storage"]
        self.storage_hash = sha1
        if not os.path.isdir(self.storage_dir):
            self.logger.warning("Creating storage path: %s" % self.storage_dir)
            os.makedirs(self.storage_dir)
        self.logger.debug("Using storage path: %s" % self.storage_dir)

    def _get_module_result_filename(self, module, parameters=None):
        """Retrieves the name of the result file corresponding to the
        given module."""
        if parameters is None:
            parameters = module.parameters
        keys = sorted(parameters.keys())
        hash = self.storage_hash()
        hash.update(",".join("%r=%r" % (key, parameters[key]) \
                             for key in keys))
        hash = hash.hexdigest()
        return os.path.join(self.storage_dir, module.name, hash)


    def get_result(self, module, parameters=None):
        """Retrieves the result object of a given module with the
        given parameters. If the parameters are not given,
        the defaults are fetched from the current configuration.
        `module` must be an instance of :class:`modula.module.Module`
        Raises :exc:`NotFoundError` if no such module or result file
        exists."""
        return pickle.load(self.get_result_stream(module, "rb"))

    def get_result_stream(self, module, parameters=None, mode="rb"):
        """Retrieves the result object of a given module with the
        given parameters. If the parameters are not given,
        the defaults are fetched from the current configuration.
        `module` must be an instance of :class:`modula.module.Module`
        Raises :exc:`NotFoundError` if no such module or result file
        exists."""
        fname = self._get_module_result_filename(module, parameters)

        if "w" in mode:
            dir = os.path.split(fname)[0]
            if not os.path.isdir(dir):
                os.makedirs(dir)

        return open(fname, mode)

    def get_result_modification_time(self, module):
        fname = self._get_module_result_filename(module)
        try:
            return os.stat(fname).st_mtime
        except OSError, ex:
            raise NotFoundError(module, fname, ex)

    def store(self, module, result):
        """Stores the result of the given module with the given parameter
        set on the disk. `module` must be an instance of
        :class:`modula.module.Module`."""
        fname = self._get_module_result_filename(module)
        dir = os.path.split(fname)[0]
        if not os.path.isdir(dir):
            os.makedirs(dir)
        pickle.dump(result, open(fname, "wb"), pickle.HIGHEST_PROTOCOL)


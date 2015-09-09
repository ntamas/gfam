"""
Modula -- a modular calculation framework for Python

Especially for scientific calculations and stuff
"""

from gfam.modula.configuration import Configuration
from gfam.modula.log import init_master_logger, ColoredConsoleHandler
from gfam.modula.module import DefaultModuleManager
from gfam.modula.storage import DiskStorageEngine

from collections import deque

import logging
import optparse
import os
import sys

__version__ = "0.1"

config = None
logger = None
module_manager = None
storage_engine = None

def init(configuration=".", \
         module_manager_factory=DefaultModuleManager,
         storage_engine_factory=DiskStorageEngine,
         debug=False):
    """Initializes the Modula engine

    `configuration` is either a `Configuration` instance or the name of the
    directory containing the module configuration file.

    `module_manager_factory` is a factory routine that constructs
    concrete `ModuleManager` instances. It is safe to leave it at its
    default value.

    `debug` is ``True`` if debug messages should be printed; otherwise
    it is ``False``.
    """
    global config, module_manager, storage_engine, logger

    if isinstance(configuration, (str, unicode)):
        config = Configuration(rootdir=configuration)
    else:
        config = Configuration(cfg=configuration)

    logger = init_master_logger(debug=debug)
    module_manager = module_manager_factory(config)
    storage_engine = storage_engine_factory(config, module_manager)


def init_project(rootdir):
    """Initializes a Modula project in a directory"""
    import stat

    if not os.path.isdir(rootdir):
        os.mkdir(rootdir)

    for dir in ["lib", "figures", "modules", "storage"]:
        d = os.path.join(rootdir, dir)
        if not os.path.isdir(d):
            os.mkdir(d)

    f = os.path.join(rootdir, "modules.cfg")
    if not os.path.isfile(f):
        f = file(f, "w")
        print >>f, """[@global]

[@inputs]
# Enter input files here in the following format:
# id1=path
# id2=path
# ...
"""
        f.close()


def main():
    """Main entry point when Modula is run from the command line"""
    parser = optparse.OptionParser(usage="%prog [options] [command]")
    parser.add_option("-f", "--force", dest="force", action="store_true",
            help="force the execution of the given command")

    # Parse command line
    options, args = parser.parse_args()
    # Initialize the Modula engine
    init()
    # Extend the Python path
    sys.path.insert(0, os.path.join(os.getcwd(), 'lib'))

    # Start the shell
    from modula.shell import Shell
    if args:
        if options.force:
            args.append("--force")
        Shell().onecmd(" ".join(args))
    else:
        Shell().cmdloop()


def run(module_name, force=False):
    """Runs the given module in the Modula framework"""
    global config, module_manager, storage_engine, logger

    # Check dependencies, run them if needed
    to_check, to_run = deque([module_name]), []
    while to_check:
        name = to_check.popleft()

        module = module_manager.get(name)
        last_updated_at = module.get_last_updated_at()
        depends = module.get_dependencies()

        should_run = False
        for dependency in depends:
            dep = module_manager.get(dependency)
            if dep.get_last_updated_at() >= last_updated_at:
                should_run = True
                if dependency not in to_check:
                    to_check.append(dependency)
            else:
                logger.debug("%s is newer than %s, not running" % (dependency, name))

        if should_run:
            to_run.append(name)

    to_run.reverse()

    if not to_run and force:
        to_run = [module_name]

    if not to_run:
        logger.info("Nothing to do")
        return

    # Run the modules that we collected, store the results
    for name in to_run:
        module = module_manager.get(name)
        result = module.run()
        if result is not None:
            storage_engine.store(module, result)


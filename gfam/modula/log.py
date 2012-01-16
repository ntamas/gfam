"""
Logging classes for Modula
"""

import copy
import logging

__all__ = ["get_logger", "ColoredConsoleHandler"]


class NullLoggingHandler(logging.Handler):
    """Null logging handler that's attached to Modula
    loggers per default"""

    __instance = None

    def __init__(self, *args, **kwds):
        if "safeguard" not in kwds:
            raise ValueError("please do not instantiate %s directly",
                    self.__class__.__name__)
        del kwds["safeguard"]
        logging.Handler.__init__(self, *args, **kwds)

    @classmethod
    def instance(cls):
        if cls.__instance is None:
            cls.__instance = cls(safeguard=True)
        return cls.__instance

    def emit(self, record):
        pass


class ColoredConsoleHandler(logging.StreamHandler):
    """Stream handler that uses a colored output using ANSI colors."""

    _has_colors = None
    def __init__(self, *args, **kwds):
        """Initializes the handler and checks whether colored output is available"""
        logging.StreamHandler.__init__(self, *args, **kwds)
        if ColoredConsoleHandler._has_colors is None:
            has_colors = True
            try:
                import curses
            except ImportError:
                has_colors = False
            if has_colors:
                try:
                    curses.initscr()
                    has_colors = curses.has_colors()
                    curses.endwin()
                except curses.error:
                    has_colors = False
            ColoredConsoleHandler._has_colors = has_colors
        self.uses_colors = has_colors

    def emit(self, record):
        """Emits the given logging message"""
        if not self.uses_colors:
            return logging.StreamHandler.emit(self, record)

        my_record = copy.copy(record)
        level = my_record.levelno
        if level >= logging.ERROR:
            color = '\x1b[31m\x1b[1m'
        elif level >= logging.WARNING:
            color = '\x1b[33m\x1b[1m'
        elif level >= logging.INFO:
            color = '\x1b[0m'
        elif level >= logging.DEBUG:
            color = '\x1b[35m'
        else:
            color = '\x1b[0m'

        msg = my_record.msg
        if msg and msg[0] == "[":
            try:
                end_pos = msg.index("]")
                msg = "'\x1b[35m%s\x1b[0m%s" % (msg[:end_pos+1], msg[end_pos+1:])
            except ValueError:
                pass
        my_record.msg = "%s%s\x1b[0m" % (color, msg)
        logging.StreamHandler.emit(self, my_record)


def get_logger(name=None):
    """Retrieves a logger with the given name and ensures that it
    is handled by a `NullLoggingHandler`. This is to avoid error
    messages being printed when the host application using `modula`
    does not use the `logging` module.
    """
    if name:
        name = "modula.%s" % name

    logger = logging.getLogger(name)
    logger.addHandler(NullLoggingHandler.instance())
    logger.setLevel(logging.DEBUG)
    return logger

def init_master_logger(datefmt="%Y/%m/%d %H:%M:%S", debug=False):
    """Initializes the master logger of modula, using the given
    date format. `debug` is ``True`` if the logger should print
    debug messages, otherwise it is ``False``."""
    logger = get_logger()
    handler = ColoredConsoleHandler()
    handler.setLevel([logging.INFO, logging.DEBUG][bool(debug)])
    if handler.uses_colors:
        formatter = logging.Formatter("[%(asctime)s] %(message)s",
                datefmt=datefmt)
    else:
        formatter = logging.Formatter("[%(asctime)s] [%(levelname)1.1s] %(message)s",
                datefmt=datefmt)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger



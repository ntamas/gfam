#!/usr/bin/env python
import sys

try:
    from gfam.scripts.download_names import DownloadNamesApp
except ImportError:
    # Insert the parent directory of the master script into the Python path
    # and try again
    from os.path import dirname, join
    sys.path.insert(0, join(dirname(sys.modules[__name__].__file__), ".."))
    from gfam.scripts.download_names import DownloadNamesApp

if __name__ == "__main__":
    sys.exit(DownloadNamesApp().run())


#!/usr/bin/env python
"""Usage: %prog [options] [result_file]

Filters the given BLAST result file (which must be in tabular
format) and drops matches not satisfying the criteria given
in the options.
"""

from gfam.scripts.blast_filter import BlastFilterApp
import sys

if __name__ == "__main__":
    sys.exit(BlastFilterApp().run())

#!/usr/bin/env python
"""
Generates a reStructuredText documentation from the configuration file
template in `gfam.scripts.master`_.
"""

import re
import string

from gfam.scripts.master import CONFIGURATION_FILE_TEMPLATE
from textwrap import dedent, TextWrapper

def main():
    regex = re.compile(r"""
        (   (?P<comment>((^\#.*)\n?)+)      # consecutive comment lines
          \n(?P<name>[a-zA-Z0-9._:]+)\s*    # variable names
          =[\t ]*(?P<value>.*)              # variable values
        ) | (                               # OR
          \[(?P<section>[a-zA-Z0-9._ :]+)\] # section names
        )""", re.VERBOSE | re.MULTILINE)
    strip_hash = re.compile("^#", re.MULTILINE)
    current_section = ""

    keys = {}
    for match in regex.finditer(CONFIGURATION_FILE_TEMPLATE):
        if match.group("section"):
            current_section = match.group('section')
            if current_section not in keys:
                keys[current_section] = []
            continue

        comment = re.sub(strip_hash, "", match.group('comment'))
        comment = dedent(comment)

        keys[current_section].append(
            (match.group('name'), match.group('value'), comment)
        )

    wrapper = TextWrapper(initial_indent     = "    ",
                          subsequent_indent  = "    ")

    for section in sorted(keys.keys()):
        if not keys[section]:
            continue

        title = "Section ``%s``" % section
        print "%s\n%s\n" % (title, "^"*len(title))
        for name, value, comment in keys[section]:
            print "``%s``" % name
            for para in comment.split("\n\n"):
                if para[-1] in string.lowercase:
                    para = para+"."
                print wrapper.fill(para)
                print "    "
            if value:
                match = re.match(r"%\(([-a-zA-Z0-9._:]+)\)s$", value)
                if match:
                    value = "same as ``%s``" % match.group(1)
                else:
                    value = "``%s``" % value
                print "    Default value: %s" % value
                print "    "

if __name__ == "__main__":
    main()
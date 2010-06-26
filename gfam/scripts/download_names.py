#!/usr/bin/env python

import re
import sys

from gfam.scripts import CommandLineApp
from gfam.utils import open_anything

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

class DownloadNamesApp(CommandLineApp):
    """\
    Usage: %prog [options]

    Downloads the ID to human readable name assignments from various data
    sources including InterPro, Pfam and Superfamily (which derives
    its IDs from SCOP sunids). Prints the assignments in tabular format
    onto the standard output. The first column of the output is the ID,
    the second is the corresponding name.
    """

    short_name = "download_names"

    def run_real(self):
        """Runs the application"""

        urls = {"interpro":
                    "ftp://ftp.ebi.ac.uk/pub/databases/interpro/names.dat",
                "pfam":
                    "http://pfam.sanger.ac.uk/families?output=text",
                "smart":
                    "http://smart.embl-heidelberg.de/smart/descriptions.pl",
                "superfamily":
                    "http://scop.mrc-lmb.cam.ac.uk/scop/parse/"
                }

        for method_name in dir(self):
            if method_name.startswith("download_"):
                getattr(self, method_name)(urls[method_name[9:]])

    def download_interpro(self, url):
        """Downloads the official InterPro ID-name mappings from the given
        URL and prints the mapping to the standard output.
        """
        self.log.info("Downloading InterPro names from %s..." % url)
        for line in open_anything(url):
            sys.stdout.write(line)

    def download_pfam(self, url):
        """Downloads the official PFam ID-name mappings from the given URL
        and prints the mapping to the standard output.
        """
        self.log.info("Downloading PFam names from %s..." % url)
        for line in open_anything(url):
            if line[0] == "#":
                continue
            parts = line.split("\t", 2)
            if len(parts) < 3 or not parts[2]:
                continue
            sys.stdout.write("%s\t%s" % (parts[0], parts[2]))

    def download_smart(self, url):
        """Downloads the official Smart ID-name mappings from the given
        URL and prints the mapping to the standard output.
        """
        self.log.info("Downloading Smart names from %s..." % url)
        for line in open_anything(url):
            parts = line.split("\t", 3)
            if len(parts) < 3 or not parts[2]:
                continue
            if parts[1] == "ACC" and parts[2] == "DEFINITION":
                continue
            sys.stdout.write("%s\t%s\n" % (parts[1], parts[2]))

    def download_superfamily(self, url):
        """Downloads the most recent mappings from SCOP sunids to human
        readable names, derives the Superfamily IDs from the SCOP sunids,
        and prints the mapping to the standard output.

        The most recent version of SCOP is identified by applying a regexp
        to the HTML output of the given page. The regexp assumes that the
        most recent description file is linked on the given page using
        an ``<a>`` tag with ``href`` equal to ``dir.des.scop.txt_X.XX``,
        where ``X.XX`` is the version number. If such an identification
        fails, this method will return without printing anything but a
        warning on the logging stream.
        """
        self.log.info("Downloading SCOP page from %s..." % url)
        contents = open_anything(url).read()

        des_link_regexp = re.compile(r"<a href=\"dir.des.scop.txt_([0-9.]+)\">",
                re.IGNORECASE)

        max_version, max_version_nosplit = None, None
        for idx, match in enumerate(des_link_regexp.finditer(contents)):
            version = [int(comp) for comp in match.group(1).split(".")]
            if version > max_version:
                max_version = version
                max_version_nosplit = match.group(1)
                max_idx = idx

        version = max_version_nosplit

        if version is None:
            self.log.warning("Cannot infer the most recent version of SCOP, "
                             "skipping Superfamily IDs")
            return

        self.log.info("Most recent SCOP version is assumed to be %s" % version)
        url = "%sdir.des.scop.txt_%s" % (url, version)

        for line in open_anything(url):
            if line[0] == '#':
                continue
            parts = line.split("\t", 4)
            if parts[1] != "sf" or not parts[4]:
                continue
            sys.stdout.write("SSF%s\t%s" % (parts[0], parts[4]))

if __name__ == "__main__":
    sys.exit(DownloadNamesApp().run())


"""FASTA parsing routines"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["Parser", "Writer"]

from gfam.sequence import SeqRecord, Sequence

from textwrap import TextWrapper


class Parser(object):
    """Parser for FASTA files"""

    def __init__(self, handle):
        """Initializes a FASTA parser that will read from the given
        file handle"""
        self.handle = handle

    def _lines(self):
        """Iterator that iterates over the interesting lines in a
        FASTA file. This will skip empty lines and trims the remaining
        ones from all unnecessary whitespace. It will also skip the
        lines before the first record enty."""
        handle = self.handle
        while True:
            line = handle.readline()
            if not line:
                return
            if line[0] == ">":
                break

        yield line
        while True:
            line = handle.readline()
            if not line:
                return
            yield line.rstrip().replace("\r", "")


    def sequences(self):
        """Returns a generator that iterates over all the sequences
        in the FASTA file."""
        seq_record, seq = None, []
        for line in self._lines():
            if line[0] == ">":
                # Starting a new sequence here
                if seq_record:
                    seq_record.seq = Sequence("".join(seq))
                    yield seq_record
                descr = line[1:]
                seq_id = descr.split()[0]
                seq_record = SeqRecord(Sequence(), id=seq_id, name=seq_id, \
                                       description=descr)
                seq = []
            else:
                # Appending to an existing sequence
                seq.extend(list(line))

        # Don't forget the last sequence
        if seq_record:
            seq_record.seq = Sequence("".join(seq))
            yield seq_record

    def __iter__(self):
        return self.sequences()

    @classmethod
    def to_dict(cls, *args, **kwds):
        """Creates a `SeqDict` out of a FASTA file.

        The resulting `SeqDict` will map sequence IDs to their corresponding
        `SeqRecord` objects. The arguments are passed on intact to the
        constructor of `FASTAParser`.
        """
        result = {}
        for seq in cls(*args, **kwds):
            result[seq.id] = seq
        return result

class Writer(object):
    """Writes `SeqRecord` objects in FASTA format to a given file handle"""

    def __init__(self, handle):
        self.handle = handle
        self.wrapper = TextWrapper(width=70)

    def write(self, seq_record):
        if seq_record.description:
            print ">%s" % seq_record.description
        else:
            print ">%s" % (seq_record.id, )
        print "\n".join(self.wrapper.wrap(seq_record.seq))


def test():
    """Short self-test routines"""
    from urllib2 import urlopen

    f = urlopen("ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Nanoarchaeum_equitans/NC_005213.ffn")
    for seq in Parser(f):
        print seq.seq

if __name__ == "__main__":
    test()

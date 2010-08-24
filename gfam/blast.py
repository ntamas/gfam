"""
Classes and functions related to BLAST files and utilities
"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

from gfam import fasta
from gfam.utils import open_anything

class BlastFilter(object):
    """Filters BLAST records, i.e. drops the ones that do not
    satisfy some criteria.
    
    You can tune the filter with the following instance variables:
        
    - ``min_sequence_identity``: the minimum sequence identity
      required by the filter
          
    - ``min_alignment_length``: the minimum alignment length
      required by the filter
          
    - ``max_e_value``: the maximum E-value required by the
      filter
          
    You can also ask the filter to normalize the alignment length to
    between zero and one by calling `set_normalize_func`.
    """

    short_name = "blast_filter"

    def __init__(self):
        self.seq_ids_to_length = {}
        self.normalize_func = None
        self.min_sequence_identity = 0.
        self.min_alignment_length = 0.
        self.max_e_value = float('inf')
        self.set_normalize_func("off")

    def accepts(self, line):
        """Returns ``True`` if the filter accepts the given line,
        ``False`` otherwise.
        """
        line = line.strip()
        if not line or line[0] == "#":
            return True

        parts = line.split("\t")
        if float(parts[2]) < self.min_sequence_identity:
            return False

        if float(parts[10]) > self.max_e_value:
            return False

        norm = self.normalize_func(parts[0], parts[1], int(parts[3]))
        if norm < self.min_alignment_length:
            return False

        return True

    def load_sequences(self, seq_generator):
        """Loads the sequences yielded by the sequence generator.

        This method iterates over the given sequences and stores their
        names and lengths in an internal dict. The dict will be used 
        later to normalize sequences by their lengths. It is
        necessary to call this method if you are using any of the
        normalizing functions."""
        for seq in seq_generator:
            self.seq_ids_to_length[seq.id] = float(len(seq.seq))

    def load_sequences_from_file(self, fname):
        """Loads the sequences from the given file. The file must
        be in FASTA format. You are allowed to pass file pointers
        or names of gzipped/bzipped files here."""
        return self.load_sequences(fasta.Parser(open_anything(fname)))

    def _normalize_smaller(self, query_id, hit_id, length):
        """Calculates a normalized alignment length by dividing the
        unnormalized length with the length of the smaller sequence"""
        max_length = min(self.seq_ids_to_length[query_id], \
                         self.seq_ids_to_length[hit_id])
        return length / max_length

    def _normalize_larger(self, query_id, hit_id, length):
        """Calculates a normalized alignment length by dividing the
        unnormalized length with the length of the larger sequence"""
        max_length = max(self.seq_ids_to_length[query_id], \
                         self.seq_ids_to_length[hit_id])
        return length / max_length

    def _normalize_query(self, query_id, _, length):
        """Calculates a normalized alignment length by dividing the
        unnormalized length with the length of the query sequence"""
        return length / self.seq_ids_to_length[query_id]

    def _normalize_hit(self, _, hit_id, length):
        """Calculates a normalized alignment length by dividing the
        unnormalized length with the length of the hit sequence"""
        return length / self.seq_ids_to_length[hit_id]

    def set_normalize_func(self, name):
        """Sets the normalizing function used by the filter when
        filtering by the length of match. The following normalizing
        functions are known:

          - ``off``: returns the unnormalized match length

          - ``smaller``: divides the match length by the length of the
            smaller sequence

          - ``larger``: divides the match length by the length of the
            larger sequence

          - ``query``: divides the match length by the length of the
            query sequence

          - ``hit``: divides the match length by the length of the hit
            sequence

        Alternatively, you may pass a callable in place of the function
        name. The callable must accept four arguments, the first is the
        `BlastFilter` object itself, the second is the length of the
        query sequence, the third is the length of the hit sequence,
        the fourth is the unnormalized match length.
        """
        if hasattr(name, "__call__"):
            self.normalize_func = name
        else:
            funcs = { \
                "off": lambda _1, _2, length: length, \
                "smaller": self._normalize_smaller, \
                "larger": self._normalize_larger, \
                "query": self._normalize_query, \
                "hit": self._normalize_hit \
            }
            self.normalize_func = funcs[name.lower()]




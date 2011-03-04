"""
Overrepresentation analysis of Gene Ontology terms
"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "MIT"

__all__ = ["OverrepresentationAnalyser"]

from collections import defaultdict
from gfam.utils import bidict
from math import exp, log
from operator import itemgetter

try:
    from scipy.special import gammaln
except ImportError:
    def gammaln(n):
        """Logarithm of Euler's gamma function for discrete values."""
        if n < 1:
            return float('inf')
        if n < 3:
            return 0.0
        c = [76.18009172947146, -86.50532032941677, \
             24.01409824083091, -1.231739572450155, \
             0.001208650973866179, -0.5395239384953 * 0.00001]
        x, y = float(n), float(n)
        tm = x + 5.5
        tm -= (x + 0.5) * log(tm)
        se = 1.0000000000000190015
        for j in range(6):
            y += 1.0
            se += c[j] / y
        return -tm + log(2.5066282746310005 * se / x)
        
def logchoose(n, k):
    """Calculates the logarithm of n-choose-k"""
    lgn1 = gammaln(n+1)
    lgk1 = gammaln(k+1)
    lgnk1 = gammaln(n-k+1)
    return lgn1 - (lgnk1 + lgk1)

def hypergeom_pmf(k, M, n, N):
    """Hypergeometric probability moment function"""
    tot, good = M, n
    bad = tot - good
    return exp(logchoose(good, k) + logchoose(bad, N-k) - logchoose(tot, N))

def hypergeom_sf(k, M, n, N):
    """Tail distribution of a hypergeometric distribution. This is
    used to calculate p-values in the overrepresentation analysis
    tests"""
    tot, good = M, n
    bad = tot - good
    den = logchoose(tot, N)
    result = 0.
    for x in xrange(k, N+1):
        a = exp(logchoose(good, x) + logchoose(bad, N-x) - den)
        result += a
    return result

class OverrepresentationAnalyser(object):
    """Performs overrepresentation analysis of Gene Ontology
    terms on sets of entities that are annotated by some
    GO terms."""

    def __init__(self, tree, mapping, confidence=0.05, min_count=5, \
            correction="fdr"):
        """Initializes the overrepresentation analysis algorithm by associating
        it to a given Gene Ontology tree and a given mapping from entities to
        their respective GO terms.

        `tree` must be an instance of `gfam.go.Tree`. `mapping` must be
        a bidirectional dictionary object (`gfam.utils.bidict`) that
        maps entities to GO terms and vice versa. For `mapping`, if an entity
        is annotated by a GO term, it is not necessary to list all the
        ancestors of that GO term for that entity, this will be taken care of
        by the class itself which always works on the copy of the given
        mapping.
        
        `confidence` is the confidence level of the test.
        
        `min_count` specifies which GO terms are to be excluded from the
        overrepresentation analysis; if a GO term occurs less than `min_count`
        times in `mapping`, it will not be considered.
        
        `correction` specifies the multiple hypothesis testing correction to be
        used and it must be one of the following:

        - ``None``, ``"none"`` or ``"off"``: no correction

        - ``"fdr"``: controlling the false discovery rate according
          to the method of Benjamini and Hochberg

        - ``"bonferroni"``: Bonferroni correction of the family-wise
          error rate

        - ``"sidak"``: Sidak correction of the family-wise error rate
        """
        self.tree = tree
        self.mapping = self._propagate_go_term_ancestors(mapping)
        self.confidence = float(confidence)
        self.min_count = max(1, int(min_count))
        self.correction = correction

    def _propagate_go_term_ancestors(self, mapping):
        """Given a mapping object which maps entities to GO terms, this
        method ensures that all the ancestor terms of each GO term
        appear at the respective entities. Returns a copy of the mapping
        which contains these modifications."""
        result = bidict(mapping)
        for entity, terms in result.iteritems_left():
            ancestors = self.tree.ancestors(*terms)
            result.add_left_multi(entity, ancestors)
        return result

    def enrichment_p(self, term_or_id, count, group_size):
        """Calculates the enrichment p-score of the given GO term or ID
        (`term_or_id`) if it occurs `count` times in a group of size
        `group_size`.
        """
        term = self.tree.ensure_term(term_or_id)
        objs = self.mapping.right[term]
        return hypergeom_sf(count, self.mapping.len_left(), \
                            len(objs), group_size)

    def test_counts(self, counts, group_size):
        """Given a dict that maps Gene Ontology terms to their
        occurrence counts and the number of entities in the group
        from which these term counts originate, calculates a
        list of overrepresented Gene Ontology terms."""
        confidence = self.confidence
        correction = self.correction
        min_count = self.min_count

        # Determine how many tests will be performed
        if correction:
            num_tests = sum(1 \
                    for term, count in counts.iteritems() \
                    if len(self.mapping.right[term]) >= min_count
            )
            if num_tests == 0:
                return []

        # If we are doing Bonferroni or Sidak correction, adjust the
        # confidence level
        if correction == "bonferroni":
            confidence /= num_tests
        elif correction == "sidak":
            confidence = 1 - (1. - confidence) ** (1. / num_tests)

        # Do the testing
        result = []
        for term, count in counts.iteritems():
            if len(self.mapping.right[term]) < min_count:
                continue
            p = self.enrichment_p(term, count, group_size)
            result.append((term, p))

        # Filter the results
        if correction == "fdr":
            result.sort(key = itemgetter(1))
            num_tests = float(num_tests)
            for k in xrange(len(result)):
                if result[k][1] > confidence * ((k+1) / num_tests):
                    result = result[0:k]
                    break
        else:
            result = [item for item in result if r[1] <= confidence]
            result.sort(key = itemgetter(1))
            if correction == "bonferroni":
                result = [(c, p * num_tests) for c, p in result]
            elif correction == "sidak":
                result = [(c, 1 - (1 - p) ** num_tests) for c, p in result]

        return result

    def test_group(self, group):
        """Overrepresentation analysis of the given group of objects.
        `group` must be an iterable yielding objects that are in
        `self.mapping.left`.
        """
        counts = defaultdict(int)
        group_size = 0
        for item in group:
            terms = self.mapping.left.get(item, [])
            for go_term in terms:
                counts[go_term] += 1
            group_size += 1
        return self.test_counts(counts, group_size)

import numpy as np
import unittest

from Bio.Seq import Seq
from scipy.sparse import csr_matrix
from strkernel.gappy_trie import gappypair_kernel as gt
from unittest import TestCase


class Test_Gappy_Kernel(TestCase):
    def test_gappy_trie(self):
        sequences = [Seq("ACGTCGATGC"), Seq("GTCGATAGC"), Seq("GTCGaaagATAGC")]
        gappy_trie = gt(sequences,k=1,t=0,g=1, gapDifferent = False)
        # Note, not all possible k-mers are considered, so k-mers which do not occur in any sequence should not be in
        # the expected results
        expected = np.array([[0,1,2,1,1,2,1,1,2,2,0,2,2],[1,1,1,1,1,1,0,1,2,2,1,1,2],[1,1,1,1,1,1,0,1,2,2,1,1,2]])
        expected = csr_matrix(expected)

        self.assertTrue(0 == (expected != csr_matrix(gappy_trie)).getnnz())

    def test_gappy_trie_flanking(self):
        sequences = [Seq("ACGTCGATGC"), Seq("GTCGATAGC"), Seq("GTCGaaagATAGC")]
        gappy_trie = gt(sequences,k=1,t=0,g=1, gapDifferent = False, include_flanking = True)
        expected = np.array([[0,1,2,1,1,2,1,1,2,2,0,2,2],[1,1,1,1,1,1,0,1,2,2,1,1,2],[5,1,3,1,1,1,0,3,2,2,1,1,2]])
        expected = csr_matrix(expected)

        self.assertTrue(0 == (expected != csr_matrix(gappy_trie)).getnnz())

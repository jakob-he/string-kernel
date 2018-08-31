from unittest import TestCase
import numpy as np
from scipy.sparse import csr_matrix

from strkernel.lib.motif import Motif
from strkernel.lib.motiftrie import MotifTrie
from strkernel.motifkernel import motifKernel


class Test_Motif_Kernel(TestCase):
    def test_motif(self):
        motif = "A[CG]T"
        test_motif = Motif(motif)
        self.assertTrue(test_motif._motif == ["A", "CG", "T"])

    def test_trie(self):
        trie = MotifTrie(["A[CG]T", "C.G", "C..G.T", "G[A][AT]", "GT.A[CA].[CT]G"])
        self.assertTrue(np.array_equal(
            trie.check_for_motifs("AGTCTGCTTGCT"), [1, 1, 1, 0, 0]))

    def test_kernel(self):
        motifs = ["A[CG]T", "C.G", "C..G.T", "G[A][AT]"
                  "GT.A[CA].[CT]G"]
        sequences = ["ACGTCGATGC", "GTCGATAGC", "GCTAGCacgtaCGC",
                     "GTAGCTgtgcGTGcgt", "CGATAGCTAGTTAGC"]

        motif_kernel = motifKernel(motifs)
        matrix_1 = motif_kernel.compute_matrix(sequences)

        row = np.array([3, 4, 4])
        col = np.array([1, 0, 2])
        data = np.array([1, 1, 1])
        matrix_2 = csr_matrix((data, (row, col)),shape = (5,4))
        self.assertTrue(np.array_equal(matrix_1.toarray(), matrix_2.toarray()))

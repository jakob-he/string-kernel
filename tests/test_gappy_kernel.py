import numpy as np
import unittest

from Bio.Seq import Seq
from scipy.sparse import csr_matrix
from strkernel.gappy_kernel import gappypair_kernel as gk
from unittest import TestCase


class Test_Gappy_Kernel(TestCase):
    def test_gappy_kernel(self):
        sequences = [Seq("ACGTCGATGC"), Seq("GTCGATAGC"), Seq("GTCGaaagATAGC")]
        gappy_kernel = gk(sequences,k=1,t=0,g=1, gapDifferent = False, sparse = False)
        expected = np.array([[0,1,2,1,1,0,2,1,1,2,0,2,0,2,2,0],[1,1,1,1,1,0,1,0,1,2,0,2,1,1,2,0],[1,1,1,1,1,0,1,0,1,2,0,2,1,1,2,0]])
        expected = expected.astype(float)

        self.assertTrue(np.array_equal(expected, gappy_kernel))

    def test_gappy_kernel_sparse(self):
        sequences = [Seq("ACGTCGATGC"), Seq("GTCGATAGC"), Seq("GTCGaaagATAGC")]
        gappy_kernel = gk(sequences,k=1,t=0,g=1, gapDifferent = False)
        expected = np.array([[0,1,2,1,1,0,2,1,1,2,0,2,0,2,2,0],[1,1,1,1,1,0,1,0,1,2,0,2,1,1,2,0],[1,1,1,1,1,0,1,0,1,2,0,2,1,1,2,0]])
        expected = expected.astype(float)
        expected = csr_matrix(expected)

        self.assertTrue(0 == (expected != gappy_kernel).getnnz())


    def test_gappy_kernel_gapDifferent(self):
        sequences = [Seq("ACGTCGATGC")]
        gappy_kernel = gk(sequences,k=1,t=0,g=1, sparse = False)
        expected = np.array([[0,1,0,1,0,0,2,0,1,1,0,1,0,1,1,0,0,0,2,0,1,0,0,1,0,1,0,1,0,1,1,0]])
        expected = expected.astype(float)

        self.assertTrue(np.array_equal(expected, gappy_kernel))

    def test_gappy_kernel_reverse(self):
        sequences = [Seq("ACGTCGATGC")]
        gappy_kernel = gk(sequences,k=1,t=0,g=1, gapDifferent = False, reverse = True, sparse = False)
        expected = np.array([[0,3,3,1,3,0,2,0,3,2,0,0,0,0,0,0]])
        expected = expected.astype(float)

        self.assertTrue(np.array_equal(expected, gappy_kernel))

    def test_gappy_kernel_flanking(self):
        sequences = [Seq("ACGTCGatgC")]
        gappy_kernel = gk(sequences,k=1,t=0,g=1, gapDifferent = False, sparse = False, include_flanking = True)
        expected = np.array([[0,1,2,1,1,0,2,1,1,2,0,2,0,2,2,0]])
        expected = expected.astype(float)

        self.assertTrue(np.array_equal(expected, gappy_kernel))

    def test_gappy_kernel_bigger_gap(self):
        sequences = [Seq("ACGTCGATGC")]
        gappy_kernel = gk(sequences,k=1,t=0,g=3, gapDifferent = False, sparse = False)
        expected = np.array([[0,3,2,2,1,1,4,2,2,3,2,2,1,2,2,1]])
        expected = expected.astype(float)

        self.assertTrue(np.array_equal(expected, gappy_kernel))

    def test_gappy_kernel_bigger_k(self):
        sequences = [Seq("ACGTCG")]
        gappy_kernel = gk(sequences,k=2,t=0,g=1, gapDifferent = False, sparse = False)
        expected = np.zeros((1,256))
        expected[0,27] = 1.0
        expected[0,29] = 1.0
        expected[0,109] = 1.0
        expected[0,102] = 1.0
        expected[0,182] = 1.0

        self.assertTrue(np.array_equal(expected, gappy_kernel))

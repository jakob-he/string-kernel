from unittest import TestCase
import unittest

from strkernel.mismatch_kernel import preprocess, MismatchKernel
import strkernel.mismatch_kernel

class Test_Mismatch_Kernel(TestCase):
  def test_preprocess(self):
    sequence = ['aACGTt', 'AACGTT']
    int_seq = preprocess(sequence, ignoreLower=False)
    self.assertEqual(int_seq, [[0,0,1,2,3,3],[0,0,1,2,3,3]])

    # When ignoring lower case, the length of first string is 4,
    # which is shorter than the 2nd string.
    # The last 2 digits are ignored in this case.
    ignore_lower = preprocess(sequence)
    self.assertEqual(ignore_lower, [[0,1,2,3],[0,0,1,2]])

  def test_kernel(self):
    sequence = ['ACGT', 'ACGT', 'CATG']
    matrix = MismatchKernel(l=4, k=3, m=1).get_kernel(preprocess(sequence))
    self.assertEqual(matrix.kernel[0,1], 1)
    self.assertLess(matrix.kernel[0,2], 1)

if __name__ == '__main__':
    unittest.main()

"""
 Module: mismatch string kernel
 Implementattion of mismatch string kernel
 Author: Meng Zhang <RainaMeng@outlook.com>
 Reference:
  Leslie C.S., Eskin E., Cohen A., Weston J., Noble W.S.
  Mismatch string kernels for discriminative protein classification.
  Bioinformatics. 2004;20:467–476. doi: 10.1093/bioinformatics/btg431.
 <https://papers.nips.cc/paper/2179-mismatch-string-kernels-for-svm-protein-classification.pdf>
"""

from strkernel.lib.mismatchTrie import MismatchTrie
import numpy as np


def integerized(sequence):
    """
    Convert the character string into numeric string.
    """

    key_dict = sorted(set(sequence))
    int_seq = []
    for char in sequence:
        to_int = key_dict.index(char)
        int_seq.append(to_int)

    return int_seq

def preprocess(sequences, ignoreLower=True):
    """
    Data preprocessing for string sequences.
    Convert lower case into upper case if 'ignoreLower' is chosen as 'False',
    else lower case is ignored(default).
    """

    upper_seq = []
    len_record = []
    for seq in sequences:
        if ignoreLower:
            seq = [x for x in seq if 'A' <= x <= 'Z']
        else:
            seq = seq.upper()
        upper_seq.append(integerized(seq))
        len_record.append(len(seq))

    length_used = min(len_record)
    post_seq = []
    for seq in upper_seq:
        seq = seq[:length_used]
        post_seq.append(seq)

    return post_seq


def normalize_kernel(kernel):
    """
    Normalizes a kernel[x, y] by doing:
    kernel[x, y] / sqrt(kernel[x, x] * kernel[y, y])
    """

    nkernel = np.copy(kernel)

    assert nkernel.ndim == 2
    assert nkernel.shape[0] == nkernel.shape[1]

    for i in range(nkernel.shape[0]):
        for j in range(i + 1, nkernel.shape[0]):
            q = np.sqrt(nkernel[i, i] * nkernel[j, j])
            if q > 0:
                nkernel[i, j] /= q
                nkernel[j, i] = nkernel[i, j]  # symmetry

    # Set diagonal elements as 1
    np.fill_diagonal(nkernel, 1.)

    return nkernel


class MismatchKernel(MismatchTrie):
    """
    Python implementation of Mismatch String Kernels.
    Parameters
    ----------
    l: int, optional (default None), size of alphabet.
       Examples of values with a natural interpretation:
       2: for binary data
       256: for data encoded as strings of bytes
       4: for DNA/RNA sequence (bioinformatics)
       20: for protein data (bioinformatics)
    k: int, optional (default None), the k in 'k-mer'.
    m: int, optional (default None)
       maximum number of mismatches for 2 k-mers to be considered 'similar'.
       Normally small values of m should work well.
       Plus, the complexity of the algorithm is exponential in m.
    **kwargs: dict, optional (default empty)
              optional parameters to pass to `tree.MismatchTrie` instantiation.

    Attributes
    ----------
    `kernel`: 2D array of shape (n_sampled, n_samples), estimated kernel.
    `n_survived_kmers`: number of leafs/k-mers that survived trie traversal.
    """

    def __init__(self, l=None, k=None, m=None, **kwargs):

        if not None in [l, k, m]:

            # invoke trie.MismatchTrie constructor
            MismatchTrie.__init__(self, **kwargs)

            # sanitize alphabet size
            if l < 2:
                raise ValueError(
                    "Alphabet too small. l must be at least 2; got %i" % l)

            # sanitize kernel parameters (k, m)
            if 2 * m > k:
                raise ValueError(
                    ("You provided k = %i and m = %i. m is too big (must"
                     "be at ""must k / 2). This doesn't make sense.") % (k, m))

            self.l = l
            self.k = k
            self.m = m

    def get_kernel(self, X, normalize = True, **kwargs):
        """
        Main calling function to get mismatch string kernel.
        """

        if isinstance(X, tuple):
            assert len(X) == 5, "Invalid model."
            self.l, self.k, self.m, self.leaf_kmers, self.kernel = X
            # sanitize the types and shapes of self.l, self.j, self.m,
            # self.leaf_kmers, and self.kernel
        else:
            # traverse/build trie proper
            for x in ['l', 'k', 'm']:
                if not hasattr(self, x):
                    raise RuntimeError(
                        ("'%s' not specified during object initialization."
                         "You must now specify complete model (tuple of l, "
                         "k, m, leafs, and, kernel).") % x)
            self.kernel, _, _ = self.traverse(
                X, self.l, self.k, self.m, **kwargs)

            if normalize:
            # normalize kernel
                self.kernel = normalize_kernel(self.kernel)

            # gather up the leafs
            self.leaf_kmers = dict((leaf.full_label,
                                    dict((index, len(kgs)) for index, kgs
                                           in leaf.kmers.items()))
                                     for leaf in self.leafs())

        return self

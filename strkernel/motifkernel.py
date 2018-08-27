#!/usr/bin/env python3
'''
Motif Kernel Module.
'''
#standard libraries
import re

# own libraries
from strkernel.lib.motiftrie import MotifTrie

# 3rd party libraries
import numpy as np
from scipy.sparse import csr_matrix


class motifKernel:
    """
    This is the main class used for the motif kernel construction. The idea is to construct a Trie from a set of Motifs and then
    use this Trie to compute the motif content of a sequnece. The motif content can then be used to compute the similarity between
    sequences and therefore is able to serve as input for machine learning algorithms.

    Example:
    An elaborate example can be found in the ``Tutorials`` sections.

    Standard Use::

        motifKernel(motifs)
    """

    def __init__(self, motifs: [str]):
        self.motif_trie = motif_trie = MotifTrie(motifs)

    def compute_matrix(self, sequences: [str], include_flanking: bool = True, return_kernel_matrix: bool = False):
        """
        Computes the motif content of a set of sequences and returns a sparse matrix which can be used as input
        for machine learning approaches. The sparse matrix has only been tested with algorithms from the python
        package *sklearn*. Optional parameters allow the computation of a kernel matrix and inclusion of flanking regions.

        Args:
            **sequences:** A list of strings that are of the same alphabet as the motifs used to construct the motif Trie.

            **include_flanking:** Option to include or disregard the flanking regions. Default is True.

            **return_kernel_matrix:** A boolean value that indicates if the function should return a sparse matrix with the similarities between sequences (True) or a sparse matrix where each row contains the motif content of a sequence (False). Default is False.

        Returns:
            **csr_matrix:** A sparse matrix object containg either the kernel matrix (*return_kernel_matrix* = True) or
            the motif content of each sequence.
        """
        if include_flanking:
            sequences = [seq.upper() for seq in sequences]
        else:
            sequences = [re.sub('[^A-Z]', '', seq) for seq in sequences]

        search_results = [self.motif_trie.check_for_motifs(sequence) for sequence in sequences]

        if return_kernel_matrix:
            kernel_matrix = csr_matrix(np.einsum('ij,kj->ik', search_results,search_results))
            return kernel_matrix
        else:
            return csr_matrix(search_results)

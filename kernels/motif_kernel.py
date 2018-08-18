#!/usr/bin/env python3
'''
Implementation of the motif kernel.
[Asa Ben-Hur and Douglas L. Brutlag. Remote homology detection: a motif based approach. Bioinformatics,
19:26–33, 2003.]
'''
# own libraries
from kernels.lib.trie import Trie

# 3rd party libraries
import numpy as np
from scipy.sparse import csr_matrix


def processMotifs(motifs: [str]) -> Trie:
    """Constructs a trie based on the given list of motifs"""
    motif_trie = Trie(motifs)
    return motif_trie


def search_sequence(trie: Trie, sequence: str) -> np.array:
    """Applies motif kernel on a sequence"""
    sequence = sequence.upper()
    matching_motifs = trie._check_for_motifs(sequence)
    return matching_motifs


def motifKernel(motifs: [str], sequences: [str], return_kernel_matrix = False) -> csr_matrix:
    """
    Main function that executes the motif kernel on a set of sequences. The idea is convert strings into numerical
    values which allow for the use of machine learning algorithms.


    Args:
        motifs: A list of strings where each string is build out three possible components:
            1. Characters of a common alphabet (e.g. {A,G,C,T}) -> this alphabet should be the same for all _motifs
            2. The wildcard character "."
            3. Subsequences that contain 2 ore more characters of the alphabet (e.g. [AG] or [CT]). "^" indicates that every character
               but those in the subsequence could be part of the motif.
        sequences: A list of strings from the same aphabet as the motifs.
    Returns:
        The function returns a sparse matrix object which can be used to compute similarities between the sequences.
    """
    motif_trie = processMotifs(motifs)
    search_results = [search_sequence(motif_trie, sequence) for sequence in sequences]
    if return_kernel_matrix:
        kernel_matrix = csr_matrix(np.einsum('ij,kj->ik', search_results,search_results))
        return kernel_matrix
    else:
        return csr_matrix(search_results)



# Test
def main():
    motifs = ["A[CG]T", "C.G", "C..G.T", "G[A][AT]"
              "GT.A[CA].[CT]G"]
    sequences = ["ACGTCGATGC", "GTCGATAGC", "GCTAGCacgtaCGC",
                 "GTAGCTgtgcGTGcgt", "CGATAGCTAGTTAGC"]

    print(motifKernel(motifs, sequences))


if __name__ == '__main__':
    main()
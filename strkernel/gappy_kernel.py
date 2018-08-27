#!/usr/bin/env python3
'''
Implementation of the gappy kernel.
'''
import numpy as np
import scipy
from Bio.Seq import Seq
from scipy.sparse import csr_matrix


sequenceTypes={'dna':0,'rna':1,'aa':2,'aa+s':3}
# DNA/RNA, Amino acids (all 20), Amino acids selenocystein
alphabets=['ACGT','ACGU','ACDEFGHIKLMNPQRSTVWY','ACDEFGHIKLMNPQRSTUVWY']

def get_numbers_for_sequence(sequence,t=0,reverse=False):
    try:
        ori=[alphabets[t].index(x) for x in sequence]
    except ValueError:
        return [-1]
    if reverse:
        rev=[alphabets[t].index(x) for x in sequence.reverse_complement()]
        if ori>rev:
            return rev
    return ori


def _extract_gappy_sequence(sequence, k, g,t=0,reverse=False):
    """Compute gappypair-spectrum for a given sequence, k-mer length k and
    gap length g. A 2*k-mer with gap is saved at the same position as a 2*k-mer
    without a gap.
    The idea is to first create a vector of the given size (4**(2*k)) and then
    transform each k-mer to a sequence of length k of numbers 0-3
    (0 = A, 1 = C, 2 = G, 3 = U).
    From there, we can multiply that sequence with a vector of length 2*k,
    containing the exponents of 4 to calculate the position in the spectrum.
    Example: AUUC -> 0331 -> 4**0*1 + 4**1*3 + 4**2*3 + 4**3*0
    """
    n = len(sequence)
    kk=2*k
    alphabet=len(alphabets[t])
    powersize=np.power(alphabet, (kk))
    multiplier = np.power(alphabet, range(kk))[::-1]
    if reverse:
        powersize=int(np.power(alphabet, (kk))/2)
    spectrum = np.zeros((powersize))
    for pos in range(n - kk + 1):
            pos_in_spectrum = np.sum(multiplier * get_numbers_for_sequence(sequence[pos:pos+(kk)],t,reverse=reverse))
            spectrum[pos_in_spectrum] += 1
            for gap in range(1,g+1):
                if (pos+gap+kk)<=n:
                    pos_gap = np.sum(multiplier * get_numbers_for_sequence(sequence[pos:pos+k] + sequence[pos+k+gap:pos+gap+kk],t,reverse=reverse))
                    spectrum[pos_gap] += 1
    return spectrum

def _extract_spectrum_sequence(sequence, k,t=0,reverse=False):
    """Compute k-spectrum for a given sequence, k-mer length k.
    This method computes the spectrum for a given sequence and k-mer-length k.
    The idea is to first create a vector of the given size (4**k) and then
    transform each k-mer to a sequence of length k of numbers 0-3
    (0 = A, 1 = C, 2 = G, 3 = U).
    From there, we can multiply that sequence with a vector of length k,
    containing the exponents of 4 to calculate the position in the spectrum.
    Example: AUUC -> 0331 -> 4**0*1 + 4**1*3 + 4**2*3 + 4**3*0
    """
    n = len(sequence)
    alphabet=len(alphabets[t])
    spectrum = np.zeros(np.power(alphabet, k))
    multiplier = np.power(alphabet, range(k))[::-1]
    for pos in range(n - k + 1):
            pos_in_spectrum = np.sum(multiplier * get_numbers_for_sequence(sequence[pos:pos+k],t,reverse))
            spectrum[pos_in_spectrum] += 1
    return spectrum

def _extract_gappy_sequence_different(sequence, k, g,t=0,reverse=False):
    """Compute gappypair-spectrum for a given sequence, k-mer length k and
    gap length g. A 2*k-mer with gap is saved at the same position as a 2*k-mer
    without a gap. A 2*k-mer with a certain gap size is saved at a different
    position than the same 2*k-mer with no gaps or another number of gaps.
    """
    n = len(sequence)
    kk=2*k
    alphabet=len(alphabets[t])
    powersize=np.power(alphabet, (kk))
    multiplier = np.power(alphabet, range(kk))[::-1]
    if reverse:
        powersize=int(np.power(alphabet, (kk))/2)
    spectrum = np.zeros((g+1)*(powersize))
    for pos in range(n - kk + 1):
            pos_in_spectrum = np.sum(multiplier * get_numbers_for_sequence(sequence[pos:pos+(kk)],t,reverse=reverse))
            spectrum[pos_in_spectrum] += 1
            if (pos+g+kk+1)<n:
                for gap in range(1,g+1):
                    pos_gap = np.sum(multiplier * get_numbers_for_sequence(sequence[pos:pos+k] + sequence[pos+k+gap:pos+gap+kk],t,reverse=reverse))
                    spectrum[(gap*(powersize))+pos_gap] += 1
    return spectrum

def gappypair_kernel(sequences, k, g=0,t=0,sparse=True, reverse=False, include_flanking=False, gapDifferent = True):
    """Compute gappypair-kernel for a set of sequences using k-mer length k
    and gap size g. The result than can be used in a linear SVM or other
    classification algorithms.
    Parameters:
    ----------
    sequences:              A list of Biopython sequences
    k:                      Integer. The length of kmers to consider
    g:                      Integer. Gapps allowed. 0 by default.
    t:                      Which alphabet according to sequenceTypes.
                            Assumes Dna (t=0).
    sparse:                 Boolean. Output as sparse matrix? True by default.
    reverse:                Boolean. Reverse complement taken into account?
                            False by default.
    include_flanking:       Boolean. Include flanking regions?
                            (the lower-case letters in the sequences given)
    gapDifferent:           Boolean. If k-mers with different gaps should be
                            threated differently or all the same.
                            True by default.
    Returns:
    -------
    A numpy array of shape (N, 4**k), containing the k-spectrum for each
    sequence. N is the number of sequences and k the length of k-mers considered.
    """
    spectrum = []
    for seq in sequences:
    # To be capable to handle string input - does that make sense?
    #seq=Seq(seq)
        if include_flanking:
            seq = seq.upper()
        else:
            seq = [x for x in seq if 'A' <= x <= 'Z']
        if (g>0) and gapDifferent:
            spectrum.append(_extract_gappy_sequence_different(seq, k, g))
        elif g>0:
            spectrum.append(_extract_gappy_sequence(seq, k, g))
        else:
            spectrum.append(_extract_spectrum_sequence(seq, k))
    if sparse:
        return csr_matrix(spectrum)
    return np.array(spectrum)

#!/usr/bin/env python3
'''
Implementation of the gappy kernel.
'''
import numpy as np
import scipy
from Bio.Seq import Seq
from scipy.sparse import csr_matrix


sequenceTypes={'dna':0,'rna':1,'aa':2}
# DNA, RNA, Amino acids (all 20 + selenocystein)
s=['ACGT','ACGU','ACDEFGHIKLMNPQRSTUVWY']

def get_numbers_for_sequence(sequence,t=0,reverse=False):
    try:
        ori=[s[t].index(x) for x in sequence]
    except ValueError:
        return [-1]
    if reverse:
        rev=[s[t].index(x) for x in sequence.reverse_complement()]
        if ori>rev:
            return rev
    return ori


def _extract_gappy_sequence(sequence, k, g,t=0,reverse=False):
    """Compute k-spectrum for a given sequence, k-mer length k and gap length g.
    This method computes the spectrum for a given sequence and k-mer-length k.
    The idea is to first create a vector of the given size (4**k) and then
    transform each k-mer to a sequence of length k of numbers 0-3
    (0 = A, 1 = C, 2 = G, 3 = U).
    From there, we can multiply that sequence with a vector of length k,
    containing the exponents of 4 to calculate the position in the spectrum.
    Example: AUUC -> 0331 -> 4**0*1 + 4**1*3 + 4**2*3 + 4**3*0
    """
    n = len(sequence)
    kk=2*k
    alphabet=len(s[t])
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
    alphabet=len(s[t])
    spectrum = np.zeros(np.power(alphabet, k))
    multiplier = np.power(alphabet, range(k))[::-1]
    for pos in range(n - k + 1):
            pos_in_spectrum = np.sum(multiplier * get_numbers_for_sequence(sequence[pos:pos+k],t,reverse))
            spectrum[pos_in_spectrum] += 1
    return spectrum

def extract_spectrum(sequences, k, m=0,g=0,t=0,sparse=True, reverse=False, include_flanking=False):
    """Compute k-spectra for a set of sequences using k-mer length k.
    Computes the k-spectra for a set of sequences. This is done such that the
    resulting vectors can be fed into a linear SVM or other classification
    algorithm. The k-spectrum of a sequence is a sparse vector containing the
    number of times that a k-mer occurs in the given sequence. It has
    as many dimensions as there are possible k-mers.
    Parameters:
    ----------
    sequences:              A list of Biopython sequences
    k:                      Integer. The length of kmers to consider
    g:                      Integer. Gapps allowed. 0 by default
    include_flanking:       Include flanking regions? (the lower-case letters
                            in the sequences given)
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
        if g>0:
            spectrum.append(_extract_gappy_sequence(seq, k, g))
        else:
            spectrum.append(_extract_spectrum_sequence(seq, k))
    if sparse:
        return csr_matrix(spectrum)
    return np.array(spectrum)

class GappyPair:
    """
    This builds a wrapper for the gappy pair calculation and SVM classifier
    to support tuning the parameters via grid search.
    """
    def __init__(self, C=1, k=3, g=1,sequenceType='dna',sparse=True,reverse=False, include_flanking=False):
        self.C = C
        self.k = k
        self.g = g
        self.t=sequenceTypes[sequenceType]
        self.sparse=sparse
        self.reverse=reverse
        self.include_flanking = include_flanking
        self.clf = SVC(C=self.C, probability=True, kernel='linear')

    def fit(self, X, y):
        spectrum_X = extract_spectrum(X, self.k,g=self.g,t=self.t,sparse=self.sparse,reverse=self.reverse, include_flanking=self.include_flanking)
        self.clf.fit(spectrum_X, y)

    def predict_proba(self, X):
        spectrum_X = extract_spectrum(X, self.k,g=self.g, t=self.t,sparse=self.sparse,reverse=self.reverse,include_flanking=self.include_flanking)
        return self.clf.predict_proba(spectrum_X)

    def predict(self, X):
        spectrum_X = extract_spectrum(X, self.k,g=self.g, t=self.t,sparse=self.sparse,reverse=self.reverse,include_flanking=self.include_flanking)
        return self.clf.predict(spectrum_X)

    def set_params(self, **parameters):
        for parameter, value in parameters.items():
            setattr(self, parameter, value)
        return self

    def get_params(self,deep=True):
        return {"C" : self.C, "k" : self.k,"gap":self.g,"sequence type":sequenceTypes[self.t],"sparse":self.sparse,"reverse_complement":self.reverse,"include_flanking" : self.include_flanking}

#!/usr/bin/env python3
'''
Trie Implementation
---
Class that includes the construction of a trie based on a set of strings
and a prefix search function.
'''
import numpy as np
import scipy
import itertools
from scipy.sparse import coo_matrix
from Bio.Seq import Seq
import time
import concurrent.futures

sequenceTypes={'dna':0,'rna':1,'aa':2,'aa+s':3}
# DNA/RNA, Amino acids (all 20), Amino acids selenocystein
alphabets=['ACGT','ACGU','ACDEFGHIKLMNPQRSTVWY','ACDEFGHIKLMNPQRSTUVWY']
s=np.array([np.arange(4),np.arange(4),np.arange(20),np.arange(21)])

class TrieNode:
    """
    Implementation of a trie node.
    """
    # Char is a number
    def __init__(self, char):
        self._char = char
        self._children = []
        # list of tuples: (seq number, pos of kmer, last seen pos in kmer)
        self._q=[]

    def __repr__(self):
        return self._char

def get_sparse(g,k,t,sequences,gap_pos,gapDifferent):
    s=[[],[],[]]
    root=TrieNode('*')
    # Initialization of all possible kmers
    root._q = np.array(list(itertools.chain.from_iterable([[(i,j,0) for j in range(sequences[i].size-k+1)] for i in range(sequences.shape[0])])))
    dfs(root,s,sequences,t,k,g,0,gap_pos,gapDifferent)
    return coo_matrix((np.int_(s[0]),(np.int_(s[1]),np.int_(s[2][1:]))))

def dfs(node,sparsem,sequences,t,k,g,i,gap_pos, gapDifferent):
    """
    Depth-first-search Implementation
    """
    if i < k:
        for letter in s[t]:
            new_q=np.array([[0,0,0]])
            # At the beginning, find positions of the current letter
            if i == 0:
                new_q=np.array([(seqn,kpos,last) for (seqn,kpos,last) in node._q if sequences[seqn][kpos] == letter])
            else:
                for (seqn,kpos,last) in node._q:
                    # Only considers gaps in certain positions
                    if i in gap_pos:
                        up=update(sequences,(seqn,kpos,last),k,g,letter)
                    else:
                        up=update(sequences,(seqn,kpos,last),k,0,letter)
                    if up.shape[0] > 0 :
                        new_q=np.append(new_q,up,axis=0)
                new_q = new_q[1:]
            # If there are still possibilities, go one step deeper
            if len(new_q)>0:
                new_node = TrieNode(letter)
                new_node._q=new_q
                node._children.append(new_node)
                dfs(new_node,sparsem,sequences,t,k,g,i+1,gap_pos,gapDifferent)
    # End reached, prepare data for conversion in sparse matrix
    elif i==k:
        # sparsem =(data,i,j)
        if len(sparsem[2]) == 0:
            sparsem[2]=[-1]
        if gapDifferent:
            for gap in range(g+1):
                used = set()
                seqns=[seq for (seq,_,a) in node._q if (seq not in used)  and (a==gap+1)]
                adding=np.zeros(len(seqns))
                for qq in ((node._q)):
                    if (qq[0] in seqns) & (qq[2]==gap+1):
                        adding[seqns.index(qq[0])]+=1
                sparsem[0]=np.append(sparsem[0],adding)
                sparsem[1]=np.append(sparsem[1],np.array(seqns))
                sparsem[2]=np.append(sparsem[2],np.array([sparsem[2][-1]+1 for s in seqns ]))
        else:
            used = set()
            seqns=[seq for (seq,_,_) in node._q if seq not in used and (used.add(seq) or True)]
            adding=np.zeros(len(seqns))
            for qq in ((node._q)):
                adding[seqns.index(qq[0])]+=1
            sparsem[0]=np.append(sparsem[0],adding)
            sparsem[1]=np.append(sparsem[1],np.array(seqns))
            sparsem[2]=np.append(sparsem[2],np.array([sparsem[2][-1]+1 for s in seqns ]))

# Checks if there are matches in the kmer
def matching(sequences,seqn,kpos,last,k,g,letter):
    return np.where(sequences[seqn][kpos+last:kpos+k+g] == letter)[0]

def update(sequences,tupleq,k,g,letter):
    new=matching(sequences,tupleq[0],tupleq[1],tupleq[2],k,g,letter)
    return np.array([(tupleq[0],tupleq[1],tupleq[2]+x) for x in new if (tupleq[2]+x <= k+g) & (x>0)])

def gapkernel(sequences,k,t,g=0,gap_pos=[], gapDifferent = False):
    """Compute gapped kernel for given sequences, k-mer length k and gap length g,
    the specific type of data and the positions where gaps can occur gap_pos.
    Parameters:
    ----------
    sequences:              A numpy array of sequences
    k:                      Integer. The length of kmers to consider
    t:                      Integer. Specifies the alphabet. See sequenceTypes.
    g:                      Integer. Gaps allowed. 0 by default
    gap_pos:                Integer list. Positions, where gaps can occur.
                            If empty, all positions are considered
    Returns:
    -------
    A sparse matrix containing the k-spectrum with g-gaps for every sequence.
    """
    if not gap_pos:
        gap_pos=[i for i in range(k)]
    return get_sparse(g,k,t,sequences,gap_pos,gapDifferent)

def prepare_data(sequences, t, include_flanking = False):
    """If sequences is not a numpy array, this function can converse them to one.
    Parameters:
    ----------
    sequences:              A list of strings or Biopython sequences
    t:                      Integer. Specifies the alphabet. See sequenceTypes.
    include_flanking:       Boolean. If true, flanks are considered. False
                            by default.
    Returns:
    -------
    A sparse matrix containing the k-spectrum with g-gaps for every sequence.
    """
    if include_flanking:
        return np.array([np.array([alphabets[t].index(p.upper()) for p in x]) for x in sequences])
    return np.array([np.array([alphabets[t].index(p) for p in x if ('A' <= p <= 'Z') & (p in alphabets[t])]) for x in sequences])

def gappypair_kernel(sequences,k,t,g=1,include_flanking=False,gapDifferent = True):
    """Compute gappypair kernel for given sequences, k-mer length k and
    gap length g, the specific type of data. If sequences are not a numpy array,
    prepare data will transform them to one.
    Parameters:
    ----------
    sequences:              A numpy array of sequences or list of strings or
                            list of Biopython sequences
    k:                      Integer. The length of kmers to consider
    t:                      Integer. Specifies the alphabet. See sequenceTypes.
    g:                      Integer. Gaps allowed. 1 by default
    include_flanking:       Boolean. If true, flanks are considered. False
                            by default.
    Returns:
    -------
    A sparse matrix containing the gappypair with g-gaps for every sequence.
    """
    if (isinstance(sequences[0], str)) | (isinstance(sequences[0], Seq)):
        sequences=prepare_data(sequences, t, include_flanking)
    return gapkernel(sequences,2*k,t,g,[k],gapDifferent)

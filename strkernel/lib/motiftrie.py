#!/usr/bin/env python3
'''
Motif Trie Module
'''
# own libraries
from strkernel.lib.motif import Motif

import numpy as np

class TrieNode:
    """
    The TrieNode class consist of an element of a motif and its children in the Motif Trie.
    This class can be considerd a helper for the MotifTrie class.
    """

    def __init__(self, char: str):
        self._char = char
        self._children = []
        self._motif_finished = False
        self._motif = None

    def __repr__(self):
        return self._char


class MotifTrie:
    """
    The main objective of this class is to construct a Trie from a set of Motifs (strings). The Trie
    can then be used to compute the motif content of a sequence.

    First, Each of the motifs which is passed to the constructer is converted into a Motif object. The
    Motif objects are then used for the Trie constrcution.

    The Trie construction can be seperated into the following steps:

    1. The Trie is initialized with a root node (TrieNode object)
    2. Each Motif object is added to the Trie by parsing the Trie and adding the parts of the Motif which are not yet present.
    3. The final Node object of each Motif is marked.

    If sequences are passed to *check_for_motifs* function a DFS is performend and a numpy array containing the motif content is returned.
    """

    def __init__(self, motifs: [str]):
        self._root = TrieNode('*')
        self._motifs = motifs
        # build the trie object based on the given motifs
        for motif in motifs:
            self.add(Motif(motif))

    def check_for_motifs(self, sequence: str) -> np.array:
        """
        Iterates over the given sequence and returns the sum of the motif content of all subsequences.

        Args:
            **sequence:** A sequence (read) that only has characters that also appear in the alphabet of the motifs used to construct the MotifTrie.

        Returns:
            Numpy array containing the motif content of the sequence.
        """

        motifdict = {motif:0 for motif in self._motifs}
        for i,c in enumerate(sequence):
            motifdict = self.dfs(sequence[i:],motifdict)

        return np.fromiter(motifdict.values(), dtype=int)


    def dfs(self, sequence: str, motifdict: dict) -> [str]:
        """
        Performs a depth first search on the input sequence and adds the motif content to the input dictionary.

        Args:
            **sequence:** A part of the sequence given to check_for_motifs. In the first iteration of check_for_motifs the complete sequence is passed to this function.
            **motifdict:** A dictionary where each entry refers to one of the motifs used to construct the MotifTrie.

        Returns:
            The motifdict with the motif content in this specific sequence.
        """

        node = self._root
        char_index = 0
        priolist = []

        if not node._children:
            return False

        matching_childs = [(child, char_index)
                           for child in node._children if child._char in sequence[char_index] or child._char == "."]
        priolist.extend(matching_childs)

        while priolist:
            # get new element and index from the priorityqueue
            node, char_index = priolist.pop()
            char_index += 1
            # check if current node is the end of a motif.
            if node._motif_finished:
                motifdict[node._motif] += 1
            # check if end of the string is reached
            elif char_index == len(sequence):
                return motifdict

            # fill priorityqueue with matching childs
            matching_childs = [(child, char_index)
                               for child in node._children if sequence[char_index] in child._char or child._char == "."]
            priolist.extend(matching_childs)

        return motifdict

    def add(self, motif: Motif):
        """
        Adds a motif to the Trie.

        Args:
            **motif:** A motif object that originates from a motif (string) passed to the constructer of the MotifTrie.
        """

        node = self._root
        for char in motif:
            found_in_child = False
            # Search for matching characters in the children
            for child in node._children:
                if child._char == char:
                    # And point the node to the child that contains this char
                    node = child
                    found_in_child = True
                    break
            # If the char is not found add a new node
            if not found_in_child:
                new_node = TrieNode(char)
                node._children.append(new_node)
                # And then point node to the new child
                node = new_node
        # Mark the end of the motif and add the motif as a string.
        node._motif_finished = True
        node._motif = motif._orginal_motif

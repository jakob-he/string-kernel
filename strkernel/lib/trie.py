#!/usr/bin/env python3
'''
Trie Implementation
---
Class that includes the construction of a trie based on a set of strings
and a prefix search function.
'''
# own libraries
from strkernel.lib.motif import Motif

import numpy as np

class TrieNode:
    """Implementation of a trie node."""

    def __init__(self, char: str):
        self._char = char
        self._children = []
        self._motif_finished = False
        self._motif = None

    def __repr__(self):
        return self._char


class Trie:
    """Main Trie class."""

    def __init__(self, motifs: [str]):
        self._root = TrieNode('*')
        self._motifs = motifs
        # build the trie object based on the given motifs
        for motif in motifs:
            self._add(Motif(motif))

    def _check_for_motifs(self, sequence: str) -> np.array:
        """Returns a numpy array containing the number of motifs in the given sequence."""

        motifdict = {motif:0 for motif in self._motifs}
        for i,c in enumerate(sequence):
            motifdict = self._dfs(sequence[i:],motifdict)

        return np.fromiter(motifdict.values(), dtype=int)


    def _dfs(self, sequence: str, motifdict: dict) -> [str]:
        """Depth-first-search Implementation"""

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

    def _add(self, motif: Motif):
        """Adds a motif to the trie."""

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

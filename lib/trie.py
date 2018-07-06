#!/usr/bin/env python3
'''
Trie Implementation
---
Class that includes the construction of a trie based on a set of strings
and a prefix search function.
'''
# own libraries
from motif import Motif

import numpy as np


class TrieNode:
    """
    Implementation of a trie node.
    """

    def __init__(self, char: str):
        self._char = char
        self._children = []
        self._motif_finished = False
        self._motif = None

    def __repr__(self):
        return self._char


class Trie:
    """
    Main Trie class.
    """

    def __init__(self, motifs: [str]):
        self._root = TrieNode('*')
        self._motifs = motifs
        # build the trie object based on the given motifs
        for motif in motifs:
            self._add(Motif(motif))

    def __repr__(self):
        return self._print_node(self._root,"")

    def __contains__(self, sequence: str) -> np.array:
        """
        Returns a numpy array containing the number of motifs in the given sequence.
        """

    def _dfs(self, motif: str) -> [str]:
        """
        Depth-first-search Implementation
        """
        node = self._root
        char_index = 0
        priolist = []
        matching_motifs = []

        if not root.children:
            return False

        while priolist:
            # check if current node is the end of a motif.
            if node._motif_finished:
                matching_motifs.append(node._motif)
            # fill priorityqueue with matching childs
            matching_childs = [(child, char_index)
                               for child in node._children if child._char in motif[char_index] or child.char == "."]
            priolist.extend(matching_childs)
            node = priolist.pop()
            char_index += 1

        return matching_motifs

    def _add(self, motif: Motif):
        """
        Adds a motif to the trie.
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


def main():
    trie = Trie(["A[CG]T", "C.G", "C..G.T", "G[A][AT]",
                 "GT.A[CA].[CT]G"])
    print(trie)


if __name__ == '__main__':
    main()

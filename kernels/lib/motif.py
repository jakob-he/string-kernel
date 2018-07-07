#!/usr/bin/env python3
'''
Motif Class
---
Implementation of the motif structure described in
[Asa Ben-Hur and Douglas L. Brutlag. Remote homology detection: a motif based approach. Bioinformatics,
19:26–33, 2003.]
'''

class Motif:

    def __init__(self,motif: str):
        self._orginal_motif = motif
        self._motif = self._process_motif()
        self.index = 0

    def __repr__(self):
        return self._orginal_motif

    def __iter__(self):
        return self

    def __next__(self):
        if self.index == len(self._motif):
            raise StopIteration
        self.index += 1
        return self._motif[self.index - 1]

    def _process_motif(self):
        """
        Converts the original motif into a list of strings
        containing characters aswell as subsequences and wildcars as elements
        """
        motif = []
        brackets = False
        alphabet = self._get_alphabet()

        for char in self._orginal_motif:
            if brackets:
                if char == "]":
                    motif.append(subsequence)
                    brackets = False
                elif char in alphabet or char == ".":
                    subsequence += char

            elif char in alphabet or char ==".":
                motif.append(char)
            elif char == "[":
                brackets = True
                subsequence = ""
            else:
                print("Motif:{} does not satisfy the format requirements".format(self._orginal_motif))
        return motif


    def _get_alphabet(self):
        stripped_motif = "".join([c for c in self._orginal_motif if c not in "[]"])
        return set(stripped_motif)

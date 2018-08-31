#!/usr/bin/env python3
"""
Motif Module.
"""


class Motif:
    """
    The input string gets transformed into an object of the Motif class. The main attribute of this class is the motif itself which is given
    as a list of strings where each string is one of the following elements:

    1. A character of a common alphabet (e.g. {A,G,C,T}) -> this alphabet should be the same for all motifs
    2. The wildcard character "."
    3. Substiution groups that contain 2 ore more characters of the alphabet (e.g. [AG] or [CT]). "^" as a leading character
       indicates that every character in the alphabet but those in the substitution group matches this part of the motif.
    """

    def __init__(self, motif: str):
        self._orginal_motif = motif
        self._motif = self.process_motif()
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

    def process_motif(self):
        """Processes the input motif which is a string and returns a list of strings."""
        motif = []
        brackets = False
        complement = False
        alphabet = self.get_alphabet()

        for char in self._orginal_motif:
            if brackets:
                if char == "^":
                    complement = True
                if char == "]":
                    if complement:
                        motif.append(alphabet - set(subsequence))
                        brackets = False
                    else:
                        motif.append(subsequence)
                        brackets = False
                elif char in alphabet or char == ".":
                    subsequence += char

            elif char in alphabet or char == ".":
                motif.append(char)
            elif char == "[":
                brackets = True
                subsequence = ""
            else:
                print("Motif:{} does not satisfy the format requirements".format(
                    self._orginal_motif))
        return motif

    def get_alphabet(self):
        """Extracts the alphabet (unique characters) from a string."""
        stripped_motif = "".join(
            [c for c in self._orginal_motif if c not in "[]"])
        return set(stripped_motif)

"""
 Module: Trie
 Python implementation of trie data structure for mismatch string kernel
 Author: Meng Zhang <RainaMeng@outlook.com>
 Based on: Dohmatob Elvis Dopgima's Trie module
 <https://github.com/dohmatob/kernels/blob/master/python/trie.py>
"""

import numpy as np


class MismatchTrie(object):
    """
    Trie implementation, specific to 'Mismatch String Kernels'.
    """

    def __init__(self, label=None, parent=None):
        """
        label: int, optional (default None), node label.
        parent: `Trie` instance, optional (default None), node's parent.
        """

        self.label = label  # label on edge connecting this node to its parent
        self.level = 0  # level of this node beyond the root node
        self.children = {}  # children of this node

        # concatenation of all labels of nodes from root node to this node
        self.full_label = ""
        # for each sample string, this dict holds pointers to it's k-mer substrings
        self.kmers = {}

        self.parent = parent

        if not parent is None:
            parent.add_child(self)

    def is_root(self):
        """
        Check whether this node is the root.
        """

        return self.parent is None

    def is_leaf(self):
        """
        Check whether this node is a leaf.
        """

        return len(self.children) == 0

    def is_empty(self):
        """
        Check whether a node has 'died'.
        """

        return len(self.kmers) == 0

    def copy_kmers(self):
        """
        Copy the kmer data for this node (not the reference pointer).
        """

        return {index: np.array(substring_pointers)
                for index, substring_pointers in self.kmers.items()}

    def add_child(self, child):
        """
        Add a new child to this node.
        """

        assert not child.label in self.children

        # initialize kmers data to that of parent
        child.kmers = self.copy_kmers()

        # child is one level beyond parent
        child.level = self.level + 1

        # parent's full label (concatenation of labels on edges leading from root node)
        # is a prefix to child's the remainder is one symbol, the child's label
        child.full_label = '%s%s' % (self.full_label, child.label)

        # let parent adopt child: commit child to parent's booklist
        self.children[child.label] = child

        # let child adopt parent
        child.parent = self

    def delete_child(self, child):
        """
        Delete a child.
        """

        # get child label
        label = child.label if isinstance(child, MismatchTrie) else child

        # check that child really exists
        assert label in self.children, "No child with label %s exists." % label

        # delete the child
        del self.children[label]


    def compute_kmers(self, training_data, k):
        """
        Compute the metadata for this node, i.e, for each input string
        training_data[index], compute the list of offsets of it's k-mers
        together with the mismatch counts (intialially zero) for this k-mers
        with the k-mer represented by this node `self`.

        Parameters
        ----------
        training_data: 2D array of shape (n_samples, n_features)
                       training data for the kernel.
        k: int, used in k-mers to compute the kernel.
        """

        # sanity checks
        if not isinstance(training_data, np.ndarray):
            training_data = np.array(training_data)

        if training_data.ndim == 1:
            training_data = np.array([training_data])

        assert training_data.ndim == 2

        # compute the len(training_data[index]) - k + 1 kmers of each input training string
        for index in range(len(training_data)):
            self.kmers[index] = np.array([(offset,
                                           0 # no mismatch yet
                                           )
                                            for offset in range(len(training_data[index])-k+1)])

    def process_node(self, training_data, k, m):
        """
        Process this node. Recompute its supported k-mers.
        Finally, determine if node survives.

        Parameters
        ----------
        training_data: 2D array of shape (n_samples, n_features)
                       training data for the kernel
        k: int, used in k-mers to compute the kernel
        m: int
           maximum number of mismatches for 2 k-mers to be considered 'similar'
           Normally small values of m should work well
           Plus, the complexity the algorithm is exponential in m
        
        Return
        -------
        True if node survives, False else
        """

        # sanity checks
        if not isinstance(training_data, np.ndarray):
            training_data = np.array(training_data)
        if training_data.ndim == 1:
            training_data = np.array([training_data])

        assert training_data.ndim == 2

        if self.is_root():
            # compute meta-data
            self.compute_kmers(training_data, k)
        else:
            # loop on all k-mers of input string training_data[index]
            for index, substring_pointers in self.kmers.items():
                # update mismatch counts
                substring_pointers[..., 1] += (training_data[index][
                        substring_pointers[..., 0] + self.level - 1
                        ] != self.label)

                # delete substring_pointers that present more than m mismatches
                self.kmers[index] = np.delete(substring_pointers,
                                               np.nonzero(substring_pointers[..., 1] > m),
                                               axis=0)

            # delete entries with empty substring_pointer list
            self.kmers = {index: substring_pointers for (
                    index, substring_pointers) in self.kmers.items(
                    ) if len(substring_pointers)}

        return not self.is_empty()

    def update_kernel(self, kernel):
        """
        Update the kernel in memory.

        Parameters
        ----------
        kernel: 2D array of shape (n_samples, n_samples)
                kernel to be updated
        full_label: mismatch kmers generated by traversing labels from root to leaf
        """

        for i in self.kmers:
            for j in self.kmers:
                kernel[i, j] += len(self.kmers[i]) * len(self.kmers[j])


    def traverse(self, training_data, l, k, m, kernel=None,
                 kernel_update_callback=None):
        """
        Traverses a node, expanding it to plausible descendants.

        Parameters
        ----------
        training_data: 2D array of shape (n_samples, n_features)
                       training data for the kernel
        l: int, size of alphabet
           Examples of values with a natural interpretation:
           2: for binary data
           256: for data encoded as strings of bytes
           4: for DNA/RNA sequence data (bioinformatics)
           20: for protein data (bioinformatics)
        k: int, we will use k-mers to compute the kernel
        m: int
           maximum number of mismatches for 2 k-mers to be considered 'similar'
           Normally small values of m should work well
           Plus, the complexity the algorithm is exponential in m
        kernel: 2D array of shape (n_samples, n_samples)
                optional (default None) kernel to be, or being, estimated

        Returns
        -------
        kernel: 2D array of shape (n_samples, n_samples), estimated kernel
        n_survived_kmers: int, number of leaf nodes that survived the traversal
        go_ahead: boolean, a flag indicating whether the node got aborted (False) or not
        """

        # initialize kernel if None
        if kernel is None:
            kernel = np.zeros((len(training_data), len(training_data)))

        # counts the number of leafs which are decendants of this node
        n_surviving_kmers = 0

        # process the node
        go_ahead = self.process_node(training_data, k, m)

        # if the node survived
        if go_ahead:
            # we've hit a leaf
            if k == 0:
                # yes, this is one more leaf/kmer
                n_surviving_kmers += 1

                # update the kernel
                self.update_kernel(kernel)

            else:
                # recursively bear and traverse child nodes
                for j in range(l):
                    # bear child
                    child = MismatchTrie(label=j, parent=self)

                    # traverse child
                    kernel, child_n_surviving_kmers, \
                        child_go_ahead = child.traverse(
                        training_data, l, k - 1, m, kernel=kernel)

                    # delete child if dead
                    if child.is_empty():
                        self.delete_child(child)

                    # update leaf counts
                    n_surviving_kmers += child_n_surviving_kmers if \
                        child_go_ahead else 0

        return kernel, n_surviving_kmers, go_ahead


    def __iter__(self):
        """
        Return an iterator on the nodes of the trie
        """
        
        yield self
        for child in self.children.values():
            for grandchild in child:
                yield grandchild

    def leafs(self):
        for leaf in self:
            if leaf.is_leaf():
                yield leaf

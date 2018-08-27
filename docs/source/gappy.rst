Gappy Kernel
============

Bio sequences such as DNA, RNA and protein sequences are widely used in SVMs, because they have some certain differences to normal text strings (e.g. a smaller alphabet) kernels fitted to their need were invented. One of them is the gappy-pair kernel [1]_.
Due to mutations and reading errors in data analysis two genome sequences can be similar even if they do not match exactly but with some gaps. When using the spectrum kernel such similarities will be lost between sequences, therefore the gappy-pair kernel was introduced. Instead of only recognising all k-mers in a sequence for a given k, the gappy-pair kernel also counts how many k-mers with a certain number of gaps(g) appear in the sequence. It does so by redefining k (called k’ from here on), which is from now on not referring to the whole length of a k-mer but to half of it, so that the gappy-pair kernel recognises every k-mer of the sort where the first k’ bases of the k-mer match the sequence, then a number of gaps from 0 to g is allowed and then the last k’ bases match the sequence again. (Therefore, in every k-mer 2k’ bases are present.)
For example for the sequence “AACG” the gappy-pair kernel for k’=1 and g=2 would count the following k-mers once: [“AA”, “A.C”, “A..G”, “AC”, “A.G”, “CG”]. Because it can make sense biologically to not distinguish between a k-mer without gaps and the same k-mer with gaps, our kernel allows the user to decide how to proceed by setting the Boolean variable gapDifferent (Default is set to true).
A simple approach of calculating the gappy-pair kernel is to go through all given sequences singularly and then scan through the sequence counting the seen k-mers. To do so it is necessary to create a matrix of the size of all possible k-mers. The number of possible k-mers with gapDifferent being true is 〖(g+1)n〗^k with n being the length of the alphabet and with gapDifferent being false n^k. Therefore, even if the output is translated into a sparse matrix, during the calculation a significant amount of memory is occupied which can cause problems.
A more memory-sufficient approach was implemented as well by using a trie structure based on the trie for gapped kernels [2]_. Instead of evaluating sequences one by one, all sequences are considered while moving along the trie with a depth-first-search.

How to use the Gappy-Pair Kernel
---------------------------

Define a k' and a g and use either the function gappypair_kernel from gappy_kernel for the simple approach or gappy_trie for the trie approach.
Both function return a sparse matrix, which can than be used in a SVM classifier.

    from gappy_kernel import gappypair_kernel as gk
    from gappy_trie import gappypair_kernel as gt

    sequences = ['ACGTAACT','GTATATAG', 'TAGAGATT']
    # Note the t in the input parameters is refering to dna.
    # Alternatively, one can choose 1 for rna, 2 for amino acids and 3 for amino acids + selenocystein
    # simple approach
    X1 = gk(sequences,k=1,t=0,g=1)
    # trie approach
    X2 = gt(sequences,k=1,t=0,g=1)

    # For usage in a SVM
    from scipy.sparse import coo_matrix, vstack
    clf = SVC(C=0.1, kernel='linear', probability=True)
    clf.fit(X1, y)
    # Or for the trie approach
    clf = SVC(C=0.1, kernel='linear', probability=True)
    clf.fit(X2, y)

Please look into the ipython tutorial Gappy_Kernel_IGF_Example.ipynb.


References
----------

.. [1] Pavel Kuksa, Pai-Hsi Huang, and Vladimir Pavlovic. A fast, large-scale learning method for protein sequence classification. In 8th Int. Workshop on Data Mining in Bioinformatics, pages 29–37, Las Vegas, NV, 2008.
.. [2] Leslie C., Kuang R., Fast String Kernels using Inexact Matching for Protein Sequences, Journal of Machine Learning Research, 2004, 1435-1455.

Modules
-------

.. automodule:: strkernel.gappy_kernel
  :members:
  :show-inheritance:

.. automodule:: strkernel.gappy_kernel
  :members:
  :show-inheritance:

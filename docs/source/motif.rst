Motif Kernel
============

The Motif Kernel depends on a user defined set of motifs. For each input sequence the number of contained motifs (**motif content**) is computed. These values can then be used as a similarity measurement. Since the features of this kernel are restricted to a user defined set, the performance in later analysis of the kernel output (e.g. classification of cell populations) depends highly on the picked motifs. While this offers a lot of flexibility, the user should carefully decide which motifs are used. The package does not yet include any method to extract significant motifs from a set of sequences. There is, however, a large variety of motif extraction applications available [1]_ [2]_ [3]_.  Information on how exactly this kernel implementation works can be found in Modules. 

A single motif is defined by a sequence of the following elements:

- A single character from an alphabet e.g. [A (Adenine), G (Guanine), T (Thymine), C (Cytosine)]. This character only matches the exact same character. 

- The wildcard character **"."** which can represent and therefore matches any of the characters from the alphabet.

- A substitution group, which is a list of characters from the alphabet enclosed in square brackets e.g. **[AG]**. This group matches every character within the brackets (e.g. A or G). If the leading character is a **"^"** the substitution group matches any character **but** those in the brackets (e.g. C or T).  

How to use the Motif Kernel
---------------------------

The collection of motifs has to be defined as a list of strings::

    motif_collection = ["A[^CG]T", "C.G", "C..G.T", "G[A][AT]", "GT.A[CA].[CT]G"]

This collection can then be used to create the motif Kernel::

    from strkernel.motifkernel import motifKernel
    motif_kernel = motifKernel(motif_collection)

Afterwards a set of sequences can be passed to the motifKernel object which returns a sparse matrix. The sparse matrix contains the motif content for each sequence (default) or the similarities between the sequences based on the motif content (with *return_kernel_matrix* enabled). If the flanking characters (lower case letters) should not be included *include_flanking* has to be disabled. Flanking characters are enabled by default since it generally improves accuracy if the kernel output is used for classification. For large sets of sequences disabling *include_flanking* characters can give a significant performance boost::

    sequences = ["ACGTCGATGC", "GTCGATAGC", "GCTAGCacgtaCGC","GTAGCTgtgcGTGcgt", "CGATAGCTAGTTAGC"] 
    #compute motif content matrix
    motif_content = motif_kernel.compute_matrix(sequences)
    #compute kernel matrix 
    kernel_matrix = motif_kernel.compute_matrix(sequences, return_kernel_matrix = True)
    #compute motif content matrix without considering the flanking characters
    matirx_without_flanking = motif_kernel.compute_matrix(sequences, include_flanking = False)

References
----------

.. [1] Pissis, S.P., Stamatakis, A. and Pavlidis, P., 2013, September. MoTeX: A word-based HPC tool for MoTif eXtraction. In Proceedings of the International Conference on Bioinformatics, Computational Biology and Biomedical Informatics (p. 13). ACM. 
.. [2] Pissis, S.P., 2014. MoTeX-II: structured MoTif eXtraction from large-scale datasets. BMC bioinformatics, 15(1), p.235.
.. [3] Zhang, Y. and Zaki, M.J., 2006. EXMOTIF: efficient structured motif extraction. Algorithms for Molecular Biology, 1(1), p.21.

Modules
-------

.. automodule:: strkernel.motifkernel
  :members:
  :show-inheritance:

Motif Trie
~~~~~~~~~~

.. automodule:: strkernel.lib.motiftrie
  :members:
  :show-inheritance:

Motif
~~~~~

.. automodule:: strkernel.lib.motif
  :members:
  :show-inheritance:

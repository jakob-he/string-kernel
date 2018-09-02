Mismatch Kernel
===============

The mismatch kernel can be used with support vector machines (SVMs) in a discriminative approach to a classification problem. It measures the sequence similarity based on shared occurrences of **k**-length subsequences, counted with up to **m** mismatches, which is what we generally state as **(k, m)-mismatch**. The kernels can be efficiently computed by using a mismatch tree data structure [1]_.

The user input consists of the following elements:

- data: 2D matrix of shape (n_samples, n_features), sequences needed to be computed as kernel

- l: int, size of alphabet. Examples of values with a natural interpretation:
     
     2: for binary data
     
     256: for data encoded as strings of bytes
     
     4: for DNA/RNA sequence data (bioinformatics)
     
     20: for protein data (bioinformatics)
     
- k: int, used in k-mers to compute the kernel
     
- m: int, maximum number of mismatches for 2 k-mers to be considered 'similar'. Normally, small values of m should work well. The complexity of the algorithm is exponential in m

The expected output is the mismatch string kernel:
     
- kernel: 2D matrix of shape (n_samples, n_samples) suggesting the similarities between sequences


How to use the Mismatch Kernel
------------------------------

The collection of computed sequences has to be defined as a list of strings, for example::

    mismatch_collection = ['aATGCg', 'ACGTTT', 'agATGC','TCACcg', 'cgTCTCGAgt']

This collection can then be used to compute the mismatch kernel::

    from strkernel.mismatch_kernel import MismatchKernel
    mismatch_kernel = MismatchKernel(l=l, k=k, m=m).get_kernel(mismatch_collection)

The sequence collection has to be preprocessed to regulate the input format and avoid subsequent errors before being used to compute the mismatch kernel, which can be achieved by calling the function *preprocess* in *mismatch_kernel*::

    from strkernel.mismatch_kernel import preprocess
    after_process = preprocess(mismatch_collection)
    mismatch_kernel = MismatchKernel(l=l, k=k, m=m).get_kernel(after_process)

The result is a sparse matrix indicating the similarities between different sequences based on the inner product of occurence counts of all (k, m)-mismatch k-mers. If the lower case letters in every string are being considered, the parameter *ignoreLower* in *preprocess* function should be set as **False**, which is **True** by default. In addition, the function *get_kernel* can also be tuned by the parameter *normalize*. Normalization is enabled by default since it generally improves accuracy. Normalization is realized by: **kernel[x, y] / sqrt(kernel[x, x] * kernel[y, y])**. That is, the diagonal elements in the kernel matrix are set to 1 since the similarity between a certain string and itself is the largest::

    sequences = ['aatgcACGTTGAgatcg','acgtgACGTTTGacggt', 'agtccATGCTGTaagtc','gttccTCACCGTcgcgt', 'gtacgTCTCGCTgtcgt']
    # preprocess
    after_process = preprocess(mismatch_collection, ignoreLower=False)
    # compute mismatch kernel
    mismatch_kernel = MismatchKernel(l=l, k=k, m=m).get_kernel(after_process, normalize = False)

We also provide a function to display the middle step of getting a kernel allowing the users to check the details before computing the final kernel. *leaf_kmers* shows us all mismatch kmers and the occurence counts of every (k, m)-mismatch k-mer in every string. n vectors of length m are the output of *leaf_kmers*, where n is number of strings, whose similarities will be computed in the susequent steps and m is the number of all (k. m)-mismatch k-mers. The similarity between string i and string j is the inner product of vector i and j, where 1<= i, j <= n::

    print(mismatch_kernel.leaf_kmers)

References
----------

.. [1] Leslie C.S., Eskin E., Cohen A., Weston J., Noble W.S. Mismatch string kernels for discriminative protein classification. Bioinformatics. 2004;20:467–476. doi: 10.1093/bioinformatics/btg431.

Modules
-------

.. automodule:: strkernel.mismatch_kernel
  :members:
  :show-inheritance:

Mismatch Trie
~~~~~~~~~~~~~

.. automodule:: strkernel.lib.MismatchTrie
  :members:
  :show-inheritance:

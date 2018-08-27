Getting Started
===============

Introduction
------------

Strkernel is a python package designed to perform a kernel based analysis of biological sequences. The implementation assumes the use of Support Vector Machines (SVMS) but does not strictly require it since each kernel can be used separately of any machine learning algorithm. Further instructions on how to combine the kernel methods with machine learning approaches can be found in the Tutorials_ section. The package provides three different kernels which can be used for sequence analysis:

- Gappy Kernel [1]_
- Mismatch Kernel [2]_
- Motif Kernel [3]_

The general use case of strkernel is the conversion of a set of sequences into a matrix with numeric similarity measurements. The expected input are either sequences or Sequence objects of the Biopython_ package. 
In addition to the default output, the Motif Kernel implementation supports the output of a similarity matrix for a set of given sequences. The default output, however, is a sparse matrix where each row represents one sequence. The meaning of the columns depends on which kernel was used and is be explained in the corresponding section_. 

.. _R: https://bioconductor.org/packages/release/bioc/html/kebabs.html
.. _Tutorials: examples.html#Tutorials
.. _Biopython: https://biopython.org/
.. _section: kernels.html#Kernels

Installation
------------

- Install strkernel via PyPi (recommended)::

     sudo pip install strkernel

- Install strkernel via github::

     git clone https://github.com/jakob-he/string-kernel
     python setup.py install

Tests
-----

You can run several tests to ensure that the package is fully working. To use the testsuite you have to download the package via github::

    git clone https://github.com/jakob-he/string-kernel
    python setup.py test


References
----------

.. [1] Pavel Kuksa, Pai-Hsi Huang, and Vladimir Pavlovic. A fast, large-scale learning method for protein sequence classification. In 8th Int. Workshop on Data Mining in Bioinformatics, pages 29–37, Las Vegas, NV, 2008.
.. [2] Christina S. Leslie, Eleazar Eskin, Adiel Cohen, Jason Weston, and William Stafford Noble. Mismatch string kernels for discriminative protein classification. Bioinformatics, 1(1):1–10, 2003.
.. [3] Asa Ben-Hur and Douglas L. Brutlag. Remote homology detection: a motif based approach. Bioinformatics, 19:26–33, 2003.


What is a String-Kernel?
-----------------
To use SVMs with strings it is necessary to have a measurement of similarity for different strings.
The idea of the most basic string kernel, the spectrum kernel is to count the appearance of k-mers in the string. As a result a string can be represented as a numerical sequence and the similarity of two strings can be calculated by using different kernels, for instance the dot product (linear kernel) of these numerical representations can be calculated.
Bio sequences such as DNA, RNA and protein sequences are widely used in SVMs, because they have some certain differences to normal text strings (e.g. a smaller alphabet) kernels fitted to their need were invented. Suchs as the gappy-pair, mismatch and motif kernel. 

Use cases
---------

Getting Started
===============

Introduction
------------
In order, to use support vector machines (SVMs) on strings a measure for similarity is necessary.
Therefore, different kernels have been developed (e.g. spectrum kernel) but while in R the use of string kernels is very simple thanks to the kebab package, a similar user friendly, easy-to-use package with different string kernels is absent in Python.
We implemented here a python3 package which is easy to install with pip and includes not only the spectrum kernel but also the gappy-pair, mismatch and motif kernel.


Installation
------------

What is a String-Kernel?
-----------------
To use SVMs with strings it is necessary to have a measurement of similarity for different strings.
The idea of the most basic string kernel, the spectrum kernel is to count the appearance of k-mers in the string. As a result a string can be represented as a numerical sequence and the similarity of two strings can be calculated by using different kernels, for instance the dot product (linear kernel) of these numerical representations can be calculated.
Bio sequences such as DNA, RNA and protein sequences are widely used in SVMs, because they have some certain differences to normal text strings (e.g. a smaller alphabet) kernels fitted to their need were invented. Suchs as the gappy-pair, mismatch and motif kernel. 

Use cases
---------

# Java_MCS_algorithms
===========
PhD project work with Sheffield University.  Variety of MCS algorithms coded in Java

The dataset
===========

Unzip the data
--------------
The data has been compressed with 7zip. It can be unzipped with p7zip on Linux, or 7zip on Windows.

Description
-----------
The benchmark data will be placed in the dataset subdirectories of the
SingleAssay and MultiAssay directories. There are 1000 files corresponding
to the 1000 repetitions. Each file contains several thousand lines of
CHEMBL IDs, where the first ID is the reference molecule, and the other
four are molecules are increasing distance (decreasing similarity) to the
reference.

How to reproduce the results
============================

Requirements
------------
1. Python 2.7
2. NumPy
3. SciPy
4. RDKit (2015.09.2)

Optional but needed to generate the graph depictions
----------------------------------------------------
1. dot (provided by GraphViz)

Get ChEMBL
----------

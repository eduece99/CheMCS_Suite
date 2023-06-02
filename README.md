MCS Suite
===========

Description
-----------
PhD project work with Sheffield University (Under supervision of Professor Peter Willett and Doctor John Holliday).  Variety of Maximum Common Substructure (MCS) algorithms for use with chemical graphs/structures.  Originally coded in Java 1.6.



Compilation
-----------
Please compile using Eclipse IDE (I used Mars.2).  Usage examples revolve around the ```ExtendedIsomorphism``` class, in ```MCSMethodsTest.java```

More recommended however is to build the project using Maven.  It is best to do this in order to create a .jar file for command line usage.  

<br />
<br />


Command Line Usage
-----------

Navigate to the directory where the jar file in question has been created (via Maven)

```
java -jar Java_MCS_algorithms-1.0.0-jar-with-dependencies.jar  ../data/input/chembl751606.smi ../data/input/aid466_decoys_5.sdf
```


**Citing** 
<br />
Please cite the article below when using this code:

Duesbury et al, 2017, Comparison of Maximum Common Subgraph Isomorphism Algorithms for the Alignment of 2D Chemical Structures, ChemMedChem.  Available at [http://onlinelibrary.wiley.com/doi/10.1002/cmdc.201700482/abstract]



Java Library Requirements
------------
1. Chemistry Development Kit (CDK) (2.8) - molecule handling
2. Ambit Core and SMARTS libraries (2.4.13) - fast SMARTS parsing and generation used in MCS representation
3. Colt (1.2) - Eigenvector calculations
4. picocli (4.7.3) - command line options
5. guava (31.1-jre) - alternative list handling
6. SMSD (core and exec) (2015 04 01) - additional graph theory tools and MCS handling
7. JUnit 5 - tests (included in Eclipse 2.9.3)


Author Information
----------------------------------------------------
Dr Edmund Duesbury (eduece99 AT g mail DOT com)

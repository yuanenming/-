# align-free-phylogentic-tree
Phylogenetic Tree Construction using Alignment Free methods.

This is a coursework of algorithm and data structure.

A simple framework written in C++ language.

Compilation: Under Linux environment, you just need to run
    $ g++ *.cpp
    then use the a.out project
    ***

Usage: ./a.out [filename]
	where filename is the file containing the organisms' nuclortide sequence which will be used for constructing the Tree. Possible file extensions: .fasta and .paml


Three options for constructing a distance matrix:
	
	Option 1: Composition Vectors (available for 3 to 10 mers) 
 	
 	Option 2: Word Context (available for 5 mers)
 	
 	Option 3: Spaced Word Frequency (available for 5 mers)
  
This project will output a result in Netwick Tree Format.

If you have installed R with ggplot2 package, you can get a phylogentic tree image in PDF format.

# phylogenetic-tree-construction
Phylogenetic Tree Construction using Alignment Free methods.

A simple framework written in C++ language.

Compilation: In the Project's folder execute make -f Makefile.win all.

Usage: ./PhyloTree [filename]
	where filename is the file containing the organisms which will be used for constructing the Tree.

Possible file extensions: .fasta and .paml

Fasta Format 

	>"Sequence/Organism Name" (e.g. AAAAAAA or GOD_Usopp)
	"The Sequence" (i.e. "ATCGTGNTACT...")

Paml Format 

		NumberOfOrganisms LengthOfSequences (These are ignored!)
      
	"Sequence/Organism Name" (e.g. AAAAAAA or GOD_Usopp)        "The Sequence" (i.e. "ATCGTGNTACT...")

Three Possible algorithms:
	
	Option 1: Composition Vectors (available for 3 to 10 mers) 
 	
 	Option 2: Word Context (available for 5 mers)
 	
 	Option 3: Spaced Word Frequency (available for 5 mers)
  
Export resulting Phylogenetic Tree in Netwick Tree Format (https://en.wikipedia.org/wiki/Newick_format) in a file given by the user.

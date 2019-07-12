readme.txt: This file explains the code and the structure of the algorithm.

Version 2.0


I. Motivation

This simulator is able to simulate a P system with active membranes, with 3 polarizations and elementary division rules. The input of the algorithm is given by a binary file, which codifies the initial configuration of the system and the rules. This version is based on the sequential algorithm, designed by I. Pérez-Hurtado et al., in the Java simulation library pLinguaCore. It is based on two phases:

1. Selection of rules to be executed.
2. Execution of rules selected in 1.


II. Installation: 

1. Install the CUDA SDK version 4.X.

2. Install the counterslib: extract the file counterslib.tar.gz inside the common folder of the CUDA SDK.

3. Copy the folder into the source folder of the CUDA SDK, and type "make".


III. Usage:

Type ./pcuda -h to list the different options.
* A sequential simulation: ./pcuda -s -i file.bin
* A fast sequential simulation: ./pcuda -f -i file.bin -m #membranes -o #objects -b #threadPerBlock -t 1
* A parallel simulation on the GPU: ./pcuda -p -i file.bin -m #membranes -o #objects -b #threadPerBlock -t 1


IV. Source:

The objective of each file is the following:

- main.cpp: Contains the main function, which reads the input parameters of the simulator, calls the parser and call the sequential algorithm.

- parserbin.h, parserbin.cpp: Contains two different classes:
  * Parserbin: Reads and returns the contained information in a binary input file.
  * Configuration: Store the configuration of a P systems (administer the set of rules, membranes, alphabet, labels, selected rules, membrane identifiers, etc.). When reading an input file, Parserbin returns a Configuration object with all the information.

- binvalues.h: Defines macros configuring aspects such as fields lengths in the binary file, etc.

- membrane.h: Defines the structure that a membrane has (first-child, next-brother structure).

- multiset.h, multiset.cpp: Implements a multiset by a linked list of objects (if an object has multiplicity 0, it is not in the list). Contains useful methods for that purpose, and one to returns, in an array of multiplicities form, the content of the multiset.

- ruleset.h, ruleset.cpp: Defines useful classes to store and loop the rules:
	* Ruleset: The rule store. Stores the rules in an ordered way: First by label, second by charge, and 
	lastly, a list of ordered rules by type (first evolution, then the rest).
	* Rulelist: An iterative class of rules. Instead of returning the list of rules associated to each label and charge, 
	it returns this class that has methods to initiate the iteration.
	* Selectedrulelist: Defines the selected rules to be used in the execution phase.
	Each node stores the selected rules, the membrane ID and the number of times. Moreover, it contains method to iteratively loop
	the list.
	
- sequential.h, sequential.cpp: Define the functions implementing the simulation algorithm (sequential_solve)
and their two phases: selection (select_rules) and execution (execute_rules).

- fast_sequential.h, fast_sequential.cpp: Define the functions implementing the fast sequential simulator.

- pcuda.cu, pcuda.h, pcudadef.cpp: Defines all the functionality implemented on CUDA. pcuda.cu contains the call to kernels and the main
	loop of the algorithm.

- seq2pcuda.cpp, seq2pcuda.h: Defines how to translate the data structures from the sequential simulator for the pcuda and fast sequential
	simulators.

- selection_kernel.cu, selection_kernel_only.cu (deprecated), dissolution_execution_kernel.cu, division_execution_kernel.cu, 		
	evolution_execution_kernel.cu (deprecated), sendin_execution_kernel.cu, sendout_execution_kernel.cu: Defines all the kernels used to 
	run the	simulator on the GPU.


/*$Id: readme.txt 2009-01-11 15:14:44 mdelamor $*/

/*
    Pcuda: Simulating P systems with active membranes on the GPU 
    This simulator is published on:
    J.M. Cecilia, J.M. García, G.D. Guerrero, M.A. Martínez-del-Amor, I. Pérez-Hurtado,
    M.J. Pérez-Jiménez. Simulation of P systems with active membranes on CUDA,
    Briefings in Bioinformatics, 11, 3 (2010), 313-322

    Pcuda is a subproject of PMCGPU (Parallel simulators for Membrane 
                                       Computing on the GPU)   
 
    Copyright (c) 2009 Miguel Á. Martínez-del-Amor (RGNC, University of Seville)
 		       Ginés D. Guerrero (GACOP, University of Murcia)
		       Chema Cecilia (GACOP, University of Murcia)
		       Ignacio Pérez-Hurtado (RGNC, University of Seville)
    
    This file is part of Pcuda.
  
    Pcuda is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Pcuda is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Pcuda.  If not, see <http://www.gnu.org/licenses/>. */

/**
  Definition file for the fields of the binary file
*/
/**********************************************/
/* Types */
#ifndef BINVALUES_H_
#define BINVALUES_H_

typedef unsigned int MembraneID;
typedef unsigned short int LabelID;
typedef unsigned int ObjectID;
typedef unsigned int RuleID;

/**********************************************/
/* Values */

#define HEADER 0xAF12FA00
#define HEADER_BITMASK 0xFFFFFF00
#define HEADER_VERSION_BITMASK 0x0000000F
#define HEADER_VARIANT_BITMASK 0x000000F0
#define END_OF_TEXT '\0'

/**********************************************/
/* Sizes */

#define MAX_FIELD_SIZE 4
#define MAX_TEXT_SIZE 256

/*   About header */
#define HEADER_SIZE 4

/*   About the multiset of objects */
#define NUMBER_INITIAL_MULTISETS_SIZE 2
#define NUMBER_DIFERENT_OBJECTS_SIZE 2
#define OBJECT_ID_SIZE_LONG 2
#define OBJECT_ID_SIZE_SHORT 1
#define OBJECT_MULTIPLICITY_SIZE 2

/*   About membranes */
#define NUMBER_LABELS_SIZE 2
#define NUMBER_MEMBRANES_SIZE 2
#define MEMBRANE_ID_SIZE_LONG 2
#define MEMBRANE_ID_SIZE_SHORT 1
#define MEMBRANE_ID_LABEL_SIZE_LONG 2
#define MEMBRANE_ID_LABEL_SIZE_SHORT 1
#define MEMBRANE_CHARGE_SIZE 1
#define CHARGE_BIT_MASK 0x03

#define NUMBER_RULES_SIZE 2

#endif /* BINVALUES_H_ */

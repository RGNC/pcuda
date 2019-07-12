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


/*
 * seq2cuda.h
 *
 *  Created on: 24-feb-2009
 *      Author: miguel
 */

#ifndef SEQ2PCUDA_H_
#define SEQ2PCUDA_H_
#include <stdlib.h>
#include <iostream>
#include "pcudadef.h"
#include "ruleset.h"
#include "multiset.h"
#include "membrane.h"

size_t get_max_memory(int num_objects, int num_labels, int max_membranes);

/* Creates the "fishing rod" rule set from the rule set of
 * the sequential version */
//Rodruleset get_rodruleset(Ruleset& ruleset, uint & size, Rodmultiset* evol_table, int alphabet_size, int label_set_size);
Rodruleset get_rodruleset(Ruleset& ruleset, uint & rod_length, uint & evol_length, Rodmultiset * evol_table, int alphabet_size, int label_set_size);

/* Creates the membrane set from the set of membranes of the sequential version */
Membraneset get_membraneset(Membrane * skin, uint & size, uint & num_membranes, uint max_membranes, int num_levels);

/* Initalizes the structure for the input multisets */
ushort * get_multiset_set(Membrane *skin, uint & size, uint max_membranes, int num_levels, int num_objects);

uint * get_skinmultiset(Membrane * skin, uint & size, int num_objects);

Selec_siodd ini_selecsiodd(uint & size, uint max_membranes);

Multiselec ini_selecevol(uint & size, uint max_membranes, uint num_objects);

#endif /* SEQ2PCUDA_H_ */

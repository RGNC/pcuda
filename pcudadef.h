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
 * pcudadef.h
 *
 *  Created on: 24-feb-2009
 *  Author: miguel, gines
 */

#ifndef PCUDADEF_H_
#define PCUDADEF_H_

#include "parserbin.h"
#include "pcuda_types.h"

/*******************/
/* Data structures */
/*******************/

struct _pcuda_configuration {
	Rodruleset rodruleset;
	Rodmultiset rodmultiset;
	Membraneset membraneset;
	ushort * multisets;
	uint *  skin_multiset;

	int skin_label;
	int skin_charge;
	ushort * environment;

	uint rodrule_len;
	uint rodmultiset_len;
	uint mmbrset_len;
	uint multisets_len;
	uint skinmultiset_len;
	
	uint current_step;
	uint max_num_membranes;

	uint block_size;

	ushort * selec_evol;
	uint * selec_siodd;
	
	int verbose;
};

typedef struct _pcuda_configuration* Pcuda_configuration;

struct _problem_size {
	uint max_membranes;
	uint num_objects;
	uint num_labels;
	uint step_limit;
	int num_levels;
};

typedef struct _problem_size* Problem_size;



/***********************/
/* Auxiliary functions */
/***********************/

int select_skin (Configuration *cfg, Pcuda_configuration pcfg);

int execute_skin (Configuration *cfg, Pcuda_configuration pcfg);

void print_multisets(Configuration * cfg, Pcuda_configuration pcfg, Problem_size ps);

void print_skin(Configuration * cfg, Pcuda_configuration pcfg, Problem_size ps);

void print_environment(Configuration * cfg, Pcuda_configuration pcfg, Problem_size ps);


#endif /* PCUDADEF_H_ */

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
 * sequential.h
 *
 *	Describe the sequential implementation of the algorithm for P systems
 *  Created on: 28-dic-2008
 *      Author: miguel
 */

#ifndef SEQUENTIAL_H_
#define SEQUENTIAL_H_

#include "parserbin.h"
#include "pcuda.h"

/**
 * Sequential_solve: performs the sequential and iterative algorithm for P systems
 * with two phases: selection and execution.
 */
void sequential_solve(Configuration *c0, int verbose_level=0, int step_limit=256, int max_membranes=1024, int num_objects=2500, int block_size=256, int threshold=1, int mode=0);

/**
 * Select_rules: selects the rules that can be executed in each membrane.
 * Returns true if there are some rule selected.
 */
bool select_rules(Configuration *cfg);

/**
 * Execute_rules: executes the rules selected by the function select_rules.
 * Returns the number of rules executed.
 */
int execute_rules(Configuration *cfg);

#endif /* SEQUENTIAL_H_ */

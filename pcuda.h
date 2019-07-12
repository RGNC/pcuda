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
 * pcuda.h
 * 
 * Describe the parallel implementation of the algorithm for P systems using CUDA
 * Created on: 4-march-2009
 */

#ifndef PCUDA_H_
#define PCUDA_H_

#include "parserbin.h"
#include "pcudadef.h"
#include "seq2pcuda.h"
#include "timestat.h"
#include "fast_sequential.h"

/*
 * pcuda_solve: performs the parallel algorithm for P systems using CUDA
 * with two phases: selection and execution.
 */
void pcuda_solve(Configuration *cfg, uint step_limit=256, uint current_step=0, uint max_membranes=1024, uint num_objects=2500, uint block_size=256, int mode=0, int verbose=0);

#endif /* PCUDA_H_ */


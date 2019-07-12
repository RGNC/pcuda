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


#if !defined(__EVOLUTION_EXECUTION_KERNEL_CU__)
#define __EVOLUTION_EXECUTION_KERNEL_CU__

#include "pcuda_types.h"
#include "p_system.h"

__global__ void evolution_execution_kernel(const Rodruleset rules, const Rodmultiset objEvo,
                                           const Membraneset membranes, const uint numMemb, 
                                           const uint numLabs, const ushort numObjects,
                                           ushort * multisets, ushort * rev) {
    __shared__ short label;
    __shared__ short charge;

    const short objsPerThread = numObjects/blockDim.x;
    const uint lid = threadIdx.x * objsPerThread;
    const uint bid = blockIdx.y * gridDim.x + blockIdx.x;
    const uint tid = bid * numObjects + lid;
    ushort multiplicity;
    uint begin;
    ushort len;
    uint max;

    if (bid >= numMemb)
	return;
    
    label = membranes[bid].label;
    charge = membranes[bid].charge;

    if (label == EMPTY_MEMBRANE)
        return;

    #pragma unroll
    for (uint obj=0; obj < objsPerThread; obj++) {
        multiplicity = rev[tid + obj];

        if (multiplicity > 0) {
            begin = rules[(lid + obj) * numLabs * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rev & 0x000FFFFF;
            len = (rules[(lid + obj) * numLabs * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rev&0xFFF00000)>>20;
            max = begin + len;

            #pragma unroll
            for (uint i=begin; i < max; i++)
                multisets[bid * numObjects + objEvo[i].obj] += objEvo[i].mult * multiplicity;
                //atomicAdd((uint *)&multisets[bid * numObjects + objEvo[i].obj], (uint)(objEvo[i].mult * multiplicity));
        }
    }
}

#endif	/* __EVOLUTION_EXECUTION_KERNEL_CU__ */

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


#if !defined(__DIVISION_EXECUTION_KERNEL_CU__)
#define __DIVISION_EXECUTION_KERNEL_CU__

#include "pcuda_types.h"
#include "p_system.h"

__global__ void division_execution_kernel(const Rodruleset rules, Membraneset membranes,
                                          uint * numMemb, uint maxNumMemb, const ushort numLabs,
                                          const ushort numObjects, ushort * multisets,
                                          uint * rsiodd) {
    __shared__ short label;
    __shared__ short charge;
    __shared__ uint newMemb;

    const short objsPerThread = numObjects/blockDim.x;
    const uint lid = threadIdx.x * objsPerThread;
    const uint bid = blockIdx.y * gridDim.x + blockIdx.x;
    const uint tid = bid * numObjects + lid;
    
    if (bid >= maxNumMemb)
        return;
    
    label = membranes[bid].label;
    charge = membranes[bid].charge;
    
    if (label == EMPTY_MEMBRANE)
        return;

    uint rule = rsiodd[bid] & 0x1FFFFFFF;
    ushort type = rsiodd[bid]>>29;

    if (type == RULE_DIVISION) {
        if (lid == 0) {
            newMemb = atomicAdd(numMemb, 1);
            membranes[bid].charge = rules[rule * numLabs * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>6 & 0x0003;
            membranes[newMemb].charge = rules[rule * numLabs * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>8 & 0x0003;
            membranes[newMemb].label = label;
        }

        __syncthreads();

        #pragma unroll
        for (uint obj=0; obj < objsPerThread; obj++)
            multisets[newMemb * numObjects + lid + obj] = multisets[tid + obj];

        __syncthreads();

        if (lid == 0) {
            ushort obj = rules[rule * numLabs * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rsodd[0];
            //printf("rule=%u, label=%u, charge=%u, obj1=%u\n",rule,label,charge,obj);
            if (obj != 0)
                multisets[bid * numObjects + obj]++;
            obj = rules[rule * numLabs * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rsodd[1];
            //printf("rule=%u, label=%u, charge=%u, obj2=%u\n",rule,label,charge,obj);
            if (obj != 0)
                multisets[newMemb * numObjects + obj]++;
        }
    }
}

#endif	/* __DIVISION_EXECUTION_KERNEL_CU__ */

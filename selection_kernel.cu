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


#if !defined(__SELECTION_KERNEL_CU__)
#define __SELECTION_KERNEL_CU__

#define NO_LOCK_SI			0
#define LOCK_SI				1

#include "pcuda_types.h"
#include "p_system.h"

__global__ void selection_kernel(const Rodruleset rules, const Rodmultiset objEvo, const Membraneset membranes,
                                 const uint numMemb, const ushort numLabs, const ushort numObjects,
                                 ushort * multisets, uint * skin, uint * rsiodd, 
                                 ushort * srsiodd) {
    __shared__ uint rulesiodd;
    __shared__ uint locksishared;
    __shared__ short label;
    __shared__ short charge;
    __shared__ ushort srsioddshared;
    const short objsPerThread = numObjects/blockDim.x;

    _rulesoddset rulessodd[MAX_OBJECTS_PER_THREAD];
    bool rulessi[MAX_OBJECTS_PER_THREAD];
    ushort multiplicity[MAX_OBJECTS_PER_THREAD];
    ushort revs[MAX_OBJECTS_PER_THREAD];
    bool select_rules[MAX_OBJECTS_PER_THREAD];

    const uint lid = threadIdx.x * objsPerThread;
    const uint bid = blockIdx.y * gridDim.x + blockIdx.x;
    const uint tid = bid * numObjects + lid;

    uint begin;
    uint len;
    uint max;

    if (bid >= numMemb)
	return;
    
    label = membranes[bid].label;
    charge = membranes[bid].charge;

    if (label == EMPTY_MEMBRANE)
        return;

    #pragma unroll
    for (uint obj=0; obj < objsPerThread; obj++) {
        multiplicity[obj] = multisets[tid + obj];
        revs[obj] = rules[(lid + obj) * numLabs * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo & 0x0001;
        rulessodd[obj].rsodd = (rules[(lid + obj) * numLabs * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>1) & 0x0001;
        if ((multiplicity[obj] > 0) && (rulessodd[obj].rsodd != EMPTY_RULE))
            rulessodd[obj].soddtype = (rules[(lid + obj) * numLabs * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>3) & 0x0007;
        else
            rulessodd[obj].soddtype = EMPTY_RULE;
        if (skin[lid + obj] > 0)
            rulessi[obj] = (rules[(lid + obj) * numLabs * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>2) & 0x0001;
        else
            rulessi[obj] = EMPTY_RULE;
        select_rules[obj] = false;
        //printf("Thread %u: mult=%u, revs=%u, rulessodd=%u, rulessodd.type=%u, rulessi=%u\n",tid,multiplicity[obj],revs[obj],rulessodd[obj].rsodd,rulessodd[obj].soddtype,rulessi[obj]);
    }
    uint skinmultold;
    uint rsioddold;
    rulesiodd = 0;
    locksishared = NO_LOCK_SI;
    uint locksi = LOCK_SI; 
    srsioddshared = 0;

    /* DISSOLUTION */
    #pragma unroll
    for (uint obj=0; ((obj < objsPerThread) && (rulesiodd == EMPTY_RULE)); obj++) {
        __syncthreads();
        if (rulessodd[obj].soddtype == RULE_DISSOLUTION) {
            rsioddold = atomicCAS(&rulesiodd, EMPTY_RULE, lid + obj);
            if (rsioddold == EMPTY_RULE) {
                multiplicity[obj]--;
                srsioddshared |= 0x0001;
                rulesiodd |= ((uint)RULE_DISSOLUTION)<<29;
                select_rules[obj] = true;
                //printf("Thread %u selecciona disolución\n\n", tid + obj);
            }
        }
        __syncthreads();
    }

    /* EVOLUTION */
    #pragma unroll
    for (uint obj=0; obj < objsPerThread; obj++)
        if ((revs[obj] != EMPTY_RULE) && (multiplicity[obj] > 0)) {
            srsioddshared |= 0x0002;
            revs[obj] = multiplicity[obj];
            multiplicity[obj] = 0;
            select_rules[obj] = true;
            //printf("Thread %u selecciona evolución\n\n", tid + obj);
        } else
            revs[obj] = 0;

    //__syncthreads();

    /* SEND_OUT */
    #pragma unroll
    for (uint obj=0; ((obj < objsPerThread) && (rulesiodd == EMPTY_RULE)); obj++) {
        __syncthreads();
        if ((multiplicity[obj] > 0) && (rulessodd[obj].soddtype == RULE_SEND_OUT)) {
            rsioddold = atomicCAS(&rulesiodd, EMPTY_RULE, lid + obj);
            if (rsioddold == EMPTY_RULE) {
                multiplicity[obj]--;
                srsioddshared |= 0x0004;
                rulesiodd |= ((uint)RULE_SEND_OUT)<<29;
                //printf("Thread %u selecciona send out\n\n", tid + obj);
                select_rules[obj] = true;
            }
        }
        __syncthreads();
    }

    //__syncthreads();

    /* SEND_IN */
    #pragma unroll
    for (uint obj=0; ((obj < objsPerThread) && (rulesiodd == EMPTY_RULE)); obj++) {
        __syncthreads();
        if (rulessi[obj] != EMPTY_RULE) {
            while ((locksi == LOCK_SI) && (rulesiodd == EMPTY_RULE)) {
                locksi = atomicCAS(&locksishared, NO_LOCK_SI, LOCK_SI);
                if (locksi == NO_LOCK_SI) {
                    skinmultold = atomicDec(&skin[lid + obj], UINT_MAX);
                    //skinmultold = skin[lid + obj]--;
                    if ((skinmultold >= 1) && (skinmultold <= (UINT_MAX - numMemb))) {
                        rulesiodd = lid + obj;
                        srsioddshared |= 0x0008;
                        rulesiodd |= ((uint)RULE_SEND_IN)<<29;
                        select_rules[obj] = true;
                        //printf("Thread %u selecciona send in\n\n", tid + obj);
                    } else {
                        skin[lid + obj] = 0;
                        locksishared = NO_LOCK_SI;
                    }
                }
            }
            rsioddold =  LOCK_SI;
        }
        __syncthreads();
    }

    //__syncthreads();

    /* DIVISION */
    #pragma unroll
    for (uint obj=0; ((obj < objsPerThread) && (rulesiodd == EMPTY_RULE)); obj++) {
        __syncthreads();
        if ((multiplicity[obj] > 0) && (rulessodd[obj].soddtype == RULE_DIVISION)) {
            rsioddold = atomicCAS(&rulesiodd, EMPTY_RULE, lid + obj);
            if (rsioddold == EMPTY_RULE) {
                multiplicity[obj]--;
                srsioddshared |= 0x0010;
                rulesiodd |= ((uint)RULE_DIVISION)<<29;
                select_rules[obj] = true;
                //printf("Thread %u selecciona división\n\n", tid + obj);
            }
        }
        __syncthreads();
    }

    if (lid == 0) {
        atomicOr((uint *)srsiodd, (uint)srsioddshared);
        rsiodd[bid] = rulesiodd;
    }

    for (uint obj=0; obj < objsPerThread; obj++)
        multisets[tid + obj] = multiplicity[obj];

    __syncthreads();

    for (uint obj=0; obj < objsPerThread; obj++) {
        //if (select_rules[obj])
            //multisets[tid + obj] = multiplicity[obj];

        if (revs[obj] != EMPTY_RULE) {
            begin = rules[(lid + obj) * numLabs * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rev & 0x000FFFFF;
            len = (rules[(lid + obj) * numLabs * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rev&0xFFF00000)>>20;
            max = begin + len;

            #pragma unroll
            for (uint i=begin; i < max; i++)
                multisets[bid * numObjects + objEvo[i].obj] += objEvo[i].mult * revs[obj];
        }
    }


    /*if (lid == 0) {
        atomicOr((uint *)srsiodd, (uint)srsioddshared);
        rsiodd[bid] = rulesiodd;
    }*/
}

#endif	/* __SELECTION_KERNEL_CU__ */

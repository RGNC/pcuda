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


// includes, system
#include <stdlib.h>
#include <math.h>
#include <iostream>

// includes, project
#include "pcuda_types.h"
//#include "p_system.h"
#include "pcudadef.h"

// includes, kernels
#include <selection_kernel_only.cu>

using namespace std;

uint dev;
cudaDeviceProp deviceProp;

extern "C"
void eraseDeviceRules(Rodruleset d_rules) {
    cutilSafeCall(cudaFree(d_rules));
}

extern "C"
Rodruleset selection_only(const Rodruleset h_rules, const Membraneset h_membranes,
                     const uint numMemb, const ushort numLabs, const ushort numObjects, const uint blockSize,
                     const ushort* h_multiselec, ushort* h_selec_evol, uint * h_skin,
                     uint * h_siodd, float& time, float& timeTransf, float& timeMalloc, uint firstTime) {
    Rodruleset d_rules;
    Membraneset d_membranes;
    ushort * d_multiselec;
    uint * d_skin;
    uint * d_siodd;
    uint blocksPerRow, rowsPerGrid;
    //cudaDeviceProp deviceProp;
    //uint dev;
    size_t rulesSize, membranesSize, numLabelsSize, multiselecSize, skinSize, sioddSize; 
    size_t deviceGlobalMem;

    if (firstTime){	
    	cudaSetDevice(dev = cutGetMaxGflopsDeviceId());

    	//dev = cutGetMaxGflopsDeviceId();
    	cutilSafeCall(cudaGetDeviceProperties(&deviceProp, dev));
    }

    // calculate sizes
    rulesSize = numLabs * NUMBER_OF_CHARGES * numObjects * sizeof(_rodruleset);
    membranesSize = numMemb * sizeof(_membraneset);
    numLabelsSize = sizeof(ushort);
    multiselecSize = numMemb * numObjects * sizeof(ushort);
    skinSize = numObjects * sizeof(uint);
    sioddSize = numMemb * sizeof(uint);
    cout << "hay " << numMemb << " membranas y " << numObjects << " objetos" << endl;
    
    deviceGlobalMem = rulesSize + membranesSize + numLabelsSize + multiselecSize +
                      skinSize + sioddSize;

    // test conditions
    cutilCondition(deviceProp.minor > 1);
    cutilCondition(numMemb <= deviceProp.maxGridSize[0] * deviceProp.maxGridSize[1]);
    cutilCondition(blockSize <= deviceProp.maxThreadsPerBlock);
    cutilCondition(deviceGlobalMem <= deviceProp.totalGlobalMem);

    // create and start timer
    unsigned int timer = 0;
    cutilCheckError(cutCreateTimer(&timer));
    cutilCheckError(cutStartTimer(timer));

    // allocate device memory
    if (firstTime)
        cutilSafeCall(cudaMalloc((void **) &d_rules, rulesSize));
    cutilSafeCall(cudaMalloc((void **) &d_membranes, membranesSize));
    cutilSafeCall(cudaMalloc((void **) &d_multiselec, multiselecSize));
    cutilSafeCall(cudaMalloc((void **) &d_skin, skinSize));
    cutilSafeCall(cudaMalloc((void **) &d_siodd, sioddSize));

    // stop and destroy timer
    cutilCheckError(cutStopTimer(timer));
    timeMalloc += cutGetTimerValue(timer);
    cutilCheckError(cutDeleteTimer(timer));

    // create and start timer
    timer = 0;
    cutilCheckError(cutCreateTimer(&timer));
    cutilCheckError(cutStartTimer(timer));

    // copy host memory to device
    if (firstTime)
        cutilSafeCall(cudaMemcpy(d_rules, h_rules, rulesSize,
                      cudaMemcpyHostToDevice));
    else
        d_rules = h_rules;
    cutilSafeCall(cudaMemcpy(d_membranes, h_membranes, membranesSize,
                  cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMemcpy(d_multiselec, h_multiselec, multiselecSize,
                  cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMemcpy(d_skin, h_skin, skinSize, cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMemcpy(d_siodd, h_siodd, sioddSize, cudaMemcpyHostToDevice));

    // stop and destroy timer
    cutilCheckError(cutStopTimer(timer));
    timeTransf += cutGetTimerValue(timer);
    cutilCheckError(cutDeleteTimer(timer));

    // setup execution parameters
    if (numMemb <= deviceProp.maxGridSize[0]) {
        // We can use a 1D Grid
        blocksPerRow = numMemb;
        rowsPerGrid  = 1;
    } else {
        // We need to use a 2D Grid
        blocksPerRow = rowsPerGrid = (uint) sqrt(numMemb);

        while ((blocksPerRow * rowsPerGrid) < numMemb)
            blocksPerRow++;
    }

    dim3 grid(blocksPerRow, rowsPerGrid);
    dim3 threads(blockSize);

    // create and start timer
    timer=0;
    cutilCheckError(cutCreateTimer(&timer));
    cutilCheckError(cutStartTimer(timer));
  
    // execute the kernel
    selection_kernel_only <<< grid, threads >>> (d_rules, d_membranes, numMemb, numLabs, 
                                            numObjects, d_multiselec, d_skin, d_siodd);
    
    // check if kernel execution generated and error
    cutilCheckMsg("Kernel execution failed");

    // stop and destroy timer
    cutilCheckError(cutStopTimer(timer));
    time = cutGetTimerValue(timer);
    cutilCheckError(cutDeleteTimer(timer));

    // create and start timer
    timer = 0;
    cutilCheckError(cutCreateTimer(&timer));
    cutilCheckError(cutStartTimer(timer));

    // copy results from device to host
    cutilSafeCall(cudaMemcpy(h_selec_evol, d_multiselec, multiselecSize,
                  cudaMemcpyDeviceToHost));
    cutilSafeCall(cudaMemcpy(h_siodd, d_siodd, sioddSize, cudaMemcpyDeviceToHost));

    // cleanup memory
    cutilSafeCall(cudaFree(d_membranes));
    cutilSafeCall(cudaFree(d_multiselec));
    cutilSafeCall(cudaFree(d_skin));
    cutilSafeCall(cudaFree(d_siodd));

    // stop and destroy timer
    cutilCheckError(cutStopTimer(timer));
    timeTransf += cutGetTimerValue(timer);
    printf("timeTransf para %d bytes = %f\n",multiselecSize+sioddSize,cutGetTimerValue(timer));
    cutilCheckError(cutDeleteTimer(timer));


    printf("timeTransf=%f\n",timeTransf);

    return d_rules;
}


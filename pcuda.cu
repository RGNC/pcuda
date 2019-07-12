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
//#include "parserbin.h"
#include "pcudadef.h"

// includes, kernels
#include <selection_kernel.cu>
#include <dissolution_execution_kernel.cu>
#include <evolution_execution_kernel.cu>
#include <sendout_execution_kernel.cu>
#include <sendin_execution_kernel.cu>
#include <division_execution_kernel.cu>


using namespace std;


/***********************/
/* AUXILIARY FUNCTIONS */

/*int select_skin (Configuration *cfg, Pcuda_configuration pcfg) {
    int mult=0;
    Rulelist * rulelist=NULL;
    Rule * rule=NULL;
    bool rule_end=false;

    cfg->new_selection();
    rulelist=cfg->get_rules()->get_rulelist(pcfg->skin_label,pcfg->skin_charge);

    rulelist->start_iteration();
    while (!rulelist->end() && !rule_end) {
        rule=rulelist->next_rule();
        if (rule==NULL) break;

        if ((rule->type == RULE_EVOLUTION) && (pcfg->skin_multiset[rule->a]>0 )) {
            mult=pcfg->skin_multiset[rule->a];
            pcfg->skin_multiset[rule->a]=0;
            if (mult>0)
                cfg->get_selection()->add_selected_rule(rule,cfg->get_skin(),mult);
        }
        else if ((rule->type == RULE_SEND_OUT) && (pcfg->skin_multiset[rule->a]>0 )) {
            mult=pcfg->skin_multiset[rule->a];
            if (mult>0) {
                pcfg->skin_multiset[rule->a]--;
                cfg->get_selection()->add_selected_rule(rule,cfg->get_skin(),1);
                rule_end=true;
            }
        }
    }

    return cfg->get_selection()->get_length();
}

int execute_skin (Configuration *cfg, Pcuda_configuration pcfg) {
    Rule* rule=NULL;
    int count=0,num_exec=0;

    if (cfg==NULL||pcfg==NULL) return -1;

    cfg->get_selection()->start_iteration();
    while (!cfg->get_selection()->end() ) {
        rule=cfg->get_selection()->get_rule();
        count=cfg->get_selection()->get_times();

        cfg->get_selection()->next_selection();

        switch (rule->type) {
            case RULE_EVOLUTION: 
                rule->new_multiset->add_to_array(pcfg->skin_multiset,pcfg->skinmultiset_len,count);
                num_exec++;
                break;
            case RULE_SEND_OUT:
                pcfg->environment[rule->b]++;
                pcfg->skin_charge=rule->new_charge;
                num_exec++;
                break;
        }
    }
    cfg->delete_selection();

    return num_exec;
}

void print_multisets(Configuration * cfg, Pcuda_configuration pcfg, Problem_size ps) {
    uint num_objects=ps->num_objects;
    string out;
    char convert[50];

    cout << "***********************************" << endl << "MULTISETS: " << endl;

    for (int i=0;i<pcfg->max_num_membranes;i++) {
	if (pcfg->membraneset[i].label==EMPTY_MEMBRANE) continue;
        cout << "MEMBRANE ID: " << i << ", Label: " << pcfg->membraneset[i].label << ", Charge: " << pcfg->membraneset[i].charge << endl;
        cout << "Multiset: ";
        out = "";
        for (int j=0; j < num_objects; j++) {     
            if (pcfg->multisets[i*num_objects+j]>0) {
                out+=cfg->get_alphabet()[j];
        
                if (pcfg->multisets[i*num_objects+j]>1) {
                    out+="*";
                    sprintf(convert,"%d",pcfg->multisets[i*num_objects+j]);
                    out+=convert;
                }
                out+=", ";
            }
        }
        cout << out << endl;
    }
    cout << endl;
}

void print_skin(Configuration * cfg, Pcuda_configuration pcfg, Problem_size ps) {
    uint num_objects=ps->num_objects;
    string out;
    char convert[50];

    cout << "***********************************" << endl << "SKIN:" << endl;

    out="";
    for (int j=0; j < num_objects; j++) {
        if (pcfg->skin_multiset[j]>0) {
            out+=cfg->get_alphabet()[j];
            
            if (pcfg->skin_multiset[j]>1) {
                out+="*";
                sprintf(convert,"%d",pcfg->skin_multiset[j]);
                out+=convert;
            }
            out+=", ";
        }
    }
    cout << out << endl;
}

void print_environment(Configuration * cfg, Pcuda_configuration pcfg, Problem_size ps) {
    uint num_objects=ps->num_objects;
    string out;
    char convert[50];

    cout << "***********************************" << endl << "ENVIRONMENT" << endl;

    out="";
    for (uint j=0; j < num_objects; j++) {
        if (pcfg->environment[j]>0) {
            out+=cfg->get_alphabet()[j];

            if (pcfg->environment[j]>1) {
                out+="*";
                sprintf(convert,"%d",pcfg->environment[j]);
                out+=convert;
            }
                out+=", ";
        }
    }
    cout << out ;
}
*/

extern "C"
void pcuda(Configuration * cfg, Pcuda_configuration pcfg, Problem_size ps) {
          /*const Rodruleset h_rules, const Rodobjevo * h_rulesEvo,
           const uint rulesEvoLen, Membraneset h_membranes, uint numMemb, 
           const uint maxNumMemb, const uint numLabs, const uint numObjects,
           const uint blockSize, ushort * h_multisets, ushort * h_skin,
           const uint stepLimit, Configuration * cfg, Pcuda_configuration pcfg) {*/
    Rodruleset d_rules;
    Rodmultiset d_rulesEvo;
    Membraneset d_membranes;
    uint * d_numMembranes;
    ushort * d_multisets;
    uint * d_skin;
    uint * d_rsiodd;
    ushort * d_srsiodd, srsiodd = 0;
    uint numMembNew;
    uint blocksPerRow, rowsPerGrid;
    cudaDeviceProp deviceProp;
    uint dev;
    uint curStep;
    size_t rulesSize, rulesEvoSize, membranesSize, maxMembranesSize;
    size_t numMembranesSize, multisetsSize, maxMultisetsSize, skinSize;
    size_t maxRevSize, maxRsioddSize, srsioddSize, deviceGlobalMem;
    dim3 grid, threads;
    uint timer = 0;
    bool end = false;
    int skinExec;
    float timeTransf, timeMalloc, timeSelect = 0.0, timeExecut = 0.0;
    float timeSelect_stp, timeExecut_stp;
    int verbose=pcfg->verbose;
    
    /* Initialize GPU */
    char * def_dev = getenv("DEFAULT_DEVICE");
    if (def_dev!=NULL)
	cudaSetDevice(dev= atoi(def_dev));
    else
	cudaSetDevice(dev = cutGetMaxGflopsDeviceId());
    
    cutilSafeCall(cudaGetDeviceProperties(&deviceProp, dev));

    // calculate sizes
    rulesSize = ps->num_labels * NUMBER_OF_CHARGES * ps->num_objects * sizeof(_rodruleset);
    rulesEvoSize = pcfg->rodmultiset_len * sizeof(_rodmultiset);
    membranesSize = pcfg->max_num_membranes * sizeof(short2);
    maxMembranesSize = ps->max_membranes * sizeof(short2);
    //cout << "hay " << ps->max_membranes << " membranas y " << ps->num_objects << " objetos" << endl;
    numMembranesSize = sizeof(uint);
    multisetsSize = pcfg->max_num_membranes * ps->num_objects * sizeof(ushort);
    maxMultisetsSize = ps->max_membranes * ps->num_objects * sizeof(ushort);
    skinSize = ps->num_objects * sizeof(uint);
    maxRevSize = ps->max_membranes * ps->num_objects * sizeof(ushort);
    maxRsioddSize = ps->max_membranes * sizeof(uint);
    srsioddSize = sizeof(ushort);
    
    deviceGlobalMem = rulesSize + rulesEvoSize + maxMembranesSize + numMembranesSize +
                      maxMultisetsSize + skinSize + maxRevSize + maxRsioddSize +
                      srsioddSize;

    // test conditions
    cutilCondition(deviceProp.minor > 1);
    cutilCondition(ps->max_membranes <= deviceProp.maxGridSize[0] * deviceProp.maxGridSize[1]);
    cutilCondition(pcfg->block_size <= deviceProp.maxThreadsPerBlock);
    cutilCondition(deviceGlobalMem <= deviceProp.totalGlobalMem);

    // create and start timer
    timer = 0;
    cutilCheckError(cutCreateTimer(&timer));
    cutilCheckError(cutStartTimer(timer));

    // allocate device memory
    cutilSafeCall(cudaMalloc((void **) &d_rules, rulesSize));
    cutilSafeCall(cudaMalloc((void **) &d_rulesEvo, rulesEvoSize));
    cutilSafeCall(cudaMalloc((void **) &d_numMembranes, numMembranesSize));
    cutilSafeCall(cudaMalloc((void **) &d_membranes, maxMembranesSize));
    cutilSafeCall(cudaMalloc((void **) &d_multisets, maxMultisetsSize));
    cutilSafeCall(cudaMalloc((void **) &d_skin, skinSize));
    cutilSafeCall(cudaMalloc((void **) &d_rsiodd, maxRsioddSize));
    cutilSafeCall(cudaMalloc((void **) &d_srsiodd, srsioddSize));

    // stop and destroy timer
    cutilCheckError(cutStopTimer(timer));
    timeMalloc = cutGetTimerValue(timer);
    cutilCheckError(cutDeleteTimer(timer));

    // create and start timer
    timer = 0;
    cutilCheckError(cutCreateTimer(&timer));
    cutilCheckError(cutStartTimer(timer));

    // copy host memory to device
    cutilSafeCall(cudaMemcpy(d_rules, pcfg->rodruleset, rulesSize,
                  cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMemcpy(d_rulesEvo, pcfg->rodmultiset, rulesEvoSize,
                  cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMemcpy(d_membranes, pcfg->membraneset, membranesSize,
                  cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMemcpy(d_numMembranes, &pcfg->max_num_membranes, numMembranesSize,
                  cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMemcpy(d_multisets, pcfg->multisets, multisetsSize,
                  cudaMemcpyHostToDevice));
    cutilSafeCall(cudaMemcpy(d_skin, pcfg->skin_multiset, skinSize, cudaMemcpyHostToDevice));
    
    // stop and destroy timer
    cutilCheckError(cutStopTimer(timer));
    timeTransf = cutGetTimerValue(timer);
    cutilCheckError(cutDeleteTimer(timer));

    for (curStep=pcfg->current_step; ((!end) && (curStep<ps->step_limit)); curStep++) {
	//cout << "hay " << pcfg->max_num_membranes << " membranas" << endl;
        numMembNew = pcfg->max_num_membranes;

        select_skin(cfg,pcfg);

        // setup execution parameters
        if (pcfg->max_num_membranes <= deviceProp.maxGridSize[0]) {
            // We can use a 1D Grid
            blocksPerRow = pcfg->max_num_membranes;
            rowsPerGrid  = 1;
        } else {
            // We need to use a 2D Grid
            blocksPerRow = rowsPerGrid = (uint) sqrt(pcfg->max_num_membranes);

            while ((blocksPerRow * rowsPerGrid) < pcfg->max_num_membranes)
                blocksPerRow++;
        }

        grid = dim3(blocksPerRow, rowsPerGrid);
        threads = dim3(pcfg->block_size);

	// create and start timer
        timer = 0;
        cutilCheckError(cutCreateTimer(&timer));
        cutilCheckError(cutStartTimer(timer));
	
        srsiodd = 0;
        // copy host memory to device
        cutilSafeCall(cudaMemcpy(d_srsiodd, &srsiodd, srsioddSize,
                      cudaMemcpyHostToDevice));
	cutilSafeCall(cudaMemcpy(d_skin, pcfg->skin_multiset, skinSize,
		      cudaMemcpyHostToDevice));
	
	    // stop and destroy timer
        cutilCheckError(cutStopTimer(timer));
        timeTransf += cutGetTimerValue(timer);
        cutilCheckError(cutDeleteTimer(timer));

        // create and start timer
        timer = 0;
        cutilCheckError(cutCreateTimer(&timer));
        cutilCheckError(cutStartTimer(timer));
  
        // execute the selection kernel
        selection_kernel <<< grid, threads >>> (d_rules, d_rulesEvo, d_membranes, pcfg->max_num_membranes, ps->num_labels, 
                                                ps->num_objects, d_multisets, d_skin,
                                                d_rsiodd, d_srsiodd);
	//cutilDeviceSynchronize();	
	CUT_DEVICE_SYNCHRONIZE() // Compatible for CUDA < 4.0
		
        // check if kernel execution generated and error
        cutilCheckMsg("Kernel execution failed");

        // stop and destroy timer
        cutilCheckError(cutStopTimer(timer));
        timeSelect_stp = cutGetTimerValue(timer);
        timeSelect += timeSelect_stp;
        cutilCheckError(cutDeleteTimer(timer));

        // create and start timer
        timer = 0;
        cutilCheckError(cutCreateTimer(&timer));
        cutilCheckError(cutStartTimer(timer));

        // copy results from device to host
        cutilSafeCall(cudaMemcpy(&srsiodd, d_srsiodd, srsioddSize,
                      cudaMemcpyDeviceToHost));

	if (srsiodd & 0x0008)
        	cutilSafeCall(cudaMemcpy(pcfg->skin_multiset, d_skin, skinSize,
                	      cudaMemcpyDeviceToHost));

        cutilCheckError(cutStopTimer(timer));
        timeTransf += cutGetTimerValue(timer);
        cutilCheckError(cutDeleteTimer(timer));

        skinExec=execute_skin(cfg,pcfg);

        if (skinExec)
            // copy host memory to device
            cutilSafeCall(cudaMemcpy(d_skin, pcfg->skin_multiset, skinSize, cudaMemcpyHostToDevice));

        timeExecut_stp = 0;
        // create and start timer
        timer = 0;
        cutilCheckError(cutCreateTimer(&timer));
        cutilCheckError(cutStartTimer(timer));

        // execute evolution rules
        if (srsiodd & 0x0002) {
	    if (verbose>0)
            cout << "Execute evolution rules" << endl;
	    // This is done now directly in the selection kernel
	    
            /* evolution_execution_kernel <<< grid, threads >>> (d_rules, d_rulesEvo, d_membranes,
                                                                        pcfg->max_num_membranes, ps->num_labels, ps->num_objects,
                                                                        d_multisets, d_rev);
            // check if kernel execution generated and error
            cutilCheckMsg("Kernel execution failed"); */
        }
        // execute dissolution rules
        if (srsiodd & 0x0001) {
	    if (verbose>0)		
            cout << "Execute dissolution rules" << endl;
	    
            dissolution_execution_kernel <<< grid, threads >>> (d_rules, d_membranes,
                                                                pcfg->max_num_membranes, ps->num_labels, ps->num_objects,
                                                                d_multisets, d_skin, d_rsiodd);
            //cutilDeviceSynchronize();	
	    CUT_DEVICE_SYNCHRONIZE() // Compatible for CUDA < 4.0
		
	    // check if kernel execution generated and error
            cutilCheckMsg("Kernel execution failed");

            // stop and destroy timer
            cutilCheckError(cutStopTimer(timer));
            timeExecut_stp += cutGetTimerValue(timer);
            cutilCheckError(cutDeleteTimer(timer));

            // create and start timer
            timer = 0;
            cutilCheckError(cutCreateTimer(&timer));
            cutilCheckError(cutStartTimer(timer));

	    cutilSafeCall(cudaMemcpy(pcfg->skin_multiset, d_skin, skinSize,
                          cudaMemcpyDeviceToHost));

            // stop and destroy timer
            cutilCheckError(cutStopTimer(timer));
            timeTransf += cutGetTimerValue(timer);
            cutilCheckError(cutDeleteTimer(timer));

            // create and start timer
            timer = 0;
            cutilCheckError(cutCreateTimer(&timer));
            cutilCheckError(cutStartTimer(timer));
        }
        // execute division rules
        if (srsiodd & 0x0010) {
	    if (verbose>0)
            cout << "Execute division rules" << endl; //",gx=" << grid.x <<",gy="<< grid.y << ",t=" << threads.x << endl;
	    
            division_execution_kernel <<< grid, threads >>> (d_rules, d_membranes,
                                                             d_numMembranes, pcfg->max_num_membranes, 
							     ps->num_labels, ps->num_objects,
                                                             d_multisets, d_rsiodd);
            //cutilDeviceSynchronize();	
	    CUT_DEVICE_SYNCHRONIZE() // Compatible for CUDA < 4.0
		    
	    // check if kernel execution generated and error
            cutilCheckMsg("Kernel execution failed");

            // stop and destroy timer
            cutilCheckError(cutStopTimer(timer));
            timeExecut_stp += cutGetTimerValue(timer);
            cutilCheckError(cutDeleteTimer(timer));

            // create and start timer
            timer = 0;
            cutilCheckError(cutCreateTimer(&timer));
            cutilCheckError(cutStartTimer(timer));

            // copy results from device to host
            cutilSafeCall(cudaMemcpy(&numMembNew, d_numMembranes, numMembranesSize,
                          cudaMemcpyDeviceToHost));

            // stop and destroy timer
            cutilCheckError(cutStopTimer(timer));
            timeTransf += cutGetTimerValue(timer);
            cutilCheckError(cutDeleteTimer(timer));

            // create and start timer
            timer = 0;
            cutilCheckError(cutCreateTimer(&timer));
            cutilCheckError(cutStartTimer(timer));
        }
        if ((srsiodd & 0x0004) || (srsiodd & 0x0008)) {
            grid = dim3(ceil((double)pcfg->max_num_membranes/(double)BLOCK_SIZE));
            threads = dim3(BLOCK_SIZE);

            // execute send_out rules
            if (srsiodd & 0x0004) {
		if (verbose>0)
                cout << "Execute send_out rules" << endl;
		
                sendout_execution_kernel <<< grid, threads >>> (d_rules, d_membranes,
                                                                pcfg->max_num_membranes, ps->num_labels,
                                                                d_skin, d_rsiodd);
		
		//cutilDeviceSynchronize();	
	        CUT_DEVICE_SYNCHRONIZE() // Compatible for CUDA < 4.0
		    
                // check if kernel execution generated and error
                cutilCheckMsg("Kernel execution failed");
            }
            // execute send_in rules
            if (srsiodd & 0x0008) {
		if (verbose>0)
                cout << "Execute send_in rules" << endl;
		
                sendin_execution_kernel <<< grid, threads >>> (d_rules, d_membranes,
                                                               pcfg->max_num_membranes, ps->num_labels, ps->num_objects,
                                                               d_multisets, d_rsiodd);
                
		//cutilDeviceSynchronize();	
	        CUT_DEVICE_SYNCHRONIZE() // Compatible for CUDA < 4.0
		
		// check if kernel execution generated and error
                cutilCheckMsg("Kernel execution failed");
            }
            // stop and destroy timer
            cutilCheckError(cutStopTimer(timer));
            timeExecut_stp += cutGetTimerValue(timer);
            cutilCheckError(cutDeleteTimer(timer));

            // create and start timer
            timer = 0;
            cutilCheckError(cutCreateTimer(&timer));
            cutilCheckError(cutStartTimer(timer));

            // copy results from device to host
            cutilSafeCall(cudaMemcpy(pcfg->skin_multiset, d_skin, skinSize,
                          cudaMemcpyDeviceToHost));

            // stop and destroy timer
            cutilCheckError(cutStopTimer(timer));
            timeTransf += cutGetTimerValue(timer);
            cutilCheckError(cutDeleteTimer(timer));

            // create and start timer
            timer = 0;
            cutilCheckError(cutCreateTimer(&timer));
            cutilCheckError(cutStartTimer(timer));
        }

	
	end = (skinExec==0 && srsiodd==0);
	uint prev_num_membranes = pcfg->max_num_membranes;
        pcfg->max_num_membranes = numMembNew;
	
        // stop and destroy timer
        cutilCheckError(cutStopTimer(timer));
        timeExecut_stp += cutGetTimerValue(timer);
        timeExecut += timeExecut_stp;
        cutilCheckError(cutDeleteTimer(timer));

	if (verbose>0) {
	cout << "**************************************" << endl;
	cout << "CONFIGURATION " << curStep << ": from " << prev_num_membranes << " to " << pcfg->max_num_membranes << " membranes:" << endl;
        cout << "TIME: sel= " << timeSelect_stp << " , exec= " << timeExecut_stp  << ", acum memory transf= " << timeTransf << ", step total= " << timeSelect_stp + timeExecut_stp << endl;
	}
        
	/*cout << endl << "Step " << curStep << ":" << endl;
        cout << "\tSelection time: " << timeSelect_stp << " ms" << endl;
        cout << "\tExecution time: " << timeExecut_stp << " ms" << endl;
        cout << "\tTotal time: " << timeSelect_stp + timeExecut_stp << " ms" << endl << endl;*/

        // copy results from device to host
        /*cutilSafeCall(cudaMemcpy(pcfg->membraneset, d_membranes, maxMembranesSize,
                     cudaMemcpyDeviceToHost));
        cutilSafeCall(cudaMemcpy(pcfg->multisets, d_multisets, maxMultisetsSize,
                     cudaMemcpyDeviceToHost));
        print_environment(cfg, pcfg, ps);
        print_skin(cfg, pcfg, ps);
        print_multisets(cfg, pcfg, ps);*/

        //printf("termina selection_only.cu con timeTransf=%f\n",timeTransf);
    }

    if (curStep>=ps->step_limit && verbose>0) 
	    cout << "Reached step limit: " << ps->step_limit-1 << endl;
    else if (verbose>0) cout << "Finished on step: " << curStep << endl;

    if (verbose>0) {	    
    print_environment(cfg, pcfg, ps);
    print_skin(cfg, pcfg, ps); }
    
    if (verbose>1) {
    // create and start timer
    timer = 0;
    cutilCheckError(cutCreateTimer(&timer));
    cutilCheckError(cutStartTimer(timer)); 

    // copy results from device to host
    cutilSafeCall(cudaMemcpy(pcfg->membraneset, d_membranes, maxMembranesSize/2, 
                  cudaMemcpyDeviceToHost));
    cutilSafeCall(cudaMemcpy(pcfg->multisets, d_multisets, maxMultisetsSize/2,
                  cudaMemcpyDeviceToHost));

    // stop and destroy timer
    cutilCheckError(cutStopTimer(timer));
    timeTransf += cutGetTimerValue(timer);
    cutilCheckError(cutDeleteTimer(timer)); 
    
    print_multisets(cfg, pcfg, ps);
    printf("timeTransf for %d bytes = %f\n",maxMembranesSize/2+maxMultisetsSize/2,cutGetTimerValue(timer));}
 

    cout << "**************************************" << endl;
    cout << "PCUDA PARALLEL TIME: " << endl;    
    cout << "\tMalloc time: " << timeMalloc << " ms" << endl;
    cout << "\tTransference time: " << timeTransf << " ms" << endl;
    cout << "\tSelection time: " << timeSelect <<  " ms" << endl;
    cout << "\tExecution time: " << timeExecut << " ms " << endl;
    cout << "\tTotal time = " << timeMalloc+timeTransf+timeSelect+timeExecut <<  " ms" << endl;

    // cleanup memory
    cutilSafeCall(cudaFree(d_rules));
    cutilSafeCall(cudaFree(d_rulesEvo));
    cutilSafeCall(cudaFree(d_membranes));
    cutilSafeCall(cudaFree(d_numMembranes));
    cutilSafeCall(cudaFree(d_multisets));
    cutilSafeCall(cudaFree(d_skin));
    cutilSafeCall(cudaFree(d_rsiodd));
    cutilSafeCall(cudaFree(d_srsiodd));
}

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
#include <timestat.h>
#include "pcudadef.h"

// includes, kernels
//#include <selection_only.cu>


////////////////////////////////////////////////////////////////////////////////
//// declaration, forward
extern "C"
void eraseDeviceRules(Rodruleset d_rules);

extern "C"
Rodruleset selection_only(const Rodruleset h_rules, const Membraneset h_membranes,
                          const uint numMemb, const ushort numLabs, const ushort numObjects, const uint blockSize,
                          const ushort* h_multiselec, ushort* h_selec_evol, uint * h_skin,
                          uint * h_siodd, float& time, float& timeTransf, float& timeMalloc, uint firstTime);


/***********************/
/* FUNCTIONS FOR PCUDA */

int pcuda_execution(Pcuda_configuration pcfg, Problem_size ps);

/*
 * pcuda_solve: performs the main loop of the algorithm
 * with two phases: selection and execution.
 */

extern "C"
void pcuda_sel(Configuration * cfg, Pcuda_configuration pcfg, Problem_size ps) {
//void pcuda_solve(Configuration *cfg, uint step_limit, uint current_step, uint max_membranes, uint num_objects, uint block_size) {

	int skin_selec, skin_exec, num_exec;
	bool end=false;	
	uint first=1;
	float time=0.0;
	float timeTransf = 0.0;
	float timeMalloc = 0.0;
	float timeSelect = 0.0;
	float timeExecut = 0.0;
	Rodruleset cr = pcfg->rodruleset;

	//time_t total_ini=0, total_exec=0,total_sel=0, aux=0, aux2=0;
	//suseconds_t utotal_ini=0, utotal_exec=0,utotal_sel=0, uaux=0, uaux2=0;
	double time_sel=0.0,time_exec=0.0;//,time_ini=0.0;
	init_time();

	pcfg->selec_evol=new ushort[ps->max_membranes*ps->num_objects];
	pcfg->selec_siodd=new uint[ps->max_membranes];

	/* Main loop */
        while ((!end) && (pcfg->current_step++<ps->step_limit)) {
		//cout << endl << "***********************************************" <<endl << endl<< "CONFIGURATION: " << current_step-1 << endl;
		/* First: select the rules for the skin */
//		start_timer();

		select_skin(cfg,pcfg);
		
		start_timer();
		/* Second: call to the parallel method that will select the rules of the membranes, in parallel */
                //num_selec=pcuda_selection(&pcfg,ps.num_objects,ps.max_membranes,cfg->get_label_set_size());
		cr=selection_only(cr,pcfg->membraneset,pcfg->max_num_membranes,ps->num_labels,ps->num_objects,pcfg->block_size,pcfg->multisets,pcfg->selec_evol,pcfg->skin_multiset,pcfg->selec_siodd, time, timeTransf, timeMalloc, first); first*=0;

		timeSelect += time;

		time_sel+=end_timer();

                /* Third: execution of the rules for the skin */
		skin_exec=execute_skin(cfg,pcfg);
		start_timer();

		/* Forth: execution of the rules for each membrane, in a sequential mode */
                num_exec=pcuda_execution(pcfg,ps);

		time_exec=end_timer();

		timeExecut += time_exec;

		cout << endl << "Step " << pcfg->current_step << ":" << endl;
		cout << "\tSelection time: " << time << " ms" << endl;
		cout << "\tExecution time: " << time_exec << " ms" << endl;
		cout << "\tTotal time: " << time + time_exec << " ms" << endl << endl;

		end = (skin_exec==0 && num_exec==0);

		print_environment(cfg, pcfg, ps);
		print_skin(cfg, pcfg, ps);
        }
	
	if (pcfg->current_step>=ps->step_limit)	cout << "Reached step limit: " << ps->step_limit-1 << endl;
	else cout << "Finished on step: " << pcfg->current_step << endl;

	print_multisets(cfg, pcfg, ps);
	//print_last_configuration(cfg,pcfg,ps);

	cout << endl << "TIME PARALLEL: " << endl;
	cout << "\tMalloc time: " << timeMalloc << " ms" << endl;
	cout << "\tTransference time: " << timeTransf << " ms" << endl;
	cout << "\tSelection time: " << timeSelect <<  " ms" << endl;
	cout << "\tExecution time: " << timeExecut << " ms " << endl;

	//cout << endl << "PARALLEL: KERNEL TIME: " << time << "ms, SELECTION TIME: " << time_sel <<  "ms, EXECUTION TIME: " << time_exec << "ms" << endl;
        //cout << "PARALLEL: INITIALIZATION TIME: " << time_ini << "ms, TRANSFERENCE TIME: " << timeTransf << "ms" << endl;
}


/*
 * Execute_rules: executes the rules selected by the function select_rules.
 * Returns the number of rules executed.
 */
int pcuda_execution(Pcuda_configuration pcfg, Problem_size ps) {
	uint rule=0,times=0;
	int num_exec=0,begin=0,len=0,max=0,label=0,charge=0;
	//uint max_membranes=ps->max_membranes;
	uint max_num_membranes=pcfg->max_num_membranes;
	uint num_objects=ps->num_objects;
	uint num_labels=ps->num_labels;
	ushort type=0,obj=0;

	/* Execution of rules */
	for (uint i=0;i<max_num_membranes;i++) {

		label=pcfg->membraneset[i].label;
		if (label==EMPTY_MEMBRANE) continue;
		charge = pcfg->membraneset[i].charge;

		/* Execution of evolution rules for membrane i */
		for (uint j=0;j<num_objects;j++) {

			/* If selected a number of times */
			if (pcfg->selec_evol[i*num_objects + j]>0) {

				times=pcfg->selec_evol[i*num_objects + j];
				pcfg->selec_evol[i*num_objects + j]=0;
				pcfg->multisets[i*num_objects + j]-=times;
				/* Index the multiset of the evolution rule */
	            		begin = pcfg->rodruleset[j * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rev & 0x000FFFFF;
		                len = (pcfg->rodruleset[j * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rev&0xFFF00000)>>20;
				max = begin + len;

				//if (j==313) cout << "beg " << begin << ", len "<<len<< ": ";
				/* Add the evolution rule to the multiset */
				for (int k=begin; k < max; k++) {
					if (pcfg->rodmultiset[k].obj==0) continue;
				//	if (j==313) cout<<pcfg->rodmultiset[k].obj<<" ";
					pcfg->multisets[i*num_objects + pcfg->rodmultiset[k].obj] += pcfg->rodmultiset[k].mult * times;
				}
				//if (j==313) cout<<endl;

				num_exec++;
			}
		}

		/* Execution of siodd rule for membrane i */
		rule = pcfg->selec_siodd[i] & 0x1FFFFFFF;
		type = pcfg->selec_siodd[i]>>29;

		pcfg->selec_siodd[i]=0;

		switch (type) {
		case RULE_SEND_OUT:
			/* Addnew element in skin */
			obj=pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rsodd[0];
//cout << "Sendout: a=" << rule << ", charge=" << charge << ", b=" << obj << endl;
			if (obj!=0) pcfg->skin_multiset[obj]++;
			/* Subs one object */
			if (pcfg->multisets[i*num_objects + rule]>0) pcfg->multisets[i*num_objects + rule]--;
			/*change charge*/
			pcfg->membraneset[i].charge=(pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>6) & 0x0003;
			num_exec++;
		break;
		case RULE_SEND_IN:
			if (pcfg->skin_multiset[rule]>0) pcfg->skin_multiset[rule]--;
			/*Add new object*/
			obj=pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rsin;
//			cout << "Sendin: a=" << rule << ", charge=" << charge << ", b=" << obj << endl;
                        if (obj!=0) pcfg->multisets[i*num_objects + obj]++;
			/*change charge*/
			pcfg->membraneset[i].charge=(pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>10) & 0x0003;
			num_exec++;
                break;
                case RULE_DISSOLUTION:
			/* subs object */
                        if (pcfg->multisets[i*num_objects + rule]>0) pcfg->multisets[i*num_objects + rule]--;
			/*add multiset to the skin */
			for (uint j=0;j<num_objects;j++)
				pcfg->skin_multiset[j]+=pcfg->multisets[i*num_objects+j];
			/*add the new object */
			obj=pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rsodd[0];
			if (obj!=0) pcfg->skin_multiset[obj]++;
			/*set to empty membrane*/
			pcfg->membraneset[i].charge=pcfg->membraneset[i].label=EMPTY_MEMBRANE;

			//pcfg->num_membranes--;
			num_exec++;
                break;
                case RULE_DIVISION:
			/*next id for a new membrane*/
			int next_id=pcfg->max_num_membranes++;
			if (next_id >= ps->max_membranes) {
				cerr << "Error, the p system needs to create more membranes than defined in the imput by -m: "<<ps->max_membranes<<endl;
				exit(0);
			}
			if (pcfg->multisets[i*num_objects+rule]>0) pcfg->multisets[i*num_objects+rule]--;
			pcfg->membraneset[next_id].label=pcfg->membraneset[i].label;
			
			/* Copy the multiset */
			for (uint j=0;j<num_objects;j++)
				pcfg->multisets[next_id*num_objects+j]=pcfg->multisets[i*num_objects+j];

//cout << "Division: a=" << rule << ", charge=" << charge;
			/* Add object b */
			obj=pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rsodd[0];
//			cout<<", b="<<obj;
			if (obj!=0) pcfg->multisets[i*num_objects+obj]++;
			/* Add object c */
			obj=pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rsodd[1];
//			cout <<", c="<<obj<<endl;
			if (obj!=0) pcfg->multisets[next_id*num_objects+obj]++;
			/* Update charges */
			pcfg->membraneset[i].charge=(pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>6) & 0x0003;
			pcfg->membraneset[next_id].charge=(pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>8) & 0x0003;
			//pcfg->num_membranes++;
			num_exec++;
//			cout <<"Despues division queda, obj="<<rule <<", " << pcfg->multisets[i*num_objects + rule]<<endl;
                break;
		};
	}

	return num_exec;
}


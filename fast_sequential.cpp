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
#include "parserbin.h"
#include "pcudadef.h"

#include "pcuda.h"
//#include "timestat.h"

#include "fast_sequential.h"

ushort * selected_times;
ushort * selected_evol;
uint selected_siodd;
uint * skin_aux;
double time_sel, time_ex;

/*
 * Execute_rules: executes the rules selected by the function select_rules.
 * Returns the number of rules executed.
 */
int execute_step(Configuration * cfg, Pcuda_configuration pcfg, Problem_size ps) {
	uint rule=0,times=0;
	int num_exec=0,begin=0,len=0,max=0,label=0,charge=0,num_selec=0;
	uint max_num_membranes=pcfg->max_num_membranes;
	uint num_objects=ps->num_objects;
	uint num_labels=ps->num_labels;
	ushort type=0,obj=0;

	ushort rsodd=0,soddtype=0,rsi=0,rev=0;
	uint mult=0,skin_mult=0;

	start_timer();

	for (uint i=0;i<num_objects;i++)
		skin_aux[i]=0;

	time_sel+=end_timer();

	for (uint i=0;i<max_num_membranes;i++) {

		start_timer();

		label=pcfg->membraneset[i].label;
		if (label==EMPTY_MEMBRANE) continue;
		charge = pcfg->membraneset[i].charge;

		selected_siodd=0;
		num_selec=0;

                /*********************************/
		/* RULE SELECTION FOR MEMBRANE i */
		/*********************************/
		for (uint j=0;j<num_objects;j++) {
		        mult=pcfg->multisets[i*num_objects + j];
		        skin_mult=pcfg->skin_multiset[j];

			if (mult==0 && skin_mult==0) continue;


			rev=pcfg->rodruleset[ j * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo & 0x0001;
			rsodd=(pcfg->rodruleset[ j * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>1) & 0x0001;
			soddtype=(pcfg->rodruleset[ j * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>3) & 0x0007;
			rsi=(pcfg->rodruleset[ j * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>2) & 0x0001;


			/* Selecting dissolution rule */
			if (soddtype == RULE_DISSOLUTION && mult>0 && selected_siodd==0) {
				mult--;
				selected_siodd=j;
				selected_siodd |= ((uint)RULE_DISSOLUTION)<<29;
			}
			if ((rev != 0) && (mult > 0)) {
				selected_evol[num_selec]=j;
				selected_times[num_selec++]=mult;
				mult = 0;
			}
                        if (rsodd!=0 && soddtype == RULE_SEND_OUT && mult>0 && selected_siodd==0) {
                                mult--;
                                selected_siodd=j;
                                selected_siodd |= ((uint)RULE_SEND_OUT)<<29;
                        }
                        if (rsi!=0 && skin_mult>0 && selected_siodd==0) {
                                skin_mult--;
                                selected_siodd=j;
				selected_siodd |= ((uint)RULE_SEND_IN)<<29;
                        }
                        if (rsodd!=0 && soddtype == RULE_DIVISION && mult>0 && selected_siodd==0) {
                                mult--;
                                selected_siodd=j;
				selected_siodd |= ((uint)RULE_DIVISION)<<29;
                        }

			pcfg->multisets[i*num_objects + j]=mult;
			pcfg->skin_multiset[j]=skin_mult;
		}

		time_sel+=end_timer();


                /*********************************/
		/* RULE EXECUTION FOR MEMBRANE i */
		/*********************************/

		start_timer();
		/* Execution of evolution rules for membrane i */
		for (uint j=0;j<num_selec;j++) {

			/* If selected a number of times */
			if (selected_evol[j]>0) {
				times=selected_times[j];
				obj=selected_evol[j];

				/* Index the multiset of the evolution rule */
	            		begin = pcfg->rodruleset[obj * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rev & 0x000FFFFF;
		                len = (pcfg->rodruleset[obj * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rev&0xFFF00000)>>20;
				max = begin + len;

				/* Add the evolution rule to the multiset */
				for (int k=begin; k < max; k++) {
					if (pcfg->rodmultiset[k].obj==0) continue;
					pcfg->multisets[i*num_objects + pcfg->rodmultiset[k].obj] += pcfg->rodmultiset[k].mult * times;
				}

				num_exec++;
			}
		}

		/* Execution of siodd rule for membrane i */
		rule = selected_siodd & 0x1FFFFFFF;
		type = selected_siodd>>29;

		switch (type) {
		case RULE_SEND_OUT:
			/* Addnew element in skin */
			obj=pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rsodd[0];
			if (obj!=0) skin_aux[obj]++;
			/*change charge*/
			pcfg->membraneset[i].charge=(pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>6) & 0x0003;
			num_exec++;
		break;
		case RULE_SEND_IN:
			/*Add new object*/
			obj=pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rsin;
                        if (obj!=0) pcfg->multisets[i*num_objects + obj]++;
			/*change charge*/
			pcfg->membraneset[i].charge=(pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>10) & 0x0003;
			num_exec++;
                break;
                case RULE_DISSOLUTION:
			/* subs object */
			/*add multiset to the skin */
			for (uint j=0;j<num_objects;j++)
				skin_aux[j]+=pcfg->multisets[i*num_objects+j];
			/*add the new object */
			obj=pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rsodd[0];
			if (obj!=0) pcfg->skin_multiset[obj]++;
			/*set to empty membrane*/
			pcfg->membraneset[i].charge=pcfg->membraneset[i].label=EMPTY_MEMBRANE;

			num_exec++;
                break;
                case RULE_DIVISION:
			/*next id for a new membrane*/
			int next_id=pcfg->max_num_membranes++;
			if (next_id >= ps->max_membranes) {
				cerr << "Error, the p system needs to create more membranes than defined in the imput by -m: "<<ps->max_membranes<<endl;
				exit(0);
			}
			pcfg->membraneset[next_id].label=pcfg->membraneset[i].label;
			
			/* Copy the multiset */
			for (uint j=0;j<num_objects;j++)
				pcfg->multisets[next_id*num_objects+j]=pcfg->multisets[i*num_objects+j];

			/* Add object b */
			obj=pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rsodd[0];
			if (obj!=0) pcfg->multisets[i*num_objects+obj]++;
			/* Add object c */
			obj=pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rsodd[1];
			if (obj!=0) pcfg->multisets[next_id*num_objects+obj]++;
			/* Update charges */
			pcfg->membraneset[i].charge=(pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>6) & 0x0003;
			pcfg->membraneset[next_id].charge=(pcfg->rodruleset[rule * num_labels * NUMBER_OF_CHARGES + label * NUMBER_OF_CHARGES + charge].rinfo>>8) & 0x0003;
			num_exec++;
                break;
		};

		time_ex+=end_timer();
	} /* End for i (membranes) */

	start_timer();

	for (int i=0;i<num_objects;i++)
		pcfg->skin_multiset[i]+=skin_aux[i];

	time_ex+=end_timer();

	return num_exec;
}

void fast_sequential(Configuration * cfg, Pcuda_configuration pcfg, Problem_size ps) {
    bool end = false;
    int skinExec,exec;
    float timetotal = 0.0, timestep = 0.0;
    float sel_total=0.0,ex_total=0.0;
    float time_malloc=0.0;
    int verbose=pcfg->verbose;

    init_time();
    
    start_timer();
    selected_evol = new ushort [ps->num_objects];
    selected_times = new ushort [ps->num_objects];
    skin_aux = new uint [ps->num_objects];
    time_malloc=end_timer();

    for (uint curStep=pcfg->current_step; ((!end) && (curStep<ps->step_limit)); curStep++) {
	time_sel=time_ex=0.0;
	uint prev_num_membranes=pcfg->max_num_membranes;
	
	if (verbose > 1) {
		cout << "**************************************" << endl;
		cout << "CONFIGURATION " << curStep << ": " << pcfg->max_num_membranes << " membranes:" << endl;

		print_environment(cfg,pcfg,ps);
		print_skin(cfg, pcfg, ps);
		print_multisets(cfg, pcfg, ps);
	}

        //start_timer();
        select_skin(cfg,pcfg);	
	//time_sel+=end_timer();

        exec=execute_step(cfg,pcfg,ps);

	//start_timer();	
        skinExec=execute_skin(cfg,pcfg);	
	//time_ex+=end_timer();

        //timestep=end_timer();

        end = (skinExec==0 && exec==0);
        sel_total+=time_sel;
        ex_total+=time_ex;
        timestep=time_sel+time_ex;
        timetotal+=timestep;

	if (verbose>0) {
		cout << "**************************************" << endl;
		cout << "CONFIGURATION " << curStep << ": from " << prev_num_membranes << " to " << pcfg->max_num_membranes << " membranes:" << endl;
		cout << "TIME: sel= " << time_sel << " , exec= " << time_ex  << " , step total= " << timestep << endl;
	}
    }

    if (verbose>0) {
    print_environment(cfg, pcfg, ps);
    print_skin(cfg, pcfg, ps);}
    //print_multisets(cfg, pcfg, ps);

    cout << "**************************************" << endl;
    cout << "FAST SEQUENTIAL TIME: " << endl;
    cout << "\tMalloc time: " << time_malloc << " ms" << endl;
    cout << "\tSelection time = "<< sel_total <<" ms" << endl;
    cout << "\tExecution time = " <<ex_total << " ms" << endl;
    cout << "\tTotal time = " << timetotal+time_malloc << " ms" << endl;

}



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
 * sequential.cpp
 *
 *  Created on: 28-dic-2008
 *      Author: miguel
 */
#include "sequential.h"
//#include "timestat/timestat.h"


/********************************************************************/
/********* PRINTING FUNCTIONS *********/

void print_ini(int verbose_level, Configuration* cfg) {
	if (verbose_level<3) return;

	cout << "ALPHABET (object, id): ";
	for (int i=0;i<cfg->get_alphabet_size();i++)
		cout<<"("<<cfg->get_alphabet()[i]<<","<<i<<"), ";
	cout << endl;

	cout << "LABELS (label, id): ";
	for (int i=0;i<cfg->get_label_set_size();i++)
		cout<<"("<<cfg->get_labels()[i]<<","<<i<<"), ";
	cout << endl;

	cout << cfg->get_rules()->to_string(cfg->get_labels(),cfg->get_alphabet());
}

void print_membrane(FILE* fout,Membrane* mmbr,Configuration * cfg) {
	if (mmbr==NULL || cfg==NULL) return;

	fprintf(fout,"\nMEMBRANE ID: %d, Label: %s, Charge: %d\n",mmbr->id,cfg->get_labels()[mmbr->label],mmbr->charge);
	fprintf(fout,"Multiset: %s\n",mmbr->multiset.to_string(cfg->get_alphabet()).c_str());
	if (mmbr->mother!=NULL)
		fprintf(fout,"Parent membrane ID: %d\n",mmbr->mother->id);


	if (mmbr->first_child!=NULL)
		print_membrane(fout,mmbr->first_child,cfg);
	if ((mmbr->mother!=NULL) && (mmbr->mother->first_child->id!=mmbr->next_sister->id))
		print_membrane(fout,mmbr->next_sister,cfg);
}

void print_selection(int verbose_level,Configuration* cfg) {
	if (verbose_level<3) return;
	cout << "----------------START SELECTED RULES---------------------" << endl;;
	cout << "STEP " << cfg->get_configuration_number() << endl;
	cout << cfg->get_selection()->to_string(cfg->get_labels(),cfg->get_alphabet());
	cout << "----------------END SELECTED RULES---------------------" << endl;
}

void print_info(int verbose_level, Configuration * cfg, bool end) {
	FILE* fout=stdout; /* TODO: For a future file output */

	if (cfg==NULL || verbose_level==0 || (verbose_level==1 && !end))
		return;

	fprintf(fout,"\n***********************************************\n\n");
	fprintf(fout,"CONFIGURATION: %d\n",cfg->get_configuration_number());
	//print_time(fout);
	//fprintf(fout,"MEMORY: %d KB\n",mem_total());
	fprintf(fout,"\n");

	if (cfg->get_environment()!=NULL) {
		fprintf(fout,"ENVIRONMENT: %s\n\n",cfg->get_environment()->to_string(cfg->get_alphabet()).c_str());
	}

	fprintf(fout,"SKIN MEMBRANE ID: %d, Label: %s, Charge: %d\n", cfg->get_skin()->id, cfg->get_labels()[cfg->get_skin()->label], cfg->get_skin()->charge);
	fprintf(fout,"Multiset: %s\n", cfg->get_skin()->multiset.to_string(cfg->get_alphabet()).c_str());
	fprintf(fout,"Internal membranes count: %d\n", cfg->get_number_membranes()-1);

	if (verbose_level==3 || end)
		print_membrane(fout,cfg->get_skin()->first_child,cfg);
	fflush(stdout);
}
/*********************************************************************/

void sequential_solve(Configuration *c0, int verbose_level, int step_limit, int max_membranes, int num_objects, int block_size, int threshold, int mode) {
	bool end=false;
	Configuration *cfg=c0;
	int i=0;
	bool onlyseq=true;
	time_t total_sec=0;
	double time_sel=0.0,time_exec=0.0;
	double time_sel_stp=0.0, time_exec_stp=0.0;
	/* Initialize the time counter */
	print_ini(verbose_level,cfg);

	init_time();

	while ((!end) && (i++<step_limit)) {
		print_info(verbose_level,cfg,end);
	  /*print_info(verbose_level,cfg,end,total_sec);*/
		/* If the number of membranes is higher of the threshold, then, call to
		 * the pcuda_solve, parallelization on the GPU */
		//cout << "Paso " << i << " con " << cfg->get_number_membranes() << " membranas" << endl;
		if (threshold < cfg->get_number_membranes()) {
			onlyseq=false;
			cout << "Threshold achieved in number of membranes: "<< threshold <<". Go to PCUDA in mode " << mode << endl;
			pcuda_solve(cfg,step_limit,i,max_membranes,num_objects,block_size,mode,verbose_level);
			break;
		}

		if (verbose_level>0)
		cout << "Step " << i << " with " << cfg->get_number_membranes() << " membranes" << endl;

		time_sel_stp=0;
		time_exec_stp=0;

		start_timer();

		end=!select_rules(cfg);

		time_sel_stp=end_timer();

		time_sel+=time_sel_stp;

		if (!end) {
			print_selection(verbose_level,cfg);

			start_timer();

			execute_rules(cfg);

			time_exec_stp=end_timer();
			time_exec+=time_exec_stp;
		}

		if (verbose_level>0)
		cout << "Time for step " << time_sel_stp+time_exec_stp << " ms" << endl;
	}

	//print_info(verbose_level,cfg,true);

	if (i-1==step_limit) {
		print_info(verbose_level,cfg,true);
		cout << endl << "Reached step limit: " << step_limit << endl;
	}
	
	if (onlyseq) {
		cout << "**************************************" << endl;
		cout << "SEQUENTIAL TIME: " << endl;
		cout << "\tSelection time = "<< time_sel <<" ms" << endl;
		cout << "\tExecution time = " << time_exec << " ms" << endl;
		cout << "\tTotal time = " << time_sel+time_exec << " ms" << endl;
	}
	
	total_sec=get_seconds();
	cout << endl << "TOTAL TIME: " << total_sec << " s" << endl;
	cout << "SEQUENTIAL: SELECTION TIME " << time_sel << "ms, EXECUTION TIME: " << time_exec << "ms" << endl;
}

void select(Membrane *current, Ruleset* rule_set, Selectedrulelist* rule_selection) {
	int mult=0;
	Rulelist * rulelist=NULL;
	Rule * rule=NULL;
	bool rule_end=false;

  	if ((current==NULL) || (rule_set==NULL) || (rule_selection==NULL))
		return;

	rulelist=rule_set->get_rulelist(current->label,current->charge);

	rulelist->start_iteration();
	while (!rulelist->end()) {
		rule=rulelist->next_rule();
		if (rule==NULL) break;

		/* [ a --> U ] */
		if ( (rule->type == RULE_EVOLUTION) && (current->multiset.check_object(rule->a) ) ) {
			mult=current->multiset.consume_all(rule->a);
			if (mult>0)
				rule_selection->add_selected_rule(rule,current,mult);
		}
		/* a[] --> [b] */
		else if ( !rule_end && (rule->type == RULE_SEND_IN) && (current->id != SKIN_ID)
				&& (current->mother->multiset.check_object(rule->a)) ) {
			mult=current->mother->multiset.consume_object(rule->a,1);
			if (mult>0) {
				rule_selection->add_selected_rule(rule,current,1);
				rule_end=true;
			}
		}
		/* [a] --> b[] */
		else if ( !rule_end && (rule->type == RULE_SEND_OUT)	&& (current->multiset.check_object(rule->a)) ) {
			mult=current->multiset.consume_object(rule->a,1);
			if (mult>0) {
				rule_selection->add_selected_rule(rule,current,1);
				rule_end=true;
			}
		}
		/* [a] --> [b][c] */
		else if ( !rule_end && (rule->type == RULE_DIVISION) && (current->id != SKIN_ID)
			&& (current->first_child == NULL) && (current->multiset.check_object(rule->a)) ) {
			mult=current->multiset.consume_object(rule->a,1);
			if (mult>0) {
				rule_selection->add_selected_rule(rule,current,1);
				rule_end=true;
			}
		}
		/* [a] --> b */
		if ( !rule_end && (rule->type == RULE_DISSOLUTION) && (current->id != SKIN_ID)
			&& (current->multiset.check_object(rule->a)) ) {
			mult=current->multiset.consume_object(rule->a,1);
			if (mult>0) {
				rule_selection->add_selected_rule(rule,current,1);
				rule_end=true;
			}
		}
	} /* end while of rulelist */

	if (rulelist!=NULL) delete rulelist;
}

bool select_rules(Configuration *cfg) {
	Membrane * current;

	if (cfg==NULL) return false;

	current=cfg->get_skin();
	cfg->new_selection();

	/* Start the selection with the skin */
	do {
		/* Do the selection */
		select(current,cfg->get_rules(),cfg->get_selection());

		/* Go first to the childrens: depth search */
		if (current->first_child!=NULL) {
			current=current->first_child;
		}
		/* If not reached the end of the list of sisters */
		else if ((current->mother != NULL) &&
				(current->next_sister->id!=current->mother->first_child->id)) {
			current=current->next_sister;
		}
		else {
			/* Go back to high levels */
			while ( (current != NULL) && (current->id != SKIN_ID) &&
					(current->next_sister->id==current->mother->first_child->id) ) {
				current=current->mother;
			}
			/* If not reached the top of the tree, go to the next sister */
			if ((current != NULL) && (current->id != SKIN_ID))
				current=current->next_sister;
		}
	} while ( (current != NULL) && (current->id != SKIN_ID) );

	return (cfg->get_selection()->get_length() > 0);
}

int execute_rules(Configuration *cfg) {
	Rule* rule=NULL;
	Membrane* membrane=NULL,* new_membrane=NULL;
	Membrane *ch1=NULL,*ch2=NULL;
	int count=0;

	if (cfg==NULL) return -1;

	cfg->get_selection()->start_iteration();
	while ( ! cfg->get_selection()->end() ) {
		rule=cfg->get_selection()->get_rule();
		membrane=cfg->get_selection()->get_membrane();
		count=cfg->get_selection()->get_times();

		cfg->get_selection()->next_selection();

		switch ( rule->type ) {
		case RULE_EVOLUTION: /* [ a --> U ] */
			membrane->multiset.add_multiset(*(rule->new_multiset),count);
		break;
		case RULE_SEND_IN: /* a[] --> [b] */
			membrane->multiset.add_object(rule->b);
			membrane->charge=rule->new_charge;
		break;
		case RULE_SEND_OUT: /* [a] --> b[] */
			if (membrane->id == SKIN_ID)
				cfg->get_environment()->add_object(rule->b);
			else
				membrane->mother->multiset.add_object(rule->b);
			membrane->charge=rule->new_charge;
		break;
		case RULE_DIVISION: /* [a] --> [b][c] */
			/* Creation of a new membrane as a new sister */
			new_membrane = new Membrane;
			new_membrane->next_sister=membrane->next_sister;
			new_membrane->prev_sister=membrane;
			membrane->next_sister->prev_sister=new_membrane;
			membrane->next_sister=new_membrane;
		      	new_membrane->mother=membrane->mother;
			new_membrane->label=membrane->label;
			new_membrane->first_child=NULL;
			new_membrane->id=cfg->get_new_membrane_id();

			/* Copy the information */
			new_membrane->multiset=membrane->multiset;

			/* Add new objects and set new charges */
			membrane->multiset.add_object(rule->b,1);
			new_membrane->multiset.add_object(rule->c,1);
			membrane->charge=rule->new_charge;
			new_membrane->charge=rule->new_charge_2;
			new_membrane=NULL;
			cfg->increment_number_membranes();
		break;
		case RULE_DISSOLUTION: /* [a] --> b */
			/* Delete the membrane */

			/* Add the multiset to his mother */
			membrane->mother->multiset.add_multiset(membrane->multiset);
			/* Add the new object to the mother */
			membrane->mother->multiset.add_object(rule->b,1);

			/* Fix the next_sister pointer */
			membrane->next_sister->prev_sister=membrane->prev_sister;
			membrane->prev_sister->next_sister=membrane->next_sister;

			/* If the membrane to be deleted has children */
			if (membrane->first_child != NULL) {
				ch1=membrane->first_child;
				ch2=ch1;
				/* Fix children's mother */
				while (ch2->next_sister->id != membrane->first_child->id) {
					ch2->mother=membrane->mother;
					ch2=ch2->next_sister;
				}
				ch2->mother=membrane->mother;
			}

			/* Check if it is the last membrane of its level */
			if (membrane->id == membrane->next_sister->id) {
				/* NULL if there is no new childrens (ch1 is NULL if no childrens)
				 * But if there are childrens, ch1 points to the first children */
				membrane->mother->first_child=ch1;
			}
			/* If there are more membrane sisters, and there are childrens, add these
			 * nodes to the level */
			else if (ch1!=NULL) {
				membrane->prev_sister->next_sister=ch1;
				ch1->prev_sister=membrane->prev_sister;
				ch2->next_sister=membrane->next_sister;
				membrane->next_sister->prev_sister=ch2;
			}
			/* If the deleted node is the first child of his mother, change it to the previous one */
			else if (membrane->id == membrane->mother->first_child->id) {
			        membrane->mother->first_child=membrane->next_sister;
			}
			ch1=ch2=NULL;

			delete membrane;

			cfg->decrement_number_membranes();

		break;
		}
	}

	cfg->delete_selection();

	/* The new configuration is ready */
	cfg->increment_configuration_number();


	return 0;
}

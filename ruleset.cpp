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
 * ruleset.cpp
 *
 *  Created on: 23-dic-2008
 *      Author: miguel
 */
#include "ruleset.h"

std::string rule_to_string (Rule* rule, int labelid, int charge, char ** label_def, char ** object_def) {
	std::string out="";
	char aux[10];
	char charges[]= {' ','+','-'};
	int i=labelid, j=charge;

	if (rule != NULL && label_def!=NULL && object_def!=NULL) {
		sprintf(aux," ID%d: ",rule->id);
		out+=aux;
		switch (rule->type) {
		case RULE_EVOLUTION:
			out+=charges[j]; out+="[";
			out+=object_def[rule->a];
			out+=" -> "; out+=rule->new_multiset->to_string(object_def);
			out+="]'"; out+=label_def[i];
			break;
		case RULE_SEND_IN:
			out+=object_def[rule->a];
			out+=charges[j]; out+="[]'"; out+=label_def[i];
			out+=" -> "; out+=charges[rule->new_charge];
			out+="["; out+=object_def[rule->b]; out+="]";
			break;
		case RULE_SEND_OUT:
			out+=charges[j]; out+="[";
			out+=object_def[rule->a]; out+="]'";
			out+=label_def[i]; out+=" -> ";
			out+=charges[rule->new_charge]; out+="[]";
			out+=object_def[rule->b];
			break;
		case RULE_DISSOLUTION:
			out+=charges[j]; out+="[";
			out+=object_def[rule->a]; out+="]'";
			out+=label_def[i]; out+=" -> ";
			out+=object_def[rule->b];
			break;
		case RULE_DIVISION:
			out+=charges[j]; out+="[";
			out+=object_def[rule->a]; out+="]'";
			out+=label_def[i]; out+=" -> ";
			out+=charges[rule->new_charge]; out+="[";
			out+=object_def[rule->b]; out+="] ";
			out+=charges[rule->new_charge_2]; out+="[";
			out+=object_def[rule->c]; out+="]";
		} /* switch */
	}
	else {
		/* No printing yet */
	} /* if char *** */
	return out;
}
/****************************/
/* Methods of class Ruleset */

Ruleset::~Ruleset() {
	Node*p,*p2;

	for (int i=0;i < this->label_set_size;i++) {
		for (int j=0; j < NUMBER_OF_CHARGES; j++) {
			p=this->rule_sort_set[i][j];
			while (p) {
				p2=p->next_rule;
				if (p->info!=NULL) {
					delete p->info->new_multiset;
					delete p->info;
				}
				delete p;
				p=p2;
			}
		}
	}
}

/* Auxiliary methods for handling rules */
bool Ruleset::set_info(Rule* rule, short int type, ObjectID a, ObjectID b, ObjectID c,
	short int new_charge, short int new_charge_2, const Multiset* new_multiset) {
		if (rule==NULL) return false;
		rule->type=type;
		rule->a=a;
		rule->b=b;
		rule->c=c;
		rule->new_charge=new_charge;
		rule->new_charge_2=new_charge_2;
		if (new_multiset!=NULL) {
			rule->new_multiset=new Multiset;
			*(rule->new_multiset)=*new_multiset;
		}
		else
			rule->new_multiset=NULL;

		rule->id=++id_count;
		return true;
}

bool Ruleset::add_sort(Rule * rule,Node* list){
	Node*p=NULL, *p2=NULL;
	//short int prio1,prio2;

	if ((list==NULL)||(rule==NULL)) return false;

	p=list->next_rule;
	p2=list;
	//prio1=rule_prio[rule->type];
	//if (p) prio2=rule_prio[p->info->type];

	/* Sort the lists by type and object */
	while ((p) && ((p->info->type > rule->type) ||
			((p->info->type == rule->type)&&(p->info->a <= rule->a)))) {
	//while ((p) && ((prio2 < prio1) ||
	//		((prio2 == prio1)&&(p->info->a <= rule->a)))) {
		p=p->next_rule;
		p2=p2->next_rule;
		//if (p) prio2=rule_prio[p->info->type];
	}

	p2->next_rule=new Node;
	p2->next_rule->info=rule;
	p2->next_rule->next_rule=p;
	p=p2=NULL;
	numrules++;
	return true;
}

/* Methods for adding rules */
bool Ruleset::add_evolution_rule(LabelID label, short int charge, ObjectID a, Multiset& multiset) {
	Rule* new_rule=new Rule;

	if (charge<0 || charge>=NUMBER_OF_CHARGES || label >= this->label_set_size) {
		return false;
	}

	evol_length+=multiset.get_length();

	if ( set_info(new_rule,RULE_EVOLUTION,a,0,0,-1,-1,&multiset) &&
	add_sort(new_rule,rule_sort_set[label][charge]) )
		return true;
	else {
		delete new_rule;
		return false;
	}
	return true;
}

bool Ruleset::add_sendin_rule(LabelID label, short int charge, short int new_charge, ObjectID a, ObjectID b){
	Rule* new_rule=new Rule;

	if (charge<0 || charge>=NUMBER_OF_CHARGES || new_charge<0 || new_charge>=NUMBER_OF_CHARGES
			|| label >= this->label_set_size) {
		return false;
	}

	if ( set_info(new_rule,RULE_SEND_IN,a,b,0,new_charge,-1,NULL) &&
	add_sort(new_rule,rule_sort_set[label][charge]) )
		return true;
	else {
		delete new_rule;
		return false;
	}
	return true;
}

bool Ruleset::add_sendout_rule(LabelID label, short int charge, short int new_charge, ObjectID a, ObjectID b){
	Rule* new_rule=new Rule;

	if (charge<0 || charge>=NUMBER_OF_CHARGES || new_charge<0 || new_charge>=NUMBER_OF_CHARGES
			|| label >= this->label_set_size) {
		return false;
	}

	if ( set_info(new_rule,RULE_SEND_OUT,a,b,0,new_charge,-1,NULL) &&
	add_sort(new_rule,rule_sort_set[label][charge]) )
		return true;
	else {
		delete new_rule;
		return false;
	}
	return true;
}

bool Ruleset::add_disolution_rule(LabelID label, short int charge, ObjectID a, ObjectID b){
	Rule* new_rule=new Rule;

	if (charge<0 || charge>=NUMBER_OF_CHARGES || label >= this->label_set_size) {
		return false;
	}

	if ( set_info(new_rule,RULE_DISSOLUTION,a,b,0,-1,-1,NULL) &&
	add_sort(new_rule,rule_sort_set[label][charge]) )
		return true;
	else {
		delete new_rule;
		return false;
	}
	return true;
}

bool Ruleset::add_division_rule(LabelID label, short int charge, short int new_charge, short int new_charge_2, ObjectID a, ObjectID b, ObjectID c){
	Rule* new_rule=new Rule;

	if (charge<0 || charge>=NUMBER_OF_CHARGES || new_charge<0 || new_charge>=NUMBER_OF_CHARGES
		|| new_charge_2<0 || new_charge_2>=NUMBER_OF_CHARGES || label >= this->label_set_size) {
		return false;
	}

	if ( set_info(new_rule,RULE_DIVISION,a,b,c,new_charge,new_charge_2,NULL) &&
	add_sort(new_rule,rule_sort_set[label][charge]) )
		return true;
	else {
		delete new_rule;
		return false;
	}
	return true;
}

Rulelist * Ruleset::get_rulelist(LabelID label, short int charge) {
	if (charge<0 || charge>=NUMBER_OF_CHARGES || label >= this->label_set_size)
		return NULL;
	return new Rulelist(this->rule_sort_set[label][charge]);
}

std::string Ruleset::to_string(char ** label_def, char ** object_def) {
	std::string out="";
	Rule * rule=NULL;
	Rulelist* rl=NULL;

	for (int i=0; i<label_set_size; i++) {
		for (int j=0; j<NUMBER_OF_CHARGES; j++) {
			rl=this->get_rulelist(i,j);
			while (!rl->end()) {
				rule=rl->next_rule();
				if (rule!=NULL) {
					out+="Info rule: ";
					out+=rule_to_string(rule,i,j,label_def,object_def);
					out+="\n";
				} /* if rule !=NULL */
			} /* while rulelist */
		} /* for j */
	} /* for i */
	return out;
}

/*****************************/
/* Methods of class Rulelist */
bool Rulelist::start_iteration(){
	iterator=head->next_rule;
	return true;
}

Rule* Rulelist::next_rule(){
	if (iterator != NULL) {
		Rule* r=iterator->info;
		iterator=iterator->next_rule;
		return r;
	}
	return NULL;
}

bool Rulelist::end(){
	return (iterator==NULL);
}


/*************************************/
/* Methods of class Selectedrulelist */

bool Selectedrulelist::add_selected_rule(Rule* rule, Membrane* membrane, int times) {
	if (iteration_mode || rule == NULL || membrane == NULL) {
		return false;
	}

	last->next_node=new Node;
	last=last->next_node;

	last->rule=rule;
	last->membrane=membrane;
	last->selected_times=times;
	last->next_node=NULL;
	length++;
	return true;
}

bool Selectedrulelist::start_iteration(){
	if (! iteration_mode) {
		iteration_mode=true;
		iterator=head->next_node;
		return true;
	}
	return false;
}

bool Selectedrulelist::next_selection(){
	if (iteration_mode && (iterator != NULL)) {
		iterator=iterator->next_node;
		return true;
	}
	return false;
}

bool Selectedrulelist::end() {
	return (iteration_mode && (iterator == NULL));
}

bool Selectedrulelist::end_iteration() {
	if (iteration_mode) {
		iteration_mode=false;
		iterator=head;
		return true;
	}
	return false;
}

Rule* Selectedrulelist::get_rule() {
	if (iteration_mode && (iterator != NULL)) {
		return iterator->rule;
	}
	return NULL;
}

Membrane* Selectedrulelist::get_membrane() {
	if (iteration_mode && (iterator != NULL)) {
		return iterator->membrane;
	}
	return NULL;
}

int Selectedrulelist::get_times() {
	if (iteration_mode && (iterator != NULL)) {
		return iterator->selected_times;
	}
	return -1;
}

std::string Selectedrulelist::to_string(char ** label_def, char ** object_def) {
	std::string out;
	char aux[50];
	Node* it=this->head->next_node;

	while (it!=NULL) {
		out+="Rule: ";
		out+=rule_to_string(it->rule,it->membrane->label,it->membrane->charge,label_def,object_def);
		out+=", Membrane id: ";
		sprintf(aux,"%d",it->membrane->id);
		out+=aux;
		out+=", Times: ";
		sprintf(aux,"%d\n",it->selected_times);
		out+=aux;
		it=it->next_node;
	}

	return out;

}

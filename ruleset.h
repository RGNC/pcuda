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
 * ruleset.h
 *
 *	Define classes for rule management: a rule sorted set, a rule selection
 *	list, and a rule iterator
 *  Created on: 17-dic-2008
 *      Author: miguel
 */

#ifndef RULESET_H_
#define RULESET_H_
#include "binvalues.h"
#include "multiset.h"
#include "membrane.h"
#include "p_system.h"

/*****************************************/

/* Abstract data structures for sequential simulation */
struct _rule {
  short int type;
  ObjectID a,b,c;
  short int new_charge, new_charge_2;
  Multiset* new_multiset;
  RuleID id;
};

typedef struct _rule Rule;

std::string rule_to_string (Rule* rule, int labelid, int charge, char ** label_def, char ** object_def);

/***********************************************************************/
/* CLASS Rulelist: Stores a list of rules sorted by type and object, for the
 * same label and charge */
class Rulelist;

/***********************************************************************/
/* CLASS Ruleset: Stores the whole set of rules sorted by label, charge,
 * type of rule and object */
class Ruleset {
	friend class Rulelist;

private:
	LabelID label_set_size;
	int numrules;
	RuleID id_count;
	int evol_length; // A counter for the total length of the evolution's multisets

	/* A sort list of the rules */
	struct _node_rule {
		struct _rule * info;
		struct _node_rule * next_rule;
	};

	typedef struct _node_rule Node;

	struct _node_rule*** rule_sort_set; /* The rules are sorted by: label,
						 charge ([3]) and then, type and object */

	bool set_info(Rule* rule, short int type, ObjectID a, ObjectID b, ObjectID c,
			short int new_charge, short int new_charge_2, const Multiset* new_multiset);

	bool add_sort(Rule * rule,Node* list);

public:
	Ruleset (){
		rule_sort_set=NULL;
		label_set_size=0;
		numrules=0;
		id_count=EMPTY_RULE;
		evol_length=0;
	}

	Ruleset(int variant,LabelID label_set_size) {
		this->label_set_size=label_set_size;
		numrules=0;
		id_count=0;
		/* TODO:
		 * - Check the variant and configure other parameters. i.e.: a global file */
		variant=0;
		variant++;
		rule_sort_set=new Node**[label_set_size];
		for (int i=0;i<label_set_size; i++) {
			rule_sort_set[i]=new Node*[NUMBER_OF_CHARGES];
			for (int j=0;j<NUMBER_OF_CHARGES;j++) {
				rule_sort_set[i][j]=new Node;
				rule_sort_set[i][j]->info=NULL;
				rule_sort_set[i][j]->next_rule=NULL;
			}
		}
	}

	~Ruleset();

	int get_num_rules() { return numrules; }
	int get_evol_length() { return evol_length; }

	bool add_evolution_rule(LabelID label, short int charge, ObjectID a, Multiset& multiset) ;
	bool add_sendin_rule(LabelID label, short int charge, short int new_charge, ObjectID a, ObjectID b);
	bool add_sendout_rule(LabelID label, short int charge, short int new_charge, ObjectID a, ObjectID b);
	bool add_disolution_rule(LabelID label, short int charge, ObjectID a, ObjectID b);
	bool add_division_rule(LabelID label, short int charge, short int new_charge, short int new_charge_2, ObjectID a, ObjectID b, ObjectID c);

	Rulelist * get_rulelist(LabelID label, short int charge);

	std::string to_string(char ** label_def=NULL, char ** object_def=NULL);
};

/****************************************************************************/
/* CLASS Rulelist: Stores a list of rules sorted by type and object, for the
 * same label and charge */
class Rulelist{
friend class Ruleset;
private:
	Ruleset::Node * head, * iterator;

	Rulelist(Ruleset::Node * head) {
		this->head=head;
		iterator=head;
	}

public:
	~Rulelist() {
		head=iterator=NULL;
	}
	bool start_iteration(); /* Initialize the iteration */
	Rule * next_rule();  /* Returns the current rule, and point to the next one
						 for the next call to this method */
	bool end();  /* True if reached the end of the list */
};


/***********************************************************************************/
/* CLASS Selectedrulelist: Stores the rules selected in each step of the algorithm */
class Selectedrulelist {
private:
	struct _node {
		Rule* rule;
		Membrane* membrane;
		int selected_times;
		struct _node* next_node;
	};

	typedef struct _node Node;

	Node* head,* last;
	Node* iterator;
	int length; /* The length of the list */
	bool iteration_mode; /* Sets the iteration mode (no adding is allowed) */

public:
	Selectedrulelist() {
		head=new Node;
		head->next_node=NULL;
		head->rule=NULL;
		head->membrane=NULL;
		iterator=last=head;
		length=0;
		iteration_mode=false;
	}

	~Selectedrulelist() {
		while (head) {
			iterator=head->next_node;
			delete head;
			head=iterator;
		}
	}

	int get_length() {
		return this->length;
	}

	/* Adding method */
	bool add_selected_rule(Rule* rule, Membrane* membrane, int times);

	/* Iteration methods */
	bool start_iteration(); /* Starts the iteration mode, the first element is available if exists */
	bool next_selection();  /* Goes to the next item */
	bool end();  /* True if reached the end of the list */
	bool end_iteration();	/* Ends the iteration mode */
	Rule* get_rule();	/* Returns the Rule of the current item */
	Membrane* get_membrane();	/* Returns the Membrane of the current item */
	int get_times();	/* Returns the number of selected times of the current item */

	std::string to_string(char ** label_def=NULL, char ** object_def=NULL);
};

#endif /* RULESET_H_ */

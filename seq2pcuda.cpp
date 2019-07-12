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
 * seq2pcuda.cpp
 *
 *  Created on: 24-feb-2009
 *      Author: miguel
 */

#include "seq2pcuda.h"

size_t get_max_memory(int num_objects, int num_labels, int max_membranes) {
	size_t rod=num_objects*num_labels*NUMBER_OF_CHARGES*sizeof(_rodruleset)/1048576;
	size_t mmbr=max_membranes*sizeof(_membraneset)/1048576;
    	size_t multiselecSize = max_membranes * num_objects * sizeof(_multiselec)/1048576;
    	size_t skinSize = num_objects * sizeof(_skin_multiset)/1048576;
	size_t env=num_objects * sizeof(_skin_multiset)/1048576;
    	size_t sioddSize = max_membranes * sizeof(_selec_siodd)/1048576;
	return rod+mmbr+2*multiselecSize+skinSize+env+sioddSize;
}

/* Creates the "fishing rod" rule set from the rule set of
 * the sequential version */
Rodruleset get_rodruleset(Ruleset& ruleset, uint & rod_length, uint & evol_length, Rodmultiset * evol_table, int alphabet_size, int label_set_size) {
	unsigned int length=alphabet_size*label_set_size*NUMBER_OF_CHARGES;
	unsigned int length_per_object=label_set_size*NUMBER_OF_CHARGES;
	unsigned int pos;
	int num_rules=0,begin=0,len=0;
	short int soddtype=0;

	/* Initializes the rod rule set*/
	Rodruleset cr=new struct _rod_rule[length];
	Rule * rule=NULL;
	Rulelist* rl=NULL;

	num_rules=ruleset.get_num_rules();
	evol_length=ruleset.get_evol_length();

	if (num_rules<0 || evol_length<0) return NULL;

	/* see that the parameter rule_table is a pointer to an array of pointers, in order
	 * to be initalized in this function*/
	if (evol_table!=NULL) {
		*evol_table=new _rodmultiset[evol_length];
	}

	rod_length=length;

	for (uint i=0;i<length;i++) {
		cr[i].rev=cr[i].rsodd[0]=cr[i].rsodd[1]=cr[i].rsin=cr[i].rinfo=0;
	}

	/* For each label and charge, it iterates the rulelist associated*/
	for (int i=0; i<label_set_size; i++) {
		for (int j=0; j<NUMBER_OF_CHARGES; j++) {
			rl=ruleset.get_rulelist(i,j);
			if (!rl->start_iteration()) continue;
			while (!rl->end()) {
				/* Extract the rule and the iterator goes to the next one */
				rule=rl->next_rule();
				if (rule!=NULL) {
					/* place of the rule in the rodruleset*/
					pos=rule->a*length_per_object + i*NUMBER_OF_CHARGES + j;
					bool isevol=cr[pos].rinfo&0x1, issodd = cr[pos].rinfo&0x2, issin=cr[pos].rinfo&0x4;

					/* if evolution rule*/
					if ((isevol)==0 && rule->type==RULE_EVOLUTION) {
						// Set the evolution rule flag
						cr[pos].rinfo|=0x1;
						// Add the multiset of the evolution rule
						len=rule->new_multiset->to_array(&((*evol_table)[begin]));
						// Add begin in the low 20 bits, and len in the 12 high bits
						cr[pos].rev=cr[pos].rev&0xFFF00000 | begin;
						cr[pos].rev=cr[pos].rev&0x000FFFFF | len << 20;
						begin+=len;
						//std::cout << std::hex << cr[pos].rinfo<<std::endl;
					}
					/* if other type but not send-in (use local object) */
					else if (rule->type!=RULE_EVOLUTION && rule->type!=RULE_SEND_IN && rule->id!=EMPTY_RULE) {
						// Extract the soddtype
						soddtype=(cr[pos].rinfo&0x0034) >> 3;
						if (soddtype<rule->type) {
							// Set the sout/div/dis rule (rsodd) flag
							cr[pos].rinfo|=0x2;
							// Add soddtype
							cr[pos].rinfo=cr[pos].rinfo&0xFFC7 | (rule->type&0x7)<<3;
							// Add charge 1
							cr[pos].rinfo=cr[pos].rinfo&0xFF3F | (rule->new_charge&0x3)<<6;
							// Add charge 2
							cr[pos].rinfo=cr[pos].rinfo&0xFCFF | (rule->new_charge_2&0x3)<<8;
							// Add objects
							cr[pos].rsodd[0]=rule->b;
							cr[pos].rsodd[1]=rule->c;
						}
					}
					/* if send-in (needs a object from the mother) */
					else if (rule->type==RULE_SEND_IN && rule->id!=EMPTY_RULE && issin == 0) {
						// Set the evolution rule flag
						cr[pos].rinfo|=0x4;
						cr[pos].rinfo=cr[pos].rinfo&0xF3FF | (rule->new_charge&0x3)<<10;
						cr[pos].rsin=rule->b;
					}

					//cr[rule->a*length_per_object + i*NUMBER_OF_CHARGES + j].id[rule->type]=rule->id;
				} /* if rule !=NULL */
			} /* while rulelist */
		} /* for j */
	} /* for i */

	return cr;
}

/* Creates the membrane set from the set of membranes of the sequential version */
Membraneset get_membraneset(Membrane* skin, uint & size, uint &num_membranes, uint max_membranes, int num_levels) {
	Membrane * current=NULL;
	num_membranes=0;
	if (num_levels>2) {
		std::cerr << "The current version cannot handle more than 2 levels of membranes"<< std::endl;
		exit(1);
	}

	Membraneset ms=new _membraneset[max_membranes];
	//ms[0].charge=skin->charge;
	//ms[0].label=skin->label;

	size=max_membranes;

	current=skin->first_child;

	do {
		ms[current->id-1].charge=current->charge;
		ms[current->id-1].label=current->label;
		current=current->next_sister;
		num_membranes++;
	}
	while (current!=NULL && current->id!=skin->first_child->id);

	return ms;
}

/* Initalizes the structure for the input multisets */
ushort * get_multiset_set(Membrane *skin, uint & size, uint max_membranes, int num_levels, int num_objects) {
	Membrane * current=NULL;
	int i=0;
	if (num_levels>2) {
		std::cerr << "The current version cannot handle more than 2 levels of membranes"<< std::endl;
		exit(1);
	}
	else if (num_objects <0){
		std::cerr << "Cannot use a negative number of objects"<< std::endl;
		exit(1);
	}
	ushort * ms=new ushort[max_membranes*num_objects];

	size=max_membranes*num_objects;	

	current=skin->first_child;

	do {
		i=current->id-1;
		current->multiset.to_array(&(ms[i*num_objects]),num_objects);
		current=current->next_sister;
	}
	while (current!=NULL && current->id!=skin->first_child->id);

	return ms;
}


uint * get_skinmultiset(Membrane * skin, uint & size, int num_objects) {
        if (num_objects <0){
                std::cerr << "Cannot use a negative number of objects"<< std::endl;
                exit(1);
        }
        uint * sms=new uint[num_objects];
	memset(sms, 0, num_objects*sizeof(uint));

	skin->multiset.to_array((ushort *)sms,num_objects);
        size=num_objects;

        return sms;

}

Multiselec ini_selecevol(uint & size, uint max_membranes, uint num_objects) {
	size=max_membranes*num_objects;

	Multiselec ms = new _multiselec[size];
	for (uint i=0;i<max_membranes;i++) {
		ms[i].rule=EMPTY_RULE;
		ms[i].multiplicity=0;
	}
	return ms;
}

Selec_siodd ini_selecsiodd(uint & size, uint max_membranes) {
	size=max_membranes;
	Selec_siodd s = new _selec_siodd[max_membranes];
	for (uint i=0; i<max_membranes; i++)
		s[i]=EMPTY_RULE;

	return s;
}



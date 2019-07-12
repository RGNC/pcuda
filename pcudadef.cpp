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


#include <stdlib.h>
#include <iostream>

#include "pcudadef.h"

using namespace std;

/***********************/
/* AUXILIARY FUNCTIONS */

int select_skin (Configuration *cfg, Pcuda_configuration pcfg) {
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

        /* [ a --> U ] */
        if ((rule->type == RULE_EVOLUTION) && (pcfg->skin_multiset[rule->a]>0 )) {
            mult=pcfg->skin_multiset[rule->a];
            pcfg->skin_multiset[rule->a]=0;
            if (mult>0)
                cfg->get_selection()->add_selected_rule(rule,cfg->get_skin(),mult);
        }
        /* [a] --> b[] */
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
            case RULE_EVOLUTION: /* [ a --> U ] */
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
        /* For each element in the multiset */
        for (int j=0; j < num_objects; j++) {     
            if (pcfg->multisets[i*num_objects+j]>0) {
                /*sprintf(convert,"%d",j);
                out+=convert;*/
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


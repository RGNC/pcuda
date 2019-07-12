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


/**
  parserbin.h: Read the binary input file and convert it into structures
  for being used on the simulator
  Define a class for the parser and the class for configurations of P systems
*/
#ifndef PARSERBIN_H_
#define PARSERBIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "binvalues.h"
#include "multiset.h"
#include "ruleset.h"
#include "membrane.h"
using namespace std;

/******************************************************************/
/* Classes: Configuration, Parserbin */

class Configuration;

/*-----------------------------------------------------------------
 * PARSERBIN: Read a binary file and returns a Configuration object
 * that stores all the information inside the file
 */
class Parserbin {
  private:
    Configuration * cfg;
    int object_id_size,label_id_size,mmbr_id_size;
    bool short_object_mode,short_label_mode,short_mmbr_mode;

    /**
      read_head: Reads the header, and returns the last 2 Bytes
    */
    unsigned int read_header(ifstream& is);
    /**
      read_alphabet: Reads the alphabet, and returns the number of objects
    */
    int read_alphabet(ifstream& is);
    /**
      read_label_set: Reads the label set, and returns the number of labels
    */
    int read_label_set(ifstream& is);
    /**
      read_membrane_set: Reads the membrane set, and returns the number of membranes
    */
    int read_membrane_set(ifstream& is);
    /**
      read_initial_multiset: Reads the multisets, and returns the number of them
    */
    int read_initial_multiset(ifstream& is,uint number_membranes);
    /**
      read_multiset: Reads a multiset, and returns the number of objects
    */
    int read_multiset(ifstream& is, Multiset& multiset, unsigned int alphabet_size);
    /**
      read_rules: Reads rules depending on the variant, and returns the number them
    */
    int read_rules (ifstream& is, int variant, LabelID number_of_labels);
    /**
      read_evolution_rule: Reads evolution rules, and returns the number them
    */
    int read_evolution_rule(ifstream& is);
    /**
      read_send-in_rule: Reads send-in rules, and returns the number them
    */
    int read_sendin_rule(ifstream& is);
    /**
      read_send-out_rule: Reads send-out rules, and returns the number them
    */
    int read_sendout_rule(ifstream& is);
    /**
      read_disolution_rule: Reads disolution rules, and returns the number them
    */
    int read_disolution_rule(ifstream& is);
    /**
      read_division_rule: Reads division rules, and returns the number them
    */
    int read_division_rule(ifstream& is);

  public:
    Configuration * readfile (const char * file);
};


/*--------------------------------------------------------
 *  CONFIGURATION: store the configuration of the P system
 *  specified in the binary file
 */
class Configuration {
  friend class Parserbin;

  /* Internal data*/
  private:
	int configuration_number;
    unsigned int alphabet_size;
    unsigned int label_set_size;
    char** alphabet; /* Alphabet of the objects (i.e. a,b,c,...) */
    char** label_set; /* All the labels available (i.e. h0,h1,h2,...) */

    Multiset* environment;

    unsigned int membrane_set_size;
    Membrane* skin; /* A pointer to the skin membrane (the head of the tree) */

    /* Handle new membrane ids */
    MembraneID _next_membrane_id;

    /* Auxiliary membrane set */
    struct __membrane_slot {
      int mother, first_child, next_sister, next_childs_sister;
      Membrane * mmbr;
    };
    typedef struct __membrane_slot Membrane_slot;
    Membrane_slot* mmbr_set;

    Ruleset *rule_set;
    Selectedrulelist *rule_selection;

  protected:

    Configuration(){
    	configuration_number=0;
    	membrane_set_size=0;
    	_next_membrane_id=0;
    	skin=NULL;
    	alphabet=NULL;
    	label_set=NULL;
    	rule_set=NULL;
    	rule_selection=NULL;
    	environment=NULL;
    }

  public:

    ~Configuration () {
      /*TODO: delete all the nodes of the membrane tree */
      if (skin != NULL)
    	  delete skin;

      if (alphabet !=NULL)
    	  delete alphabet;
      if (label_set != NULL)
    	  delete label_set;
      if (rule_set != NULL)
    	  delete rule_set;
      if (rule_selection != NULL)
    	  delete rule_selection;
      if (environment != NULL)
    	  delete environment;
    }

    /* Get values of the configuration */
    char** get_alphabet() {
      return alphabet;
    }
    int get_alphabet_size () {
    	return this->alphabet_size;
    }
    char** get_labels() {
      return label_set;
    }
    int get_label_set_size() {
    	return this->label_set_size;
    }
    Membrane* get_skin() {
      return skin;
    }
    Ruleset* get_rules() {
      return rule_set;
    }
    /************************************/

    /* Sets and handle a selection set of rules */
    void new_selection() {
    	if (rule_selection == NULL)
    		delete rule_selection;
    	rule_selection=new Selectedrulelist();
    }
    Selectedrulelist* get_selection() {
      return rule_selection;
    }
    void delete_selection() {
    	if (rule_selection != NULL)
    		delete rule_selection;
    	rule_selection=NULL;
    }
    /********************************************/

    /* Administration of the environment of the configuration */
    void new_environment() {
    	if (environment == NULL)
    		delete environment;
    	environment = new Multiset();
    }
    Multiset * get_environment() {
    	return environment;
    }
    void delete_environment() {
    	if (environment!=NULL)
    		delete environment;
    	environment=NULL;
    }
    /***********************************************/

    /* Administration of the membrane ids */
    void ini_membrane_id() {
    	_next_membrane_id=SKIN_ID;
    }
    MembraneID get_new_membrane_id () {
    	/*TODO: reuse membrane ids of deleted membranes */
    	return _next_membrane_id++;
    }
    int get_number_membranes () {
    	return membrane_set_size;
    }
    int increment_number_membranes() {
    	return ++membrane_set_size;
    }
    int decrement_number_membranes() {
        return --membrane_set_size;
    }

    /***************************/

    /* Administration of the configuration number */
    void set_configuration_number (int new_configuration_number) {
    	configuration_number=new_configuration_number;
    }
    int get_configuration_number () {
    	return configuration_number;
    }
    void increment_configuration_number () {
    	configuration_number++;
    }
    /*********************************************/

};

#endif /* PARSERBIN_H_ */

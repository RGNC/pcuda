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
 * multiset.h
 *
 *	Implementation of a multiset
 *  Created on: 16-dic-2008
 *      Author: miguel
 */

#ifndef MULTISET_H_
#define MULTISET_H_
#include <stdio.h>
#include <string.h>
#include <string>
#include "binvalues.h"
#include "pcuda_types.h"

#define EMPTY_OBJECT_ID 0

class Multiset {
private:
	struct _object {
		ObjectID id;
		unsigned int multiplicity;
		struct _object* next_object;
	};

	typedef struct _object Object;

	Object * head;
	unsigned int length,multiset_size;

	/* Create new and empty internal structures */
	void initialize () {
		head=new Object;
		head->next_object=NULL;
		multiset_size=0;
		length=0;
	}

	/* Delete the internal structures */
	void deleteAll() {
		Object * aux;
		while (head!=NULL) {
			aux=head->next_object;
			delete head;
			head=aux;
		}
		head=NULL;
		length=multiset_size=0;
	}

public:
	Multiset() {
		initialize();
	}
	Multiset(unsigned int multiset_size) {
		initialize();
		this->multiset_size=multiset_size;
	}
	~Multiset() {
		deleteAll();
	}

	/*
	 * Set the size of the multiset (number of different objects in the alphabet)
	 * Returns the size (only change if the size was 0 before)
	 */
	int set_size (unsigned int size) {
		if (multiset_size==0) {
			multiset_size=size;
		}
		return multiset_size;
	}

	unsigned int get_length(); 
	/*
	 * Add an object to the multiset, with the specified multiplicity.
	 * Returns the number of objects added (0 if the object is not from
	 * the multiset).
	 * If the object is the empty object (#, with id 0), doesn't add anything
	 */
	int add_object(ObjectID id, unsigned int multiplicity=1);
	/*
	 * Consume an object from the multiset, with the specified multiplicity.
	 * Returns the number of objects consumed (0 if couldn't consume)
	 */
	int consume_object(ObjectID id, unsigned int multiplicity=1);
	/*
	 * Consume all the instances of the object
	 * Returns the multiplicity that the object had (0 if couldn't consume)
	 */
	int consume_all(ObjectID id);
	/*
	 * Check if there is an object in the multiset
	 * Returns the multiplicity of it (0 if there is not)
	 */
	int check_object(ObjectID id);
        /*
         * Translate the multiset to an array of multiplicities
         * Returns false if num_obj is negative or it's an empty multiset
         */
        bool to_array(uint * mult, uint num_obj);
	/*
	 * Translate the multiset to an array of multiplicities
	 * Returns false if num_obj is negative or it's an empty multiset
	 */
	bool to_array(ushort * mult, uint num_obj);
	/* Translate the multiset to an array of (object, multiplicity) pairs 
	 * Returns the length of the array */
	int to_array(Rodmultiset evol_mult);
         /*
         * Adds the content of the multiset to the array of multiplicites
         * Returns false if num_obj is negative or it's an empty multiset
         */
        bool add_to_array(uint *mult, uint num_obj, uint times);
        /*
	 * Adds the content of the multiset to the array of multiplicites
	 * Returns false if num_obj is negative or it's an empty multiset
	 */
        bool add_to_array(ushort *mult, uint num_obj, uint times);
	/*
	* Adds the content of the multiset to the array of multiplicites
	* Returns false if num_obj is negative or it's an empty multiset
	*/
	bool add_to_array(Multiselec mult, uint num_obj, uint times);
	/*
	 * Add the elements of an input multiset to the multiset, several times
	 */
	void add_multiset(const Multiset &multiset,int times=1);
	/*
	 * Copy one multiset to another
	 */
	Multiset& operator=(const Multiset &multiset);
	/*
	 * Returns a string with the information of the multiset, in order to be printed out
	 */
	std::string to_string(char ** object_description=NULL);
};

#endif /* MULTISET_H_ */

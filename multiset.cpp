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
 * multiset.cpp
 *
 *  Created on: 16-dic-2008
 *      Author: miguel
 */
#include "multiset.h"

int Multiset::add_object(ObjectID id, unsigned int multiplicity) {
	Object * before,*after;
	bool found=false;

	if (id >= this->multiset_size)
		return 0;
	else if (id==0) /* #, empty object */
		return multiplicity;

	before=head;
	after=head->next_object;
	while (!found) {
		/* If we are not at the end of the list, and the current
		 * node is the one that we are looking for */
		if ((after!=NULL)&&(after->id == id)) {
			after->multiplicity+=multiplicity;
			multiplicity=after->multiplicity;
			found=true;
		}
		/* If we are at the end of the list, or if the current object
		 * is higher than the input object (sort list) */
		else if ((after==NULL)||(after->id > id)) {
			before->next_object=new Object;
			before->next_object->next_object=after;
			before->next_object->id=id;
			before->next_object->multiplicity=multiplicity;
			length++;
			found=true;
		}
		else {
			before=after;
			after=after->next_object;
		}
	}

	return multiplicity;

}


int Multiset::consume_object(ObjectID id, unsigned int multiplicity) {
	Object * before,*p;

	if (id >= multiset_size)
		return 0;
	else if (id==0) /* #, empty object */
		return multiplicity;

	before=head;
	p=head->next_object;
	while (p!=NULL) {
		/* If the current object is the input object, and the
		 * multiplicity is enough for the input multiplicity */
		if ((p->id==id)&&(p->multiplicity>multiplicity)) {
			p->multiplicity-=multiplicity;
			return multiplicity;
		}
		/* If the current object is the input object, and the
		 * multiplicity is equal to the input multiplicity */
		else if ((p->id==id)&&(p->multiplicity==multiplicity)) {
			before->next_object=p->next_object;
			delete p;
			length--;
			return multiplicity;
		}
		/* If the current object is the input object, and the
		 * multiplicity is not enough for the input multiplicity */
		else if ((p->id==id)&&(p->multiplicity<multiplicity)) {
			return 0;
		}
		before=p;
		p=p->next_object;
	}
	return 0;
}

int Multiset::consume_all(ObjectID id) {
	Object * before,*p;
	unsigned int multiplicity=0;

	if (id >= multiset_size)
		return 0;
	else if (id==0) /* #, empty object */
		return 1;

	before=head;
	p=head->next_object;
	while (p!=NULL) {
		/* If found the input object */
		if (p->id==id) {
			multiplicity=p->multiplicity;
			before->next_object=p->next_object;
			delete p;
			length--;
			return multiplicity;
		}
		before=p;
		p=p->next_object;
	}
	return 0;
}

int Multiset::check_object(ObjectID id) {
	Object* p=head->next_object;

	if (id >= multiset_size)
		return 0;

	while (p!=NULL) {
		if (p->id==id)
			return p->multiplicity;
		p=p->next_object;
	}
	return 0;
}

bool Multiset::to_array(uint * mult, uint num_obj) {
        Object * p=head->next_object;

        if (multiset_size!=num_obj)
                return false;

        /* For each element in the multiset */
        for (uint i=0; i<multiset_size; i++) {
                mult[i]=0;
                if ((p!=NULL)&&(p->id==i)) {
                        mult[i]=p->multiplicity;
                        p=p->next_object;
                }
        }

        return true;
}

bool Multiset::to_array(ushort * mult, uint num_obj) {
	Object * p=head->next_object;

	if (multiset_size!=num_obj)
		return false;

	/* For each element in the multiset */
	for (uint i=0; i<multiset_size; i++) {
		mult[i]=0;
		if ((p!=NULL)&&(p->id==i)) {
			mult[i]=p->multiplicity;
			p=p->next_object;
		}
	}

	return true;
}

/* Translate the multiset to an array of (object, multiplicity) pairs 
 * Returns the length of the array */
int Multiset::to_array(Rodmultiset evol_mult) {
        Object * p=head->next_object;
	uint i=0;


        /* For each element in the multiset */
        for (; ((i<multiset_size) && (p!=NULL)); p=p->next_object) {
                if (p!=NULL) {
                        evol_mult[i].mult=p->multiplicity;
                        evol_mult[i++].obj=p->id;
                }
        }
        return i;

}

bool Multiset::add_to_array(uint *mult, uint num_obj, uint times) {
        Object * p=head->next_object;

        if (multiset_size!=num_obj)
                return false;

        /* For each element in the multiset */
	while (p!=NULL) {
		if (p->id>0 && p->id < multiset_size) mult[p->id]+=p->multiplicity*times;
		p=p->next_object;
	}

        return true;
}

bool Multiset::add_to_array(ushort *mult, uint num_obj, uint times) {
        Object * p=head->next_object;

        if (multiset_size!=num_obj)
		return false;

        /* For each element in the multiset */
        while (p!=NULL) {
                if (p->id>0 && p->id < multiset_size) mult[p->id]+=p->multiplicity*times;
                p=p->next_object;
	}						        

        return true;
}

bool Multiset::add_to_array(Multiselec mult, uint num_obj, uint times) {
	Object * p=head->next_object;

        if (multiset_size!=num_obj)
                return false;

	while (p!=NULL) {
		if (p->id>0 && p->id < multiset_size) mult[p->id].multiplicity+=p->multiplicity*times;
		p=p->next_object;
	}

        return true;
}

void Multiset::add_multiset(const Multiset &multiset,int times) {
	Object * before,*after,*current;
	ObjectID id;
	bool found=false;
	if (multiset.length>1)
		found=false;

	if (times <= 0)
		return;

	before=head;
	after=head->next_object;
	current=multiset.head->next_object;

	/* Adding each element of the other multiset */
	while (current != NULL) {
		id=current->id;

		while (!found) {
			/* If we are not at the end of the list, and the current
			 * node is the one that we are looking for */
			if ((after!=NULL)&&(after->id == id)) {
				after->multiplicity+=current->multiplicity*times;
				found=true;
			}
			/* If we are at the end of the list, or if the current object
			 * is higher than the input object (sort list) */
			else if ((after==NULL)||(after->id > id)) {
				before->next_object=new Object;
				before->next_object->next_object=after; /* NULL or next object */
				before->next_object->id=id;
				before->next_object->multiplicity=current->multiplicity*times;
				length++;
				found=true;
				before=before->next_object;
			}
			else {
				before=after;
				after=after->next_object;
			}
		}
		found=false;

		current=current->next_object;
	}
}

Multiset& Multiset::operator=(const Multiset &multiset) {
	Object *p1, *p2;
	if (this->length>0) {
		deleteAll();
		initialize();
	}
	this->length=multiset.length;
	this->multiset_size=multiset.multiset_size;

	p1=this->head;
	p2=multiset.head;

	while (p2->next_object!=NULL) {
		p2=p2->next_object;
		p1->next_object=new Object;
		p1=p1->next_object;

		p1->id=p2->id;
		p1->multiplicity=p2->multiplicity;
		p1->next_object=NULL;
	}

	return *this;
}

std::string Multiset::to_string(char ** object_description) {
	std::string out="";
	Object * p=head->next_object;
	char convert[50];

	if (multiset_size==0) {
		return "#";
	}
	/* For each element in the multiset */
	for (uint i=0; ((i<multiset_size) && (p!=NULL)); i++) {
		if ((p!=NULL)&&(p->id==i)) {
			if (object_description!=NULL)
				out+=object_description[i];
			else {
				sprintf(convert,"%d",i);
				out+=convert;
			}
			if (p->multiplicity!=1) {
				out+="*";
				sprintf(convert,"%d",p->multiplicity);
				out+=convert;
			}
			if (p->next_object!=NULL) out+=", ";
			p=p->next_object;
		}
	}

	return out;
}

unsigned int Multiset::get_length() {
	return length;
}

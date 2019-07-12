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


#include "parserbin.h"

/********************************************/
/* AUXILIAR FUNCTIONS */

/* Auxiliary function:
    Read bytes and returns them into a integer (better than a char buffer for bitmasks */
unsigned int read_bytes (ifstream& is, int size) {
  unsigned int out=0x0;
  char buffer[MAX_FIELD_SIZE];

  is.read(buffer,size);
  // check for I/O errors
  if (is.fail()||is.eof()) {
    is.close();
    if ( (is.rdstate() & ifstream::eofbit ) != 0 )
      cerr << "Error while reading file: In read_bytes, cannot read bytes because reached an unexpected end of file." << endl;
    else if ( (is.rdstate() & ifstream::failbit ) != 0 )
      cerr << "Error while reading file: In read_bytes, cannot read bytes because of an internal error." << endl;
    else if ( (is.rdstate() & ifstream::badbit ) != 0 )
      cerr << "Error while reading file: In read_bytes, cannot read bytes because of an internal I/O error." << endl;

    exit (1);
  }

  for (int i=0; i<size; i++) {
    out = out | ((unsigned char)buffer[i])<<((size-i-1)*8);
  }

  //is >> out;
  return out;
}

/* Auxiliary function:
    Read a text ended on end_text and returns the size of it */
int read_text (ifstream& is, char text_buffer[], int text_buffer_size=256, char end_text='\0') {
  char c='p';
  int i=0;
  for (i=0; (i<text_buffer_size) && (c!=end_text); i++) {
    is>>c;
    // check for I/O errors
    if (is.fail()||is.eof()) {
      is.close();
      if ( (is.rdstate() & ifstream::eofbit ) != 0 )
        cerr << "Error while reading file: In read_text, cannot read char because reached an unexpected end of file." << endl;
      else if ( (is.rdstate() & ifstream::failbit ) != 0 )
        cerr << "Error while reading file: In read_text, cannot read char because of an internal error." << endl;
      else if ( (is.rdstate() & ifstream::badbit ) != 0 )
        cerr << "Error while reading file: In read_text, cannot read char because of an internal I/O error." << endl;

      exit (1);
    }

    if (i==text_buffer_size-1)
      text_buffer[i]=end_text;
    else
      text_buffer[i]=c;
  }

  return i-1;
}

/* Auxiliary function:
    If condition is true, there is an error: close the binary file, print the message and returns true. */
int check_file_error (bool condition, string message, ifstream& is) {
  if (condition) {
      is.close();
      cerr << "Error while reading file: " << message << endl;
      exit (1); // TODO: handle the error in other way
      return 1;
  }
  else {
    return 0;
  }
}


/*********************************/
/* MEMBER FUNCTIONS OF PARSERBIN */

Configuration * Parserbin::readfile (const char * file) {
  int number_objects=0, number_labels=0, number_membranes=0, number_multiset=0, number_rules=0;
  unsigned int recbyte;

  object_id_size=OBJECT_ID_SIZE_LONG;
  label_id_size=MEMBRANE_ID_LABEL_SIZE_LONG;
  mmbr_id_size=MEMBRANE_ID_SIZE_LONG;
  short_object_mode=false;
  short_label_mode=false;
  short_mmbr_mode=false;

  /************************************/
  /* Open the input file as binary */
  ifstream is;
  is.open (file, ios::binary);

  /* If couldn't open it, error*/
  check_file_error(is.fail(), string("Couldn't open binary file ") + file, is);
  /*************************************/

  /************************/
  /* Read the file header */
  recbyte = read_header(is);

  cout << "Using binary file version " << (recbyte & HEADER_VERSION_BITMASK) << endl;

  // TODO: rec & HEADER_VARIANT_BITMASK --> get the variant mode (type of rules)

  cfg=new Configuration();
  /************************/

  /********************************/
  /* Read objects of the alphabet */
  number_objects = read_alphabet(is);
  cout << "Number of objects: " << number_objects << endl;
  /********************************/

  /***********************************/
  /* Read the set of possible labels */
  number_labels = read_label_set(is);
  cout << "Number of labels: " << number_labels << endl;
  /***********************************/

  /*************************************/
  /* Read the initial set of membranes */
  number_membranes = read_membrane_set(is);
  cout << "Number of membranes: " << number_membranes << endl;
  /*************************************/

  /******************************/
  /* Read the initial multisets */
  number_multiset = read_initial_multiset(is,number_membranes);
  cout << "Number of initial multisets: " << number_multiset << endl;
  /******************************/

  /**********************************/
  /* Read the total amount of rules */
  cout << "Using variant: " << (unsigned int)((recbyte & HEADER_VARIANT_BITMASK)>4) << endl;

  number_rules=read_rules(is,recbyte & HEADER_VARIANT_BITMASK, number_labels);
  cout << "Total number of rules: " << number_rules << endl;
  /**********************************/

  cout << "Good binary file, let continue!" << endl;

  is.close();
  //delete cfg;

  return cfg;

}

/*  read_head: Reads the header, and returns the type of P system variant */
unsigned int Parserbin::read_header(ifstream& is) {
  unsigned int recbyte;

  recbyte=read_bytes(is,HEADER_SIZE);

  check_file_error((recbyte & HEADER_BITMASK) != HEADER, "Header is not correct", is);
  return (int)recbyte;
}

/*  read_alphabet: Reads the alphabet, and returns the number of objects */
int Parserbin::read_alphabet(ifstream& is){
  unsigned int recbyte;
  int len=0,number_objects;
  char text_buffer [MAX_TEXT_SIZE];

  recbyte=read_bytes(is,NUMBER_DIFERENT_OBJECTS_SIZE);

  if (recbyte <= 0xFF) {
    short_object_mode=true;
    object_id_size=OBJECT_ID_SIZE_SHORT;
  }

  number_objects=recbyte;

  check_file_error(number_objects<=0, "Number of objects is not correct.", is);

  cfg->alphabet_size=number_objects;
  cfg->alphabet=new char*[cfg->alphabet_size];

  for (int i=0; i<number_objects; i++) {
    len=read_text(is,text_buffer,MAX_TEXT_SIZE,END_OF_TEXT);
    check_file_error(len<=0, "Text size is not correct.", is);

    cfg->alphabet[i]=new char[len+1];
    strcpy(cfg->alphabet[i],text_buffer);
  }

  return number_objects;
}

/*  read_label_set: Reads the label set, and returns the number of labels */
int Parserbin::read_label_set(ifstream& is){
  unsigned int recbyte;
  int len=0,number_labels;
  char text_buffer [MAX_TEXT_SIZE];

  recbyte=read_bytes(is,NUMBER_LABELS_SIZE);

  if (recbyte <= 0xFF) {
    short_label_mode=true;
    label_id_size=MEMBRANE_ID_LABEL_SIZE_SHORT;
  }
  number_labels=recbyte;

  check_file_error(number_labels<=0, "Number of labels is not correct.", is);

  cfg->label_set_size=number_labels;
  cfg->label_set=new char*[cfg->label_set_size];

  for (int i=0; i<number_labels; i++) {
    len=read_text(is,text_buffer,MAX_TEXT_SIZE,END_OF_TEXT);
    check_file_error(len<=0, "Text size is not correct.", is);

    cfg->label_set[i]=new char[len+1];
    strcpy(cfg->label_set[i],text_buffer);
  }
  return number_labels;
}

/*  read_membrane_set: Reads the membrane set, and returns the number of membranes */
int Parserbin::read_membrane_set(ifstream& is){
  unsigned int recbyte;
  int number_membranes;
  unsigned int rmother=0, rlabel=0, rcharge=0;

  recbyte=read_bytes(is,NUMBER_MEMBRANES_SIZE);

  if (recbyte <= 0xFF) {
    short_mmbr_mode=true;
    mmbr_id_size=MEMBRANE_ID_SIZE_SHORT;
  }
  number_membranes=recbyte;

  check_file_error(number_membranes<=0, "Number of membranes is not correct.", is);

  cfg->membrane_set_size=number_membranes;
  cfg->mmbr_set=new Configuration::Membrane_slot[cfg->membrane_set_size];

  /* Initialization of the set*/
  for (int i=0; i<number_membranes; i++) {
    cfg->mmbr_set[i].mother=-1;
    cfg->mmbr_set[i].first_child=-1;
    cfg->mmbr_set[i].next_sister=-1;
    cfg->mmbr_set[i].next_childs_sister=-1;
    cfg->mmbr_set[i].mmbr=NULL;
  }

  /* Information about the skin */
  rmother=read_bytes(is,mmbr_id_size);
  rlabel=read_bytes(is,label_id_size);
  rcharge=read_bytes(is,MEMBRANE_CHARGE_SIZE);
  /* The mother of the skin has not any meaning */
  cfg->mmbr_set[0].mother=rmother;
  cfg->mmbr_set[0].mmbr=new Membrane;
  cfg->mmbr_set[0].mmbr->multiset.set_size(cfg->get_alphabet_size());
  cfg->ini_membrane_id(); /* sets 0 to skin membrane */
  cfg->mmbr_set[0].mmbr->id=cfg->get_new_membrane_id();
  cfg->mmbr_set[0].mmbr->label=rlabel;
  cfg->mmbr_set[0].mmbr->charge=rcharge;

  cfg->skin = cfg->mmbr_set[0].mmbr;

  /* Information about the rest of membranes */
  for (int i=1; i<number_membranes; i++) {
    rmother=read_bytes(is,mmbr_id_size);
    rlabel=read_bytes(is,label_id_size);
    rcharge=read_bytes(is,MEMBRANE_CHARGE_SIZE);
    /* The mother of the current membrane */
    cfg->mmbr_set[i].mother=rmother;
    /* Create new membrane */
    cfg->mmbr_set[i].mmbr=new Membrane;
    cfg->mmbr_set[i].mmbr->id=cfg->get_new_membrane_id();
    cfg->mmbr_set[i].mmbr->label=rlabel;
    cfg->mmbr_set[i].mmbr->charge=rcharge;

    /* If the mother has not any frist child, now it is the current */
    if (cfg->mmbr_set[rmother].first_child<0)
      cfg->mmbr_set[rmother].first_child=i;
    /* Next sister of the current membrane is the one prepared by the mother */
    cfg->mmbr_set[i].next_sister=cfg->mmbr_set[rmother].next_childs_sister;
    /* Update the next sister for childrens from the mother as the current
       (The next sister will take the current as next sister) */
    cfg->mmbr_set[rmother].next_childs_sister=i;
    /* Update next sister of the first children of the mother to the current (round list) */
    cfg->mmbr_set[cfg->mmbr_set[rmother].first_child].next_sister=i;
  }

  return number_membranes;
}

/*  read_multiset: Reads the multisets, and returns the number of them */
int Parserbin::read_initial_multiset(ifstream& is, uint number_membranes){
  unsigned int recbyte,mmbr;
  uint num_multisets;

  /* Read the number of multisets */
  recbyte=read_bytes(is,NUMBER_INITIAL_MULTISETS_SIZE);
  check_file_error((recbyte>number_membranes)||(recbyte==0), string("Number of multisets is not correct."), is);
  num_multisets=recbyte;

  /* Read the membrane id which owns the first multiset */
  mmbr=read_bytes(is,mmbr_id_size);
  check_file_error(mmbr>number_membranes, "Membrane id is not correct.", is);
  /* Fix the pointers of all the membranes and, in the same loop, read the multisets */
  for (uint i=0,j=0;i<number_membranes;i++) {
    /* Set the pointer to the mother */
    if (cfg->mmbr_set[i].mother >= 0 && cfg->mmbr_set[i].mother != i) {
      cfg->mmbr_set[i].mmbr->mother=cfg->mmbr_set[cfg->mmbr_set[i].mother].mmbr; }
    else {
      cfg->mmbr_set[i].mmbr->mother=NULL; }
    /* Set the pointer to the first children */
    if (cfg->mmbr_set[i].first_child > 0) {
      cfg->mmbr_set[i].mmbr->first_child=cfg->mmbr_set[cfg->mmbr_set[i].first_child].mmbr;
    }
    else {
      cfg->mmbr_set[i].mmbr->first_child=NULL;
    }
    /* Set the pointer to the next and prev sister */
    if (cfg->mmbr_set[i].next_sister > 0) {
      cfg->mmbr_set[i].mmbr->next_sister=cfg->mmbr_set[cfg->mmbr_set[i].next_sister].mmbr;
      cfg->mmbr_set[cfg->mmbr_set[i].next_sister].mmbr->prev_sister=cfg->mmbr_set[i].mmbr;
    }
    else {
      cfg->mmbr_set[i].mmbr->next_sister=NULL;
      cfg->mmbr_set[i].mmbr->prev_sister=NULL;
    }

    if (i==mmbr) {
      /* Read membrane multiset */
      read_multiset(is,cfg->mmbr_set[i].mmbr->multiset,cfg->alphabet_size);
      /* Next multiset to read */
      if (++j<num_multisets) {
        mmbr=read_bytes(is,mmbr_id_size);
        check_file_error(mmbr>number_membranes, "Membrane id is not correct.", is);
      }
    }
  }

  cfg->new_environment();
  cfg->get_environment()->set_size(cfg->get_alphabet_size());

  return num_multisets;
}

/*  read_multiset: Reads a multiset, and returns the number of objects */
int Parserbin::read_multiset(ifstream& is, Multiset& multiset, unsigned int alphabet_size) {
  unsigned int recbyte,obj,mult;
  int num_obj;

  multiset.set_size(alphabet_size);

  recbyte=read_bytes(is,NUMBER_DIFERENT_OBJECTS_SIZE);
  num_obj=recbyte;
  check_file_error(num_obj<0,"Number of objects in the multiset is not correct",is);

  for (int i=0;i<num_obj;i++) {
    obj=read_bytes(is,object_id_size);
    mult=read_bytes(is,OBJECT_MULTIPLICITY_SIZE);
    multiset.add_object(obj,mult);
  }
  return num_obj;
}

/**
  read_rules: Reads rules depending on the variant, and returns the number them
*/
int Parserbin::read_rules (ifstream& is, int variant, LabelID number_of_labels) {
	int number_evolution=0, number_sendin=0, number_sendout=0, number_disolution=0, number_division=0;

	cfg->rule_set=new Ruleset(variant,number_of_labels);
	// TODO: This part depends on the P system's variant
	/************************/
	/* Read evolution rules */
	number_evolution = read_evolution_rule(is);
	cout << "Number of evolution rules: " << number_evolution << endl;
	/************************/

	/**********************/
	/* Read send in rules */
	number_sendin = read_sendin_rule(is);
	cout << "Number of sendin rules: " << number_sendin << endl;
	/**********************/

	/***********************/
	/* Read send out rules */
	number_sendout = read_sendout_rule(is);
	cout << "Number of sendout rules: " << number_sendout << endl;
	/***********************/

	/*************************/
	/* Read disolution rules */
	number_disolution = read_disolution_rule(is);
	cout << "Number of disolution rules: " << number_disolution << endl;
	/*************************/

	/***********************/
	/* Read division rules */
	number_division = read_division_rule(is);
	cout << "Number of division rules: " << number_division << endl;
	/***********************/

	return number_evolution+number_sendin+number_sendout+number_disolution+number_division;
}

/*  read_evolution_rule: Reads evolution rules, and returns the number them */
int Parserbin::read_evolution_rule(ifstream& is){
  unsigned int recbyte,lab,chrg,obj;
  int num_rules;
  Multiset *multiset;
  recbyte=read_bytes(is,NUMBER_RULES_SIZE);
  num_rules=recbyte;
  check_file_error(num_rules<0,"Number of evolution rules is not correct",is);

  for (int i=0;i<num_rules;i++) {
    lab=read_bytes(is,label_id_size);
    chrg=read_bytes(is,MEMBRANE_CHARGE_SIZE);
    obj=read_bytes(is,object_id_size);
    multiset = new Multiset;
    read_multiset(is,*multiset,cfg->alphabet_size);
    cfg->rule_set->add_evolution_rule(lab,chrg,obj,*multiset);
    delete multiset;
    multiset=NULL;
  }
  return num_rules;
}

/*  read_send-in_rule: Reads send-in rules, and returns the number them */
int Parserbin::read_sendin_rule(ifstream& is){
  unsigned int recbyte,lab,chrg,nchrg,obj_l,obj_r;
  int num_rules;
  recbyte=read_bytes(is,NUMBER_RULES_SIZE);
  num_rules=recbyte;
  check_file_error(num_rules<0,"Number of sendin rules is not correct",is);

  for (int i=0;i<num_rules;i++) {
    lab=read_bytes(is,label_id_size);
    chrg=read_bytes(is,MEMBRANE_CHARGE_SIZE);
    nchrg=read_bytes(is,MEMBRANE_CHARGE_SIZE);
    obj_l=read_bytes(is,object_id_size);
    obj_r=read_bytes(is,object_id_size);
    cfg->rule_set->add_sendin_rule(lab,chrg,nchrg,obj_l,obj_r);
  }
  return num_rules;
}

/*  read_send-out_rule: Reads send-out rules, and returns the number them */
int Parserbin::read_sendout_rule(ifstream& is){
  unsigned int recbyte,lab,chrg,nchrg,obj_l,obj_r;
  int num_rules;
  recbyte=read_bytes(is,NUMBER_RULES_SIZE);
  num_rules=recbyte;
  check_file_error(num_rules<0,"Number of sendout rules is not correct",is);

  for (int i=0;i<num_rules;i++) {
    lab=read_bytes(is,label_id_size);
    chrg=read_bytes(is,MEMBRANE_CHARGE_SIZE);
    nchrg=read_bytes(is,MEMBRANE_CHARGE_SIZE);
    obj_l=read_bytes(is,object_id_size);
    obj_r=read_bytes(is,object_id_size);
    cfg->rule_set->add_sendout_rule(lab,chrg,nchrg,obj_l,obj_r);
  }
  return num_rules;
}

/*  read_disolution_rule: Reads disolution rules, and returns the number them */
int Parserbin::read_disolution_rule(ifstream& is){
  unsigned int recbyte,lab,chrg,obj_l,obj_r;
  int num_rules;
  recbyte=read_bytes(is,NUMBER_RULES_SIZE);
  num_rules=recbyte;
  check_file_error(num_rules<0,"Number of disolution rules is not correct",is);

  for (int i=0;i<num_rules;i++) {
    lab=read_bytes(is,label_id_size);
    chrg=read_bytes(is,MEMBRANE_CHARGE_SIZE);
    obj_l=read_bytes(is,object_id_size);
    obj_r=read_bytes(is,object_id_size);
    cfg->rule_set->add_disolution_rule(lab,chrg,obj_l,obj_r);
  }
  return num_rules;
}

/*  read_division_rule: Reads division rules, and returns the number them */
int Parserbin::read_division_rule(ifstream& is){
  unsigned int recbyte,lab,chrg,nchrg1,nchrg2,obj_l,obj_r1,obj_r2;
  int num_rules;
  recbyte=read_bytes(is,NUMBER_RULES_SIZE);
  num_rules=recbyte;
  check_file_error(num_rules<0,"Number of division rules is not correct",is);

  for (int i=0;i<num_rules;i++) {
    lab=read_bytes(is,label_id_size);
    chrg=read_bytes(is,MEMBRANE_CHARGE_SIZE);
    nchrg1=read_bytes(is,MEMBRANE_CHARGE_SIZE);
    nchrg2=read_bytes(is,MEMBRANE_CHARGE_SIZE);
    obj_l=read_bytes(is,object_id_size);
    obj_r1=read_bytes(is,object_id_size);
    obj_r2=read_bytes(is,object_id_size);
    cfg->rule_set->add_division_rule(lab,chrg,nchrg1,nchrg2,obj_l,obj_r1,obj_r2);
  }
  return num_rules;
}

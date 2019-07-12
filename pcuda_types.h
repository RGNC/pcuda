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


#if !defined(__PCUDA_TYPES_H__)
#define __PCUDA_TYPES_H__

#include <cutil_inline.h>

/********************/
/* INPUT STRUCTURES */
/********************/

/* Data structure for the input rules */
struct __align__ (16) _rod_rule {
    ushort rinfo;
    uint rev;	 		        /* Evolution rule */
    ushort rsodd[2];			/* Rule that can be: dissolution, send-out, division */
    ushort rsin;			/* Rule for send-in, if exists (it doesn't need a local object) */
};

typedef struct _rod_rule _rodruleset;
typedef struct _rod_rule *Rodruleset;

/* Data structure for the input rules */
struct __align__ (8) _rod_object {
    ushort obj;
    ushort mult;      
};

typedef struct _rod_object _rodmultiset;
typedef struct _rod_object* Rodmultiset;

/* Data structure of the rest membrane information (label and charge) */
struct __align__ (4) _pcuda_membrane {
    short label;
    short charge;
};

typedef struct _pcuda_membrane _membraneset;
typedef struct _pcuda_membrane *Membraneset;


/***************************/
/* INPUT/OUTPUT STRUCTURES */
/***************************/

/* Where the multiset of each membrane will be stored, and where the selected evolution are output */
struct __align__ (8) _multi_selec {
    uint multiplicity;		/* Can be the input multiplicity of the object in the multiset, and the output number of selection times */
    uint rule;			/* The selected evolution rule */
};

typedef struct _multi_selec _multiselec;
typedef struct _multi_selec *Multiselec;

typedef uint _skin_multiset;
typedef uint *Skin_multiset;


/*********************/
/* OUTPUT STRUCTURES */
/*********************/

typedef uint _selec_siodd;
typedef uint *Selec_siodd;


/**************************************/
/* AUXILIARY STRUCTURES IN SHARED MEM */
/**************************************/

struct __align__ (8) _rule_sodd {
	bool rsodd;
	ushort soddtype;
};
typedef struct _rule_sodd _rulesoddset;


#endif	/* !__PCUDA_TYPES_H__ */

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


#include "pcuda.h"
//#include "timestat.h"
//#include "pcuda_sel.cpp"

////////////////////////////////////////////////////////////////////////////////
//// declaration, forward
extern "C"
void pcuda(Configuration * cfg, Pcuda_configuration pcfg, Problem_size ps);

extern "C"
void pcuda_sel(Configuration * cfg, Pcuda_configuration pcfg, Problem_size ps);

/***********************/
/* FUNCTIONS FOR PCUDA */
void pcuda_solve(Configuration *cfg, uint step_limit, uint current_step, uint max_membranes, uint num_objects, uint block_size, int mode, int verbose) {
        struct _pcuda_configuration pcfg;
        struct _problem_size ps;

	ps.step_limit=step_limit;
	ps.max_membranes=max_membranes;
	ps.num_objects=num_objects;
	ps.num_labels=cfg->get_label_set_size();
	ps.num_levels=2;

        /* Initialization of structures */
	pcfg.max_num_membranes=cfg->get_number_membranes()-1;
	pcfg.rodruleset=get_rodruleset(*(cfg->get_rules()),pcfg.rodrule_len,pcfg.rodmultiset_len,&pcfg.rodmultiset, cfg->get_alphabet_size(),cfg->get_label_set_size());
        pcfg.membraneset=get_membraneset(cfg->get_skin(),pcfg.mmbrset_len,pcfg.max_num_membranes,ps.max_membranes,ps.num_levels);
        pcfg.multisets=get_multiset_set(cfg->get_skin(),pcfg.multisets_len,ps.max_membranes,ps.num_levels,ps.num_objects);
        pcfg.skin_multiset=get_skinmultiset(cfg->get_skin(),pcfg.skinmultiset_len,ps.num_objects);
        pcfg.skin_label=cfg->get_skin()->label;
        pcfg.skin_charge=cfg->get_skin()->charge;
        pcfg.environment=new ushort[ps.num_objects];
        cfg->get_environment()->to_array(pcfg.environment,ps.num_objects);
	pcfg.current_step=current_step;
	pcfg.block_size=block_size;
	pcfg.verbose=verbose;

	/* Fast sequential */
	if (mode == 1)
		fast_sequential(cfg,&pcfg, &ps);
	/* Only selection on GPU*/
	else if (mode == 2)
		pcuda_sel(cfg, &pcfg, &ps);
	/* Selection and execution on GPU */
	else
		pcuda(cfg, &pcfg, &ps);
}

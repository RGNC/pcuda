#!/bin/bash

################################################################################
#    Pcuda: Simulating P systems with active membranes on the GPU 
#    This simulator is published on:
#    J.M. Cecilia, J.M. García, G.D. Guerrero, M.A. Martínez-del-Amor, I. Pérez-Hurtado,
#    M.J. Pérez-Jiménez. Simulation of P systems with active membranes on CUDA,
#    Briefings in Bioinformatics, 11, 3 (2010), 313-322
#
#    Pcuda is a subproject of PMCGPU (Parallel simulators for Membrane 
#                                       Computing on the GPU)   
# 
#    Copyright (c) 2009 Miguel Á. Martínez-del-Amor (RGNC, University of Seville)
# 		       Ginés D. Guerrero (GACOP, University of Murcia)
#		       Chema Cecilia (GACOP, University of Murcia)
#		       Ignacio Pérez-Hurtado (RGNC, University of Seville)
#    
#    This file is part of Pcuda.
#  
#    Pcuda is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Pcuda is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Pcuda.  If not, see <http://www.gnu.org/licenses/>.

mkdir -p out
mkdir -p out/test

m=65536
l=17
v=1

for(( objs=4; objs <= 32768; objs = objs*2 ))
do
    if [ $objs -le 256 ]; 
    then
	../../bin/linux/release/pcuda -v $v -i test_data/test/test$objs.bin -l $l -m $m -o $objs -b $objs -t 1 -f > out/test/out_fseq_$objs

        ../../bin/linux/release/pcuda -v $v -i test_data/test/test$objs.bin -l $l -m $m -o $objs -b $objs -t 1 -p 3 > out/test/out_par_$objs
    else
	../../bin/linux/release/pcuda -v $v -i test_data/test/test$objs.bin -l $l -m $m -o $objs -b 256 -t 1 -f > out/test/out_fseq_$objs

        ../../bin/linux/release/pcuda -v $v -i test_data/test/test$objs.bin -l $l -m $m -o $objs -b 256 -t 1 -p 3 > out/test/out_par_$objs
    fi
done


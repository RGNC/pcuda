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

mkdir -p out/test_sat

v=1

echo For n=6

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n6.bin -m 64 -o 920 -b 460 -t 1 -f > out/test_sat/out_fseq_6

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n6.bin -m 64 -o 920 -b 460 -t 1 -p 3 > out/test_sat/out_par_6

echo For n=8

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n8.bin -m 256 -o 1155 -b 385 -t 1 -f > out/test_sat/out_fseq_8

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n8.bin -m 256 -o 1155 -b 385 -t 1 -p 3 > out/test_sat/out_par_8

echo For n=10

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n10.bin -m 1024 -o 729 -b 243 -t 1 -f > out/test_sat/out_fseq_10

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n10.bin -m 1024 -o 729 -b 243 -t 1 -p 3 > out/test_sat/out_par_10

echo For n=12

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n12.bin -m 4096 -o 914 -b 457 -t 1 -f > out/test_sat/out_fseq_12

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n12.bin -m 4096 -o 914 -b 457 -t 1 -p 3 > out/test_sat/out_par_12

echo For n=14

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n14.bin -m 16384 -o 1056 -b 352 -t 1 -f > out/test_sat/out_fseq_14

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n14.bin -m 16384 -o 1056 -b 352 -t 1 -p 3 > out/test_sat/out_par_14

echo For n=16

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n16.bin -m 65536 -o 1131 -b 377 -t 1 -f > out/test_sat/out_fseq_16

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n16.bin -m 65536 -o 1131 -b 377 -t 1 -p 3 > out/test_sat/out_par_16

echo For n=18

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n18.bin -m 262144 -o 2465 -b 493 -t 1 -f > out/test_sat/out_fseq_18

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n18.bin -m 262144 -o 2465 -b 493 -t 1 -p 3 > out/test_sat/out_par_18


echo And now, the sequential version

echo For n=6

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n6.bin > out/test_sat/out_seq_6

echo For n=8

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n8.bin > out/test_sat/out_seq_8

echo For n=10

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n10.bin > out/test_sat/out_seq_10

echo For n=12

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n12.bin > out/test_sat/out_seq_12

echo For n=14

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n14.bin > out/test_sat/out_seq_14

echo For n=16

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n16.bin > out/test_sat/out_seq_16

echo For n=18

../../bin/linux/release/pcuda -v $v -i test_data/sat/sat_n18.bin > out/test_sat/out_seq_18





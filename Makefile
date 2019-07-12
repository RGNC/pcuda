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


################################################################################
#
# Build script for project
#
################################################################################

# Add source files here
EXECUTABLE	:= pcuda
# Cuda source files (compiled with cudacc)
CUFILES_sm_12	:= pcuda.cu selection_only.cu
# C/C++ source files (compiled with gcc / c++)
CCFILES		:= main.cpp multiset.cpp parserbin.cpp ruleset.cpp sequential.cpp seq2pcuda.cpp fast_sequential.cpp pcuda.cpp pcudadef.cpp pcuda_sel.cpp

# Uncomment to debug
#dbg = 1

INCLUDES	+=  -I../../common/counterslib/inc -I../../common/inc -I../../common/include

################################################################################
# Rules and targets

include ../../common/common.mk

LIB		+= -L../../common/lib -L../../common/counterslib/lib -ltimestat

# Uncomment to debug
#CXXFLAGS	+= -g -O0


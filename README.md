# Parallel simulation of P systems with active membranes (PCUDA) #

----------
## 1. Simulation of P systems with active membranes ##

The simulators we have developed are based on a two staged simulation algorithm. This scheme was adopted by the sequential simulator for P systems with active membranes developed in [pLinguaCore](https://www.p-lingua.org). In this design, the simulation process is divided into two stages: *selection stage* and *execution stage*. The selection stage consists of the search for the rules to be executed in each membrane in a given configuration. The rules selected are executed at the execution stage. At the end of the execution stage, the simulation process restarts the selection stage in an iterative way until a halting configuration is reached. This stop condition is 2-fold: a certain number of iterations or a final configuration is reached. 

On one hand, at the beginning of the simulation, we define the maximum number of iterations. On the other hand, a halting configuration is obtained when there are no more rules to select at selection stage. As previously explained, the halting configuration is always reached since it is a simulator for recognizer P systems. The input data for the selection stage consists of a description of the membranes with their multisets (strings over the working alphabet O, labels associated with the membrane in H, etc.) and the set of rules R to be selected. The output data of this stage is the set of selected rules per membrane which will be executed on the execution stage.

The execution stage applies the rules previously selected on the selection stage. During the execution stage, membranes can vary by including new objects, dissolving membranes, dividing membranes, etc, obtaining a new configuration of the simulated P system. This new configuration will be the input data for the selection stage of the next iteration.

----------
## 2. Design of the simulator based on CUDA ##

In the simulators developed in CUDA we have designed a solution based on mapping the double parallelism of P systems over the double parallelism of CUDA. In this sense, a *Block of Threads* (CUDA) is assigned to each *Membrane* of the system, and each *Thread* to each *Object* in the membrane. In this case, each thread corresponds to each object defined in the alphabet.

----------
## 3. Publications ##

### 3.1. Journals ###

* Jose M. Cecilia, José M. García, Ginés D. Guerrero, Miguel A. Martínez-del-Amor, Ignacio Pérez-Hurtado, Mario J. Pérez-Jiménez. **Simulation of P systems with active membranes on CUDA**, *Briefings in Bioinformatics*, 11, 3 (2010), 313-322.
* Jose M. Cecilia, José M. García, Ginés D. Guerrero, Miguel A. Martínez-del-Amor, Ignacio Pérez-Hurtado, Mario J. Pérez-Jiménez. **Implementing P systems parallelism by means of GPUs**, *Lecture Notes in Computer Science*, 5957 (2010), 227-241.

### 3.2. Conference Contributions ###

* Miguel A. Martínez-del-Amor, Jose M. Cecilia, Ginés D. Guerrero, Ignacio Pérez-Hurtado. **An Overview of P System Simulation on GPUs**, *I Jornadas Jóvenes Investigadores*, 17 April 2010, Cáceres, Spain, p.2-7, (2010).
* Ginés D. Guerrero, José M. Cecilia, José M. García, Miguel A. Martínez-del-Amor, Ignacio Pérez-Hurtado, Mario J. Pérez-Jiménez. **Analysis of P systems simulation on CUDA**, *XX Jornadas de Paralelismo*, September 2009, A coruña, Spain, 289-294, (2009). 
* José M. Cecilia, Ginés D. Guerrero, José M. García, Miguel Á. Martínez-del-Amor, Ignacio Pérez-Hurtado, Mario J. Pérez-Jiménez. **Simulation of P Systems with active membranes on CUDA**, *2009 International Workshop on High Performance Computational Systems Biology*, October 2009, Trento, Italy, 61-71, (2009).
* José M. Cecilia, Ginés D. Guerrero, José M. García, Miguel Á. Martínez-del-Amor, Ignacio Pérez-Hurtado, Mario J. Pérez-Jiménez. **A massively parallel framework using P systems and GPUs**, *Symposium on Application Accelerators in High Performance Computing*, July 2009, Illinois, USA, (2009)
* Miguel A. Martínez-del-Amor, Ignacio Pérez-Hurtado, Mario J. Pérez-Jiménez, Jose M. Cecilia, Ginés D. Guerrero, José M. García. **Simulating active membrane systems using GPUs**, *10th Workshop on Membrane Computing*, August 2009, Curtea de Arges, Rumania, 369-384, (2009).
* Miguel A. Martínez-del-Amor, Ignacio Pérez-Hurtado, Mario J. Pérez-Jiménez, José M. Cecilia, Ginés D. Guerrero, José M. García. **Simulation of recognizer P systems by using manycore GPUs**, *7th Brainstorming Week on Membrane Computing*, February 2009, Volume II, Seville, Spain, 45-58, (2009).

### 3.3 Ph.D. Thesis ###

* Miguel Á. Martínez-del-Amor. [Accelerating Membrane Systems Simulators using High Performance Computing with GPU.](http://www.cs.us.es/~mdelamor/research.html#thesis) May 2013, University of Seville. Advised by Mario J. Pérez-Jiménez and Ignacio Pérez-Hurtado.
* José M. Cecilia. The GPU as a Processor for Novel Computation: Analysis and Contributions. 2011, University of Murcia. Advised by J.M. García and M. Ujaldón.

----------
## 4. Downloads ##

[Link to PCUDA Files](http://sourceforge.net/projects/pmcgpu/files/PCUDA/)

[Required Counterslib library](http://sourceforge.net/projects/pmcgpu/files/counterslib)

Read the howto.pdf (extract from Miguel A. Martínez-del-Amor's thesis) for futher information about the simulators and a getting started guide. You can find it in the [root folder of files of PMCGPU](http://sourceforge.net/projects/pmcgpu/files).

**(Warning)**
*Results showed in our publications were obtained not using the fast simulator as it is provided here. This simulator was developed and improved after that time. Please, read Martínez-del-Amor's [doctoral dissertation](http://www.cs.us.es/~mdelamor/research.html#thesis) for further information, and new results.*

----------
## 5. How to acknowledge ##

If you intend to create a branch of PCUDA, or use its produced results, please consider citing the following publication:

*Jose M. Cecilia, José M. García, Ginés D. Guerrero, Miguel A. Martínez-del-Amor, Ignacio Pérez-Hurtado, Mario J. Pérez-Jiménez. Simulation of P systems with active membranes on CUDA, Briefings in Bioinformatics, 11, 3 (2010), 313-322.*

----------
## 6. Funding ##

TIN2006–13425 of the Ministerio de Educación y Ciencia of Spain, cofinanced by FEDER funds, and ‘‘Proyecto de Excelencia con Investigador de Reconocida Valía’’ of the Junta de Andalucía under grant P08-TIC04200 to M.A.M., I.P.-H. and M.J.P.-J.

Fundación Séneca (Agencia Regional de Ciencia y Tecnología, Región de Murcia) under grant 00001/CS/2007, and also by the Spanish MEC and European Commission FEDER under grant CSD2006-00046 to J.M.C., J.M.G., G.D.G.

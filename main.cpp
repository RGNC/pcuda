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
 * main.cpp
 *
 *	The main routine of the simulator
 *  Created on: 28-dic-2008
 *      Author: miguel
 */

#include "parserbin.h"
#include "sequential.h"

/* An auxiliary function for creating a test binary file */
void build_binary_file();
void read_config_file(const char * file);

int main (int argc, char* argv[]) {

    Parserbin p;
    Configuration * cfg = NULL;
    int verbose_mode=0, step_limit=256,mode=0;
    char c;
    bool config=false,seq_mode=false;
    string input_file="",config_file="";
    int numparams=0,max_membranes=0,num_objects=0,block_size=0,threshold=0;

    while ((c = getopt (argc, argv, "v:c:l:m:o:b:t:i:hsfp:")) != -1) {
        switch (c) {
        case 'v':
            verbose_mode = atoi(optarg);
            break;
        case 'c':
	    config=true;
            config_file=optarg;
            read_config_file(config_file.c_str());
            break;
        case 'l':
	    numparams++;
            step_limit=atoi(optarg);
            break;
        case 'i':
	    numparams++;
            input_file = optarg;
            break;
        case 'm':
            numparams++;
            max_membranes= atoi(optarg);
            break;
        case 'o':
            numparams++;
            num_objects = atoi(optarg);
            break;
        case 'b':
            numparams++;
            block_size = atoi(optarg);
            break;
        case 't':
            numparams++;
            threshold = atoi(optarg);
            break;
	case 's':
	    seq_mode=true;
	    break;
	case 'f':
	    mode=1;
	    break;
	case 'p':
	    mode=atoi(optarg);
	    break;
        default:
        case 'h':
        case '?':
	    cout << "Copyright (C) 2009 M.A. Martínez-del-Amor, G.D. Guerrero, J.M. Cecilia, I. Pérez-Hurtado" << endl <<
	            "This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome to redistribute it under certain conditions" << endl << endl;
            cout<< "Usage: pcuda <params>. Where <params> can be:" << endl;
            cout<< "-vX: Indicates the verbosity level. No verbose activated by default."<< endl << "The verbosity levels are: " << endl;
            cout <<"* -v1: Print only the last configuration." << endl;
            cout <<"* -v2: Print only the configuration of the skin in every step, and at the end, the configuration of the remaining membranes"<<endl;
            cout <<"* -v3: Print all the information of all the membranes in each configuration, the alphabet and the set of rules"<<endl;
            cout <<"-c: Define a configuration file (no implemented yet). Default: config.cfg"<< endl;
            cout <<"-l: Define the maximum number of steps to do. Default: 256"<< endl;
            cout <<"-i: Define the input binary file" << endl;
            cout <<"-s: Executes only the sequential algorithm. Activated by default until defining some of the next params" << endl;
	    cout <<"-f: Launch fast sequential simulation (also -p 1)" << endl;
	    cout <<"-p: Set the algorithm to be used in fast mode: 1=Fast sequential (also -f), 2=Only selection in parallel, 3(default)=Selection and execution in parallel" << endl << endl;
	    cout <<"The next params are mandatory when configuring the fast algorithms:" << endl;
	    cout <<"-m: Define the maximum number of membranes that the P system will create" << endl;
	    cout <<"-o: Define the maximum number of objects that one membrane can have in one step" << endl;
	    cout <<"-b: Define the number of threads to execute per block" << endl;
	    cout <<"-t: Define the threshold to achieve for executing the parallel algorithm (number of membranes)"<<endl << endl;
            cout << "Pcuda, a simulator for recognizer P systems with active membranes." << endl;
            cout << "Version 1.0." << endl;

	    return 0;
        }

    }

    if (!seq_mode && numparams<5) {
	    cout << "Error with the input params. Please, use the help typing -h" << endl;
	    return -1;
    }
    else if (seq_mode && threshold==0) {
	    threshold=INT_MAX;
    }

    /*
     * Read and parse an input binary file. Exits if there is
     * any syntax or semantic mistake inside the binary file */
    cfg = p.readfile(input_file.c_str());

    if (cfg == NULL) {
        cerr << "Error while reading binary file " << input_file << endl;
        return 1;
    }

    if (verbose_mode>1) {
    Rulelist* r=NULL;
    int rl[2]={0,0};
    int rc[3]={0,0,0};
    int rt=0;
    for (int i=0; i< cfg->get_label_set_size(); i++)
	    for (int j=0; j<3; j++) {
		    r=cfg->get_rules()->get_rulelist(i,j);
		    int k=0;
		    for (r->start_iteration();!r->end();r->next_rule())
			    k++;
		    cout << "Rules for (label,charge): (" << i << "," << j << "): are" << k << endl; 
		    rl[i]+=k;
		    rc[j]+=k;
		    rt+=k;
	    }
		 
    cout << "Number of rules per label: 1=" << rl[0] << ", 2=" << rl[1]<<endl;
    cout << "Number of rules per charge: 0=" << rc[0] << ", +=" << rc[1] << ", -=" << rc[2] << endl;
    cout << "Total number: " << rt << endl;
    }    
    
    /* Call to the sequential solution */
    sequential_solve(cfg,verbose_mode,step_limit,max_membranes,num_objects,block_size,threshold,mode);

    delete cfg;
}



/******************/
void read_config_file(const char * file) {
    char buffer[512], c='\0';
    int i=0;

    FILE * f=NULL;
    f=fopen(file,"ro");
    if (f==NULL)
        perror("Couldn't open config file");

    while (!feof(f)&&!ferror(f)) {
        fscanf(f,"%s",(char *) &buffer);

        if (buffer[0]=='#') { /* Discard comments */
            do {
                c=fgetc(f);
            } while (c!='\n' && c!=EOF);
        } else {				/* Read the variable value */
            fscanf(f,"%d", &i);

            cout<<buffer<<" "<<i<<endl;
        }
    }

    fclose(f);

}


void build_binary_file() {


    ofstream outfile ("small.bin",ofstream::binary);
    char buffer [4];

// Header
    buffer[0]=0xAF;
    buffer[1]=0x12;
    buffer[2]=0xFA;
    buffer[3]=0x11;

    outfile.write(buffer,4);

// Number of objects
    buffer[0]=0x00;
    buffer[1]=0x02;
    outfile.write(buffer,2);
// Objects ids
    buffer[0]='o';
    buffer[1]='b';
    buffer[2]='1';
    buffer[3]='\0';
    outfile.write(buffer,4);
    buffer[2]='2';
    outfile.write(buffer,4);

// Number of labels
    buffer[0]=0x00;
    buffer[1]=0x02;
    outfile.write(buffer,2);
// Labels ids
    buffer[0]='m';
    buffer[1]='e';
    buffer[2]='1';
    buffer[3]='\0';
    outfile.write(buffer,4);
    buffer[2]='2';
    outfile.write(buffer,4);

// Number of membranes
    buffer[0]=0x00;
    buffer[1]=0x04;
    outfile.write(buffer,2);
// Membrane 0
    outfile.write(buffer,1); // Id father
    outfile.write(buffer,1); // Label id
    outfile.write(buffer,1); // Charge
// Membrane 1
    buffer[0]=0x00;
    outfile.write(buffer,1); // Id father
    buffer[0]=0x01;
    outfile.write(buffer,1); // Label id
    buffer[0]=0x00;
    outfile.write(buffer,1); // Charge
// Membrane 2
    buffer[0]=0x00;
    outfile.write(buffer,1); // Id father
    buffer[0]=0x01;
    outfile.write(buffer,1); // Label id
    buffer[0]=0x00;
    outfile.write(buffer,1); // Charge
// Membrane 3
    buffer[0]=0x01;
    outfile.write(buffer,1); // Id father
    buffer[0]=0x00;
    outfile.write(buffer,1); // Label id
    buffer[0]=0x01;
    outfile.write(buffer,1); // Charge

// Number of multisets
    buffer[0]=0x00;
    buffer[1]=0x03;
    outfile.write(buffer,2);
// Membrane 0
    buffer[0]=0x00;
    outfile.write(buffer,1); // Membrane ID
    buffer[0]=0x00;
    buffer[1]=0x01;
    outfile.write(buffer,2); // Number of objects
    buffer[0]=0x00;
    outfile.write(buffer,1); // Object ID
    buffer[0]=0x00;
    buffer[1]=0x01;
    outfile.write(buffer,2); // Multiplicity
// Membrane 1
    buffer[0]=0x01;
    outfile.write(buffer,1); // Membrane ID
    buffer[0]=0x00;
    buffer[1]=0x01;
    outfile.write(buffer,2); // Number of objects
    buffer[0]=0x01;
    outfile.write(buffer,1); // Object ID
    buffer[0]=0x00;
    buffer[1]=0x02;
    outfile.write(buffer,2); // Multiplicity
// Membrane 2
    buffer[0]=0x02;
    outfile.write(buffer,1); // Membrane ID
    buffer[0]=0x00;
    buffer[1]=0x02;
    outfile.write(buffer,2); // Number of objects
    buffer[0]=0x00;
    outfile.write(buffer,1); // Object ID
    buffer[0]=0x00;
    buffer[1]=0x04;
    outfile.write(buffer,2); // Multiplicity
    buffer[0]=0x01;
    outfile.write(buffer,1); // Object ID
    buffer[0]=0x00;
    buffer[1]=0x05;
    outfile.write(buffer,2); // Multiplicity

// Number of evolution rules
    buffer[0]=0x00;
    buffer[1]=0x02;
    outfile.write(buffer,2);
// Rule 0
    buffer[0]=0x00;
    outfile.write(buffer,1); // Label ID
    buffer[0]=0x00;
    outfile.write(buffer,1); // Charge
    buffer[0]=0x00;
    outfile.write(buffer,1); // Object ID left
    buffer[0]=0x00;
    buffer[1]=0x01;
    outfile.write(buffer,2); // Number of objects
    buffer[0]=0x01;
    outfile.write(buffer,1); // Object ID
    buffer[0]=0x00;
    buffer[1]=0x05;
    outfile.write(buffer,2); // Multiplicity
// Rule 1
    buffer[0]=0x01;
    outfile.write(buffer,1); // Label ID
    buffer[0]=0x00;
    outfile.write(buffer,1); // Charge
    buffer[0]=0x01;
    outfile.write(buffer,1); // Object ID left
    buffer[0]=0x00;
    buffer[1]=0x02;
    outfile.write(buffer,2); // Number of objects
    buffer[0]=0x00;
    outfile.write(buffer,1); // Object ID
    buffer[0]=0x00;
    buffer[1]=0x01;
    outfile.write(buffer,2); // Multiplicity
    buffer[0]=0x01;
    outfile.write(buffer,1); // Object ID
    buffer[0]=0x00;
    buffer[1]=0x01;
    outfile.write(buffer,2); // Multiplicity

// Number of sendin rules
    buffer[0]=0x00;
    buffer[1]=0x01;
    outfile.write(buffer,2);
// Rule 0
    buffer[0]=0x00;
    outfile.write(buffer,1); // Label ID
    buffer[0]=0x01;
    outfile.write(buffer,1); // Charge
    buffer[0]=0x02;
    outfile.write(buffer,1); // New charge
    buffer[0]=0x00;
    outfile.write(buffer,1); // Object ID left
    buffer[0]=0x01;
    outfile.write(buffer,1); // Object ID right

// Number of sendout rules
    buffer[0]=0x00;
    buffer[1]=0x00;
    outfile.write(buffer,2);

// Number of disolution rules
    buffer[0]=0x00;
    buffer[1]=0x00;
    outfile.write(buffer,2);

// Number of division rules
    buffer[0]=0x00;
    buffer[1]=0x00;
    outfile.write(buffer,2);

    outfile.close();

}

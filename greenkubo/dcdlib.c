/*
 * 
 * 
 * 
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

float PDBVELFACTOR=20.45482706;

//int DIM=3;

/*
    read_dcd_header reads the header of a dcd file returning the number of atoms and number of steps
*/
void read_dcd_header(FILE *in, int *nAtoms, int *nSteps) {
	int magicnumber;
	char hdr[5];
	int junk;
	int i;
	int ntitle;
	char title[80];
	float junkf;
	int deltaStep;
	float deltaTime;
	hdr[4]='\0';
	// magic number
	fread(&magicnumber,sizeof(int),1,in);
	// header
	fread(hdr,4*sizeof(char),1,in);
	//read a junk int
	// number of steps
	fread(nSteps,sizeof(int),1,in);
	printf("Number of steps in dcd file: %d\n",*nSteps);
	fread(&junk,sizeof(int),1,in);
	// number of steps between writes
	fread(&deltaStep,sizeof(int),1,in);
	// a bunch of crap integers I dont care about
	for(i=0;i<6;i++) {
		fread(&junk,sizeof(int),1,in);
	}
	// read simulation time between steps
	fread(&deltaTime,sizeof(float),1,in);
	// a bunch of crap integers I dont care about
	for(i=0;i<12;i++) {
		fread(&junk,sizeof(int),1,in);
	}
	// the number of title strings
	fread(&ntitle,sizeof(int),1,in);
	// the title
	for(i=0;i<ntitle;i++) {
		fread(title,80*sizeof(char),1,in);
	}
	// a bunch of crap integers I dont care about
	for(i=0;i<2;i++) {
		fread(&junk,sizeof(int),1,in);
	}
	// number of atoms
	fread(nAtoms,sizeof(int),1,in);
	printf("Number of atoms in dcd file: %d\n",*nAtoms);
	// magic number
	fread(&magicnumber,sizeof(int),1,in);
	//  the end of the header

}



/* 
 * a subroutine to read the coordinates from a dcd file
 * this subroutine takes the file pointer, number of stepsm, and number of atoms as input.
 * the filled out axis and coordinate arrays are returned.
 */

void read_dcd_coord(FILE *d, int nAtoms, int nSteps, double **axis, float ***coord) {
	int   atom;
	int   k;
	int   junk;
	float temp;
	int step;
	long int filePosition;

	filePosition = ftell(d);
	for (step=0;step<nSteps;step++) {

		if(step%100==0) printf("Reading step %d from trajectory file\n",step);
/*
		// skip the first input of this data section
		fseek(d,filePosition+1*sizeof(int),0);
		// read box size informtaion
		fread(&axis[step][0],sizeof(double),1,d);
		fseek(d,ftell(d)+1*sizeof(double),0); // skip one double
		fread(&axis[step][1],sizeof(double),1,d); 
		fseek(d,ftell(d)+2*sizeof(double),0); // skip two doubles
		fread(&axis[step][2],sizeof(double),1,d);
       	 	// skip the last line of this data section
		fseek(d,ftell(d)+sizeof(int),0);
*/
        	// read coordinates
		for(k=0;k<3;k++) {
        	   	fread(&junk,sizeof(int),1,d);
           		for (atom=0;atom<nAtoms;atom++) {
              			fread(&temp,sizeof(float),1,d);
	      			coord[step][atom][k] = temp*PDBVELFACTOR;
	   		}
           		fread(&junk,sizeof(int),1,d);
		}	

		filePosition = ftell(d);
	}

}



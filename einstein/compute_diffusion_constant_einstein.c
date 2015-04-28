// Program to compute the diffussion constant of water from MD simulation using the einstein relationship
// USAGE :: ./compute_diffusion_constant_einstein.x -dcd 
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dcdlib.h"

// declare constants
int MAX_STRING_SIZE=1024;
int DIM = 3;
int buffer = 100;

// declare subroutines
void parse_command_line(int, char **, char*, char*, float *);
void compute_pos2_autocorr(float ***, double **, int, int, float, FILE *);
double dist2_box(float *, float *, double *);

int main(int argc, char *argv[]) {


	char dcdFileName[MAX_STRING_SIZE];
	FILE *dcdFile;

	char outputFileName[MAX_STRING_SIZE];
	FILE *outFile;

	int nAtoms;
	int nSteps;
	float ***coord;
	double **box;
	float deltaTime;

	int junk;
	int i, j, k;

	// read the name of the dcd file from the command line
	parse_command_line(argc, argv, dcdFileName, outputFileName,&deltaTime);

	// open dcd file
	dcdFile = fopen(dcdFileName,"r");

	// read dcd header info
	read_dcd_header(dcdFile,&nAtoms,&nSteps);

	// allocate coordinate and box array
	coord = (float***) malloc(nSteps*sizeof(float**));
	box = (double**) malloc(nSteps*sizeof(double*));
	for (i=0;i<nSteps;i++) {
		coord[i] = (float**) malloc(nAtoms*sizeof(float*));
		box[i] = (double*) malloc(DIM*sizeof(double));
		for (j=0;j<nAtoms;j++) {
			coord[i][j] = (float*) malloc(DIM*sizeof(float));
		}
	}

	// read all coordinates from DCD file
	read_dcd_coord(dcdFile, nAtoms, nSteps, box, coord);
	printf("Finished reading coordinates\n");
	// close dcd file
	fclose(dcdFile);

	// open output data file
	outFile = fopen(outputFileName,"w");
	// compute square position autocorrelation function
	compute_pos2_autocorr(coord,box,nAtoms,nSteps,deltaTime,outFile);
	// close output file
	fclose(outFile);

}
/*
 *
 *				SUBROUTINES
 *
 */

// parse command line (look for options and place values into variables)
void parse_command_line(int argc, char *argv[], char *dcdFileName, char *outputFileName, float *deltaTime) {

	int i;

	// initialize strings
	memset(dcdFileName, '\0', sizeof(dcdFileName));
	memset(outputFileName, '\0', sizeof(outputFileName));

	// Check to see that there is more than one command line argument
	if (argc==1) {
		printf("ERROR: too few command line arguments.\n");
		printf("USAGE: compute_diffusion_constant_einstein.x -dcd [dcd file name]\n");
		exit(1);
	}

	// Loop through command line arguments and parse appropriate files
	for (i=1; i<argc;i++) {
		if (strcmp(argv[i],"-dcd") ==0) {
			strcpy(dcdFileName,argv[++i]);
			printf("dcd file:%s\n",dcdFileName);
		} else if (strcmp(argv[i],"-out") ==0) {
			strcpy(outputFileName,argv[++i]);
			printf("output file:%s\n",outputFileName);
		} else if (strcmp(argv[i],"-dt") ==0) {
			*deltaTime = atof(argv[++i]);
			printf("delta time:%10.5f\n",*deltaTime);
		} else {
			printf("ERROR: unrecognized command line argument: %s\n:", argv[i]);
			printf("USAGE: compute_diffusion_constant_einstein.x -dcd [dcd file name]\n");
			exit(2);
		}
	}
}

// function to compute square distance autocorrelation
void compute_pos2_autocorr(float ***coord, double **box, int nAtoms, int nSteps, float deltaTime, FILE *outFile) {

	int atom;
	int step1;
	int step2;
	int i;
	int deltaStep;
	int maxDeltaStep;
	int count;
	double dist2AutoCorr;
	double temp;

	maxDeltaStep = nSteps-1-buffer;

	for (deltaStep=1;deltaStep<maxDeltaStep;deltaStep++) {
		if (deltaStep%100==0) {
			printf("Working on correlation for delta T=%10d\n", deltaStep);
		}
		count=0;
		dist2AutoCorr=0;
		for(step1=0;step1<(nSteps-deltaStep);step1++) {
			step2 = step1+deltaStep;
			for (atom=0;atom<nAtoms;atom+=3) { // skip hydrogens
//				dist2AutoCorr += dist2_box(coord[step1][atom],coord[step2][atom],box[step1]);
				temp = dist2_box(coord[step1][atom],coord[step2][atom],box[step1]);
				if (temp<0) {
					printf("dist2 < 0 - step1: %5d, step2: %5d, atom: %5d\n",step1,step2,atom);
				}
				dist2AutoCorr += temp;
				count++;
			}
		}
		dist2AutoCorr /= (double) (count);
		fprintf(outFile,"%10d %30.15f\n", deltaStep,dist2AutoCorr);
		fflush(outFile);
	}


}

// function to compute square distance between two vectors in R3 given box size info
double dist2_box (float *coord1, float *coord2, double *box) {

	double temp;
	int i;
	double dist2;

	dist2=0;
	for (i=0;i<3;i++) {
		temp = (double) (coord1[i]-coord2[i]);
		// check periodic boundaries
		if (temp< -box[i]/2.0) {
			temp += box[i];
		} else if (temp > box[i]/2.0) {
			temp -= box[i];
		}
		dist2 += temp*temp;
	}

	return dist2;

}



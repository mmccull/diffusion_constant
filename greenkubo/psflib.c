#include "psflib.h"

void read_psf_header( FILE *psfFile, int *nAtoms)
{
	char buffer[1024];
	char key[8];
	char temp[10];
	
	
	
	/*
		Read the number of atoms from the header of the psf file
	*/
	while( fgets( buffer, 1024, psfFile ) != NULL )
	{
		strncpy(key,buffer+9,6);
//		key[6]='\0';
		if( strncmp(key,"!NATOM",6)==0 )
		{
			strncpy(temp,buffer,9);
			temp[9]='\0';
			*nAtoms = atoi( temp );
			break;
		}
	}
}

char** read_psf_atom_data(FILE  *psfFile, int nAtoms, int *nonHatoms, int *numNonHatoms, int *nAtomTypes, int *uniqueAtomTypeNum) {

	char buffer[1024];
	int i, j;
	char temp[6];
	int flag;
	int type;
	char atomType[6];
	char atomName[6];
	char uniqueAtomType[nAtoms][5];
	char **smallUniqueType;
	int uniqueTypeCount;


	uniqueTypeCount = 0;
	*numNonHatoms = 0;
	for (i=0;i<nAtoms;i++) {
		fgets(buffer,1024,psfFile);
		strncpy(temp,buffer+24,4);
		temp[4]='\0';
		strncpy(atomType,buffer+29,4);
		atomType[4]='\0';
		if (temp[0]=='H') {
			nonHatoms[i] = -99;
			uniqueAtomTypeNum[i] = -99;
		} else {
			if (*numNonHatoms==0) {
				strncpy(uniqueAtomType[uniqueTypeCount],atomType,5);
				uniqueAtomTypeNum[i] = uniqueTypeCount;
				uniqueTypeCount++;
			} else {
				flag = 0;
				for (type=0;type<uniqueTypeCount;type++) {
					if (strncmp(uniqueAtomType[type],atomType,4)==0) {
						flag = 1;
						uniqueAtomTypeNum[i] = type;
						break;
					} 
				}
				if (flag == 0) {
					strncpy(uniqueAtomType[uniqueTypeCount],atomType,5);
					uniqueAtomTypeNum[i] = uniqueTypeCount;
					uniqueTypeCount++;
				}
			}
			nonHatoms[i] = *numNonHatoms;
			*numNonHatoms = *numNonHatoms + 1;
		}

	}

	printf ("Number of non-hydrogens = %d\n",*numNonHatoms);

	*nAtomTypes = uniqueTypeCount;

	smallUniqueType = new char* [*nAtomTypes];

	for (type=0;type<*nAtomTypes;type++) {
		smallUniqueType[type] = new char [5];
		strncpy(smallUniqueType[type],uniqueAtomType[type],5);
	}


	return smallUniqueType;

	delete [] smallUniqueType;

}



void read_psf_bond_data( FILE *psfFile, int nAtoms, int numNonHatoms, int *nonHatoms, double **sasaPij)
{
	double pij_bonded=0.8875;
	double pij_nonbonded=0.3516;
	int bondsPerLine = 4;
	char buffer[1024];
	char key[8];
	char temp[10];
	int nBonds;
	int nLines;
	int line;
	int bond;
	int bondLineCount;
	int atom1;
	int atom2;
	int stringPos;
	
	
	// assume nonbonded	
	for (atom1=0;atom1<numNonHatoms;atom1++) {
		for (atom2=0;atom2<numNonHatoms;atom2++) {
			sasaPij[atom1][atom2] = pij_nonbonded;
		}
	}

	/*
		Read the number of atoms from the header of the psf file
	*/
	while( fgets( buffer, 1024, psfFile ) != NULL )
	{
		strncpy(key,buffer+9,6);
		if( strncmp(key,"!NBOND",6)==0 )
		{
			strncpy(temp,buffer,9);
			temp[9]='\0';
			nBonds = atoi( temp );
			printf("number of bonds: %d\n",nBonds);
			break;
		}
	}
	nLines = (int) ( ((double) nBonds) / ((double) bondsPerLine) +0.75);

	bond = 0;
	for (line=0;line<nLines;line++) {
		
		fgets( buffer, 1024, psfFile );
		stringPos = 0;
		for(bondLineCount=0;bondLineCount<bondsPerLine;bondLineCount++) {
			if (bond>=nBonds) break;
			bond++;
			strncpy(temp,buffer+stringPos,8);
			temp[8] = '\0';
			atom1 = nonHatoms[atoi(temp)-1];
			strncpy(temp,buffer+stringPos+8,8);
			atom2 = nonHatoms[atoi(temp)-1];
			if (atom1 !=-99 && atom2 !=-99) {
				sasaPij[atom1][atom2] = pij_bonded;
				sasaPij[atom2][atom1] = pij_bonded;
			}
			stringPos+=16;
		}	
		

	}

}


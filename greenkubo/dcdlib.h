
/*
 *  A set of subroutines used for read, writing and manipulating coordinates from a DCD file
 * 
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void read_dcd_header(FILE *, int *, int *);

void read_dcd_coord(FILE *, int, int, double **, float ***);


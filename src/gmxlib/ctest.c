/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Great Red Oystrich Makes All Chemists Sane
 */
#include <limits.h>
#include <math.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>

#include "xdrf.h"

/*________________________________________________________________________
 |
 | read_frame - ad hoc implementation to read ascii data
 |
 | given a open file handle and the number of coordinates, this routine
 | reads sets of three floating point numbers and stores the result in
 | coord. the coord array should be preallocated to hold 3 * num_of_coord
 | coordinates. Minval and maxval will return the minumum and maximum
 | values of the coordinates.
 | Return value will be -1 when eof is reached, and 0 on success
*/

int read_frame(FILE *in, int num_of_coord, float coord[]) {
    int i;
    double d[3];
    char line[1024];
    for (i = 0; i < num_of_coord; i++) {
	if (fgets(line, 1024 -1, in) == NULL) {
	    if (i == 0) return -1;
	    fprintf(stderr,"incomplete coordinate frame during read %d", i);
	    exit (1);
	}
	sscanf(line," %lf %lf %lf", &d[0], &d[1], &d[2]);
	*coord++ = (float) d[0];
	*coord++ = (float) d[1];
	*coord++ = (float) d[2];
    }
    return 0;
}



int main() {
    XDR xd, xd2;
    int i, j;
    float prec = 80000000.0;
    float *coord, *coord2;
    int num_of_coord;
    char *line;
    int framecnt = 0;
    double maxdiff = 0;
    float d0, d1, d2, d3;
    FILE *fgmx, *fout;
    
    line = (char *) malloc(1024);

    /* open file containing 3d coordinates */
    fgmx = fopen("test.gmx","r");
    if (fgmx == NULL) {
	perror("could open gmx test data");
	exit(1);
    }
    if (fgets(line, 1024 -1, fgmx) == NULL) {
	perror("cannot read gmx test data");
	exit (1);
    }
    
    /* read first line which contains number of coordinates */
    sscanf(line," %d %f %f %f %f", &num_of_coord, &d0, &d1, &d2, &d3);
    
    coord = (float *)malloc(num_of_coord * 3 * sizeof(float));
    coord2 = (float *)malloc(num_of_coord * 3 * sizeof(float));

    /* open xdr file to which compressed coordinates will be written */
    if (xdropen(&xd, "test.xdr","w") == 0) {
	fprintf(stderr,"failed to open file\n");
    }
    
    /* just as test write the first line using normal xdr routine */
    xdr_int(&xd, &num_of_coord);
    xdr_float(&xd, &d0);
    xdr_float(&xd, &d1);
    xdr_float(&xd, &d2);
    xdr_float(&xd, &d3);
    /* read all frames from the data and write compressed coordinates */
    while ( read_frame(fgmx, num_of_coord, coord) == 0 ) {
	framecnt++;
	if (xdr3dfcoord(&xd, coord, &num_of_coord, &prec) == 0) {
	    fprintf(stderr,"error while writing coordinates\n");
	}
    }
    xdrclose(&xd);
    fclose(fgmx);
    
    
    
    /* Now do the inverse ! */
    
    /* open file to write decompressed data */
    fout = fopen("test.out", "w+");
    if (fout == NULL) {
	perror("could not open test.out to write data\n");
	exit(1);
    }
    if (xdropen(&xd2, "test.xdr","r") == 0) {
	fprintf(stderr,"error while opening test.xdr for read\n");
	exit(1);
    }
    *line = '\0';
    xdr_int(&xd2, &num_of_coord);
    xdr_float(&xd2, &d0);
    xdr_float(&xd2, &d1);
    xdr_float(&xd2, &d2);
    xdr_float(&xd2, &d3);
    fprintf(fout, "%5d%8.3f%8.3f%8.3f%8.3f\n", num_of_coord, d0, d1, d2, d3);
    for (i = 0; i < framecnt ; i++) {
	if (xdr3dfcoord(&xd2, (float *)coord2, &num_of_coord, &prec) == 0) {
	    fprintf(stderr, "error while reading coordinates\n");
	}
	for (j=0; j < num_of_coord * 3; j += 3) {
	    fprintf(fout, "%8.3f %8.3f %8.3f\n", coord2[j], 
		coord2[j+1] ,coord2[j+2]);
	}
    }
    xdrclose(&xd2);
    fclose(fout);
    
    
    
    /* And compare the result */
    
    fgmx = fopen("test.gmx", "r");
    fout = fopen("test.out", "r");
    maxdiff = 0;
    fgets(line, 1024 -1, fgmx);
    fgets(line, 1024 -1, fout);
    for (i = 0; i < framecnt ; i++) {
        read_frame(fgmx, num_of_coord, coord);
        read_frame(fout, num_of_coord, coord2);
	for (j=0; j < num_of_coord * 3; j++) {
	    if (fabs(coord[j] - coord2[j]) > maxdiff) 
	    maxdiff = fabs(coord[j] - coord2[j]) ;
	}
    }
    fprintf(stderr,"\nmaxdiff  = %f\n", maxdiff);    
    return 0;
}

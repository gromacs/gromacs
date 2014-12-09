/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */


#ifndef LZ77_H
#define LZ77_H

void Ptngc_comp_to_lz77(unsigned int *vals, const int nvals,
		  unsigned int *data, int *ndata,
		  unsigned int *len, int *nlens,
		  unsigned int *offsets, int *noffsets);

void Ptngc_comp_from_lz77(unsigned int *data, const int ndata,
		    unsigned int *len, const int nlens,
		    unsigned int *offsets, const int noffsets,
		    unsigned int *vals, const int nvals);

#endif

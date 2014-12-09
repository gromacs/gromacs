/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */


#ifndef RLE_H
#define RLE_H

void Ptngc_comp_conv_to_rle(unsigned int *vals, const int nvals,
		                    unsigned int *rle, int *nrle,
		                    const int min_rle);

void Ptngc_comp_conv_from_rle(unsigned int *rle,
			                  unsigned int *vals, const int nvals);

#endif

/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */


#ifndef VALS16_H
#define VALS16_H

void Ptngc_comp_conv_to_vals16(unsigned int *vals, const int nvals,
                               unsigned int *vals16, int *nvals16);

void Ptngc_comp_conv_from_vals16(unsigned int *vals16, const int nvals16,
                                 unsigned int *vals, int *nvals);

#endif

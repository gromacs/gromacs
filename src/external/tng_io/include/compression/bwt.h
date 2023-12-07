/*
 * This code is part of the tng binary trajectory format.
 *
 * Copyright (c) 2010,2013, The GROMACS development team.
 * Copyright (c) 2020, by the GROMACS development team.
 * TNG was orginally written by Magnus Lundborg, Daniel Spångberg and
 * Rossen Apostolov. The API is implemented mainly by Magnus Lundborg,
 * Daniel Spångberg and Anders Gärdenäs.
 *
 * Please see the AUTHORS file for more information.
 *
 * The TNG library is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 *
 * To help us fund future development, we humbly ask that you cite
 * the research papers on the package.
 *
 * Check out http://www.gromacs.org for more information.
 */

/* This code is part of the tng compression routines
 * Written by Daniel Spangberg
 */

#ifndef BWT_H
#define BWT_H

void Ptngc_comp_to_bwt(unsigned int* vals, int nvals, unsigned int* output, int* index);

void Ptngc_comp_from_bwt(const unsigned int* input, int nvals, int index, unsigned int* vals);

void Ptngc_bwt_merge_sort_inner(int*          indices,
                                int           nvals,
                                unsigned int* vals,
                                int           start,
                                int           end,
                                unsigned int* nrepeat,
                                int*          workarray);

#endif

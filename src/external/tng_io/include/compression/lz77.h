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

#ifndef LZ77_H
#define LZ77_H

void Ptngc_comp_to_lz77(unsigned int* vals,
                        int           nvals,
                        unsigned int* data,
                        int*          ndata,
                        unsigned int* len,
                        int*          nlens,
                        unsigned int* offsets,
                        int*          noffsets);

void Ptngc_comp_from_lz77(const unsigned int* data,
                          int                 ndata,
                          const unsigned int* len,
                          int                 nlens,
                          const unsigned int* offsets,
                          int                 noffsets,
                          unsigned int*       vals,
                          int                 nvals);

#endif

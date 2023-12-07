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

#ifndef MTF_H
#define MTF_H

void Ptngc_comp_conv_to_mtf(const unsigned int* vals,
                            int                 nvals,
                            const unsigned int* dict,
                            int                 ndict,
                            unsigned int*       valsmtf);

void Ptngc_comp_conv_from_mtf(const unsigned int* valsmtf,
                              int                 nvals,
                              const unsigned int* dict,
                              int                 ndict,
                              unsigned int*       vals);

void Ptngc_comp_conv_to_mtf_partial(const unsigned int* vals, int nvals, unsigned int* valsmtf);

void Ptngc_comp_conv_from_mtf_partial(const unsigned int* valsmtf, int nvals, unsigned int* vals);

void Ptngc_comp_conv_to_mtf_partial3(const unsigned int* vals, int nvals, unsigned char* valsmtf);

void Ptngc_comp_conv_from_mtf_partial3(unsigned char* valsmtf, int nvals, unsigned int* vals);

#endif

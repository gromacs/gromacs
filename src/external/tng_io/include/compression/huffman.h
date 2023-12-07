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

#ifndef HUFFMAN_H
#define HUFFMAN_H

void Ptngc_comp_conv_to_huffman(const unsigned int* vals,
                                int                 nvals,
                                const unsigned int* dict,
                                int                 ndict,
                                unsigned int*       prob,
                                unsigned char*      huffman,
                                int*                huffman_len,
                                unsigned char*      huffman_dict,
                                int*                huffman_dictlen,
                                unsigned int*       huffman_dict_unpacked,
                                int*                huffman_dict_unpackedlen);

void Ptngc_comp_conv_from_huffman(unsigned char*      huffman,
                                  unsigned int*       vals,
                                  int                 nvals,
                                  int                 ndict,
                                  unsigned char*      huffman_dict,
                                  int                 huffman_dictlen,
                                  const unsigned int* huffman_dict_unpacked,
                                  int                 huffman_dict_unpackedlen);

#endif

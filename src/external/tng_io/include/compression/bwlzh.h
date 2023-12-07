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

#ifndef BWLZH_H
#define BWLZH_H

/* Compress the integers (positive, small integers are preferable)
   using bwlzh compression.  The unsigned char *output should be
   allocated to be able to hold worst case. You can obtain this length
   conveniently by calling comp_get_buflen()
*/
void DECLSPECDLLEXPORT bwlzh_compress(unsigned int* vals, int nvals, unsigned char* output, int* output_len);

void DECLSPECDLLEXPORT bwlzh_compress_no_lz77(unsigned int* vals, int nvals, unsigned char* output, int* output_len);

int DECLSPECDLLEXPORT bwlzh_get_buflen(int nvals);

void DECLSPECDLLEXPORT bwlzh_decompress(unsigned char* input, int nvals, unsigned int* vals);


/* The routines below are mostly useful for testing, and for internal
   use by the library. */

void DECLSPECDLLEXPORT bwlzh_compress_verbose(unsigned int* vals, int nvals, unsigned char* output, int* output_len);

void DECLSPECDLLEXPORT bwlzh_compress_no_lz77_verbose(unsigned int*  vals,
                                                      int            nvals,
                                                      unsigned char* output,
                                                      int*           output_len);

void DECLSPECDLLEXPORT bwlzh_decompress_verbose(unsigned char* input, int nvals, unsigned int* vals);

/* Compress the integers (positive, small integers are preferable)
   using huffman coding, with automatic selection of how to handle the
   huffman dictionary.  The unsigned char *huffman should be allocated
   to be able to hold worst case. You can obtain this length
   conveniently by calling comp_huff_buflen()
*/
void Ptngc_comp_huff_compress(unsigned int* vals, int nvals, unsigned char* huffman, int* huffman_len);

int Ptngc_comp_huff_buflen(int nvals);

void Ptngc_comp_huff_decompress(unsigned char* huffman, int huffman_len, unsigned int* vals);


/* the value pointed to by chosen_algo should be
   sent as -1 for autodetect. */
void Ptngc_comp_huff_compress_verbose(unsigned int*  vals,
                                      int            nvals,
                                      unsigned char* huffman,
                                      int*           huffman_len,
                                      int*           huffdatalen,
                                      int*           huffman_lengths,
                                      int*           chosen_algo,
                                      int            isvals16);

#define N_HUFFMAN_ALGO 3
char* Ptngc_comp_get_huff_algo_name(int algo);
char* Ptngc_comp_get_algo_name(int algo);


#endif

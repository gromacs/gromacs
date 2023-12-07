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

#ifndef CODER_H
#define CODER_H

#ifndef DECLSPECDLLEXPORT
#    ifdef USE_WINDOWS
#        define DECLSPECDLLEXPORT __declspec(dllexport)
#    else /* USE_WINDOWS */
#        define DECLSPECDLLEXPORT
#    endif /* USE_WINDOWS */
#endif     /* DECLSPECDLLEXPORT */

struct coder
{
    unsigned int pack_temporary;
    int          pack_temporary_bits;
    int          stat_overflow;
    int          stat_numval;
};

struct coder DECLSPECDLLEXPORT* Ptngc_coder_init(void);
void DECLSPECDLLEXPORT Ptngc_coder_deinit(struct coder* coder);
unsigned char DECLSPECDLLEXPORT* Ptngc_pack_array(struct coder* coder,
                                                  int*          input,
                                                  int*          length,
                                                  int           coding,
                                                  int           coding_parameter,
                                                  int           natoms,
                                                  int           speed);
int DECLSPECDLLEXPORT Ptngc_unpack_array(struct coder*  coder,
                                         unsigned char* packed,
                                         int*           output,
                                         int            length,
                                         int            coding,
                                         int            coding_parameter,
                                         int            natoms);
unsigned char DECLSPECDLLEXPORT* Ptngc_pack_array_xtc2(struct coder* coder, int* input, int* length);
int DECLSPECDLLEXPORT Ptngc_unpack_array_xtc2(struct coder* coder, unsigned char* packed, int* output, int length);
unsigned char DECLSPECDLLEXPORT* Ptngc_pack_array_xtc3(int* input, int* length, int natoms, int speed);
int DECLSPECDLLEXPORT Ptngc_unpack_array_xtc3(unsigned char* packed, int* output, int length, int natoms);

void DECLSPECDLLEXPORT Ptngc_out8bits(struct coder* coder, unsigned char** output);
void DECLSPECDLLEXPORT Ptngc_pack_flush(struct coder* coder, unsigned char** output);
void DECLSPECDLLEXPORT Ptngc_write_pattern(struct coder*   coder,
                                           unsigned int    pattern,
                                           int             nbits,
                                           unsigned char** output);

void DECLSPECDLLEXPORT Ptngc_writebits(struct coder* coder, unsigned int value, int nbits, unsigned char** output_ptr);
void DECLSPECDLLEXPORT Ptngc_write32bits(struct coder*   coder,
                                         unsigned int    value,
                                         int             nbits,
                                         unsigned char** output_ptr);
void DECLSPECDLLEXPORT Ptngc_writemanybits(struct coder*   coder,
                                           unsigned char*  value,
                                           int             nbits,
                                           unsigned char** output_ptr);


#endif

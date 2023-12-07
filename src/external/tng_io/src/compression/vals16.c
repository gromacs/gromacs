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

#include "../../include/compression/vals16.h"

/* Coding 32 bit ints in sequences of 16 bit ints. Worst case
   the output is 3*nvals long. */
void Ptngc_comp_conv_to_vals16(const unsigned int* vals, const int nvals, unsigned int* vals16, int* nvals16)
{
    int i;
    int j = 0;
    for (i = 0; i < nvals; i++)
    {
        if (vals[i] <= 0x7FFFU)
        {
            vals16[j++] = vals[i];
        }
        else
        {
            unsigned int lo = (vals[i] & 0x7FFFU) | 0x8000U;
            unsigned int hi = vals[i] >> 15;
            vals16[j++]     = lo;
            if (hi <= 0x7FFFU)
            {
                vals16[j++] = hi;
            }
            else
            {
                unsigned int lohi = (hi & 0x7FFFU) | 0x8000U;
                unsigned int hihi = hi >> 15;
                vals16[j++]       = lohi;
                vals16[j++]       = hihi;
            }
        }
    }
#if 0
  /* Test that things that detect that this is bad really works. */
  vals16[0]=0;
#endif
    *nvals16 = j;
}

void Ptngc_comp_conv_from_vals16(const unsigned int* vals16, const int nvals16, unsigned int* vals, int* nvals)
{
    int i = 0;
    int j = 0;
    while (i < nvals16)
    {
        if (vals16[i] <= 0x7FFFU)
        {
            vals[j++] = vals16[i++];
        }
        else
        {
            unsigned int lo = vals16[i++];
            unsigned int hi = vals16[i++];
            if (hi <= 0x7FFFU)
            {
                vals[j++] = (lo & 0x7FFFU) | (hi << 15);
            }
            else
            {
                unsigned int hihi = vals16[i++];
                vals[j++]         = (lo & 0x7FFFU) | ((hi & 0x7FFFU) << 15) | (hihi << 30);
            }
        }
    }
    *nvals = j;
}

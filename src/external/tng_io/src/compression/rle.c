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

#include "../../include/compression/rle.h"

static void add_rle(unsigned int* rle, const int v, int nsim, int* j, const int min_rle)
{
    if (nsim > min_rle)
    {
        /* Insert run-length */
        unsigned int run = ((unsigned int)nsim);
        while (run > 1)
        {
            if (run & 0x1U)
            {
                rle[(*j)++] = 1;
            }
            else
            {
                rle[(*j)++] = 0;
            }
            run >>= 1;
        }
        nsim = 1;
    }
    while (nsim--)
    {
        rle[(*j)++] = v + 2;
    }
}

/* Run length encoding.
   Acceptable inputs are about 16 bits (0-0xFFFF)
   If input is 0-N output will be be values of 0-(N+2) */
void Ptngc_comp_conv_to_rle(const unsigned int* vals, const int nvals, unsigned int* rle, int* nrle, const int min_rle)
{
    int i;
    int j    = 0;
    int nsim = 0;
    int v    = -1;
    for (i = 0; i < nvals; i++)
    {
        if (!nsim)
        {
            v    = vals[i];
            nsim = 1;
        }
        else
        {
            if (v == vals[i])
            {
                nsim++;
            }
            else
            {
                add_rle(rle, v, nsim, &j, min_rle);
                nsim = 1;
                v    = vals[i];
            }
        }
    }
    if (nsim != 0)
    {
        add_rle(rle, v, nsim, &j, min_rle);
    }
    *nrle = j;
}

void Ptngc_comp_conv_from_rle(const unsigned int* rle, unsigned int* vals, const int nvals)
{
    int i = 0;
    int j = 0;
    while (i < nvals)
    {
        int          k;
        unsigned int len    = 0;
        unsigned int mask   = 0x1;
        unsigned int v      = rle[j++];
        unsigned int hasrle = 0;
        while (v < 2)
        {
            if (v)
            {
                len |= mask;
            }
            mask <<= 1;
            hasrle = 1;
            v      = rle[j++];
        }
        if (!hasrle)
        {
            len = 1;
        }
        else
        {
            len |= mask;
        }
        for (k = 0; k < (int)len; k++)
        {
            vals[i++] = v - 2;
        }
    }
}

/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "nbnxn_kernel_common.h"

static void
clear_f_all(const nbnxn_atomdata_t *nbat, real *f)
{
    int i;

    for (i = 0; i < nbat->natoms*nbat->fstride; i++)
    {
        f[i] = 0;
    }
}

static void
clear_f_flagged(const nbnxn_atomdata_t *nbat, int output_index, real *f)
{
    const nbnxn_buffer_flags_t *flags;
    unsigned                    our_flag;
    int g, b, a0, a1, i;

    flags = &nbat->buffer_flags;

    our_flag = (1U << output_index);

    for (b = 0; b < flags->nflag; b++)
    {
        if (flags->flag[b] & our_flag)
        {
            a0 = b*NBNXN_BUFFERFLAG_SIZE;
            a1 = a0 + NBNXN_BUFFERFLAG_SIZE;
            for (i = a0*nbat->fstride; i < a1*nbat->fstride; i++)
            {
                f[i] = 0;
            }
        }
    }
}

void
clear_f(const nbnxn_atomdata_t *nbat, int output_index, real *f)
{
    if (nbat->bUseBufferFlags)
    {
        clear_f_flagged(nbat, output_index, f);
    }
    else
    {
        clear_f_all(nbat, f);
    }
}

void
clear_fshift(real *fshift)
{
    int i;

    for (i = 0; i < SHIFTS*DIM; i++)
    {
        fshift[i] = 0;
    }
}

void
reduce_energies_over_lists(const nbnxn_atomdata_t     *nbat,
                           int                         nlist,
                           real                       *Vvdw,
                           real                       *Vc)
{
    int nb;
    int i, j, ind, indr;

    for (nb = 0; nb < nlist; nb++)
    {
        for (i = 0; i < nbat->nenergrp; i++)
        {
            /* Reduce the diagonal terms */
            ind        = i*nbat->nenergrp + i;
            Vvdw[ind] += nbat->out[nb].Vvdw[ind];
            Vc[ind]   += nbat->out[nb].Vc[ind];

            /* Reduce the off-diagonal terms */
            for (j = i+1; j < nbat->nenergrp; j++)
            {
                /* The output should contain only one off-diagonal part */
                ind        = i*nbat->nenergrp + j;
                indr       = j*nbat->nenergrp + i;
                Vvdw[ind] += nbat->out[nb].Vvdw[ind] + nbat->out[nb].Vvdw[indr];
                Vc[ind]   += nbat->out[nb].Vc[ind]   + nbat->out[nb].Vc[indr];
            }
        }
    }
}

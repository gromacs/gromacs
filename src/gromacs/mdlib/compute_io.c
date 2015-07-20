/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include "compute_io.h"

#include <signal.h>
#include <stdlib.h>

#include "gromacs/legacyheaders/typedefs.h"

static int div_nsteps(int nsteps, int nst)
{
    if (nst > 0)
    {
        return (1 + nsteps + nst - 1)/nst;
    }
    else
    {
        return 0;
    }
}

double compute_io(t_inputrec *ir, int natoms, gmx_groups_t *groups,
                  int nrener, int nrepl)
{

    int    nsteps = ir->nsteps;
    int    i, nxtcatoms = 0;
    int    nstx, nstv, nstf, nste, nstlog, nstxtc, nfep = 0;
    double cio;

    nstx   = div_nsteps(nsteps, ir->nstxout);
    nstv   = div_nsteps(nsteps, ir->nstvout);
    nstf   = div_nsteps(nsteps, ir->nstfout);
    nstxtc = div_nsteps(nsteps, ir->nstxout_compressed);
    if (ir->nstxout_compressed > 0)
    {
        for (i = 0; i < natoms; i++)
        {
            if (groups->grpnr[egcCompressedX] == NULL || groups->grpnr[egcCompressedX][i] == 0)
            {
                nxtcatoms++;
            }
        }
    }
    nstlog = div_nsteps(nsteps, ir->nstlog);
    /* We add 2 for the header */
    nste   = div_nsteps(2+nsteps, ir->nstenergy);

    cio  = 80*natoms;
    cio += (nstx+nstf+nstv)*sizeof(real)*(3.0*natoms);
    cio += nstxtc*(14*4 + nxtcatoms*5.0); /* roughly 5 bytes per atom */
    cio += nstlog*(nrener*16*2.0);        /* 16 bytes per energy term plus header */
    /* t_energy contains doubles, but real is written to edr */
    cio += (1.0*nste)*nrener*3*sizeof(real);

    if ((ir->efep != efepNO || ir->bSimTemp) && (ir->fepvals->nstdhdl > 0))
    {
        int ndh    = ir->fepvals->n_lambda;
        int ndhdl  = 0;
        int nchars = 0;

        for (i = 0; i < efptNR; i++)
        {
            if (ir->fepvals->separate_dvdl[i])
            {
                ndhdl += 1;
            }
        }

        if (ir->fepvals->separate_dhdl_file == esepdhdlfileYES)
        {
            nchars = 8 + ndhdl*8 + ndh*10; /* time data ~8 chars/entry, dH data ~10 chars/entry */
            if (ir->expandedvals->elmcmove > elmcmoveNO)
            {
                nchars += 5;   /* alchemical state */
            }

            if (ir->fepvals->edHdLPrintEnergy != edHdLPrintEnergyNO)
            {
                nchars += 12; /* energy for dhdl */
            }
            cio += div_nsteps(nsteps, ir->fepvals->nstdhdl)*nchars;
        }
        else
        {
            /* dH output to ener.edr: */
            if (ir->fepvals->dh_hist_size <= 0)
            {
                int ndh_tot = ndh+ndhdl;
                if (ir->expandedvals->elmcmove > elmcmoveNO)
                {
                    ndh_tot += 1;
                }
                if (ir->fepvals->edHdLPrintEnergy != edHdLPrintEnergyNO)
                {
                    ndh_tot += 1;
                }
                /* as data blocks: 1 real per dH point */
                cio += div_nsteps(nsteps, ir->fepvals->nstdhdl)*ndh_tot*sizeof(real);
            }
            else
            {
                /* as histograms: dh_hist_size ints per histogram */
                cio += div_nsteps(nsteps, ir->nstenergy)*
                    sizeof(int)*ir->fepvals->dh_hist_size*ndh;
            }
        }
    }
    if (ir->pull != NULL)
    {
        cio += div_nsteps(nsteps, ir->pull->nstxout)*20; /* roughly 20 chars per line */
        cio += div_nsteps(nsteps, ir->pull->nstfout)*20; /* roughly 20 chars per line */
    }

    return cio*nrepl/(1024*1024);
}

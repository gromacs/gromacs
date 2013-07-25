/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#include "typedefs.h"
#include "xmdrun.h"
#include "vec.h"

real mol_dipole(int k0, int k1, rvec x[], real q[])
{
    int  k, m;
    rvec mu;

    clear_rvec(mu);
    for (k = k0; (k < k1); k++)
    {
        for (m = 0; (m < DIM); m++)
        {
            mu[m] += q[k]*x[k][m];
        }
    }
    return norm(mu); /* Dipole moment of this molecule in e nm */
}

real calc_mu_aver(rvec x[], real q[], t_block *mols, t_mdatoms *md, int gnx, atom_id grpindex[])
{
    int     i, start, end;
    real    mu_ave;

    start = md->start;
    end   = md->homenr + start;

    /*
       clear_rvec(mu);
       for(i=start; (i<end); i++)
       for(m=0; (m<DIM); m++)
        mu[m] += q[i]*x[i][m];
       if (PAR(cr)) {
       gmx_sum(DIM,mu,cr);
       }
     */
    /* I guess we have to parallelise this one! */

    if (gnx > 0)
    {
        mu_ave = 0.0;
        for (i = 0; (i < gnx); i++)
        {
            int gi = grpindex[i];
            mu_ave += mol_dipole(mols->index[gi], mols->index[gi+1], x, q);
        }

        return(mu_ave/gnx);
    }
    else
    {
        return 0;
    }
}


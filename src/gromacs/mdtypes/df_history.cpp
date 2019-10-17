/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2018,2019, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "df_history.h"

#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/smalloc.h"

void init_df_history(df_history_t* dfhist, int nlambda)
{
    int i;

    dfhist->nlambda  = nlambda;
    dfhist->bEquil   = false;
    dfhist->wl_delta = 0;

    if (nlambda > 0)
    {
        snew(dfhist->sum_weights, dfhist->nlambda);
        snew(dfhist->sum_dg, dfhist->nlambda);
        snew(dfhist->sum_minvar, dfhist->nlambda);
        snew(dfhist->sum_variance, dfhist->nlambda);
        snew(dfhist->n_at_lam, dfhist->nlambda);
        snew(dfhist->wl_histo, dfhist->nlambda);

        /* allocate transition matrices here */
        snew(dfhist->Tij, dfhist->nlambda);
        snew(dfhist->Tij_empirical, dfhist->nlambda);

        /* allocate accumulators for various transition matrix
           free energy methods here */
        snew(dfhist->accum_p, dfhist->nlambda);
        snew(dfhist->accum_m, dfhist->nlambda);
        snew(dfhist->accum_p2, dfhist->nlambda);
        snew(dfhist->accum_m2, dfhist->nlambda);

        for (i = 0; i < dfhist->nlambda; i++)
        {
            snew(dfhist->Tij[i], dfhist->nlambda);
            snew(dfhist->Tij_empirical[i], dfhist->nlambda);
            snew((dfhist->accum_p)[i], dfhist->nlambda);
            snew((dfhist->accum_m)[i], dfhist->nlambda);
            snew((dfhist->accum_p2)[i], dfhist->nlambda);
            snew((dfhist->accum_m2)[i], dfhist->nlambda);
        }
    }
}

void copy_df_history(df_history_t* df_dest, df_history_t* df_source)
{
    int i, j;

    /* Currently, there should not be any difference in nlambda between the two,
       but this is included for completeness for potential later functionality */
    df_dest->nlambda  = df_source->nlambda;
    df_dest->bEquil   = df_source->bEquil;
    df_dest->wl_delta = df_source->wl_delta;

    for (i = 0; i < df_dest->nlambda; i++)
    {
        df_dest->sum_weights[i]  = df_source->sum_weights[i];
        df_dest->sum_dg[i]       = df_source->sum_dg[i];
        df_dest->sum_minvar[i]   = df_source->sum_minvar[i];
        df_dest->sum_variance[i] = df_source->sum_variance[i];
        df_dest->n_at_lam[i]     = df_source->n_at_lam[i];
        df_dest->wl_histo[i]     = df_source->wl_histo[i];
    }

    for (i = 0; i < df_dest->nlambda; i++)
    {
        for (j = 0; j < df_dest->nlambda; j++)
        {
            df_dest->accum_p[i][j]       = df_source->accum_p[i][j];
            df_dest->accum_m[i][j]       = df_source->accum_m[i][j];
            df_dest->accum_p2[i][j]      = df_source->accum_p2[i][j];
            df_dest->accum_m2[i][j]      = df_source->accum_m2[i][j];
            df_dest->Tij[i][j]           = df_source->Tij[i][j];
            df_dest->Tij_empirical[i][j] = df_source->Tij_empirical[i][j];
        }
    }
}

void done_df_history(df_history_t* dfhist)
{
    int i;

    if (dfhist->nlambda > 0)
    {
        sfree(dfhist->n_at_lam);
        sfree(dfhist->wl_histo);
        sfree(dfhist->sum_weights);
        sfree(dfhist->sum_dg);
        sfree(dfhist->sum_minvar);
        sfree(dfhist->sum_variance);

        for (i = 0; i < dfhist->nlambda; i++)
        {
            sfree(dfhist->Tij[i]);
            sfree(dfhist->Tij_empirical[i]);
            sfree(dfhist->accum_p[i]);
            sfree(dfhist->accum_m[i]);
            sfree(dfhist->accum_p2[i]);
            sfree(dfhist->accum_m2[i]);
        }
    }
    dfhist->bEquil   = false;
    dfhist->nlambda  = 0;
    dfhist->wl_delta = 0;
}

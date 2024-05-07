/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "mdebin_bar.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstring>

#include <memory>

#include "gromacs/fileio/enxio.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

/* Number of entries in subblock_d preceding the lambda components */
constexpr int c_subblockDNumPreEntries = 5;
/* Number of entries in subblock_i preceding the lambda components */
constexpr int c_subblockINumPreEntries = 2;

/* reset the delta_h list to prepare it for new values */
static void mde_delta_h_reset(t_mde_delta_h* dh)
{
    dh->ndh     = 0;
    dh->written = FALSE;
}

/* initialize the delta_h list */
static void mde_delta_h_init(t_mde_delta_h* dh,
                             int            nbins,
                             double         dx,
                             unsigned int   ndhmax,
                             int            type,
                             int            derivative,
                             int            nlambda,
                             const double*  lambda)
{
    int i;

    dh->type       = type;
    dh->derivative = derivative;
    dh->nlambda    = nlambda;

    dh->lambda.resize(nlambda);
    for (i = 0; i < nlambda; i++)
    {
        assert(lambda);
        dh->lambda[i] = lambda[i];
    }


    dh->subblock_meta_d.resize(dh->nlambda + 1);

    dh->ndhmax = ndhmax + 2;

    dh->dh.resize(dh->ndhmax);
    dh->dhf.resize(dh->ndhmax);

    if (nbins <= 0 || dx < GMX_REAL_EPS * 10)
    {
        dh->nhist = 0;
    }
    else
    {
        int i;
        /* pre-allocate the histogram */
        dh->nhist = 2; /* energies and derivatives histogram */
        dh->dx    = dx;
        dh->nbins = nbins;
        for (i = 0; i < dh->nhist; i++)
        {
            dh->bin[i].resize(dh->nbins);
        }
    }
    mde_delta_h_reset(dh);
}

/* Add a value to the delta_h list */
static void mde_delta_h_add_dh(t_mde_delta_h* dh, double delta_h)
{
    if (dh->ndh >= dh->ndhmax)
    {
        gmx_incons("delta_h array not big enough!");
    }
    dh->dh[dh->ndh] = delta_h;
    dh->ndh++;
}

/* construct histogram with index hi */
static void mde_delta_h_make_hist(t_mde_delta_h* dh, int hi, gmx_bool invert)
{
    double       min_dh = FLT_MAX;
    double       max_dh = -FLT_MAX;
    unsigned int i;
    double       max_dh_hist; /* maximum binnable dh value */
    double       min_dh_hist; /* minimum binnable dh value */
    double       dx = dh->dx;
    double       f; /* energy mult. factor */

    /* by applying a -1 scaling factor on the energies we get the same as
       having a negative dx, but we don't need to fix the min/max values
       beyond inverting x0 */
    f = invert ? -1 : 1;

    /* first find min and max */
    for (i = 0; i < dh->ndh; i++)
    {
        if (f * dh->dh[i] < min_dh)
        {
            min_dh = f * dh->dh[i];
        }
        if (f * dh->dh[i] > max_dh)
        {
            max_dh = f * dh->dh[i];
        }
    }

    /* reset the histogram */
    for (i = 0; i < dh->nbins; i++)
    {
        dh->bin[hi][i] = 0;
    }
    dh->maxbin[hi] = 0;

    /* The starting point of the histogram is the lowest value found:
       that value has the highest contribution to the free energy.

       Get this start value in number of histogram dxs from zero,
       as an integer.*/

    dh->x0[hi] = static_cast<int64_t>(floor(min_dh / dx));

    min_dh_hist = (dh->x0[hi]) * dx;
    max_dh_hist = (dh->x0[hi] + dh->nbins + 1) * dx;

    /* and fill the histogram*/
    for (i = 0; i < dh->ndh; i++)
    {
        unsigned int bin;

        /* Determine the bin number. If it doesn't fit into the histogram,
           add it to the last bin.
           We check the max_dh_int range because converting to integers
           might lead to overflow with unpredictable results.*/
        if ((f * dh->dh[i] >= min_dh_hist) && (f * dh->dh[i] <= max_dh_hist))
        {
            bin = static_cast<unsigned int>((f * dh->dh[i] - min_dh_hist) / dx);
        }
        else
        {
            bin = dh->nbins - 1;
        }
        /* double-check here because of possible round-off errors*/
        if (bin >= dh->nbins)
        {
            bin = dh->nbins - 1;
        }
        if (bin > dh->maxbin[hi])
        {
            dh->maxbin[hi] = bin;
        }

        dh->bin[hi][bin]++;
    }

    /* make sure we include a bin with 0 if we didn't use the full
       histogram width. This can then be used as an indication that
       all the data was binned. */
    if (dh->maxbin[hi] < dh->nbins - 1)
    {
        dh->maxbin[hi] += 1;
    }
}


static void mde_delta_h_handle_block(t_mde_delta_h* dh, t_enxblock* blk)
{
    /* first check which type we should use: histogram or raw data */
    if (dh->nhist == 0)
    {
        int i;

        /* We write raw data.
           Raw data consists of 3 subblocks: an int metadata block
           with type and derivative index, a foreign lambda block
           and and the data itself */
        add_subblocks_enxblock(blk, 3);

        blk->id = enxDH;

        /* subblock 1 */
        dh->subblock_meta_i[0] = dh->type;       /* block data type */
        dh->subblock_meta_i[1] = dh->derivative; /* derivative direction if
                                                    applicable (in indices
                                                    starting from first coord in
                                                    the main delta_h_coll) */
        blk->sub[0].nr   = 2;
        blk->sub[0].type = XdrDataType::Int;
        blk->sub[0].ival = dh->subblock_meta_i.data();

        /* subblock 2 */
        for (i = 0; i < dh->nlambda; i++)
        {
            dh->subblock_meta_d[i] = dh->lambda[i];
        }
        blk->sub[1].nr   = dh->nlambda;
        blk->sub[1].type = XdrDataType::Double;
        blk->sub[1].dval = dh->subblock_meta_d.data();

        /* subblock 3 */
        /* check if there's actual data to be written. */
        /*if (dh->ndh > 1)*/
        if (dh->ndh > 0)
        {
            unsigned int i;

            blk->sub[2].nr = dh->ndh;
            /* Michael commented in 2012 that this use of explicit
               XdrDataType::Float was good for F@H for now.
               Apparently it's still good enough. */
            blk->sub[2].type = XdrDataType::Float;
            for (i = 0; i < dh->ndh; i++)
            {
                dh->dhf[i] = static_cast<float>(dh->dh[i]);
            }
            blk->sub[2].fval = dh->dhf.data();
            dh->written      = TRUE;
        }
        else
        {
            blk->sub[2].nr   = 0;
            blk->sub[2].type = XdrDataType::Float;
            blk->sub[2].fval = nullptr;
        }
    }
    else
    {
        int nhist_written = 0;
        int i;
        int k;

        /* TODO histogram metadata */
        /* check if there's actual data to be written. */
        if (dh->ndh > 1)
        {
            gmx_bool prev_complete = FALSE;
            /* Make the histogram(s) */
            for (i = 0; i < dh->nhist; i++)
            {
                if (!prev_complete)
                {
                    /* the first histogram is always normal, and the
                       second one is always reverse */
                    mde_delta_h_make_hist(dh, i, i == 1);
                    nhist_written++;
                    /* check whether this histogram contains all data: if the
                       last bin is 0, it does */
                    if (dh->bin[i][dh->nbins - 1] == 0)
                    {
                        prev_complete = TRUE;
                    }
                    if (!dh->derivative)
                    {
                        prev_complete = TRUE;
                    }
                }
            }
            dh->written = TRUE;
        }

        /* A histogram consists of 2, 3 or 4 subblocks:
           the foreign lambda value + histogram spacing, the starting point,
           and the histogram data (0, 1 or 2 blocks). */
        add_subblocks_enxblock(blk, nhist_written + 2);
        blk->id = enxDHHIST;

        /* subblock 1: the lambda value + the histogram spacing */
        if (dh->nlambda == 1)
        {
            /* for backward compatibility */
            dh->subblock_meta_d[0] = dh->lambda[0];
        }
        else
        {
            dh->subblock_meta_d[0] = -1;
            for (i = 0; i < dh->nlambda; i++)
            {
                dh->subblock_meta_d[2 + i] = dh->lambda[i];
            }
        }
        dh->subblock_meta_d[1] = dh->dx;
        blk->sub[0].nr         = 2 + ((dh->nlambda > 1) ? dh->nlambda : 0);
        blk->sub[0].type       = XdrDataType::Double;
        blk->sub[0].dval       = dh->subblock_meta_d.data();

        /* subblock 2: the starting point(s) as a long integer */
        dh->subblock_meta_l[0] = nhist_written;
        dh->subblock_meta_l[1] = dh->type; /*dh->derivative ? 1 : 0;*/
        k                      = 2;
        for (i = 0; i < nhist_written; i++)
        {
            dh->subblock_meta_l[k++] = dh->x0[i];
        }
        /* append the derivative data */
        dh->subblock_meta_l[k++] = dh->derivative;

        blk->sub[1].nr   = nhist_written + 3;
        blk->sub[1].type = XdrDataType::Int64;
        blk->sub[1].lval = dh->subblock_meta_l.data();

        /* subblock 3 + 4 : the histogram data */
        for (i = 0; i < nhist_written; i++)
        {
            blk->sub[i + 2].nr = dh->maxbin[i] + 1; /* it's +1 because size=index+1
                                                       in C */
            blk->sub[i + 2].type = XdrDataType::Int;
            blk->sub[i + 2].ival = dh->bin[i].data();
        }
    }
}

/* initialize the collection*/
t_mde_delta_h_coll::t_mde_delta_h_coll(const t_inputrec& inputrec)
{
    int       i, j, n;
    double*   lambda_vec;
    int       ndhmax = inputrec.nstenergy / inputrec.nstcalcenergy;
    t_lambda* fep    = inputrec.fepvals.get();

    /* only store system temperature */
    temperature = (haveConstantEnsembleTemperature(inputrec) ? constantEnsembleTemperature(inputrec)
                                                             : 0.0_real);
    start_time  = 0.;
    delta_time  = inputrec.delta_t * inputrec.fepvals->nstdhdl;
    start_time_set = FALSE;

    /* this is the compatibility lambda value. If it is >=0, it is valid,
       and there is either an old-style lambda or a slow growth simulation. */
    start_lambda = inputrec.fepvals->init_lambda;
    /* for continuous change of lambda values */
    delta_lambda = inputrec.fepvals->delta_lambda * inputrec.fepvals->nstdhdl;

    if (start_lambda < 0)
    {
        /* create the native lambda vectors */
        lambda_index = fep->init_fep_state;
        n_lambda_vec = 0;
        for (auto i : keysOf(fep->separate_dvdl))
        {
            if (fep->separate_dvdl[i])
            {
                n_lambda_vec++;
            }
        }
        native_lambda_vec.resize(n_lambda_vec);
        native_lambda_components.resize(n_lambda_vec);
        j = 0;
        for (auto i : keysOf(fep->separate_dvdl))
        {
            if (fep->separate_dvdl[i])
            {
                native_lambda_components[j] = static_cast<int>(i);
                if (fep->init_fep_state >= 0 && fep->init_fep_state < fep->n_lambda)
                {
                    native_lambda_vec[j] = fep->all_lambda[i][fep->init_fep_state];
                }
                else
                {
                    native_lambda_vec[j] = -1;
                }
                j++;
            }
        }
    }
    else
    {
        /* don't allocate the meta-data subblocks for lambda vectors */
        n_lambda_vec = 0;
        lambda_index = -1;
    }
    /* allocate metadata subblocks */
    subblock_d.resize(c_subblockDNumPreEntries + n_lambda_vec);
    subblock_i.resize(c_subblockINumPreEntries + n_lambda_vec);

    /* now decide which data to write out */
    nlambda           = 0;
    ndhdl             = 0;
    dh_expanded_index = -1;
    dh_energy_index   = -1;
    dh_pv_index       = -1;

    /* total number of raw data point collections in the sample */
    ndh = 0;

    {
        gmx_bool bExpanded           = FALSE;
        gmx_bool bEnergy             = FALSE;
        gmx_bool bPV                 = FALSE;
        int      n_lambda_components = 0;

        /* first count the number of states */

        /* add the dhdl's */
        if (fep->dhdl_derivatives == DhDlDerivativeCalculation::Yes)
        {
            for (auto i : keysOf(fep->separate_dvdl))
            {
                if (inputrec.fepvals->separate_dvdl[i])
                {
                    ndh += 1;
                    ndhdl += 1;
                }
            }
        }
        /* add the lambdas */
        nlambda = inputrec.fepvals->lambda_stop_n - inputrec.fepvals->lambda_start_n;
        ndh += nlambda;
        /* another compatibility check */
        if (start_lambda < 0)
        {
            /* include one more for the specification of the state, by lambda or
               fep_state*/
            if (inputrec.expandedvals->elmcmove > LambdaMoveCalculation::No)
            {
                ndh += 1;
                bExpanded = TRUE;
            }
            /* whether to print energies */
            if (inputrec.fepvals->edHdLPrintEnergy != FreeEnergyPrintEnergy::No)
            {
                ndh += 1;
                bEnergy = TRUE;
            }
            if (inputrec.pressureCouplingOptions.epc > PressureCoupling::No)
            {
                ndh += 1; /* include pressure-volume work */
                bPV = TRUE;
            }
        }
        /* allocate them */
        dh.resize(ndh);

        /* now initialize them */
        /* the order, for now, must match that of the dhdl.xvg file because of
           how gmx energy -odh is implemented */
        n = 0;
        if (bExpanded)
        {
            dh_expanded_index = n;
            mde_delta_h_init(&dh[n],
                             inputrec.fepvals->dh_hist_size,
                             inputrec.fepvals->dh_hist_spacing,
                             ndhmax,
                             dhbtEXPANDED,
                             0,
                             0,
                             nullptr);
            n++;
        }
        if (bEnergy)
        {
            dh_energy_index = n;
            mde_delta_h_init(&dh[n],
                             inputrec.fepvals->dh_hist_size,
                             inputrec.fepvals->dh_hist_spacing,
                             ndhmax,
                             dhbtEN,
                             0,
                             0,
                             nullptr);
            n++;
        }
        /* add the dhdl's */
        n_lambda_components = 0;
        if (fep->dhdl_derivatives == DhDlDerivativeCalculation::Yes)
        {
            dh_dhdl_index = n;
            for (auto i : keysOf(fep->separate_dvdl))
            {
                if (inputrec.fepvals->separate_dvdl[i])
                {
                    /* we give it init_lambda for compatibility */
                    mde_delta_h_init(&dh[n],
                                     inputrec.fepvals->dh_hist_size,
                                     inputrec.fepvals->dh_hist_spacing,
                                     ndhmax,
                                     dhbtDHDL,
                                     n_lambda_components,
                                     1,
                                     &(fep->init_lambda));
                    n++;
                    n_lambda_components++;
                }
            }
        }
        else
        {
            for (auto i : keysOf(fep->separate_dvdl))
            {
                if (fep->separate_dvdl[i])
                {
                    n_lambda_components++; /* count the components */
                }
            }
        }
        /* add the lambdas */
        dh_du_index = n;
        snew(lambda_vec, n_lambda_components);
        for (i = inputrec.fepvals->lambda_start_n; i < inputrec.fepvals->lambda_stop_n; i++)
        {
            int k = 0;

            for (auto j : keysOf(fep->separate_dvdl))
            {
                if (fep->separate_dvdl[j])
                {
                    lambda_vec[k++] = fep->all_lambda[j][i];
                }
            }

            mde_delta_h_init(&dh[n],
                             inputrec.fepvals->dh_hist_size,
                             inputrec.fepvals->dh_hist_spacing,
                             ndhmax,
                             dhbtDH,
                             0,
                             n_lambda_components,
                             lambda_vec);
            n++;
        }
        sfree(lambda_vec);
        if (bPV)
        {
            dh_pv_index = n;
            mde_delta_h_init(&dh[n],
                             inputrec.fepvals->dh_hist_size,
                             inputrec.fepvals->dh_hist_spacing,
                             ndhmax,
                             dhbtPV,
                             0,
                             0,
                             nullptr);
            n++;
        }
    }
}

/* add a bunch of samples - note fep_state is double to allow for better data storage */
void mde_delta_h_coll_add_dh(t_mde_delta_h_coll*   dhc,
                             double                fep_state,
                             double                energy,
                             double                pV,
                             gmx::ArrayRef<double> dhdl,
                             double*               foreign_dU,
                             double                time)
{
    int i;

    if (!dhc->start_time_set)
    {
        dhc->start_time_set = TRUE;
        dhc->start_time     = time;
    }

    for (i = 0; i < dhc->ndhdl; i++)
    {
        mde_delta_h_add_dh(&dhc->dh[dhc->dh_dhdl_index + i], dhdl[i]);
    }
    for (i = 0; i < dhc->nlambda; i++)
    {
        mde_delta_h_add_dh(&dhc->dh[dhc->dh_du_index + i], foreign_dU[i]);
    }
    if (dhc->dh_pv_index >= 0)
    {
        mde_delta_h_add_dh(&dhc->dh[dhc->dh_pv_index], pV);
    }
    if (dhc->dh_energy_index >= 0)
    {
        mde_delta_h_add_dh(&dhc->dh[dhc->dh_energy_index], energy);
    }
    if (dhc->dh_expanded_index >= 0)
    {
        mde_delta_h_add_dh(&dhc->dh[dhc->dh_expanded_index], fep_state);
    }
}

/* write the metadata associated with all the du blocks, and call
   handle_block to write out all the du blocks */
void mde_delta_h_coll_handle_block(t_mde_delta_h_coll* dhc, t_enxframe* fr, int nblock)
{
    int         i;
    t_enxblock* blk;

    /* add one block with one subblock as the collection's own data */
    nblock++;
    add_blocks_enxframe(fr, nblock);
    blk = fr->block + (nblock - 1);

    /* only allocate lambda vector component blocks if they must be written out
       for backward compatibility */
    if (!dhc->native_lambda_components.empty())
    {
        add_subblocks_enxblock(blk, 2);
    }
    else
    {
        add_subblocks_enxblock(blk, 1);
    }

    dhc->subblock_d[0] = dhc->temperature;  /* temperature */
    dhc->subblock_d[1] = dhc->start_time;   /* time of first sample */
    dhc->subblock_d[2] = dhc->delta_time;   /* time difference between samples */
    dhc->subblock_d[3] = dhc->start_lambda; /* old-style lambda at starttime */
    dhc->subblock_d[4] = dhc->delta_lambda; /* lambda diff. between samples */
    /* set the lambda vector components if they exist */
    if (!dhc->native_lambda_components.empty())
    {
        for (i = 0; i < dhc->n_lambda_vec; i++)
        {
            dhc->subblock_d[c_subblockDNumPreEntries + i] = dhc->native_lambda_vec[i];
        }
    }
    blk->id          = enxDHCOLL;
    blk->sub[0].nr   = c_subblockDNumPreEntries + dhc->n_lambda_vec;
    blk->sub[0].type = XdrDataType::Double;
    blk->sub[0].dval = dhc->subblock_d.data();

    if (!dhc->native_lambda_components.empty())
    {
        dhc->subblock_i[0] = dhc->lambda_index;
        /* set the lambda vector component IDs if they exist */
        dhc->subblock_i[1] = dhc->n_lambda_vec;
        for (i = 0; i < dhc->n_lambda_vec; i++)
        {
            dhc->subblock_i[c_subblockINumPreEntries + i] = dhc->native_lambda_components[i];
        }
        blk->sub[1].nr   = c_subblockINumPreEntries + dhc->n_lambda_vec;
        blk->sub[1].type = XdrDataType::Int;
        blk->sub[1].ival = dhc->subblock_i.data();
    }

    for (i = 0; i < dhc->ndh; i++)
    {
        nblock++;
        add_blocks_enxframe(fr, nblock);
        blk = fr->block + (nblock - 1);

        mde_delta_h_handle_block(&dhc->dh[i], blk);
    }
}

/* reset the data for a new round */
void mde_delta_h_coll_reset(t_mde_delta_h_coll* dhc)
{
    int i;
    for (i = 0; i < dhc->ndh; i++)
    {
        if (dhc->dh[i].written)
        {
            /* we can now throw away the data */
            mde_delta_h_reset(&dhc->dh[i]);
        }
    }
    dhc->start_time_set = FALSE;
}

/* set the energyhistory variables to save state */
void mde_delta_h_coll_update_energyhistory(const t_mde_delta_h_coll* dhc, energyhistory_t* enerhist)
{
    if (enerhist->deltaHForeignLambdas == nullptr)
    {
        enerhist->deltaHForeignLambdas = std::make_unique<delta_h_history_t>();
        enerhist->deltaHForeignLambdas->dh.resize(dhc->ndh);
    }

    delta_h_history_t* const deltaH = enerhist->deltaHForeignLambdas.get();

    GMX_RELEASE_ASSERT(
            deltaH->dh.size() == static_cast<size_t>(dhc->ndh),
            "energy history number of delta_h histograms should match inputrec's number");

    for (int i = 0; i < dhc->ndh; i++)
    {
        std::vector<real>& dh = deltaH->dh[i];
        for (unsigned int j = 0; j < dhc->dh[i].ndh; j++)
        {
            dh.emplace_back(dhc->dh[i].dh[j]);
        }
    }
    deltaH->start_time   = dhc->start_time;
    deltaH->start_lambda = dhc->start_lambda;
}


/* restore the variables from an energyhistory */
void mde_delta_h_coll_restore_energyhistory(t_mde_delta_h_coll* dhc, const delta_h_history_t* deltaH)
{
    GMX_RELEASE_ASSERT(dhc, "Should have delta_h histograms");
    GMX_RELEASE_ASSERT(deltaH, "Should have delta_h histograms in energy history");
    GMX_RELEASE_ASSERT(
            deltaH->dh.size() == static_cast<size_t>(dhc->ndh),
            "energy history number of delta_h histograms should match inputrec's number");

    for (gmx::Index i = 0; i < gmx::ssize(deltaH->dh); i++)
    {
        dhc->dh[i].ndh = deltaH->dh[i].size();
        for (unsigned int j = 0; j < dhc->dh[i].ndh; j++)
        {
            dhc->dh[i].dh[j] = deltaH->dh[i][j];
        }
    }
    dhc->start_time = deltaH->start_time;
    if (deltaH->start_lambda_set)
    {
        dhc->start_lambda = deltaH->start_lambda;
    }
    dhc->start_time_set = (dhc->dh[0].ndh > 0);
}

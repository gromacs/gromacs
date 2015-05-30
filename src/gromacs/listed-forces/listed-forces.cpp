/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief This file defines high-level functions for mdrun to compute
 * energies and forces for listed interactions.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_listed-forces
 */
#include "gmxpre.h"

#include "listed-forces.h"

#include "config.h"

#include <assert.h>

#include <algorithm>

#include "gromacs/legacyheaders/disre.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/orires.h"
#include "gromacs/legacyheaders/types/force_flags.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/listed-forces/position-restraints.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/forcerec-threading.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/smalloc.h"

#include "pairs.h"

namespace
{

/*! \brief Return true if ftype is an explicit pair-listed LJ or
 * COULOMB interaction type: bonded LJ (usually 1-4), or special
 * listed non-bonded for FEP. */
bool
isPairInteraction(int ftype)
{
    return ((ftype) >= F_LJ14 && (ftype) <= F_LJC_PAIRS_NB);
}

/*! \brief Zero thread-local force-output buffers */
void
zero_thread_forces(f_thread_t *f_t, int n,
                   int nblock, int blocksize)
{
    int b, a0, a1, a, i, j;

    if (n > f_t->f_nalloc)
    {
        f_t->f_nalloc = over_alloc_large(n);
        srenew(f_t->f, f_t->f_nalloc);
    }

    if (!bitmask_is_zero(f_t->red_mask))
    {
        for (b = 0; b < nblock; b++)
        {
            if (bitmask_is_set(f_t->red_mask, b))
            {
                a0 = b*blocksize;
                a1 = std::min((b+1)*blocksize, n);
                for (a = a0; a < a1; a++)
                {
                    clear_rvec(f_t->f[a]);
                }
            }
        }
    }
    for (i = 0; i < SHIFTS; i++)
    {
        clear_rvec(f_t->fshift[i]);
    }
    for (i = 0; i < F_NRE; i++)
    {
        f_t->ener[i] = 0;
    }
    for (i = 0; i < egNR; i++)
    {
        for (j = 0; j < f_t->grpp.nener; j++)
        {
            f_t->grpp.ener[i][j] = 0;
        }
    }
    for (i = 0; i < efptNR; i++)
    {
        f_t->dvdl[i] = 0;
    }
}

/*! \brief The max thread number is arbitrary, we used a fixed number
 * to avoid memory management.  Using more than 16 threads is probably
 * never useful performance wise. */
#define MAX_BONDED_THREADS 256

/*! \brief Reduce thread-local force buffers */
void
reduce_thread_force_buffer(int n, rvec *f,
                           int nthreads, f_thread_t *f_t,
                           int nblock, int block_size)
{
    int b;

    if (nthreads > MAX_BONDED_THREADS)
    {
        gmx_fatal(FARGS, "Can not reduce bonded forces on more than %d threads",
                  MAX_BONDED_THREADS);
    }

    /* This reduction can run on any number of threads,
     * independently of nthreads.
     */
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (b = 0; b < nblock; b++)
    {
        rvec *fp[MAX_BONDED_THREADS];
        int   nfb, ft, fb;
        int   a0, a1, a;

        /* Determine which threads contribute to this block */
        nfb = 0;
        for (ft = 1; ft < nthreads; ft++)
        {
            if (bitmask_is_set(f_t[ft].red_mask, b))
            {
                fp[nfb++] = f_t[ft].f;
            }
        }
        if (nfb > 0)
        {
            /* Reduce force buffers for threads that contribute */
            a0 =  b   *block_size;
            a1 = (b+1)*block_size;
            a1 = std::min(a1, n);
            for (a = a0; a < a1; a++)
            {
                for (fb = 0; fb < nfb; fb++)
                {
                    rvec_inc(f[a], fp[fb][a]);
                }
            }
        }
    }
}

/*! \brief Reduce thread-local forces */
void
reduce_thread_forces(int n, rvec *f, rvec *fshift,
                     real *ener, gmx_grppairener_t *grpp, real *dvdl,
                     int nthreads, f_thread_t *f_t,
                     int nblock, int block_size,
                     gmx_bool bCalcEnerVir,
                     gmx_bool bDHDL)
{
    if (nblock > 0)
    {
        /* Reduce the bonded force buffer */
        reduce_thread_force_buffer(n, f, nthreads, f_t, nblock, block_size);
    }

    /* When necessary, reduce energy and virial using one thread only */
    if (bCalcEnerVir)
    {
        int t, i, j;

        for (i = 0; i < SHIFTS; i++)
        {
            for (t = 1; t < nthreads; t++)
            {
                rvec_inc(fshift[i], f_t[t].fshift[i]);
            }
        }
        for (i = 0; i < F_NRE; i++)
        {
            for (t = 1; t < nthreads; t++)
            {
                ener[i] += f_t[t].ener[i];
            }
        }
        for (i = 0; i < egNR; i++)
        {
            for (j = 0; j < f_t[1].grpp.nener; j++)
            {
                for (t = 1; t < nthreads; t++)
                {

                    grpp->ener[i][j] += f_t[t].grpp.ener[i][j];
                }
            }
        }
        if (bDHDL)
        {
            for (i = 0; i < efptNR; i++)
            {

                for (t = 1; t < nthreads; t++)
                {
                    dvdl[i] += f_t[t].dvdl[i];
                }
            }
        }
    }
}

/*! \brief Calculate one element of the list of bonded interactions
    for this thread */
real
calc_one_bond(int thread,
              int ftype, const t_idef *idef,
              const rvec x[], rvec f[], rvec fshift[],
              t_forcerec *fr,
              const t_pbc *pbc, const t_graph *g,
              gmx_grppairener_t *grpp,
              t_nrnb *nrnb,
              real *lambda, real *dvdl,
              const t_mdatoms *md, t_fcdata *fcd,
              gmx_bool bCalcEnerVir,
              int *global_atom_index)
{
#ifdef GMX_SIMD_HAVE_REAL
    gmx_bool bUseSIMD;
    /* MSVC 2010 produces buggy SIMD PBC code, disable SIMD for MSVC <= 2010 */
#if defined _MSC_VER && _MSC_VER < 1700 && !defined(__ICL)
    bUseSIMD = FALSE;
#else
    bUseSIMD = fr->use_simd_kernels;
#endif
#endif

    int      nat1, nbonds, efptFTYPE;
    real     v = 0;
    t_iatom *iatoms;
    int      nb0, nbn;

    if (IS_RESTRAINT_TYPE(ftype))
    {
        efptFTYPE = efptRESTRAINT;
    }
    else
    {
        efptFTYPE = efptBONDED;
    }

    nat1      = interaction_function[ftype].nratoms + 1;
    nbonds    = idef->il[ftype].nr/nat1;
    iatoms    = idef->il[ftype].iatoms;

    nb0 = idef->il_thread_division[ftype*(idef->nthreads+1)+thread];
    nbn = idef->il_thread_division[ftype*(idef->nthreads+1)+thread+1] - nb0;

    if (!isPairInteraction(ftype))
    {
        if (ftype == F_CMAP)
        {
            /* TODO The execution time for CMAP dihedrals might be
               nice to account to its own subtimer, but first
               wallcycle needs to be extended to support calling from
               multiple threads. */
            v = cmap_dihs(nbn, iatoms+nb0,
                          idef->iparams, &idef->cmap_grid,
                          x, f, fshift,
                          pbc, g, lambda[efptFTYPE], &(dvdl[efptFTYPE]),
                          md, fcd, global_atom_index);
        }
#ifdef GMX_SIMD_HAVE_REAL
        else if (ftype == F_ANGLES && bUseSIMD &&
                 !bCalcEnerVir && fr->efep == efepNO)
        {
            /* No energies, shift forces, dvdl */
            angles_noener_simd(nbn, idef->il[ftype].iatoms+nb0,
                               idef->iparams,
                               x, f,
                               pbc, g, lambda[efptFTYPE], md, fcd,
                               global_atom_index);
            v = 0;
        }
#endif
        else if (ftype == F_PDIHS &&
                 !bCalcEnerVir && fr->efep == efepNO)
        {
            /* No energies, shift forces, dvdl */
#ifdef GMX_SIMD_HAVE_REAL
            if (bUseSIMD)
            {
                pdihs_noener_simd(nbn, idef->il[ftype].iatoms+nb0,
                                  idef->iparams,
                                  x, f,
                                  pbc, g, lambda[efptFTYPE], md, fcd,
                                  global_atom_index);
            }
            else
#endif
            {
                pdihs_noener(nbn, idef->il[ftype].iatoms+nb0,
                             idef->iparams,
                             x, f,
                             pbc, g, lambda[efptFTYPE], md, fcd,
                             global_atom_index);
            }
            v = 0;
        }
#ifdef GMX_SIMD_HAVE_REAL
        else if (ftype == F_RBDIHS && bUseSIMD &&
                 !bCalcEnerVir && fr->efep == efepNO)
        {
            /* No energies, shift forces, dvdl */
            rbdihs_noener_simd(nbn, idef->il[ftype].iatoms+nb0,
                               idef->iparams,
                               (const rvec*)x, f,
                               pbc, g, lambda[efptFTYPE], md, fcd,
                               global_atom_index);
            v = 0;
        }
#endif
        else
        {
            v = interaction_function[ftype].ifunc(nbn, iatoms+nb0,
                                                  idef->iparams,
                                                  x, f, fshift,
                                                  pbc, g, lambda[efptFTYPE], &(dvdl[efptFTYPE]),
                                                  md, fcd, global_atom_index);
        }
    }
    else
    {
        /* TODO The execution time for pairs might be nice to account
           to its own subtimer, but first wallcycle needs to be
           extended to support calling from multiple threads. */
        v = do_pairs(ftype, nbn, iatoms+nb0, idef->iparams, x, f, fshift,
                     pbc, g, lambda, dvdl, md, fr, grpp, global_atom_index);
    }

    if (thread == 0)
    {
        inc_nrnb(nrnb, interaction_function[ftype].nrnb_ind, nbonds);
    }

    return v;
}

} // namespace

gmx_bool
ftype_is_bonded_potential(int ftype)
{
    return
        (interaction_function[ftype].flags & IF_BOND) &&
        !(ftype == F_CONNBONDS || ftype == F_POSRES || ftype == F_FBPOSRES) &&
        (ftype < F_GB12 || ftype > F_GB14);
}

void calc_listed(const gmx_multisim_t *ms,
                 gmx_wallcycle        *wcycle,
                 const t_idef *idef,
                 const rvec x[], history_t *hist,
                 rvec f[], t_forcerec *fr,
                 const struct t_pbc *pbc,
                 const struct t_pbc *pbc_full,
                 const struct t_graph *g,
                 gmx_enerdata_t *enerd, t_nrnb *nrnb,
                 real *lambda,
                 const t_mdatoms *md,
                 t_fcdata *fcd, int *global_atom_index,
                 int force_flags)
{
    gmx_bool      bCalcEnerVir;
    int           i;
    real          dvdl[efptNR]; /* The dummy array is to have a place to store the dhdl at other values
                                                        of lambda, which will be thrown away in the end*/
    const  t_pbc *pbc_null;
    int           thread;

    assert(fr->nthreads == idef->nthreads);

    bCalcEnerVir = (force_flags & (GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY));

    for (i = 0; i < efptNR; i++)
    {
        dvdl[i] = 0.0;
    }
    if (fr->bMolPBC)
    {
        pbc_null = pbc;
    }
    else
    {
        pbc_null = NULL;
    }

#ifdef DEBUG
    if (g && debug)
    {
        p_graph(debug, "Bondage is fun", g);
    }
#endif

    if ((idef->il[F_POSRES].nr > 0) ||
        (idef->il[F_FBPOSRES].nr > 0) ||
        (idef->il[F_ORIRES].nr > 0) ||
        (idef->il[F_DISRES].nr > 0))
    {
        /* TODO Use of restraints triggers further function calls
           inside the loop over calc_one_bond(), but those are too
           awkward to account to this subtimer properly in the present
           code. We don't test / care much about performance with
           restraints, anyway. */
        wallcycle_sub_start(wcycle, ewcsRESTRAINTS);

        if (idef->il[F_POSRES].nr > 0)
        {
            posres_wrapper(nrnb, idef, pbc_full, x, enerd, lambda, fr);
        }

        if (idef->il[F_FBPOSRES].nr > 0)
        {
            fbposres_wrapper(nrnb, idef, pbc_full, x, enerd, fr);
        }

        /* Do pre force calculation stuff which might require communication */
        if (idef->il[F_ORIRES].nr > 0)
        {
            enerd->term[F_ORIRESDEV] =
                calc_orires_dev(ms, idef->il[F_ORIRES].nr,
                                idef->il[F_ORIRES].iatoms,
                                idef->iparams, md, x,
                                pbc_null, fcd, hist);
        }
        if (idef->il[F_DISRES].nr)
        {
            calc_disres_R_6(idef->il[F_DISRES].nr,
                            idef->il[F_DISRES].iatoms,
                            idef->iparams, x, pbc_null,
                            fcd, hist);
#ifdef GMX_MPI
            if (fcd->disres.nsystems > 1)
            {
                gmx_sum_sim(2*fcd->disres.nres, fcd->disres.Rt_6, ms);
            }
#endif
        }

        wallcycle_sub_stop(wcycle, ewcsRESTRAINTS);
    }

    wallcycle_sub_start(wcycle, ewcsLISTED);
#pragma omp parallel for num_threads(fr->nthreads) schedule(static)
    for (thread = 0; thread < fr->nthreads; thread++)
    {
        int                ftype;
        real              *epot, v;
        /* thread stuff */
        rvec              *ft, *fshift;
        real              *dvdlt;
        gmx_grppairener_t *grpp;

        if (thread == 0)
        {
            ft     = f;
            fshift = fr->fshift;
            epot   = enerd->term;
            grpp   = &enerd->grpp;
            dvdlt  = dvdl;
        }
        else
        {
            zero_thread_forces(&fr->f_t[thread], fr->natoms_force,
                               fr->red_nblock, 1<<fr->red_ashift);

            ft     = fr->f_t[thread].f;
            fshift = fr->f_t[thread].fshift;
            epot   = fr->f_t[thread].ener;
            grpp   = &fr->f_t[thread].grpp;
            dvdlt  = fr->f_t[thread].dvdl;
        }
        /* Loop over all bonded force types to calculate the bonded forces */
        for (ftype = 0; (ftype < F_NRE); ftype++)
        {
            if (idef->il[ftype].nr > 0 && ftype_is_bonded_potential(ftype))
            {
                v = calc_one_bond(thread, ftype, idef, x,
                                  ft, fshift, fr, pbc_null, g, grpp,
                                  nrnb, lambda, dvdlt,
                                  md, fcd, bCalcEnerVir,
                                  global_atom_index);
                epot[ftype] += v;
            }
        }
    }
    wallcycle_sub_stop(wcycle, ewcsLISTED);

    if (fr->nthreads > 1)
    {
        wallcycle_sub_start(wcycle, ewcsLISTED_BUF_OPS);
        reduce_thread_forces(fr->natoms_force, f, fr->fshift,
                             enerd->term, &enerd->grpp, dvdl,
                             fr->nthreads, fr->f_t,
                             fr->red_nblock, 1<<fr->red_ashift,
                             bCalcEnerVir,
                             force_flags & GMX_FORCE_DHDL);
        wallcycle_sub_stop(wcycle, ewcsLISTED_BUF_OPS);
    }

    /* Remaining code does not have enough flops to bother counting */
    if (force_flags & GMX_FORCE_DHDL)
    {
        for (i = 0; i < efptNR; i++)
        {
            enerd->dvdl_nonlin[i] += dvdl[i];
        }
    }

    /* Copy the sum of violations for the distance restraints from fcd */
    if (fcd)
    {
        enerd->term[F_DISRESVIOL] = fcd->disres.sumviol;

    }
}

void calc_listed_lambda(const t_idef *idef,
                        const rvec x[],
                        t_forcerec *fr,
                        const struct t_pbc *pbc, const struct t_graph *g,
                        gmx_grppairener_t *grpp, real *epot, t_nrnb *nrnb,
                        real *lambda,
                        const t_mdatoms *md,
                        t_fcdata *fcd,
                        int *global_atom_index)
{
    int           ftype, nr_nonperturbed, nr;
    real          v;
    real          dvdl_dum[efptNR] = {0};
    rvec         *f, *fshift;
    const  t_pbc *pbc_null;
    t_idef        idef_fe;

    if (fr->bMolPBC)
    {
        pbc_null = pbc;
    }
    else
    {
        pbc_null = NULL;
    }

    /* Copy the whole idef, so we can modify the contents locally */
    idef_fe          = *idef;
    idef_fe.nthreads = 1;
    snew(idef_fe.il_thread_division, F_NRE*(idef_fe.nthreads+1));

    /* We already have the forces, so we use temp buffers here */
    snew(f, fr->natoms_force);
    snew(fshift, SHIFTS);

    /* Loop over all bonded force types to calculate the bonded energies */
    for (ftype = 0; (ftype < F_NRE); ftype++)
    {
        if (ftype_is_bonded_potential(ftype))
        {
            /* Set the work range of thread 0 to the perturbed bondeds only */
            nr_nonperturbed                       = idef->il[ftype].nr_nonperturbed;
            nr                                    = idef->il[ftype].nr;
            idef_fe.il_thread_division[ftype*2+0] = nr_nonperturbed;
            idef_fe.il_thread_division[ftype*2+1] = nr;

            /* This is only to get the flop count correct */
            idef_fe.il[ftype].nr = nr - nr_nonperturbed;

            if (nr - nr_nonperturbed > 0)
            {
                v = calc_one_bond(0, ftype, &idef_fe,
                                  x, f, fshift, fr, pbc_null, g,
                                  grpp, nrnb, lambda, dvdl_dum,
                                  md, fcd, TRUE,
                                  global_atom_index);
                epot[ftype] += v;
            }
        }
    }

    sfree(fshift);
    sfree(f);

    sfree(idef_fe.il_thread_division);
}

void
do_force_listed(gmx_wallcycle        *wcycle,
                matrix                box,
                const t_lambda       *fepvals,
                const gmx_multisim_t *ms,
                const t_idef         *idef,
                const rvec            x[],
                history_t            *hist,
                rvec                  f[],
                t_forcerec           *fr,
                const struct t_pbc   *pbc,
                const struct t_graph *graph,
                gmx_enerdata_t       *enerd,
                t_nrnb               *nrnb,
                real                 *lambda,
                const t_mdatoms      *md,
                t_fcdata             *fcd,
                int                  *global_atom_index,
                int                   flags)
{
    t_pbc pbc_full; /* Full PBC is needed for position restraints */

    if (!(flags & GMX_FORCE_LISTED))
    {
        return;
    }

    if ((idef->il[F_POSRES].nr > 0) ||
        (idef->il[F_FBPOSRES].nr > 0))
    {
        /* Not enough flops to bother counting */
        set_pbc(&pbc_full, fr->ePBC, box);
    }
    calc_listed(ms, wcycle, idef, x, hist, f, fr, pbc, &pbc_full,
                graph, enerd, nrnb, lambda, md, fcd,
                global_atom_index, flags);

    /* Check if we have to determine energy differences
     * at foreign lambda's.
     */
    if (fepvals->n_lambda > 0 && (flags & GMX_FORCE_DHDL))
    {
        posres_wrapper_lambda(wcycle, fepvals, idef, &pbc_full, x, enerd, lambda, fr);

        if (idef->ilsort != ilsortNO_FE)
        {
            wallcycle_sub_start(wcycle, ewcsLISTED_FEP);
            if (idef->ilsort != ilsortFE_SORTED)
            {
                gmx_incons("The bonded interactions are not sorted for free energy");
            }
            for (int i = 0; i < enerd->n_lambda; i++)
            {
                real lam_i[efptNR];

                reset_foreign_enerdata(enerd);
                for (int j = 0; j < efptNR; j++)
                {
                    lam_i[j] = (i == 0 ? lambda[j] : fepvals->all_lambda[j][i-1]);
                }
                calc_listed_lambda(idef, x, fr, pbc, graph, &(enerd->foreign_grpp), enerd->foreign_term, nrnb, lam_i, md,
                                   fcd, global_atom_index);
                sum_epot(&(enerd->foreign_grpp), enerd->foreign_term);
                enerd->enerpart_lambda[i] += enerd->foreign_term[F_EPOT];
            }
            wallcycle_sub_stop(wcycle, ewcsLISTED_FEP);
        }
    }
    debug_gmx();
}

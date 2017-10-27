/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <assert.h>

#include <algorithm>

#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/listed-forces/disre.h"
#include "gromacs/listed-forces/orires.h"
#include "gromacs/listed-forces/pairs.h"
#include "gromacs/listed-forces/position-restraints.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/simd/simd.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "listed-internal.h"

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

/*! \brief Zero thread-local output buffers */
static void
zero_thread_output(struct bonded_threading_t *bt, int thread)
{
    if (!bt->haveBondeds)
    {
        return;
    }

    f_thread_t *f_t      = &bt->f_t[thread];
    const int   nelem_fa = sizeof(*f_t->f)/sizeof(real);

    for (int i = 0; i < f_t->nblock_used; i++)
    {
        int a0 = f_t->block_index[i]*reduction_block_size;
        int a1 = a0 + reduction_block_size;
        for (int a = a0; a < a1; a++)
        {
            for (int d = 0; d < nelem_fa; d++)
            {
                f_t->f[a][d] = 0;
            }
        }
    }

    for (int i = 0; i < SHIFTS; i++)
    {
        clear_rvec(f_t->fshift[i]);
    }
    for (int i = 0; i < F_NRE; i++)
    {
        f_t->ener[i] = 0;
    }
    for (int i = 0; i < egNR; i++)
    {
        for (int j = 0; j < f_t->grpp.nener; j++)
        {
            f_t->grpp.ener[i][j] = 0;
        }
    }
    for (int i = 0; i < efptNR; i++)
    {
        f_t->dvdl[i] = 0;
    }
}

/*! \brief The max thread number is arbitrary, we used a fixed number
 * to avoid memory management.  Using more than 16 threads is probably
 * never useful performance wise. */
#define MAX_BONDED_THREADS 256

/*! \brief Reduce thread-local force buffers */
static void
reduce_thread_forces(int n, rvec *f,
                     struct bonded_threading_t *bt,
                     int nthreads)
{
    if (nthreads > MAX_BONDED_THREADS)
    {
        gmx_fatal(FARGS, "Can not reduce bonded forces on more than %d threads",
                  MAX_BONDED_THREADS);
    }

    /* This reduction can run on any number of threads,
     * independently of bt->nthreads.
     * But if nthreads matches bt->nthreads (which it currently does)
     * the uniform distribution of the touched blocks over nthreads will
     * match the distribution of bonded over threads well in most cases,
     * which means that threads mostly reduce their own data which increases
     * the number of cache hits.
     */
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int b = 0; b < bt->nblock_used; b++)
    {
        try
        {
            int    ind = bt->block_index[b];
            rvec4 *fp[MAX_BONDED_THREADS];

            /* Determine which threads contribute to this block */
            int nfb = 0;
            for (int ft = 0; ft < bt->nthreads; ft++)
            {
                if (bitmask_is_set(bt->mask[ind], ft))
                {
                    fp[nfb++] = bt->f_t[ft].f;
                }
            }
            if (nfb > 0)
            {
                /* Reduce force buffers for threads that contribute */
                int a0 =  ind     *reduction_block_size;
                int a1 = (ind + 1)*reduction_block_size;
                /* It would be nice if we could pad f to avoid this min */
                a1     = std::min(a1, n);
                for (int a = a0; a < a1; a++)
                {
                    for (int fb = 0; fb < nfb; fb++)
                    {
                        rvec_inc(f[a], fp[fb][a]);
                    }
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
}

/*! \brief Reduce thread-local forces, shift forces and energies */
static void
reduce_thread_output(int n, rvec *f, rvec *fshift,
                     real *ener, gmx_grppairener_t *grpp, real *dvdl,
                     struct bonded_threading_t *bt,
                     gmx_bool bCalcEnerVir,
                     gmx_bool bDHDL)
{
    assert(bt->haveBondeds);

    if (bt->nblock_used > 0)
    {
        /* Reduce the bonded force buffer */
        reduce_thread_forces(n, f, bt, bt->nthreads);
    }

    /* When necessary, reduce energy and virial using one thread only */
    if (bCalcEnerVir && bt->nthreads > 1)
    {
        f_thread_t *f_t = bt->f_t;

        for (int i = 0; i < SHIFTS; i++)
        {
            for (int t = 1; t < bt->nthreads; t++)
            {
                rvec_inc(fshift[i], f_t[t].fshift[i]);
            }
        }
        for (int i = 0; i < F_NRE; i++)
        {
            for (int t = 1; t < bt->nthreads; t++)
            {
                ener[i] += f_t[t].ener[i];
            }
        }
        for (int i = 0; i < egNR; i++)
        {
            for (int j = 0; j < f_t[1].grpp.nener; j++)
            {
                for (int t = 1; t < bt->nthreads; t++)
                {
                    grpp->ener[i][j] += f_t[t].grpp.ener[i][j];
                }
            }
        }
        if (bDHDL)
        {
            for (int i = 0; i < efptNR; i++)
            {

                for (int t = 1; t < bt->nthreads; t++)
                {
                    dvdl[i] += f_t[t].dvdl[i];
                }
            }
        }
    }
}

/*! \brief Calculate one element of the list of bonded interactions
    for this thread */
static real
calc_one_bond(int thread,
              int ftype, const t_idef *idef,
              const rvec x[], rvec4 f[], rvec fshift[],
              const t_forcerec *fr,
              const t_pbc *pbc, const t_graph *g,
              gmx_grppairener_t *grpp,
              t_nrnb *nrnb,
              const real *lambda, real *dvdl,
              const t_mdatoms *md, t_fcdata *fcd,
              gmx_bool bCalcEnerVir,
              int *global_atom_index)
{
#if GMX_SIMD_HAVE_REAL
    bool bUseSIMD = fr->use_simd_kernels;
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
                          md, fcd, global_atom_index
#ifdef BUILD_WITH_FDA
                          , fr->fda
#endif
                         );
        }
#if GMX_SIMD_HAVE_REAL
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

        else if (ftype == F_UREY_BRADLEY && bUseSIMD &&
                 !bCalcEnerVir && fr->efep == efepNO)
        {
            /* No energies, shift forces, dvdl */
            urey_bradley_noener_simd(nbn, idef->il[ftype].iatoms+nb0,
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
#if GMX_SIMD_HAVE_REAL
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
#if GMX_SIMD_HAVE_REAL
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
                                                  md, fcd, global_atom_index
#ifdef BUILD_WITH_FDA
                                                  , fr->fda
#endif
                                                 );
        }
    }
    else
    {
        /* TODO The execution time for pairs might be nice to account
           to its own subtimer, but first wallcycle needs to be
           extended to support calling from multiple threads. */
        do_pairs(ftype, nbn, iatoms+nb0, idef->iparams, x, f, fshift,
                 pbc, g, lambda, dvdl, md, fr,
                 bCalcEnerVir, grpp, global_atom_index);
        v = 0;
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

/*! \brief Compute the bonded part of the listed forces, parallelized over threads
 */
static void
calcBondedForces(const t_idef     *idef,
                 const rvec        x[],
                 const t_forcerec *fr,
                 const t_pbc      *pbc_null,
                 const t_graph    *g,
                 gmx_enerdata_t   *enerd,
                 t_nrnb           *nrnb,
                 const real       *lambda,
                 real             *dvdl,
                 const t_mdatoms  *md,
                 t_fcdata         *fcd,
                 gmx_bool          bCalcEnerVir,
                 int              *global_atom_index)
{
    struct bonded_threading_t *bt = fr->bonded_threading;

#pragma omp parallel for num_threads(bt->nthreads) schedule(static)
    for (int thread = 0; thread < bt->nthreads; thread++)
    {
        try
        {
            int                ftype;
            real              *epot, v;
            /* thread stuff */
            rvec4             *ft;
            rvec              *fshift;
            real              *dvdlt;
            gmx_grppairener_t *grpp;

            zero_thread_output(bt, thread);

            ft = bt->f_t[thread].f;

            if (thread == 0)
            {
                fshift = fr->fshift;
                epot   = enerd->term;
                grpp   = &enerd->grpp;
                dvdlt  = dvdl;
            }
            else
            {
                fshift = bt->f_t[thread].fshift;
                epot   = bt->f_t[thread].ener;
                grpp   = &bt->f_t[thread].grpp;
                dvdlt  = bt->f_t[thread].dvdl;
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
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
}

void calc_listed(const t_commrec             *cr,
                 struct gmx_wallcycle        *wcycle,
                 const t_idef *idef,
                 const rvec x[], history_t *hist,
                 rvec f[],
                 gmx::ForceWithVirial *forceWithVirial,
                 const t_forcerec *fr,
                 const struct t_pbc *pbc,
                 const struct t_pbc *pbc_full,
                 const struct t_graph *g,
                 gmx_enerdata_t *enerd, t_nrnb *nrnb,
                 const real *lambda,
                 const t_mdatoms *md,
                 t_fcdata *fcd, int *global_atom_index,
                 int force_flags)
{
    struct bonded_threading_t *bt;
    gmx_bool                   bCalcEnerVir;
    const  t_pbc              *pbc_null;

    bt = fr->bonded_threading;

    assert(bt->nthreads == idef->nthreads);

    bCalcEnerVir = (force_flags & (GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY));

    if (fr->bMolPBC)
    {
        pbc_null = pbc;
    }
    else
    {
        pbc_null = nullptr;
    }

    if ((idef->il[F_POSRES].nr > 0) ||
        (idef->il[F_FBPOSRES].nr > 0) ||
        fcd->orires.nr > 0 ||
        fcd->disres.nres > 0)
    {
        /* TODO Use of restraints triggers further function calls
           inside the loop over calc_one_bond(), but those are too
           awkward to account to this subtimer properly in the present
           code. We don't test / care much about performance with
           restraints, anyway. */
        wallcycle_sub_start(wcycle, ewcsRESTRAINTS);

        if (idef->il[F_POSRES].nr > 0)
        {
            posres_wrapper(nrnb, idef, pbc_full, x, enerd, lambda, fr,
                           forceWithVirial);
        }

        if (idef->il[F_FBPOSRES].nr > 0)
        {
            fbposres_wrapper(nrnb, idef, pbc_full, x, enerd, fr,
                             forceWithVirial);
        }

        /* Do pre force calculation stuff which might require communication */
        if (fcd->orires.nr > 0)
        {
            /* This assertion is to ensure we have whole molecules.
             * Unfortunately we do not have an mdrun state variable that tells
             * us if molecules in x are not broken over PBC, so we have to make
             * do with checking graph!=nullptr, which should tell us if we made
             * molecules whole before calling the current function.
             */
            GMX_RELEASE_ASSERT(fr->ePBC == epbcNONE || g != nullptr, "With orientation restraints molecules should be whole");
            enerd->term[F_ORIRESDEV] =
                calc_orires_dev(cr->ms, idef->il[F_ORIRES].nr,
                                idef->il[F_ORIRES].iatoms,
                                idef->iparams, md, x,
                                pbc_null, fcd, hist);
        }
        if (fcd->disres.nres > 0)
        {
            calc_disres_R_6(cr,
                            idef->il[F_DISRES].nr,
                            idef->il[F_DISRES].iatoms,
                            x, pbc_null,
                            fcd, hist);
        }

        wallcycle_sub_stop(wcycle, ewcsRESTRAINTS);
    }

    if (bt->haveBondeds)
    {
        wallcycle_sub_start(wcycle, ewcsLISTED);
        /* The dummy array is to have a place to store the dhdl at other values
           of lambda, which will be thrown away in the end */
        real dvdl[efptNR] = {0};
        calcBondedForces(idef, x, fr, pbc_null, g, enerd, nrnb, lambda, dvdl, md,
                         fcd, bCalcEnerVir, global_atom_index);
        wallcycle_sub_stop(wcycle, ewcsLISTED);

        wallcycle_sub_start(wcycle, ewcsLISTED_BUF_OPS);
        reduce_thread_output(fr->natoms_force, f, fr->fshift,
                             enerd->term, &enerd->grpp, dvdl,
                             bt,
                             bCalcEnerVir,
                             force_flags & GMX_FORCE_DHDL);

        if (force_flags & GMX_FORCE_DHDL)
        {
            for (int i = 0; i < efptNR; i++)
            {
                enerd->dvdl_nonlin[i] += dvdl[i];
            }
        }
        wallcycle_sub_stop(wcycle, ewcsLISTED_BUF_OPS);
    }

    /* Copy the sum of violations for the distance restraints from fcd */
    if (fcd)
    {
        enerd->term[F_DISRESVIOL] = fcd->disres.sumviol;
    }
}

void calc_listed_lambda(const t_idef *idef,
                        const rvec x[],
                        const t_forcerec *fr,
                        const struct t_pbc *pbc, const struct t_graph *g,
                        gmx_grppairener_t *grpp, real *epot, t_nrnb *nrnb,
                        const real *lambda,
                        const t_mdatoms *md,
                        t_fcdata *fcd,
                        int *global_atom_index)
{
    int           ftype, nr_nonperturbed, nr;
    real          v;
    real          dvdl_dum[efptNR] = {0};
    rvec4        *f;
    rvec         *fshift;
    const  t_pbc *pbc_null;
    t_idef        idef_fe;

    if (fr->bMolPBC)
    {
        pbc_null = pbc;
    }
    else
    {
        pbc_null = nullptr;
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
do_force_listed(struct gmx_wallcycle        *wcycle,
                matrix                       box,
                const t_lambda              *fepvals,
                const t_commrec             *cr,
                const t_idef                *idef,
                const rvec                   x[],
                history_t                   *hist,
                rvec                        *forceForUseWithShiftForces,
                gmx::ForceWithVirial        *forceWithVirial,
                const t_forcerec            *fr,
                const struct t_pbc          *pbc,
                const struct t_graph        *graph,
                gmx_enerdata_t              *enerd,
                t_nrnb                      *nrnb,
                const real                  *lambda,
                const t_mdatoms             *md,
                t_fcdata                    *fcd,
                int                         *global_atom_index,
                int                          flags)
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
    calc_listed(cr, wcycle, idef, x, hist,
                forceForUseWithShiftForces, forceWithVirial,
                fr, pbc, &pbc_full,
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
}

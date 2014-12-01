/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2006,2007,2008,2009,2010,2012,2013,2014, by the GROMACS development team, led by
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
 * \brief This file implements functions for domdec to use
 * while managing inter-atomic constraints.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "domdec_constraints.h"

#include <assert.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/legacyheaders/constr.h"
#include "gromacs/legacyheaders/gmx_ga2la.h"
#include "gromacs/legacyheaders/gmx_hash.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "domdec_specatomcomm.h"

/*! \brief Struct used during constraint setup with domain decomposition */
typedef struct gmx_domdec_constraints {
    //! @cond Doxygen_Suppress
    int       *molb_con_offset; /**< Offset in the constraint array for each molblock */
    int       *molb_ncon_mol;   /**< The number of constraints per molecule for each molblock */

    int        ncon;            /**< The fully local and conneced constraints */
    /* The global constraint number, only required for clearing gc_req */
    int       *con_gl;          /**< Global constraint indices for local constraints */
    int       *con_nlocat;      /**< Number of local atoms (2/1/0) for each constraint */
    int        con_nalloc;      /**< Allocation size for \p con_gl and \p con_nlocat */

    char      *gc_req;          /**< Boolean that tells if a global constraint index has been requested; note: size global #constraints */
    gmx_hash_t ga2la;           /**< Global to local communicated constraint atom only index */

    /* Multi-threading stuff */
    int      nthread;           /**< Number of threads used for DD constraint setup */
    t_ilist *ils;               /**< Constraint ilist working arrays, size \p nthread */
    //! @endcond
} gmx_domdec_constraints_t;

void dd_move_x_constraints(gmx_domdec_t *dd, matrix box,
                           rvec *x0, rvec *x1, gmx_bool bX1IsCoord)
{
    if (dd->constraint_comm)
    {
        dd_move_x_specat(dd, dd->constraint_comm, box, x0, x1, bX1IsCoord);
    }
}

int *dd_constraints_nlocalatoms(gmx_domdec_t *dd)
{
    if (dd->constraints)
    {
        return dd->constraints->con_nlocat;
    }
    else
    {
        return NULL;
    }
}

void dd_clear_local_constraint_indices(gmx_domdec_t *dd)
{
    gmx_domdec_constraints_t *dc;
    int i;

    dc = dd->constraints;

    for (i = 0; i < dc->ncon; i++)
    {
        dc->gc_req[dc->con_gl[i]] = 0;
    }

    if (dd->constraint_comm)
    {
        gmx_hash_clear_and_optimize(dc->ga2la);
    }
}

void dd_clear_local_vsite_indices(gmx_domdec_t *dd)
{
    if (dd->vsite_comm)
    {
        gmx_hash_clear_and_optimize(dd->ga2la_vsite);
    }
}

/*! \brief Walks over the constraints out from the local atoms into the non-local atoms and adds them to a list */
static void walk_out(int con, int con_offset, int a, int offset, int nrec,
                     int ncon1, const t_iatom *ia1, const t_iatom *ia2,
                     const t_blocka *at2con,
                     const gmx_ga2la_t ga2la, gmx_bool bHomeConnect,
                     gmx_domdec_constraints_t *dc,
                     gmx_domdec_specat_comm_t *dcc,
                     t_ilist *il_local,
                     ind_req_t *ireq)
{
    int            a1_gl, a2_gl, a_loc, i, coni, b;
    const t_iatom *iap;

    if (dc->gc_req[con_offset+con] == 0)
    {
        /* Add this non-home constraint to the list */
        if (dc->ncon+1 > dc->con_nalloc)
        {
            dc->con_nalloc = over_alloc_large(dc->ncon+1);
            srenew(dc->con_gl, dc->con_nalloc);
            srenew(dc->con_nlocat, dc->con_nalloc);
        }
        dc->con_gl[dc->ncon]       = con_offset + con;
        dc->con_nlocat[dc->ncon]   = (bHomeConnect ? 1 : 0);
        dc->gc_req[con_offset+con] = 1;
        if (il_local->nr + 3 > il_local->nalloc)
        {
            il_local->nalloc = over_alloc_dd(il_local->nr+3);
            srenew(il_local->iatoms, il_local->nalloc);
        }
        iap = constr_iatomptr(ncon1, ia1, ia2, con);
        il_local->iatoms[il_local->nr++] = iap[0];
        a1_gl = offset + iap[1];
        a2_gl = offset + iap[2];
        /* The following indexing code can probably be optizimed */
        if (ga2la_get_home(ga2la, a1_gl, &a_loc))
        {
            il_local->iatoms[il_local->nr++] = a_loc;
        }
        else
        {
            /* We set this index later */
            il_local->iatoms[il_local->nr++] = -a1_gl - 1;
        }
        if (ga2la_get_home(ga2la, a2_gl, &a_loc))
        {
            il_local->iatoms[il_local->nr++] = a_loc;
        }
        else
        {
            /* We set this index later */
            il_local->iatoms[il_local->nr++] = -a2_gl - 1;
        }
        dc->ncon++;
    }
    /* Check to not ask for the same atom more than once */
    if (gmx_hash_get_minone(dc->ga2la, offset+a) == -1)
    {
        assert(dcc);
        /* Add this non-home atom to the list */
        if (ireq->n+1 > ireq->nalloc)
        {
            ireq->nalloc = over_alloc_large(ireq->n+1);
            srenew(ireq->ind, ireq->nalloc);
        }
        ireq->ind[ireq->n++] = offset + a;
        /* Temporarily mark with -2, we get the index later */
        gmx_hash_set(dc->ga2la, offset+a, -2);
    }

    if (nrec > 0)
    {
        for (i = at2con->index[a]; i < at2con->index[a+1]; i++)
        {
            coni = at2con->a[i];
            if (coni != con)
            {
                /* Walk further */
                iap = constr_iatomptr(ncon1, ia1, ia2, coni);
                if (a == iap[1])
                {
                    b = iap[2];
                }
                else
                {
                    b = iap[1];
                }
                if (!ga2la_get_home(ga2la, offset+b, &a_loc))
                {
                    walk_out(coni, con_offset, b, offset, nrec-1,
                             ncon1, ia1, ia2, at2con,
                             ga2la, FALSE, dc, dcc, il_local, ireq);
                }
            }
        }
    }
}

/*! \brief Looks up SETTLE constraints for a range of charge-groups */
static void atoms_to_settles(gmx_domdec_t *dd,
                             const gmx_mtop_t *mtop,
                             const int *cginfo,
                             const int **at2settle_mt,
                             int cg_start, int cg_end,
                             t_ilist *ils_local,
                             ind_req_t *ireq)
{
    gmx_ga2la_t           ga2la;
    gmx_mtop_atomlookup_t alook;
    int                   settle;
    int                   nral, sa;
    int                   cg, a, a_gl, a_glsa, a_gls[3], a_locs[3];
    int                   mb, molnr, a_mol, offset;
    const gmx_molblock_t *molb;
    const t_iatom        *ia1;
    gmx_bool              a_home[3];
    int                   nlocal;
    gmx_bool              bAssign;

    ga2la  = dd->ga2la;

    alook = gmx_mtop_atomlookup_settle_init(mtop);

    nral = NRAL(F_SETTLE);

    for (cg = cg_start; cg < cg_end; cg++)
    {
        if (GET_CGINFO_SETTLE(cginfo[cg]))
        {
            for (a = dd->cgindex[cg]; a < dd->cgindex[cg+1]; a++)
            {
                a_gl = dd->gatindex[a];

                gmx_mtop_atomnr_to_molblock_ind(alook, a_gl, &mb, &molnr, &a_mol);
                molb = &mtop->molblock[mb];

                settle = at2settle_mt[molb->type][a_mol];

                if (settle >= 0)
                {
                    offset = a_gl - a_mol;

                    ia1 = mtop->moltype[molb->type].ilist[F_SETTLE].iatoms;

                    bAssign = FALSE;
                    nlocal  = 0;
                    for (sa = 0; sa < nral; sa++)
                    {
                        a_glsa     = offset + ia1[settle*(1+nral)+1+sa];
                        a_gls[sa]  = a_glsa;
                        a_home[sa] = ga2la_get_home(ga2la, a_glsa, &a_locs[sa]);
                        if (a_home[sa])
                        {
                            if (nlocal == 0 && a_gl == a_glsa)
                            {
                                bAssign = TRUE;
                            }
                            nlocal++;
                        }
                    }

                    if (bAssign)
                    {
                        if (ils_local->nr+1+nral > ils_local->nalloc)
                        {
                            ils_local->nalloc = over_alloc_dd(ils_local->nr+1+nral);
                            srenew(ils_local->iatoms, ils_local->nalloc);
                        }

                        ils_local->iatoms[ils_local->nr++] = ia1[settle*4];

                        for (sa = 0; sa < nral; sa++)
                        {
                            if (ga2la_get_home(ga2la, a_gls[sa], &a_locs[sa]))
                            {
                                ils_local->iatoms[ils_local->nr++] = a_locs[sa];
                            }
                            else
                            {
                                ils_local->iatoms[ils_local->nr++] = -a_gls[sa] - 1;
                                /* Add this non-home atom to the list */
                                if (ireq->n+1 > ireq->nalloc)
                                {
                                    ireq->nalloc = over_alloc_large(ireq->n+1);
                                    srenew(ireq->ind, ireq->nalloc);
                                }
                                ireq->ind[ireq->n++] = a_gls[sa];
                                /* A check on double atom requests is
                                 * not required for settle.
                                 */
                            }
                        }
                    }
                }
            }
        }
    }

    gmx_mtop_atomlookup_destroy(alook);
}

/*! \brief Looks up constraint for the local atoms */
static void atoms_to_constraints(gmx_domdec_t *dd,
                                 const gmx_mtop_t *mtop,
                                 const int *cginfo,
                                 const t_blocka *at2con_mt, int nrec,
                                 t_ilist *ilc_local,
                                 ind_req_t *ireq)
{
    const t_blocka           *at2con;
    gmx_ga2la_t               ga2la;
    gmx_mtop_atomlookup_t     alook;
    int                       ncon1;
    gmx_molblock_t           *molb;
    t_iatom                  *ia1, *ia2, *iap;
    int                       nhome, cg, a, a_gl, a_mol, a_loc, b_lo, offset, mb, molnr, b_mol, i, con, con_offset;
    gmx_domdec_constraints_t *dc;
    gmx_domdec_specat_comm_t *dcc;

    dc  = dd->constraints;
    dcc = dd->constraint_comm;

    ga2la  = dd->ga2la;

    alook = gmx_mtop_atomlookup_init(mtop);

    nhome = 0;
    for (cg = 0; cg < dd->ncg_home; cg++)
    {
        if (GET_CGINFO_CONSTR(cginfo[cg]))
        {
            for (a = dd->cgindex[cg]; a < dd->cgindex[cg+1]; a++)
            {
                a_gl = dd->gatindex[a];

                gmx_mtop_atomnr_to_molblock_ind(alook, a_gl, &mb, &molnr, &a_mol);
                molb = &mtop->molblock[mb];

                ncon1 = mtop->moltype[molb->type].ilist[F_CONSTR].nr/NRAL(F_SETTLE);

                ia1 = mtop->moltype[molb->type].ilist[F_CONSTR].iatoms;
                ia2 = mtop->moltype[molb->type].ilist[F_CONSTRNC].iatoms;

                /* Calculate the global constraint number offset for the molecule.
                 * This is only required for the global index to make sure
                 * that we use each constraint only once.
                 */
                con_offset =
                    dc->molb_con_offset[mb] + molnr*dc->molb_ncon_mol[mb];

                /* The global atom number offset for this molecule */
                offset = a_gl - a_mol;
                at2con = &at2con_mt[molb->type];
                for (i = at2con->index[a_mol]; i < at2con->index[a_mol+1]; i++)
                {
                    con = at2con->a[i];
                    iap = constr_iatomptr(ncon1, ia1, ia2, con);
                    if (a_mol == iap[1])
                    {
                        b_mol = iap[2];
                    }
                    else
                    {
                        b_mol = iap[1];
                    }
                    if (ga2la_get_home(ga2la, offset+b_mol, &a_loc))
                    {
                        /* Add this fully home constraint at the first atom */
                        if (a_mol < b_mol)
                        {
                            if (dc->ncon+1 > dc->con_nalloc)
                            {
                                dc->con_nalloc = over_alloc_large(dc->ncon+1);
                                srenew(dc->con_gl, dc->con_nalloc);
                                srenew(dc->con_nlocat, dc->con_nalloc);
                            }
                            dc->con_gl[dc->ncon]     = con_offset + con;
                            dc->con_nlocat[dc->ncon] = 2;
                            if (ilc_local->nr + 3 > ilc_local->nalloc)
                            {
                                ilc_local->nalloc = over_alloc_dd(ilc_local->nr + 3);
                                srenew(ilc_local->iatoms, ilc_local->nalloc);
                            }
                            b_lo = a_loc;
                            ilc_local->iatoms[ilc_local->nr++] = iap[0];
                            ilc_local->iatoms[ilc_local->nr++] = (a_gl == iap[1] ? a    : b_lo);
                            ilc_local->iatoms[ilc_local->nr++] = (a_gl == iap[1] ? b_lo : a   );
                            dc->ncon++;
                            nhome++;
                        }
                    }
                    else
                    {
                        /* We need the nrec constraints coupled to this constraint,
                         * so we need to walk out of the home cell by nrec+1 atoms,
                         * since already atom bg is not locally present.
                         * Therefore we call walk_out with nrec recursions to go
                         * after this first call.
                         */
                        walk_out(con, con_offset, b_mol, offset, nrec,
                                 ncon1, ia1, ia2, at2con,
                                 dd->ga2la, TRUE, dc, dcc, ilc_local, ireq);
                    }
                }
            }
        }
    }

    gmx_mtop_atomlookup_destroy(alook);

    if (debug)
    {
        fprintf(debug,
                "Constraints: home %3d border %3d atoms: %3d\n",
                nhome, dc->ncon-nhome,
                dd->constraint_comm ? ireq->n : 0);
    }
}

int dd_make_local_constraints(gmx_domdec_t *dd, int at_start,
                              const gmx_mtop_t *mtop,
                              const int *cginfo,
                              gmx_constr_t constr, int nrec,
                              t_ilist *il_local)
{
    gmx_domdec_constraints_t *dc;
    t_ilist                  *ilc_local, *ils_local;
    ind_req_t                *ireq;
    const t_blocka           *at2con_mt;
    const int               **at2settle_mt;
    gmx_hash_t                ga2la_specat;
    int at_end, i, j;
    t_iatom                  *iap;

    // This code should not be called unless this condition is true,
    // because that's the only time init_domdec_constraints is
    // called...
    GMX_RELEASE_ASSERT(dd->bInterCGcons || dd->bInterCGsettles, "dd_make_local_constraints called when there are no local constraints");
    // ... and init_domdec_constraints always sets
    // dd->constraint_comm...
    GMX_RELEASE_ASSERT(dd->constraint_comm, "Invalid use of dd_make_local_constraints before construction of constraint_comm");
    // ... which static analysis needs to be reassured about, because
    // otherwise, when dd->bInterCGsettles is
    // true. dd->constraint_comm is unilaterally dereferenced before
    // the call to atoms_to_settles.

    dc = dd->constraints;

    ilc_local = &il_local[F_CONSTR];
    ils_local = &il_local[F_SETTLE];

    dc->ncon      = 0;
    ilc_local->nr = 0;
    if (dd->constraint_comm)
    {
        at2con_mt = atom2constraints_moltype(constr);
        ireq      = &dd->constraint_comm->ireq[0];
        ireq->n   = 0;
    }
    else
    {
        // Currently unreachable
        at2con_mt = NULL;
        ireq      = NULL;
    }

    if (dd->bInterCGsettles)
    {
        at2settle_mt  = atom2settle_moltype(constr);
        ils_local->nr = 0;
    }
    else
    {
        /* Settle works inside charge groups, we assigned them already */
        at2settle_mt = NULL;
    }

    if (at2settle_mt == NULL)
    {
        atoms_to_constraints(dd, mtop, cginfo, at2con_mt, nrec,
                             ilc_local, ireq);
    }
    else
    {
        int t0_set;
        int thread;

        /* Do the constraints, if present, on the first thread.
         * Do the settles on all other threads.
         */
        t0_set = ((at2con_mt != NULL && dc->nthread > 1) ? 1 : 0);

#pragma omp parallel for num_threads(dc->nthread) schedule(static)
        for (thread = 0; thread < dc->nthread; thread++)
        {
            if (at2con_mt && thread == 0)
            {
                atoms_to_constraints(dd, mtop, cginfo, at2con_mt, nrec,
                                     ilc_local, ireq);
            }

            if (thread >= t0_set)
            {
                int        cg0, cg1;
                t_ilist   *ilst;
                ind_req_t *ireqt;

                /* Distribute the settle check+assignments over
                 * dc->nthread or dc->nthread-1 threads.
                 */
                cg0 = (dd->ncg_home*(thread-t0_set  ))/(dc->nthread-t0_set);
                cg1 = (dd->ncg_home*(thread-t0_set+1))/(dc->nthread-t0_set);

                if (thread == t0_set)
                {
                    ilst = ils_local;
                }
                else
                {
                    ilst = &dc->ils[thread];
                }
                ilst->nr = 0;

                ireqt = &dd->constraint_comm->ireq[thread];
                if (thread > 0)
                {
                    ireqt->n = 0;
                }

                atoms_to_settles(dd, mtop, cginfo, at2settle_mt,
                                 cg0, cg1,
                                 ilst, ireqt);
            }
        }

        /* Combine the generate settles and requested indices */
        for (thread = 1; thread < dc->nthread; thread++)
        {
            t_ilist   *ilst;
            ind_req_t *ireqt;
            int        ia;

            if (thread > t0_set)
            {
                ilst = &dc->ils[thread];
                if (ils_local->nr + ilst->nr > ils_local->nalloc)
                {
                    ils_local->nalloc = over_alloc_large(ils_local->nr + ilst->nr);
                    srenew(ils_local->iatoms, ils_local->nalloc);
                }
                for (ia = 0; ia < ilst->nr; ia++)
                {
                    ils_local->iatoms[ils_local->nr+ia] = ilst->iatoms[ia];
                }
                ils_local->nr += ilst->nr;
            }

            ireqt = &dd->constraint_comm->ireq[thread];
            if (ireq->n+ireqt->n > ireq->nalloc)
            {
                ireq->nalloc = over_alloc_large(ireq->n+ireqt->n);
                srenew(ireq->ind, ireq->nalloc);
            }
            for (ia = 0; ia < ireqt->n; ia++)
            {
                ireq->ind[ireq->n+ia] = ireqt->ind[ia];
            }
            ireq->n += ireqt->n;
        }

        if (debug)
        {
            fprintf(debug, "Settles: total %3d\n", ils_local->nr/4);
        }
    }

    if (dd->constraint_comm)
    {
        int nral1;

        at_end =
            setup_specat_communication(dd, ireq, dd->constraint_comm,
                                       dd->constraints->ga2la,
                                       at_start, 2,
                                       "constraint", " or lincs-order");

        /* Fill in the missing indices */
        ga2la_specat = dd->constraints->ga2la;

        nral1 = 1 + NRAL(F_CONSTR);
        for (i = 0; i < ilc_local->nr; i += nral1)
        {
            iap = ilc_local->iatoms + i;
            for (j = 1; j < nral1; j++)
            {
                if (iap[j] < 0)
                {
                    iap[j] = gmx_hash_get_minone(ga2la_specat, -iap[j]-1);
                }
            }
        }

        nral1 = 1 + NRAL(F_SETTLE);
        for (i = 0; i < ils_local->nr; i += nral1)
        {
            iap = ils_local->iatoms + i;
            for (j = 1; j < nral1; j++)
            {
                if (iap[j] < 0)
                {
                    iap[j] = gmx_hash_get_minone(ga2la_specat, -iap[j]-1);
                }
            }
        }
    }
    else
    {
        // Currently unreachable
        at_end = at_start;
    }

    return at_end;
}

void init_domdec_constraints(gmx_domdec_t *dd,
                             gmx_mtop_t   *mtop)
{
    gmx_domdec_constraints_t *dc;
    gmx_molblock_t           *molb;
    int mb, ncon, c;

    if (debug)
    {
        fprintf(debug, "Begin init_domdec_constraints\n");
    }

    snew(dd->constraints, 1);
    dc = dd->constraints;

    snew(dc->molb_con_offset, mtop->nmolblock);
    snew(dc->molb_ncon_mol, mtop->nmolblock);

    ncon = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb                    = &mtop->molblock[mb];
        dc->molb_con_offset[mb] = ncon;
        dc->molb_ncon_mol[mb]   =
            mtop->moltype[molb->type].ilist[F_CONSTR].nr/3 +
            mtop->moltype[molb->type].ilist[F_CONSTRNC].nr/3;
        ncon += molb->nmol*dc->molb_ncon_mol[mb];
    }

    if (ncon > 0)
    {
        snew(dc->gc_req, ncon);
        for (c = 0; c < ncon; c++)
        {
            dc->gc_req[c] = 0;
        }
    }

    /* Use a hash table for the global to local index.
     * The number of keys is a rough estimate, it will be optimized later.
     */
    dc->ga2la = gmx_hash_init(std::min(mtop->natoms/20,
                                       mtop->natoms/(2*dd->nnodes)));

    dc->nthread = gmx_omp_nthreads_get(emntDomdec);
    snew(dc->ils, dc->nthread);

    dd->constraint_comm = specat_comm_init(dc->nthread);
}

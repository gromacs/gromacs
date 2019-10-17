/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2006,2007,2008,2009,2010,2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include <cassert>

#include <algorithm>
#include <memory>

#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/domdec/hashedmap.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h" // only for GET_CGINFO_*
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "domdec_internal.h"
#include "domdec_specatomcomm.h"

/*! \brief Struct used during constraint setup with domain decomposition */
struct gmx_domdec_constraints_t
{
    //! @cond Doxygen_Suppress
    std::vector<int> molb_con_offset; /**< Offset in the constraint array for each molblock */
    std::vector<int> molb_ncon_mol; /**< The number of constraints per molecule for each molblock */

    int ncon; /**< The fully local and conneced constraints */
    /* The global constraint number, only required for clearing gc_req */
    std::vector<int> con_gl;     /**< Global constraint indices for local constraints */
    std::vector<int> con_nlocat; /**< Number of local atoms (2/1/0) for each constraint */

    std::vector<bool> gc_req; /**< Boolean that tells if a global constraint index has been requested; note: size global #constraints */

    /* Hash table for keeping track of requests */
    std::unique_ptr<gmx::HashedMap<int>> ga2la; /**< Global to local communicated constraint atom only index */

    /* Multi-threading stuff */
    int                  nthread; /**< Number of threads used for DD constraint setup */
    std::vector<t_ilist> ils;     /**< Constraint ilist working arrays, size \p nthread */

    /* Buffers for requesting atoms */
    std::vector<std::vector<int>> requestedGlobalAtomIndices; /**< Buffers for requesting global atom indices, one per thread */

    //! @endcond
};

void dd_move_x_constraints(gmx_domdec_t* dd, const matrix box, rvec* x0, rvec* x1, gmx_bool bX1IsCoord)
{
    if (dd->constraint_comm)
    {
        dd_move_x_specat(dd, dd->constraint_comm, box, x0, x1, bX1IsCoord);

        ddReopenBalanceRegionCpu(dd);
    }
}

gmx::ArrayRef<const int> dd_constraints_nlocalatoms(const gmx_domdec_t* dd)
{
    if (dd && dd->constraints)
    {
        return dd->constraints->con_nlocat;
    }
    else
    {
        return {};
    }
}

void dd_clear_local_constraint_indices(gmx_domdec_t* dd)
{
    gmx_domdec_constraints_t* dc = dd->constraints;

    std::fill(dc->gc_req.begin(), dc->gc_req.end(), false);

    if (dd->constraint_comm)
    {
        dc->ga2la->clear();
    }
}

/*! \brief Walks over the constraints out from the local atoms into the non-local atoms and adds them to a list */
static void walk_out(int                       con,
                     int                       con_offset,
                     int                       a,
                     int                       offset,
                     int                       nrec,
                     gmx::ArrayRef<const int>  ia1,
                     gmx::ArrayRef<const int>  ia2,
                     const t_blocka*           at2con,
                     const gmx_ga2la_t&        ga2la,
                     gmx_bool                  bHomeConnect,
                     gmx_domdec_constraints_t* dc,
                     gmx_domdec_specat_comm_t* dcc,
                     t_ilist*                  il_local,
                     std::vector<int>*         ireq)
{
    int            a1_gl, a2_gl, i, coni, b;
    const t_iatom* iap;

    if (!dc->gc_req[con_offset + con])
    {
        /* Add this non-home constraint to the list */
        dc->con_gl.push_back(con_offset + con);
        dc->con_nlocat.push_back(bHomeConnect ? 1 : 0);
        dc->gc_req[con_offset + con] = true;
        if (il_local->nr + 3 > il_local->nalloc)
        {
            il_local->nalloc = over_alloc_dd(il_local->nr + 3);
            srenew(il_local->iatoms, il_local->nalloc);
        }
        iap                              = constr_iatomptr(ia1, ia2, con);
        il_local->iatoms[il_local->nr++] = iap[0];
        a1_gl                            = offset + iap[1];
        a2_gl                            = offset + iap[2];
        /* The following indexing code can probably be optizimed */
        if (const int* a_loc = ga2la.findHome(a1_gl))
        {
            il_local->iatoms[il_local->nr++] = *a_loc;
        }
        else
        {
            /* We set this index later */
            il_local->iatoms[il_local->nr++] = -a1_gl - 1;
        }
        if (const int* a_loc = ga2la.findHome(a2_gl))
        {
            il_local->iatoms[il_local->nr++] = *a_loc;
        }
        else
        {
            /* We set this index later */
            il_local->iatoms[il_local->nr++] = -a2_gl - 1;
        }
        dc->ncon++;
    }
    /* Check to not ask for the same atom more than once */
    if (!dc->ga2la->find(offset + a))
    {
        assert(dcc);
        /* Add this non-home atom to the list */
        ireq->push_back(offset + a);
        /* Temporarily mark with -2, we get the index later */
        dc->ga2la->insert(offset + a, -2);
    }

    if (nrec > 0)
    {
        for (i = at2con->index[a]; i < at2con->index[a + 1]; i++)
        {
            coni = at2con->a[i];
            if (coni != con)
            {
                /* Walk further */
                iap = constr_iatomptr(ia1, ia2, coni);
                if (a == iap[1])
                {
                    b = iap[2];
                }
                else
                {
                    b = iap[1];
                }
                if (!ga2la.findHome(offset + b))
                {
                    walk_out(coni, con_offset, b, offset, nrec - 1, ia1, ia2, at2con, ga2la, FALSE,
                             dc, dcc, il_local, ireq);
                }
            }
        }
    }
}

/*! \brief Looks up SETTLE constraints for a range of charge-groups */
static void atoms_to_settles(gmx_domdec_t*                         dd,
                             const gmx_mtop_t*                     mtop,
                             const int*                            cginfo,
                             gmx::ArrayRef<const std::vector<int>> at2settle_mt,
                             int                                   cg_start,
                             int                                   cg_end,
                             t_ilist*                              ils_local,
                             std::vector<int>*                     ireq)
{
    const gmx_ga2la_t& ga2la = *dd->ga2la;
    int                nral  = NRAL(F_SETTLE);

    int mb = 0;
    for (int a = cg_start; a < cg_end; a++)
    {
        if (GET_CGINFO_SETTLE(cginfo[a]))
        {
            int a_gl = dd->globalAtomIndices[a];
            int a_mol;
            mtopGetMolblockIndex(mtop, a_gl, &mb, nullptr, &a_mol);

            const gmx_molblock_t* molb   = &mtop->molblock[mb];
            int                   settle = at2settle_mt[molb->type][a_mol];

            if (settle >= 0)
            {
                int offset = a_gl - a_mol;

                const int* ia1 = mtop->moltype[molb->type].ilist[F_SETTLE].iatoms.data();

                int      a_gls[3];
                gmx_bool bAssign = FALSE;
                int      nlocal  = 0;
                for (int sa = 0; sa < nral; sa++)
                {
                    int a_glsa = offset + ia1[settle * (1 + nral) + 1 + sa];
                    a_gls[sa]  = a_glsa;
                    if (ga2la.findHome(a_glsa))
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
                    if (ils_local->nr + 1 + nral > ils_local->nalloc)
                    {
                        ils_local->nalloc = over_alloc_dd(ils_local->nr + 1 + nral);
                        srenew(ils_local->iatoms, ils_local->nalloc);
                    }

                    ils_local->iatoms[ils_local->nr++] = ia1[settle * 4];

                    for (int sa = 0; sa < nral; sa++)
                    {
                        if (const int* a_loc = ga2la.findHome(a_gls[sa]))
                        {
                            ils_local->iatoms[ils_local->nr++] = *a_loc;
                        }
                        else
                        {
                            ils_local->iatoms[ils_local->nr++] = -a_gls[sa] - 1;
                            /* Add this non-home atom to the list */
                            ireq->push_back(a_gls[sa]);
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

/*! \brief Looks up constraint for the local atoms */
static void atoms_to_constraints(gmx_domdec_t*                 dd,
                                 const gmx_mtop_t*             mtop,
                                 const int*                    cginfo,
                                 gmx::ArrayRef<const t_blocka> at2con_mt,
                                 int                           nrec,
                                 t_ilist*                      ilc_local,
                                 std::vector<int>*             ireq)
{
    const t_blocka* at2con;
    int             b_lo, offset, b_mol, i, con, con_offset;

    gmx_domdec_constraints_t* dc  = dd->constraints;
    gmx_domdec_specat_comm_t* dcc = dd->constraint_comm;

    const gmx_ga2la_t& ga2la = *dd->ga2la;

    dc->con_gl.clear();
    dc->con_nlocat.clear();

    int mb    = 0;
    int nhome = 0;
    for (int a = 0; a < dd->ncg_home; a++)
    {
        if (GET_CGINFO_CONSTR(cginfo[a]))
        {
            int a_gl = dd->globalAtomIndices[a];
            int molnr, a_mol;
            mtopGetMolblockIndex(mtop, a_gl, &mb, &molnr, &a_mol);

            const gmx_molblock_t& molb = mtop->molblock[mb];

            gmx::ArrayRef<const int> ia1 = mtop->moltype[molb.type].ilist[F_CONSTR].iatoms;
            gmx::ArrayRef<const int> ia2 = mtop->moltype[molb.type].ilist[F_CONSTRNC].iatoms;

            /* Calculate the global constraint number offset for the molecule.
             * This is only required for the global index to make sure
             * that we use each constraint only once.
             */
            con_offset = dc->molb_con_offset[mb] + molnr * dc->molb_ncon_mol[mb];

            /* The global atom number offset for this molecule */
            offset = a_gl - a_mol;
            at2con = &at2con_mt[molb.type];
            for (i = at2con->index[a_mol]; i < at2con->index[a_mol + 1]; i++)
            {
                con            = at2con->a[i];
                const int* iap = constr_iatomptr(ia1, ia2, con);
                if (a_mol == iap[1])
                {
                    b_mol = iap[2];
                }
                else
                {
                    b_mol = iap[1];
                }
                if (const int* a_loc = ga2la.findHome(offset + b_mol))
                {
                    /* Add this fully home constraint at the first atom */
                    if (a_mol < b_mol)
                    {
                        dc->con_gl.push_back(con_offset + con);
                        dc->con_nlocat.push_back(2);
                        if (ilc_local->nr + 3 > ilc_local->nalloc)
                        {
                            ilc_local->nalloc = over_alloc_dd(ilc_local->nr + 3);
                            srenew(ilc_local->iatoms, ilc_local->nalloc);
                        }
                        b_lo                               = *a_loc;
                        ilc_local->iatoms[ilc_local->nr++] = iap[0];
                        ilc_local->iatoms[ilc_local->nr++] = (a_gl == iap[1] ? a : b_lo);
                        ilc_local->iatoms[ilc_local->nr++] = (a_gl == iap[1] ? b_lo : a);
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
                    walk_out(con, con_offset, b_mol, offset, nrec, ia1, ia2, at2con, ga2la, TRUE,
                             dc, dcc, ilc_local, ireq);
                }
            }
        }
    }

    GMX_ASSERT(dc->con_gl.size() == static_cast<size_t>(dc->ncon),
               "con_gl size should match the number of constraints");
    GMX_ASSERT(dc->con_nlocat.size() == static_cast<size_t>(dc->ncon),
               "con_nlocat size should match the number of constraints");

    if (debug)
    {
        fprintf(debug, "Constraints: home %3d border %3d atoms: %3zu\n", nhome, dc->ncon - nhome,
                dd->constraint_comm ? ireq->size() : 0);
    }
}

int dd_make_local_constraints(gmx_domdec_t*            dd,
                              int                      at_start,
                              const struct gmx_mtop_t* mtop,
                              const int*               cginfo,
                              gmx::Constraints*        constr,
                              int                      nrec,
                              t_ilist*                 il_local)
{
    gmx_domdec_constraints_t*     dc;
    t_ilist *                     ilc_local, *ils_local;
    std::vector<int>*             ireq;
    gmx::ArrayRef<const t_blocka> at2con_mt;
    gmx::HashedMap<int>*          ga2la_specat;
    int                           at_end, i, j;
    t_iatom*                      iap;

    // This code should not be called unless this condition is true,
    // because that's the only time init_domdec_constraints is
    // called...
    GMX_RELEASE_ASSERT(dd->comm->systemInfo.haveSplitConstraints || dd->comm->systemInfo.haveSplitSettles,
                       "dd_make_local_constraints called when there are no local constraints");
    // ... and init_domdec_constraints always sets
    // dd->constraint_comm...
    GMX_RELEASE_ASSERT(
            dd->constraint_comm,
            "Invalid use of dd_make_local_constraints before construction of constraint_comm");
    // ... which static analysis needs to be reassured about, because
    // otherwise, when dd->splitSettles is
    // true. dd->constraint_comm is unilaterally dereferenced before
    // the call to atoms_to_settles.

    dc = dd->constraints;

    ilc_local = &il_local[F_CONSTR];
    ils_local = &il_local[F_SETTLE];

    dc->ncon      = 0;
    ilc_local->nr = 0;
    if (dd->constraint_comm)
    {
        // TODO Perhaps gmx_domdec_constraints_t should keep a valid constr?
        GMX_RELEASE_ASSERT(constr != nullptr, "Must have valid constraints object");
        at2con_mt = constr->atom2constraints_moltype();
        ireq      = &dc->requestedGlobalAtomIndices[0];
        ireq->clear();
    }
    else
    {
        // Currently unreachable
        at2con_mt = {};
        ireq      = nullptr;
    }

    gmx::ArrayRef<const std::vector<int>> at2settle_mt;
    /* When settle works inside charge groups, we assigned them already */
    if (dd->comm->systemInfo.haveSplitSettles)
    {
        // TODO Perhaps gmx_domdec_constraints_t should keep a valid constr?
        GMX_RELEASE_ASSERT(constr != nullptr, "Must have valid constraints object");
        at2settle_mt  = constr->atom2settle_moltype();
        ils_local->nr = 0;
    }

    if (at2settle_mt.empty())
    {
        atoms_to_constraints(dd, mtop, cginfo, at2con_mt, nrec, ilc_local, ireq);
    }
    else
    {
        int t0_set;

        /* Do the constraints, if present, on the first thread.
         * Do the settles on all other threads.
         */
        t0_set = ((!at2con_mt.empty() && dc->nthread > 1) ? 1 : 0);

#pragma omp parallel for num_threads(dc->nthread) schedule(static)
        for (int thread = 0; thread < dc->nthread; thread++)
        {
            try
            {
                if (!at2con_mt.empty() && thread == 0)
                {
                    atoms_to_constraints(dd, mtop, cginfo, at2con_mt, nrec, ilc_local, ireq);
                }

                if (thread >= t0_set)
                {
                    int      cg0, cg1;
                    t_ilist* ilst;

                    /* Distribute the settle check+assignments over
                     * dc->nthread or dc->nthread-1 threads.
                     */
                    cg0 = (dd->ncg_home * (thread - t0_set)) / (dc->nthread - t0_set);
                    cg1 = (dd->ncg_home * (thread - t0_set + 1)) / (dc->nthread - t0_set);

                    if (thread == t0_set)
                    {
                        ilst = ils_local;
                    }
                    else
                    {
                        ilst = &dc->ils[thread];
                    }
                    ilst->nr = 0;

                    std::vector<int>& ireqt = dc->requestedGlobalAtomIndices[thread];
                    if (thread > 0)
                    {
                        ireqt.clear();
                    }

                    atoms_to_settles(dd, mtop, cginfo, at2settle_mt, cg0, cg1, ilst, &ireqt);
                }
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }

        /* Combine the generate settles and requested indices */
        for (int thread = 1; thread < dc->nthread; thread++)
        {
            t_ilist* ilst;
            int      ia;

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
                    ils_local->iatoms[ils_local->nr + ia] = ilst->iatoms[ia];
                }
                ils_local->nr += ilst->nr;
            }

            const std::vector<int>& ireqt = dc->requestedGlobalAtomIndices[thread];
            ireq->insert(ireq->end(), ireqt.begin(), ireqt.end());
        }

        if (debug)
        {
            fprintf(debug, "Settles: total %3d\n", ils_local->nr / 4);
        }
    }

    if (dd->constraint_comm)
    {
        int nral1;

        at_end = setup_specat_communication(dd, ireq, dd->constraint_comm, dd->constraints->ga2la.get(),
                                            at_start, 2, "constraint", " or lincs-order");

        /* Fill in the missing indices */
        ga2la_specat = dd->constraints->ga2la.get();

        nral1 = 1 + NRAL(F_CONSTR);
        for (i = 0; i < ilc_local->nr; i += nral1)
        {
            iap = ilc_local->iatoms + i;
            for (j = 1; j < nral1; j++)
            {
                if (iap[j] < 0)
                {
                    const int* a = ga2la_specat->find(-iap[j] - 1);
                    GMX_ASSERT(a, "We have checked before that this atom index has been set");
                    iap[j] = *a;
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
                    const int* a = ga2la_specat->find(-iap[j] - 1);
                    GMX_ASSERT(a, "We have checked before that this atom index has been set");
                    iap[j] = *a;
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

void init_domdec_constraints(gmx_domdec_t* dd, const gmx_mtop_t* mtop)
{
    gmx_domdec_constraints_t* dc;
    const gmx_molblock_t*     molb;

    if (debug)
    {
        fprintf(debug, "Begin init_domdec_constraints\n");
    }

    dd->constraints = new gmx_domdec_constraints_t;
    dc              = dd->constraints;

    dc->molb_con_offset.resize(mtop->molblock.size());
    dc->molb_ncon_mol.resize(mtop->molblock.size());

    int ncon = 0;
    for (size_t mb = 0; mb < mtop->molblock.size(); mb++)
    {
        molb                    = &mtop->molblock[mb];
        dc->molb_con_offset[mb] = ncon;
        dc->molb_ncon_mol[mb]   = mtop->moltype[molb->type].ilist[F_CONSTR].size() / 3
                                + mtop->moltype[molb->type].ilist[F_CONSTRNC].size() / 3;
        ncon += molb->nmol * dc->molb_ncon_mol[mb];
    }

    if (ncon > 0)
    {
        dc->gc_req.resize(ncon);
    }

    /* Use a hash table for the global to local index.
     * The number of keys is a rough estimate, it will be optimized later.
     */
    int numKeysEstimate = std::min(mtop->natoms / 20, mtop->natoms / (2 * dd->nnodes));
    dc->ga2la           = std::make_unique<gmx::HashedMap<int>>(numKeysEstimate);

    dc->nthread = gmx_omp_nthreads_get(emntDomdec);
    dc->ils.resize(dc->nthread);

    dd->constraint_comm = new gmx_domdec_specat_comm_t;

    dc->requestedGlobalAtomIndices.resize(dc->nthread);
}

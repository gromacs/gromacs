/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "vsite.h"

#include <stdio.h>

#include <algorithm>
#include <vector>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"


/* The strategy used here for assigning virtual sites to (thread-)tasks
 * is as follows:
 *
 * We divide the atom range that vsites operate on (natoms_local with DD,
 * 0 - last atom involved in vsites without DD) equally over all threads.
 *
 * Vsites in the local range constructed from atoms in the local range
 * and/or other vsites that are fully local are assigned to a simple,
 * independent task.
 *
 * Vsites that are not assigned after using the above criterion get assigned
 * to a so called "interdependent" thread task when none of the constructing
 * atoms is a vsite. These tasks are called interdependent, because one task
 * accesses atoms assigned to a different task/thread.
 * Note that this option is turned off with large (local) atom counts
 * to avoid high memory usage.
 *
 * Any remaining vsites are assigned to a separate master thread task.
 */

using gmx::RVec;

static void init_ilist(t_ilist *ilist)
{
    for (int i = 0; i < F_NRE; i++)
    {
        ilist[i].nr     = 0;
        ilist[i].nalloc = 0;
        ilist[i].iatoms = nullptr;
    }
}

/*! \brief List of atom indices belonging to a task */
struct AtomIndex {
    //! List of atom indices
    std::vector<int> atom;
};

/*! \brief Data structure for thread tasks that use constructing atoms outside their own atom range */
struct InterdependentTask
{
    //! The interaction lists, only vsite entries are used
    t_ilist                ilist[F_NRE];
    //! Thread/task-local force buffer
    std::vector<RVec>      force;
    //! The atom indices of the vsites of our task
    std::vector<int>       vsite;
    //! Flags if elements in force are spread to or not
    std::vector<bool>      use;
    //! The number of entries set to true in use
    int                    nuse;
    //! Array of atoms indices, size nthreads, covering all nuse set elements in use
    std::vector<AtomIndex> atomIndex;
    //! List of tasks (force blocks) this task spread forces to
    std::vector<int>       spreadTask;
    //! List of tasks that write to this tasks force block range
    std::vector<int>       reduceTask;

    InterdependentTask()
    {
        init_ilist(ilist);
        nuse = 0;
    }
};

/*! \brief Vsite thread task data structure */
struct VsiteThread {
    //! Start of atom range of this task
    int                rangeStart;
    //! End of atom range of this task
    int                rangeEnd;
    //! The interaction lists, only vsite entries are used
    t_ilist            ilist[F_NRE];
    //! Local fshift accumulation buffer
    rvec               fshift[SHIFTS];
    //! Local virial dx*df accumulation buffer
    matrix             dxdf;
    //! Tells if interdependent task idTask should be used (in addition to the rest of this task), this bool has the same value on all threads
    bool               useInterdependentTask;
    //! Data for vsites that involve constructing atoms in the atom range of other threads/tasks
    InterdependentTask idTask;

    /*! \brief Constructor */
    VsiteThread()
    {
        rangeStart            = -1;
        rangeEnd              = -1;
        init_ilist(ilist);
        clear_rvecs(SHIFTS, fshift);
        clear_mat(dxdf);
        useInterdependentTask = false;
    }
};


/* The start and end values of for the vsite indices in the ftype enum.
 * The validity of these values is checked in init_vsite.
 * This is used to avoid loops over all ftypes just to get the vsite entries.
 * (We should replace the fixed ilist array by only the used entries.)
 */
static const int c_ftypeVsiteStart = F_VSITE2;
static const int c_ftypeVsiteEnd   = F_VSITEN + 1;


/*! \brief Returns the sum of the vsite ilist sizes over all vsite types
 *
 * \param[in] ilist  The interaction list
 */
static int vsiteIlistNrCount(const t_ilist *ilist)
{
    int nr = 0;
    for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
    {
        nr += ilist[ftype].nr;
    }

    return nr;
}

static int pbc_rvec_sub(const t_pbc *pbc, const rvec xi, const rvec xj, rvec dx)
{
    if (pbc)
    {
        return pbc_dx_aiuc(pbc, xi, xj, dx);
    }
    else
    {
        rvec_sub(xi, xj, dx);
        return CENTRAL;
    }
}

/* Vsite construction routines */

static void constr_vsite2(const rvec xi, const rvec xj, rvec x, real a, const t_pbc *pbc)
{
    real b = 1 - a;
    /* 1 flop */

    if (pbc)
    {
        rvec dx;
        pbc_dx_aiuc(pbc, xj, xi, dx);
        x[XX] = xi[XX] + a*dx[XX];
        x[YY] = xi[YY] + a*dx[YY];
        x[ZZ] = xi[ZZ] + a*dx[ZZ];
    }
    else
    {
        x[XX] = b*xi[XX] + a*xj[XX];
        x[YY] = b*xi[YY] + a*xj[YY];
        x[ZZ] = b*xi[ZZ] + a*xj[ZZ];
        /* 9 Flops */
    }

    /* TOTAL: 10 flops */
}

static void constr_vsite3(const rvec xi, const rvec xj, const rvec xk, rvec x, real a, real b,
                          const t_pbc *pbc)
{
    real c = 1 - a - b;
    /* 2 flops */

    if (pbc)
    {
        rvec dxj, dxk;

        pbc_dx_aiuc(pbc, xj, xi, dxj);
        pbc_dx_aiuc(pbc, xk, xi, dxk);
        x[XX] = xi[XX] + a*dxj[XX] + b*dxk[XX];
        x[YY] = xi[YY] + a*dxj[YY] + b*dxk[YY];
        x[ZZ] = xi[ZZ] + a*dxj[ZZ] + b*dxk[ZZ];
    }
    else
    {
        x[XX] = c*xi[XX] + a*xj[XX] + b*xk[XX];
        x[YY] = c*xi[YY] + a*xj[YY] + b*xk[YY];
        x[ZZ] = c*xi[ZZ] + a*xj[ZZ] + b*xk[ZZ];
        /* 15 Flops */
    }

    /* TOTAL: 17 flops */
}

static void constr_vsite3FD(const rvec xi, const rvec xj, const rvec xk, rvec x, real a, real b,
                            const t_pbc *pbc)
{
    rvec xij, xjk, temp;
    real c;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xj, xjk);
    /* 6 flops */

    /* temp goes from i to a point on the line jk */
    temp[XX] = xij[XX] + a*xjk[XX];
    temp[YY] = xij[YY] + a*xjk[YY];
    temp[ZZ] = xij[ZZ] + a*xjk[ZZ];
    /* 6 flops */

    c = b*gmx::invsqrt(iprod(temp, temp));
    /* 6 + 10 flops */

    x[XX] = xi[XX] + c*temp[XX];
    x[YY] = xi[YY] + c*temp[YY];
    x[ZZ] = xi[ZZ] + c*temp[ZZ];
    /* 6 Flops */

    /* TOTAL: 34 flops */
}

static void constr_vsite3FAD(const rvec xi, const rvec xj, const rvec xk, rvec x, real a, real b, const t_pbc *pbc)
{
    rvec xij, xjk, xp;
    real a1, b1, c1, invdij;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xj, xjk);
    /* 6 flops */

    invdij = gmx::invsqrt(iprod(xij, xij));
    c1     = invdij * invdij * iprod(xij, xjk);
    xp[XX] = xjk[XX] - c1*xij[XX];
    xp[YY] = xjk[YY] - c1*xij[YY];
    xp[ZZ] = xjk[ZZ] - c1*xij[ZZ];
    a1     = a*invdij;
    b1     = b*gmx::invsqrt(iprod(xp, xp));
    /* 45 */

    x[XX] = xi[XX] + a1*xij[XX] + b1*xp[XX];
    x[YY] = xi[YY] + a1*xij[YY] + b1*xp[YY];
    x[ZZ] = xi[ZZ] + a1*xij[ZZ] + b1*xp[ZZ];
    /* 12 Flops */

    /* TOTAL: 63 flops */
}

static void constr_vsite3OUT(const rvec xi, const rvec xj, const rvec xk, rvec x,
                             real a, real b, real c, const t_pbc *pbc)
{
    rvec xij, xik, temp;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xi, xik);
    cprod(xij, xik, temp);
    /* 15 Flops */

    x[XX] = xi[XX] + a*xij[XX] + b*xik[XX] + c*temp[XX];
    x[YY] = xi[YY] + a*xij[YY] + b*xik[YY] + c*temp[YY];
    x[ZZ] = xi[ZZ] + a*xij[ZZ] + b*xik[ZZ] + c*temp[ZZ];
    /* 18 Flops */

    /* TOTAL: 33 flops */
}

static void constr_vsite4FD(const rvec xi, const rvec xj, const rvec xk, const rvec xl, rvec x,
                            real a, real b, real c, const t_pbc *pbc)
{
    rvec xij, xjk, xjl, temp;
    real d;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xj, xjk);
    pbc_rvec_sub(pbc, xl, xj, xjl);
    /* 9 flops */

    /* temp goes from i to a point on the plane jkl */
    temp[XX] = xij[XX] + a*xjk[XX] + b*xjl[XX];
    temp[YY] = xij[YY] + a*xjk[YY] + b*xjl[YY];
    temp[ZZ] = xij[ZZ] + a*xjk[ZZ] + b*xjl[ZZ];
    /* 12 flops */

    d = c*gmx::invsqrt(iprod(temp, temp));
    /* 6 + 10 flops */

    x[XX] = xi[XX] + d*temp[XX];
    x[YY] = xi[YY] + d*temp[YY];
    x[ZZ] = xi[ZZ] + d*temp[ZZ];
    /* 6 Flops */

    /* TOTAL: 43 flops */
}

static void constr_vsite4FDN(const rvec xi, const rvec xj, const rvec xk, const rvec xl, rvec x,
                             real a, real b, real c, const t_pbc *pbc)
{
    rvec xij, xik, xil, ra, rb, rja, rjb, rm;
    real d;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xi, xik);
    pbc_rvec_sub(pbc, xl, xi, xil);
    /* 9 flops */

    ra[XX] = a*xik[XX];
    ra[YY] = a*xik[YY];
    ra[ZZ] = a*xik[ZZ];

    rb[XX] = b*xil[XX];
    rb[YY] = b*xil[YY];
    rb[ZZ] = b*xil[ZZ];

    /* 6 flops */

    rvec_sub(ra, xij, rja);
    rvec_sub(rb, xij, rjb);
    /* 6 flops */

    cprod(rja, rjb, rm);
    /* 9 flops */

    d = c*gmx::invsqrt(norm2(rm));
    /* 5+5+1 flops */

    x[XX] = xi[XX] + d*rm[XX];
    x[YY] = xi[YY] + d*rm[YY];
    x[ZZ] = xi[ZZ] + d*rm[ZZ];
    /* 6 Flops */

    /* TOTAL: 47 flops */
}


static int constr_vsiten(const t_iatom *ia, const t_iparams ip[],
                         rvec *x, const t_pbc *pbc)
{
    rvec x1, dx;
    dvec dsum;
    int  n3, av, ai;
    real a;

    n3 = 3*ip[ia[0]].vsiten.n;
    av = ia[1];
    ai = ia[2];
    copy_rvec(x[ai], x1);
    clear_dvec(dsum);
    for (int i = 3; i < n3; i += 3)
    {
        ai = ia[i+2];
        a  = ip[ia[i]].vsiten.a;
        if (pbc)
        {
            pbc_dx_aiuc(pbc, x[ai], x1, dx);
        }
        else
        {
            rvec_sub(x[ai], x1, dx);
        }
        dsum[XX] += a*dx[XX];
        dsum[YY] += a*dx[YY];
        dsum[ZZ] += a*dx[ZZ];
        /* 9 Flops */
    }

    x[av][XX] = x1[XX] + dsum[XX];
    x[av][YY] = x1[YY] + dsum[YY];
    x[av][ZZ] = x1[ZZ] + dsum[ZZ];

    return n3;
}

/*! \brief PBC modes for vsite construction and spreading */
enum class PbcMode
{
    all,         // Apply normal, simple PBC for all vsites
    chargeGroup, // Keep vsite in the same periodic image as the rest of it's charge group
    none         // No PBC treatment needed
};

/*! \brief Returns the PBC mode based on the system PBC and vsite properties
 *
 * \param[in] pbcPtr  A pointer to a PBC struct or nullptr when no PBC treatment is required
 * \param[in] vsite   A pointer to the vsite struct, can be nullptr
 */
static PbcMode getPbcMode(const t_pbc       *pbcPtr,
                          const gmx_vsite_t *vsite)
{
    if (pbcPtr == nullptr)
    {
        return PbcMode::none;
    }
    else if (vsite != nullptr && vsite->bHaveChargeGroups)
    {
        return PbcMode::chargeGroup;
    }
    else
    {
        return PbcMode::all;
    }
}

static void construct_vsites_thread(const gmx_vsite_t *vsite,
                                    rvec x[],
                                    real dt, rvec *v,
                                    const t_iparams ip[], const t_ilist ilist[],
                                    const t_pbc *pbc_null)
{
    real         inv_dt;
    if (v != nullptr)
    {
        inv_dt = 1.0/dt;
    }
    else
    {
        inv_dt = 1.0;
    }

    const PbcMode  pbcMode   = getPbcMode(pbc_null, vsite);
    /* We need another pbc pointer, as with charge groups we switch per vsite */
    const t_pbc   *pbc_null2 = pbc_null;
    const int     *vsite_pbc = nullptr;

    for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
    {
        if (ilist[ftype].nr == 0)
        {
            continue;
        }

        {   // TODO remove me
            int            nra = interaction_function[ftype].nratoms;
            int            inc = 1 + nra;
            int            nr  = ilist[ftype].nr;

            const t_iatom *ia = ilist[ftype].iatoms;

            if (pbcMode == PbcMode::chargeGroup)
            {
                vsite_pbc = vsite->vsite_pbc_loc[ftype - c_ftypeVsiteStart];
            }

            for (int i = 0; i < nr; )
            {
                int  tp     = ia[0];
                /* The vsite and constructing atoms */
                int  avsite = ia[1];
                int  ai     = ia[2];
                /* Constants for constructing vsites */
                real a1     = ip[tp].vsite.a;
                /* Check what kind of pbc we need to use */
                int  pbc_atom;
                rvec xpbc;
                if (pbcMode == PbcMode::all)
                {
                    /* No charge groups, vsite follows its own pbc */
                    pbc_atom = avsite;
                    copy_rvec(x[avsite], xpbc);
                }
                else if (pbcMode == PbcMode::chargeGroup)
                {
                    pbc_atom = vsite_pbc[i/(1 + nra)];
                    if (pbc_atom > -2)
                    {
                        if (pbc_atom >= 0)
                        {
                            /* We need to copy the coordinates here,
                             * single for single atom cg's pbc_atom
                             * is the vsite itself.
                             */
                            copy_rvec(x[pbc_atom], xpbc);
                        }
                        pbc_null2 = pbc_null;
                    }
                    else
                    {
                        pbc_null2 = nullptr;
                    }
                }
                else
                {
                    pbc_atom = -2;
                }
                /* Copy the old position */
                rvec xv;
                copy_rvec(x[avsite], xv);

                /* Construct the vsite depending on type */
                int  aj, ak, al;
                real b1, c1;
                switch (ftype)
                {
                    case F_VSITE2:
                        aj = ia[3];
                        constr_vsite2(x[ai], x[aj], x[avsite], a1, pbc_null2);
                        break;
                    case F_VSITE3:
                        aj = ia[3];
                        ak = ia[4];
                        b1 = ip[tp].vsite.b;
                        constr_vsite3(x[ai], x[aj], x[ak], x[avsite], a1, b1, pbc_null2);
                        break;
                    case F_VSITE3FD:
                        aj = ia[3];
                        ak = ia[4];
                        b1 = ip[tp].vsite.b;
                        constr_vsite3FD(x[ai], x[aj], x[ak], x[avsite], a1, b1, pbc_null2);
                        break;
                    case F_VSITE3FAD:
                        aj = ia[3];
                        ak = ia[4];
                        b1 = ip[tp].vsite.b;
                        constr_vsite3FAD(x[ai], x[aj], x[ak], x[avsite], a1, b1, pbc_null2);
                        break;
                    case F_VSITE3OUT:
                        aj = ia[3];
                        ak = ia[4];
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        constr_vsite3OUT(x[ai], x[aj], x[ak], x[avsite], a1, b1, c1, pbc_null2);
                        break;
                    case F_VSITE4FD:
                        aj = ia[3];
                        ak = ia[4];
                        al = ia[5];
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        constr_vsite4FD(x[ai], x[aj], x[ak], x[al], x[avsite], a1, b1, c1,
                                        pbc_null2);
                        break;
                    case F_VSITE4FDN:
                        aj = ia[3];
                        ak = ia[4];
                        al = ia[5];
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        constr_vsite4FDN(x[ai], x[aj], x[ak], x[al], x[avsite], a1, b1, c1,
                                         pbc_null2);
                        break;
                    case F_VSITEN:
                        inc = constr_vsiten(ia, ip, x, pbc_null2);
                        break;
                    default:
                        gmx_fatal(FARGS, "No such vsite type %d in %s, line %d",
                                  ftype, __FILE__, __LINE__);
                }

                if (pbc_atom >= 0)
                {
                    /* Match the pbc of this vsite to the rest of its charge group */
                    rvec dx;
                    int  ishift = pbc_dx_aiuc(pbc_null, x[avsite], xpbc, dx);
                    if (ishift != CENTRAL)
                    {
                        rvec_add(xpbc, dx, x[avsite]);
                    }
                }
                if (v != nullptr)
                {
                    /* Calculate velocity of vsite... */
                    rvec vv;
                    rvec_sub(x[avsite], xv, vv);
                    svmul(inv_dt, vv, v[avsite]);
                }

                /* Increment loop variables */
                i  += inc;
                ia += inc;
            }
        }
    }
}

void construct_vsites(const gmx_vsite_t *vsite,
                      rvec x[],
                      real dt, rvec *v,
                      const t_iparams ip[], const t_ilist ilist[],
                      int ePBC, gmx_bool bMolPBC,
                      t_commrec *cr,
                      const matrix box)
{
    const bool useDomdec = (vsite != nullptr && vsite->useDomdec);
    GMX_ASSERT(!useDomdec || (cr != nullptr && DOMAINDECOMP(cr)), "When vsites are set up with domain decomposition, we need a valid commrec");
    // TODO: Remove this assertion when we remove charge groups
    GMX_ASSERT(vsite != nullptr || ePBC == epbcNONE, "Without a vsite struct we can not do PBC (in case we have charge groups)");

    t_pbc     pbc, *pbc_null;

    /* We only need to do pbc when we have inter-cg vsites.
     * Note that with domain decomposition we do not need to apply PBC here
     * when we have at least 3 domains along each dimension. Currently we
     * do not optimize this case.
     */
    if (ePBC != epbcNONE && (useDomdec || bMolPBC) &&
        !(vsite != nullptr && vsite->n_intercg_vsite == 0))
    {
        /* This is wasting some CPU time as we now do this multiple times
         * per MD step.
         */
        ivec null_ivec;
        clear_ivec(null_ivec);
        pbc_null = set_pbc_dd(&pbc, ePBC,
                              useDomdec ? cr->dd->nc : null_ivec,
                              FALSE, box);
    }
    else
    {
        pbc_null = nullptr;
    }

    if (useDomdec)
    {
        dd_move_x_vsites(cr->dd, box, x);
    }

    // cppcheck-suppress nullPointerRedundantCheck
    if (vsite == nullptr || vsite->nthreads == 1)
    {
        construct_vsites_thread(vsite,
                                x, dt, v,
                                ip, ilist,
                                pbc_null);
    }
    else
    {
#pragma omp parallel num_threads(vsite->nthreads)
        {
            try
            {
                const int          th    = gmx_omp_get_thread_num();
                const VsiteThread &tData = *vsite->tData[th];
                GMX_ASSERT(tData.rangeStart >= 0, "The thread data should be initialized before calling construct_vsites");

                construct_vsites_thread(vsite,
                                        x, dt, v,
                                        ip, tData.ilist,
                                        pbc_null);
                if (tData.useInterdependentTask)
                {
                    /* Here we don't need a barrier (unlike the spreading),
                     * since both tasks only construct vsites from particles,
                     * or local vsites, not from non-local vsites.
                     */
                    construct_vsites_thread(vsite,
                                            x, dt, v,
                                            ip, tData.idTask.ilist,
                                            pbc_null);
                }
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        }
        /* Now we can construct the vsites that might depend on other vsites */
        construct_vsites_thread(vsite,
                                x, dt, v,
                                ip, vsite->tData[vsite->nthreads]->ilist,
                                pbc_null);
    }
}

static void spread_vsite2(const t_iatom ia[], real a,
                          const rvec x[],
                          rvec f[], rvec fshift[],
                          const t_pbc *pbc, const t_graph *g)
{
    rvec    fi, fj, dx;
    t_iatom av, ai, aj;
    ivec    di;
    int     siv, sij;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];

    svmul(1 - a, f[av], fi);
    svmul(    a, f[av], fj);
    /* 7 flop */

    rvec_inc(f[ai], fi);
    rvec_inc(f[aj], fj);
    /* 6 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, av), di);
        siv = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), di);
        sij = IVEC2IS(di);
    }
    else if (pbc)
    {
        siv = pbc_dx_aiuc(pbc, x[ai], x[av], dx);
        sij = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
    }
    else
    {
        siv = CENTRAL;
        sij = CENTRAL;
    }

    if (fshift && (siv != CENTRAL || sij != CENTRAL))
    {
        rvec_inc(fshift[siv], f[av]);
        rvec_dec(fshift[CENTRAL], fi);
        rvec_dec(fshift[sij], fj);
    }

    /* TOTAL: 13 flops */
}

void constructVsitesGlobal(const gmx_mtop_t         &mtop,
                           gmx::ArrayRef<gmx::RVec>  x)
{
    GMX_ASSERT(x.size() >= static_cast<size_t>(mtop.natoms), "x should contain the whole system");

    for (int mb = 0; mb < mtop.nmolblock; mb++)
    {
        const gmx_molblock_t &molb = mtop.molblock[mb];
        const gmx_moltype_t  &molt = mtop.moltype[molb.type];
        if (vsiteIlistNrCount(molt.ilist) > 0)
        {
            int atomOffset = molb.globalAtomStart;
            for (int mol = 0; mol < molb.nmol; mol++)
            {
                construct_vsites(nullptr, as_rvec_array(x.data()) + atomOffset,
                                 0.0, nullptr,
                                 mtop.ffparams.iparams, molt.ilist,
                                 epbcNONE, TRUE, nullptr, nullptr);
                atomOffset += molt.atoms.nr;
            }
        }
    }
}

static void spread_vsite3(const t_iatom ia[], real a, real b,
                          const rvec x[],
                          rvec f[], rvec fshift[],
                          const t_pbc *pbc, const t_graph *g)
{
    rvec    fi, fj, fk, dx;
    int     av, ai, aj, ak;
    ivec    di;
    int     siv, sij, sik;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];

    svmul(1 - a - b, f[av], fi);
    svmul(        a, f[av], fj);
    svmul(        b, f[av], fk);
    /* 11 flops */

    rvec_inc(f[ai], fi);
    rvec_inc(f[aj], fj);
    rvec_inc(f[ak], fk);
    /* 9 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, ia[1]), di);
        siv = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), di);
        sij = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, ak), di);
        sik = IVEC2IS(di);
    }
    else if (pbc)
    {
        siv = pbc_dx_aiuc(pbc, x[ai], x[av], dx);
        sij = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        sik = pbc_dx_aiuc(pbc, x[ai], x[ak], dx);
    }
    else
    {
        siv = CENTRAL;
        sij = CENTRAL;
        sik = CENTRAL;
    }

    if (fshift && (siv != CENTRAL || sij != CENTRAL || sik != CENTRAL))
    {
        rvec_inc(fshift[siv], f[av]);
        rvec_dec(fshift[CENTRAL], fi);
        rvec_dec(fshift[sij], fj);
        rvec_dec(fshift[sik], fk);
    }

    /* TOTAL: 20 flops */
}

static void spread_vsite3FD(const t_iatom ia[], real a, real b,
                            const rvec x[],
                            rvec f[], rvec fshift[],
                            gmx_bool VirCorr, matrix dxdf,
                            const t_pbc *pbc, const t_graph *g)
{
    real    c, invl, fproj, a1;
    rvec    xvi, xij, xjk, xix, fv, temp;
    t_iatom av, ai, aj, ak;
    int     svi, sji, skj;
    ivec    di;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];
    copy_rvec(f[av], fv);

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    skj = pbc_rvec_sub(pbc, x[ak], x[aj], xjk);
    /* 6 flops */

    /* xix goes from i to point x on the line jk */
    xix[XX] = xij[XX]+a*xjk[XX];
    xix[YY] = xij[YY]+a*xjk[YY];
    xix[ZZ] = xij[ZZ]+a*xjk[ZZ];
    /* 6 flops */

    invl = gmx::invsqrt(iprod(xix, xix));
    c    = b*invl;
    /* 4 + ?10? flops */

    fproj = iprod(xix, fv)*invl*invl; /* = (xix . f)/(xix . xix) */

    temp[XX] = c*(fv[XX]-fproj*xix[XX]);
    temp[YY] = c*(fv[YY]-fproj*xix[YY]);
    temp[ZZ] = c*(fv[ZZ]-fproj*xix[ZZ]);
    /* 16 */

    /* c is already calculated in constr_vsite3FD
       storing c somewhere will save 26 flops!     */

    a1         = 1 - a;
    f[ai][XX] += fv[XX] - temp[XX];
    f[ai][YY] += fv[YY] - temp[YY];
    f[ai][ZZ] += fv[ZZ] - temp[ZZ];
    f[aj][XX] += a1*temp[XX];
    f[aj][YY] += a1*temp[YY];
    f[aj][ZZ] += a1*temp[ZZ];
    f[ak][XX] += a*temp[XX];
    f[ak][YY] += a*temp[YY];
    f[ak][ZZ] += a*temp[ZZ];
    /* 19 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ia[1]), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sji = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, aj), di);
        skj = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift && (svi != CENTRAL || sji != CENTRAL || skj != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        fshift[CENTRAL][XX] += fv[XX] - (1 + a)*temp[XX];
        fshift[CENTRAL][YY] += fv[YY] - (1 + a)*temp[YY];
        fshift[CENTRAL][ZZ] += fv[ZZ] - (1 + a)*temp[ZZ];
        fshift[    sji][XX] += temp[XX];
        fshift[    sji][YY] += temp[YY];
        fshift[    sji][ZZ] += temp[ZZ];
        fshift[    skj][XX] += a*temp[XX];
        fshift[    skj][YY] += a*temp[YY];
        fshift[    skj][ZZ] += a*temp[ZZ];
    }

    if (VirCorr)
    {
        /* When VirCorr=TRUE, the virial for the current forces is not
         * calculated from the redistributed forces. This means that
         * the effect of non-linear virtual site constructions on the virial
         * needs to be added separately. This contribution can be calculated
         * in many ways, but the simplest and cheapest way is to use
         * the first constructing atom ai as a reference position in space:
         * subtract (xv-xi)*fv and add (xj-xi)*fj + (xk-xi)*fk.
         */
        rvec xiv;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                /* As xix is a linear combination of j and k, use that here */
                dxdf[i][j] += -xiv[i]*fv[j] + xix[i]*temp[j];
            }
        }
    }

    /* TOTAL: 61 flops */
}

static void spread_vsite3FAD(const t_iatom ia[], real a, real b,
                             const rvec x[],
                             rvec f[], rvec fshift[],
                             gmx_bool VirCorr, matrix dxdf,
                             const t_pbc *pbc, const t_graph *g)
{
    rvec    xvi, xij, xjk, xperp, Fpij, Fppp, fv, f1, f2, f3;
    real    a1, b1, c1, c2, invdij, invdij2, invdp, fproj;
    t_iatom av, ai, aj, ak;
    int     svi, sji, skj, d;
    ivec    di;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];
    copy_rvec(f[ia[1]], fv);

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    skj = pbc_rvec_sub(pbc, x[ak], x[aj], xjk);
    /* 6 flops */

    invdij    = gmx::invsqrt(iprod(xij, xij));
    invdij2   = invdij * invdij;
    c1        = iprod(xij, xjk) * invdij2;
    xperp[XX] = xjk[XX] - c1*xij[XX];
    xperp[YY] = xjk[YY] - c1*xij[YY];
    xperp[ZZ] = xjk[ZZ] - c1*xij[ZZ];
    /* xperp in plane ijk, perp. to ij */
    invdp = gmx::invsqrt(iprod(xperp, xperp));
    a1    = a*invdij;
    b1    = b*invdp;
    /* 45 flops */

    /* a1, b1 and c1 are already calculated in constr_vsite3FAD
       storing them somewhere will save 45 flops!     */

    fproj = iprod(xij, fv)*invdij2;
    svmul(fproj,                      xij,  Fpij);    /* proj. f on xij */
    svmul(iprod(xperp, fv)*invdp*invdp, xperp, Fppp); /* proj. f on xperp */
    svmul(b1*fproj,                   xperp, f3);
    /* 23 flops */

    rvec_sub(fv, Fpij, f1); /* f1 = f - Fpij */
    rvec_sub(f1, Fppp, f2); /* f2 = f - Fpij - Fppp */
    for (d = 0; (d < DIM); d++)
    {
        f1[d] *= a1;
        f2[d] *= b1;
    }
    /* 12 flops */

    c2         = 1 + c1;
    f[ai][XX] += fv[XX] - f1[XX] + c1*f2[XX] + f3[XX];
    f[ai][YY] += fv[YY] - f1[YY] + c1*f2[YY] + f3[YY];
    f[ai][ZZ] += fv[ZZ] - f1[ZZ] + c1*f2[ZZ] + f3[ZZ];
    f[aj][XX] +=          f1[XX] - c2*f2[XX] - f3[XX];
    f[aj][YY] +=          f1[YY] - c2*f2[YY] - f3[YY];
    f[aj][ZZ] +=          f1[ZZ] - c2*f2[ZZ] - f3[ZZ];
    f[ak][XX] +=                      f2[XX];
    f[ak][YY] +=                      f2[YY];
    f[ak][ZZ] +=                      f2[ZZ];
    /* 30 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ia[1]), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sji = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, aj), di);
        skj = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift && (svi != CENTRAL || sji != CENTRAL || skj != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        fshift[CENTRAL][XX] += fv[XX] - f1[XX] - (1-c1)*f2[XX] + f3[XX];
        fshift[CENTRAL][YY] += fv[YY] - f1[YY] - (1-c1)*f2[YY] + f3[YY];
        fshift[CENTRAL][ZZ] += fv[ZZ] - f1[ZZ] - (1-c1)*f2[ZZ] + f3[ZZ];
        fshift[    sji][XX] +=          f1[XX] -    c1 *f2[XX] - f3[XX];
        fshift[    sji][YY] +=          f1[YY] -    c1 *f2[YY] - f3[YY];
        fshift[    sji][ZZ] +=          f1[ZZ] -    c1 *f2[ZZ] - f3[ZZ];
        fshift[    skj][XX] +=                          f2[XX];
        fshift[    skj][YY] +=                          f2[YY];
        fshift[    skj][ZZ] +=                          f2[ZZ];
    }

    if (VirCorr)
    {
        rvec xiv;
        int  i, j;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                /* Note that xik=xij+xjk, so we have to add xij*f2 */
                dxdf[i][j] +=
                    -xiv[i]*fv[j]
                    + xij[i]*(f1[j] + (1 - c2)*f2[j] - f3[j])
                    + xjk[i]*f2[j];
            }
        }
    }

    /* TOTAL: 113 flops */
}

static void spread_vsite3OUT(const t_iatom ia[], real a, real b, real c,
                             const rvec x[],
                             rvec f[], rvec fshift[],
                             gmx_bool VirCorr, matrix dxdf,
                             const t_pbc *pbc, const t_graph *g)
{
    rvec    xvi, xij, xik, fv, fj, fk;
    real    cfx, cfy, cfz;
    int     av, ai, aj, ak;
    ivec    di;
    int     svi, sji, ski;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    ski = pbc_rvec_sub(pbc, x[ak], x[ai], xik);
    /* 6 Flops */

    copy_rvec(f[av], fv);

    cfx = c*fv[XX];
    cfy = c*fv[YY];
    cfz = c*fv[ZZ];
    /* 3 Flops */

    fj[XX] = a*fv[XX]     -  xik[ZZ]*cfy +  xik[YY]*cfz;
    fj[YY] =  xik[ZZ]*cfx + a*fv[YY]     -  xik[XX]*cfz;
    fj[ZZ] = -xik[YY]*cfx +  xik[XX]*cfy + a*fv[ZZ];

    fk[XX] = b*fv[XX]     +  xij[ZZ]*cfy -  xij[YY]*cfz;
    fk[YY] = -xij[ZZ]*cfx + b*fv[YY]     +  xij[XX]*cfz;
    fk[ZZ] =  xij[YY]*cfx -  xij[XX]*cfy + b*fv[ZZ];
    /* 30 Flops */

    f[ai][XX] += fv[XX] - fj[XX] - fk[XX];
    f[ai][YY] += fv[YY] - fj[YY] - fk[YY];
    f[ai][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ];
    rvec_inc(f[aj], fj);
    rvec_inc(f[ak], fk);
    /* 15 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ia[1]), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sji = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, ai), di);
        ski = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift && (svi != CENTRAL || sji != CENTRAL || ski != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        fshift[CENTRAL][XX] += fv[XX] - fj[XX] - fk[XX];
        fshift[CENTRAL][YY] += fv[YY] - fj[YY] - fk[YY];
        fshift[CENTRAL][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ];
        rvec_inc(fshift[sji], fj);
        rvec_inc(fshift[ski], fk);
    }

    if (VirCorr)
    {
        rvec xiv;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                dxdf[i][j] += -xiv[i]*fv[j] + xij[i]*fj[j] + xik[i]*fk[j];
            }
        }
    }

    /* TOTAL: 54 flops */
}

static void spread_vsite4FD(const t_iatom ia[], real a, real b, real c,
                            const rvec x[],
                            rvec f[], rvec fshift[],
                            gmx_bool VirCorr, matrix dxdf,
                            const t_pbc *pbc, const t_graph *g)
{
    real    d, invl, fproj, a1;
    rvec    xvi, xij, xjk, xjl, xix, fv, temp;
    int     av, ai, aj, ak, al;
    ivec    di;
    int     svi, sji, skj, slj, m;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];
    al = ia[5];

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    skj = pbc_rvec_sub(pbc, x[ak], x[aj], xjk);
    slj = pbc_rvec_sub(pbc, x[al], x[aj], xjl);
    /* 9 flops */

    /* xix goes from i to point x on the plane jkl */
    for (m = 0; m < DIM; m++)
    {
        xix[m] = xij[m] + a*xjk[m] + b*xjl[m];
    }
    /* 12 flops */

    invl = gmx::invsqrt(iprod(xix, xix));
    d    = c*invl;
    /* 4 + ?10? flops */

    copy_rvec(f[av], fv);

    fproj = iprod(xix, fv)*invl*invl; /* = (xix . f)/(xix . xix) */

    for (m = 0; m < DIM; m++)
    {
        temp[m] = d*(fv[m] - fproj*xix[m]);
    }
    /* 16 */

    /* c is already calculated in constr_vsite3FD
       storing c somewhere will save 35 flops!     */

    a1 = 1 - a - b;
    for (m = 0; m < DIM; m++)
    {
        f[ai][m] += fv[m] - temp[m];
        f[aj][m] += a1*temp[m];
        f[ak][m] += a*temp[m];
        f[al][m] += b*temp[m];
    }
    /* 26 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ia[1]), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sji = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, aj), di);
        skj = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, al), SHIFT_IVEC(g, aj), di);
        slj = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift &&
        (svi != CENTRAL || sji != CENTRAL || skj != CENTRAL || slj != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        for (m = 0; m < DIM; m++)
        {
            fshift[CENTRAL][m] += fv[m] - (1 + a + b)*temp[m];
            fshift[    sji][m] += temp[m];
            fshift[    skj][m] += a*temp[m];
            fshift[    slj][m] += b*temp[m];
        }
    }

    if (VirCorr)
    {
        rvec xiv;
        int  i, j;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                dxdf[i][j] += -xiv[i]*fv[j] + xix[i]*temp[j];
            }
        }
    }

    /* TOTAL: 77 flops */
}


static void spread_vsite4FDN(const t_iatom ia[], real a, real b, real c,
                             const rvec x[],
                             rvec f[], rvec fshift[],
                             gmx_bool VirCorr, matrix dxdf,
                             const t_pbc *pbc, const t_graph *g)
{
    rvec xvi, xij, xik, xil, ra, rb, rja, rjb, rab, rm, rt;
    rvec fv, fj, fk, fl;
    real invrm, denom;
    real cfx, cfy, cfz;
    ivec di;
    int  av, ai, aj, ak, al;
    int  svi, sij, sik, sil;

    /* DEBUG: check atom indices */
    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];
    al = ia[5];

    copy_rvec(f[av], fv);

    sij = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    sik = pbc_rvec_sub(pbc, x[ak], x[ai], xik);
    sil = pbc_rvec_sub(pbc, x[al], x[ai], xil);
    /* 9 flops */

    ra[XX] = a*xik[XX];
    ra[YY] = a*xik[YY];
    ra[ZZ] = a*xik[ZZ];

    rb[XX] = b*xil[XX];
    rb[YY] = b*xil[YY];
    rb[ZZ] = b*xil[ZZ];

    /* 6 flops */

    rvec_sub(ra, xij, rja);
    rvec_sub(rb, xij, rjb);
    rvec_sub(rb, ra, rab);
    /* 9 flops */

    cprod(rja, rjb, rm);
    /* 9 flops */

    invrm = gmx::invsqrt(norm2(rm));
    denom = invrm*invrm;
    /* 5+5+2 flops */

    cfx = c*invrm*fv[XX];
    cfy = c*invrm*fv[YY];
    cfz = c*invrm*fv[ZZ];
    /* 6 Flops */

    cprod(rm, rab, rt);
    /* 9 flops */

    rt[XX] *= denom;
    rt[YY] *= denom;
    rt[ZZ] *= denom;
    /* 3flops */

    fj[XX] = (        -rm[XX]*rt[XX]) * cfx + ( rab[ZZ]-rm[YY]*rt[XX]) * cfy + (-rab[YY]-rm[ZZ]*rt[XX]) * cfz;
    fj[YY] = (-rab[ZZ]-rm[XX]*rt[YY]) * cfx + (        -rm[YY]*rt[YY]) * cfy + ( rab[XX]-rm[ZZ]*rt[YY]) * cfz;
    fj[ZZ] = ( rab[YY]-rm[XX]*rt[ZZ]) * cfx + (-rab[XX]-rm[YY]*rt[ZZ]) * cfy + (        -rm[ZZ]*rt[ZZ]) * cfz;
    /* 30 flops */

    cprod(rjb, rm, rt);
    /* 9 flops */

    rt[XX] *= denom*a;
    rt[YY] *= denom*a;
    rt[ZZ] *= denom*a;
    /* 3flops */

    fk[XX] = (          -rm[XX]*rt[XX]) * cfx + (-a*rjb[ZZ]-rm[YY]*rt[XX]) * cfy + ( a*rjb[YY]-rm[ZZ]*rt[XX]) * cfz;
    fk[YY] = ( a*rjb[ZZ]-rm[XX]*rt[YY]) * cfx + (          -rm[YY]*rt[YY]) * cfy + (-a*rjb[XX]-rm[ZZ]*rt[YY]) * cfz;
    fk[ZZ] = (-a*rjb[YY]-rm[XX]*rt[ZZ]) * cfx + ( a*rjb[XX]-rm[YY]*rt[ZZ]) * cfy + (          -rm[ZZ]*rt[ZZ]) * cfz;
    /* 36 flops */

    cprod(rm, rja, rt);
    /* 9 flops */

    rt[XX] *= denom*b;
    rt[YY] *= denom*b;
    rt[ZZ] *= denom*b;
    /* 3flops */

    fl[XX] = (          -rm[XX]*rt[XX]) * cfx + ( b*rja[ZZ]-rm[YY]*rt[XX]) * cfy + (-b*rja[YY]-rm[ZZ]*rt[XX]) * cfz;
    fl[YY] = (-b*rja[ZZ]-rm[XX]*rt[YY]) * cfx + (          -rm[YY]*rt[YY]) * cfy + ( b*rja[XX]-rm[ZZ]*rt[YY]) * cfz;
    fl[ZZ] = ( b*rja[YY]-rm[XX]*rt[ZZ]) * cfx + (-b*rja[XX]-rm[YY]*rt[ZZ]) * cfy + (          -rm[ZZ]*rt[ZZ]) * cfz;
    /* 36 flops */

    f[ai][XX] += fv[XX] - fj[XX] - fk[XX] - fl[XX];
    f[ai][YY] += fv[YY] - fj[YY] - fk[YY] - fl[YY];
    f[ai][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ] - fl[ZZ];
    rvec_inc(f[aj], fj);
    rvec_inc(f[ak], fk);
    rvec_inc(f[al], fl);
    /* 21 flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, av), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sij = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, ai), di);
        sik = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, al), SHIFT_IVEC(g, ai), di);
        sil = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift && (svi != CENTRAL || sij != CENTRAL || sik != CENTRAL || sil != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        fshift[CENTRAL][XX] += fv[XX] - fj[XX] - fk[XX] - fl[XX];
        fshift[CENTRAL][YY] += fv[YY] - fj[YY] - fk[YY] - fl[YY];
        fshift[CENTRAL][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ] - fl[ZZ];
        rvec_inc(fshift[sij], fj);
        rvec_inc(fshift[sik], fk);
        rvec_inc(fshift[sil], fl);
    }

    if (VirCorr)
    {
        rvec xiv;
        int  i, j;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                dxdf[i][j] += -xiv[i]*fv[j] + xij[i]*fj[j] + xik[i]*fk[j] + xil[i]*fl[j];
            }
        }
    }

    /* Total: 207 flops (Yuck!) */
}


static int spread_vsiten(const t_iatom ia[], const t_iparams ip[],
                         const rvec x[],
                         rvec f[], rvec fshift[],
                         const t_pbc *pbc, const t_graph *g)
{
    rvec xv, dx, fi;
    int  n3, av, i, ai;
    real a;
    ivec di;
    int  siv;

    n3 = 3*ip[ia[0]].vsiten.n;
    av = ia[1];
    copy_rvec(x[av], xv);

    for (i = 0; i < n3; i += 3)
    {
        ai = ia[i+2];
        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, av), di);
            siv = IVEC2IS(di);
        }
        else if (pbc)
        {
            siv = pbc_dx_aiuc(pbc, x[ai], xv, dx);
        }
        else
        {
            siv = CENTRAL;
        }
        a = ip[ia[i]].vsiten.a;
        svmul(a, f[av], fi);
        rvec_inc(f[ai], fi);
        if (fshift && siv != CENTRAL)
        {
            rvec_inc(fshift[siv], fi);
            rvec_dec(fshift[CENTRAL], fi);
        }
        /* 6 Flops */
    }

    return n3;
}


static int vsite_count(const t_ilist *ilist, int ftype)
{
    if (ftype == F_VSITEN)
    {
        return ilist[ftype].nr/3;
    }
    else
    {
        return ilist[ftype].nr/(1 + interaction_function[ftype].nratoms);
    }
}

static void spread_vsite_f_thread(const gmx_vsite_t *vsite,
                                  const rvec x[],
                                  rvec f[], rvec *fshift,
                                  gmx_bool VirCorr, matrix dxdf,
                                  t_iparams ip[], const t_ilist ilist[],
                                  const t_graph *g, const t_pbc *pbc_null)
{
    const PbcMode  pbcMode   = getPbcMode(pbc_null, vsite);
    /* We need another pbc pointer, as with charge groups we switch per vsite */
    const t_pbc   *pbc_null2 = pbc_null;
    const int     *vsite_pbc = nullptr;

    /* this loop goes backwards to be able to build *
     * higher type vsites from lower types         */
    for (int ftype = c_ftypeVsiteEnd - 1; ftype >= c_ftypeVsiteStart; ftype--)
    {
        if (ilist[ftype].nr == 0)
        {
            continue;
        }

        {   // TODO remove me
            int            nra = interaction_function[ftype].nratoms;
            int            inc = 1 + nra;
            int            nr  = ilist[ftype].nr;

            const t_iatom *ia = ilist[ftype].iatoms;

            if (pbcMode == PbcMode::all)
            {
                pbc_null2 = pbc_null;
            }
            else if (pbcMode == PbcMode::chargeGroup)
            {
                vsite_pbc = vsite->vsite_pbc_loc[ftype - c_ftypeVsiteStart];
            }

            for (int i = 0; i < nr; )
            {
                if (vsite_pbc != nullptr)
                {
                    if (vsite_pbc[i/(1 + nra)] > -2)
                    {
                        pbc_null2 = pbc_null;
                    }
                    else
                    {
                        pbc_null2 = nullptr;
                    }
                }

                int tp = ia[0];

                /* Constants for constructing */
                real a1, b1, c1;
                a1 = ip[tp].vsite.a;
                /* Construct the vsite depending on type */
                switch (ftype)
                {
                    case F_VSITE2:
                        spread_vsite2(ia, a1, x, f, fshift, pbc_null2, g);
                        break;
                    case F_VSITE3:
                        b1 = ip[tp].vsite.b;
                        spread_vsite3(ia, a1, b1, x, f, fshift, pbc_null2, g);
                        break;
                    case F_VSITE3FD:
                        b1 = ip[tp].vsite.b;
                        spread_vsite3FD(ia, a1, b1, x, f, fshift, VirCorr, dxdf, pbc_null2, g);
                        break;
                    case F_VSITE3FAD:
                        b1 = ip[tp].vsite.b;
                        spread_vsite3FAD(ia, a1, b1, x, f, fshift, VirCorr, dxdf, pbc_null2, g);
                        break;
                    case F_VSITE3OUT:
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        spread_vsite3OUT(ia, a1, b1, c1, x, f, fshift, VirCorr, dxdf, pbc_null2, g);
                        break;
                    case F_VSITE4FD:
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        spread_vsite4FD(ia, a1, b1, c1, x, f, fshift, VirCorr, dxdf, pbc_null2, g);
                        break;
                    case F_VSITE4FDN:
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        spread_vsite4FDN(ia, a1, b1, c1, x, f, fshift, VirCorr, dxdf, pbc_null2, g);
                        break;
                    case F_VSITEN:
                        inc = spread_vsiten(ia, ip, x, f, fshift, pbc_null2, g);
                        break;
                    default:
                        gmx_fatal(FARGS, "No such vsite type %d in %s, line %d",
                                  ftype, __FILE__, __LINE__);
                }
                clear_rvec(f[ia[1]]);

                /* Increment loop variables */
                i  += inc;
                ia += inc;
            }
        }
    }
}

/*! \brief Clears the task force buffer elements that are written by task idTask */
static void clearTaskForceBufferUsedElements(InterdependentTask *idTask)
{
    int ntask = idTask->spreadTask.size();
    for (int ti = 0; ti < ntask; ti++)
    {
        const AtomIndex *atomList = &idTask->atomIndex[idTask->spreadTask[ti]];
        int              natom    = atomList->atom.size();
        RVec            *force    = idTask->force.data();
        for (int i = 0; i < natom; i++)
        {
            clear_rvec(force[atomList->atom[i]]);
        }
    }
}

void spread_vsite_f(const gmx_vsite_t *vsite,
                    const rvec * gmx_restrict x,
                    rvec * gmx_restrict f, rvec * gmx_restrict fshift,
                    gmx_bool VirCorr, matrix vir,
                    t_nrnb *nrnb, const t_idef *idef,
                    int ePBC, gmx_bool bMolPBC, const t_graph *g, const matrix box,
                    t_commrec *cr)
{
    const bool useDomdec = vsite->useDomdec;
    GMX_ASSERT(!useDomdec || (cr != nullptr && DOMAINDECOMP(cr)), "When vsites are set up with domain decomposition, we need a valid commrec");

    t_pbc pbc, *pbc_null;

    /* We only need to do pbc when we have inter-cg vsites */
    if ((useDomdec || bMolPBC) && vsite->n_intercg_vsite)
    {
        /* This is wasting some CPU time as we now do this multiple times
         * per MD step.
         */
        pbc_null = set_pbc_dd(&pbc, ePBC, useDomdec ? cr->dd->nc : nullptr, FALSE, box);
    }
    else
    {
        pbc_null = nullptr;
    }

    if (useDomdec)
    {
        dd_clear_f_vsites(cr->dd, f);
    }

    if (vsite->nthreads == 1)
    {
        matrix dxdf;
        if (VirCorr)
        {
            clear_mat(dxdf);
        }
        spread_vsite_f_thread(vsite,
                              x, f, fshift,
                              VirCorr, dxdf,
                              idef->iparams, idef->il,
                              g, pbc_null);

        if (VirCorr)
        {
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    vir[i][j] += -0.5*dxdf[i][j];
                }
            }
        }
    }
    else
    {
        /* First spread the vsites that might depend on non-local vsites */
        if (VirCorr)
        {
            clear_mat(vsite->tData[vsite->nthreads]->dxdf);
        }
        spread_vsite_f_thread(vsite,
                              x, f, fshift,
                              VirCorr, vsite->tData[vsite->nthreads]->dxdf,
                              idef->iparams,
                              vsite->tData[vsite->nthreads]->ilist,
                              g, pbc_null);

#pragma omp parallel num_threads(vsite->nthreads)
        {
            try
            {
                int          thread = gmx_omp_get_thread_num();
                VsiteThread *tData  = vsite->tData[thread];

                rvec        *fshift_t;
                if (thread == 0 || fshift == nullptr)
                {
                    fshift_t = fshift;
                }
                else
                {
                    fshift_t = tData->fshift;

                    for (int i = 0; i < SHIFTS; i++)
                    {
                        clear_rvec(fshift_t[i]);
                    }
                }
                if (VirCorr)
                {
                    clear_mat(tData->dxdf);
                }

                if (tData->useInterdependentTask)
                {
                    /* Spread the vsites that spread outside our local range.
                     * This is done using a thread-local force buffer force.
                     * First we need to copy the input vsite forces to force.
                     */
                    InterdependentTask *idTask = &tData->idTask;

                    /* Clear the buffer elements set by our task during
                     * the last call to spread_vsite_f.
                     */
                    clearTaskForceBufferUsedElements(idTask);

                    int nvsite = idTask->vsite.size();
                    for (int i = 0; i < nvsite; i++)
                    {
                        copy_rvec(f[idTask->vsite[i]],
                                  idTask->force[idTask->vsite[i]]);
                    }
                    spread_vsite_f_thread(vsite,
                                          x, as_rvec_array(idTask->force.data()), fshift_t,
                                          VirCorr, tData->dxdf,
                                          idef->iparams,
                                          tData->idTask.ilist,
                                          g, pbc_null);

                    /* We need a barrier before reducing forces below
                     * that have been produced by a different thread above.
                     */
#pragma omp barrier

                    /* Loop over all thread task and reduce forces they
                     * produced on atoms that fall in our range.
                     * Note that atomic reduction would be a simpler solution,
                     * but that might not have good support on all platforms.
                     */
                    int ntask = idTask->reduceTask.size();
                    for (int ti = 0; ti < ntask; ti++)
                    {
                        const InterdependentTask *idt_foreign = &vsite->tData[idTask->reduceTask[ti]]->idTask;
                        const AtomIndex          *atomList    = &idt_foreign->atomIndex[thread];
                        const RVec               *f_foreign   = idt_foreign->force.data();

                        int natom = atomList->atom.size();
                        for (int i = 0; i < natom; i++)
                        {
                            int ind = atomList->atom[i];
                            rvec_inc(f[ind], f_foreign[ind]);
                            /* Clearing of f_foreign is done at the next step */
                        }
                    }
                    /* Clear the vsite forces, both in f and force */
                    for (int i = 0; i < nvsite; i++)
                    {
                        int ind = tData->idTask.vsite[i];
                        clear_rvec(f[ind]);
                        clear_rvec(tData->idTask.force[ind]);
                    }
                }

                /* Spread the vsites that spread locally only */
                spread_vsite_f_thread(vsite,
                                      x, f, fshift_t,
                                      VirCorr, tData->dxdf,
                                      idef->iparams,
                                      tData->ilist,
                                      g, pbc_null);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        }

        if (fshift != nullptr)
        {
            for (int th = 1; th < vsite->nthreads; th++)
            {
                for (int i = 0; i < SHIFTS; i++)
                {
                    rvec_inc(fshift[i], vsite->tData[th]->fshift[i]);
                }
            }
        }

        if (VirCorr)
        {
            for (int th = 0; th < vsite->nthreads + 1; th++)
            {
                /* MSVC doesn't like matrix references, so we use a pointer */
                const matrix *dxdf = &vsite->tData[th]->dxdf;

                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        vir[i][j] += -0.5*(*dxdf)[i][j];
                    }
                }
            }
        }
    }

    if (useDomdec)
    {
        dd_move_f_vsites(cr->dd, f, fshift);
    }

    inc_nrnb(nrnb, eNR_VSITE2,   vsite_count(idef->il, F_VSITE2));
    inc_nrnb(nrnb, eNR_VSITE3,   vsite_count(idef->il, F_VSITE3));
    inc_nrnb(nrnb, eNR_VSITE3FD, vsite_count(idef->il, F_VSITE3FD));
    inc_nrnb(nrnb, eNR_VSITE3FAD, vsite_count(idef->il, F_VSITE3FAD));
    inc_nrnb(nrnb, eNR_VSITE3OUT, vsite_count(idef->il, F_VSITE3OUT));
    inc_nrnb(nrnb, eNR_VSITE4FD, vsite_count(idef->il, F_VSITE4FD));
    inc_nrnb(nrnb, eNR_VSITE4FDN, vsite_count(idef->il, F_VSITE4FDN));
    inc_nrnb(nrnb, eNR_VSITEN,   vsite_count(idef->il, F_VSITEN));
}

/*! \brief Returns the an array with charge-group indices for each atom
 *
 * \param[in] chargeGroups  The charge group block struct
 */
static std::vector<int> atom2cg(const t_block &chargeGroups)
{
    std::vector<int> a2cg(chargeGroups.index[chargeGroups.nr], 0);

    for (int chargeGroup = 0; chargeGroup < chargeGroups.nr; chargeGroup++)
    {
        std::fill(a2cg.begin() + chargeGroups.index[chargeGroup],
                  a2cg.begin() + chargeGroups.index[chargeGroup + 1],
                  chargeGroup);
    }

    return a2cg;
}

int count_intercg_vsites(const gmx_mtop_t *mtop)
{
    gmx_molblock_t *molb;
    gmx_moltype_t  *molt;
    int             n_intercg_vsite;

    n_intercg_vsite = 0;
    for (int mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        molt = &mtop->moltype[molb->type];

        std::vector<int> a2cg = atom2cg(molt->cgs);
        for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
        {
            int            nral = NRAL(ftype);
            t_ilist       *il   = &molt->ilist[ftype];
            const t_iatom *ia   = il->iatoms;
            for (int i = 0; i < il->nr; i += 1 + nral)
            {
                int cg = a2cg[ia[1+i]];
                for (int a = 1; a < nral; a++)
                {
                    if (a2cg[ia[1+a]] != cg)
                    {
                        n_intercg_vsite += molb->nmol;
                        break;
                    }
                }
            }
        }
    }

    return n_intercg_vsite;
}

static int **get_vsite_pbc(const t_iparams *iparams, const t_ilist *ilist,
                           const t_atom *atom, const t_mdatoms *md,
                           const t_block &cgs)
{
    /* Make an atom to charge group index */
    std::vector<int> a2cg = atom2cg(cgs);

    /* Make an array that tells if the pbc of an atom is set */
    std::vector<bool> pbc_set(cgs.index[cgs.nr], false);
    /* PBC is set for all non vsites */
    for (int a = 0; a < cgs.index[cgs.nr]; a++)
    {
        if ((atom && atom[a].ptype != eptVSite) ||
            (md   && md->ptype[a]  != eptVSite))
        {
            pbc_set[a] = true;
        }
    }

    int **vsite_pbc;
    snew(vsite_pbc, F_VSITEN-F_VSITE2+1);

    for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
    {
        {   // TODO remove me
            int            nral = NRAL(ftype);
            const t_ilist *il   = &ilist[ftype];
            const t_iatom *ia   = il->iatoms;
            int           *vsite_pbc_f;

            snew(vsite_pbc[ftype-F_VSITE2], il->nr/(1 + nral));
            vsite_pbc_f = vsite_pbc[ftype-F_VSITE2];

            int i = 0;
            while (i < il->nr)
            {
                int     vsi   = i/(1 + nral);
                t_iatom vsite = ia[i+1];
                int     cg_v  = a2cg[vsite];
                /* A value of -2 signals that this vsite and its contructing
                 * atoms are all within the same cg, so no pbc is required.
                 */
                vsite_pbc_f[vsi] = -2;
                /* Check if constructing atoms are outside the vsite's cg */
                int nc3 = 0;
                if (ftype == F_VSITEN)
                {
                    nc3 = 3*iparams[ia[i]].vsiten.n;
                    for (int j = 0; j < nc3; j += 3)
                    {
                        if (a2cg[ia[i+j+2]] != cg_v)
                        {
                            vsite_pbc_f[vsi] = -1;
                        }
                    }
                }
                else
                {
                    for (int a = 1; a < nral; a++)
                    {
                        if (a2cg[ia[i+1+a]] != cg_v)
                        {
                            vsite_pbc_f[vsi] = -1;
                        }
                    }
                }
                if (vsite_pbc_f[vsi] == -1)
                {
                    /* Check if this is the first processed atom of a vsite only cg */
                    gmx_bool bViteOnlyCG_and_FirstAtom = TRUE;
                    for (int a = cgs.index[cg_v]; a < cgs.index[cg_v + 1]; a++)
                    {
                        /* Non-vsites already have pbc set, so simply check for pbc_set */
                        if (pbc_set[a])
                        {
                            bViteOnlyCG_and_FirstAtom = FALSE;
                            break;
                        }
                    }
                    if (bViteOnlyCG_and_FirstAtom)
                    {
                        /* First processed atom of a vsite only charge group.
                         * The pbc of the input coordinates to construct_vsites
                         * should be preserved.
                         */
                        vsite_pbc_f[vsi] = vsite;
                    }
                    else if (cg_v != a2cg[ia[1+i+1]])
                    {
                        /* This vsite has a different charge group index
                         * than it's first constructing atom
                         * and the charge group has more than one atom,
                         * search for the first normal particle
                         * or vsite that already had its pbc defined.
                         * If nothing is found, use full pbc for this vsite.
                         */
                        for (int a = cgs.index[cg_v]; a < cgs.index[cg_v + 1]; a++)
                        {
                            if (a != vsite && pbc_set[a])
                            {
                                vsite_pbc_f[vsi] = a;
                                if (gmx_debug_at)
                                {
                                    fprintf(debug, "vsite %d match pbc with atom %d\n",
                                            vsite+1, a+1);
                                }
                                break;
                            }
                        }
                        if (gmx_debug_at)
                        {
                            fprintf(debug, "vsite atom %d  cg %d - %d pbc atom %d\n",
                                    vsite+1, cgs.index[cg_v] + 1, cgs.index[cg_v + 1],
                                    vsite_pbc_f[vsi] + 1);
                        }
                    }
                }
                if (ftype == F_VSITEN)
                {
                    /* The other entries in vsite_pbc_f are not used for center vsites */
                    i += nc3;
                }
                else
                {
                    i += 1 + nral;
                }

                /* This vsite now has its pbc defined */
                pbc_set[vsite] = true;
            }
        }
    }

    return vsite_pbc;
}


gmx_vsite_t *initVsite(const gmx_mtop_t &mtop,
                       t_commrec        *cr)
{
    GMX_RELEASE_ASSERT(cr != nullptr, "We need a valid commrec");

    /* check if there are vsites */
    int nvsite = 0;
    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            GMX_ASSERT(ftype >= c_ftypeVsiteStart && ftype < c_ftypeVsiteEnd, "c_ftypeVsiteStart and/or c_ftypeVsiteEnd do not have correct values");

            nvsite += gmx_mtop_ftype_count(&mtop, ftype);
        }
        else
        {
            GMX_ASSERT(ftype < c_ftypeVsiteStart || ftype >= c_ftypeVsiteEnd, "c_ftypeVsiteStart and/or c_ftypeVsiteEnd do not have correct values");
        }
    }

    if (nvsite == 0)
    {
        return nullptr;
    }

    gmx_vsite_t *vsite = new(gmx_vsite_t);

    vsite->n_intercg_vsite   = count_intercg_vsites(&mtop);

    vsite->bHaveChargeGroups = (ncg_mtop(&mtop) < mtop.natoms);

    vsite->useDomdec         = (DOMAINDECOMP(cr) && cr->dd->nnodes > 1);

    /* If we don't have charge groups, the vsite follows its own pbc.
     *
     * With charge groups, each vsite needs to follow the pbc of the charge
     * group. Thus we need to keep track of PBC. Here we assume that without
     * domain decomposition all molecules are whole (which will not be
     * the case with periodic molecules).
     */
    if (vsite->bHaveChargeGroups &&
        vsite->n_intercg_vsite > 0 &&
        DOMAINDECOMP(cr))
    {
        vsite->nvsite_pbc_molt = mtop.nmoltype;
        snew(vsite->vsite_pbc_molt, vsite->nvsite_pbc_molt);
        for (int mt = 0; mt < mtop.nmoltype; mt++)
        {
            const gmx_moltype_t &molt = mtop.moltype[mt];
            vsite->vsite_pbc_molt[mt] = get_vsite_pbc(mtop.ffparams.iparams,
                                                      molt.ilist,
                                                      molt.atoms.atom, nullptr,
                                                      molt.cgs);
        }

        snew(vsite->vsite_pbc_loc_nalloc, c_ftypeVsiteEnd - c_ftypeVsiteStart);
        snew(vsite->vsite_pbc_loc,        c_ftypeVsiteEnd - c_ftypeVsiteStart);
    }
    else
    {
        vsite->vsite_pbc_molt = nullptr;
        vsite->vsite_pbc_loc  = nullptr;
    }

    vsite->nthreads = gmx_omp_nthreads_get(emntVSITE);

    if (vsite->nthreads > 1)
    {
        /* We need one extra thread data structure for the overlap vsites */
        snew(vsite->tData, vsite->nthreads + 1);
#pragma omp parallel for num_threads(vsite->nthreads) schedule(static)
        for (int thread = 0; thread < vsite->nthreads; thread++)
        {
            try
            {
                vsite->tData[thread] = new VsiteThread;

                InterdependentTask *idTask = &vsite->tData[thread]->idTask;
                idTask->nuse               = 0;
                idTask->atomIndex.resize(vsite->nthreads);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        }
        if (vsite->nthreads > 1)
        {
            vsite->tData[vsite->nthreads] = new VsiteThread;
        }
    }

    vsite->taskIndex       = nullptr;
    vsite->taskIndexNalloc = 0;

    return vsite;
}

static gmx_inline void flagAtom(InterdependentTask *idTask, int atom,
                                int thread, int nthread, int natperthread)
{
    if (!idTask->use[atom])
    {
        idTask->use[atom] = true;
        thread            = atom/natperthread;
        /* Assign all non-local atom force writes to thread 0 */
        if (thread >= nthread)
        {
            thread        = 0;
        }
        idTask->atomIndex[thread].atom.push_back(atom);
    }
}

/*\brief Here we try to assign all vsites that are in our local range.
 *
 * Our task local atom range is tData->rangeStart - tData->rangeEnd.
 * Vsites that depend only on local atoms, as indicated by taskIndex[]==thread,
 * are assigned to task tData->ilist. Vsites that depend on non-local atoms
 * but not on other vsites are assigned to task tData->id_task.ilist.
 * taskIndex[] is set for all vsites in our range, either to our local tasks
 * or to the single last task as taskIndex[]=2*nthreads.
 */
static void assignVsitesToThread(VsiteThread           *tData,
                                 int                    thread,
                                 int                    nthread,
                                 int                    natperthread,
                                 int                   *taskIndex,
                                 const t_ilist         *ilist,
                                 const t_iparams       *ip,
                                 const unsigned short  *ptype)
{
    for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
    {
        tData->ilist[ftype].nr        = 0;
        tData->idTask.ilist[ftype].nr = 0;

        int      nral1 = 1 + NRAL(ftype);
        int      inc   = nral1;
        t_iatom *iat   = ilist[ftype].iatoms;
        for (int i = 0; i < ilist[ftype].nr; )
        {
            if (ftype == F_VSITEN)
            {
                /* The 3 below is from 1+NRAL(ftype)=3 */
                inc = ip[iat[i]].vsiten.n*3;
            }

            if (iat[1 + i] <  tData->rangeStart ||
                iat[1 + i] >= tData->rangeEnd)
            {
                /* This vsite belongs to a different thread */
                i += inc;
                continue;
            }

            /* We would like to assign this vsite to task thread,
             * but it might depend on atoms outside the atom range of thread
             * or on another vsite not assigned to task thread.
             */
            int task = thread;
            if (ftype != F_VSITEN)
            {
                for (int j = i + 2; j < i + nral1; j++)
                {
                    /* Do a range check to avoid a harmless race on taskIndex */
                    if (iat[j] <  tData->rangeStart ||
                        iat[j] >= tData->rangeEnd ||
                        taskIndex[iat[j]] != thread)
                    {
                        if (!tData->useInterdependentTask ||
                            ptype[iat[j]] == eptVSite)
                        {
                            /* At least one constructing atom is a vsite
                             * that is not assigned to the same thread.
                             * Put this vsite into a separate task.
                             */
                            task = 2*nthread;
                            break;
                        }

                        /* There are constructing atoms outside our range,
                         * put this vsite into a second task to be executed
                         * on the same thread. During construction no barrier
                         * is needed between the two tasks on the same thread.
                         * During spreading we need to run this task with
                         * an additional thread-local intermediate force buffer
                         * (or atomic reduction) and a barrier between the two
                         * tasks.
                         */
                        task = nthread + thread;
                    }
                }
            }
            else
            {
                for (int j = i + 2; j < i + inc; j += 3)
                {
                    /* Do a range check to avoid a harmless race on taskIndex */
                    if (iat[j] <  tData->rangeStart ||
                        iat[j] >= tData->rangeEnd ||
                        taskIndex[iat[j]] != thread)
                    {
                        GMX_ASSERT(ptype[iat[j]] != eptVSite, "A vsite to be assigned in assignVsitesToThread has a vsite as a constructing atom that does not belong to our task, such vsites should be assigned to the single 'master' task");

                        task = nthread + thread;
                    }
                }
            }

            /* Update this vsite's thread index entry */
            taskIndex[iat[1+i]] = task;

            if (task == thread || task == nthread + thread)
            {
                /* Copy this vsite to the thread data struct of thread */
                t_ilist *il_task;
                if (task == thread)
                {
                    il_task = &tData->ilist[ftype];
                }
                else
                {
                    il_task = &tData->idTask.ilist[ftype];
                }
                /* Ensure we have sufficient memory allocated */
                if (il_task->nr + inc > il_task->nalloc)
                {
                    il_task->nalloc = over_alloc_large(il_task->nr + inc);
                    srenew(il_task->iatoms, il_task->nalloc);
                }
                /* Copy the vsite data to the thread-task local array */
                for (int j = i; j < i + inc; j++)
                {
                    il_task->iatoms[il_task->nr++] = iat[j];
                }
                if (task == nthread + thread)
                {
                    /* This vsite write outside our own task force block.
                     * Put it into the interdependent task list and flag
                     * the atoms involved for reduction.
                     */
                    tData->idTask.vsite.push_back(iat[i + 1]);
                    if (ftype != F_VSITEN)
                    {
                        for (int j = i + 2; j < i + nral1; j++)
                        {
                            flagAtom(&tData->idTask, iat[j],
                                     thread, nthread, natperthread);
                        }
                    }
                    else
                    {
                        for (int j = i + 2; j < i + inc; j += 3)
                        {
                            flagAtom(&tData->idTask, iat[j],
                                     thread, nthread, natperthread);
                        }
                    }
                }
            }

            i += inc;
        }
    }
}

/*! \brief Assign all vsites with taskIndex[]==task to task tData */
static void assignVsitesToSingleTask(VsiteThread     *tData,
                                     int              task,
                                     const int       *taskIndex,
                                     const t_ilist   *ilist,
                                     const t_iparams *ip)
{
    for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
    {
        tData->ilist[ftype].nr        = 0;
        tData->idTask.ilist[ftype].nr = 0;

        int      nral1   = 1 + NRAL(ftype);
        int      inc     = nral1;
        t_iatom *iat     = ilist[ftype].iatoms;
        t_ilist *il_task = &tData->ilist[ftype];

        for (int i = 0; i < ilist[ftype].nr; )
        {
            if (ftype == F_VSITEN)
            {
                /* The 3 below is from 1+NRAL(ftype)=3 */
                inc = ip[iat[i]].vsiten.n*3;
            }
            /* Check if the vsite is assigned to our task */
            if (taskIndex[iat[1 + i]] == task)
            {
                /* Ensure we have sufficient memory allocated */
                if (il_task->nr + inc > il_task->nalloc)
                {
                    il_task->nalloc = over_alloc_large(il_task->nr + inc);
                    srenew(il_task->iatoms, il_task->nalloc);
                }
                /* Copy the vsite data to the thread-task local array */
                for (int j = i; j < i + inc; j++)
                {
                    il_task->iatoms[il_task->nr++] = iat[j];
                }
            }

            i += inc;
        }
    }
}

void split_vsites_over_threads(const t_ilist   *ilist,
                               const t_iparams *ip,
                               const t_mdatoms *mdatoms,
                               gmx_vsite_t     *vsite)
{
    int      vsite_atom_range, natperthread;

    if (vsite->nthreads == 1)
    {
        /* Nothing to do */
        return;
    }

    /* The current way of distributing the vsites over threads in primitive.
     * We divide the atom range 0 - natoms_in_vsite uniformly over threads,
     * without taking into account how the vsites are distributed.
     * Without domain decomposition we at least tighten the upper bound
     * of the range (useful for common systems such as a vsite-protein
     * in 3-site water).
     * With domain decomposition, as long as the vsites are distributed
     * uniformly in each domain along the major dimension, usually x,
     * it will also perform well.
     */
    if (!vsite->useDomdec)
    {
        vsite_atom_range = -1;
        for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
        {
            {   // TODO remove me
                if (ftype != F_VSITEN)
                {
                    int            nral1 = 1 + NRAL(ftype);
                    const t_iatom *iat   = ilist[ftype].iatoms;
                    for (int i = 0; i < ilist[ftype].nr; i += nral1)
                    {
                        for (int j = i + 1; j < i + nral1; j++)
                        {
                            vsite_atom_range = std::max(vsite_atom_range, iat[j]);
                        }
                    }
                }
                else
                {
                    int            vs_ind_end;

                    const t_iatom *iat = ilist[ftype].iatoms;

                    int            i = 0;
                    while (i < ilist[ftype].nr)
                    {
                        /* The 3 below is from 1+NRAL(ftype)=3 */
                        vs_ind_end = i + ip[iat[i]].vsiten.n*3;

                        vsite_atom_range = std::max(vsite_atom_range, iat[i+1]);
                        while (i < vs_ind_end)
                        {
                            vsite_atom_range = std::max(vsite_atom_range, iat[i+2]);
                            i               += 3;
                        }
                    }
                }
            }
        }
        vsite_atom_range++;
        natperthread     = (vsite_atom_range + vsite->nthreads - 1)/vsite->nthreads;
    }
    else
    {
        /* Any local or not local atom could be involved in virtual sites.
         * But since we usually have very few non-local virtual sites
         * (only non-local vsites that depend on local vsites),
         * we distribute the local atom range equally over the threads.
         * When assigning vsites to threads, we should take care that the last
         * threads also covers the non-local range.
         */
        vsite_atom_range = mdatoms->nr;
        natperthread     = (mdatoms->homenr + vsite->nthreads - 1)/vsite->nthreads;
    }

    if (debug)
    {
        fprintf(debug, "virtual site thread dist: natoms %d, range %d, natperthread %d\n", mdatoms->nr, vsite_atom_range, natperthread);
    }

    /* To simplify the vsite assignment, we make an index which tells us
     * to which task particles, both non-vsites and vsites, are assigned.
     */
    if (mdatoms->nr > vsite->taskIndexNalloc)
    {
        vsite->taskIndexNalloc = over_alloc_large(mdatoms->nr);
        srenew(vsite->taskIndex, vsite->taskIndexNalloc);
    }

    /* Initialize the task index array. Here we assign the non-vsite
     * particles to task=thread, so we easily figure out if vsites
     * depend on local and/or non-local particles in assignVsitesToThread.
     */
    int *taskIndex = vsite->taskIndex;
    {
        int  thread = 0;
        for (int i = 0; i < mdatoms->nr; i++)
        {
            if (mdatoms->ptype[i] == eptVSite)
            {
                /* vsites are not assigned to a task yet */
                taskIndex[i] = -1;
            }
            else
            {
                /* assign non-vsite particles to task thread */
                taskIndex[i] = thread;
            }
            if (i == (thread + 1)*natperthread && thread < vsite->nthreads)
            {
                thread++;
            }
        }
    }

#pragma omp parallel num_threads(vsite->nthreads)
    {
        try
        {
            int          thread = gmx_omp_get_thread_num();
            VsiteThread *tData  = vsite->tData[thread];

            /* Clear the buffer use flags that were set before */
            if (tData->useInterdependentTask)
            {
                InterdependentTask *idTask = &tData->idTask;

                /* To avoid an extra OpenMP barrier in spread_vsite_f,
                 * we clear the force buffer at the next step,
                 * so we need to do it here as well.
                 */
                clearTaskForceBufferUsedElements(idTask);

                idTask->vsite.resize(0);
                for (int t = 0; t < vsite->nthreads; t++)
                {
                    AtomIndex *atomIndex = &idTask->atomIndex[t];
                    int        natom     = atomIndex->atom.size();
                    for (int i = 0; i < natom; i++)
                    {
                        idTask->use[atomIndex->atom[i]] = false;
                    }
                    atomIndex->atom.resize(0);
                }
                idTask->nuse = 0;
            }

            /* To avoid large f_buf allocations of #threads*vsite_atom_range
             * we don't use task2 with more than 200000 atoms. This doesn't
             * affect performance, since with such a large range relatively few
             * vsites will end up in the separate task.
             * Note that useTask2 should be the same for all threads.
             */
            tData->useInterdependentTask = (vsite_atom_range <= 200000);
            if (tData->useInterdependentTask)
            {
                size_t              natoms_use_in_vsites = vsite_atom_range;
                InterdependentTask *idTask               = &tData->idTask;
                /* To avoid resizing and re-clearing every nstlist steps,
                 * we never down size the force buffer.
                 */
                if (natoms_use_in_vsites > idTask->force.size() ||
                    natoms_use_in_vsites > idTask->use.size())
                {
                    idTask->force.resize(natoms_use_in_vsites, { 0, 0, 0 });
                    idTask->use.resize(natoms_use_in_vsites, false);
                }
            }

            /* Assign all vsites that can execute independently on threads */
            tData->rangeStart     =  thread     *natperthread;
            if (thread < vsite->nthreads - 1)
            {
                tData->rangeEnd   = (thread + 1)*natperthread;
            }
            else
            {
                /* The last thread should cover up to the end of the range */
                tData->rangeEnd   = mdatoms->nr;
            }
            assignVsitesToThread(tData,
                                 thread, vsite->nthreads,
                                 natperthread,
                                 taskIndex,
                                 ilist, ip, mdatoms->ptype);

            if (tData->useInterdependentTask)
            {
                /* In the worst case, all tasks write to force ranges of
                 * all other tasks, leading to #tasks^2 scaling (this is only
                 * the overhead, the actual flops remain constant).
                 * But in most cases there is far less coupling. To improve
                 * scaling at high thread counts we therefore construct
                 * an index to only loop over the actually affected tasks.
                 */
                InterdependentTask *idTask = &tData->idTask;

                /* Ensure assignVsitesToThread finished on other threads */
#pragma omp barrier

                idTask->spreadTask.resize(0);
                idTask->reduceTask.resize(0);
                for (int t = 0; t < vsite->nthreads; t++)
                {
                    /* Do we write to the force buffer of task t? */
                    if (idTask->atomIndex[t].atom.size() > 0)
                    {
                        idTask->spreadTask.push_back(t);
                    }
                    /* Does task t write to our force buffer? */
                    if (vsite->tData[t]->idTask.atomIndex[thread].atom.size() > 0)
                    {
                        idTask->reduceTask.push_back(t);
                    }
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
    /* Assign all remaining vsites, that will have taskIndex[]=2*vsite->nthreads,
     * to a single task that will not run in parallel with other tasks.
     */
    assignVsitesToSingleTask(vsite->tData[vsite->nthreads],
                             2*vsite->nthreads,
                             taskIndex,
                             ilist, ip);

    if (debug && vsite->nthreads > 1)
    {
        fprintf(debug, "virtual site useInterdependentTask %d, nuse:\n",
                vsite->tData[0]->useInterdependentTask);
        for (int th = 0; th < vsite->nthreads + 1; th++)
        {
            fprintf(debug, " %4d", vsite->tData[th]->idTask.nuse);
        }
        fprintf(debug, "\n");

        for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
        {
            if (ilist[ftype].nr > 0)
            {
                fprintf(debug, "%-20s thread dist:",
                        interaction_function[ftype].longname);
                for (int th = 0; th < vsite->nthreads + 1; th++)
                {
                    fprintf(debug, " %4d %4d ",
                            vsite->tData[th]->ilist[ftype].nr,
                            vsite->tData[th]->idTask.ilist[ftype].nr);
                }
                fprintf(debug, "\n");
            }
        }
    }

#ifndef NDEBUG
    int nrOrig     = vsiteIlistNrCount(ilist);
    int nrThreaded = 0;
    for (int th = 0; th < vsite->nthreads + 1; th++)
    {
        nrThreaded +=
            vsiteIlistNrCount(vsite->tData[th]->ilist) +
            vsiteIlistNrCount(vsite->tData[th]->idTask.ilist);
    }
    GMX_ASSERT(nrThreaded == nrOrig, "The number of virtual sites assigned to all thread task has to match the total number of virtual sites");
#endif
}

void set_vsite_top(gmx_vsite_t          *vsite,
                   const gmx_localtop_t *top,
                   const t_mdatoms      *md)
{
    if (vsite->n_intercg_vsite > 0 && vsite->bHaveChargeGroups)
    {
        vsite->vsite_pbc_loc = get_vsite_pbc(top->idef.iparams,
                                             top->idef.il, nullptr, md,
                                             top->cgs);
    }

    if (vsite->nthreads > 1)
    {
        if (vsite->bHaveChargeGroups)
        {
            gmx_fatal(FARGS, "The combination of threading, virtual sites and charge groups is not implemented");
        }

        split_vsites_over_threads(top->idef.il, top->idef.iparams,
                                  md, vsite);
    }
}

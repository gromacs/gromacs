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
/*! \internal \file
 * \brief Implements the VirtualSitesHandler class and vsite standalone functions
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "vsite.h"

#include <cstdio>

#include <algorithm>
#include <memory>
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
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxomp.h"

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
 * Any remaining vsites are assigned to a separate main thread task.
 */
namespace gmx
{

//! VirialHandling is often used outside VirtualSitesHandler class members
using VirialHandling = VirtualSitesHandler::VirialHandling;

/*! \brief Information on PBC and domain decomposition for virtual sites
 */
struct DomainInfo
{
public:
    //! Constructs without PBC and DD
    DomainInfo() = default;

    //! Constructs with PBC and DD, if !=nullptr
    DomainInfo(PbcType pbcType, bool haveInterUpdateGroupVirtualSites, gmx_domdec_t* domdec) :
        pbcType_(pbcType),
        useMolPbc_(pbcType != PbcType::No && haveInterUpdateGroupVirtualSites),
        domdec_(domdec)
    {
    }

    //! Returns whether we are using domain decomposition with more than 1 DD rank
    bool useDomdec() const { return (domdec_ != nullptr); }

    //! The pbc type
    const PbcType pbcType_ = PbcType::No;
    //! Whether molecules are broken over PBC
    const bool useMolPbc_ = false;
    //! Pointer to the domain decomposition struct, nullptr without PP DD
    const gmx_domdec_t* domdec_ = nullptr;
};

/*! \brief List of atom indices belonging to a task
 */
struct AtomIndex
{
    //! List of atom indices
    std::vector<int> atom;
};

/*! \brief Data structure for thread tasks that use constructing atoms outside their own atom range
 */
struct InterdependentTask
{
    //! The interaction lists, only vsite entries are used
    InteractionLists ilist;
    //! Thread/task-local force buffer
    std::vector<RVec> force;
    //! The atom indices of the vsites of our task
    std::vector<int> vsite;
    //! Flags if elements in force are spread to or not
    std::vector<bool> use;
    //! The number of entries set to true in use
    int nuse = 0;
    //! Array of atoms indices, size nthreads, covering all nuse set elements in use
    std::vector<AtomIndex> atomIndex;
    //! List of tasks (force blocks) this task spread forces to
    std::vector<int> spreadTask;
    //! List of tasks that write to this tasks force block range
    std::vector<int> reduceTask;
};

/*! \brief Vsite thread task data structure
 */
struct VsiteThread
{
    //! Start of atom range of this task
    int rangeStart;
    //! End of atom range of this task
    int rangeEnd;
    //! The interaction lists, only vsite entries are used
    std::array<InteractionList, F_NRE> ilist;
    //! Local fshift accumulation buffer
    std::array<RVec, c_numShiftVectors> fshift;
    //! Local virial dx*df accumulation buffer
    matrix dxdf;
    //! Tells if interdependent task idTask should be used (in addition to the rest of this task), this bool has the same value on all threads
    bool useInterdependentTask;
    //! Data for vsites that involve constructing atoms in the atom range of other threads/tasks
    InterdependentTask idTask;

    /*! \brief Constructor */
    VsiteThread()
    {
        rangeStart = -1;
        rangeEnd   = -1;
        for (auto& elem : fshift)
        {
            elem = { 0.0_real, 0.0_real, 0.0_real };
        }
        clear_mat(dxdf);
        useInterdependentTask = false;
    }
};


/*! \brief Information on how the virtual site work is divided over thread tasks
 */
class ThreadingInfo
{
public:
    //! Constructor, retrieves the number of threads to use from gmx_omp_nthreads.h
    ThreadingInfo();

    //! Returns the number of threads to use for vsite operations
    int numThreads() const { return numThreads_; }

    //! Returns the thread data for the given thread
    const VsiteThread& threadData(int threadIndex) const { return *tData_[threadIndex]; }

    //! Returns the thread data for the given thread
    VsiteThread& threadData(int threadIndex) { return *tData_[threadIndex]; }

    //! Returns the thread data for vsites that depend on non-local vsites
    const VsiteThread& threadDataNonLocalDependent() const { return *tData_[numThreads_]; }

    //! Returns the thread data for vsites that depend on non-local vsites
    VsiteThread& threadDataNonLocalDependent() { return *tData_[numThreads_]; }

    //! Set VSites and distribute VSite work over threads, should be called after DD partitioning
    void setVirtualSites(ArrayRef<const InteractionList> ilist,
                         ArrayRef<const t_iparams>       iparams,
                         int                             numAtoms,
                         int                             homenr,
                         ArrayRef<const ParticleType>    ptype,
                         bool                            useDomdec);

private:
    //! Number of threads used for vsite operations
    const int numThreads_;
    //! Thread local vsites and work structs
    std::vector<std::unique_ptr<VsiteThread>> tData_;
    //! Work array for dividing vsites over threads
    std::vector<int> taskIndex_;
};

/*! \brief Impl class for VirtualSitesHandler
 */
class VirtualSitesHandler::Impl
{
public:
    //! Constructor, domdec should be nullptr without DD
    Impl(const gmx_mtop_t&                 mtop,
         gmx_domdec_t*                     domdec,
         PbcType                           pbcType,
         ArrayRef<const RangePartitioning> updateGroupingPerMoleculeType);

    //! Returns the number of virtual sites acting over multiple update groups
    int numInterUpdategroupVirtualSites() const { return numInterUpdategroupVirtualSites_; }

    //! Set VSites and distribute VSite work over threads, should be called after DD partitioning
    void setVirtualSites(ArrayRef<const InteractionList> ilist,
                         int                             numAtoms,
                         int                             homenr,
                         ArrayRef<const ParticleType>    ptype);

    /*! \brief Create positions of vsite atoms based for the local system
     *
     * \param[in,out] x          The coordinates
     * \param[in,out] v          The velocities, needed if operation requires it
     * \param[in]     box        The box
     * \param[in]     operation  Whether we calculate positions, velocities, or both
     */
    void construct(ArrayRef<RVec> x, ArrayRef<RVec> v, const matrix box, VSiteOperation operation) const;

    /*! \brief Spread the force operating on the vsite atoms on the surrounding atoms.
     *
     * vsite should point to a valid object.
     * The virialHandling parameter determines how virial contributions are handled.
     * If this is set to Linear, shift forces are accumulated into fshift.
     * If this is set to NonLinear, non-linear contributions are added to virial.
     * This non-linear correction is required when the virial is not calculated
     * afterwards from the particle position and forces, but in a different way,
     * as for instance for the PME mesh contribution.
     */
    void spreadForces(ArrayRef<const RVec> x,
                      ArrayRef<RVec>       f,
                      VirialHandling       virialHandling,
                      ArrayRef<RVec>       fshift,
                      matrix               virial,
                      t_nrnb*              nrnb,
                      const matrix         box,
                      gmx_wallcycle*       wcycle);

private:
    //! The number of vsites that cross update groups, when =0 no PBC treatment is needed
    const int numInterUpdategroupVirtualSites_;
    //! PBC and DD information
    const DomainInfo domainInfo_;
    //! The interaction parameters
    const ArrayRef<const t_iparams> iparams_;
    //! The interaction lists
    ArrayRef<const InteractionList> ilists_;
    //! Information for handling vsite threading
    ThreadingInfo threadingInfo_;
};

VirtualSitesHandler::~VirtualSitesHandler() = default;

int VirtualSitesHandler::numInterUpdategroupVirtualSites() const
{
    return impl_->numInterUpdategroupVirtualSites();
}

/*! \brief Returns the sum of the vsite ilist sizes over all vsite types
 *
 * \param[in] ilist  The interaction list
 */
static int vsiteIlistNrCount(ArrayRef<const InteractionList> ilist)
{
    int nr = 0;
    for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
    {
        nr += ilist[ftype].size();
    }

    return nr;
}

//! Computes the distance between xi and xj, pbc is used when pbc!=nullptr
static int pbc_rvec_sub(const t_pbc* pbc, const rvec xi, const rvec xj, rvec dx)
{
    if (pbc)
    {
        return pbc_dx_aiuc(pbc, xi, xj, dx);
    }
    else
    {
        rvec_sub(xi, xj, dx);
        return c_centralShiftIndex;
    }
}

//! Returns the 1/norm(x)
static inline real inverseNorm(const rvec x)
{
    return gmx::invsqrt(iprod(x, x));
}

//! Whether we're calculating the virtual site position
enum class VSiteCalculatePosition
{
    Yes,
    No
};
//! Whether we're calculating the virtual site velocity
enum class VSiteCalculateVelocity
{
    Yes,
    No
};

#ifndef DOXYGEN
/* Vsite construction routines */

template<VSiteCalculatePosition calculatePosition, VSiteCalculateVelocity calculateVelocity>
static void constr_vsite1(const rvec xi, rvec x, const rvec vi, rvec v)
{
    if (calculatePosition == VSiteCalculatePosition::Yes)
    {
        copy_rvec(xi, x);
        /* TOTAL: 0 flops */
    }
    if (calculateVelocity == VSiteCalculateVelocity::Yes)
    {
        copy_rvec(vi, v);
    }
}

template<VSiteCalculatePosition calculatePosition, VSiteCalculateVelocity calculateVelocity>
static void
constr_vsite2(const rvec xi, const rvec xj, rvec x, real a, const t_pbc* pbc, const rvec vi, const rvec vj, rvec v)
{
    const real b = 1 - a;
    /* 1 flop */

    if (calculatePosition == VSiteCalculatePosition::Yes)
    {
        if (pbc)
        {
            rvec dx;
            pbc_dx_aiuc(pbc, xj, xi, dx);
            x[XX] = xi[XX] + a * dx[XX];
            x[YY] = xi[YY] + a * dx[YY];
            x[ZZ] = xi[ZZ] + a * dx[ZZ];
        }
        else
        {
            x[XX] = b * xi[XX] + a * xj[XX];
            x[YY] = b * xi[YY] + a * xj[YY];
            x[ZZ] = b * xi[ZZ] + a * xj[ZZ];
            /* 9 Flops */
        }
        /* TOTAL: 10 flops */
    }
    if (calculateVelocity == VSiteCalculateVelocity::Yes)
    {
        v[XX] = b * vi[XX] + a * vj[XX];
        v[YY] = b * vi[YY] + a * vj[YY];
        v[ZZ] = b * vi[ZZ] + a * vj[ZZ];
    }
}

template<VSiteCalculatePosition calculatePosition, VSiteCalculateVelocity calculateVelocity>
static void
constr_vsite2FD(const rvec xi, const rvec xj, rvec x, real a, const t_pbc* pbc, const rvec vi, const rvec vj, rvec v)
{
    rvec xij = { 0 };
    pbc_rvec_sub(pbc, xj, xi, xij);
    /* 3 flops */

    const real invNormXij = inverseNorm(xij);
    const real b          = a * invNormXij;
    /* 6 + 10 flops */

    if (calculatePosition == VSiteCalculatePosition::Yes)
    {
        x[XX] = xi[XX] + b * xij[XX];
        x[YY] = xi[YY] + b * xij[YY];
        x[ZZ] = xi[ZZ] + b * xij[ZZ];
        /* 6 Flops */
        /* TOTAL: 25 flops */
    }
    if (calculateVelocity == VSiteCalculateVelocity::Yes)
    {
        rvec vij = { 0 };
        rvec_sub(vj, vi, vij);
        const real vijDotXij = iprod(vij, xij);

        v[XX] = vi[XX] + b * (vij[XX] - xij[XX] * vijDotXij * invNormXij * invNormXij);
        v[YY] = vi[YY] + b * (vij[YY] - xij[YY] * vijDotXij * invNormXij * invNormXij);
        v[ZZ] = vi[ZZ] + b * (vij[ZZ] - xij[ZZ] * vijDotXij * invNormXij * invNormXij);
    }
}

template<VSiteCalculatePosition calculatePosition, VSiteCalculateVelocity calculateVelocity>
static void constr_vsite3(const rvec   xi,
                          const rvec   xj,
                          const rvec   xk,
                          rvec         x,
                          real         a,
                          real         b,
                          const t_pbc* pbc,
                          const rvec   vi,
                          const rvec   vj,
                          const rvec   vk,
                          rvec         v)
{
    const real c = 1 - a - b;
    /* 2 flops */

    if (calculatePosition == VSiteCalculatePosition::Yes)
    {
        if (pbc)
        {
            rvec dxj, dxk;

            pbc_dx_aiuc(pbc, xj, xi, dxj);
            pbc_dx_aiuc(pbc, xk, xi, dxk);
            x[XX] = xi[XX] + a * dxj[XX] + b * dxk[XX];
            x[YY] = xi[YY] + a * dxj[YY] + b * dxk[YY];
            x[ZZ] = xi[ZZ] + a * dxj[ZZ] + b * dxk[ZZ];
        }
        else
        {
            x[XX] = c * xi[XX] + a * xj[XX] + b * xk[XX];
            x[YY] = c * xi[YY] + a * xj[YY] + b * xk[YY];
            x[ZZ] = c * xi[ZZ] + a * xj[ZZ] + b * xk[ZZ];
            /* 15 Flops */
        }
        /* TOTAL: 17 flops */
    }
    if (calculateVelocity == VSiteCalculateVelocity::Yes)
    {
        v[XX] = c * vi[XX] + a * vj[XX] + b * vk[XX];
        v[YY] = c * vi[YY] + a * vj[YY] + b * vk[YY];
        v[ZZ] = c * vi[ZZ] + a * vj[ZZ] + b * vk[ZZ];
    }
}

template<VSiteCalculatePosition calculatePosition, VSiteCalculateVelocity calculateVelocity>
static void constr_vsite3FD(const rvec   xi,
                            const rvec   xj,
                            const rvec   xk,
                            rvec         x,
                            real         a,
                            real         b,
                            const t_pbc* pbc,
                            const rvec   vi,
                            const rvec   vj,
                            const rvec   vk,
                            rvec         v)
{
    rvec xij, xjk, temp;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xj, xjk);
    /* 6 flops */

    /* temp goes from i to a point on the line jk */
    temp[XX] = xij[XX] + a * xjk[XX];
    temp[YY] = xij[YY] + a * xjk[YY];
    temp[ZZ] = xij[ZZ] + a * xjk[ZZ];
    /* 6 flops */

    const real invNormTemp = inverseNorm(temp);
    const real c           = b * invNormTemp;
    /* 6 + 10 flops */

    if (calculatePosition == VSiteCalculatePosition::Yes)
    {
        x[XX] = xi[XX] + c * temp[XX];
        x[YY] = xi[YY] + c * temp[YY];
        x[ZZ] = xi[ZZ] + c * temp[ZZ];
        /* 6 Flops */
        /* TOTAL: 34 flops */
    }
    if (calculateVelocity == VSiteCalculateVelocity::Yes)
    {
        rvec vij = { 0 };
        rvec vjk = { 0 };
        rvec_sub(vj, vi, vij);
        rvec_sub(vk, vj, vjk);
        const rvec tempV = { vij[XX] + a * vjk[XX], vij[YY] + a * vjk[YY], vij[ZZ] + a * vjk[ZZ] };
        const real tempDotTempV = iprod(temp, tempV);

        v[XX] = vi[XX] + c * (tempV[XX] - temp[XX] * tempDotTempV * invNormTemp * invNormTemp);
        v[YY] = vi[YY] + c * (tempV[YY] - temp[YY] * tempDotTempV * invNormTemp * invNormTemp);
        v[ZZ] = vi[ZZ] + c * (tempV[ZZ] - temp[ZZ] * tempDotTempV * invNormTemp * invNormTemp);
    }
}

template<VSiteCalculatePosition calculatePosition, VSiteCalculateVelocity calculateVelocity>
static void constr_vsite3FAD(const rvec   xi,
                             const rvec   xj,
                             const rvec   xk,
                             rvec         x,
                             real         a,
                             real         b,
                             const t_pbc* pbc,
                             const rvec   vi,
                             const rvec   vj,
                             const rvec   vk,
                             rvec         v)
{ // Note: a = d * cos(theta)
    //       b = d * sin(theta)
    rvec xij, xjk, xp;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xj, xjk);
    /* 6 flops */

    const real invdij    = inverseNorm(xij);
    const real xijDotXjk = iprod(xij, xjk);
    const real c1        = invdij * invdij * xijDotXjk;
    xp[XX]               = xjk[XX] - c1 * xij[XX];
    xp[YY]               = xjk[YY] - c1 * xij[YY];
    xp[ZZ]               = xjk[ZZ] - c1 * xij[ZZ];
    const real a1        = a * invdij;
    const real invNormXp = inverseNorm(xp);
    const real b1        = b * invNormXp;
    /* 45 */

    if (calculatePosition == VSiteCalculatePosition::Yes)
    {
        x[XX] = xi[XX] + a1 * xij[XX] + b1 * xp[XX];
        x[YY] = xi[YY] + a1 * xij[YY] + b1 * xp[YY];
        x[ZZ] = xi[ZZ] + a1 * xij[ZZ] + b1 * xp[ZZ];
        /* 12 Flops */
        /* TOTAL: 63 flops */
    }

    if (calculateVelocity == VSiteCalculateVelocity::Yes)
    {
        rvec vij = { 0 };
        rvec vjk = { 0 };
        rvec_sub(vj, vi, vij);
        rvec_sub(vk, vj, vjk);

        const real vijDotXjkPlusXijDotVjk = iprod(vij, xjk) + iprod(xij, vjk);
        const real xijDotVij              = iprod(xij, vij);
        const real invNormXij2            = invdij * invdij;

        rvec vp = { 0 };
        vp[XX]  = vjk[XX]
                 - xij[XX] * invNormXij2
                           * (vijDotXjkPlusXijDotVjk - invNormXij2 * xijDotXjk * xijDotVij * 2)
                 - vij[XX] * xijDotXjk * invNormXij2;
        vp[YY] = vjk[YY]
                 - xij[YY] * invNormXij2
                           * (vijDotXjkPlusXijDotVjk - invNormXij2 * xijDotXjk * xijDotVij * 2)
                 - vij[YY] * xijDotXjk * invNormXij2;
        vp[ZZ] = vjk[ZZ]
                 - xij[ZZ] * invNormXij2
                           * (vijDotXjkPlusXijDotVjk - invNormXij2 * xijDotXjk * xijDotVij * 2)
                 - vij[ZZ] * xijDotXjk * invNormXij2;

        const real xpDotVp = iprod(xp, vp);

        v[XX] = vi[XX] + a1 * (vij[XX] - xij[XX] * xijDotVij * invdij * invdij)
                + b1 * (vp[XX] - xp[XX] * xpDotVp * invNormXp * invNormXp);
        v[YY] = vi[YY] + a1 * (vij[YY] - xij[YY] * xijDotVij * invdij * invdij)
                + b1 * (vp[YY] - xp[YY] * xpDotVp * invNormXp * invNormXp);
        v[ZZ] = vi[ZZ] + a1 * (vij[ZZ] - xij[ZZ] * xijDotVij * invdij * invdij)
                + b1 * (vp[ZZ] - xp[ZZ] * xpDotVp * invNormXp * invNormXp);
    }
}

template<VSiteCalculatePosition calculatePosition, VSiteCalculateVelocity calculateVelocity>
static void constr_vsite3OUT(const rvec   xi,
                             const rvec   xj,
                             const rvec   xk,
                             rvec         x,
                             real         a,
                             real         b,
                             real         c,
                             const t_pbc* pbc,
                             const rvec   vi,
                             const rvec   vj,
                             const rvec   vk,
                             rvec         v)
{
    rvec xij, xik, temp;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xi, xik);
    cprod(xij, xik, temp);
    /* 15 Flops */

    if (calculatePosition == VSiteCalculatePosition::Yes)
    {
        x[XX] = xi[XX] + a * xij[XX] + b * xik[XX] + c * temp[XX];
        x[YY] = xi[YY] + a * xij[YY] + b * xik[YY] + c * temp[YY];
        x[ZZ] = xi[ZZ] + a * xij[ZZ] + b * xik[ZZ] + c * temp[ZZ];
        /* 18 Flops */
        /* TOTAL: 33 flops */
    }

    if (calculateVelocity == VSiteCalculateVelocity::Yes)
    {
        rvec vij = { 0 };
        rvec vik = { 0 };
        rvec_sub(vj, vi, vij);
        rvec_sub(vk, vi, vik);

        rvec temp1 = { 0 };
        rvec temp2 = { 0 };
        cprod(vij, xik, temp1);
        cprod(xij, vik, temp2);

        v[XX] = vi[XX] + a * vij[XX] + b * vik[XX] + c * (temp1[XX] + temp2[XX]);
        v[YY] = vi[YY] + a * vij[YY] + b * vik[YY] + c * (temp1[YY] + temp2[YY]);
        v[ZZ] = vi[ZZ] + a * vij[ZZ] + b * vik[ZZ] + c * (temp1[ZZ] + temp2[ZZ]);
    }
}

template<VSiteCalculatePosition calculatePosition, VSiteCalculateVelocity calculateVelocity>
static void constr_vsite4FD(const rvec   xi,
                            const rvec   xj,
                            const rvec   xk,
                            const rvec   xl,
                            rvec         x,
                            real         a,
                            real         b,
                            real         c,
                            const t_pbc* pbc,
                            const rvec   vi,
                            const rvec   vj,
                            const rvec   vk,
                            const rvec   vl,
                            rvec         v)
{
    rvec xij, xjk, xjl, temp;
    real d;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xj, xjk);
    pbc_rvec_sub(pbc, xl, xj, xjl);
    /* 9 flops */

    /* temp goes from i to a point on the plane jkl */
    temp[XX] = xij[XX] + a * xjk[XX] + b * xjl[XX];
    temp[YY] = xij[YY] + a * xjk[YY] + b * xjl[YY];
    temp[ZZ] = xij[ZZ] + a * xjk[ZZ] + b * xjl[ZZ];
    /* 12 flops */

    const real invRm = inverseNorm(temp);
    d                = c * invRm;
    /* 6 + 10 flops */

    if (calculatePosition == VSiteCalculatePosition::Yes)
    {
        x[XX] = xi[XX] + d * temp[XX];
        x[YY] = xi[YY] + d * temp[YY];
        x[ZZ] = xi[ZZ] + d * temp[ZZ];
        /* 6 Flops */
        /* TOTAL: 43 flops */
    }
    if (calculateVelocity == VSiteCalculateVelocity::Yes)
    {
        rvec vij = { 0 };
        rvec vjk = { 0 };
        rvec vjl = { 0 };

        rvec_sub(vj, vi, vij);
        rvec_sub(vk, vj, vjk);
        rvec_sub(vl, vj, vjl);

        rvec vm = { 0 };
        vm[XX]  = vij[XX] + a * vjk[XX] + b * vjl[XX];
        vm[YY]  = vij[YY] + a * vjk[YY] + b * vjl[YY];
        vm[ZZ]  = vij[ZZ] + a * vjk[ZZ] + b * vjl[ZZ];

        const real vmDotRm = iprod(vm, temp);
        v[XX]              = vi[XX] + d * (vm[XX] - temp[XX] * vmDotRm * invRm * invRm);
        v[YY]              = vi[YY] + d * (vm[YY] - temp[YY] * vmDotRm * invRm * invRm);
        v[ZZ]              = vi[ZZ] + d * (vm[ZZ] - temp[ZZ] * vmDotRm * invRm * invRm);
    }
}

template<VSiteCalculatePosition calculatePosition, VSiteCalculateVelocity calculateVelocity>
static void constr_vsite4FDN(const rvec   xi,
                             const rvec   xj,
                             const rvec   xk,
                             const rvec   xl,
                             rvec         x,
                             real         a,
                             real         b,
                             real         c,
                             const t_pbc* pbc,
                             const rvec   vi,
                             const rvec   vj,
                             const rvec   vk,
                             const rvec   vl,
                             rvec         v)
{
    rvec xij, xik, xil, ra, rb, rja, rjb, rm;
    real d;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xi, xik);
    pbc_rvec_sub(pbc, xl, xi, xil);
    /* 9 flops */

    ra[XX] = a * xik[XX];
    ra[YY] = a * xik[YY];
    ra[ZZ] = a * xik[ZZ];

    rb[XX] = b * xil[XX];
    rb[YY] = b * xil[YY];
    rb[ZZ] = b * xil[ZZ];

    /* 6 flops */

    rvec_sub(ra, xij, rja);
    rvec_sub(rb, xij, rjb);
    /* 6 flops */

    cprod(rja, rjb, rm);
    /* 9 flops */

    const real invNormRm = inverseNorm(rm);
    d                    = c * invNormRm;
    /* 5+5+1 flops */

    if (calculatePosition == VSiteCalculatePosition::Yes)
    {
        x[XX] = xi[XX] + d * rm[XX];
        x[YY] = xi[YY] + d * rm[YY];
        x[ZZ] = xi[ZZ] + d * rm[ZZ];
        /* 6 Flops */
        /* TOTAL: 47 flops */
    }

    if (calculateVelocity == VSiteCalculateVelocity::Yes)
    {
        rvec vij = { 0 };
        rvec vik = { 0 };
        rvec vil = { 0 };
        rvec_sub(vj, vi, vij);
        rvec_sub(vk, vi, vik);
        rvec_sub(vl, vi, vil);

        rvec vja = { 0 };
        rvec vjb = { 0 };

        vja[XX] = a * vik[XX] - vij[XX];
        vja[YY] = a * vik[YY] - vij[YY];
        vja[ZZ] = a * vik[ZZ] - vij[ZZ];
        vjb[XX] = b * vil[XX] - vij[XX];
        vjb[YY] = b * vil[YY] - vij[YY];
        vjb[ZZ] = b * vil[ZZ] - vij[ZZ];

        rvec temp1 = { 0 };
        rvec temp2 = { 0 };
        cprod(vja, rjb, temp1);
        cprod(rja, vjb, temp2);

        rvec vm = { 0 };
        vm[XX]  = temp1[XX] + temp2[XX];
        vm[YY]  = temp1[YY] + temp2[YY];
        vm[ZZ]  = temp1[ZZ] + temp2[ZZ];

        const real rmDotVm = iprod(rm, vm);
        v[XX]              = vi[XX] + d * (vm[XX] - rm[XX] * rmDotVm * invNormRm * invNormRm);
        v[YY]              = vi[YY] + d * (vm[YY] - rm[YY] * rmDotVm * invNormRm * invNormRm);
        v[ZZ]              = vi[ZZ] + d * (vm[ZZ] - rm[ZZ] * rmDotVm * invNormRm * invNormRm);
    }
}

template<VSiteCalculatePosition calculatePosition, VSiteCalculateVelocity calculateVelocity>
static int constr_vsiten(const t_iatom*            ia,
                         ArrayRef<const t_iparams> ip,
                         ArrayRef<RVec>            x,
                         const t_pbc*              pbc,
                         ArrayRef<RVec>            v)
{
    rvec x1, dx;
    dvec dsum;
    real a;
    dvec dvsum = { 0 };
    rvec v1    = { 0 };

    const int n3 = 3 * ip[ia[0]].vsiten.n;
    const int av = ia[1];
    int       ai = ia[2];
    if (calculatePosition == VSiteCalculatePosition::Yes)
    {
        copy_rvec(x[ai], x1);
    }
    if (calculateVelocity == VSiteCalculateVelocity::Yes)
    {
        copy_rvec(v[ai], v1);
    }
    clear_dvec(dsum);
    for (int i = 3; i < n3; i += 3)
    {
        ai = ia[i + 2];
        a  = ip[ia[i]].vsiten.a;
        if (calculatePosition == VSiteCalculatePosition::Yes)
        {
            if (pbc)
            {
                pbc_dx_aiuc(pbc, x[ai], x1, dx);
            }
            else
            {
                rvec_sub(x[ai], x1, dx);
            }
            dsum[XX] += a * dx[XX];
            dsum[YY] += a * dx[YY];
            dsum[ZZ] += a * dx[ZZ];
            /* 9 Flops */
        }
        if (calculateVelocity == VSiteCalculateVelocity::Yes)
        {
            rvec_sub(v[ai], v1, dx);
            dvsum[XX] += a * dx[XX];
            dvsum[YY] += a * dx[YY];
            dvsum[ZZ] += a * dx[ZZ];
            /* 9 Flops */
        }
    }

    if (calculatePosition == VSiteCalculatePosition::Yes)
    {
        x[av][XX] = x1[XX] + dsum[XX];
        x[av][YY] = x1[YY] + dsum[YY];
        x[av][ZZ] = x1[ZZ] + dsum[ZZ];
    }

    if (calculateVelocity == VSiteCalculateVelocity::Yes)
    {
        v[av][XX] = v1[XX] + dvsum[XX];
        v[av][YY] = v1[YY] + dvsum[YY];
        v[av][ZZ] = v1[ZZ] + dvsum[ZZ];
    }

    return n3;
}

#endif // DOXYGEN

//! PBC modes for vsite construction and spreading
enum class PbcMode
{
    all, //!< Apply normal, simple PBC for all vsites
    none //!< No PBC treatment needed
};

/*! \brief Returns the PBC mode based on the system PBC and vsite properties
 *
 * \param[in] pbcPtr  A pointer to a PBC struct or nullptr when no PBC treatment is required
 */
static PbcMode getPbcMode(const t_pbc* pbcPtr)
{
    if (pbcPtr == nullptr)
    {
        return PbcMode::none;
    }
    else
    {
        return PbcMode::all;
    }
}

/*! \brief Executes the vsite construction task for a single thread
 *
 * \tparam        operation  Whether we are calculating positions, velocities, or both
 * \param[in,out] x   Coordinates to construct vsites for
 * \param[in,out] v   Velocities are generated for virtual sites if `operation` requires it
 * \param[in]     ip  Interaction parameters for all interaction, only vsite parameters are used
 * \param[in]     ilist  The interaction lists, only vsites are usesd
 * \param[in]     pbc_null  PBC struct, used for PBC distance calculations when !=nullptr
 */
template<VSiteCalculatePosition calculatePosition, VSiteCalculateVelocity calculateVelocity>
static void construct_vsites_thread(ArrayRef<RVec>                  x,
                                    ArrayRef<RVec>                  v,
                                    ArrayRef<const t_iparams>       ip,
                                    ArrayRef<const InteractionList> ilist,
                                    const t_pbc*                    pbc_null)
{
    if (calculateVelocity == VSiteCalculateVelocity::Yes)
    {
        GMX_RELEASE_ASSERT(x.empty() || !v.empty(),
                           "Can't calculate velocities without access to velocity vector.");
    }

    // Work around clang bug (unfixed as of Feb 2021)
    // https://bugs.llvm.org/show_bug.cgi?id=35450
    // clang-format off
    CLANG_DIAGNOSTIC_IGNORE(-Wunused-lambda-capture)
    // clang-format on
    // getVOrNull returns a velocity rvec if we need it, nullptr otherwise.
    auto getVOrNull = [v](int idx) -> real* {
        if (calculateVelocity == VSiteCalculateVelocity::Yes)
        {
            return v[idx].as_vec();
        }
        else
        {
            return nullptr;
        }
    };
    CLANG_DIAGNOSTIC_RESET

    const PbcMode pbcMode = getPbcMode(pbc_null);
    /* We need another pbc pointer, as with charge groups we switch per vsite */
    const t_pbc* pbc_null2 = pbc_null;

    for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
    {
        if (ilist[ftype].empty())
        {
            continue;
        }

        { // TODO remove me
            int nra = interaction_function[ftype].nratoms;
            int inc = 1 + nra;
            int nr  = ilist[ftype].size();

            const t_iatom* ia = ilist[ftype].iatoms.data();

            for (int i = 0; i < nr;)
            {
                int tp = ia[0];
                /* The vsite and constructing atoms */
                int avsite = ia[1];
                int ai     = ia[2];
                /* Constants for constructing vsites */
                real a1 = ip[tp].vsite.a;
                /* Copy the old position */
                rvec xv;
                copy_rvec(x[avsite], xv);

                /* Construct the vsite depending on type */
                int  aj, ak, al;
                real b1, c1;
                switch (ftype)
                {
                    case F_VSITE1:
                        constr_vsite1<calculatePosition, calculateVelocity>(
                                x[ai], x[avsite], getVOrNull(ai), getVOrNull(avsite));
                        break;
                    case F_VSITE2:
                        aj = ia[3];
                        constr_vsite2<calculatePosition, calculateVelocity>(x[ai],
                                                                            x[aj],
                                                                            x[avsite],
                                                                            a1,
                                                                            pbc_null2,
                                                                            getVOrNull(ai),
                                                                            getVOrNull(aj),
                                                                            getVOrNull(avsite));
                        break;
                    case F_VSITE2FD:
                        aj = ia[3];
                        constr_vsite2FD<calculatePosition, calculateVelocity>(x[ai],
                                                                              x[aj],
                                                                              x[avsite],
                                                                              a1,
                                                                              pbc_null2,
                                                                              getVOrNull(ai),
                                                                              getVOrNull(aj),
                                                                              getVOrNull(avsite));
                        break;
                    case F_VSITE3:
                        aj = ia[3];
                        ak = ia[4];
                        b1 = ip[tp].vsite.b;
                        constr_vsite3<calculatePosition, calculateVelocity>(x[ai],
                                                                            x[aj],
                                                                            x[ak],
                                                                            x[avsite],
                                                                            a1,
                                                                            b1,
                                                                            pbc_null2,
                                                                            getVOrNull(ai),
                                                                            getVOrNull(aj),
                                                                            getVOrNull(ak),
                                                                            getVOrNull(avsite));
                        break;
                    case F_VSITE3FD:
                        aj = ia[3];
                        ak = ia[4];
                        b1 = ip[tp].vsite.b;
                        constr_vsite3FD<calculatePosition, calculateVelocity>(x[ai],
                                                                              x[aj],
                                                                              x[ak],
                                                                              x[avsite],
                                                                              a1,
                                                                              b1,
                                                                              pbc_null2,
                                                                              getVOrNull(ai),
                                                                              getVOrNull(aj),
                                                                              getVOrNull(ak),
                                                                              getVOrNull(avsite));
                        break;
                    case F_VSITE3FAD:
                        aj = ia[3];
                        ak = ia[4];
                        b1 = ip[tp].vsite.b;
                        constr_vsite3FAD<calculatePosition, calculateVelocity>(x[ai],
                                                                               x[aj],
                                                                               x[ak],
                                                                               x[avsite],
                                                                               a1,
                                                                               b1,
                                                                               pbc_null2,
                                                                               getVOrNull(ai),
                                                                               getVOrNull(aj),
                                                                               getVOrNull(ak),
                                                                               getVOrNull(avsite));
                        break;
                    case F_VSITE3OUT:
                        aj = ia[3];
                        ak = ia[4];
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        constr_vsite3OUT<calculatePosition, calculateVelocity>(x[ai],
                                                                               x[aj],
                                                                               x[ak],
                                                                               x[avsite],
                                                                               a1,
                                                                               b1,
                                                                               c1,
                                                                               pbc_null2,
                                                                               getVOrNull(ai),
                                                                               getVOrNull(aj),
                                                                               getVOrNull(ak),
                                                                               getVOrNull(avsite));
                        break;
                    case F_VSITE4FD:
                        aj = ia[3];
                        ak = ia[4];
                        al = ia[5];
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        constr_vsite4FD<calculatePosition, calculateVelocity>(x[ai],
                                                                              x[aj],
                                                                              x[ak],
                                                                              x[al],
                                                                              x[avsite],
                                                                              a1,
                                                                              b1,
                                                                              c1,
                                                                              pbc_null2,
                                                                              getVOrNull(ai),
                                                                              getVOrNull(aj),
                                                                              getVOrNull(ak),
                                                                              getVOrNull(al),
                                                                              getVOrNull(avsite));
                        break;
                    case F_VSITE4FDN:
                        aj = ia[3];
                        ak = ia[4];
                        al = ia[5];
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        constr_vsite4FDN<calculatePosition, calculateVelocity>(x[ai],
                                                                               x[aj],
                                                                               x[ak],
                                                                               x[al],
                                                                               x[avsite],
                                                                               a1,
                                                                               b1,
                                                                               c1,
                                                                               pbc_null2,
                                                                               getVOrNull(ai),
                                                                               getVOrNull(aj),
                                                                               getVOrNull(ak),
                                                                               getVOrNull(al),
                                                                               getVOrNull(avsite));
                        break;
                    case F_VSITEN:
                        inc = constr_vsiten<calculatePosition, calculateVelocity>(ia, ip, x, pbc_null2, v);
                        break;
                    default:
                        gmx_fatal(FARGS, "No such vsite type %d in %s, line %d", ftype, __FILE__, __LINE__);
                }

                if (pbcMode == PbcMode::all)
                {
                    /* Keep the vsite in the same periodic image as before */
                    rvec dx;
                    int  ishift = pbc_dx_aiuc(pbc_null, x[avsite], xv, dx);
                    if (ishift != c_centralShiftIndex)
                    {
                        rvec_add(xv, dx, x[avsite]);
                    }
                }

                /* Increment loop variables */
                i += inc;
                ia += inc;
            }
        }
    }
}

/*! \brief Dispatch the vsite construction tasks for all threads
 *
 * \param[in]     threadingInfo  Used to divide work over threads when != nullptr
 * \param[in,out] x   Coordinates to construct vsites for
 * \param[in,out] v   When not empty, velocities are generated for virtual sites
 * \param[in]     ip  Interaction parameters for all interaction, only vsite parameters are used
 * \param[in]     ilist  The interaction lists, only vsites are usesd
 * \param[in]     domainInfo  Information about PBC and DD
 * \param[in]     box  Used for PBC when PBC is set in domainInfo
 */
template<VSiteCalculatePosition calculatePosition, VSiteCalculateVelocity calculateVelocity>
static void construct_vsites(const ThreadingInfo*            threadingInfo,
                             ArrayRef<RVec>                  x,
                             ArrayRef<RVec>                  v,
                             ArrayRef<const t_iparams>       ip,
                             ArrayRef<const InteractionList> ilist,
                             const DomainInfo&               domainInfo,
                             const matrix                    box)
{
    const bool useDomdec = domainInfo.useDomdec();

    t_pbc pbc, *pbc_null;

    /* We only need to do pbc when we have inter update-group vsites.
     * Note that with domain decomposition we do not need to apply PBC here
     * when we have at least 3 domains along each dimension. Currently we
     * do not optimize this case.
     */
    if (domainInfo.pbcType_ != PbcType::No && domainInfo.useMolPbc_)
    {
        /* This is wasting some CPU time as we now do this multiple times
         * per MD step.
         */
        ivec null_ivec;
        clear_ivec(null_ivec);
        pbc_null = set_pbc_dd(
                &pbc, domainInfo.pbcType_, useDomdec ? domainInfo.domdec_->numCells : null_ivec, FALSE, box);
    }
    else
    {
        pbc_null = nullptr;
    }

    if (useDomdec)
    {
        if (calculateVelocity == VSiteCalculateVelocity::Yes)
        {
            dd_move_x_and_v_vsites(*domainInfo.domdec_, box, x, v);
        }
        else
        {
            dd_move_x_vsites(*domainInfo.domdec_, box, x);
        }
    }

    if (threadingInfo == nullptr || threadingInfo->numThreads() == 1)
    {
        construct_vsites_thread<calculatePosition, calculateVelocity>(x, v, ip, ilist, pbc_null);
    }
    else
    {
#pragma omp parallel num_threads(threadingInfo->numThreads())
        {
            try
            {
                const int          th    = gmx_omp_get_thread_num();
                const VsiteThread& tData = threadingInfo->threadData(th);
                GMX_ASSERT(tData.rangeStart >= 0,
                           "The thread data should be initialized before calling construct_vsites");

                construct_vsites_thread<calculatePosition, calculateVelocity>(
                        x, v, ip, tData.ilist, pbc_null);
                if (tData.useInterdependentTask)
                {
                    /* Here we don't need a barrier (unlike the spreading),
                     * since both tasks only construct vsites from particles,
                     * or local vsites, not from non-local vsites.
                     */
                    construct_vsites_thread<calculatePosition, calculateVelocity>(
                            x, v, ip, tData.idTask.ilist, pbc_null);
                }
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }
        /* Now we can construct the vsites that might depend on other vsites */
        construct_vsites_thread<calculatePosition, calculateVelocity>(
                x, v, ip, threadingInfo->threadDataNonLocalDependent().ilist, pbc_null);
    }
}

void VirtualSitesHandler::Impl::construct(ArrayRef<RVec> x,
                                          ArrayRef<RVec> v,
                                          const matrix   box,
                                          VSiteOperation operation) const
{
    switch (operation)
    {
        case VSiteOperation::Positions:
            construct_vsites<VSiteCalculatePosition::Yes, VSiteCalculateVelocity::No>(
                    &threadingInfo_, x, v, iparams_, ilists_, domainInfo_, box);
            break;
        case VSiteOperation::Velocities:
            construct_vsites<VSiteCalculatePosition::No, VSiteCalculateVelocity::Yes>(
                    &threadingInfo_, x, v, iparams_, ilists_, domainInfo_, box);
            break;
        case VSiteOperation::PositionsAndVelocities:
            construct_vsites<VSiteCalculatePosition::Yes, VSiteCalculateVelocity::Yes>(
                    &threadingInfo_, x, v, iparams_, ilists_, domainInfo_, box);
            break;
        default: gmx_fatal(FARGS, "Unknown virtual site operation");
    }
}

void VirtualSitesHandler::construct(ArrayRef<RVec> x, ArrayRef<RVec> v, const matrix box, VSiteOperation operation) const
{
    impl_->construct(x, v, box, operation);
}

void constructVirtualSites(ArrayRef<RVec> x, ArrayRef<const t_iparams> ip, ArrayRef<const InteractionList> ilist)

{
    // No PBC, no DD
    const DomainInfo domainInfo;
    construct_vsites<VSiteCalculatePosition::Yes, VSiteCalculateVelocity::No>(
            nullptr, x, {}, ip, ilist, domainInfo, nullptr);
}

#ifndef DOXYGEN
/* Force spreading routines */

static void spread_vsite1(const t_iatom ia[], ArrayRef<RVec> f)
{
    const int av = ia[1];
    const int ai = ia[2];

    f[ai] += f[av];
}

template<VirialHandling virialHandling>
static void spread_vsite2(const t_iatom        ia[],
                          real                 a,
                          ArrayRef<const RVec> x,
                          ArrayRef<RVec>       f,
                          ArrayRef<RVec>       fshift,
                          const t_pbc*         pbc)
{
    rvec    fi, fj, dx;
    t_iatom av, ai, aj;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];

    svmul(1 - a, f[av], fi);
    svmul(a, f[av], fj);
    /* 7 flop */

    rvec_inc(f[ai], fi);
    rvec_inc(f[aj], fj);
    /* 6 Flops */

    if (virialHandling == VirialHandling::Pbc)
    {
        int siv;
        int sij;
        if (pbc)
        {
            siv = pbc_dx_aiuc(pbc, x[ai], x[av], dx);
            sij = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        }
        else
        {
            siv = c_centralShiftIndex;
            sij = c_centralShiftIndex;
        }

        if (siv != c_centralShiftIndex || sij != c_centralShiftIndex)
        {
            rvec_inc(fshift[siv], f[av]);
            rvec_dec(fshift[c_centralShiftIndex], fi);
            rvec_dec(fshift[sij], fj);
        }
    }

    /* TOTAL: 13 flops */
}

void constructVirtualSitesGlobal(const gmx_mtop_t& mtop, gmx::ArrayRef<gmx::RVec> x)
{
    GMX_ASSERT(x.ssize() >= mtop.natoms, "x should contain the whole system");
    GMX_ASSERT(!mtop.moleculeBlockIndices.empty(),
               "molblock indices are needed in constructVsitesGlobal");

    for (size_t mb = 0; mb < mtop.molblock.size(); mb++)
    {
        const gmx_molblock_t& molb = mtop.molblock[mb];
        const gmx_moltype_t&  molt = mtop.moltype[molb.type];
        if (vsiteIlistNrCount(molt.ilist) > 0)
        {
            int atomOffset = mtop.moleculeBlockIndices[mb].globalAtomStart;
            for (int mol = 0; mol < molb.nmol; mol++)
            {
                constructVirtualSites(
                        x.subArray(atomOffset, molt.atoms.nr), mtop.ffparams.iparams, molt.ilist);
                atomOffset += molt.atoms.nr;
            }
        }
    }
}

template<VirialHandling virialHandling>
static void spread_vsite2FD(const t_iatom        ia[],
                            real                 a,
                            ArrayRef<const RVec> x,
                            ArrayRef<RVec>       f,
                            ArrayRef<RVec>       fshift,
                            matrix               dxdf,
                            const t_pbc*         pbc)
{
    const int av = ia[1];
    const int ai = ia[2];
    const int aj = ia[3];
    rvec      fv;
    copy_rvec(f[av], fv);

    rvec xij;
    int  sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    /* 6 flops */

    const real invDistance = inverseNorm(xij);
    const real b           = a * invDistance;
    /* 4 + ?10? flops */

    const real fproj = iprod(xij, fv) * invDistance * invDistance;

    rvec fj;
    fj[XX] = b * (fv[XX] - fproj * xij[XX]);
    fj[YY] = b * (fv[YY] - fproj * xij[YY]);
    fj[ZZ] = b * (fv[ZZ] - fproj * xij[ZZ]);
    /* 9 */

    /* b is already calculated in constr_vsite2FD
       storing b somewhere will save flops.     */

    f[ai][XX] += fv[XX] - fj[XX];
    f[ai][YY] += fv[YY] - fj[YY];
    f[ai][ZZ] += fv[ZZ] - fj[ZZ];
    f[aj][XX] += fj[XX];
    f[aj][YY] += fj[YY];
    f[aj][ZZ] += fj[ZZ];
    /* 9 Flops */

    if (virialHandling == VirialHandling::Pbc)
    {
        int svi;
        if (pbc)
        {
            rvec xvi;
            svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
        }
        else
        {
            svi = c_centralShiftIndex;
        }

        if (svi != c_centralShiftIndex || sji != c_centralShiftIndex)
        {
            rvec_dec(fshift[svi], fv);
            fshift[c_centralShiftIndex][XX] += fv[XX] - fj[XX];
            fshift[c_centralShiftIndex][YY] += fv[YY] - fj[YY];
            fshift[c_centralShiftIndex][ZZ] += fv[ZZ] - fj[ZZ];
            fshift[sji][XX] += fj[XX];
            fshift[sji][YY] += fj[YY];
            fshift[sji][ZZ] += fj[ZZ];
        }
    }

    if (virialHandling == VirialHandling::NonLinear)
    {
        /* Under this condition, the virial for the current forces is not
         * calculated from the redistributed forces. This means that
         * the effect of non-linear virtual site constructions on the virial
         * needs to be added separately. This contribution can be calculated
         * in many ways, but the simplest and cheapest way is to use
         * the first constructing atom ai as a reference position in space:
         * subtract (xv-xi)*fv and add (xj-xi)*fj.
         */
        rvec xiv;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                /* As xix is a linear combination of j and k, use that here */
                dxdf[i][j] += -xiv[i] * fv[j] + xij[i] * fj[j];
            }
        }
    }

    /* TOTAL: 38 flops */
}

template<VirialHandling virialHandling>
static void spread_vsite3(const t_iatom        ia[],
                          real                 a,
                          real                 b,
                          ArrayRef<const RVec> x,
                          ArrayRef<RVec>       f,
                          ArrayRef<RVec>       fshift,
                          const t_pbc*         pbc)
{
    rvec fi, fj, fk, dx;
    int  av, ai, aj, ak;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];

    svmul(1 - a - b, f[av], fi);
    svmul(a, f[av], fj);
    svmul(b, f[av], fk);
    /* 11 flops */

    rvec_inc(f[ai], fi);
    rvec_inc(f[aj], fj);
    rvec_inc(f[ak], fk);
    /* 9 Flops */

    if (virialHandling == VirialHandling::Pbc)
    {
        int siv;
        int sij;
        int sik;
        if (pbc)
        {
            siv = pbc_dx_aiuc(pbc, x[ai], x[av], dx);
            sij = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
            sik = pbc_dx_aiuc(pbc, x[ai], x[ak], dx);
        }
        else
        {
            siv = c_centralShiftIndex;
            sij = c_centralShiftIndex;
            sik = c_centralShiftIndex;
        }

        if (siv != c_centralShiftIndex || sij != c_centralShiftIndex || sik != c_centralShiftIndex)
        {
            rvec_inc(fshift[siv], f[av]);
            rvec_dec(fshift[c_centralShiftIndex], fi);
            rvec_dec(fshift[sij], fj);
            rvec_dec(fshift[sik], fk);
        }
    }

    /* TOTAL: 20 flops */
}

template<VirialHandling virialHandling>
static void spread_vsite3FD(const t_iatom        ia[],
                            real                 a,
                            real                 b,
                            ArrayRef<const RVec> x,
                            ArrayRef<RVec>       f,
                            ArrayRef<RVec>       fshift,
                            matrix               dxdf,
                            const t_pbc*         pbc)
{
    real    fproj, a1;
    rvec    xvi, xij, xjk, xix, fv, temp;
    t_iatom av, ai, aj, ak;
    int     sji, skj;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];
    copy_rvec(f[av], fv);

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    skj = pbc_rvec_sub(pbc, x[ak], x[aj], xjk);
    /* 6 flops */

    /* xix goes from i to point x on the line jk */
    xix[XX] = xij[XX] + a * xjk[XX];
    xix[YY] = xij[YY] + a * xjk[YY];
    xix[ZZ] = xij[ZZ] + a * xjk[ZZ];
    /* 6 flops */

    const real invDistance = inverseNorm(xix);
    const real c           = b * invDistance;
    /* 4 + ?10? flops */

    fproj = iprod(xix, fv) * invDistance * invDistance; /* = (xix . f)/(xix . xix) */

    temp[XX] = c * (fv[XX] - fproj * xix[XX]);
    temp[YY] = c * (fv[YY] - fproj * xix[YY]);
    temp[ZZ] = c * (fv[ZZ] - fproj * xix[ZZ]);
    /* 16 */

    /* c is already calculated in constr_vsite3FD
       storing c somewhere will save 26 flops!     */

    a1 = 1 - a;
    f[ai][XX] += fv[XX] - temp[XX];
    f[ai][YY] += fv[YY] - temp[YY];
    f[ai][ZZ] += fv[ZZ] - temp[ZZ];
    f[aj][XX] += a1 * temp[XX];
    f[aj][YY] += a1 * temp[YY];
    f[aj][ZZ] += a1 * temp[ZZ];
    f[ak][XX] += a * temp[XX];
    f[ak][YY] += a * temp[YY];
    f[ak][ZZ] += a * temp[ZZ];
    /* 19 Flops */

    if (virialHandling == VirialHandling::Pbc)
    {
        int svi;
        if (pbc)
        {
            svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
        }
        else
        {
            svi = c_centralShiftIndex;
        }

        if (svi != c_centralShiftIndex || sji != c_centralShiftIndex || skj != c_centralShiftIndex)
        {
            rvec_dec(fshift[svi], fv);
            fshift[c_centralShiftIndex][XX] += fv[XX] - (1 + a) * temp[XX];
            fshift[c_centralShiftIndex][YY] += fv[YY] - (1 + a) * temp[YY];
            fshift[c_centralShiftIndex][ZZ] += fv[ZZ] - (1 + a) * temp[ZZ];
            fshift[sji][XX] += temp[XX];
            fshift[sji][YY] += temp[YY];
            fshift[sji][ZZ] += temp[ZZ];
            fshift[skj][XX] += a * temp[XX];
            fshift[skj][YY] += a * temp[YY];
            fshift[skj][ZZ] += a * temp[ZZ];
        }
    }

    if (virialHandling == VirialHandling::NonLinear)
    {
        /* Under this condition, the virial for the current forces is not
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
                dxdf[i][j] += -xiv[i] * fv[j] + xix[i] * temp[j];
            }
        }
    }

    /* TOTAL: 61 flops */
}

template<VirialHandling virialHandling>
static void spread_vsite3FAD(const t_iatom        ia[],
                             real                 a,
                             real                 b,
                             ArrayRef<const RVec> x,
                             ArrayRef<RVec>       f,
                             ArrayRef<RVec>       fshift,
                             matrix               dxdf,
                             const t_pbc*         pbc)
{
    rvec    xvi, xij, xjk, xperp, Fpij, Fppp, fv, f1, f2, f3;
    real    a1, b1, c1, c2, invdij, invdij2, invdp, fproj;
    t_iatom av, ai, aj, ak;
    int     sji, skj;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];
    copy_rvec(f[ia[1]], fv);

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    skj = pbc_rvec_sub(pbc, x[ak], x[aj], xjk);
    /* 6 flops */

    invdij    = inverseNorm(xij);
    invdij2   = invdij * invdij;
    c1        = iprod(xij, xjk) * invdij2;
    xperp[XX] = xjk[XX] - c1 * xij[XX];
    xperp[YY] = xjk[YY] - c1 * xij[YY];
    xperp[ZZ] = xjk[ZZ] - c1 * xij[ZZ];
    /* xperp in plane ijk, perp. to ij */
    invdp = inverseNorm(xperp);
    a1    = a * invdij;
    b1    = b * invdp;
    /* 45 flops */

    /* a1, b1 and c1 are already calculated in constr_vsite3FAD
       storing them somewhere will save 45 flops!     */

    fproj = iprod(xij, fv) * invdij2;
    svmul(fproj, xij, Fpij);                              /* proj. f on xij */
    svmul(iprod(xperp, fv) * invdp * invdp, xperp, Fppp); /* proj. f on xperp */
    svmul(b1 * fproj, xperp, f3);
    /* 23 flops */

    rvec_sub(fv, Fpij, f1); /* f1 = f - Fpij */
    rvec_sub(f1, Fppp, f2); /* f2 = f - Fpij - Fppp */
    for (int d = 0; d < DIM; d++)
    {
        f1[d] *= a1;
        f2[d] *= b1;
    }
    /* 12 flops */

    c2 = 1 + c1;
    f[ai][XX] += fv[XX] - f1[XX] + c1 * f2[XX] + f3[XX];
    f[ai][YY] += fv[YY] - f1[YY] + c1 * f2[YY] + f3[YY];
    f[ai][ZZ] += fv[ZZ] - f1[ZZ] + c1 * f2[ZZ] + f3[ZZ];
    f[aj][XX] += f1[XX] - c2 * f2[XX] - f3[XX];
    f[aj][YY] += f1[YY] - c2 * f2[YY] - f3[YY];
    f[aj][ZZ] += f1[ZZ] - c2 * f2[ZZ] - f3[ZZ];
    f[ak][XX] += f2[XX];
    f[ak][YY] += f2[YY];
    f[ak][ZZ] += f2[ZZ];
    /* 30 Flops */

    if (virialHandling == VirialHandling::Pbc)
    {
        int svi;

        if (pbc)
        {
            svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
        }
        else
        {
            svi = c_centralShiftIndex;
        }

        if (svi != c_centralShiftIndex || sji != c_centralShiftIndex || skj != c_centralShiftIndex)
        {
            rvec_dec(fshift[svi], fv);
            fshift[c_centralShiftIndex][XX] += fv[XX] - f1[XX] - (1 - c1) * f2[XX] + f3[XX];
            fshift[c_centralShiftIndex][YY] += fv[YY] - f1[YY] - (1 - c1) * f2[YY] + f3[YY];
            fshift[c_centralShiftIndex][ZZ] += fv[ZZ] - f1[ZZ] - (1 - c1) * f2[ZZ] + f3[ZZ];
            fshift[sji][XX] += f1[XX] - c1 * f2[XX] - f3[XX];
            fshift[sji][YY] += f1[YY] - c1 * f2[YY] - f3[YY];
            fshift[sji][ZZ] += f1[ZZ] - c1 * f2[ZZ] - f3[ZZ];
            fshift[skj][XX] += f2[XX];
            fshift[skj][YY] += f2[YY];
            fshift[skj][ZZ] += f2[ZZ];
        }
    }

    if (virialHandling == VirialHandling::NonLinear)
    {
        rvec xiv;
        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                /* Note that xik=xij+xjk, so we have to add xij*f2 */
                dxdf[i][j] += -xiv[i] * fv[j] + xij[i] * (f1[j] + (1 - c2) * f2[j] - f3[j])
                              + xjk[i] * f2[j];
            }
        }
    }

    /* TOTAL: 113 flops */
}

template<VirialHandling virialHandling>
static void spread_vsite3OUT(const t_iatom        ia[],
                             real                 a,
                             real                 b,
                             real                 c,
                             ArrayRef<const RVec> x,
                             ArrayRef<RVec>       f,
                             ArrayRef<RVec>       fshift,
                             matrix               dxdf,
                             const t_pbc*         pbc)
{
    rvec xvi, xij, xik, fv, fj, fk;
    real cfx, cfy, cfz;
    int  av, ai, aj, ak;
    int  sji, ski;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    ski = pbc_rvec_sub(pbc, x[ak], x[ai], xik);
    /* 6 Flops */

    copy_rvec(f[av], fv);

    cfx = c * fv[XX];
    cfy = c * fv[YY];
    cfz = c * fv[ZZ];
    /* 3 Flops */

    fj[XX] = a * fv[XX] - xik[ZZ] * cfy + xik[YY] * cfz;
    fj[YY] = xik[ZZ] * cfx + a * fv[YY] - xik[XX] * cfz;
    fj[ZZ] = -xik[YY] * cfx + xik[XX] * cfy + a * fv[ZZ];

    fk[XX] = b * fv[XX] + xij[ZZ] * cfy - xij[YY] * cfz;
    fk[YY] = -xij[ZZ] * cfx + b * fv[YY] + xij[XX] * cfz;
    fk[ZZ] = xij[YY] * cfx - xij[XX] * cfy + b * fv[ZZ];
    /* 30 Flops */

    f[ai][XX] += fv[XX] - fj[XX] - fk[XX];
    f[ai][YY] += fv[YY] - fj[YY] - fk[YY];
    f[ai][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ];
    rvec_inc(f[aj], fj);
    rvec_inc(f[ak], fk);
    /* 15 Flops */

    if (virialHandling == VirialHandling::Pbc)
    {
        int svi;
        if (pbc)
        {
            svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
        }
        else
        {
            svi = c_centralShiftIndex;
        }

        if (svi != c_centralShiftIndex || sji != c_centralShiftIndex || ski != c_centralShiftIndex)
        {
            rvec_dec(fshift[svi], fv);
            fshift[c_centralShiftIndex][XX] += fv[XX] - fj[XX] - fk[XX];
            fshift[c_centralShiftIndex][YY] += fv[YY] - fj[YY] - fk[YY];
            fshift[c_centralShiftIndex][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ];
            rvec_inc(fshift[sji], fj);
            rvec_inc(fshift[ski], fk);
        }
    }

    if (virialHandling == VirialHandling::NonLinear)
    {
        rvec xiv;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                dxdf[i][j] += -xiv[i] * fv[j] + xij[i] * fj[j] + xik[i] * fk[j];
            }
        }
    }

    /* TOTAL: 54 flops */
}

template<VirialHandling virialHandling>
static void spread_vsite4FD(const t_iatom        ia[],
                            real                 a,
                            real                 b,
                            real                 c,
                            ArrayRef<const RVec> x,
                            ArrayRef<RVec>       f,
                            ArrayRef<RVec>       fshift,
                            matrix               dxdf,
                            const t_pbc*         pbc)
{
    real fproj, a1;
    rvec xvi, xij, xjk, xjl, xix, fv, temp;
    int  av, ai, aj, ak, al;
    int  sji, skj, slj, m;

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
        xix[m] = xij[m] + a * xjk[m] + b * xjl[m];
    }
    /* 12 flops */

    const real invDistance = inverseNorm(xix);
    const real d           = c * invDistance;
    /* 4 + ?10? flops */

    copy_rvec(f[av], fv);

    fproj = iprod(xix, fv) * invDistance * invDistance; /* = (xix . f)/(xix . xix) */

    for (m = 0; m < DIM; m++)
    {
        temp[m] = d * (fv[m] - fproj * xix[m]);
    }
    /* 16 */

    /* c is already calculated in constr_vsite3FD
       storing c somewhere will save 35 flops!     */

    a1 = 1 - a - b;
    for (m = 0; m < DIM; m++)
    {
        f[ai][m] += fv[m] - temp[m];
        f[aj][m] += a1 * temp[m];
        f[ak][m] += a * temp[m];
        f[al][m] += b * temp[m];
    }
    /* 26 Flops */

    if (virialHandling == VirialHandling::Pbc)
    {
        int svi;
        if (pbc)
        {
            svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
        }
        else
        {
            svi = c_centralShiftIndex;
        }

        if (svi != c_centralShiftIndex || sji != c_centralShiftIndex || skj != c_centralShiftIndex
            || slj != c_centralShiftIndex)
        {
            rvec_dec(fshift[svi], fv);
            for (m = 0; m < DIM; m++)
            {
                fshift[c_centralShiftIndex][m] += fv[m] - (1 + a + b) * temp[m];
                fshift[sji][m] += temp[m];
                fshift[skj][m] += a * temp[m];
                fshift[slj][m] += b * temp[m];
            }
        }
    }

    if (virialHandling == VirialHandling::NonLinear)
    {
        rvec xiv;
        int  i, j;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                dxdf[i][j] += -xiv[i] * fv[j] + xix[i] * temp[j];
            }
        }
    }

    /* TOTAL: 77 flops */
}

template<VirialHandling virialHandling>
static void spread_vsite4FDN(const t_iatom        ia[],
                             real                 a,
                             real                 b,
                             real                 c,
                             ArrayRef<const RVec> x,
                             ArrayRef<RVec>       f,
                             ArrayRef<RVec>       fshift,
                             matrix               dxdf,
                             const t_pbc*         pbc)
{
    rvec xvi, xij, xik, xil, ra, rb, rja, rjb, rab, rm, rt;
    rvec fv, fj, fk, fl;
    real invrm, denom;
    real cfx, cfy, cfz;
    int  av, ai, aj, ak, al;
    int  sij, sik, sil;

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

    ra[XX] = a * xik[XX];
    ra[YY] = a * xik[YY];
    ra[ZZ] = a * xik[ZZ];

    rb[XX] = b * xil[XX];
    rb[YY] = b * xil[YY];
    rb[ZZ] = b * xil[ZZ];

    /* 6 flops */

    rvec_sub(ra, xij, rja);
    rvec_sub(rb, xij, rjb);
    rvec_sub(rb, ra, rab);
    /* 9 flops */

    cprod(rja, rjb, rm);
    /* 9 flops */

    invrm = inverseNorm(rm);
    denom = invrm * invrm;
    /* 5+5+2 flops */

    cfx = c * invrm * fv[XX];
    cfy = c * invrm * fv[YY];
    cfz = c * invrm * fv[ZZ];
    /* 6 Flops */

    cprod(rm, rab, rt);
    /* 9 flops */

    rt[XX] *= denom;
    rt[YY] *= denom;
    rt[ZZ] *= denom;
    /* 3flops */

    fj[XX] = (-rm[XX] * rt[XX]) * cfx + (rab[ZZ] - rm[YY] * rt[XX]) * cfy
             + (-rab[YY] - rm[ZZ] * rt[XX]) * cfz;
    fj[YY] = (-rab[ZZ] - rm[XX] * rt[YY]) * cfx + (-rm[YY] * rt[YY]) * cfy
             + (rab[XX] - rm[ZZ] * rt[YY]) * cfz;
    fj[ZZ] = (rab[YY] - rm[XX] * rt[ZZ]) * cfx + (-rab[XX] - rm[YY] * rt[ZZ]) * cfy
             + (-rm[ZZ] * rt[ZZ]) * cfz;
    /* 30 flops */

    cprod(rjb, rm, rt);
    /* 9 flops */

    rt[XX] *= denom * a;
    rt[YY] *= denom * a;
    rt[ZZ] *= denom * a;
    /* 3flops */

    fk[XX] = (-rm[XX] * rt[XX]) * cfx + (-a * rjb[ZZ] - rm[YY] * rt[XX]) * cfy
             + (a * rjb[YY] - rm[ZZ] * rt[XX]) * cfz;
    fk[YY] = (a * rjb[ZZ] - rm[XX] * rt[YY]) * cfx + (-rm[YY] * rt[YY]) * cfy
             + (-a * rjb[XX] - rm[ZZ] * rt[YY]) * cfz;
    fk[ZZ] = (-a * rjb[YY] - rm[XX] * rt[ZZ]) * cfx + (a * rjb[XX] - rm[YY] * rt[ZZ]) * cfy
             + (-rm[ZZ] * rt[ZZ]) * cfz;
    /* 36 flops */

    cprod(rm, rja, rt);
    /* 9 flops */

    rt[XX] *= denom * b;
    rt[YY] *= denom * b;
    rt[ZZ] *= denom * b;
    /* 3flops */

    fl[XX] = (-rm[XX] * rt[XX]) * cfx + (b * rja[ZZ] - rm[YY] * rt[XX]) * cfy
             + (-b * rja[YY] - rm[ZZ] * rt[XX]) * cfz;
    fl[YY] = (-b * rja[ZZ] - rm[XX] * rt[YY]) * cfx + (-rm[YY] * rt[YY]) * cfy
             + (b * rja[XX] - rm[ZZ] * rt[YY]) * cfz;
    fl[ZZ] = (b * rja[YY] - rm[XX] * rt[ZZ]) * cfx + (-b * rja[XX] - rm[YY] * rt[ZZ]) * cfy
             + (-rm[ZZ] * rt[ZZ]) * cfz;
    /* 36 flops */

    f[ai][XX] += fv[XX] - fj[XX] - fk[XX] - fl[XX];
    f[ai][YY] += fv[YY] - fj[YY] - fk[YY] - fl[YY];
    f[ai][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ] - fl[ZZ];
    rvec_inc(f[aj], fj);
    rvec_inc(f[ak], fk);
    rvec_inc(f[al], fl);
    /* 21 flops */

    if (virialHandling == VirialHandling::Pbc)
    {
        int svi;
        if (pbc)
        {
            svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
        }
        else
        {
            svi = c_centralShiftIndex;
        }

        if (svi != c_centralShiftIndex || sij != c_centralShiftIndex || sik != c_centralShiftIndex
            || sil != c_centralShiftIndex)
        {
            rvec_dec(fshift[svi], fv);
            fshift[c_centralShiftIndex][XX] += fv[XX] - fj[XX] - fk[XX] - fl[XX];
            fshift[c_centralShiftIndex][YY] += fv[YY] - fj[YY] - fk[YY] - fl[YY];
            fshift[c_centralShiftIndex][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ] - fl[ZZ];
            rvec_inc(fshift[sij], fj);
            rvec_inc(fshift[sik], fk);
            rvec_inc(fshift[sil], fl);
        }
    }

    if (virialHandling == VirialHandling::NonLinear)
    {
        rvec xiv;
        int  i, j;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                dxdf[i][j] += -xiv[i] * fv[j] + xij[i] * fj[j] + xik[i] * fk[j] + xil[i] * fl[j];
            }
        }
    }

    /* Total: 207 flops (Yuck!) */
}

template<VirialHandling virialHandling>
static int spread_vsiten(const t_iatom             ia[],
                         ArrayRef<const t_iparams> ip,
                         ArrayRef<const RVec>      x,
                         ArrayRef<RVec>            f,
                         ArrayRef<RVec>            fshift,
                         const t_pbc*              pbc)
{
    rvec xv, dx, fi;
    int  n3, av, i, ai;
    real a;
    int  siv;

    n3 = 3 * ip[ia[0]].vsiten.n;
    av = ia[1];
    copy_rvec(x[av], xv);

    for (i = 0; i < n3; i += 3)
    {
        ai = ia[i + 2];
        if (pbc)
        {
            siv = pbc_dx_aiuc(pbc, x[ai], xv, dx);
        }
        else
        {
            siv = c_centralShiftIndex;
        }
        a = ip[ia[i]].vsiten.a;
        svmul(a, f[av], fi);
        rvec_inc(f[ai], fi);

        if (virialHandling == VirialHandling::Pbc && siv != c_centralShiftIndex)
        {
            rvec_inc(fshift[siv], fi);
            rvec_dec(fshift[c_centralShiftIndex], fi);
        }
        /* 6 Flops */
    }

    return n3;
}

#endif // DOXYGEN

//! Returns the number of virtual sites in the interaction list, for VSITEN the number of atoms
static int vsite_count(ArrayRef<const InteractionList> ilist, int ftype)
{
    if (ftype == F_VSITEN)
    {
        return ilist[ftype].size() / 3;
    }
    else
    {
        return ilist[ftype].size() / (1 + interaction_function[ftype].nratoms);
    }
}

//! Executes the force spreading task for a single thread
template<VirialHandling virialHandling>
static void spreadForceForThread(ArrayRef<const RVec>            x,
                                 ArrayRef<RVec>                  f,
                                 ArrayRef<RVec>                  fshift,
                                 matrix                          dxdf,
                                 ArrayRef<const t_iparams>       ip,
                                 ArrayRef<const InteractionList> ilist,
                                 const t_pbc*                    pbc_null)
{
    const PbcMode pbcMode = getPbcMode(pbc_null);
    /* We need another pbc pointer, as with charge groups we switch per vsite */
    const t_pbc*             pbc_null2 = pbc_null;
    gmx::ArrayRef<const int> vsite_pbc;

    /* this loop goes backwards to be able to build *
     * higher type vsites from lower types         */
    for (int ftype = c_ftypeVsiteEnd - 1; ftype >= c_ftypeVsiteStart; ftype--)
    {
        if (ilist[ftype].empty())
        {
            continue;
        }

        { // TODO remove me
            int nra = interaction_function[ftype].nratoms;
            int inc = 1 + nra;
            int nr  = ilist[ftype].size();

            const t_iatom* ia = ilist[ftype].iatoms.data();

            if (pbcMode == PbcMode::all)
            {
                pbc_null2 = pbc_null;
            }

            for (int i = 0; i < nr;)
            {
                int tp = ia[0];

                /* Constants for constructing */
                real a1, b1, c1;
                a1 = ip[tp].vsite.a;
                /* Construct the vsite depending on type */
                switch (ftype)
                {
                    case F_VSITE1: spread_vsite1(ia, f); break;
                    case F_VSITE2:
                        spread_vsite2<virialHandling>(ia, a1, x, f, fshift, pbc_null2);
                        break;
                    case F_VSITE2FD:
                        spread_vsite2FD<virialHandling>(ia, a1, x, f, fshift, dxdf, pbc_null2);
                        break;
                    case F_VSITE3:
                        b1 = ip[tp].vsite.b;
                        spread_vsite3<virialHandling>(ia, a1, b1, x, f, fshift, pbc_null2);
                        break;
                    case F_VSITE3FD:
                        b1 = ip[tp].vsite.b;
                        spread_vsite3FD<virialHandling>(ia, a1, b1, x, f, fshift, dxdf, pbc_null2);
                        break;
                    case F_VSITE3FAD:
                        b1 = ip[tp].vsite.b;
                        spread_vsite3FAD<virialHandling>(ia, a1, b1, x, f, fshift, dxdf, pbc_null2);
                        break;
                    case F_VSITE3OUT:
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        spread_vsite3OUT<virialHandling>(ia, a1, b1, c1, x, f, fshift, dxdf, pbc_null2);
                        break;
                    case F_VSITE4FD:
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        spread_vsite4FD<virialHandling>(ia, a1, b1, c1, x, f, fshift, dxdf, pbc_null2);
                        break;
                    case F_VSITE4FDN:
                        b1 = ip[tp].vsite.b;
                        c1 = ip[tp].vsite.c;
                        spread_vsite4FDN<virialHandling>(ia, a1, b1, c1, x, f, fshift, dxdf, pbc_null2);
                        break;
                    case F_VSITEN:
                        inc = spread_vsiten<virialHandling>(ia, ip, x, f, fshift, pbc_null2);
                        break;
                    default:
                        gmx_fatal(FARGS, "No such vsite type %d in %s, line %d", ftype, __FILE__, __LINE__);
                }
                clear_rvec(f[ia[1]]);

                /* Increment loop variables */
                i += inc;
                ia += inc;
            }
        }
    }
}

//! Wrapper function for calling the templated thread-local spread function
static void spreadForceWrapper(ArrayRef<const RVec>            x,
                               ArrayRef<RVec>                  f,
                               const VirialHandling            virialHandling,
                               ArrayRef<RVec>                  fshift,
                               matrix                          dxdf,
                               const bool                      clearDxdf,
                               ArrayRef<const t_iparams>       ip,
                               ArrayRef<const InteractionList> ilist,
                               const t_pbc*                    pbc_null)
{
    if (virialHandling == VirialHandling::NonLinear && clearDxdf)
    {
        clear_mat(dxdf);
    }

    switch (virialHandling)
    {
        case VirialHandling::None:
            spreadForceForThread<VirialHandling::None>(x, f, fshift, dxdf, ip, ilist, pbc_null);
            break;
        case VirialHandling::Pbc:
            spreadForceForThread<VirialHandling::Pbc>(x, f, fshift, dxdf, ip, ilist, pbc_null);
            break;
        case VirialHandling::NonLinear:
            spreadForceForThread<VirialHandling::NonLinear>(x, f, fshift, dxdf, ip, ilist, pbc_null);
            break;
    }
}

//! Clears the task force buffer elements that are written by task idTask
static void clearTaskForceBufferUsedElements(InterdependentTask* idTask)
{
    int ntask = idTask->spreadTask.size();
    for (int ti = 0; ti < ntask; ti++)
    {
        const AtomIndex* atomList = &idTask->atomIndex[idTask->spreadTask[ti]];
        int              natom    = atomList->atom.size();
        RVec*            force    = idTask->force.data();
        for (int i = 0; i < natom; i++)
        {
            clear_rvec(force[atomList->atom[i]]);
        }
    }
}

void VirtualSitesHandler::Impl::spreadForces(ArrayRef<const RVec> x,
                                             ArrayRef<RVec>       f,
                                             const VirialHandling virialHandling,
                                             ArrayRef<RVec>       fshift,
                                             matrix               virial,
                                             t_nrnb*              nrnb,
                                             const matrix         box,
                                             gmx_wallcycle*       wcycle)
{
    wallcycle_start(wcycle, WallCycleCounter::VsiteSpread);

    const bool useDomdec = domainInfo_.useDomdec();

    t_pbc pbc, *pbc_null;

    if (domainInfo_.useMolPbc_)
    {
        /* This is wasting some CPU time as we now do this multiple times
         * per MD step.
         */
        pbc_null = set_pbc_dd(
                &pbc, domainInfo_.pbcType_, useDomdec ? domainInfo_.domdec_->numCells : nullptr, FALSE, box);
    }
    else
    {
        pbc_null = nullptr;
    }

    if (useDomdec)
    {
        dd_clear_f_vsites(*domainInfo_.domdec_, f);
    }

    const int numThreads = threadingInfo_.numThreads();

    if (numThreads == 1)
    {
        matrix dxdf;
        spreadForceWrapper(x, f, virialHandling, fshift, dxdf, true, iparams_, ilists_, pbc_null);

        if (virialHandling == VirialHandling::NonLinear)
        {
            for (int i = 0; i < DIM; i++)
            {
                for (int j = 0; j < DIM; j++)
                {
                    virial[i][j] += -0.5 * dxdf[i][j];
                }
            }
        }
    }
    else
    {
        /* First spread the vsites that might depend on non-local vsites */
        auto& nlDependentVSites = threadingInfo_.threadDataNonLocalDependent();
        spreadForceWrapper(x,
                           f,
                           virialHandling,
                           fshift,
                           nlDependentVSites.dxdf,
                           true,
                           iparams_,
                           nlDependentVSites.ilist,
                           pbc_null);

#pragma omp parallel num_threads(numThreads)
        {
            try
            {
                int          thread = gmx_omp_get_thread_num();
                VsiteThread& tData  = threadingInfo_.threadData(thread);

                ArrayRef<RVec> fshift_t;
                if (virialHandling == VirialHandling::Pbc)
                {
                    if (thread == 0)
                    {
                        fshift_t = fshift;
                    }
                    else
                    {
                        fshift_t = tData.fshift;

                        for (int i = 0; i < c_numShiftVectors; i++)
                        {
                            clear_rvec(fshift_t[i]);
                        }
                    }
                }

                if (tData.useInterdependentTask)
                {
                    /* Spread the vsites that spread outside our local range.
                     * This is done using a thread-local force buffer force.
                     * First we need to copy the input vsite forces to force.
                     */
                    InterdependentTask* idTask = &tData.idTask;

                    /* Clear the buffer elements set by our task during
                     * the last call to spread_vsite_f.
                     */
                    clearTaskForceBufferUsedElements(idTask);

                    int nvsite = idTask->vsite.size();
                    for (int i = 0; i < nvsite; i++)
                    {
                        copy_rvec(f[idTask->vsite[i]], idTask->force[idTask->vsite[i]]);
                    }
                    spreadForceWrapper(x,
                                       idTask->force,
                                       virialHandling,
                                       fshift_t,
                                       tData.dxdf,
                                       true,
                                       iparams_,
                                       tData.idTask.ilist,
                                       pbc_null);

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
                        const InterdependentTask& idt_foreign =
                                threadingInfo_.threadData(idTask->reduceTask[ti]).idTask;
                        const AtomIndex& atomList  = idt_foreign.atomIndex[thread];
                        const RVec*      f_foreign = idt_foreign.force.data();

                        for (int ind : atomList.atom)
                        {
                            rvec_inc(f[ind], f_foreign[ind]);
                            /* Clearing of f_foreign is done at the next step */
                        }
                    }
                    /* Clear the vsite forces, both in f and force */
                    for (int i = 0; i < nvsite; i++)
                    {
                        int ind = tData.idTask.vsite[i];
                        clear_rvec(f[ind]);
                        clear_rvec(tData.idTask.force[ind]);
                    }
                }

                /* Spread the vsites that spread locally only */
                spreadForceWrapper(
                        x, f, virialHandling, fshift_t, tData.dxdf, false, iparams_, tData.ilist, pbc_null);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }

        if (virialHandling == VirialHandling::Pbc)
        {
            for (int th = 1; th < numThreads; th++)
            {
                for (int i = 0; i < c_numShiftVectors; i++)
                {
                    rvec_inc(fshift[i], threadingInfo_.threadData(th).fshift[i]);
                }
            }
        }

        if (virialHandling == VirialHandling::NonLinear)
        {
            for (int th = 0; th < numThreads + 1; th++)
            {
                /* MSVC doesn't like matrix references, so we use a pointer */
                const matrix& dxdf = threadingInfo_.threadData(th).dxdf;

                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        virial[i][j] += -0.5 * dxdf[i][j];
                    }
                }
            }
        }
    }

    if (useDomdec)
    {
        dd_move_f_vsites(*domainInfo_.domdec_, f, fshift);
    }

    inc_nrnb(nrnb, eNR_VSITE1, vsite_count(ilists_, F_VSITE1));
    inc_nrnb(nrnb, eNR_VSITE2, vsite_count(ilists_, F_VSITE2));
    inc_nrnb(nrnb, eNR_VSITE2FD, vsite_count(ilists_, F_VSITE2FD));
    inc_nrnb(nrnb, eNR_VSITE3, vsite_count(ilists_, F_VSITE3));
    inc_nrnb(nrnb, eNR_VSITE3FD, vsite_count(ilists_, F_VSITE3FD));
    inc_nrnb(nrnb, eNR_VSITE3FAD, vsite_count(ilists_, F_VSITE3FAD));
    inc_nrnb(nrnb, eNR_VSITE3OUT, vsite_count(ilists_, F_VSITE3OUT));
    inc_nrnb(nrnb, eNR_VSITE4FD, vsite_count(ilists_, F_VSITE4FD));
    inc_nrnb(nrnb, eNR_VSITE4FDN, vsite_count(ilists_, F_VSITE4FDN));
    inc_nrnb(nrnb, eNR_VSITEN, vsite_count(ilists_, F_VSITEN));

    wallcycle_stop(wcycle, WallCycleCounter::VsiteSpread);
}

/*! \brief Returns the an array with group indices for each atom
 *
 * \param[in] grouping  The partitioning of the atom range into atom groups
 */
static std::vector<int> makeAtomToGroupMapping(const gmx::RangePartitioning& grouping)
{
    std::vector<int> atomToGroup(grouping.fullRange().end(), 0);

    for (int group = 0; group < grouping.numBlocks(); group++)
    {
        auto block = grouping.block(group);
        std::fill(atomToGroup.begin() + block.begin(), atomToGroup.begin() + block.end(), group);
    }

    return atomToGroup;
}

int countNonlinearVsites(const gmx_mtop_t& mtop)
{
    int numNonlinearVsites = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        const gmx_moltype_t& molt = mtop.moltype[molb.type];

        for (const auto& ilist : extractILists(molt.ilist, IF_VSITE))
        {
            if (ilist.functionType != F_VSITE2 && ilist.functionType != F_VSITE3
                && ilist.functionType != F_VSITEN)
            {
                numNonlinearVsites += molb.nmol * ilist.iatoms.size() / (1 + NRAL(ilist.functionType));
            }
        }
    }

    return numNonlinearVsites;
}

void VirtualSitesHandler::spreadForces(ArrayRef<const RVec> x,
                                       ArrayRef<RVec>       f,
                                       const VirialHandling virialHandling,
                                       ArrayRef<RVec>       fshift,
                                       matrix               virial,
                                       t_nrnb*              nrnb,
                                       const matrix         box,
                                       gmx_wallcycle*       wcycle)
{
    impl_->spreadForces(x, f, virialHandling, fshift, virial, nrnb, box, wcycle);
}

int countInterUpdategroupVsites(const gmx_mtop_t&                           mtop,
                                gmx::ArrayRef<const gmx::RangePartitioning> updateGroupingsPerMoleculeType)
{
    int n_intercg_vsite = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        const gmx_moltype_t& molt = mtop.moltype[molb.type];

        std::vector<int> atomToGroup;
        if (!updateGroupingsPerMoleculeType.empty())
        {
            atomToGroup = makeAtomToGroupMapping(updateGroupingsPerMoleculeType[molb.type]);
        }
        for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
        {
            const int              nral = NRAL(ftype);
            const InteractionList& il   = molt.ilist[ftype];
            for (int i = 0; i < il.size(); i += 1 + nral)
            {
                bool isInterGroup = atomToGroup.empty();
                if (!isInterGroup)
                {
                    const int group = atomToGroup[il.iatoms[1 + i]];
                    for (int a = 1; a < nral; a++)
                    {
                        if (atomToGroup[il.iatoms[1 + a]] != group)
                        {
                            isInterGroup = true;
                            break;
                        }
                    }
                }
                if (isInterGroup)
                {
                    n_intercg_vsite += molb.nmol;
                }
            }
        }
    }

    return n_intercg_vsite;
}

std::unique_ptr<VirtualSitesHandler> makeVirtualSitesHandler(const gmx_mtop_t& mtop,
                                                             const t_commrec*  cr,
                                                             PbcType           pbcType,
                                                             ArrayRef<const RangePartitioning> updateGroupingPerMoleculeType)
{
    GMX_RELEASE_ASSERT(cr != nullptr, "We need a valid commrec");

    std::unique_ptr<VirtualSitesHandler> vsite;

    /* check if there are vsites */
    int nvsite = 0;
    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            GMX_ASSERT(ftype >= c_ftypeVsiteStart && ftype < c_ftypeVsiteEnd,
                       "c_ftypeVsiteStart and/or c_ftypeVsiteEnd do not have correct values");

            nvsite += gmx_mtop_ftype_count(mtop, ftype);
        }
        else
        {
            GMX_ASSERT(ftype < c_ftypeVsiteStart || ftype >= c_ftypeVsiteEnd,
                       "c_ftypeVsiteStart and/or c_ftypeVsiteEnd do not have correct values");
        }
    }

    if (nvsite == 0)
    {
        return vsite;
    }

    return std::make_unique<VirtualSitesHandler>(mtop, cr->dd, pbcType, updateGroupingPerMoleculeType);
}

ThreadingInfo::ThreadingInfo() : numThreads_(gmx_omp_nthreads_get(ModuleMultiThread::VirtualSite))
{
    if (numThreads_ > 1)
    {
        /* We need one extra thread data structure for the overlap vsites */
        tData_.resize(numThreads_ + 1);
#pragma omp parallel for num_threads(numThreads_) schedule(static)
        for (int thread = 0; thread < numThreads_; thread++)
        {
            try
            {
                tData_[thread] = std::make_unique<VsiteThread>();

                InterdependentTask& idTask = tData_[thread]->idTask;
                idTask.nuse                = 0;
                idTask.atomIndex.resize(numThreads_);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }
        if (numThreads_ > 1)
        {
            tData_[numThreads_] = std::make_unique<VsiteThread>();
        }
    }
}

VirtualSitesHandler::Impl::Impl(const gmx_mtop_t&                       mtop,
                                gmx_domdec_t*                           domdec,
                                const PbcType                           pbcType,
                                const ArrayRef<const RangePartitioning> updateGroupingPerMoleculeType) :
    numInterUpdategroupVirtualSites_(countInterUpdategroupVsites(mtop, updateGroupingPerMoleculeType)),
    domainInfo_({ pbcType, pbcType != PbcType::No && numInterUpdategroupVirtualSites_ > 0, domdec }),
    iparams_(mtop.ffparams.iparams)
{
}

VirtualSitesHandler::VirtualSitesHandler(const gmx_mtop_t&                       mtop,
                                         gmx_domdec_t*                           domdec,
                                         const PbcType                           pbcType,
                                         const ArrayRef<const RangePartitioning> updateGroupingPerMoleculeType) :
    impl_(new Impl(mtop, domdec, pbcType, updateGroupingPerMoleculeType))
{
}

//! Flag that atom \p atom which is home in another task, if it has not already been added before
static inline void flagAtom(InterdependentTask* idTask, const int atom, const int numThreads, const int numAtomsPerThread)
{
    if (!idTask->use[atom])
    {
        idTask->use[atom] = true;
        int thread        = atom / numAtomsPerThread;
        /* Assign all non-local atom force writes to thread 0 */
        if (thread >= numThreads)
        {
            thread = 0;
        }
        idTask->atomIndex[thread].atom.push_back(atom);
    }
}

/*! \brief Here we try to assign all vsites that are in our local range.
 *
 * Our task local atom range is tData->rangeStart - tData->rangeEnd.
 * Vsites that depend only on local atoms, as indicated by taskIndex[]==thread,
 * are assigned to task tData->ilist. Vsites that depend on non-local atoms
 * but not on other vsites are assigned to task tData->id_task.ilist.
 * taskIndex[] is set for all vsites in our range, either to our local tasks
 * or to the single last task as taskIndex[]=2*nthreads.
 */
static void assignVsitesToThread(VsiteThread*                    tData,
                                 int                             thread,
                                 int                             nthread,
                                 int                             natperthread,
                                 gmx::ArrayRef<int>              taskIndex,
                                 ArrayRef<const InteractionList> ilist,
                                 ArrayRef<const t_iparams>       ip,
                                 ArrayRef<const ParticleType>    ptype)
{
    for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
    {
        tData->ilist[ftype].clear();
        tData->idTask.ilist[ftype].clear();

        const int  nral1 = 1 + NRAL(ftype);
        const int* iat   = ilist[ftype].iatoms.data();
        for (int i = 0; i < ilist[ftype].size();)
        {
            /* Get the number of iatom entries in this virtual site.
             * The 3 below for F_VSITEN is from 1+NRAL(ftype)=3
             */
            const int numIAtoms = (ftype == F_VSITEN ? ip[iat[i]].vsiten.n * 3 : nral1);

            if (iat[1 + i] < tData->rangeStart || iat[1 + i] >= tData->rangeEnd)
            {
                /* This vsite belongs to a different thread */
                i += numIAtoms;
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
                    if (iat[j] < tData->rangeStart || iat[j] >= tData->rangeEnd || taskIndex[iat[j]] != thread)
                    {
                        if (!tData->useInterdependentTask || ptype[iat[j]] == ParticleType::VSite)
                        {
                            /* At least one constructing atom is a vsite
                             * that is not assigned to the same thread.
                             * Put this vsite into a separate task.
                             */
                            task = 2 * nthread;
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
                for (int j = i + 2; j < i + numIAtoms; j += 3)
                {
                    /* Do a range check to avoid a harmless race on taskIndex */
                    if (iat[j] < tData->rangeStart || iat[j] >= tData->rangeEnd || taskIndex[iat[j]] != thread)
                    {
                        GMX_ASSERT(ptype[iat[j]] != ParticleType::VSite,
                                   "A vsite to be assigned in assignVsitesToThread has a vsite as "
                                   "a constructing atom that does not belong to our task, such "
                                   "vsites should be assigned to the single 'main' task");

                        if (tData->useInterdependentTask)
                        {
                            // Assign to the interdependent task
                            task = nthread + thread;
                        }
                        else
                        {
                            // Assign to the separate, non-parallel task
                            task = 2 * nthread;
                        }
                    }
                }
            }

            /* Update this vsite's thread index entry */
            taskIndex[iat[1 + i]] = task;

            if (task == thread || task == nthread + thread)
            {
                /* Copy this vsite to the thread data struct of thread */
                InteractionList* il_task;
                if (task == thread)
                {
                    il_task = &tData->ilist[ftype];
                }
                else
                {
                    il_task = &tData->idTask.ilist[ftype];
                }
                /* Copy the vsite data to the thread-task local array */
                il_task->push_back(iat[i], numIAtoms - 1, iat + i + 1);
                if (task == nthread + thread)
                {
                    /* This vsite writes outside our own task force block.
                     * Put it into the interdependent task list and flag
                     * the atoms involved for reduction.
                     */
                    tData->idTask.vsite.push_back(iat[i + 1]);
                    if (ftype != F_VSITEN)
                    {
                        for (int j = i + 2; j < i + nral1; j++)
                        {
                            flagAtom(&tData->idTask, iat[j], nthread, natperthread);
                        }
                    }
                    else
                    {
                        for (int j = i + 2; j < i + numIAtoms; j += 3)
                        {
                            flagAtom(&tData->idTask, iat[j], nthread, natperthread);
                        }
                    }
                }
            }

            i += numIAtoms;
        }
    }
}

/*! \brief Assign all vsites with taskIndex[]==task to task tData */
static void assignVsitesToSingleTask(VsiteThread*                    tData,
                                     int                             task,
                                     gmx::ArrayRef<const int>        taskIndex,
                                     ArrayRef<const InteractionList> ilist,
                                     ArrayRef<const t_iparams>       ip)
{
    for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
    {
        tData->ilist[ftype].clear();
        tData->idTask.ilist[ftype].clear();

        int              nral1   = 1 + NRAL(ftype);
        int              inc     = nral1;
        const int*       iat     = ilist[ftype].iatoms.data();
        InteractionList* il_task = &tData->ilist[ftype];

        for (int i = 0; i < ilist[ftype].size();)
        {
            if (ftype == F_VSITEN)
            {
                /* The 3 below is from 1+NRAL(ftype)=3 */
                inc = ip[iat[i]].vsiten.n * 3;
            }
            /* Check if the vsite is assigned to our task */
            if (taskIndex[iat[1 + i]] == task)
            {
                /* Copy the vsite data to the thread-task local array */
                il_task->push_back(iat[i], inc - 1, iat + i + 1);
            }

            i += inc;
        }
    }
}

void ThreadingInfo::setVirtualSites(ArrayRef<const InteractionList> ilists,
                                    ArrayRef<const t_iparams>       iparams,
                                    const int                       numAtoms,
                                    const int                       homenr,
                                    ArrayRef<const ParticleType>    ptype,
                                    const bool                      useDomdec)
{
    if (numThreads_ <= 1)
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
    int vsite_atom_range;
    int natperthread;
    if (!useDomdec)
    {
        vsite_atom_range = -1;
        for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
        {
            { // TODO remove me
                if (ftype != F_VSITEN)
                {
                    int                 nral1 = 1 + NRAL(ftype);
                    ArrayRef<const int> iat   = ilists[ftype].iatoms;
                    for (int i = 0; i < ilists[ftype].size(); i += nral1)
                    {
                        for (int j = i + 1; j < i + nral1; j++)
                        {
                            vsite_atom_range = std::max(vsite_atom_range, iat[j]);
                        }
                    }
                }
                else
                {
                    int vs_ind_end;

                    ArrayRef<const int> iat = ilists[ftype].iatoms;

                    int i = 0;
                    while (i < ilists[ftype].size())
                    {
                        /* The 3 below is from 1+NRAL(ftype)=3 */
                        vs_ind_end = i + iparams[iat[i]].vsiten.n * 3;

                        vsite_atom_range = std::max(vsite_atom_range, iat[i + 1]);
                        while (i < vs_ind_end)
                        {
                            vsite_atom_range = std::max(vsite_atom_range, iat[i + 2]);
                            i += 3;
                        }
                    }
                }
            }
        }
        vsite_atom_range++;
        natperthread = (vsite_atom_range + numThreads_ - 1) / numThreads_;
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
        vsite_atom_range = numAtoms;
        natperthread     = (homenr + numThreads_ - 1) / numThreads_;
    }

    if (debug)
    {
        fprintf(debug,
                "virtual site thread dist: natoms %d, range %d, natperthread %d\n",
                numAtoms,
                vsite_atom_range,
                natperthread);
    }

    /* To simplify the vsite assignment, we make an index which tells us
     * to which task particles, both non-vsites and vsites, are assigned.
     */
    taskIndex_.resize(numAtoms);

    /* Initialize the task index array. Here we assign the non-vsite
     * particles to task=thread, so we easily figure out if vsites
     * depend on local and/or non-local particles in assignVsitesToThread.
     */
    {
        int thread = 0;
        for (int i = 0; i < numAtoms; i++)
        {
            if (ptype[i] == ParticleType::VSite)
            {
                /* vsites are not assigned to a task yet */
                taskIndex_[i] = -1;
            }
            else
            {
                /* assign non-vsite particles to task thread */
                taskIndex_[i] = thread;
            }
            if (i == (thread + 1) * natperthread && thread < numThreads_)
            {
                thread++;
            }
        }
    }

#pragma omp parallel num_threads(numThreads_)
    {
        try
        {
            int          thread = gmx_omp_get_thread_num();
            VsiteThread& tData  = *tData_[thread];

            /* Clear the buffer use flags that were set before */
            if (tData.useInterdependentTask)
            {
                InterdependentTask& idTask = tData.idTask;

                /* To avoid an extra OpenMP barrier in spread_vsite_f,
                 * we clear the force buffer at the next step,
                 * so we need to do it here as well.
                 */
                clearTaskForceBufferUsedElements(&idTask);

                idTask.vsite.resize(0);
                for (int t = 0; t < numThreads_; t++)
                {
                    AtomIndex& atomIndex = idTask.atomIndex[t];
                    int        natom     = atomIndex.atom.size();
                    for (int i = 0; i < natom; i++)
                    {
                        idTask.use[atomIndex.atom[i]] = false;
                    }
                    atomIndex.atom.resize(0);
                }
                idTask.nuse = 0;
            }

            /* To avoid large f_buf allocations of #threads*vsite_atom_range
             * we don't use the interdependent tasks with more than 200000 atoms.
             * This doesn't affect performance, since with such a large range
             * relatively few vsites will end up in the separate task.
             * Note that useInterdependentTask should be the same for all threads.
             */
            const int c_maxNumLocalAtomsForInterdependentTask = 200000;
            tData.useInterdependentTask = (vsite_atom_range <= c_maxNumLocalAtomsForInterdependentTask);
            if (tData.useInterdependentTask)
            {
                size_t              natoms_use_in_vsites = vsite_atom_range;
                InterdependentTask& idTask               = tData.idTask;
                /* To avoid resizing and re-clearing every nstlist steps,
                 * we never down size the force buffer.
                 */
                if (natoms_use_in_vsites > idTask.force.size() || natoms_use_in_vsites > idTask.use.size())
                {
                    idTask.force.resize(natoms_use_in_vsites, { 0, 0, 0 });
                    idTask.use.resize(natoms_use_in_vsites, false);
                }
            }

            /* Assign all vsites that can execute independently on threads */
            tData.rangeStart = thread * natperthread;
            if (thread < numThreads_ - 1)
            {
                tData.rangeEnd = (thread + 1) * natperthread;
            }
            else
            {
                /* The last thread should cover up to the end of the range */
                tData.rangeEnd = numAtoms;
            }
            assignVsitesToThread(
                    &tData, thread, numThreads_, natperthread, taskIndex_, ilists, iparams, ptype);

            if (tData.useInterdependentTask)
            {
                /* In the worst case, all tasks write to force ranges of
                 * all other tasks, leading to #tasks^2 scaling (this is only
                 * the overhead, the actual flops remain constant).
                 * But in most cases there is far less coupling. To improve
                 * scaling at high thread counts we therefore construct
                 * an index to only loop over the actually affected tasks.
                 */
                InterdependentTask& idTask = tData.idTask;

                /* Ensure assignVsitesToThread finished on other threads */
#pragma omp barrier

                idTask.spreadTask.resize(0);
                idTask.reduceTask.resize(0);
                for (int t = 0; t < numThreads_; t++)
                {
                    /* Do we write to the force buffer of task t? */
                    if (!idTask.atomIndex[t].atom.empty())
                    {
                        idTask.spreadTask.push_back(t);
                    }
                    /* Does task t write to our force buffer? */
                    if (!tData_[t]->idTask.atomIndex[thread].atom.empty())
                    {
                        idTask.reduceTask.push_back(t);
                    }
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    /* Assign all remaining vsites, that will have taskIndex[]=2*vsite->nthreads,
     * to a single task that will not run in parallel with other tasks.
     */
    assignVsitesToSingleTask(tData_[numThreads_].get(), 2 * numThreads_, taskIndex_, ilists, iparams);

    if (debug && numThreads_ > 1)
    {
        fprintf(debug,
                "virtual site useInterdependentTask %d, nuse:\n",
                static_cast<int>(tData_[0]->useInterdependentTask));
        for (int th = 0; th < numThreads_ + 1; th++)
        {
            fprintf(debug, " %4d", tData_[th]->idTask.nuse);
        }
        fprintf(debug, "\n");

        for (int ftype = c_ftypeVsiteStart; ftype < c_ftypeVsiteEnd; ftype++)
        {
            if (!ilists[ftype].empty())
            {
                fprintf(debug, "%-20s thread dist:", interaction_function[ftype].longname);
                for (int th = 0; th < numThreads_ + 1; th++)
                {
                    fprintf(debug,
                            " %4d %4d ",
                            tData_[th]->ilist[ftype].size(),
                            tData_[th]->idTask.ilist[ftype].size());
                }
                fprintf(debug, "\n");
            }
        }
    }

#ifndef NDEBUG
    int nrOrig     = vsiteIlistNrCount(ilists);
    int nrThreaded = 0;
    for (int th = 0; th < numThreads_ + 1; th++)
    {
        nrThreaded += vsiteIlistNrCount(tData_[th]->ilist) + vsiteIlistNrCount(tData_[th]->idTask.ilist);
    }
    GMX_ASSERT(nrThreaded == nrOrig,
               "The number of virtual sites assigned to all thread task has to match the total "
               "number of virtual sites");
#endif
}

void VirtualSitesHandler::Impl::setVirtualSites(ArrayRef<const InteractionList> ilists,
                                                const int                       numAtoms,
                                                const int                       homenr,
                                                ArrayRef<const ParticleType>    ptype)
{
    ilists_ = ilists;

    threadingInfo_.setVirtualSites(ilists, iparams_, numAtoms, homenr, ptype, domainInfo_.useDomdec());
}

void VirtualSitesHandler::setVirtualSites(ArrayRef<const InteractionList> ilists,
                                          const int                       numAtoms,
                                          const int                       homenr,
                                          ArrayRef<const ParticleType>    ptype)
{
    impl_->setVirtualSites(ilists, numAtoms, homenr, ptype);
}

} // namespace gmx

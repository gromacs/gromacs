/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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

/*! \libinternal \file
 *
 * \brief This file declares functions for timing the load imbalance due to domain decomposition.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_DLBTIMING_H
#define GMX_DOMDEC_DLBTIMING_H

#include <memory>

#include "gromacs/domdec/domdec_struct.h"

struct gmx_domdec_t;
struct t_nrnb;

/*! \brief Tells if we should open the balancing region */
enum class DdAllowBalanceRegionReopen
{
    no, //!< Do not allow opening an already open region
    yes //!< Allow opening an already open region
};

/*! \brief Tells if we had to wait for a GPU to finish computation */
enum class DdBalanceRegionWaitedForGpu
{
    no, //!< The GPU finished computation before the CPU needed the result
    yes //!< We had to wait for the GPU to finish computation
};

/*! \brief Re-open the, already opened, load balance timing region
 *
 * This function should be called after every MPI communication that occurs
 * in the main MD loop.
 * Note that the current setup assumes that all MPI communication acts like
 * a global barrier. But if some ranks don't participate in communication
 * or if some ranks communicate faster with neighbors than others,
 * the obtained timings might not accurately reflect the computation time.
 *
 * \param[in,out] dd  The domain decomposition struct
 */
void ddReopenBalanceRegionCpu(const gmx_domdec_t* dd);

/*! \libinternal
 * \brief Manager for starting and stopping the dynamic load balancing region
 */
class DDBalanceRegionHandler
{
public:
    //! Constructor, pass a pointer to gmx_domdec_t or nullptr when not using domain decomposition
    DDBalanceRegionHandler(gmx_domdec_t* dd) :
        useBalancingRegion_(dd != nullptr ? havePPDomainDecomposition(dd) : false), dd_(dd)
    {
    }

    /*! \brief Returns whether were are actually using the balancing region
     */
    bool useBalancingRegion() const { return useBalancingRegion_; }

    /*! \brief Open the load balance timing region on the CPU
     *
     * Opens the balancing region for timing how much time it takes to perform
     * the (balancable part of) the MD step. This should be called right after
     * the last communication during the previous step to maximize the region.
     * In practice this means right after the force communication finished
     * or just before neighbor search at search steps.
     * It is assumed that computation done in the region either scales along
     * with the domain size or takes constant time.
     *
     * \param[in] allowReopen  Allows calling with a potentially already opened region
     */
    void openBeforeForceComputationCpu(DdAllowBalanceRegionReopen allowReopen) const
    {
        if (useBalancingRegion_)
        {
            openRegionCpuImpl(allowReopen);
        }
    }

    /*! \brief Open the load balance timing region for the CPU
     *
     * This can only be called within a region that is open on the CPU side.
     */
    void openBeforeForceComputationGpu() const
    {
        if (useBalancingRegion_)
        {
            openRegionGpuImpl();
        }
    }

    /*! \brief Re-open the, already opened, load balance timing region
     *
     * This function should be called after every MPI communication that occurs
     * in the main MD loop.
     * Note that the current setup assumes that all MPI communication acts like
     * a global barrier. But if some ranks don't participate in communication
     * or if some ranks communicate faster with neighbors than others,
     * the obtained timings might not accurately reflect the computation time.
     */
    void reopenRegionCpu() const
    {
        if (useBalancingRegion_)
        {
            ddReopenBalanceRegionCpu(dd_);
        }
    }

    /*! \brief Close the load balance timing region on the CPU side
     */
    void closeAfterForceComputationCpu() const
    {
        if (useBalancingRegion_)
        {
            closeRegionCpuImpl();
        }
    }

    /*! \brief Close the load balance timing region on the GPU side
     *
     * This should be called after the CPU receives the last (local) results
     * from the GPU. The wait time for these results is estimated, depending
     * on the \p waitedForGpu parameter.
     * If called on an already closed region, this call does nothing.
     *
     * \param[in] waitCyclesGpuInCpuRegion  The time we waited for the GPU earlier, overlapping completely with the open CPU region
     * \param[in] waitedForGpu              Tells if we waited for the GPU to finish now
     */
    void closeAfterForceComputationGpu(float                       waitCyclesGpuInCpuRegion,
                                       DdBalanceRegionWaitedForGpu waitedForGpu) const
    {
        if (useBalancingRegion_)
        {
            closeRegionGpuImpl(waitCyclesGpuInCpuRegion, waitedForGpu);
        }
    }

private:
    /*! \brief Open the load balance timing region on the CPU
     *
     * \param[in] allowReopen  Allows calling with a potentially already opened region
     */
    void openRegionCpuImpl(DdAllowBalanceRegionReopen allowReopen) const;

    /*! \brief Open the load balance timing region for the GPU
     *
     * This can only be called within a region that is open on the CPU side.
     */
    void openRegionGpuImpl() const;

    /*! \brief Close the load balance timing region on the CPU side
     */
    void closeRegionCpuImpl() const;

    /*! \brief Close the load balance timing region on the GPU side
     *
     * \param[in] waitCyclesGpuInCpuRegion  The time we waited for the GPU earlier, overlapping completely with the open CPU region
     * \param[in] waitedForGpu              Tells if we waited for the GPU to finish now
     */
    void closeRegionGpuImpl(float waitCyclesGpuInCpuRegion, DdBalanceRegionWaitedForGpu waitedForGpu) const;

    //! Tells whether the balancing region should be active
    bool useBalancingRegion_;
    //! A pointer to the DD struct, only valid with useBalancingRegion_=true
    gmx_domdec_t* dd_;
};

//! Object that describes a DLB balancing region
class BalanceRegion
{
public:
    BalanceRegion();
    ~BalanceRegion();
    // Not private because used by DDBalanceRegionHandler
    // private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

/*! \brief Start the force flop count */
void dd_force_flop_start(struct gmx_domdec_t* dd, t_nrnb* nrnb);

/*! \brief Stop the force flop count */
void dd_force_flop_stop(struct gmx_domdec_t* dd, t_nrnb* nrnb);

//! Clear the cycle counts used for tuning.
void clear_dd_cycle_counts(gmx_domdec_t* dd);

#endif

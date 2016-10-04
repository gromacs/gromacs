/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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

struct BalanceRegion;
struct gmx_domdec_t;

/*! \brief Tells if we should open the balancing region */
enum class DdOpenBalanceRegionBeforeForceComputation
{
    no,  //!< Do not open a balancing region
    yes  //!< Open the balancing region before update or after pair-search
};

/*! \brief Tells if we should close the balancing region after the force computation has completed */
enum class DdCloseBalanceRegionAfterForceComputation
{
    no,  //!< Do not close a balancing region
    yes  //!< Close the balancing region after computation completed
};

/*! \brief Tells if we should open the balancing region */
enum class DdAllowBalanceRegionReopen
{
    no,  //!< Do not allow opening an already open region
    yes  //!< Allow opening an already open region
};

/*! \brief Tells if we had to wait for a GPU to finish computation */
enum class DdBalanceRegionWaitedForGpu
{
    no,  //!< The GPU finished computation before the CPU needed the result
    yes  //!< We had to wait for the GPU to finish computation
};

/*! \brief Returns a pointer to a constructed \p BalanceRegion struct
 *
 * Should be replaced by a proper constructor once BalanceRegion is a proper
 * class (requires restructering in domdec.cpp).
 */
BalanceRegion *ddBalanceRegionAllocate();

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
 * \param[in,out] dd           The domain decomposition struct
 * \param[in]     allowReopen  Allows calling with a potentially already opened region
 */
void ddOpenBalanceRegionCpu(const gmx_domdec_t         *dd,
                            DdAllowBalanceRegionReopen  allowReopen);

/*! \brief Open the load balance timing region for the CPU
 *
 * This can only be called within a region that is open on the CPU side.
 */
void ddOpenBalanceRegionGpu(const gmx_domdec_t *dd);

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
void ddReopenBalanceRegionCpu(const gmx_domdec_t *dd);

/*! \brief Close the load balance timing region on the CPU side
 *
 * \param[in,out] dd  The domain decomposition struct
 */
void ddCloseBalanceRegionCpu(const gmx_domdec_t *dd);

/*! \brief Close the load balance timing region on the GPU side
 *
 * This should be called after the CPU receives the last (local) results
 * from the GPU. The wait time for these results is estimated, depending
 * on the \p waitedForGpu parameter.
 * If called on an already closed region, this call does nothing.
 *
 * \param[in,out] dd                        The domain decomposition struct
 * \param[in]     waitCyclesGpuInCpuRegion  The time we waited for the GPU earlier, overlapping completely with the open CPU region
 * \param[in]     waitedForGpu              Tells if we waited for the GPU to finish now
 */
void ddCloseBalanceRegionGpu(const gmx_domdec_t          *dd,
                             float                        waitCyclesGpuInCpuRegion,
                             DdBalanceRegionWaitedForGpu  waitedForGpu);

#endif

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 *
 * \brief This file contains declarations for BondedThreading
 * class.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_BONDED_THREADING_H
#define GMX_LISTED_FORCES_BONDED_THREADING_H

#include <vector>

#include "gromacs/mdtypes/threaded_force_buffer.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/vectypes.h"

#include "listed_internal.h"

namespace gmx
{
class MDLogger;

/*! \internal \brief The division of bonded interactions of the threads */
class WorkDivision
{
public:
    //! Constructor
    WorkDivision(int numThreads) :
        stride_(numThreads + 1), packedBounds_(static_cast<int>(InteractionFunction::Count) * stride_)
    {
    }

    //! Sets the bound between threads \p boundIndex-1 and \p boundIndex to \p count
    void setBound(InteractionFunction functionType, int boundIndex, int count)
    {
        packedBounds_[static_cast<int>(functionType) * stride_ + boundIndex] = count;
    }

    int bound(InteractionFunction functionType, int boundIndex) const
    {
        return bound(static_cast<int>(functionType), boundIndex);
    }

    //! Returns the bound between threads \p boundIndex-1 and \p boundIndex
    int bound(int functionType, int boundIndex) const
    {
        return packedBounds_[functionType * stride_ + boundIndex];
    }

    //! Returns the last bound
    int end(int ftype) const { return bound(ftype, stride_ - 1); }

private:
    //! The stride_ between and size of the entries for a function type
    int stride_;
    //! The bounds stored as a flat array for fast access
    std::vector<int> packedBounds_;
};

struct CentersOfMassScaledBuffers
{
    //! Buffer for computing scaled centers of mass for position restraints
    std::vector<RVec> comA_;
    //! Buffer for computing scaled centers of mass for topology B for position restraints
    std::vector<RVec> comB_;
};

/*! \brief struct for passing all data required for a function type */
struct IListData
{
    const InteractionList* il;       /**< pointer to t_ilist entry corresponding to ftype */
    InteractionFunction    ftype;    /**< the function type index */
    int                    numAtoms; /**< Number of atoms involved in a single ftype interaction */
};


/*! \internal \brief struct contain all data for bonded force threading */
class BondedThreading
{
public:
    //! Constructor
    BondedThreading(int numThreads, int numEnergyGroups, int numComGroups, const MDLogger& mdlog);

    //! Number of threads to be used for bondeds
    int numThreads_ = 0;

    //! The thread parallel force and energy buffers
    ThreadedForceBuffer<rvec4> threadedForceBuffers_;
    //! true if we have and thus need to reduce bonded forces
    bool haveBondeds_ = false;

private:
    /* There are two different ways to distribute the bonded force calculation
     * over the threads. We decide which to use based on the number of threads.
     */
    //! Maximum thread count for uniform distribution of bondeds over threads
    int maxNthreadUniform_ = 0;

public:
    //! The division of work in the t_list over threads.
    WorkDivision workDivision_;

    //! Work division for free-energy foreign lambda calculations, always uses 1 thread
    WorkDivision foreignLambdaWorkDivision_;

    //! Buffers for each thread for computing scaled centers of mass for position restraints
    std::vector<CentersOfMassScaledBuffers> centersOfMassScaledBuffers_;

    //! Update data structures after possible DD repartition
    void updateAfterPartition(int numAtomsForce, bool useGpuForBondeds, const InteractionDefinitions& idef);
    GMX_DISALLOW_COPY_MOVE_AND_ASSIGN(BondedThreading);

private:
    //! Divides bonded interactions over threads and the GPU
    void divideWork(bool useGpuForBondeds, const InteractionDefinitions& idef);
    /*! \brief Divides listed interactions over threads
     *
     * This routine attempts to divide all interactions of the numType bondeds
     * types stored in ild over the threads such that each thread has roughly
     * equal load and different threads avoid touching the same atoms as much
     * as possible so that reduction cost is minimized.
     *
     * Note that with few threads, a simpler division that focuses on
     * balancing load is used.
     */
    void divideByLocality(int numType, const EnumerationArray<InteractionFunction, IListData>& ild);
    //! Construct a reduction mask for the blocks of the force array touched by this thread
    void calculateReductionMask(int natoms, const InteractionDefinitions& idef, int thread);
};

} // namespace gmx

#endif

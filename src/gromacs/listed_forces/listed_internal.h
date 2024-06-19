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
 * \brief This file contains declarations for functions needed
 * internally by the module.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_LISTED_INTERNAL_H
#define GMX_LISTED_FORCES_LISTED_INTERNAL_H

#include <cstdio>

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/threaded_force_buffer.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/bitmask.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/enumerationhelpers.h"

/*! \internal \brief The division of bonded interactions of the threads */
class WorkDivision
{
public:
    //! Constructor
    WorkDivision(int numThreads) : stride_(numThreads + 1), packedBounds_(F_NRE * stride_) {}

    //! Sets the bound between threads \p boundIndex-1 and \p boundIndex to \p count
    void setBound(int functionType, int boundIndex, int count)
    {
        packedBounds_[functionType * stride_ + boundIndex] = count;
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

/*! \internal \brief struct contain all data for bonded force threading */
struct bonded_threading_t
{
    //! Constructor
    bonded_threading_t(int numThreads, int numEnergyGroups, FILE* fplog);

    //! Number of threads to be used for bondeds
    int nthreads = 0;
    //! The thread parallel force and energy buffers
    gmx::ThreadedForceBuffer<rvec4> threadedForceBuffer;
    //! true if we have and thus need to reduce bonded forces
    bool haveBondeds = false;

    /* There are two different ways to distribute the bonded force calculation
     * over the threads. We decide which to use based on the number of threads.
     */
    //! Maximum thread count for uniform distribution of bondeds over threads
    int max_nthread_uniform = 0;

    //! The division of work in the t_list over threads.
    WorkDivision workDivision;

    //! Work division for free-energy foreign lambda calculations, always uses 1 thread
    WorkDivision foreignLambdaWorkDivision;

    GMX_DISALLOW_COPY_MOVE_AND_ASSIGN(bonded_threading_t);
};


/*! \brief Returns the global topology atom number belonging to local
 * atom index i.
 *
 * This function is intended for writing ascii output and returns atom
 * numbers starting at 1.  When global_atom_index=NULL returns i+1.
 */
int glatnr(const int* global_atom_index, int i);

#endif

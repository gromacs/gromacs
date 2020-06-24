/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017,2018 by the GROMACS development team.
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
 * \brief This file contains declarations for functions needed
 * internally by the module.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_LISTED_INTERNAL_H
#define GMX_LISTED_FORCES_LISTED_INTERNAL_H

#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/bitmask.h"
#include "gromacs/utility/classhelpers.h"

/* We reduce the force array in blocks of 32 atoms. This is large enough
 * to not cause overhead and 32*sizeof(rvec) is a multiple of the cache-line
 * size on all systems.
 */
static const int reduction_block_size = 32; /**< Force buffer block size in atoms*/
static const int reduction_block_bits = 5;  /**< log2(reduction_block_size) */

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

/*! \internal \brief struct with output for bonded forces, used per thread */
struct f_thread_t
{
    //! Constructor
    f_thread_t(int numEnergyGroups);

    ~f_thread_t() = default;

    //! Force array pointer, equals fBuffer.data(), needed because rvec4 is not a C++ type
    rvec4* f = nullptr;
    //! Force array buffer
    std::vector<real, gmx::AlignedAllocator<real>> fBuffer;
    //! Mask for marking which parts of f are filled, working array for constructing mask in bonded_threading_t
    std::vector<gmx_bitmask_t> mask;
    //! Number of blocks touched by our thread
    int nblock_used = 0;
    //! Index to touched blocks
    std::vector<int> block_index;

    //! Shift force array, size SHIFTS
    std::vector<gmx::RVec> fshift;
    //! Energy array
    real ener[F_NRE];
    //! Group pair energy data for pairs
    gmx_grppairener_t grpp;
    //! Free-energy dV/dl output
    real dvdl[efptNR];

    GMX_DISALLOW_COPY_MOVE_AND_ASSIGN(f_thread_t);
};

/*! \internal \brief struct contain all data for bonded force threading */
struct bonded_threading_t
{
    //! Constructor
    bonded_threading_t(int numThreads, int numEnergyGroups, FILE* fplog);

    //! Number of threads to be used for bondeds
    int nthreads = 0;
    //! Force/energy data per thread, size nthreads, stored in unique_ptr to allow thread local allocation
    std::vector<std::unique_ptr<f_thread_t>> f_t;
    //! The number of force blocks to reduce
    int nblock_used = 0;
    //! Index of size nblock_used into mask
    std::vector<int> block_index;
    //! Mask array, one element corresponds to a block of reduction_block_size atoms of the force array, bit corresponding to thread indices set if a thread writes to that block
    std::vector<gmx_bitmask_t> mask;
    //! true if we have and thus need to reduce bonded forces
    bool haveBondeds = false;
    //! The number of atoms forces are computed for
    int numAtomsForce = 0;

    /* There are two different ways to distribute the bonded force calculation
     * over the threads. We dedice which to use based on the number of threads.
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

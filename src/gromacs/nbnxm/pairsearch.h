/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016 by the GROMACS development team.
 * Copyright (c) 2017,2018,2019,2020,2021, by the GROMACS development team, led by
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
 * \brief Declares the PairSearch class and helper structs
 *
 * The PairSearch class holds the domain setup, the search grids
 * and helper object for the pair search. It manages the search work.
 * The actual gridding and pairlist generation is performed by the
 * GridSet/Grid and PairlistSet/Pairlist classes, respectively.
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_PAIRSEARCH_H
#define GMX_NBNXM_PAIRSEARCH_H

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "gridset.h"
#include "pairlist.h"

struct gmx_domdec_zones_t;
struct PairsearchWork;


/*! \brief Convenience declaration for an std::vector with aligned memory */
template<class T>
using AlignedVector = std::vector<T, gmx::AlignedAllocator<T>>;


//! Local cycle count struct for profiling \internal
class nbnxn_cycle_t
{
public:
    //! Start counting cycles
    void start() { start_ = gmx_cycles_read(); }
    //! Stop counting cycles
    void stop()
    {
        cycles_ += gmx_cycles_read() - start_;
        count_++;
    }
    //! Return the number of periods of cycle counting
    int count() const { return count_; }

    //! Return the average number of million cycles per counting period
    double averageMCycles() const
    {
        if (count_ > 0)
        {
            return static_cast<double>(cycles_) * 1e-6 / count_;
        }
        else
        {
            return 0;
        }
    }

private:
    //! Number of counting periods
    int count_ = 0;
    //! Total cycles in all counting periods
    gmx_cycles_t cycles_ = 0;
    //! Cycle count at the most recent start
    gmx_cycles_t start_ = 0;
};

//! Local cycle count enum for profiling different parts of search
enum
{
    enbsCCgrid,
    enbsCCsearch,
    enbsCCcombine,
    enbsCCnr
};

/*! \internal
 * \brief Struct for collecting detailed cycle counts for the search
 */
struct SearchCycleCounting
{
    //! Start a pair search cycle counter
    void start(const int enbsCC) { cc_[enbsCC].start(); }

    //! Stop a pair search cycle counter
    void stop(const int enbsCC) { cc_[enbsCC].stop(); }

    //! Print the cycle counts to \p fp
    void printCycles(FILE* fp, gmx::ArrayRef<const PairsearchWork> work) const;

    //! Tells whether we record cycles
    bool recordCycles_ = false;
    //! The number of times pairsearching has been performed, local+non-local count as 1
    int searchCount_ = 0;
    //! The set of cycle counters
    nbnxn_cycle_t cc_[enbsCCnr];
};

// TODO: Move nbnxn_search_work_t definition to its own file

//! Thread-local work struct, contains working data for Grid \internal
struct PairsearchWork
{
    PairsearchWork();

    ~PairsearchWork();

    //! Buffer to avoid cache pollution
    gmx_cache_protect_t cp0;

    //! Temporary buffer for sorting atoms within a grid column
    std::vector<int> sortBuffer;

    //! Flags for force buffer access
    std::vector<gmx_bitmask_t> buffer_flags;

    //! Number of distance checks for flop counting
    int ndistc;


    //! Temporary FEP list for load balancing
    std::unique_ptr<t_nblist> nbl_fep;

    //! Counter for thread-local cycles
    nbnxn_cycle_t cycleCounter;

    //! Buffer to avoid cache polution
    gmx_cache_protect_t cp1;
};

//! Main pair-search struct, contains the grid(s), not the pair-list(s) \internal
class PairSearch
{
public:
    //! Puts the atoms in \p ddZone on the grid and copies the coordinates to \p nbat
    void putOnGrid(const matrix                   box,
                   int                            ddZone,
                   const rvec                     lowerCorner,
                   const rvec                     upperCorner,
                   const gmx::UpdateGroupsCog*    updateGroupsCog,
                   gmx::Range<int>                atomRange,
                   real                           atomDensity,
                   gmx::ArrayRef<const int64_t>   atomInfo,
                   gmx::ArrayRef<const gmx::RVec> x,
                   int                            numAtomsMoved,
                   const int*                     move,
                   nbnxn_atomdata_t*              nbat)
    {
        cycleCounting_.start(enbsCCgrid);

        gridSet_.putOnGrid(box,
                           ddZone,
                           lowerCorner,
                           upperCorner,
                           updateGroupsCog,
                           atomRange,
                           atomDensity,
                           atomInfo,
                           x,
                           numAtomsMoved,
                           move,
                           nbat);

        cycleCounting_.stop(enbsCCgrid);
    }

    /*! \brief Constructor
     *
     * \param[in] pbcType                  The periodic boundary conditions
     * \param[in] doTestParticleInsertion  Whether test-particle insertion is active
     * \param[in] numDDCells               The number of domain decomposition cells per dimension, without DD nullptr should be passed
     * \param[in] zones                    The domain decomposition zone setup, without DD nullptr should be passed
     * \param[in] pairlistType             The type of tte pair list
     * \param[in] haveFep                  Tells whether non-bonded interactions are perturbed
     * \param[in] maxNumThreads            The maximum number of threads used in the search
     * \param[in] pinningPolicy            Sets the pinning policy for all buffers used on the GPU
     */
    PairSearch(PbcType                   pbcType,
               bool                      doTestParticleInsertion,
               const ivec*               numDDCells,
               const gmx_domdec_zones_t* zones,
               PairlistType              pairlistType,
               bool                      haveFep,
               int                       maxNumThreads,
               gmx::PinningPolicy        pinningPolicy);

    //! Sets the order of the local atoms to the order grid atom ordering
    void setLocalAtomOrder() { gridSet_.setLocalAtomOrder(); }

    //! Returns the set of search grids
    const Nbnxm::GridSet& gridSet() const { return gridSet_; }

    //! Returns the list of thread-local work objects
    gmx::ArrayRef<const PairsearchWork> work() const { return work_; }

    //! Returns the list of thread-local work objects
    gmx::ArrayRef<PairsearchWork> work() { return work_; }

private:
    //! The set of search grids
    Nbnxm::GridSet gridSet_;
    //! Work objects, one entry for each thread
    std::vector<PairsearchWork> work_;

public:
    //! Cycle counting for measuring components of the search
    SearchCycleCounting cycleCounting_;
};

#endif

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
 * \brief Declares internal nbnxm module details
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_INTERNAL_H
#define GMX_NBNXM_INTERNAL_H

#include <memory>
#include <vector>

#include "gromacs/domdec/domdec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/simd/simd.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "gridset.h"

struct gmx_domdec_zones_t;


/*! \brief Convenience declaration for an std::vector with aligned memory */
template <class T>
using AlignedVector = std::vector < T, gmx::AlignedAllocator < T>>;


/* Local cycle count struct for profiling */
typedef struct {
    int          count;
    gmx_cycles_t c;
    gmx_cycles_t start;
} nbnxn_cycle_t;

/* Local cycle count enum for profiling */
enum {
    enbsCCgrid, enbsCCsearch, enbsCCcombine, enbsCCreducef, enbsCCnr
};

/* Thread-local work struct, contains working data for Grid */
struct nbnxn_search_work_t
{
    nbnxn_search_work_t();

    ~nbnxn_search_work_t();

    gmx_cache_protect_t       cp0;          /* Buffer to avoid cache polution */

    std::vector<int>          sortBuffer;   /* Temporary buffer for sorting atoms within a grid column */

    nbnxn_buffer_flags_t      buffer_flags; /* Flags for force buffer access */

    int                       ndistc;       /* Number of distance checks for flop counting */


    std::unique_ptr<t_nblist> nbl_fep;      /* Temporary FEP list for load balancing */

    nbnxn_cycle_t             cc[enbsCCnr]; /* Counters for thread-local cycles */

    gmx_cache_protect_t       cp1;          /* Buffer to avoid cache polution */
};

/* Main pair-search struct, contains the grid(s), not the pair-list(s) */
struct nbnxn_search
{
    /* \brief Constructor
     *
     * \param[in] ePBC            The periodic boundary conditions
     * \param[in] n_dd_cells      The number of domain decomposition cells per dimension, without DD nullptr should be passed
     * \param[in] zones           The domain decomposition zone setup, without DD nullptr should be passed
     * \param[in] haveFep         Tells whether non-bonded interactions are perturbed
     * \param[in] maxNumThreads   The maximum number of threads used in the search
     */
    nbnxn_search(int                       ePBC,
                 const ivec               *n_dd_cells,
                 const gmx_domdec_zones_t *zones,
                 PairlistType              pairlistType,
                 bool                      haveFep,
                 int                       maxNumthreads);

    const Nbnxm::GridSet &gridSet() const
    {
        return gridSet_;
    }

    int                        ePBC;            /* PBC type enum                              */
    gmx_bool                   DomDec;          /* Are we doing domain decomposition?         */
    ivec                       dd_dim;          /* Are we doing DD in x,y,z?                  */
    const gmx_domdec_zones_t  *zones;           /* The domain decomposition zones        */

    //! The local and non-local grids
    Nbnxm::GridSet       gridSet_;

    gmx_bool             print_cycles;
    int                  search_count;
    nbnxn_cycle_t        cc[enbsCCnr];

    /* Thread-local work data */
    mutable std::vector<nbnxn_search_work_t> work; /* Work array, one entry for each thread */
};


/*! \brief Start an nbnxn cycle counter */
static inline void nbs_cycle_start(nbnxn_cycle_t *cc)
{
    cc->start = gmx_cycles_read();
}

/*! \brief Stop an nbnxn cycle counter */
static inline void nbs_cycle_stop(nbnxn_cycle_t *cc)
{
    cc->c += gmx_cycles_read() - cc->start;
    cc->count++;
}

#endif

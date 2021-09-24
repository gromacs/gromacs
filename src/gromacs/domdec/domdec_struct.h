/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2018,2019,2020,2021, by the GROMACS development team, led by
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
 * \brief Declares structures related to domain decomposition.
 *
 * \author Berk Hess <hess@kth.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inlibraryapi
 * \ingroup module_domdec
 */
#ifndef GMX_DOMDEC_DOMDEC_STRUCT_H
#define GMX_DOMDEC_DOMDEC_STRUCT_H

#include <cstddef>

#include <array>
#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/real.h"

//! Max number of zones in domain decomposition
#define DD_MAXZONE 8
//! Max number of izones in domain decomposition
#define DD_MAXIZONE 4

struct AtomDistribution;
struct gmx_domdec_comm_t;
struct gmx_domdec_constraints_t;
struct gmx_domdec_specat_comm_t;
class gmx_ga2la_t;
struct gmx_pme_comm_n_box_t;
struct t_inputrec;
class gmx_reverse_top_t;
struct gmx_mtop_t;
struct ReverseTopOptions;

namespace gmx
{
template<typename T>
class HashedMap;
class LocalAtomSetManager;
class LocalTopologyChecker;
class GpuHaloExchange;
} // namespace gmx

/*! \internal
 * \brief Pair interaction zone and atom range for an i-zone
 */
struct DDPairInteractionRanges
{
    //! The index of this i-zone in the i-zone list
    int iZoneIndex = -1;
    //! The range of j-zones
    gmx::Range<int> jZoneRange;
    //! The i-atom range
    gmx::Range<int> iAtomRange;
    //! The j-atom range
    gmx::Range<int> jAtomRange;
    //! Minimum shifts to consider
    gmx::IVec shift0 = { 0, 0, 0 };
    //! Maximum shifts to consider
    gmx::IVec shift1 = { 0, 0, 0 };
};

typedef struct gmx_domdec_zone_size
{
    /* Zone lower corner in triclinic coordinates         */
    gmx::RVec x0 = { 0, 0, 0 };
    /* Zone upper corner in triclinic coordinates         */
    gmx::RVec x1 = { 0, 0, 0 };
    /* Zone bounding box lower corner in Cartesian coords */
    gmx::RVec bb_x0 = { 0, 0, 0 };
    /* Zone bounding box upper corner in Cartesian coords */
    gmx::RVec bb_x1 = { 0, 0, 0 };
} gmx_domdec_zone_size_t;

struct gmx_domdec_zones_t
{
    /* The number of zones including the home zone */
    int n = 0;
    /* The shift of the zones with respect to the home zone */
    std::array<ivec, DD_MAXZONE> shift;
    /* The charge group boundaries for the zones */
    std::array<int, DD_MAXZONE + 1> cg_range;
    /* The pair interaction zone and atom ranges per each i-zone */
    std::vector<DDPairInteractionRanges> iZones;
    /* Boundaries of the zones */
    std::array<gmx_domdec_zone_size_t, DD_MAXZONE> size;
    /* The cg density of the home zone */
    real dens_zone0 = 0;
};

struct gmx_ddbox_t
{
    int       npbcdim;
    int       nboundeddim;
    gmx::RVec box0     = { 0, 0, 0 };
    gmx::RVec box_size = { 0, 0, 0 };
    /* Tells if the box is skewed for each of the three cartesian directions */
    gmx::IVec tric_dir = { 0, 0, 0 };
    gmx::RVec skew_fac = { 0, 0, 0 };
    /* Orthogonal vectors for triclinic cells, Cartesian index */
    rvec v[DIM][DIM];
    /* Normal vectors for the cells walls */
    rvec normal[DIM];
};

/*! \internal \brief Provides information about properties of the unit cell */
struct UnitCellInfo
{
    //! Constructor
    UnitCellInfo(const t_inputrec& ir);

    //! We have PBC from dim 0 (X) up to npbcdim
    int npbcdim;
    //! The system is bounded from 0 (X) to numBoundedDimensions
    int numBoundedDimensions;
    //! Tells whether the box bounding the atoms is dynamic
    bool ddBoxIsDynamic;
    //! Screw PBC?
    bool haveScrewPBC;
};

struct gmx_domdec_t
{ //NOLINT(clang-analyzer-optin.performance.Padding)
    //! Constructor, only partial for now
    gmx_domdec_t(const t_inputrec& ir);
    ~gmx_domdec_t();

    /* The DD particle-particle nodes only */
    /* The communication setup within the communicator all
     * defined in dd->comm in domdec.c
     */
    int      nnodes       = 1;
    MPI_Comm mpi_comm_all = MPI_COMM_NULL;
    /* The local DD cell index and rank */
    gmx::IVec ci         = { 0, 0, 0 };
    int       rank       = 0;
    gmx::IVec master_ci  = { 0, 0, 0 };
    int       masterrank = 0;
    /* Communication with the PME only nodes */
    int                   pme_nodeid           = 0;
    gmx_bool              pme_receive_vir_ener = false;
    gmx_pme_comm_n_box_t* cnb                  = nullptr;
    int                   nreq_pme             = 0;
    MPI_Request           req_pme[8];

    /* Properties of the unit cell */
    UnitCellInfo unitCellInfo;

    /* The communication setup, identical for each cell, cartesian index */
    //! Todo: refactor nbnxm to not rely on this sometimes being a nullptr so this can be IVec
    ivec      numCells = { 0, 0, 0 };
    int       ndim     = 0;
    gmx::IVec dim      = { 0, 0, 0 }; /* indexed by 0 to ndim */

    /* Forward and backward neighboring cells, indexed by 0 to ndim */
    int neighbor[DIM][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } };

    /* Only available on the master node */
    std::unique_ptr<AtomDistribution> ma;

    /* Global atom number to interaction list */
    std::unique_ptr<gmx_reverse_top_t> reverse_top;

    /* Whether we have non-self exclusion */
    bool haveExclusions = false;

    /* Vsite stuff */
    gmx::HashedMap<int>*      ga2la_vsite = nullptr;
    gmx_domdec_specat_comm_t* vsite_comm  = nullptr;
    std::vector<int>          vsite_requestedGlobalAtomIndices;

    /* Constraint stuff */
    gmx_domdec_constraints_t* constraints     = nullptr;
    gmx_domdec_specat_comm_t* constraint_comm = nullptr;

    /* The number of home atoms */
    int numHomeAtoms = 0;
    /* Global atom group indices for the home and all non-home groups */
    std::vector<int> globalAtomGroupIndices;

    /* Index from the local atoms to the global atoms, covers home and received zones */
    std::vector<int> globalAtomIndices;

    /* Global atom number to local atom number list */
    gmx_ga2la_t* ga2la = nullptr;

    /* Communication stuff */
    gmx_domdec_comm_t* comm = nullptr;

    /* The partioning count, to keep track of the state */
    int64_t ddp_count = 0;

    /* The managed atom sets that are updated in domain decomposition */
    gmx::LocalAtomSetManager* atomSets = nullptr;

    //! The handler for checking whether the local topology is missing interactions
    std::unique_ptr<gmx::LocalTopologyChecker> localTopologyChecker;

    /* gmx_pme_recv_f buffer */
    std::vector<gmx::RVec> pmeForceReceiveBuffer;

    /* GPU halo exchange objects: this structure supports a vector of pulses for each dimension */
    std::vector<std::unique_ptr<gmx::GpuHaloExchange>> gpuHaloExchange[DIM];
};

//! Are we the master node for domain decomposition
static inline bool DDMASTER(const gmx_domdec_t& dd)
{
    return dd.rank == dd.masterrank;
};

//! Are we the master node for domain decomposition, deprecated
static inline bool DDMASTER(const gmx_domdec_t* dd)
{
    return dd->rank == dd->masterrank;
};

#endif

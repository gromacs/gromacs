/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2018,2019, by the GROMACS development team, led by
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

#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/real.h"

//! Max number of zones in domain decomposition
#define DD_MAXZONE  8
//! Max number of izones in domain decomposition
#define DD_MAXIZONE 4

struct AtomDistribution;
struct gmx_domdec_comm_t;
struct gmx_domdec_constraints_t;
struct gmx_domdec_specat_comm_t;
class gmx_ga2la_t;
struct gmx_pme_comm_n_box_t;
struct gmx_reverse_top_t;

namespace gmx
{
template <typename T> class HashedMap;
class LocalAtomSetManager;
}

typedef struct {
    int  j0;     /* j-zone start               */
    int  j1;     /* j-zone end                 */
    int  cg1;    /* i-charge-group end         */
    int  jcg0;   /* j-charge-group start       */
    int  jcg1;   /* j-charge-group end         */
    ivec shift0; /* Minimum shifts to consider */
    ivec shift1; /* Maximum shifts to consider */
} gmx_domdec_ns_ranges_t;

typedef struct {
    rvec x0;     /* Zone lower corner in triclinic coordinates         */
    rvec x1;     /* Zone upper corner in triclinic coordinates         */
    rvec bb_x0;  /* Zone bounding box lower corner in Cartesian coords */
    rvec bb_x1;  /* Zone bounding box upper corner in Cartesian coords */
} gmx_domdec_zone_size_t;

struct gmx_domdec_zones_t {
    /* The number of zones including the home zone */
    int                    n;
    /* The shift of the zones with respect to the home zone */
    ivec                   shift[DD_MAXZONE];
    /* The charge group boundaries for the zones */
    int                    cg_range[DD_MAXZONE+1];
    /* The number of neighbor search zones with i-particles */
    int                    nizone;
    /* The neighbor search charge group ranges for each i-zone */
    gmx_domdec_ns_ranges_t izone[DD_MAXIZONE];
    /* Boundaries of the zones */
    gmx_domdec_zone_size_t size[DD_MAXZONE];
    /* The cg density of the home zone */
    real                   dens_zone0;
};

struct gmx_ddbox_t {
    int  npbcdim;
    int  nboundeddim;
    rvec box0;
    rvec box_size;
    /* Tells if the box is skewed for each of the three cartesian directions */
    ivec tric_dir;
    rvec skew_fac;
    /* Orthogonal vectors for triclinic cells, Cartesian index */
    rvec v[DIM][DIM];
    /* Normal vectors for the cells walls */
    rvec normal[DIM];
};


struct gmx_domdec_t { //NOLINT(clang-analyzer-optin.performance.Padding)
    /* The DD particle-particle nodes only */
    /* The communication setup within the communicator all
     * defined in dd->comm in domdec.c
     */
    int                    nnodes;
    MPI_Comm               mpi_comm_all;
    /* Use MPI_Sendrecv communication instead of non-blocking calls */
    gmx_bool               bSendRecv2;
    /* The local DD cell index and rank */
    ivec                   ci;
    int                    rank;
    ivec                   master_ci;
    int                    masterrank;
    /* Communication with the PME only nodes */
    int                    pme_nodeid;
    gmx_bool               pme_receive_vir_ener;
    gmx_pme_comm_n_box_t  *cnb      = nullptr;
    int                    nreq_pme = 0;
    MPI_Request            req_pme[8];


    /* The communication setup, identical for each cell, cartesian index */
    ivec     nc;
    int      ndim;
    ivec     dim; /* indexed by 0 to ndim */

    /* TODO: Move the next 4, and more from domdec_internal.h, to a simulation system */

    /* PBC from dim 0 (X) to npbcdim */
    int  npbcdim;
    /* The system is bounded from 0 (X) to numBoundedDimensions */
    int  numBoundedDimensions;
    /* Does the box size change during the simulaton? */
    bool haveDynamicBox;

    /* Screw PBC? */
    gmx_bool bScrewPBC;

    /* Forward and backward neighboring cells, indexed by 0 to ndim */
    int  neighbor[DIM][2];

    /* Only available on the master node */
    std::unique_ptr<AtomDistribution> ma;

    /* Can atoms connected by constraints be assigned to different domains? */
    bool splitConstraints;
    /* Can atoms connected by settles be assigned to different domains? */
    bool splitSettles;

    /* Global atom number to interaction list */
    gmx_reverse_top_t  *reverse_top;
    int                 nbonded_global;
    int                 nbonded_local;

    /* The number of inter charge-group exclusions */
    int  n_intercg_excl;

    /* Vsite stuff */
    gmx::HashedMap<int>       *ga2la_vsite = nullptr;
    gmx_domdec_specat_comm_t  *vsite_comm  = nullptr;
    std::vector<int>           vsite_requestedGlobalAtomIndices;

    /* Constraint stuff */
    gmx_domdec_constraints_t *constraints     = nullptr;
    gmx_domdec_specat_comm_t *constraint_comm = nullptr;

    /* The number of home atom groups */
    int                           ncg_home = 0;
    /* Global atom group indices for the home and all non-home groups */
    std::vector<int>              globalAtomGroupIndices;

    /* Index from the local atoms to the global atoms, covers home and received zones */
    std::vector<int> globalAtomIndices;

    /* Global atom number to local atom number list */
    gmx_ga2la_t  *ga2la = nullptr;

    /* Communication stuff */
    gmx_domdec_comm_t *comm;

    /* The partioning count, to keep track of the state */
    int64_t ddp_count;

    /* The managed atom sets that are updated in domain decomposition */
    gmx::LocalAtomSetManager * atomSets;

    /* gmx_pme_recv_f buffer */
    int   pme_recv_f_alloc = 0;
    rvec *pme_recv_f_buf   = nullptr;
};

//! Are we the master node for domain decomposition
static inline bool DDMASTER(const gmx_domdec_t &dd)
{
    return dd.rank == dd.masterrank;
};

//! Are we the master node for domain decomposition, deprecated
static inline bool DDMASTER(const gmx_domdec_t *dd)
{
    return dd->rank == dd->masterrank;
};

#endif

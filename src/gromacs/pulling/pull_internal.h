/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2018 by the GROMACS development team.
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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
 *
 * \brief
 * This file contains datatypes and function declarations for internal
   use in the pull code.
 *
 * \author Berk Hess
 */

#ifndef GMX_PULLING_PULL_INTERNAL_H
#define GMX_PULLING_PULL_INTERNAL_H

#include "config.h"

#include <memory>
#include <vector>

#include "gromacs/domdec/localatomset.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/utility/gmxmpi.h"

#include "pullcoordexpressionparser.h"

/*! \brief Determines up to what local atom count a pull group gets processed single-threaded.
 *
 * We set this limit to 1 with debug to catch bugs.
 * On Haswell with GCC 5 the cross-over point is around 400 atoms,
 * independent of thread count and hyper-threading.
 */
#ifdef NDEBUG
static const int c_pullMaxNumLocalAtomsSingleThreaded = 100;
#else
static const int c_pullMaxNumLocalAtomsSingleThreaded = 1;
#endif

class PullHistory;
enum class PbcType : int;

class t_state;

enum
{
    epgrppbcNONE,
    epgrppbcREFAT,
    epgrppbcCOS,
    epgrppbcPREVSTEPCOM
};

/*! \internal
 * \brief Pull group data used during pulling
 */
struct pull_group_work_t
{
    /*! \brief Constructor
     *
     * \param[in] params                  The group parameters set by the user
     * \param[in] atomSet                 The global to local atom set manager
     * \param[in] setPbcRefToPrevStepCOM Does this pull group use the COM from the previous step as reference position?
     * \param[in] maxNumThreads           Use either this number of threads of 1 for operations on x and f
     */
    pull_group_work_t(const t_pull_group& params,
                      gmx::LocalAtomSet   atomSet,
                      bool                setPbcRefToPrevStepCOM,
                      int                 maxNumThreads);

    //! Returns the number of threads to use for local atom operations based on the local atom count
    int numThreads() const
    {
        return atomSet.numAtomsLocal() <= c_pullMaxNumLocalAtomsSingleThreaded ? 1 : maxNumThreads;
    }

    /* Data only modified at initialization */
    const t_pull_group params;   /**< The pull group parameters */
    const int          epgrppbc; /**< The type of pbc for this pull group, see enum above */
    const int maxNumThreads; /**< The maximum number of threads to use for operations on x and f */
    bool      needToCalcCom; /**< Do we need to calculate the COM? (Not for group 0 or if only used as cylinder group) */
    std::vector<real> globalWeights; /**< Weights per atom set by the user and/or mass/friction coefficients, if empty all weights are equal */

    /* Data modified only at init or at domain decomposition */
    gmx::LocalAtomSet                  atomSet;      /**< Global to local atom set mapper */
    std::vector<real>                  localWeights; /**< Weights for the local atoms */
    std::unique_ptr<gmx::LocalAtomSet> pbcAtomSet;   /**< Keeps index of the pbc reference atom.
                                                          The stored LocalAtomSet consists of exactly   one atom when pbc reference atom is required.
                                                          When no pbc refence atom is used, this   pointer   shall be null. */

    /* Data, potentially, changed at every pull call */
    real mwscale; /**< mass*weight scaling factor 1/sum w m */
    real wscale;  /**< scaling factor for the weights: sum w m/sum w w m */
    real invtm;   /**< inverse total mass of the group: 1/wscale sum w m */
    std::vector<gmx::BasicVector<double>> mdw; /**< mass*gradient(weight) for atoms */
    std::vector<double>                   dv;  /**< distance to the other group(s) along vec */
    dvec                                  x;   /**< COM before update */
    dvec                                  xp;  /**< COM after update before constraining */
    dvec                                  x_prev_step; /**< center of mass of the previous step */
};

/* Struct describing the instantaneous spatial layout of a pull coordinate */
struct PullCoordSpatialData
{
    dvec   dr01;       /* The direction vector of group 1 relative to group 0 */
    dvec   dr23;       /* The direction vector of group 3 relative to group 2 */
    dvec   dr45;       /* The direction vector of group 5 relative to group 4 */
    dvec   vec;        /* The pull direction */
    double vec_len;    /* Length of vec for direction-relative */
    dvec   ffrad;      /* conversion factor from vec to radial force */
    double cyl_dev;    /* The deviation from the reference position */
    dvec   planevec_m; /* Normal of plane for groups 0, 1, 2, 3 for geometry dihedral */
    dvec   planevec_n; /* Normal of plane for groups 2, 3, 4, 5 for geometry dihedral */

    double value; /* The current value of the coordinate, units of nm or rad */
};

//! \brief Struct with parameters and force evaluation local data for a pull coordinate
struct pull_coord_work_t
{
    //! Constructor
    pull_coord_work_t(const t_pull_coord& params) :
        params(params),
        value_ref(0),
        spatialData(),
        scalarForce(0),
        bExternalPotentialProviderHasBeenRegistered(false),
        expressionParser(params.eGeom == PullGroupGeometry::Transformation ? params.expression : "",
                         params.coordIndex),
        transformationVariables(params.eGeom == PullGroupGeometry::Transformation ? params.coordIndex : 0)
    {
    }

    //! Pull coordinate parameters
    const t_pull_coord params;

    //! Dynamic pull group 0 for this coordinate with dynamic weights, only present when needed */
    std::unique_ptr<pull_group_work_t> dynamicGroup0;
    //! The reference value, usually init+rate*t, units of nm or rad.
    double value_ref;

    //! Data defining the current geometry
    PullCoordSpatialData spatialData;

    //! Scalar force for this cooordinate
    double scalarForce;

    //! For external-potential coordinates only, for checking if a provider has been registered
    bool bExternalPotentialProviderHasBeenRegistered;

    //! The expression parser for a transformation coordinate
    gmx::PullCoordExpressionParser expressionParser;
    //! Variables from other pull coordinates for a transformation coordinate
    std::vector<double> transformationVariables;
};

/* Struct for storing vectorial forces for a pull coordinate */
struct PullCoordVectorForces
{
    dvec force01; /* Force due to the pulling/constraining for groups 0, 1 */
    dvec force23; /* Force for groups 2 and 3 */
    dvec force45; /* Force for groups 4 and 5 */
};

/* Struct for sums over (local) atoms in a pull group */
struct ComSums
{
    /* For normal weighting */
    double sum_wm;   /* Sum of weight*mass        */
    double sum_wwm;  /* Sum of weight*weight*mass */
    dvec   sum_wmx;  /* Sum of weight*mass*x      */
    dvec   sum_wmxp; /* Sum of weight*mass*xp     */

    /* For cosine weighting */
    double sum_cm;  /* Sum of cos(x)*mass          */
    double sum_sm;  /* Sum of sin(x)*mass          */
    double sum_ccm; /* Sum of cos(x)*cos(x)*mass   */
    double sum_csm; /* Sum of cos(x)*sin(x)*mass   */
    double sum_ssm; /* Sum of sin(x)*sin(x)*mass   */
    double sum_cmp; /* Sum of cos(xp)*sin(xp)*mass */
    double sum_smp; /* Sum of sin(xp)*sin(xp)*mass */

    /* Dummy data to ensure adjacent elements in an array are separated
     * by a cache line size, max 128 bytes.
     * TODO: Replace this by some automated mechanism.
     */
    int dummy[32];
};

/*! \brief The normal COM buffer needs 3 elements per group */
static constexpr int c_comBufferStride = 3;

/*! \brief The cylinder buffer needs 9 elements per group */
static constexpr int c_cylinderBufferStride = 9;

struct pull_comm_t
{
    gmx_bool bParticipateAll; /* Do all ranks always participate in pulling? */
    gmx_bool bParticipate;    /* Does our rank participate in pulling? */
#if GMX_MPI
    MPI_Comm mpi_comm_com; /* Communicator for pulling */
#endif
    int nparticipate; /* The number of ranks participating */
    bool isMasterRank; /* Tells whether our rank is the master rank and thus should add the pull virial */

    int64_t setup_count; /* The number of decomposition calls */
    int64_t must_count;  /* The last count our rank needed to be part */

    /* Buffers for parallel reductions */
    std::vector<gmx::RVec>                pbcAtomBuffer; /* COM calculation buffer */
    std::vector<gmx::BasicVector<double>> comBuffer;     /* COM calculation buffer */
    std::vector<double> cylinderBuffer; /* cylinder ref. groups calculation buffer */
};

// The COM pull force calculation data structure
// TODO Convert this into a ForceProvider
struct pull_t
{
    /* Global parameters */
    pull_params_t params; /* The pull parameters, from inputrec */

    gmx_bool bPotential;  /* Are there coordinates with potential? */
    gmx_bool bConstraint; /* Are there constrained coordinates? */
    gmx_bool bAngle;      /* Are there angle geometry coordinates? */

    PbcType  pbcType;   /* the boundary conditions */
    int      npbcdim;   /* do pbc in dims 0 <= dim < npbcdim */
    gmx_bool bRefAt;    /* do we need reference atoms for a group COM ? */
    int      cosdim;    /* dimension for cosine weighting, -1 if none */
    gmx_bool bCylinder; /* Is group 0 a cylinder group? */

    /* Parameters + dynamic data for groups */
    std::vector<pull_group_work_t> group; /* The pull group param and work data */

    /* Parameters + dynamic data for coordinates */
    std::vector<pull_coord_work_t> coord; /* The pull group param and work data */

    /* Global dynamic data */
    gmx_bool bSetPBCatoms; /* Do we need to set x_pbc for the groups? */

    std::vector<ComSums> comSums; /* Work array for summing for COM, 1 entry per thread */

    pull_comm_t comm; /* Communication parameters, communicator and buffers */

    FILE* out_x; /* Output file for pull data */
    FILE* out_f; /* Output file for pull data */

    bool bXOutAverage; /* Output average pull coordinates */
    bool bFOutAverage; /* Output average pull forces */

    PullHistory* coordForceHistory; /* Pull coordinate and force history */

    /* The number of coordinates using an external potential */
    int numCoordinatesWithExternalPotential;
    /* Counter for checking external potential registration */
    int numUnregisteredExternalPotentials;
    /* */
    int numExternalPotentialsStillToBeAppliedThisStep;
};

/*! \brief Copies the pull group COM of the previous step from the checkpoint state to the pull state
 *
 * \param[in]   pull  The COM pull force calculation data structure
 * \param[in]   state The global state container
 */
void setPrevStepPullComFromState(struct pull_t* pull, const t_state* state);

/*! \brief Resizes the vector, in the state container, containing the COMs from the previous step
 *
 * \param[in]   state The global state container
 * \param[in]   pull  The COM pull force calculation data structure
 */
void allocStatePrevStepPullCom(t_state* state, const pull_t* pull);


#endif

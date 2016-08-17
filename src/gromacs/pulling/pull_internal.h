/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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
 *
 * \brief
 * This file contains datatypes and function declarations for internal
   use in the pull code.
 *
 * \author Berk Hess
 *
 * \inlibraryapi
 */

#ifndef GMX_PULLING_PULL_INTERNAL_H
#define GMX_PULLING_PULL_INTERNAL_H

#include "config.h"

#include "gromacs/mdtypes/pull-params.h"
#include "gromacs/utility/gmxmpi.h"

/*! \cond INTERNAL */

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

enum {
    epgrppbcNONE, epgrppbcREFAT, epgrppbcCOS
};

typedef struct
{
    t_pull_group  params;

    gmx_bool      bCalcCOM;   /* Calculate COM? Not if only used as cylinder group */
    int           epgrppbc;   /* The type of pbc for this pull group, see enum above */

    int           nat_loc;    /* Number of local pull atoms */
    int           nalloc_loc; /* Allocation size for ind_loc and weight_loc */
    int          *ind_loc;    /* Local pull indices */
    real         *weight_loc; /* Weights for the local indices */

    real          mwscale;    /* mass*weight scaling factor 1/sum w m */
    real          wscale;     /* scaling factor for the weights: sum w m/sum w w m */
    real          invtm;      /* inverse total mass of the group: 1/wscale sum w m */
    dvec         *mdw;        /* mass*gradient(weight) for atoms */
    double       *dv;         /* distance to the other group along vec */
    dvec          x;          /* center of mass before update */
    dvec          xp;         /* center of mass after update before constraining */
}
pull_group_work_t;

typedef struct
{
    t_pull_coord  params;     /* Pull coordinate (constant) parameters */

    double        value_ref;  /* The reference value, usually init+rate*t, units of nm or rad */
    double        value;      /* The current value of the coordinate, units of nm or rad */
    dvec          dr01;       /* The direction vector of group 1 relative to group 0 */
    dvec          dr23;       /* The direction vector of group 3 relative to group 2 */
    dvec          dr45;       /* The direction vector of group 5 relative to group 4 */
    dvec          vec;        /* The pull direction */
    double        vec_len;    /* Length of vec for direction-relative */
    dvec          ffrad;      /* conversion factor from vec to radial force */
    double        cyl_dev;    /* The deviation from the reference position */
    double        f_scal;     /* Scalar force for directional pulling */
    dvec          f01;        /* Force due to the pulling/constraining for groups 0, 1 */
    dvec          f23;        /* Force for groups 2 and 3 */
    dvec          f45;        /* Force for groups 4 and 5 */
    dvec          planevec_m; /* Normal of plane for groups 0, 1, 2, 3 for geometry dihedral */
    dvec          planevec_n; /* Normal of plane for groups 2, 3, 4, 5 for geometry dihedral */

    /* For external-potential coordinates only, for checking if a provider has been registered */
    bool          bExternalPotentialProviderHasBeenRegistered;
}
pull_coord_work_t;

/* Struct for sums over (local) atoms in a pull group */
struct pull_sum_com_t {
    /* For normal weighting */
    double sum_wm;    /* Sum of weight*mass        */
    double sum_wwm;   /* Sum of weight*weight*mass */
    dvec   sum_wmx;   /* Sum of weight*mass*x      */
    dvec   sum_wmxp;  /* Sum of weight*mass*xp     */

    /* For cosine weighting */
    double sum_cm;    /* Sum of cos(x)*mass          */
    double sum_sm;    /* Sum of sin(x)*mass          */
    double sum_ccm;   /* Sum of cos(x)*cos(x)*mass   */
    double sum_csm;   /* Sum of cos(x)*sin(x)*mass   */
    double sum_ssm;   /* Sum of sin(x)*sin(x)*mass   */
    double sum_cmp;   /* Sum of cos(xp)*sin(xp)*mass */
    double sum_smp;   /* Sum of sin(xp)*sin(xp)*mass */

    /* Dummy data to ensure adjacent elements in an array are separated
     * by a cache line size, max 128 bytes.
     * TODO: Replace this by some automated mechanism.
     */
    int    dummy[32];
};

typedef struct {
    gmx_bool    bParticipateAll; /* Do all ranks always participate in pulling? */
    gmx_bool    bParticipate;    /* Does our rank participate in pulling? */
#if GMX_MPI
    MPI_Comm    mpi_comm_com;    /* Communicator for pulling */
#endif
    int         nparticipate;    /* The number of ranks participating */

    gmx_int64_t setup_count;     /* The number of decomposition calls */
    gmx_int64_t must_count;      /* The last count our rank needed to be part */

    rvec       *rbuf;            /* COM calculation buffer */
    dvec       *dbuf;            /* COM calculation buffer */
    double     *dbuf_cyl;        /* cylinder ref. groups calculation buffer */
}
pull_comm_t;

struct pull_t
{
    pull_params_t      params;       /* The pull parameters, from inputrec */

    gmx_bool           bPotential;   /* Are there coordinates with potential? */
    gmx_bool           bConstraint;  /* Are there constrained coordinates? */
    gmx_bool           bAngle;       /* Are there angle geometry coordinates? */

    int                ePBC;         /* the boundary conditions */
    int                npbcdim;      /* do pbc in dims 0 <= dim < npbcdim */
    gmx_bool           bRefAt;       /* do we need reference atoms for a group COM ? */
    int                cosdim;       /* dimension for cosine weighting, -1 if none */

    int                ngroup;       /* Number of pull groups */
    int                ncoord;       /* Number of pull coordinates */
    pull_group_work_t *group;        /* The pull group param and work data */
    pull_group_work_t *dyna;         /* Dynamic groups for geom=cylinder */
    pull_coord_work_t *coord;        /* The pull group param and work data */

    gmx_bool           bCylinder;    /* Is group 0 a cylinder group? */

    gmx_bool           bSetPBCatoms; /* Do we need to set x_pbc for the groups? */

    int                nthreads;     /* Number of threads used by the pull code */
    pull_sum_com_t    *sum_com;      /* Work array for summing for COM, 1 entry per thread */

    pull_comm_t        comm;         /* Communication parameters, communicator and buffers */

    FILE              *out_x;        /* Output file for pull data */
    FILE              *out_f;        /* Output file for pull data */

    /* The number of coordinates using an external potential */
    int                numCoordinatesWithExternalPotential;
    /* Counter for checking external potential registration */
    int                numUnregisteredExternalPotentials;
    /* */
    int                numExternalPotentialsStillToBeAppliedThisStep;
};

/*! \endcond */

#endif

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
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
    atom_id      *ind_loc;    /* Local pull indices */
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
    t_pull_coord  params;    /* Pull coordinate (constant) parameters */

    double        value_ref; /* The reference value, usually init+rate*t */
    double        value;     /* The current value of the coordinate */
    dvec          dr;        /* The distance from the reference group */
    rvec          vec;       /* The pull direction */
    double        vec_len;   /* Length of vec for direction-relative */
    dvec          ffrad;     /* conversion factor from vec to radial force */
    double        cyl_dev;   /* The deviation from the reference position */
    double        f_scal;    /* Scalar force for directional pulling */
    dvec          f;         /* force due to the pulling/constraining */
}
pull_coord_work_t;

typedef struct {
    gmx_bool    bParticipateAll; /* Do all ranks always participate in pulling? */
    gmx_bool    bParticipate;    /* Does our rank participate in pulling? */
#ifdef GMX_MPI
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

    pull_comm_t        comm;         /* Communication parameters, communicator and buffers */

    FILE              *out_x;        /* Output file for pull data */
    FILE              *out_f;        /* Output file for pull data */
};

#ifdef __cplusplus
}
#endif

#endif

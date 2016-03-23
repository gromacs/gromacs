/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
 * This file contains datatypes for the mdp options used by the pull code.
 *
 * \author Berk Hess
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

#ifndef GMX_MDTYPES_PULL_PARAMS_H
#define GMX_MDTYPES_PULL_PARAMS_H

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

/*! \cond INTERNAL */

/*! \brief Struct that defines a pull group */
typedef struct {
    int   nat;     /**< Number of atoms in the pull group */
    int  *ind;     /**< The global atoms numbers */
    int   nweight; /**< The number of weights (0 or nat) */
    real *weight;  /**< Weights (use all 1 when weight==NULL) */
    int   pbcatom; /**< The reference atom for pbc (global number) */
} t_pull_group;

/*! Maximum number of pull groups that can be used in a pull coordinate */
static const int c_pullCoordNgroupMax = 6;

/*! \brief Struct that defines a pull coordinate */
typedef struct {
    int      eType;                        /**< The pull type: umbrella, constraint, ... */
    char    *externalPotentialProvider;    /**< Name of the module providing the external potential, only used with eType==epullEXTERNAL */
    int      eGeom;                        /**< The pull geometry */
    int      ngroup;                       /**< The number of groups, depends on eGeom */
    int      group[c_pullCoordNgroupMax];  /**< The pull groups: indices into the group arrays in pull_t and pull_params_t, ngroup indices are used */
    ivec     dim;                          /**< Used to select components for constraint */
    rvec     origin;                       /**< The origin for the absolute reference */
    rvec     vec;                          /**< The pull vector, direction or position */
    gmx_bool bStart;                       /**< Set init based on the initial structure */
    real     init;                         /**< Initial reference displacement (nm) or (deg) */
    real     rate;                         /**< Rate of motion (nm/ps) or (deg/ps) */
    real     k;                            /**< Force constant (kJ/(mol nm^2) or kJ/(mol rad^2) for umbrella pull type, or kJ/(mol nm) or kJ/(mol rad) for constant force pull type */
    real     kB;                           /**< Force constant for state B */
} t_pull_coord;

/*! \brief Struct containing all pull parameters */
typedef struct pull_params_t {
    int            ngroup;         /**< Number of pull groups */
    int            ncoord;         /**< Number of pull coordinates */
    real           cylinder_r;     /**< Radius of cylinder for dynamic COM (nm) */
    real           constr_tol;     /**< Absolute tolerance for constraints in (nm) */
    gmx_bool       bPrintCOM;      /**< Print coordinates of COM for each coord */
    gmx_bool       bPrintRefValue; /**< Print the reference value for each coord */
    gmx_bool       bPrintComp;     /**< Print cartesian components for each coord with geometry=distance */
    int            nstxout;        /**< Output interval for pull x */
    int            nstfout;        /**< Output interval for pull f */

    t_pull_group  *group;          /**< groups to pull/restrain/etc/ */
    t_pull_coord  *coord;          /**< the pull coordinates */
} pull_params_t;

/*! \endcond */

#endif /* GMX_MDTYPES_PULL_PARAMS_H */

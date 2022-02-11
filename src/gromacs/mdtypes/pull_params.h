/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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

#include <array>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

/*! \cond INTERNAL */

/*! \brief Struct that defines a pull group */
struct t_pull_group
{
    std::vector<int>  ind;     /**< The global atoms numbers */
    std::vector<real> weight;  /**< Weights (use all 1 when weight==NULL) */
    int               pbcatom; /**< The reference atom for pbc (global number) */
    int               pbcatom_input; /**< The reference atom for pbc (global number) as specified in the input parameters */
};

/*! Maximum number of pull groups that can be used in a pull coordinate */
static constexpr int c_pullCoordNgroupMax = 6;

/*! \brief Struct that defines a pull coordinate */
struct t_pull_coord
{
    //! The pull type: umbrella, constraint, ...
    PullingAlgorithm eType = PullingAlgorithm::Umbrella;
    //! Name of the module providing   the external potential, only used with eType==epullEXTERNAL
    std::string externalPotentialProvider;
    //! The pull geometry
    PullGroupGeometry eGeom = PullGroupGeometry::Distance;
    //! Mathematical expression evaluated by the pull code for transformation coordinates.
    std::string expression;
    //! The finite difference to use in numerical derivation of mathematical expressions
    double dx = 1e-9;
    //! The number of groups, depends on eGeom
    int ngroup = 0;
    /*! \brief The pull groups:
     *
     *  indices into the group arrays in pull_t and pull_params_t,
     *   ngroup indices are used
     */
    std::array<int, c_pullCoordNgroupMax> group;
    //! Used to select components for constraint
    gmx::IVec dim = { 0, 0, 0 };
    //! The origin for the absolute reference
    gmx::RVec origin = { 0, 0, 0 };
    //! The pull vector, direction or position
    gmx::RVec vec = { 0, 0, 0 };
    //! Set init based on the initial structure
    bool bStart = false;
    //! Initial reference displacement (nm) or (deg)
    real init = 0.0;
    //! Rate of motion (nm/ps) or (deg/ps)
    real rate = 0.0;
    /*! \brief Force constant
     *
     * For umbrella pull type this is (kJ/(mol nm^2) or kJ/(mol rad^2).
     * For constant force pull type it is kJ/(mol nm) or kJ/(mol rad).
     */
    real k = 0.0;
    //! Force constant for state B
    real kB = 0.0;
    //! The index of this coordinate in the list of coordinates
    int coordIndex = -1;
};

/*! \brief Struct containing all pull parameters */
struct pull_params_t
{
    //! Number of pull groups
    int ngroup = 0;
    //! Number of pull coordinates
    int ncoord = 0;
    //! Radius of cylinder for dynamic COM (nm)
    real cylinder_r = 0.0;
    //! Absolute tolerance for constraints in (nm)
    real constr_tol = 0.0;
    //! Print coordinates of COM for each coord
    bool bPrintCOM = false;
    //! Print the reference value for each coord
    bool bPrintRefValue = false;
    //! Print cartesian components for each coord with geometry=distance
    bool bPrintComp = false;
    //! Use the COM of each group from the previous step as reference
    bool bSetPbcRefToPrevStepCOM = false;
    //! Output interval for pull x
    int nstxout = 0;
    //! Output interval for pull f
    int nstfout = 0;
    //! Write the average coordinate during the output interval
    bool bXOutAverage = false;
    //! Write the average force during the output interval
    bool bFOutAverage = false;
    //! groups to pull/restrain/etc/
    std::vector<t_pull_group> group;
    //! the pull coordinates
    std::vector<t_pull_coord> coord;
};

/*! \endcond */

#endif /* GMX_MDTYPES_PULL_PARAMS_H */

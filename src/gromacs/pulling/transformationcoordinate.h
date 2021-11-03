/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 * \brief
 * Declares function for compute transformation coordinate values and forces
 *
 * \author Oliver Fleetwood <oliver.fleetwood@gmail.com>
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Berk Hess <hess@kth.se>
 *
 */
#ifndef GMX_PULL_TRANSFORMATIONCOORDINATE_H
#define GMX_PULL_TRANSFORMATIONCOORDINATE_H

struct pull_t;
struct pull_coord_work_t;

namespace gmx
{
template<typename>
class ArrayRef;

/*! \brief Calculates pull->coord[coord_ind].spatialData.value for a transformation pull coordinate
 *
 * This requires the values of the pull coordinates of lower indices to be set
 * \param[in] coord  The (transformation) coordinate to compute the value for
 * \param[in] variableCoords  Pull coordinates used as variables, entries 0 to coord->coordIndex
 *                            will be used
 * \returns Transformation value for pull coordinate.
 */
double getTransformationPullCoordinateValue(pull_coord_work_t*                coord,
                                            ArrayRef<const pull_coord_work_t> variableCoords);

/*! \brief Applies a force of a transformation pull coordinate and distributes it to pull coordinates of lower rank
 *
 * \param[in,out] pcrd            The transformation pull coordinate to act on
 * \param[in,out] variableCoords  List of variable coords up to the coord index of \p pcrd
 * \param[in] transformationCoordForce  The force working on coord \p pcrd
 */
void applyTransformationPullCoordForce(pull_coord_work_t*               pcrd,
                                       gmx::ArrayRef<pull_coord_work_t> variableCoords,
                                       double                           transformationCoordForce);

} // namespace gmx

#endif // GMX_PULL_TRANSFORMATIONCOORDINATE_H

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
/*! \internal \file
 * \brief Declares CP2K support information and version-dependent MPI requirements.
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_QMMMVERSION_H
#define GMX_APPLIED_FORCES_QMMMVERSION_H

#include <string>
#include <string_view>

namespace gmx
{

/*! \brief Returns a description of the CP2K QM/MM support.
 *
 * \return "disabled" if CP2K support is not compiled in, "enabled" if CP2K
 *         provides no version string, or "enabled (<version>)" containing
 *         the version string reported by CP2K.
 */
std::string qmmmDescription();

/*! \brief Determines whether the QM/MM support requires MPI_THREAD_MULTIPLE.
 *
 * Inspects the description without initializing CP2K or MPI.
 *
 * \param[in] description QM/MM support description returned by qmmmDescription().
 * \retval true For CP2K versions from 2027 onward, development versions of
 *              CP2K 2026.2, or unrecognized descriptions or versions, which
 *              conservatively require MPI_THREAD_MULTIPLE.
 * \retval false If CP2K support is disabled or the recognized version does
 *               not require MPI_THREAD_MULTIPLE.
 */
bool qmmmRequiresMpiThreadMultiple(std::string_view description);

} // namespace gmx

#endif // GMX_APPLIED_FORCES_QMMMVERSION_H

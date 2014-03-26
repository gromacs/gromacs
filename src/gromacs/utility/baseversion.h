/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 * \brief
 * Declares functions to get basic version information.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_BASEVERSION_H
#define GMX_UTILITY_BASEVERSION_H

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief
 * Version string, containing the version, date, and abbreviated hash.
 *
 * This can be a plain version if git version info was disabled during the
 * build.
 * The returned string starts with a literal word `VERSION`.
 *
 * \ingroup module_utility
 */
const char *gmx_version(void);
/*! \brief
 * Full git hash of the latest commit.
 *
 * If git version info was disabled during the build, returns an empty string.
 *
 * \ingroup module_utility
 */
const char *gmx_version_git_full_hash(void);
/*! \brief
 * Full git hash of the latest commit in a central \Gromacs repository.
 *
 * If git version info was disabled during the build, returns an empty string.
 * Also, if the latest commit was from a central repository, the return value
 * is an empty string.
 *
 * \ingroup module_utility
 */
const char *gmx_version_git_central_base_hash(void);

#ifdef __cplusplus
}
#endif

#endif

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
 * Declares functions for obtaining information about the operating environment
 * and the current process.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_SYSINFO_H
#define GMX_UTILITY_SYSINFO_H

#include <stddef.h>
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! \addtogroup module_utility
 * \{
 */

/*! \brief
 * Gets the hostname as given by gethostname(), if available.
 *
 * \param[out] buf  Buffer to receive the hostname.
 * \param[in]  len  Length of buffer \p buf (must be >= 8).
 * \returns 0 on success, -1 on error.
 *
 * If the value is not available, "unknown" is returned.
 * \p name should have at least size \p len.
 *
 * Does not throw.
 */
int gmx_gethostname(char *buf, size_t len);

/*! \brief
 * Returns the process ID of the current process.
 *
 * Does not throw.
 */
int gmx_getpid(void);
/*! \brief
 * Returns the current user ID, or -1 if not available.
 *
 * Does not throw.
 */
int gmx_getuid(void);
/*! \brief
 * Gets the current user name, if available.
 *
 * \param[out] buf  Buffer to receive the username.
 * \param[in]  len  Length of buffer \p buf (must be >= 8).
 * \returns 0 on success, -1 on error.
 *
 * Does not throw.
 */
int gmx_getusername(char *buf, size_t len);

/*! \brief
 * Portable version of ctime_r.
 *
 * Does not throw.
 */
char *gmx_ctime_r(const time_t *clock, char *buf, size_t len);
/*! \brief
 * Gets the current time as a string.
 *
 * \param[out] buf  Buffer to receive the string.
 * \param[in]  len  Length of buffer \p buf (26 characters should be sufficient).
 *
 * Does not throw.
 */
void gmx_format_current_time(char *buf, size_t len);

/*! \brief
 * Wrapper for nice().
 *
 * Does not throw.
 */
int gmx_set_nice(int level);

/*! \} */

#ifdef __cplusplus
}
#endif

#endif

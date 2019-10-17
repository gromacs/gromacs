/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief Declares functions that write timestamps to e.g. log files.
 *
 * \ingroup module_mdrunutility
 * \inlibraryapi
 */
#ifndef GMX_MDRUNUTILITY_PRINTTIME_H
#define GMX_MDRUNUTILITY_PRINTTIME_H

#include <cstdint>
#include <cstdio>
#include <cstdlib>

struct gmx_walltime_accounting;
struct t_commrec;
struct t_inputrec;

//! Print time to \c out.
void print_time(FILE*                    out,
                gmx_walltime_accounting* walltime_accounting,
                int64_t                  step,
                const t_inputrec*        ir,
                const t_commrec*         cr);

/*! \brief Print date, time, MPI rank and a description of this point
 * in time.
 *
 * \param[in] log       logfile, or NULL to suppress output
 * \param[in] rank      MPI rank to include in the output
 * \param[in] title     Description to include in the output
 * \param[in] the_time  Seconds since the epoch, e.g. as reported by gmx_gettime
 */
void print_date_and_time(FILE* log, int rank, const char* title, double the_time);

//! Print start time to \c fplog.
void print_start(FILE* fplog, const t_commrec* cr, gmx_walltime_accounting* walltime_accounting, const char* name);

#endif

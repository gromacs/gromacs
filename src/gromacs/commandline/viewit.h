/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
/*! \file
 * \brief
 * Provides function to open output files automatically (with some X11
 * programs).
 *
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_VIEWIT_H
#define GMX_COMMANDLINE_VIEWIT_H

struct gmx_output_env_t;
struct t_filenm;

/*! \brief
 * Executes an external (X11) command to view a file.
 *
 * Currently eps, xpm, xvg and pdb are supported.
 * Default programs are provided, can be overriden with environment vars
 * (but note that if the caller provides program-specific \p opts, setting the
 * environment variable most likely breaks things).
 */
void do_view(const gmx_output_env_t* oenv, const char* fn, const char* opts);

/*! \brief
 * Calls do_view() for all viewable output files.
 */
void view_all(const gmx_output_env_t* oenv, int nf, t_filenm fnm[]);

#endif

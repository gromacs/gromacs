/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
 *
 * \brief This file declares C-style entrypoints for mdrun
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_mdrun
 */
#ifndef GMX_PROGRAMS_MDRUN_MDRUN_H
#define GMX_PROGRAMS_MDRUN_MDRUN_H

#include "gromacs/utility/gmxmpi.h"

struct gmx_hw_info_t;

namespace gmx
{

/*! \brief Implements C-style main function for mdrun
 *
 * This implementation detects hardware itself, as suits
 * the gmx wrapper binary.
 *
 * \param[in]  argc          Number of C-style command-line arguments
 * \param[in]  argv          C-style command-line argument strings
 */
int gmx_mdrun(int argc, char* argv[]);

/*! \brief Implements C-style main function for mdrun
 *
 * This implementation facilitates reuse of infrastructure. This
 * includes the information about the hardware detected across the
 * given \c communicator. That suits e.g. efficient implementation of
 * test fixtures.
 *
 * \param[in]  communicator  The communicator to use for the simulation
 * \param[in]  hwinfo        Describes the hardware detected on the physical nodes of the communicator
 * \param[in]  argc          Number of C-style command-line arguments
 * \param[in]  argv          C-style command-line argument strings
 *
 * \todo Progress on https://gitlab.com/gromacs/gromacs/-/issues/3774
 * will remove the need of test binaries to call gmx_mdrun in a way
 * that is different from the command-line and gmxapi.
 */
int gmx_mdrun(MPI_Comm communicator, const gmx_hw_info_t& hwinfo, int argc, char* argv[]);

} // namespace gmx

#endif

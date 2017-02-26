/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Declares gmx::IMDOutputProvider.
 *
 * See \ref page_mdmodules for an overview of this and associated interfaces.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_IMDOUTPUTPROVIDER_H
#define GMX_MDTYPES_IMDOUTPUTPROVIDER_H

#include <cstdio>

struct gmx_output_env_t;
struct t_filenm;

namespace gmx
{

/*! \libinternal \brief
 * Interface for handling additional output files during a simulation.
 *
 * This interface provides a mechanism for additional modules to initialize
 * and finalize output files during the simulation.  Writing values to the
 * output files is currently handled elsewhere (e.g., when the module has
 * computed its forces).
 *
 * The interface is not very generic, as it has been written purely based on
 * extraction of existing functions related to electric field handling.
 * Also, the command-line parameters to specify the output files cannot be
 * specified by the module, but are hard-coded in mdrun.
 * This needs to be generalized when more modules are moved to use the
 * interface.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
class IMDOutputProvider
{
    public:
        /*! \brief
         * Initializes file output from a simulation run.
         *
         * \param[in] fplog File pointer for log messages
         * \param[in] nfile Number of files
         * \param[in] fnm   Array of filenames and properties
         * \param[in] bAppendFiles Whether or not we should append to files
         * \param[in] oenv  The output environment for xvg files
         */
        virtual void initOutput(FILE *fplog, int nfile, const t_filenm fnm[],
                                bool bAppendFiles, const gmx_output_env_t *oenv) = 0;

        //! Finalizes output from a simulation run.
        virtual void finishOutput() = 0;

    protected:
        ~IMDOutputProvider() {}
};

} // namespace gmx

#endif

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef GROMACS_SIMULATIONCONTEXT_H
#define GROMACS_SIMULATIONCONTEXT_H

/*! \libinternal
 * \file
 * \brief Provide ways for client code to own simulation resources.
 *
 * For `gmx mdrun` to be implemented as a client program, public API needs to
 * provide a way to create and manipulate the SimulationContext.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup module_mdrun
 * \inlibraryapi
 */

#include <memory>

#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdrun/mdfilenames.h"

struct t_filenm;
struct t_commrec;
struct gmx_output_env_t;

namespace gmx
{

class SimulationContext final
{
    public:
        /*!
         * \brief Object must be initialized with non-default constructor.
         */
        SimulationContext() = delete;

        /*!
         * \brief Initializate with borrowed values.
         *
         * \param communicator pointer to communicator handle (non-owning pointer to pointer).
         * \param hardwareOptions options to copy.
         * \param filenames preconfigured MD filename options container to copy.
         * \param outputEnvironment pointer to output environment handle (non-owning pointer to pointer).
         * \param logFile non-owning pointer to log file handle.
         *
         * \internal
         * SimulationContext should be the owner of these objects and these implementation
         * details are subject to change as ownership semantics are clarified in future
         * development.
         */
        SimulationContext(t_commrec        **communicator,
                          const gmx_hw_opt_t &hardwareOptions,
                          const MdFilenames  &filenames,
                          gmx_output_env_t **outputEnvironment,
                          FILE             **logFile);

        /*!
         * \brief         Non-owning communicator handle.
         *
         * \todo Context should own communicator.
         */
        t_commrec**           communicator_;

        //! \brief Parallelism information.
        gmx_hw_opt_t hardwareOptions_;

        //! filename options for simulation.
        MdFilenames filenames_;

        /*! \brief Non-owning handle to output environment.
         *
         * \todo Context should own output environment for client.
         */
        gmx_output_env_t**    outputEnvironment_;

        /*! \brief Non-owning handle to MD log file.
         *
         * \todo Context should own output facilities for client.
         */
        FILE**                logFile_;
};

}      // end namespace gmx

#endif //GROMACS_SIMULATIONCONTEXT_H

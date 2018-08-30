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
/*!\libinternal \file
 * \brief
 * Factory function to build OutputManager and helper object.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 * \inlibraryapi
 */
#ifndef GMX_COORDINATEIO_OUTPUTMANAGERBUILDER_H
#define GMX_COORDINATEIO_OUTPUTMANAGERBUILDER_H

#include <algorithm>
#include <utility>

#include "gromacs/coordinateio/outputmanager.h"

namespace gmx
{
/*! \brief
 * Factory function for OutputManager.
 *
 * Used to initialize a new instance of OutputManager with the user supplied information
 * for writing trajectory data to disk. Information needed is the file type, file name
 * corresponding to the type, if available topology information and selection information.
 *
 * If supplied, the modules contained within \p adapters are registered on the OutputManager
 * if possible.
 *
 * The factory function is responsible for the initial santity checks concerning file types and
 * availability of topology information, with the registration of modules being the second part.
 *
 * \param[in] mtop                   Pointer to full topology or null.
 * \param[in] sel                    Reference to global selection used to construct the object.
 * \param[in] filename               Name of new output file, used to deduce file type.
 * \param[in] adapters               Container for ICoordinateOutput modules that should be
 *                                   registered on the OutputManager during construction.
 * \throws    InconsistentInputError When user input and requirements don't match.
 */
OutputManagerPointer createOutputManager(const gmx_mtop_t        *mtop,
                                         const Selection         &sel,
                                         const std::string       &filename,
                                         CoordinateOutputAdapters adapters);


/*! \internal
 * \brief Helper struct for the factory function that builds a new OutputManager.
 *
 * Inherits from OutputManager to help build the object in the factory function of
 * OutputManagerBuilder. Because it inherits it can access the private constructor and
 * get things created. Shamelessly copied from the approach in make_constraints.h.
 */
struct OutputManager::OutputManagerBuildHelper : public OutputManager
{
    public:
        //! Constructor that can call the private OutputManager constructor.
        OutputManagerBuildHelper(std::string              name,
                                 int                      filetype,
                                 unsigned long            flag,
                                 const Selection         &sel,
                                 const gmx_mtop_t        *mtop,
                                 CoordinateOutputAdapters adapters) :
            OutputManager(std::move(name), filetype, flag, sel, mtop)
        {
            for (auto &adapter : adapters)
            {
                this->addOutputAdapter(std::move(adapter.module_));
            }
        }
};

} // namespace gmx

#endif

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
/*! \file
 * \brief
 * Declares gmx::OutputAdapterContainer, a storage object for
 * multiple outputadapters derived from the IOutputadaper interface.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_OUTPUTADAPTERCONTAINER_H
#define GMX_COORDINATEIO_OUTPUTADAPTERCONTAINER_H

#include <memory>
#include <vector>

#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/enumerationhelpers.h"

namespace gmx
{

/*! \libinternal \brief
 * Storage for output adapters that modify the state of a t_trxframe object.
 *
 * The OutputAdapterContainer is responsible for storing the number of
 * OutputAdapters, as well as the bitmask representing the current requirements
 * for constructing an CoordinateFile object with the modules registered. It is responsible
 * for ensuring that no module can be registered multiple times, and that the
 * correct order for some modifications is observed (e.g. we can not reduce the
 * number of coordinates written to a file before we have set all the other flags).
 * Any other behaviour indicates a programming error and triggers an assertion.
 *
 * The abilities that need to be cross checked for the container are usually constrained
 * by the file format the coordinate data will be written to. When declaring new abilities,
 * these must match the file type for the output.
 *
 * \todo Keeping track of those abilities has to be the responsibility of an object
 *       implementing and interface that declares it capabilities and will execute the
 *       the function of writing to a file.
 * \todo This could be changed to instead construct the container with a pointer to an
 *       ICoordinateOutputWriter that can be passed to the IOutputAdapter modules to check
 *       their cross-dependencies.
 */
class OutputAdapterContainer
{
public:
    //! Only allow constructing the container with defined output abilities.
    explicit OutputAdapterContainer(unsigned long abilities) : abilities_(abilities) {}
    //! Allow abilities to be also defined using the enum class.
    explicit OutputAdapterContainer(CoordinateFileFlags abilities) :
        abilities_(convertFlag(abilities))
    {
    }

    /*! \brief
     * Add an adapter of a type not previously added.
     *
     * Only one adapter of each type can be registered, and the order of adapters
     * is predefined in the underlying storage object.
     * Calls internal checks to make sure that the new adapter does not violate
     * any of the preconditions set to make an CoordinateFile object containing
     * the registered modules.
     *
     * \param[in] adapter unique_ptr to adapter, with container taking ownership here.
     * \param[in] type What kind of adapter is being added.
     * \throws InternalError When registering an adapter of a type already registered .
     * \throws InconsistentInputError When incompatible modules are added.
     */
    void addAdapter(OutputAdapterPointer adapter, CoordinateFileFlags type);

    //! Get vector of all registered adapters.
    ArrayRef<const OutputAdapterPointer> getAdapters() { return outputAdapters_; }
    //! Get info if we have any registered adapters.
    bool isEmpty() const;

private:
    //! Array of registered modules.
    EnumerationArray<CoordinateFileFlags, OutputAdapterPointer> outputAdapters_;
    //! Construction time bitmask declaring what the OutputManager can do.
    unsigned long abilities_ = convertFlag(CoordinateFileFlags::Base);
};

} // namespace gmx

#endif

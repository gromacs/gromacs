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
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

/*! \libinternal \brief
 * Storage for output adapters that modify the state of a t_trxframe object.
 *
 * The OutputAdapterContainer is responsible for storing the number of
 * OutputAdapters, as well as the bitmask representing the current requirements
 * for constructing an OutputManager with the modules registered. It is responsible
 * for ensuring that no module can be registered multiple times, and that the
 * correct order for some modifications is observed (e.g. we can not reduce the
 * number of coordinates written to a file before we have set all the other flags).
 * Any other behaviour indicates a programming error and leads to an assertion
 * being triggered.
 */
class OutputAdapterContainer
{
    public:
        explicit OutputAdapterContainer(unsigned long abilities) : abilities_(abilities)
        {}

        /*! \brief
         * Add new adapter to be processed.
         *
         * Adds \p adapter to the vector of modules being processed.
         * Calls internal checks to make sure that the new adapter does not violate
         * any of the preconditions set to make an OutputManager containing
         * the registered modules.
         *
         * \param[in] adapter unique_ptr to adapter, with container taking ownership here.
         * \throws InternalError When double registering or violating registration order.
         */
        void addAdapter(OutputAdapterPointer adapter);

        //! Get vector of all registered adapters.
        std::vector<OutputAdapterPointer> &getAdapters() { return outputAdapters_; }
        //! Get info if we have any registered adapters.
        bool isEmpty() { return outputAdapters_.empty(); }

    private:
        //! Vector of registered modules.
        std::vector<OutputAdapterPointer> outputAdapters_;
        //! Bitmask storing the local requirements.
        unsigned long                     requirements_ = 0;
        //! Bitmask storing the registered modules.
        unsigned long                     registered_ = 0;
        //! Construction time bitmask declaring what the OutputManager can do.
        unsigned long                     abilities_ = efBaseOutputManager;
};

} // namespace gmx

#endif

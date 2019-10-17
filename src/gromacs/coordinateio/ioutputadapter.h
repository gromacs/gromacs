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
 * Declares gmx::IOutputAdapter interface for modifying coordinate
 * file structures before writing them to disk.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_IOUTPUTADAPTER_H
#define GMX_COORDINATEIO_IOUTPUTADAPTER_H

#include <memory>

#include "gromacs/utility/classhelpers.h"

#include "enums.h"

struct t_trxframe;

namespace gmx
{

/*!\brief
 * OutputAdapter class for handling trajectory file flag setting and processing.
 *
 * This interface provides the base point upon which modules that modify trajectory frame
 * datastructures should be build. The interface itself does not provide any direct means
 * to modify the data, but only gives the virtual method to perform work on a
 * t_trxframe object. Classes that modify trajectory frames should implement this interface.
 *
 * \inlibraryapi
 * \ingroup module_coordinateio
 *
 */
class IOutputAdapter
{
public:
    /*! \brief
     * Default constructor for IOutputAdapter interface.
     */
    IOutputAdapter() {}
    virtual ~IOutputAdapter() {}
    //! Move constructor for old object.
    explicit IOutputAdapter(IOutputAdapter&& old) noexcept = default;

    /*! \brief
     * Change t_trxframe according to user input.
     *
     * \param[in] framenumber Frame number as reported from the
     *                        trajectoryanalysis framework or set by user.
     * \param[in,out] input   Pointer to trajectory analysis frame that will
     *                        be worked on.
     */
    virtual void processFrame(int framenumber, t_trxframe* input) = 0;

    /*! \brief
     * Checks that the abilities of the output writer are sufficient for this adapter.
     *
     * It can happen that a method to write coordinate files does not match with
     * a requested operation on the input data (e.g. the user requires velocities or
     * forces to be written to a PDB file).
     * To check those dependencies, derived classes need to implement a version of this
     * function to make sure that only matching methods can be used.
     *
     * \param[in] abilities The abilities of an output method that need to be checked against
     *                      the dependencies created by using the derived method.
     * \throws InconsistentInputError If dependencies can not be matched to abilities.
     */
    virtual void checkAbilityDependencies(unsigned long abilities) const = 0;

    GMX_DISALLOW_COPY_AND_ASSIGN(IOutputAdapter);
};

//! Smart pointer to manage the frame adapter object.
using OutputAdapterPointer = std::unique_ptr<IOutputAdapter>;

} // namespace gmx

#endif

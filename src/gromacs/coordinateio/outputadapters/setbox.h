/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Declares gmx::SetBox.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_SETBOX_H
#define GMX_COORDINATEIO_SETBOX_H

#include <memory>

#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/math/vec.h"
#include "gromacs/trajectory/trajectoryframe.h"

namespace gmx
{

/*!\brief
 * Allows changing box information when writing a coordinate file.
 *
 * \inpublicapi
 * \ingroup module_coordinateio
 */
class SetBox : public IOutputAdapter
{
public:
    /*! \brief
     * Construct SetBox object with a new user defined box.
     */
    explicit SetBox(const matrix box) { copy_mat(box, box_); }
    /*! \brief
     *  Move constructor for SetBox.
     */
    SetBox(SetBox&& old) noexcept { copy_mat(old.box_, box_); }

    ~SetBox() override { clear_mat(box_); }

    /*! \brief
     * Change coordinate frame information for output.
     *
     * In this case, box information is added to the \p t_trxframe object
     * depending on the user input.
     *
     * \param[in] input Coordinate frame to be modified later.
     */
    void processFrame(int /*framenumner*/, t_trxframe* input) override
    {
        copy_mat(box_, input->box);
    }

    void checkAbilityDependencies(unsigned long /* abilities */) const override {}

private:
    //! New box information from the user.
    matrix box_;
};

//! Smart pointer to manage the object.
using SetBoxPointer = std::unique_ptr<SetBox>;

} // namespace gmx

#endif

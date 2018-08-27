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
 * Declares gmx::SetBox.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_FILEIO_SETBOX_H
#define GMX_FILEIO_SETBOX_H

#include <algorithm>

#include "gromacs/coordinateio/coordinateoutput.h"
#include "gromacs/math/vec.h"

namespace gmx
{

/*!\brief
 * SetBox class allows changing writing of velocities to file.
 *
 * This class allows the user to define if velocities should be written
 * to the output coordinate file, and checks if they are available from the
 * currently processed data.
 *
 * \inpublicapi
 * \ingroup module_coordinateio
 */
class SetBox : public ICoordinateOutput
{
    public:
        /*! \brief
         * Construct SetBox object with choice for boolean value.
         *
         * Can be used to initialize SetBox from outside of trajectoryanalysis
         * with the user specified option to write coordinate velocities or not.
         */
        explicit SetBox(matrix box) : ICoordinateOutput(efAnyOutputSupported)
        {
            copy_mat(box, box_);
        }
        /*! \brief
         * Move constructor for SetBox.
         */
        SetBox &operator=(SetBox &&old)
        {
            copy_mat(old.box_, box_);
            return *this;
        }
        /*! \brief
         *  Move constructor for SetBox.
         */
        SetBox(SetBox &&old) : ICoordinateOutput(old.moduleFlags_)
        {
            copy_mat(old.box_, box_);
        }

        ~SetBox() {}

        /*! \brief
         * Change coordinate frame information for output.
         *
         * In this case, box information is added to the \p t_trxframe object
         * depending on the user input.
         *
         * \param[in] input Coordinate frame to be modified later.
         */
        virtual void processFrame(const int /*framenumner*/, t_trxframe *input);

    private:
        //! New box information from the user.
        matrix                            box_;
};

//! Smart pointer to manage the outputselector object.
typedef std::unique_ptr<SetBox>
    SetBoxPointer;

} // namespace gmx

#endif

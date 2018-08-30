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
 * Declares gmx::SetPrecision.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_FILEIO_SETPRECISION_H
#define GMX_FILEIO_SETPRECISION_H

#include <algorithm>

#include "gromacs/coordinateio/coordinateoutput.h"

namespace gmx
{

/*!\brief
 * SetPrecision class allows changing file writing precision.
 *
 * This class allows the user to define the precision for writing
 * coordinate data to output files.
 *
 * \inpublicapi
 * \ingroup module_coordinateio
 *
 */
class SetPrecision : public ICoordinateOutput
{
    public:
        /*! \brief
         * Construct SetPrecision object with choice for boolean value.
         *
         * Can be used to initialize SetPrecision from outside of trajectoryanalysis
         * with the user specified option to write coordinate velocities or not.
         * framework.
         */
        explicit SetPrecision(int precision) : ICoordinateOutput(efAnyOutputSupported),
                                               precision_(precision)
        {
        }
        /*! \brief
         * Move constructor for SetPrecision.
         */
        SetPrecision &operator=(SetPrecision &&old) noexcept
        {
            precision_ = old.precision_;
            return *this;
        }
        /*! \brief
         *  Move constructor for SetPrecision.
         */
        SetPrecision(SetPrecision &&old) noexcept : ICoordinateOutput(old.moduleFlags_),
                                                    precision_(old.precision_)
        {
        }

        ~SetPrecision() {}

        /*! \brief
         * Change coordinate frame information for output.
         *
         * In this case, the correct value for the coordinate precision
         * is applied for file writing.
         *
         * \param[in] input Coordinate frame to be modified later.
         */
        virtual void processFrame(int /*framenumber*/, t_trxframe *input);

    private:
        /*! \brief
         * Internal method to convert user provided precision value to internal value.
         *
         * \param[in] ndec User provided precision as interger.
         */
        real setFramePrecision(int ndec);
        //! User specified changes to default precision.
        int                             precision_;
};

//! Smart pointer to manage the outputselector object.
typedef std::unique_ptr<SetPrecision>
    SetPrecisionPointer;

} // namespace gmx

#endif

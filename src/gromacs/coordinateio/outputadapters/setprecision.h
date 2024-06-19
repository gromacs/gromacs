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
 * Declares gmx::SetPrecision.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inpublicapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_FILEIO_SETPRECISION_H
#define GMX_FILEIO_SETPRECISION_H

#include <algorithm>
#include <memory>

#include "gromacs/coordinateio/coordinatefileenums.h"
#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/utility/real.h"

struct t_trxframe;

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
class SetPrecision : public IOutputAdapter
{
public:
    /*! \brief
     * Construct SetPrecision object with user defined value.
     *
     * Can be used to initialize SetPrecision from outside of trajectoryanalysis
     * with the user specified option to change precision or not.
     *
     * \param[in] precision User defined value for output precision in file types that support it.
     */
    explicit SetPrecision(int precision) : precision_(precision)
    {
        // Only request special treatment if precision is not the default.
        if (precision == 3)
        {
            moduleRequirements_ = CoordinateFileFlags::Base;
        }
        else
        {
            moduleRequirements_ = CoordinateFileFlags::RequireChangedOutputPrecision;
        }
    }
    /*! \brief
     *  Move constructor for SetPrecision.
     */
    SetPrecision(SetPrecision&& old) noexcept = default;

    ~SetPrecision() override {}

    void processFrame(int /*framenumber*/, t_trxframe* input) override;

    void checkAbilityDependencies(unsigned long abilities) const override;

private:
    //! User specified changes to default precision.
    int precision_;
    //! Module requirements dependent on user input.
    CoordinateFileFlags moduleRequirements_;
};

//! Smart pointer to manage the outputselector object.
using SetPrecisionPointer = std::unique_ptr<SetPrecision>;

} // namespace gmx

#endif

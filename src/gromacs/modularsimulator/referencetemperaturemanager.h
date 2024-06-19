/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
/*! \internal \file
 * \brief Declares the a helper struct managing reference temperature changes
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_REFERENCETEMPERATUREMANAGER_H
#define GMX_MODULARSIMULATOR_REFERENCETEMPERATUREMANAGER_H

#include <functional>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "modularsimulatorinterfaces.h"

class gmx_ekindata_t;

namespace gmx
{

/*! \internal
 * \brief The algorithm changing the reference temperature
 *
 * In the legacy implementation, reference temperature changes by
 * different algorithms are not handled identically. This enum is
 * used to inform clients what algorithm caused the temperature
 * change, allowing them to customize their response.
 */
enum class ReferenceTemperatureChangeAlgorithm
{
};

/*! \internal
 * \brief Object managing reference temperature changes
 *
 * The ReferenceTemperatureManager allows to change the reference
 * temperatures of the temperature groups. Elements can register a callback
 * if they need to be informed about changes.
 *
 * The ReferenceTemperatureManager updates the inputrec. Elements
 * might, however, have a copy of the reference temperature they
 * need updated, or perform another action upon change of the
 * reference temperature (e.g. velocity scaling or recalculating
 * a temperature coupling integral).
 */
class ReferenceTemperatureManager final
{
public:
    //! Constructor
    ReferenceTemperatureManager(gmx_ekindata_t* ekindata);
    //! Register a callback for reference temperature update
    void registerUpdateCallback(ReferenceTemperatureCallback referenceTemperatureCallback);
    //! Set reference temperatures (one per temperature group)
    void setReferenceTemperature(ArrayRef<const real>                newReferenceTemperatures,
                                 ReferenceTemperatureChangeAlgorithm algorithm);

private:
    //! List of callbacks
    std::vector<ReferenceTemperatureCallback> callbacks_;
    //! Pointer to the kinetic energy data
    gmx_ekindata_t* ekindata_;
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_REFERENCETEMPERATUREMANAGER_H

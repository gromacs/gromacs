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
/*! \internal \file
 * \brief Declares the composite element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */
#ifndef GROMACS_MDTYPES_COMPOSITESIMULATORELEMENT_H
#define GROMACS_MDTYPES_COMPOSITESIMULATORELEMENT_H

#include <memory>
#include <vector>

#include "gromacs/compat/pointers.h"

#include "modularsimulatorinterfaces.h"

namespace gmx
{

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Composite simulator element
 *
 * The composite simulator element takes a call list of elements and implements
 * the ISimulatorElement interface, making a group of elements effectively
 * behave as one. This simplifies building algorithms.
 *
 * The CompositeSimulatorElement can optionally also own the elements, but does
 * not require this. The owner of a CompositeSimulatorElement object can hence
 * decide to either pass the ownership to CompositeSimulatorElement, or keep
 * the ownership (and guarantee that they remain valid during the life time
 * of the CompositeSimulatorElement object). CompositeSimulatorElement will only
 * call the setup and teardown methods on the owned elements, thereby avoiding
 * to call them more than once. Consequently, the owner of the elements not
 * owned by CompositeSimulatorElement is responsible to call setup and teardown
 * methods on these elements.
 */
class CompositeSimulatorElement final : public ISimulatorElement
{
public:
    //! Constructor
    explicit CompositeSimulatorElement(std::vector<compat::not_null<ISimulatorElement*>> elementCallList,
                                       std::vector<std::unique_ptr<ISimulatorElement>> elements,
                                       int                                             frequency);

    /*! \brief Register run function for step / time
     *
     * Lets every member of the composite simulator register run functions
     * for the given step.
     *
     * \param step                 The step number
     * \param time                 The time
     * \param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction) override;

    /*! \brief Element setup
     *
     * Calls the setup functions of the single elements.
     */
    void elementSetup() override;

    /*! \brief Element teardown
     *
     * Calls the teardown functions of the single elements.
     */
    void elementTeardown() override;

private:
    //! The call list of elements forming the composite element
    std::vector<compat::not_null<ISimulatorElement*>> elementCallList_;
    //! List of elements owned by composite element
    std::vector<std::unique_ptr<ISimulatorElement>> elementOwnershipList_;
    //! The frequency at which the composite element is running
    const int frequency_;
};

} // namespace gmx

#endif // GROMACS_MDTYPES_COMPOSITESIMULATORELEMENT_H

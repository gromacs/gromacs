/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

#ifndef GROMACS_RESTRAINTMDMODULE_H
#define GROMACS_RESTRAINTMDMODULE_H

/*! \libinternal \file
 * \brief Library interface for RestraintMDModule
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_restraint
 */

#include <memory>
#include <vector>

#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/restraint/restraintpotential.h"

namespace gmx
{

// Forward declaration to allow opaque pointer to library internal class.
class RestraintMDModuleImpl;
struct MDModulesNotifiers;
class ForceProviders;
class IMDOutputProvider;
class IMdpOptionProvider;

/*! \libinternal \ingroup module_restraint
 * \brief MDModule wrapper for Restraint implementations.
 *
 * Shares ownership of an object implementing the IRestraintPotential interface.
 * Provides the IMDModule interfaces.
 */
class RestraintMDModule final : public gmx::IMDModule
{
public:
    RestraintMDModule() = delete;

    /*!
     * \brief Constructor used by static create() method.
     */
    explicit RestraintMDModule(std::unique_ptr<RestraintMDModuleImpl> restraint);

    ~RestraintMDModule() override;

    /*!
     * \brief Wrap a restraint potential as an MDModule
     *
     * Consumers of the interfaces provided by an IMDModule do not extend the lifetime
     * of the interface objects returned by mdpOptionProvider(), outputProvider(), or
     * registered via initForceProviders(). Calling code must keep this object alive
     * as long as those interfaces are needed (probably the duration of an MD run).
     *
     * \param restraint shared ownership of an object for calculating restraint forces
     * \param sites list of sites for the framework to pass to the restraint
     * \return new wrapper object sharing ownership of restraint
     */
    static std::unique_ptr<RestraintMDModule> create(std::shared_ptr<gmx::IRestraintPotential> restraint,
                                                     const std::vector<int>& sites);

    /*!
     * \brief Implement IMDModule interface
     *
     * Unused.
     *
     * \return nullptr.
     */
    IMdpOptionProvider* mdpOptionProvider() override;

    /*!
     * \brief Implement IMDModule interface
     *
     * Unused.
     *
     * \return nullptr.
     */
    IMDOutputProvider* outputProvider() override;

    /*!
     * \brief Implement IMDModule interface.
     *
     * See gmx::IMDModule::initForceProviders()
     * \param forceProviders manager in the force record.
     */
    void initForceProviders(ForceProviders* forceProviders) override;

    //! Subscribe to simulation setup notifications
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifiers) override;
    //! Subscribe to pre processing notifications
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* notifiers) override;

private:
    /*!
     * \brief Private implementation opaque pointer.
     */
    std::unique_ptr<RestraintMDModuleImpl> impl_;
};

} // end namespace gmx

#endif // GROMACS_RESTRAINTMDMODULE_H

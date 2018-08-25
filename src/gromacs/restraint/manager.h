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

#ifndef GROMACS_RESTRAINT_MANAGER_H
#define GROMACS_RESTRAINT_MANAGER_H

/*! \libinternal \file
 * \brief Declare the Manager for restraint potentials.
 *
 * \inlibraryapi
 * \ingroup module_restraint
 */

#include <memory>
#include <mutex>
#include <string>

#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/basedefinitions.h"

struct t_commrec;
struct t_mdatoms;
struct pull_t;

namespace gmx
{
class LegacyPuller;

/*! \libinternal
 * \brief Implementation details for MD restraints
 *
 * \ingroup module_restraint
 */
namespace restraint
{

class ManagerImpl;
class ICalculation;

/*! \libinternal \ingroup module_restraint
 * \brief Manage the Restraint potentials available for Molecular Dynamics.
 *
 * Until further factoring of the MD integrators and force calculations, we use a singleton
 * to reduce coupling between rapidly changing GROMACS components. Ultimately, this manager
 * should either not be necessary or can be used in more tightly scoped instances.
 *
 * The manager takes ownership of the "pull groups" (or atomic selections) and of
 * the various restraints and constraints applied for a given simulation.
 *
 * Calling code provides the manager with a means to access the various required input data
 * to be used when restraints are computed.
 *
 * \internal
 * When calculate(t) is called, the various input and output data sources are provided to
 * a CalculationBuilder to produce a Calculation object for the point in time, t.
 * Constructing the Calculation object triggers updates to the simulation state force array
 * and virial tensor. After construction, the Calculation object can be queried for calculated
 * data such as energy or pulling work.
 */
class Manager final
{
    public:
        ~Manager();

        /// Get a shared reference to the global manager.
        static std::shared_ptr<Manager> instance();

        Manager(const Manager&) = delete;
        Manager              &operator=(const Manager &) = delete;
        Manager(Manager &&) = delete;
        Manager              &operator=(Manager &&) = delete;

        /*!
         * \brief Clear registered restraints and reset the manager.
         */
        void clear() noexcept;

        /*!
         * \brief Get the number of currently managed restraints.
         *
         * \return number of restraints.
         *
         * \internal
         * Only considers the IRestraintPotential objects
         */
        unsigned long countRestraints() noexcept;

        /*! \brief Obtain the ability to create a restraint MDModule
         *
         * Though the name is reminiscent of the evolving idea of a work specification, the
         * Spec here is just a list of restraint modules.
         *
         * \param puller shared ownership of a restraint potential interface.
         * \param name key by which to reference the restraint.
         */
        void addToSpec(std::shared_ptr<gmx::IRestraintPotential> puller,
                       std::string                               name);

        /*!
         * \brief Get a copy of the current set of restraints to be applied.
         *
         * \return a copy of the list of restraint potentials.
         */
        std::vector < std::shared_ptr < IRestraintPotential>> getSpec() const;

    private:
        /// Private constructor enforces singleton pattern.
        Manager();

        /// Regulate initialization of the manager instance when the singleton is first accessed.
        static std::mutex initializationMutex_;

        /// Ownership of the shared reference to the global manager.
        static std::shared_ptr<Manager> instance_;

        /// Opaque implementation pointer.
        std::unique_ptr<ManagerImpl> impl_;
};


}      // end namespace restraint

}      // end namespace gmx

#endif //GROMACS_RESTRAINT_MANAGER_H

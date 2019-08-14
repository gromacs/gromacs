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
/*! \internal \file
 * \brief
 * Declares data structure and utilities for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "densityfitting.h"

#include "gromacs/mdrunutility/mdmodulenotification.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

#include "densityfittingforceprovider.h"
#include "densityfittingoptions.h"
#include "densityfittingoutputprovider.h"


namespace gmx
{

class IMdpOptionProvider;
class DensityFittingForceProvider;

namespace
{

/*! \internal
 * \brief Density fitting
 *
 * Class that implements the density fitting forces and potential
 * \note the virial calculation is not yet implemented
 */
class DensityFitting final : public IMDModule
{
    public:
        /*! \brief Construct the density fitting module.
         *
         * \param[in] notifier allows the module to subscribe to notifications from MdModules.
         *
         * The density fitting code subscribes to these notifications:
         *   - setting atom group indices in the densityFittingOptions_ by
         *     taking a parmeter const IndexGroupsAndNames &
         *   - storing its internal parameters in a tpr file by writing to a
         *     key-value-tree during pre-processing by a function taking a
         *     KeyValueTreeObjectBuilder as parameter
         *   - reading its internal parameters from a key-value-tree during
         *     simulation setup by taking a const KeyValueTreeObject & parameter
         */
        explicit DensityFitting(MdModulesNotifier *notifier)
        {
            // Callbacks for several kinds of MdModuleNotification are created
            // and subscribed, and will be dispatched correctly at run time
            // based on the type of the parameter required by the lambda.

            // Setting atom group indices
            const auto setFitGroupIndicesFunction = [this](const IndexGroupsAndNames &indexGroupsAndNames) {
                    densityFittingOptions_.setFitGroupIndices(indexGroupsAndNames);
                };
            notifier->notifier_.subscribe(setFitGroupIndicesFunction);

            // Writing internal parameters during pre-processing
            const auto writeInternalParametersFunction = [this](KeyValueTreeObjectBuilder treeBuilder) {
                    densityFittingOptions_.writeInternalParametersToKvt(treeBuilder);
                };
            notifier->notifier_.subscribe(writeInternalParametersFunction);

            // Reading internal parameters during simulation setup
            const auto readInternalParametersFunction = [this](const KeyValueTreeObject &tree) {
                    densityFittingOptions_.readInternalParametersFromKvt(tree);
                };
            notifier->notifier_.subscribe(readInternalParametersFunction);
        }

        //! From IMDModule; this class provides the mdpOptions itself
        IMdpOptionProvider *mdpOptionProvider() override { return &densityFittingOptions_; }

        //! Add this module to the force providers if active
        void initForceProviders(ForceProviders *forceProviders) override
        {
            if (densityFittingOptions_.active())
            {
                const auto &parameters = densityFittingOptions_.buildParameters();
                forceProvider_ = std::make_unique<DensityFittingForceProvider>(parameters);
                forceProviders->addForceProvider(forceProvider_.get());
            }
        }

        //! This MDModule provides its own output
        IMDOutputProvider *outputProvider() override { return &densityFittingOutputProvider_; }

    private:
        //! The output provider
        DensityFittingOutputProvider                 densityFittingOutputProvider_;
        //! The options provided for density fitting
        DensityFittingOptions                        densityFittingOptions_;
        //! Object that evaluates the forces
        std::unique_ptr<DensityFittingForceProvider> forceProvider_;
};

}   // namespace

std::unique_ptr<IMDModule> DensityFittingModuleInfo::create(MdModulesNotifier * notifier)
{
    return std::make_unique<DensityFitting>(notifier);
}

const std::string DensityFittingModuleInfo::name_ = "density-guided-simulation";

} // namespace gmx

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

#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/stringutil.h"

#include "densityfittingforceprovider.h"
#include "densityfittingparameters.h"

namespace gmx
{

namespace
{

/*! \internal
 * \brief Input data storage for density fitting
 */
class DensityFittingOptions : public IMdpOptionProvider
{
    public:
        DensityFittingOptions()
        { }

        //! From IMdpOptionProvider
        void initMdpTransform(IKeyValueTreeTransformRules * /*transform*/ ) override
        {}

        /*! \brief
         * Build mdp parameters for density fitting to be output after pre-processing.
         * \param[in, out] builder the builder for the mdp options output KV-tree.
         * \note This should be symmetrical to option initialization without
         *       employing manual prefixing with the section name string once
         *       the legacy code blocking this design is removed.
         */
        void buildMdpOutput(KeyValueTreeObjectBuilder *builder) const override
        {
            builder->addValue<bool>(inputSectionName_ + "-" + c_activeTag_,
                                    active_);
        }

        /*! \brief
         * Connect option name and data.
         */
        void initMdpOptions(IOptionsContainerWithSections *options) override
        {
            auto section = options->addSection(OptionSection(inputSectionName_.c_str()));
            section.addOption(
                    BooleanOption(c_activeTag_.c_str()).store(&active_));
        }

        //! Report if this set of options is active
        bool active() const
        {
            return active_;
        }

        //! Process input options to parameters, including input file reading.
        DensityFittingParameters buildParameters()
        {
            // read map to mrc-header and data vector
            // convert mrcheader and data vector to grid data
            return {};
        }

    private:
        //! The name of the density-fitting module
        const std::string     inputSectionName_ = "density-guided-simulation";

        const std::string     c_activeTag_             = "active";
        bool                  active_                  = false;
};

/*! \internal
 * \brief Density fitting
 *
 * Class that implements the density fitting forces and potential
 * \note the virial calculation is not yet implemented
 */
class DensityFitting final : public IMDModule,
                             public IMDOutputProvider
{
    public:
        DensityFitting() = default;

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
        IMDOutputProvider *outputProvider() override { return this; }

        //! Initialize output
        void initOutput(FILE * /*fplog*/, int /*nfile*/, const t_filenm /*fnm*/[],
                        bool /*bAppendFiles*/, const gmx_output_env_t * /*oenv*/) override
        {}

        //! Finalizes output from a simulation run.
        void finishOutput() override {}

    private:
        //! The options provided for density fitting
        DensityFittingOptions densityFittingOptions_;
        //! Object that evaluates the forces
        std::unique_ptr<DensityFittingForceProvider> forceProvider_;
};

}   // namespace

std::unique_ptr<IMDModule> createDensityFittingModule()
{
    return std::unique_ptr<IMDModule>(new DensityFitting());
}

} // namespace gmx

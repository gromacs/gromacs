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


namespace gmx
{

namespace
{

/*! \internal
 * \brief Input data storage for density fitting
 */
class DensityFittingData
{
    public:
        DensityFittingData() : forceConstant_ {0},
        referenceMapSource_ {"map.mrc"}
        { }

        /*! \brief
         * Connect option name and data.
         */
        void initMdpOptions(IOptionsContainerWithSections *options)
        {
            options->addOption(
                    StringOption(mdpReferenceMapSourceTag_c.c_str()).store(&referenceMapSource_));
            options->addOption(
                    DoubleOption(mdpForceConstantTag_c.c_str()).store(&forceConstant_));
        }
        /*! \brief
         * Creates mdp parameters for density fitting
         */
        void buildMdpOutput(KeyValueTreeObjectBuilder *builder) const
        {
            builder->addValue<std::string>("density-fitting-" +
                                           mdpReferenceMapSourceTag_c, referenceMapSource_);
            builder->addValue<double>("density-fitting-" +
                                      mdpForceConstantTag_c, forceConstant_);
        }

        //! Return the force constant for density fitting
        double forceConstant() const
        {
            return forceConstant_;
        }

    private:
        const std::string mdpReferenceMapSourceTag_c = "reference-map-source";
        const std::string mdpForceConstantTag_c      = "force-constant";
        double            forceConstant_;
        std::string       referenceMapSource_;
};

/*! \internal
 * \brief Density fitting
 *
 * Class that implements the density fitting forces
 */
class DensityFitting final : public IMDModule,
                             public IMdpOptionProvider,
                             public IMDOutputProvider,
                             public IForceProvider
{
    public:
        DensityFitting() = default;

        // From IMDModule; this class provides the mdpOptions itself
        IMdpOptionProvider *mdpOptionProvider() override { return this; }

        //! Add this module to the force providers if active
        void initForceProviders(ForceProviders *forceProviders) override
        {
            if (isActive())
            {
                forceProviders->addForceProvider(this);
            }
        }

        //! From IMdpOptionProvider
        void initMdpTransform(IKeyValueTreeTransformRules * /*transform*/ ) override
        {}

        //! Initialize the mdp options
        void initMdpOptions(IOptionsContainerWithSections *options) override
        {
            auto section = options->addSection(OptionSection("density-fitting"));
            densityFittingData_.initMdpOptions(&section);
        }

        //! Generate the mdp file output for file that is written after pre-processsing the input mdp
        void buildMdpOutput(KeyValueTreeObjectBuilder *builder) const override
        {
            const std::vector<std::string> comment = {
                "",
                "; Density fitting ",
                "; Apply additional forces to a simulation that maximise similarity of",
                "; a simulated three-dimensional density map with a given reference map.",
                "; "
            };
            builder->addValue<std::string>("comment-density-fitting", joinStrings(comment, "\n"));
            densityFittingData_.buildMdpOutput(builder);
        }

        //! Calculate the density fitting forces
        void calculateForces(const ForceProviderInput & /*forceProviderInput*/,
                             ForceProviderOutput      * /*forceProviderOutput*/) override
        {}

        // This MDModule provides its own output
        IMDOutputProvider *outputProvider() override { return this; }

        //! Initialize output
        void initOutput(FILE * /*fplog*/, int /*nfile*/, const t_filenm /*fnm*/[],
                        bool /*bAppendFiles*/, const gmx_output_env_t * /*oenv*/) override
        {}

        //! Finalizes output from a simulation run.
        void finishOutput() override {}

    private:

        //! DensityFitting is active if the force constant for the fitting is not zero
        bool isActive() const
        {
            return (densityFittingData_.forceConstant() != 0);
        }
        //! the data needed to evaluate forces for density fitting
        DensityFittingData densityFittingData_;
};

}   // namespace

std::unique_ptr<IMDModule> createDensityFittingModule()
{
    return std::unique_ptr<IMDModule>(new DensityFitting());
}

} // namespace gmx

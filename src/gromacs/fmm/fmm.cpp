/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 *
 * \brief
 * Interface of \Gromacs to the fast multipole method (FMM).
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \ingroup module_fmm
 *
 */

#include "gmxpre.h"

#include "fmm.h"
#ifdef GMX_WITH_FMM
#include "gromacs/fmm/fmm-impl.h"
#endif

#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/strconvert.h"

struct FmmImpl;

namespace gmx
{


/*! \brief
 *  Class providing basic infrastructure allowing \Gromacs to interface a FMM.
 *
 *  This class handles FMM-related .mdp input file parameters, including forwarding
 *  them to the underlying FMM implementation class. The actual force evaluation call
 *  is also simply passed on to the FMM implementation.
 */
class Fmm final : public IMDModule,
                  public IMdpOptionProvider,
                  public IForceProvider
{
    public:
        Fmm() { };

        // From IMDModule
        IMdpOptionProvider *mdpOptionProvider() override { return this; }
        IMDOutputProvider  *outputProvider()    override { return nullptr; }
        void initForceProviders(gmx_unused ForceProviders *forceProviders, gmx_unused const t_inputrec *ir, gmx_unused const gmx_mtop_t *mtop) override
        {
            if (EEL_FMM(ir->coulombtype))
            {
#ifdef GMX_WITH_FMM
                this->impl_ = createFmmImpl(ir, &fmmParameters, mtop);
                forceProviders->addForceProvider(this);
#else
                gmx_fatal(FARGS, "FMM electrostatics requested, but GROMACS was compiled without FMM.\n"
                          "Set -DCMAKE_CXX_FLAGS=-DGMX_WITH_FMM when building GROMACS.");
#endif
            }
        }

        // From IMdpOptionProvider
        void initMdpTransform(IKeyValueTreeTransformRules *transform) override;
        void initMdpOptions(IOptionsContainerWithSections *options) override;

        // From IForceProvider
        //! \copydoc IForceProvider::calculateForces()
        void calculateForces(const t_commrec *cr,
                             const t_mdatoms *mdatoms,
                             const matrix     box,
                             double           t,
                             const rvec      *x,
                             ArrayRef<RVec>   force,
                             gmx_enerdata_t  *enerd) override;

    private:
        fmmInputParameters fmmParameters;

#ifdef GMX_WITH_FMM
        // Pointer to the underlying FMM implementation (could for example be the ExaFmm or the Juelich FmSolvr library):
        std::unique_ptr<FmmImpl> impl_;
#endif
};


std::unique_ptr<IMDModule> createFastMultipoleModule()
{
    return std::unique_ptr<IMDModule>(new Fmm());
}


void Fmm::initMdpTransform(IKeyValueTreeTransformRules *rules)
{
    std::string prefix = "/";
    prefix.append(nameOfFmmOptionsSection);

    rules->addRule().from<std::string>("/fmm-precision"      ).to<real>(prefix + "/fmm-precision"      ).transformWith(&fromStdString<real>);
    rules->addRule().from<std::string>("/fmm-multipole-order").to<int> (prefix + "/fmm-multipole-order").transformWith(&fromStdString<int>);
    rules->addRule().from<std::string>("/fmm-tree-depth"     ).to<int> (prefix + "/fmm-tree-depth"     ).transformWith(&fromStdString<int>);
    rules->addRule().from<std::string>("/fmm-separation"     ).to<int> (prefix + "/fmm-separation"     ).transformWith(&fromStdString<int>);
};

void Fmm::initMdpOptions(IOptionsContainerWithSections *options)
{
    // Create a section for FMM input
    auto section = options->addSection(OptionSection(nameOfFmmOptionsSection));

    // Create the FMM-related options
    section.addOption(RealOption   ("fmm-precision"      ).store(&fmmParameters.precision));
    section.addOption(IntegerOption("fmm-multipole-order").store(&fmmParameters.order));
    section.addOption(IntegerOption("fmm-tree-depth"     ).store(&fmmParameters.depth));
    section.addOption(IntegerOption("fmm-separation"     ).store(&fmmParameters.separation));
}


void Fmm::calculateForces(gmx_unused const t_commrec *cr,
                          gmx_unused const t_mdatoms *mdatoms,
                          gmx_unused const matrix     box,
                          gmx_unused double           t,
                          gmx_unused const rvec      *x,
                          gmx_unused ArrayRef<RVec>   force,
                          gmx_unused gmx_enerdata_t  *enerd)
{
#ifdef GMX_WITH_FMM
    double coulombEnergy = 0.0;
    impl_->calculateForcesAndEnergies(cr, box, mdatoms->chargeA, x, force, &coulombEnergy);

    // For now, add the FMM Coulomb energy to group with index 0
    enerd->grpp.ener[egCOULSR][0] += coulombEnergy;
#endif
}

} // namespace gmx

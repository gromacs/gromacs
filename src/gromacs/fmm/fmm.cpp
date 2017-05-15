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
#include "fmsolver.h"
#endif

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/strconvert.h"

namespace gmx
{

const char *nameOfFmmOptionsSection = "fast-multipole-method";

/*! \brief
 *  Identification string that can be prepended to FMM-related output
 */
const char *fmmStr                  = "FMM:";


class FmmSolverImpl;

/*! \brief
 *  Class providing functionality for Coulomb force evaluations with the FMM
 */
class Fmm : public IMDModule,
            public IMdpOptionProvider,
            public IMDOutputProvider,
            public IForceProvider
{
    public:
        Fmm() : depth_(-1), order_(-1), separation_(-1), precision_(0.001) { };

        // From IMDModule
        virtual IMdpOptionProvider *mdpOptionProvider() override { return this; }
        virtual IMDOutputProvider  *outputProvider()    override { return this; }
        virtual IForceProvider     *forceProvider()     override { return this; }

        // From IMdpOptionProvider
        virtual void initMdpTransform(IKeyValueTreeTransformRules *transform) override;
        virtual void initMdpOptions(IOptionsContainerWithSections *options) override;

        // From IMDOutputProvider
        virtual void initOutput(FILE*, int, const t_filenm [], bool, const gmx_output_env_t*) override { };
        virtual void finishOutput() override { };

        // From IForceProvider
        virtual void initForcerec(gmx_unused t_forcerec *fr) override;

        //! \copydoc IForceProvider::calculateForces()
        virtual void calculateForces(gmx_unused const t_commrec  * /* cr      */,
                                     gmx_unused const t_mdatoms  * /* mdatoms */,
                                     gmx_unused PaddedRVecVector * /* force   */,
                                     gmx_unused double             /* t       */) override { }; // not used

        /*! \brief Main FMM routine, here the Coulomb forces and energies are calculated.
         *
         * Should be called at every time step.
         *
         * \param[in]  cr            Communication record, keeps information needed for parallel runs.
         * \param[in]  box           The simulation box vectors.
         * \param[in]  q             Array of charges.
         * \param[in]  x             Array of the atomic positions.
         * \param[out] f_fmm         Array of resulting Coulomb forces on the atoms.
         * \param[out] coulombEnergy The Coulomb energy of the charge configuration.
         */
        virtual void calculateForces(const t_commrec *cr,
                                     const matrix    &box,
                                     const real      *q,
                                     const rvec      *x,
                                     rvec            *f_fmm,
                                     double          *coulombEnergy) override;

        /*! \brief
         * Initialize the Fast Multipole Solver.
         *
         * \param[in] mtop   The topology information.
         * \param[in] hw_opt Threading and GPU options that were set automatically or by the user.
         * \param[in] useGpu Use GPU acceleration?
         */
        virtual void initialize(const gmx_mtop_t   *mtop,
                                const gmx_hw_opt_t *hw_opt,
                                gmx_bool            useGpu) override;

    private:
        int  depth_;      //! tree depth
        int  order_;      //! order of the multipole expansion
        int  separation_; //! separation criterion
        real precision_;  //! Max. accepted relative error in Coulomb potential energy when using FMM

#ifdef GMX_WITH_FMM
        // Pointer to the FMM solver implementation:
        std::unique_ptr<FmmSolverImpl> fmsolvr_;
#endif
};

#ifdef GMX_WITH_FMM
std::unique_ptr<FmmSolverImpl> createFmSolvr(t_forcerec *fr);
#endif


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
    options->addOption(RealOption   ("fmm-precision"      ).store(&precision_));
    options->addOption(IntegerOption("fmm-multipole-order").store(&order_));
    options->addOption(IntegerOption("fmm-tree-depth"     ).store(&depth_));
    options->addOption(IntegerOption("fmm-separation"     ).store(&separation_));
}


void Fmm::initForcerec(gmx_unused t_forcerec *fr)
{
#ifdef GMX_WITH_FMM
    if (fr->eeltype == eelFMM)
    {
        fr->fmm  = this;
        fmsolvr_ = createFmSolvr(fr);
    }
#endif
}

void Fmm::initialize(gmx_unused const gmx_mtop_t *mtop, gmx_unused const gmx_hw_opt_t *hw_opt, gmx_unused gmx_bool useGpu)
{
#ifdef GMX_WITH_FMM
    fmsolvr_->n_                  = mtop->natoms;
    fmsolvr_->precision_          = precision_;
    fmsolvr_->override_d_         = depth_;
    fmsolvr_->override_p_         = order_;
    fmsolvr_->override_ws_        = separation_;
    fmsolvr_->useGpuAcceleration_ = useGpu;
    fmsolvr_->hw_opt_             = hw_opt;

    // Generate a naive list of atom IDs and excluded IDs.
    fmsolvr_->makeExclusionList(mtop);
#endif
}


void Fmm::calculateForces(gmx_unused const t_commrec *cr,
                          gmx_unused const matrix    &box,
                          gmx_unused const real      *q,
                          gmx_unused const rvec      *x,
                          gmx_unused rvec            *f_fmm,
                          gmx_unused double          *coulombEnergy)
{
#ifdef GMX_WITH_FMM
    fprintf(stderr, "%s Entering main FMM routine on rank %d\n", fmmStr, cr->sim_nodeid);

    if (DOMAINDECOMP(cr))
    {
        fprintf(stderr, "%s rank %d, local atoms %d\n", fmmStr, cr->sim_nodeid, cr->dd->nat_home);
    }
    fmsolvr_->calculateForcesAndEnergies(cr, box, q, x, f_fmm, coulombEnergy);
#endif
}

} // namespace gmx

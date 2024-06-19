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
 * \brief Implements test of 1-4 interactions
 *
 * This test is copied from the bonded interactions test and slightly
 * modified since 'do_pairs' takes a different set of arguments than
 * 'calculateSimpleBond'. To keep the test setup uncluttered this test is
 * therefore not merged into the bonded test but implemented standalone.
 *
 * The test setup consists of 2 atom pairs that are tested in an fep setting
 * (vanishing charge and lennard-jones parameters of one atom) and without
 * fep. Placement of the atoms in the box is such that shift-forces and pbc
 * paths in do_pairs are covered.
 *
 * \author Sebastian Kehl <sebastian.kehl@mpcdf.mpg.de>
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include "gromacs/listed_forces/pairs.h"

#include <cmath>

#include <filesystem>
#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/listed_forces/bonded.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/booltype.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

//! Number of atoms used in these tests.
constexpr int c_numAtoms = 3;

/*! \brief Output from pairs kernels
 *
 */
struct OutputQuantities
{
    OutputQuantities(int energyGroup) :
        energy(energyGroup), dvdLambda(static_cast<int>(FreeEnergyPerturbationCouplingType::Count), 0.0)
    {
    }

    //! Energy of this interaction
    gmx_grppairener_t energy;
    //! Derivative with respect to lambda
    std::vector<real> dvdLambda;
    //! Shift vectors
    rvec fShift[gmx::detail::c_numIvecs] = { { 0 } };
    //! Forces
    alignas(GMX_REAL_MAX_SIMD_WIDTH * sizeof(real)) rvec4 f[c_numAtoms] = { { 0 } };
};

/*! \brief Utility to check the output from pairs tests
 *
 * \param[in] checker Reference checker
 * \param[in] output  The output from the test to check
 * \param[in] bondedKernelFlavor  Flavor for determining what output to check
 * \param[in] functionType type of the interaction
 */
void checkOutput(TestReferenceChecker*    checker,
                 const OutputQuantities&  output,
                 const BondedKernelFlavor bondedKernelFlavor,
                 const int                functionType)
{
    if (computeEnergy(bondedKernelFlavor))
    {
        switch (functionType)
        {
            case F_LJ14:
            case F_LJC14_Q:
                checker->checkReal(output.energy.energyGroupPairTerms[NonBondedEnergyTerms::Coulomb14][0],
                                   "Epot Coulomb14");
                checker->checkReal(output.energy.energyGroupPairTerms[NonBondedEnergyTerms::LJ14][0],
                                   "Epot LJ14");
                break;
            case F_LJC_PAIRS_NB:
                checker->checkReal(output.energy.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR][0],
                                   "Epot Coulomb14");
                checker->checkReal(output.energy.energyGroupPairTerms[NonBondedEnergyTerms::LJSR][0],
                                   "Epot LJ14");
                break;
            default: gmx_fatal(FARGS, "Unknown function type %d in do_nonbonded14", functionType);
        }
        checker->checkReal(output.dvdLambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)],
                           "dVdlCoul ");
        checker->checkReal(output.dvdLambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)],
                           "dVdlVdw ");
    }
    checker->checkSequence(std::begin(output.f), std::end(output.f), "Forces");
}

/* \brief Utility class to setup forcerec and interaction parameters
 *
 * Data is only initialized as necessary for the 1-4 interactions to work!
 */
class ForcerecHelper
{
public:
    ForcerecHelper()
    {
        fepVals_.sc_alpha                = 0.3;
        fepVals_.sc_power                = 1;
        fepVals_.sc_r_power              = 6.0;
        fepVals_.sc_sigma                = 0.3;
        fepVals_.sc_sigma_min            = 0.3;
        fepVals_.bScCoul                 = true;
        fepVals_.scGapsysScaleLinpointLJ = 0.85;
        fepVals_.scGapsysScaleLinpointQ  = 0.3;
        fepVals_.scGapsysSigmaLJ         = 0.3;
        fepVals_.softcoreFunction        = SoftcoreType::Beutler;

        fr_.ic = std::make_unique<interaction_const_t>();
        // set data in ic
        fr_.ic->softCoreParameters = std::make_unique<interaction_const_t::SoftCoreParameters>(fepVals_);

        // set data in fr
        real tableRange = 2.9;
        fr_.pairsTable = make_tables(nullptr, fr_.ic.get(), nullptr, tableRange, GMX_MAKETABLES_14ONLY);
        fr_.efep       = haveFep_;
        fr_.fudgeQQ          = 0.5;
        fr_.use_simd_kernels = useSimd_;
        fr_.bMolPBC          = haveMolPBC_;
    }

    t_forcerec* get() { return &fr_; }

private:
    FreeEnergyPerturbationType haveFep_    = FreeEnergyPerturbationType::No;
    bool                       useSimd_    = false;
    bool                       haveMolPBC_ = false;

    t_lambda   fepVals_;
    t_forcerec fr_;
};

ForcerecHelper frHelper;


/*! \brief Input structure for listed forces tests
 */
struct ListInput
{
public:
    //! Function type
    int fType = -1;
    //! do fep
    bool fep = false;
    //! Tolerance for float evaluation
    float floatToler = 1e-6;
    //! Tolerance for double evaluation
    double doubleToler = 1e-8;
    //! Interaction parameters
    t_iparams iparams = { { 0 } };

    //! Constructor
    ListInput() {}

    /*! \brief Constructor with tolerance
     *
     * \param[in] ftol Single precision tolerance
     * \param[in] dtol Double precision tolerance
     */
    ListInput(float ftol, double dtol)
    {
        floatToler  = ftol;
        doubleToler = dtol;
    }

    /*! \brief Set parameters for lj14 interaction
     *
     * Fep is used if either c6A != c6B or c12A != c12B.
     *
     * \param[in] c6A  lj-c6 of state A
     * \param[in] c12A lj-c12 of state A
     * \param[in] c6B  lj-c6 of state B
     * \param[in] c12B lj-c12 of state B
     */
    ListInput setLj14Interaction(real c6A, real c12A, real c6B, real c12B)
    {
        fType             = F_LJ14;
        fep               = (c6A != c6B || c12A != c12B);
        iparams.lj14.c6A  = c6A;
        iparams.lj14.c12A = c12A;
        iparams.lj14.c6B  = c6B;
        iparams.lj14.c12B = c12B;

        return *this;
    }

    /*! \brief Set parameters for ljc14 interaction
     *
     * \param[in] qi     charge i
     * \param[in] qj     charge j
     * \param[in] c6     lj-c6
     * \param[in] c12    lj-c12
     * \param[in] fudgeQ fudge factor
     */
    ListInput setLjc14Interaction(real qi, real qj, real c6, real c12, real fudgeQ)
    {
        fType             = F_LJC14_Q;
        iparams.ljc14.qi  = qi;
        iparams.ljc14.qj  = qj;
        iparams.ljc14.c6  = c6;
        iparams.ljc14.c12 = c12;
        iparams.ljc14.fqq = fudgeQ;

        return *this;
    }

    /*! \brief Set parameters for ljcnb interaction
     *
     * \param[in] qi  charge i
     * \param[in] qj  charge j
     * \param[in] c6  lj-c6
     * \param[in] c12 lj-c12
     */
    ListInput setLjcnbInteraction(real qi, real qj, real c6, real c12)
    {
        fType             = F_LJC_PAIRS_NB;
        iparams.ljcnb.qi  = qi;
        iparams.ljcnb.qj  = qj;
        iparams.ljcnb.c6  = c6;
        iparams.ljcnb.c12 = c12;

        return *this;
    }
};

class ListedForcesPairsTest :
    public ::testing::TestWithParam<std::tuple<ListInput, PaddedVector<RVec>, PbcType>>
{
protected:
    matrix               box_;
    t_pbc                pbc_;
    PaddedVector<RVec>   x_;
    PbcType              pbcType_;
    ListInput            input_;
    TestReferenceData    refData_;
    TestReferenceChecker checker_;

    ListedForcesPairsTest() : checker_(refData_.rootChecker())
    {
        input_   = std::get<0>(GetParam());
        x_       = std::get<1>(GetParam());
        pbcType_ = std::get<2>(GetParam());
        clear_mat(box_);
        box_[0][0] = box_[1][1] = box_[2][2] = 1.0;
        set_pbc(&pbc_, pbcType_, box_);

        FloatingPointTolerance tolerance = relativeToleranceAsPrecisionDependentFloatingPoint(
                1.0, input_.floatToler, input_.doubleToler);

        checker_.setDefaultTolerance(tolerance);
    }

    void testOneIfunc(TestReferenceChecker* checker, const real lambda, const SoftcoreType softcoreType)
    {
        SCOPED_TRACE(std::string("Testing PBC type: ") + c_pbcTypeNames[pbcType_]);

        // 'definition of pairs' is a concatenation of #npairs (here 2)
        // 'nAtomsPerPair+1'-tuples (fType a_0 a_i ... a_nAtomsPerPair)
        std::vector<t_iatom> iatoms = { 0, 1, 2, 0, 0, 2 };

        std::vector<int>            ddgatindex = { 0, 1, 2 };
        std::vector<real>           chargeA    = { 1.0, -0.5, -0.5 };
        std::vector<real>           chargeB    = { 0.0, 0.0, 0.0 };
        std::vector<BoolType>       perturbed  = { true, true, true };
        std::vector<unsigned short> egrp       = { 0, 0, 0 };
        t_mdatoms                   mdatoms    = { 0 };

        mdatoms.chargeA    = chargeA;
        mdatoms.chargeB    = chargeB;
        mdatoms.bPerturbed = perturbed;
        mdatoms.cENER      = egrp;
        mdatoms.nPerturbed = 3;

        t_forcerec* fr = frHelper.get();
        fr->efep = input_.fep ? FreeEnergyPerturbationType::Yes : FreeEnergyPerturbationType::No;
        fr->ic->softCoreParameters->softcoreType = softcoreType;
        if (pbcType_ != PbcType::No)
        {
            fr->bMolPBC = true;
        }

        std::vector<BondedKernelFlavor> flavors = { BondedKernelFlavor::ForcesAndVirialAndEnergy };

        if (!input_.fep || lambda == 0)
        {
            fr->use_simd_kernels = true;
            flavors.push_back(BondedKernelFlavor::ForcesSimdWhenAvailable);
        }

        for (const auto flavor : flavors)
        {
            SCOPED_TRACE("Testing bonded kernel flavor: " + c_bondedKernelFlavorStrings[flavor]);

            StepWorkload stepWork;
            stepWork.computeVirial = computeVirial(flavor);
            stepWork.computeEnergy = computeEnergy(flavor);

            bool havePerturbedInteractions = input_.fep;
            if (flavor == BondedKernelFlavor::ForcesSimdWhenAvailable)
            {
                havePerturbedInteractions = false;
            }

            int numEnergyTerms      = static_cast<int>(NonBondedEnergyTerms::Count);
            int numFepCouplingTerms = static_cast<int>(FreeEnergyPerturbationCouplingType::Count);
            OutputQuantities  output(numEnergyTerms);
            std::vector<real> lambdas(numFepCouplingTerms, lambda);

            do_pairs(input_.fType,
                     iatoms.size(),
                     iatoms.data(),
                     &input_.iparams,
                     as_rvec_array(x_.data()),
                     output.f,
                     output.fShift,
                     &pbc_,
                     lambdas.data(),
                     output.dvdLambda.data(),
                     mdatoms.chargeA,
                     mdatoms.chargeB,
                     makeArrayRef(mdatoms.bPerturbed),
                     mdatoms.cENER,
                     mdatoms.nPerturbed,
                     fr,
                     havePerturbedInteractions,
                     stepWork,
                     &output.energy,
                     ddgatindex.data());

            checkOutput(checker, output, flavor, input_.fType);
            auto shiftForcesChecker = checker->checkCompound("Shift-Forces", "Shift-forces");

            if (computeVirial(flavor))
            {
                shiftForcesChecker.checkVector(output.fShift[gmx::c_centralShiftIndex], "Central");
            }
            else
            {
                // Permit omitting to compare shift forces with
                // reference data when that is useless.
                shiftForcesChecker.disableUnusedEntriesCheck();
            }
        }
    }

    void testIfunc()
    {
        TestReferenceChecker thisChecker =
                checker_.checkCompound("FunctionType", interaction_function[input_.fType].name)
                        .checkCompound("FEP", (input_.fep ? "Yes" : "No"));

        if (input_.fep)
        {
            const int numLambdas = 3;
            for (int i = 0; i < numLambdas; ++i)
            {
                const real lambda = i / (numLambdas - 1.0);
                for (SoftcoreType c : EnumerationWrapper<SoftcoreType>{})
                {
                    auto lambdaChecker = thisChecker.checkCompound("Lambda", toString(lambda));
                    auto softcoreChecker = lambdaChecker.checkCompound("Sofcore", enumValueToString(c));
                    testOneIfunc(&softcoreChecker, lambda, c);
                }
            }
        }
        else
        {
            testOneIfunc(&thisChecker, 0.0, SoftcoreType::Beutler);
        }
    }
};

TEST_P(ListedForcesPairsTest, Ifunc)
{
    testIfunc();
}

//! Function types for testing 1-4 interaction. Add new terms at the end.
std::vector<ListInput> c_14Interaction = {
    { ListInput(1e-5, 1e-7).setLj14Interaction(0.001458, 1.0062882e-6, 0.0, 0.0) },
    { ListInput(1e-5, 1e-7).setLj14Interaction(0.001458, 1.0062882e-6, 0.001458, 1.0062882e-6) },
    { ListInput(1e-5, 1e-7).setLjc14Interaction(1.0, -1.0, 0.001458, 1.0062882e-6, 0.5) },
    { ListInput(1e-5, 1e-7).setLjcnbInteraction(1.0, -1.0, 0.001458, 1.0062882e-6) }
};

//! PBC values for testing
std::vector<PbcType> c_pbcForTests = { PbcType::No, PbcType::XY, PbcType::Xyz };

/*! \brief Coordinates for testing 1-4 interaction
 *
 * Define coordinates for 3 atoms here, which will be used in 2 interactions.
 */
std::vector<PaddedVector<RVec>> c_coordinatesFor14Interaction = {
    { { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, { 1.1, 1.2, 1.3 } }
};

INSTANTIATE_TEST_SUITE_P(14Interaction,
                         ListedForcesPairsTest,
                         ::testing::Combine(::testing::ValuesIn(c_14Interaction),
                                            ::testing::ValuesIn(c_coordinatesFor14Interaction),
                                            ::testing::ValuesIn(c_pbcForTests)));

} // namespace

} // namespace test

} // namespace gmx

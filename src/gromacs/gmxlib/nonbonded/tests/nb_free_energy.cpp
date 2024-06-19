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
 * \brief Implements test of nonbonded fep kernel
 *
 * Implements the test logic from the bonded interactions also for the
 * nonbonded fep kernel. This requires setting up some more input
 * structures that in the bonded case.
 *
 * The test setup consists of an atom pair that is evaluated in an fep setting
 * (vanishing charge and lennard-jones parameters of atom #2) with and without
 * softcore Potentials.
 *
 * \author Sebastian Kehl <sebastian.kehl@mpcdf.mpg.de>
 * \ingroup module_gmxlib_nonbonded
 */
#include "gmxpre.h"

#include "gromacs/gmxlib/nonbonded/nb_free_energy.h"

#include "config.h"

#include <cmath>

#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
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
constexpr int c_numAtoms     = 4;
constexpr int c_numAtomTypes = 3;

/*! \brief Output from nonbonded fep kernel
 *
 */
struct OutputQuantities
{
    OutputQuantities() :
        energy(static_cast<int>(NonBondedEnergyTerms::Count)),
        dvdLambda(static_cast<int>(FreeEnergyPerturbationCouplingType::Count), 0.0),
        fShift(1, { 0.0, 0.0, 0.0 }),
        f(c_numAtoms, { 0.0, 0.0, 0.0 })
    {
    }

    //! Energies of this interaction (size EgNR)
    gmx_grppairener_t energy;
    //! Derivative with respect to lambda (size efptNR)
    std::vector<real> dvdLambda;
    //! Shift force vectors (size N_IVEC but in this test only 1)
    std::vector<RVec> fShift;
    //! Forces (size c_numAtoms)
    PaddedVector<RVec> f;
};

/*! \brief Utility to check the output from nonbonded test
 *
 * \param[in] checker Reference checker
 * \param[in] output  The output from the test to check
 */
void checkOutput(TestReferenceChecker* checker, const OutputQuantities& output)
{
    checker->checkReal(output.energy.energyGroupPairTerms[NonBondedEnergyTerms::LJSR][0], "EVdw ");
    checker->checkReal(output.energy.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR][0],
                       "ECoul ");
    checker->checkReal(output.dvdLambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)],
                       "dVdlCoul ");
    checker->checkReal(output.dvdLambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)],
                       "dVdlVdw ");

    checker->checkSequence(std::begin(output.f), std::end(output.f), "Forces");

    auto shiftForcesChecker = checker->checkCompound("Shift-Forces", "Shift-forces");
    shiftForcesChecker.checkVector(output.fShift[0], "Central");
}

class InteractionConstHelper
{
public:
    InteractionConstHelper() {}

    //! init data to construct interaction_const
    void initInteractionConst(CoulombInteractionType coulType, VanDerWaalsType vdwType, InteractionModifiers vdwMod)
    {
        coulType_ = coulType;
        vdwType_  = vdwType;
        vdwMod_   = vdwMod;

        // initialize correction tables
        interaction_const_t tmp;
        tmp.ewaldcoeff_q       = calc_ewaldcoeff_q(1.0, 1.0e-5);
        coulEwaldCoeff_        = tmp.ewaldcoeff_q;
        tmp.ewaldcoeff_lj      = calc_ewaldcoeff_lj(1.0, 1.0e-5);
        vdwEwaldCoeff_         = tmp.ewaldcoeff_lj;
        tmp.eeltype            = coulType;
        tmp.vdwtype            = vdwType;
        tmp.coulombEwaldTables = std::make_unique<EwaldCorrectionTables>();
        tmp.vdwEwaldTables     = std::make_unique<EwaldCorrectionTables>();

        init_interaction_const_tables(nullptr, &tmp, 1.0, 0.0);
        coulombTables_ = *tmp.coulombEwaldTables;
        vdwTables_     = *tmp.vdwEwaldTables;
    }

    /*! \brief Setup interaction_const_t
     *
     * \param[in]  fepVals t_lambda struct of fep values
     * \param[out] ic      interaction_const_t pointer with data
     */
    void getInteractionConst(const t_lambda& fepVals, interaction_const_t* ic)
    {
        ic->softCoreParameters = std::unique_ptr<interaction_const_t::SoftCoreParameters>(
                new interaction_const_t::SoftCoreParameters(fepVals));

        ic->coulombEwaldTables  = std::unique_ptr<EwaldCorrectionTables>(new EwaldCorrectionTables);
        *ic->coulombEwaldTables = coulombTables_;

        ic->vdwEwaldTables  = std::unique_ptr<EwaldCorrectionTables>(new EwaldCorrectionTables);
        *ic->vdwEwaldTables = vdwTables_;

        // set coulomb and vdw types
        ic->eeltype      = coulType_;
        ic->vdwtype      = vdwType_;
        ic->vdw_modifier = vdwMod_;

        // some non default parameters used in this testcase
        ic->epsfac                   = gmx::c_one4PiEps0 * 0.25;
        ic->reactionFieldCoefficient = 0.0; // former k_rf
        ic->reactionFieldShift       = 1.0; // former c_rf
        ic->ewaldcoeff_q             = coulEwaldCoeff_;
        ic->ewaldcoeff_lj            = vdwEwaldCoeff_;
        ic->sh_ewald                 = 1.0e-5;
        ic->sh_lj_ewald              = -1.0;
        ic->dispersion_shift.cpot    = -1.0;
        ic->repulsion_shift.cpot     = -1.0;
    }

private:
    //! correction tables
    EwaldCorrectionTables coulombTables_;
    EwaldCorrectionTables vdwTables_;

    //! coulomb and vdw type specifiers
    CoulombInteractionType coulType_;
    real                   coulEwaldCoeff_;
    VanDerWaalsType        vdwType_;
    InteractionModifiers   vdwMod_;
    real                   vdwEwaldCoeff_;
};


/* \brief Utility class to setup forcerec
 *
 * This helper takes care of handling the neccessary data that are kept
 * at various places and in various forms in the forcerec hierarchy such
 * that this class can safely be used.
 *
 * Data is only initialized as necessary for the nonbonded kernel to work!
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
    }

    //! initialize data structure to construct forcerec
    void initForcerec(const gmx_ffparams_t&  idef,
                      CoulombInteractionType coulType,
                      VanDerWaalsType        vdwType,
                      InteractionModifiers   vdwMod)
    {
        icHelper_.initInteractionConst(coulType, vdwType, vdwMod);
        nbfp_        = makeNonBondedParameterLists(idef.atnr, idef.iparams, false);
        ljPmeC6Grid_ = makeLJPmeC6GridCorrectionParameters(idef.atnr, idef.iparams, LongRangeVdW::Geom);
    }

    void setSoftcoreAlpha(const real scBeutlerAlphaOrGapsysLinpointScaling)
    {
        fepVals_.sc_alpha                = scBeutlerAlphaOrGapsysLinpointScaling;
        fepVals_.scGapsysScaleLinpointLJ = scBeutlerAlphaOrGapsysLinpointScaling;
        fepVals_.scGapsysScaleLinpointQ  = scBeutlerAlphaOrGapsysLinpointScaling;
    }
    void setSoftcoreCoulomb(const bool scCoulomb) { fepVals_.bScCoul = scCoulomb; }
    void setSoftcoreType(const SoftcoreType softcoreType)
    {
        fepVals_.softcoreFunction = softcoreType;
    }


    //! get forcerec data as wanted by the nonbonded kernel
    void getForcerec(t_forcerec* fr)
    {
        fr->ic = std::make_unique<interaction_const_t>();

        // set data in ic
        icHelper_.getInteractionConst(fepVals_, fr->ic.get());

        // set data in fr
        fr->ljpme_c6grid = ljPmeC6Grid_;
        fr->nbfp         = nbfp_;
        fr->shift_vec    = { { 0.0, 0.0, 0.0 } };
        fr->rlist        = rlist_;
        fr->ntype        = c_numAtomTypes;

        // simd
        fr->use_simd_kernels = GMX_USE_SIMD_KERNELS;
    }

private:
    InteractionConstHelper icHelper_;
    std::vector<real>      ljPmeC6Grid_;
    std::vector<real>      nbfp_;
    t_lambda               fepVals_;
    real                   rlist_ = 1.0;
};

/*! \brief Utility structure to hold atoms data
 *
 * A system having 4 interactions with the following perturbation pattern:
 *  - no perturbation
 *  - vdw- and coulomb-perturbation
 *  - coulomb-perturbation only
 *  - vdw-perturbation only
 *
 * This is realized by defining 3 different atom types that control
 * the vdw-perturbation. The coulomb-perturbation is controlled by directly
 * setting the charge of the atoms at the lambda states A/B.
 */
struct AtomData
{
    AtomData()
    {
        idef.atnr = c_numAtomTypes;
        idef.iparams.resize(2 * c_numAtomTypes * c_numAtomTypes);

        // set interaction parameters for different combinations of types
        idef.iparams[0].lj = { 0.001458, 1.0062882e-6 }; // 0-0
        idef.iparams[1].lj = { 0.0, 0.0 };               // 0-1
        idef.iparams[2].lj = { 0.001458, 1.0062882e-6 }; // 0-2
        idef.iparams[3].lj = { 0.0, 0.0 };               // 1-0
        idef.iparams[4].lj = { 0.0, 0.0 };               // 1-1
        idef.iparams[5].lj = { 0.0, 0.0 };               // 1-2
        idef.iparams[6].lj = { 0.001458, 1.0062882e-6 }; // 2-0
        idef.iparams[7].lj = { 0.0, 0.0 };               // 2-1
        idef.iparams[8].lj = { 0.001458, 1.0062882e-6 }; // 2-2

        GMX_ASSERT(chargeA.size() == c_numAtoms, "This test wants 4 atoms");
        GMX_ASSERT(chargeB.size() == c_numAtoms, "This test wants 4 atoms");
        GMX_ASSERT(typeA.size() == c_numAtoms, "This test wants 4 atoms");
        GMX_ASSERT(typeB.size() == c_numAtoms, "This test wants 4 atoms");
    }

    // forcefield parameters
    gmx_ffparams_t idef;

    // atom data
    std::vector<real> chargeA = { 1.0, -1.0, -1.0, 1.0 };
    std::vector<real> chargeB = { 1.0, 0.0, 0.0, 1.0 };
    std::vector<int>  typeA   = { 0, 0, 0, 0 };
    std::vector<int>  typeB   = { 0, 1, 2, 1 };
    // perturbation pattern: {no-pert, vdw- and coul-pert, coul-pert, vdw-pert}

    // neighbourhood information
    std::vector<int> iAtoms  = { 0 };
    std::vector<int> jAtoms  = { 0, 1, 2, 3 };
    std::vector<int> jIndex  = { 0, 4 };
    std::vector<int> shift   = { 0 };
    std::vector<int> gid     = { 0 };
    std::vector<int> exclFep = { 0, 1, 1, 1 };

    // construct t_nblist
    t_nblist getNbList()
    {
        t_nblist nbl;
        nbl.nri      = 1;
        nbl.nrj      = 4;
        nbl.iinr     = iAtoms;
        nbl.jindex   = jIndex;
        nbl.jjnr     = jAtoms;
        nbl.shift    = shift;
        nbl.gid      = gid;
        nbl.excl_fep = exclFep;
        return nbl;
    }
};

/*! \brief Input structure for nonbonded fep kernel
 */
struct ListInput
{
public:
    //! Function type
    int fType = F_LJ;
    //! Tolerance for float evaluation
    float floatToler = 1e-6;
    //! Tolerance for double evaluation
    double doubleToler = 1e-8;
    //! atom parameters
    AtomData atoms;
    //! forcerec helper
    ForcerecHelper frHelper;

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

    /*! \brief Set parameters for nonbonded interaction
     *
     * \param[in] coulType coulomb type
     * \param[in] vdwType  vdw type
     * \param[in] vdwMod   vdw potential modifier
     */
    ListInput setInteraction(CoulombInteractionType coulType, VanDerWaalsType vdwType, InteractionModifiers vdwMod)
    {
        frHelper.initForcerec(atoms.idef, coulType, vdwType, vdwMod);
        return *this;
    }
};

class NonbondedFepTest :
    public ::testing::TestWithParam<std::tuple<SoftcoreType, ListInput, PaddedVector<RVec>, real, real, bool>>
{
protected:
    PaddedVector<RVec>   x_;
    ListInput            input_;
    real                 lambda_;
    real                 softcoreAlpha_;
    bool                 softcoreCoulomb_;
    SoftcoreType         softcoreType_;
    TestReferenceData    refData_;
    TestReferenceChecker checker_;

    NonbondedFepTest() : checker_(refData_.rootChecker())
    {
        softcoreType_    = std::get<0>(GetParam());
        input_           = std::get<1>(GetParam());
        x_               = std::get<2>(GetParam());
        lambda_          = std::get<3>(GetParam());
        softcoreAlpha_   = std::get<4>(GetParam());
        softcoreCoulomb_ = std::get<5>(GetParam());

        // Note that the reference data for Ewald type interactions has been generated
        // with accurate analytical approximations for the long-range corrections.
        // When the free-energy kernel switches from tabulated to analytical corrections,
        // the double precision tolerance can be tightend to 1e-11.
        test::FloatingPointTolerance tolerance(
                input_.floatToler, input_.doubleToler, 1.0e-6, 1.0e-11, 10000, 100, false);
        checker_.setDefaultTolerance(tolerance);
    }

    void testKernel()
    {
        input_.frHelper.setSoftcoreAlpha(softcoreAlpha_);
        input_.frHelper.setSoftcoreCoulomb(softcoreCoulomb_);
        input_.frHelper.setSoftcoreType(softcoreType_);

        // get forcerec and interaction_const
        t_forcerec fr;
        input_.frHelper.getForcerec(&fr);

        // t_nblist
        t_nblist nbl = input_.atoms.getNbList();

        // output buffers
        OutputQuantities output;

        // lambda vector
        int numFepCouplingTerms = static_cast<int>(FreeEnergyPerturbationCouplingType::Count);
        std::vector<real> lambdas(numFepCouplingTerms, lambda_);

        // fep kernel data
        int doNBFlags = 0;
        doNBFlags |= GMX_NONBONDED_DO_FORCE;
        doNBFlags |= GMX_NONBONDED_DO_SHIFTFORCE;
        doNBFlags |= GMX_NONBONDED_DO_POTENTIAL;

        // force buffers
        bool                      unusedBool = true; // this bool has no effect in the kernel
        gmx::ForceWithShiftForces forces(output.f.arrayRefWithPadding(), unusedBool, output.fShift);

        // dummy counter
        t_nrnb nrnb;

        // run fep kernel
        gmx_nb_free_energy_kernel(nbl,
                                  x_.arrayRefWithPadding(),
                                  fr.use_simd_kernels,
                                  fr.ntype,
                                  *fr.ic,
                                  fr.shift_vec,
                                  fr.nbfp,
                                  fr.ljpme_c6grid,
                                  input_.atoms.chargeA,
                                  input_.atoms.chargeB,
                                  input_.atoms.typeA,
                                  input_.atoms.typeB,
                                  doNBFlags,
                                  lambdas,
                                  &nrnb,
                                  output.f.arrayRefWithPadding(),
                                  as_rvec_array(output.fShift.data()),
                                  output.energy.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR],
                                  output.energy.energyGroupPairTerms[NonBondedEnergyTerms::LJSR],
                                  output.dvdLambda);

        checkOutput(&checker_, output);
    }
};

TEST_P(NonbondedFepTest, testKernel)
{
    testKernel();
}

//! configurations to test
std::vector<ListInput> c_interaction = {
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::Cut, VanDerWaalsType::Cut, InteractionModifiers::None) },
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::Cut, VanDerWaalsType::Cut, InteractionModifiers::PotSwitch) },
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::Pme, VanDerWaalsType::Pme, InteractionModifiers::None) }
};

//! test parameters
std::vector<real>         c_fepLambdas                                  = { 0.0, 0.5, 1.0 };
std::vector<real>         c_softcoreBeutlerAlphaOrGapsysLinpointScaling = { 0.0, 0.3 };
std::vector<bool>         c_softcoreCoulomb                             = { true, false };
std::vector<SoftcoreType> c_softcoreType = { SoftcoreType::Beutler, SoftcoreType::Gapsys };

//! Coordinates for testing
std::vector<PaddedVector<RVec>> c_coordinates = {
    { { 1.0, 1.0, 1.0 }, { 1.1, 1.15, 1.2 }, { 0.9, 0.85, 0.8 }, { 1.1, 1.15, 0.8 } }
};

/*! \brief Coordinates for testing with near atomic overlap
 *
 * The distance between atoms 0 and 1 is set up such that the LJ-PME
 * exclusion correction, in combination with the parameters used,
 * uses the full expression in double precision and the approximation
 * in single precision, so these are both checked.
 */
std::vector<PaddedVector<RVec>> c_coordinatesShortDistance = {
    { { 1.0, 1.0, 1.0 }, { 1.014, 1.005, 1.008 }, { 0.9, 0.85, 0.8 }, { 1.1, 1.15, 0.8 } }
};

INSTANTIATE_TEST_SUITE_P(NBInteraction,
                         NonbondedFepTest,
                         ::testing::Combine(::testing::ValuesIn(c_softcoreType),
                                            ::testing::ValuesIn(c_interaction),
                                            ::testing::ValuesIn(c_coordinates),
                                            ::testing::ValuesIn(c_fepLambdas),
                                            ::testing::ValuesIn(c_softcoreBeutlerAlphaOrGapsysLinpointScaling),
                                            ::testing::ValuesIn(c_softcoreCoulomb)));

// Tests whether nearly overlapping atoms are handled correctly
INSTANTIATE_TEST_SUITE_P(NBInteractionShortDistance,
                         NonbondedFepTest,
                         ::testing::Combine(::testing::ValuesIn(c_softcoreType),
                                            ::testing::ValuesIn(c_interaction),
                                            ::testing::ValuesIn(c_coordinatesShortDistance),
                                            ::testing::Values(1.0),
                                            ::testing::Values(0.3),
                                            ::testing::Values(true)));

} // namespace

} // namespace test

} // namespace gmx

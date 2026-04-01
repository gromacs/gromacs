/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * \brief Implements test of nonbonded fep kernel on gpu
 *
 * Implements the test logic from the bonded interactions also for the
 * nonbonded fep kernel on gpu. This requires setting up some more input
 * structures that in the bonded case.
 *
 * The test setup consists of an atom pair that is evaluated in an fep setting
 * (vanishing charge and lennard-jones parameters of atom #2) with and without
 * softcore Potentials.
 *
 * \author Yiqi Chen <yiqi.chen@metax-tech.com>
 * \ingroup module_gmxlib_nonbonded
 */
#include "gmxpre.h"

#include "config.h"

#include <cmath>

#include <array>
#include <iterator>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

// Currently FEP-on-GPU calculations are only implemented on CUDA platform
#if GMX_GPU_CUDA
#    include "gromacs/ewald/ewald_utils.h"
#    include "gromacs/gpu_utils/device_stream_manager.h"
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/gpu_utils.h"
#    include "gromacs/gpu_utils/gputraits.h"
#    include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#    include "gromacs/hardware/device_information.h"
#    include "gromacs/hardware/device_management.h"
#    include "gromacs/math/arrayrefwithpadding.h"
#    include "gromacs/math/paddedvector.h"
#    include "gromacs/math/units.h"
#    include "gromacs/mdlib/forcerec.h"
#    include "gromacs/mdtypes/enerdata.h"
#    include "gromacs/mdtypes/forceoutput.h"
#    include "gromacs/mdtypes/forcerec.h"
#    include "gromacs/mdtypes/interaction_const.h"
#    include "gromacs/mdtypes/md_enums.h"
#    include "gromacs/mdtypes/mdatom.h"
#    include "gromacs/mdtypes/simulation_workload.h"
#    include "gromacs/nbnxm/atomdata.h"
#    include "gromacs/nbnxm/atompairlist.h"
#    if GMX_GPU_CUDA
#        include "gromacs/nbnxm/cuda/nbnxm_cuda_types.h"
#    endif
#    include "gromacs/nbnxm/gpu_data_mgmt.h"
#    include "gromacs/nbnxm/gpu_types_common.h"
#    include "gromacs/nbnxm/nbnxm.h"
#    include "gromacs/nbnxm/nbnxm_gpu.h"
#    include "gromacs/nbnxm/pairlist.h"
#    include "gromacs/nbnxm/pairlistparams.h"
#    include "gromacs/pbcutil/ishift.h"
#    include "gromacs/pbcutil/pbc.h"
#    include "gromacs/tables/forcetable.h"
#    include "gromacs/topology/forcefieldparameters.h"
#    include "gromacs/topology/idef.h"
#    include "gromacs/topology/ifunc.h"
#    include "gromacs/utility/arrayref.h"
#    include "gromacs/utility/enumerationhelpers.h"
#    include "gromacs/utility/gmxassert.h"
#    include "gromacs/utility/logger.h"
#    include "gromacs/utility/real.h"
#    include "gromacs/utility/stringutil.h"
#    include "gromacs/utility/template_mp.h"
#    include "gromacs/utility/vec.h"
#    include "gromacs/utility/vectypes.h"

#    include "testutils/refdata.h"
#    include "testutils/test_hardware_environment.h"
#    include "testutils/testasserts.h"

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
        tmp.coulomb.ewaldCoeff = calc_ewaldcoeff_q(1.0, 1.0e-5);
        coulEwaldCoeff_        = tmp.coulomb.ewaldCoeff;
        tmp.vdw.ewaldCoeff     = calc_ewaldcoeff_lj(1.0, 1.0e-5);
        vdwEwaldCoeff_         = tmp.vdw.ewaldCoeff;
        tmp.coulomb.type       = coulType;
        tmp.vdw.type           = vdwType;
        tmp.vdw.modifier       = vdwMod;
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
        ic->softCoreParameters = std::make_unique<interaction_const_t::SoftCoreParameters>(fepVals);

        ic->coulombEwaldTables  = std::make_unique<EwaldCorrectionTables>();
        *ic->coulombEwaldTables = coulombTables_;

        ic->vdwEwaldTables  = std::make_unique<EwaldCorrectionTables>();
        *ic->vdwEwaldTables = vdwTables_;

        // set coulomb and vdw types
        ic->coulomb.type = coulType_;
        ic->vdw.type     = vdwType_;
        ic->vdw.modifier = vdwMod_;

        // some non default parameters used in this testcase
        ic->coulomb.epsfac                   = gmx::c_one4PiEps0 * 0.25;
        ic->coulomb.reactionFieldCoefficient = 0.0; // former k_rf
        ic->coulomb.reactionFieldShift       = 1.0; // former c_rf
        ic->coulomb.ewaldCoeff               = coulEwaldCoeff_;
        ic->vdw.ewaldCoeff                   = vdwEwaldCoeff_;
        ic->coulomb.ewaldShift               = 1.0e-5;
        ic->vdw.ewaldShift                   = -1.0;
        ic->vdw.dispersionShift.cpot         = -1.0;
        ic->vdw.repulsionShift.cpot          = -1.0;
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
        nbfp_ = makeNonBondedParameterLists(idef.atnr, false, idef.iparams, false);
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
    std::vector<int>  iAtoms  = { 0 };
    std::vector<int>  jAtoms  = { 0, 1, 2, 3 };
    std::vector<int>  shift   = { 0 };
    std::vector<int>  gid     = { 0 };
    std::vector<bool> exclFep = { false, true, true, true };

    // construct AtomPairlist
    AtomPairlist getNbList()
    {
        AtomPairlist nbl;
        nbl.addIEntry({ iAtoms[0], shift[0], gid[0] }, jAtoms.size());
        for (size_t j = 0; j < jAtoms.size(); j++)
        {
            nbl.addJEntry({ jAtoms[j], exclFep[j] });
        }
        return nbl;
    }
};

/*! \brief Input structure for nonbonded fep kernel
 */
struct ListInput
{
public:
    //! Function type
    InteractionFunction fType = InteractionFunction::LennardJonesShortRange;
    //! Tolerance for float evaluation
    float floatToler = 1e-6;
    //! Tolerance for double evaluation
    double doubleToler = 1e-8;
    //! atom parameters
    AtomData atoms;
    //! forcerec helper
    ForcerecHelper frHelper;
    //! LJ combination rule (for Cut VdW; determines Cut vs CutCombGeom vs CutCombLB)
    LJCombinationRule ljCombinationRule = LJCombinationRule::None;
    //! Description for test naming
    std::string description;

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
     * \param[in] ljComb   LJ combination rule (for VanDerWaalsType::Cut with None/PotShift)
     */
    ListInput setInteraction(CoulombInteractionType coulType,
                             VanDerWaalsType        vdwType,
                             InteractionModifiers   vdwMod,
                             LJCombinationRule      ljComb = LJCombinationRule::None)
    {
        ljCombinationRule = ljComb;
        frHelper.initForcerec(atoms.idef, coulType, vdwType, vdwMod);
        description = formatString("coul_%s_vdw_%s_vdwmod_%s_ljrule_%s",
                                   enumValueToString(coulType),
                                   enumValueToString(vdwType),
                                   enumValueToString(vdwMod),
                                   enumValueToString(ljComb));
        return *this;
    }
};

//! Build nbnxn_atomdata_t for the 4-atom FEP test using production code path
static std::unique_ptr<nbnxn_atomdata_t> makeNbatForFepTest(const AtomData&           atoms,
                                                            const PaddedVector<RVec>& x,
                                                            const t_forcerec&         fr,
                                                            ArrayRef<const real>      lambdas,
                                                            LJCombinationRule ljCombinationRule)
{
    const int  numAtoms   = c_numAtoms;
    const real lambdaCoul = lambdas[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)];
    const real lambdaVdw  = lambdas[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)];

    std::optional<LJCombinationRule> ljRule = (ljCombinationRule != LJCombinationRule::None)
                                                      ? std::optional<LJCombinationRule>(ljCombinationRule)
                                                      : std::nullopt;

    auto nbat = std::make_unique<nbnxn_atomdata_t>(PinningPolicy::PinnedIfSupported,
                                                   MDLogger(),
                                                   NbnxmKernelType::Gpu8x8x8,
                                                   ljRule,
                                                   LJCombinationRule::None,
                                                   fr.nbfp,
                                                   true,
                                                   1,
                                                   1);

    nbnxn_atomdata_t::Params& params = nbat->paramsDeprecated();
    params.typeA.resize(numAtoms);
    params.typeB.resize(numAtoms);
    params.qA.resize(numAtoms);
    params.qB.resize(numAtoms);
    params.type.resize(numAtoms);
    for (int i = 0; i < numAtoms; i++)
    {
        params.typeA[i] = atoms.typeA[i];
        params.typeB[i] = atoms.typeB[i];
        params.qA[i]    = atoms.chargeA[i];
        params.qB[i]    = atoms.chargeB[i];
        params.type[i]  = atoms.typeA[i];
    }

    if (params.ljCombinationRule != LJCombinationRule::None)
    {
        params.lj_comb.resize(numAtoms * 2);
        params.ljCombA.resize(numAtoms * 2);
        params.ljCombB.resize(numAtoms * 2);
        for (int i = 0; i < numAtoms; i++)
        {
            const int tA = params.typeA[i], tB = params.typeB[i];
            params.lj_comb[2 * i]     = params.nbfp_comb[tA * 2];
            params.lj_comb[2 * i + 1] = params.nbfp_comb[tA * 2 + 1];
            params.ljCombA[2 * i]     = params.nbfp_comb[tA * 2];
            params.ljCombA[2 * i + 1] = params.nbfp_comb[tA * 2 + 1];
            params.ljCombB[2 * i]     = params.nbfp_comb[tB * 2];
            params.ljCombB[2 * i + 1] = params.nbfp_comb[tB * 2 + 1];
        }
    }

    nbat->resizeCoordinateBuffer(numAtoms, 0);
    std::copy(fr.shift_vec.begin(), fr.shift_vec.end(), nbat->shift_vec.begin());
    nbat->bDynamicBox = false;

    // Fill x with x,y,z,q (nbatXYZQ format); q = interpolated charge
    ArrayRef<real> xbat = nbat->x();
    for (int i = 0; i < numAtoms; i++)
    {
        const real q    = (1 - lambdaCoul) * atoms.chargeA[i] + lambdaCoul * atoms.chargeB[i];
        xbat[4 * i]     = x[i][XX];
        xbat[4 * i + 1] = x[i][YY];
        xbat[4 * i + 2] = x[i][ZZ];
        xbat[4 * i + 3] = q;
    }

    return nbat;
}

//! Identity atom index mapping for the 4-atom test (grid-local index == global index)
static const std::array<int, c_numAtoms> c_identityAtomIndices = { 0, 1, 2, 3 };

// Copy results back to CPU
void gpuLaunchCpyback(NbnxmGpu* nb, struct OutputQuantities* output)
{
    /* extract the data */
    NBAtomDataGpu*      atdat        = nb->atdat;
    const DeviceStream& deviceStream = *nb->deviceStreams[0];

    // copy f back to ouput.f
    copyFromDeviceBuffer(
            output->f.data(), &atdat->f, 0, c_numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    // copy fShift back to ouput.fShift
    copyFromDeviceBuffer(
            output->fShift.data(), &atdat->fShift, 0, 1, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    copyFromDeviceBuffer(output->energy.energyGroupPairTerms[NonBondedEnergyTerms::LJSR].data(),
                         &atdat->eLJ,
                         0,
                         1,
                         deviceStream,
                         GpuApiCallBehavior::Sync,
                         nullptr);

    copyFromDeviceBuffer(output->energy.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR].data(),
                         &atdat->eElec,
                         0,
                         1,
                         deviceStream,
                         GpuApiCallBehavior::Sync,
                         nullptr);

    copyFromDeviceBuffer(&output->dvdLambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)],
                         &atdat->dvdlLJ,
                         0,
                         1,
                         deviceStream,
                         GpuApiCallBehavior::Sync,
                         nullptr);

    copyFromDeviceBuffer(&output->dvdLambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)],
                         &atdat->dvdlElec,
                         0,
                         1,
                         deviceStream,
                         GpuApiCallBehavior::Sync,
                         nullptr);
}

struct CoordinateSet
{
    std::string        name;
    PaddedVector<RVec> x;
};

using NonbondedFepGpuTestParameters = std::tuple<SoftcoreType, ListInput, CoordinateSet, real, real, bool>;

/*! \brief Help GoogleTest name our test cases */
std::string nameOfTest(const testing::TestParamInfo<NonbondedFepGpuTestParameters>& info)
{
    std::string testName = formatString("softcore_%s_list_%s_coords_%s_%g_%g_scCoulomb_%s",
                                        enumValueToString(std::get<0>(info.param)),
                                        std::get<1>(info.param).description.c_str(),
                                        std::get<2>(info.param).name.c_str(),
                                        std::get<3>(info.param),
                                        std::get<4>(info.param),
                                        std::get<5>(info.param) ? "Yes" : "No");

    // Note that the returned names must be unique and may use only
    // alphanumeric ASCII characters. It's not supposed to contain
    // underscores (see the GoogleTest FAQ
    // why-should-test-suite-names-and-test-names-not-contain-underscore),
    // but doing so works for now, is likely to remain so, and makes
    // such test names much more readable.
    testName = replaceAll(testName, "-", "_");
    testName = replaceAll(testName, ".", "_");
    testName = replaceAll(testName, " ", "_");
    return testName;
}

class NonbondedFepGpuTest : public ::testing::TestWithParam<NonbondedFepGpuTestParameters>
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

    NonbondedFepGpuTest() : checker_(refData_.rootChecker())
    {
        softcoreType_    = std::get<0>(GetParam());
        input_           = std::get<1>(GetParam());
        x_               = std::get<2>(GetParam()).x;
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
        // Skip the check for unused entries to pass this test when running GPU builds on CPU-only machines.
        checker_.disableUnusedEntriesCheck();
    }

    void testGpuKernel()
    {
        input_.frHelper.setSoftcoreAlpha(softcoreAlpha_);
        input_.frHelper.setSoftcoreCoulomb(softcoreCoulomb_);
        input_.frHelper.setSoftcoreType(softcoreType_);

        // get forcerec and interaction_const
        t_forcerec fr{ false };
        input_.frHelper.getForcerec(&fr);

        // AtomPairlist
        AtomPairlist nbl = input_.atoms.getNbList();

        // output buffers
        OutputQuantities output;

        // lambda vector
        int numFepCouplingTerms = static_cast<int>(FreeEnergyPerturbationCouplingType::Count);
        std::vector<real> lambdas(numFepCouplingTerms, lambda_);

        // fep kernel data
        StepWorkload stepWork;
        stepWork.computeForces = true;
        stepWork.computeVirial = true;
        stepWork.computeEnergy = true;

        // dummy counter
        // t_nrnb nrnb;
        for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
        {
            testDevice->deviceContext().activate();
            SimulationWorkload simulationWorkload;
            auto               deviceStreamManager = std::make_shared<DeviceStreamManager>(
                    testDevice->deviceInfo(), simulationWorkload, false);

            auto nbat = makeNbatForFepTest(input_.atoms, x_, fr, lambdas, input_.ljCombinationRule);

            PairlistParams pairlistParams(NbnxmKernelType::Gpu8x8x8, sc_layoutType, true, fr.rlist, false);
            pairlistParams.haveNonbondedFEGpu_ = true;

            NbnxmGpu* nb = gpu_init(
                    *deviceStreamManager, fr.ic.get(), pairlistParams, nbat.get(), false, std::optional<int>(1));

            EnumerationArray<FreeEnergyPerturbationCouplingType, std::vector<double>> all_lambda;
            all_lambda[FreeEnergyPerturbationCouplingType::Coul] = {
                lambdas[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)]
            };
            all_lambda[FreeEnergyPerturbationCouplingType::Vdw] = {
                lambdas[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)]
            };
            copy_gpu_fepparams(nb,
                               true,
                               fr.ic->softCoreParameters->alphaCoulomb,
                               fr.ic->softCoreParameters->alphaVdw,
                               fr.ic->softCoreParameters->lambdaPower,
                               fr.ic->softCoreParameters->sigma6WithInvalidSigma,
                               fr.ic->softCoreParameters->sigma6Minimum,
                               lambdas[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)],
                               lambdas[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)],
                               1,
                               all_lambda);

            gpu_init_atomdata(nb, nbat.get());
            gpu_upload_shiftvec(nb, nbat.get());
            gpu_copy_xq_to_gpu(nb, nbat.get(), AtomLocality::Local);
            gpu_init_feppairlist(nb, nbl, InteractionLocality::Local, c_identityAtomIndices);

            // Clear force buffer before kernel (matches production step flow)
            gpu_clear_outputs(nb, stepWork.computeVirial);

            gpu_launch_free_energy_kernel(nb, simulationWorkload, stepWork, InteractionLocality::Local);

            gpuLaunchCpyback(nb, &output);

            gpu_free(nb);

            checkOutput(&checker_, output);
        }
    }
};

TEST_P(NonbondedFepGpuTest, testGpuKernel)
{
    testGpuKernel();
}

//! configurations to test (each (elec, vdw, modifier) exercises a different GPU kernel flavor)
std::vector<ListInput> c_interaction = {
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::Cut, VanDerWaalsType::Cut, InteractionModifiers::None) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::Cut,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::None,
                              LJCombinationRule::Geometric) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::Cut,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::None,
                              LJCombinationRule::LorentzBerthelot) },
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::Cut, VanDerWaalsType::Cut, InteractionModifiers::PotShift) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::Cut,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::PotShift,
                              LJCombinationRule::Geometric) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::Cut,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::PotShift,
                              LJCombinationRule::LorentzBerthelot) },
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::Cut, VanDerWaalsType::Cut, InteractionModifiers::ForceSwitch) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::Cut,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::ForceSwitch,
                              LJCombinationRule::Geometric) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::Cut,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::ForceSwitch,
                              LJCombinationRule::LorentzBerthelot) },
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::RF, VanDerWaalsType::Cut, InteractionModifiers::None) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::RF,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::None,
                              LJCombinationRule::Geometric) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::RF,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::None,
                              LJCombinationRule::LorentzBerthelot) },
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::RF, VanDerWaalsType::Cut, InteractionModifiers::PotShift) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::RF,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::PotShift,
                              LJCombinationRule::Geometric) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::RF,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::PotShift,
                              LJCombinationRule::LorentzBerthelot) },
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::RF, VanDerWaalsType::Cut, InteractionModifiers::ForceSwitch) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::RF,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::ForceSwitch,
                              LJCombinationRule::Geometric) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::RF,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::ForceSwitch,
                              LJCombinationRule::LorentzBerthelot) },
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::Pme, VanDerWaalsType::Cut, InteractionModifiers::None) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::Pme,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::None,
                              LJCombinationRule::Geometric) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::Pme,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::None,
                              LJCombinationRule::LorentzBerthelot) },
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::Pme, VanDerWaalsType::Cut, InteractionModifiers::PotShift) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::Pme,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::PotShift,
                              LJCombinationRule::Geometric) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::Pme,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::PotShift,
                              LJCombinationRule::LorentzBerthelot) },
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::Pme, VanDerWaalsType::Cut, InteractionModifiers::ForceSwitch) },
    { ListInput(1e-6, 1e-8)
              .setInteraction(CoulombInteractionType::Pme,
                              VanDerWaalsType::Cut,
                              InteractionModifiers::ForceSwitch,
                              LJCombinationRule::Geometric) }
};

//! test parameters
std::vector<real> c_fepLambdas                                  = { 0.0, 0.5, 1.0 };
std::vector<real> c_softcoreBeutlerAlphaOrGapsysLinpointScaling = { 0.0, 0.3 };
std::vector<bool> c_softcoreCoulomb                             = { true, false };
//  No Gapsys Softcore function support on GPU yet
std::vector<SoftcoreType> c_softcoreType = { SoftcoreType::Beutler };

//! Named sets of coordinates for testing
std::vector<CoordinateSet> c_coordinateSets = {
    { "A", { { 1.0, 1.0, 1.0 }, { 1.1, 1.15, 1.2 }, { 0.9, 0.85, 0.8 }, { 1.1, 1.15, 0.8 } } },
};

INSTANTIATE_TEST_SUITE_P(NBInteraction,
                         NonbondedFepGpuTest,
                         ::testing::Combine(::testing::ValuesIn(c_softcoreType),
                                            ::testing::ValuesIn(c_interaction),
                                            ::testing::ValuesIn(c_coordinateSets),
                                            ::testing::ValuesIn(c_fepLambdas),
                                            ::testing::ValuesIn(c_softcoreBeutlerAlphaOrGapsysLinpointScaling),
                                            ::testing::ValuesIn(c_softcoreCoulomb)),
                         nameOfTest);


} // namespace

} // namespace test

} // namespace gmx
#endif

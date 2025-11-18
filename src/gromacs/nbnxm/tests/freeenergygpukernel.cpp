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

#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/atompairlist.h"
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
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

#if GMX_GPU_CUDA
#    include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#    include "gromacs/gpu_utils/cuda_kernel_utils.cuh"
#    include "gromacs/gpu_utils/cudautils.cuh"
#    include "gromacs/gpu_utils/device_stream_manager.h"
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/gpu_utils.h"
#    include "gromacs/gpu_utils/gputraits.h"
#    include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#    include "gromacs/gpu_utils/vectype_ops_cuda.h"
#    include "gromacs/hardware/device_information.h"
#    include "gromacs/hardware/device_management.h"
#    include "gromacs/mdtypes/simulation_workload.h"
#    include "gromacs/nbnxm/cuda/nbnxm_cuda.h"
#    include "gromacs/nbnxm/cuda/nbnxm_cuda_kernel_utils.cuh"
#    include "gromacs/nbnxm/cuda/nbnxm_cuda_types.h"
#    include "gromacs/nbnxm/gpu_types_common.h"
#    include "gromacs/nbnxm/nbnxm.h"
#    include "gromacs/nbnxm/nbnxm_gpu.h"

/** Force & energy **/
#    define CALC_ENERGIES
#    include "gromacs/nbnxm/cuda/nbfe_cuda_kernels.cuh"
#    undef CALC_ENERGIES

#endif

#include "testutils/refdata.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{
// Currently FEP-on-GPU calculations are only implemented on CUDA platform
#if GMX_GPU_CUDA
/*! Nonbonded FEP kernel function pointer type */
typedef void (*nbfe_cu_kfunc_ptr_t)(const NBAtomDataGpu, const NBParamGpu, const GpuFeplist, bool);

/*! Force + energy fep kernel function pointers. */
static const nbfe_cu_kfunc_ptr_t nbfe_kfunc_ener_ptr[c_numElecTypes][c_numVdwTypes] = {
    { nbfe_kernel_ElecCut_VdwLJ_VF_cuda,
      nbfe_kernel_ElecCut_VdwLJCombGeom_VF_cuda,
      nbfe_kernel_ElecCut_VdwLJCombLB_VF_cuda,
      nbfe_kernel_ElecCut_VdwLJFsw_VF_cuda,
      nbfe_kernel_ElecCut_VdwLJPsw_VF_cuda,
      nbfe_kernel_ElecCut_VdwLJEwCombGeom_VF_cuda,
      nbfe_kernel_ElecCut_VdwLJEwCombLB_VF_cuda },
    { nbfe_kernel_ElecRF_VdwLJ_VF_cuda,
      nbfe_kernel_ElecRF_VdwLJCombGeom_VF_cuda,
      nbfe_kernel_ElecRF_VdwLJCombLB_VF_cuda,
      nbfe_kernel_ElecRF_VdwLJFsw_VF_cuda,
      nbfe_kernel_ElecRF_VdwLJPsw_VF_cuda,
      nbfe_kernel_ElecRF_VdwLJEwCombGeom_VF_cuda,
      nbfe_kernel_ElecRF_VdwLJEwCombLB_VF_cuda },
    { nbfe_kernel_ElecEwQSTab_VdwLJ_VF_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJCombGeom_VF_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJCombLB_VF_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJFsw_VF_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJPsw_VF_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJEwCombGeom_VF_cuda,
      nbfe_kernel_ElecEwQSTab_VdwLJEwCombLB_VF_cuda },
    { nbfe_kernel_ElecEwQSTabTwinCut_VdwLJ_VF_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_VF_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_VF_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJFsw_VF_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJPsw_VF_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_VF_cuda,
      nbfe_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_VF_cuda },
    { nbfe_kernel_ElecEw_VdwLJ_VF_cuda,
      nbfe_kernel_ElecEw_VdwLJCombGeom_VF_cuda,
      nbfe_kernel_ElecEw_VdwLJCombLB_VF_cuda,
      nbfe_kernel_ElecEw_VdwLJFsw_VF_cuda,
      nbfe_kernel_ElecEw_VdwLJPsw_VF_cuda,
      nbfe_kernel_ElecEw_VdwLJEwCombGeom_VF_cuda,
      nbfe_kernel_ElecEw_VdwLJEwCombLB_VF_cuda },
    { nbfe_kernel_ElecEwTwinCut_VdwLJ_VF_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJCombGeom_VF_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJCombLB_VF_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJFsw_VF_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJPsw_VF_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJEwCombGeom_VF_cuda,
      nbfe_kernel_ElecEwTwinCut_VdwLJEwCombLB_VF_cuda }
};

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

// Function to calculate potential switch constants
static void potential_switch_constants(real rsw, real rc, switch_consts_t* sc)
{
    /* The switch function is 1 at rsw and 0 at rc.
     * The derivative and second derivate are zero at both ends.
     * rsw        = max(r - r_switch, 0)
     * sw         = 1 + c3*rsw^3 + c4*rsw^4 + c5*rsw^5
     * dsw        = 3*c3*rsw^2 + 4*c4*rsw^3 + 5*c5*rsw^4
     * force      = force*dsw - potential*sw
     * potential *= sw
     */
    sc->c3 = -10 / gmx::power3(rc - rsw);
    sc->c4 = 15 / gmx::power4(rc - rsw);
    sc->c5 = -6 / gmx::power5(rc - rsw);
}

// Brief function to set up FEP GPU parameters
void setGpuNbParams(NBParamGpu*                nbp,
                    const t_forcerec&          fr,
                    const interaction_const_t& ic,
                    const std::vector<real>&   lambdas,
                    const DeviceContext&       deviceContext,
                    const DeviceStream&        localStream)
{

    nbp->rvdw_sq     = fr.rlist * fr.rlist;
    nbp->rcoulomb_sq = nbp->rvdw_sq;

    nbp->epsfac        = ic.coulomb.epsfac;
    nbp->c_rf          = ic.coulomb.reactionFieldShift;
    nbp->two_k_rf      = 2.0 * ic.coulomb.reactionFieldCoefficient;
    nbp->sh_ewald      = ic.coulomb.ewaldShift;
    nbp->sh_lj_ewald   = ic.vdw.ewaldShift;
    nbp->ewaldcoeff_lj = ic.vdw.ewaldCoeff;
    nbp->ewald_beta    = ic.coulomb.ewaldCoeff;


    VdwType vdwType = ic.vdw.modifier == InteractionModifiers::PotSwitch ? VdwType::PSwitch : VdwType::Cut;
    nbp->vdwType  = vdwType;
    nbp->elecType = static_cast<ElecType>(ic.coulomb.type);

    nbp->dispersion_shift = ic.vdw.dispersionShift;
    nbp->repulsion_shift  = ic.vdw.repulsionShift;
    nbp->rvdw_switch      = ic.vdw.switchDistance;
    nbp->vdw_switch       = ic.vdw.switchConstants;

    if (nbp->vdwType == VdwType::PSwitch)
    {
        potential_switch_constants(ic.vdw.switchDistance, ic.vdw.cutoff, &(nbp->vdw_switch));
    }

    static_assert(sizeof(decltype(nbp->nbfp)) == 2 * sizeof(decltype(*fr.nbfp.data())),
                  "Mismatch in the size of host / device data types");

    allocateDeviceBuffer(&nbp->nbfp, c_numAtomTypes * c_numAtomTypes, deviceContext);
    copyToDeviceBuffer(&nbp->nbfp,
                       reinterpret_cast<const Float2*>(fr.nbfp.data()),
                       0,
                       c_numAtomTypes * c_numAtomTypes,
                       localStream,
                       GpuApiCallBehavior::Sync,
                       nullptr);

    if (!c_disableCudaTextures)
    {
        cudaResourceDesc rd;
        cudaTextureDesc  td;

        memset(&rd, 0, sizeof(rd));
        rd.resType                = cudaResourceTypeLinear;
        rd.res.linear.devPtr      = nbp->nbfp;
        rd.res.linear.desc        = cudaCreateChannelDesc<Float2>();
        rd.res.linear.sizeInBytes = c_numAtomTypes * sizeof(Float2);

        memset(&td, 0, sizeof(td));
        td.readMode      = cudaReadModeElementType;
        cudaError_t stat = cudaCreateTextureObject(&nbp->nbfp_texobj, &rd, &td, nullptr);
        gmx::checkDeviceError(stat, "Binding of the texture object failed.");
    }

    // FEP params
    nbp->bFepGpuNonBonded       = true;
    nbp->lambdaPower            = ic.softCoreParameters->lambdaPower;
    nbp->sigma6WithInvalidSigma = ic.softCoreParameters->sigma6WithInvalidSigma;
    nbp->sigma6Minimum          = ic.softCoreParameters->sigma6Minimum;

    nbp->alphaCoul  = ic.softCoreParameters->alphaCoulomb;
    nbp->alphaVdw   = ic.softCoreParameters->alphaVdw;
    nbp->lambdaCoul = lambdas[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)];
    nbp->lambdaVdw  = lambdas[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)];
}

// Brief function to copy fep pair list to GPU
void constructFeppairlist(GpuFeplist*          d_feplist,
                          const AtomPairlist*  h_feplist,
                          const DeviceContext& deviceContext,
                          const DeviceStream&  localStream)
{

    // Extract pair list data
    ArrayRef<const AtomPairlist::IEntry> iList     = h_feplist->iList();
    ArrayRef<const AtomPairlist::JEntry> jList     = h_feplist->flatJList();
    const int                            numiAtoms = iList.ssize();
    const int                            numjAtoms = jList.ssize();

    HostVector<int> iinr(numiAtoms, { PinningPolicy::PinnedIfSupported });
    for (int i = 0; i < numiAtoms; i++)
    {
        iinr[i] = iList[i].atom;
    }

    HostVector<int> jjnr(numjAtoms, { PinningPolicy::PinnedIfSupported });
    for (int j = 0; j < numjAtoms; j++)
    {
        jjnr[j] = jList[j].atom;
    }

    HostVector<int> jIndex(numiAtoms + 1, { PinningPolicy::PinnedIfSupported });
    int             count = 0;
    for (int i = 0; i < numiAtoms; i++)
    {
        int incr = h_feplist->jList(i).size();
        count += incr;
        jIndex[i + 1] = count;
    }

    HostVector<int> shift(numiAtoms, { PinningPolicy::PinnedIfSupported });
    for (int i = 0; i < numiAtoms; i++)
    {
        shift[i] = iList[i].shiftIndex;
    }

    HostVector<int> exclFep(numjAtoms, { PinningPolicy::PinnedIfSupported });
    for (int j = 0; j < numjAtoms; j++)
    {
        exclFep[j] = jList[j].interacts;
    }

    reallocateDeviceBuffer(
            &d_feplist->iinr, numiAtoms, &d_feplist->numiAtoms, &d_feplist->maxNumiAtoms, deviceContext);
    copyToDeviceBuffer(
            &d_feplist->iinr, iinr.data(), 0, numiAtoms, localStream, GpuApiCallBehavior::Async, nullptr);

    reallocateDeviceBuffer(
            &d_feplist->shift, numiAtoms, &d_feplist->numShift, &d_feplist->maxNumShift, deviceContext);
    copyToDeviceBuffer(
            &d_feplist->shift, shift.data(), 0, numiAtoms, localStream, GpuApiCallBehavior::Async, nullptr);

    reallocateDeviceBuffer(
            &d_feplist->jIndex, numiAtoms + 1, &d_feplist->numjIndex, &d_feplist->maxNumjIndex, deviceContext);
    copyToDeviceBuffer(
            &d_feplist->jIndex, jIndex.data(), 0, numiAtoms + 1, localStream, GpuApiCallBehavior::Async, nullptr);

    reallocateDeviceBuffer(
            &d_feplist->jjnr, numjAtoms, &d_feplist->numjAtoms, &d_feplist->maxNumjAtoms, deviceContext);
    copyToDeviceBuffer(
            &d_feplist->jjnr, jjnr.data(), 0, numjAtoms, localStream, GpuApiCallBehavior::Async, nullptr);

    reallocateDeviceBuffer(
            &d_feplist->exclFep, numjAtoms, &d_feplist->numExcl, &d_feplist->maxNumExcl, deviceContext);
    copyToDeviceBuffer(
            &d_feplist->exclFep, exclFep.data(), 0, numjAtoms, localStream, GpuApiCallBehavior::Async, nullptr);
}

// A simpler version of gpu_init where only fep relevant things are initialized
NbnxmGpu* gpu_fep_init(const DeviceStreamManager& deviceStreamManager,
                       const AtomData             atoms,
                       const AtomPairlist         nblist,
                       const PaddedVector<RVec>&  x,
                       const t_forcerec&          fr,
                       const std::vector<real>&   lambdas)
{
    auto* nb           = new NbnxmGpu();
    nb->deviceContext_ = &deviceStreamManager.context();
    nb->atdat          = new NBAtomDataGpu;
    nb->nbparam        = new NBParamGpu;

    nb->feplist[0] = std::make_unique<GpuFeplist>();

    nb->bUseTwoStreams = false;

    const DeviceStream& localStream = deviceStreamManager.stream(DeviceStreamType::NonBondedLocal);
    nb->deviceStreams[InteractionLocality::Local] = &localStream;

    // initialize atdat and nbparam
    nb->atdat->numTypes = c_numAtomTypes;
    nb->atdat->numAtoms = c_numAtoms;

    allocateDeviceBuffer(&nb->atdat->shiftVec, 1, *nb->deviceContext_);
    copyToDeviceBuffer(&nb->atdat->shiftVec,
                       asGenericFloat3Pointer(fr.shift_vec),
                       0,
                       1,
                       localStream,
                       GpuApiCallBehavior::Sync,
                       nullptr);
    nb->atdat->shiftVecUploaded = true;

    allocateDeviceBuffer(&nb->atdat->fShift, 1, *nb->deviceContext_);
    allocateDeviceBuffer(&nb->atdat->eLJ, 1, *nb->deviceContext_);
    allocateDeviceBuffer(&nb->atdat->eElec, 1, *nb->deviceContext_);
    allocateDeviceBuffer(&nb->atdat->dvdlLJ, 1, *nb->deviceContext_);
    allocateDeviceBuffer(&nb->atdat->dvdlElec, 1, *nb->deviceContext_);

    clearDeviceBufferAsync(&nb->atdat->fShift, 0, 1, localStream);
    clearDeviceBufferAsync(&nb->atdat->eElec, 0, 1, localStream);
    clearDeviceBufferAsync(&nb->atdat->eLJ, 0, 1, localStream);
    clearDeviceBufferAsync(&nb->atdat->dvdlElec, 0, 1, localStream);
    clearDeviceBufferAsync(&nb->atdat->dvdlLJ, 0, 1, localStream);

    HostVector<RVec> f_h(c_numAtoms, { 0.0, 0.0, 0.0 }, { PinningPolicy::PinnedIfSupported });
    allocateDeviceBuffer(&nb->atdat->f, c_numAtoms, *nb->deviceContext_);
    copyToDeviceBuffer(
            &nb->atdat->f, f_h.data(), 0, c_numAtoms, localStream, GpuApiCallBehavior::Async, nullptr);

    allocateDeviceBuffer(&nb->atdat->xq, c_numAtoms, *nb->deviceContext_);
    HostVector<float> xq_h(c_numAtoms * 4, { PinningPolicy::PinnedIfSupported });
    for (int k = 0; k < c_numAtoms; k++)
    {
        xq_h[4 * k]     = x[k][0];
        xq_h[4 * k + 1] = x[k][1];
        xq_h[4 * k + 2] = x[k][2];
    }
    static_assert(sizeof(nb->atdat->xq[0]) == sizeof(Float4),
                  "Size of the xq parameters element should be equal to the size of float4.");
    copyToDeviceBuffer(&nb->atdat->xq,
                       reinterpret_cast<Float4*>(xq_h.data()),
                       0,
                       c_numAtoms,
                       localStream,
                       GpuApiCallBehavior::Async,
                       nullptr);

    allocateDeviceBuffer(&nb->atdat->q4, c_numAtoms, *nb->deviceContext_);
    HostVector<float> q4_h(c_numAtoms * 4, { PinningPolicy::PinnedIfSupported });
    for (int k = 0; k < c_numAtoms; k++)
    {
        q4_h[4 * k]     = (float)atoms.chargeA[k];
        q4_h[4 * k + 1] = (float)atoms.chargeB[k];
    }
    static_assert(sizeof(nb->atdat->q4[0]) == sizeof(Float4),
                  "Size of the q4 parameters element should be equal to the size of float4.");
    copyToDeviceBuffer(&nb->atdat->q4,
                       reinterpret_cast<Float4*>(q4_h.data()),
                       0,
                       c_numAtoms,
                       localStream,
                       GpuApiCallBehavior::Async,
                       nullptr);

    allocateDeviceBuffer(&nb->atdat->atomTypes4, c_numAtoms, *nb->deviceContext_);
    HostVector<int> atomTypes4_h(c_numAtoms * 4, { PinningPolicy::PinnedIfSupported });
    for (int k = 0; k < c_numAtoms; k++)
    {
        atomTypes4_h[4 * k]     = (int)atoms.typeA[k];
        atomTypes4_h[4 * k + 1] = (int)atoms.typeB[k];
    }
    static_assert(sizeof(nb->atdat->atomTypes4[0]) == sizeof(Int4),
                  "Size of the atomTypes4 parameters element should be equal to the size of Int4.");
    copyToDeviceBuffer(&nb->atdat->atomTypes4,
                       reinterpret_cast<Int4*>(atomTypes4_h.data()),
                       0,
                       c_numAtoms,
                       localStream,
                       GpuApiCallBehavior::Async,
                       nullptr);

    // set up nb fep parameters
    setGpuNbParams(nb->nbparam, fr, *fr.ic, lambdas, *nb->deviceContext_, localStream);
    // construct the gpu fep list
    constructFeppairlist(nb->feplist[0].get(), &nblist, *nb->deviceContext_, localStream);

    return nb;
}

// Brief function to launch FEP kernel on GPU
void gpuLaunchFepKernel(NbnxmGpu* nb, StepWorkload& stepWork)
{
    NBAtomDataGpu*      atdat        = nb->atdat;
    NBParamGpu*         nbp          = nb->nbparam;
    const DeviceStream& deviceStream = *nb->deviceStreams[0];

    auto* feplist = nb->feplist[0].get();

    KernelLaunchConfig fepConfig;
    fepConfig.blockSize[0] = 64;
    fepConfig.blockSize[1] = 1;
    fepConfig.blockSize[2] = 1;

    const int bsize = fepConfig.blockSize[0] * fepConfig.blockSize[1] * fepConfig.blockSize[2];
    // one warp per nri
    const int nriPerBlock      = bsize / warp_size;
    int       nblock           = (feplist->numiAtoms + nriPerBlock - 1) / nriPerBlock;
    fepConfig.gridSize[0]      = nblock;
    fepConfig.gridSize[1]      = 1;
    fepConfig.gridSize[2]      = 1;
    fepConfig.sharedMemorySize = 0;

    const int elecTypeIdx = static_cast<int>(nbp->elecType);
    const int vdwTypeIdx  = static_cast<int>(nbp->vdwType);

    const auto fepKernel     = nbfe_kfunc_ener_ptr[elecTypeIdx][vdwTypeIdx];
    const auto fepKernelArgs = prepareGpuKernelArguments(
            fepKernel, fepConfig, atdat, nbp, feplist, &stepWork.computeVirial);

    launchGpuKernel(fepKernel, fepConfig, deviceStream, nullptr, "k_calc_nb_fep", fepKernelArgs);
}

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

class NonbondedFepGpuTest :
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

    NonbondedFepGpuTest() : checker_(refData_.rootChecker())
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

    void testGpuKernel()
    {
        input_.frHelper.setSoftcoreAlpha(softcoreAlpha_);
        input_.frHelper.setSoftcoreCoulomb(softcoreCoulomb_);
        input_.frHelper.setSoftcoreType(softcoreType_);

        // get forcerec and interaction_const
        t_forcerec fr;
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
            // Currently only supports Nvidia device
            if (testDevice->deviceInfo().deviceVendor == DeviceVendor::Nvidia)
            {
                testDevice->deviceContext().activate();
                SimulationWorkload simulationWorkload;
                auto               deviceStreamManager = std::make_shared<DeviceStreamManager>(
                        testDevice->deviceInfo(), simulationWorkload, false);
                // Construct NbnxmGpu object nb
                NbnxmGpu* nb = gpu_fep_init(*deviceStreamManager, input_.atoms, nbl, x_, fr, lambdas);

                // Run fep gpu kernel
                gpuLaunchFepKernel(nb, stepWork);

                // Transfer results to host buffer
                gpuLaunchCpyback(nb, &output);

                /* Free atdat */
                freeDeviceBuffer(&nb->atdat->xq);
                freeDeviceBuffer(&nb->atdat->f);
                freeDeviceBuffer(&nb->atdat->eLJ);
                freeDeviceBuffer(&nb->atdat->eElec);
                freeDeviceBuffer(&nb->atdat->dvdlLJ);
                freeDeviceBuffer(&nb->atdat->dvdlElec);
                freeDeviceBuffer(&nb->atdat->fShift);
                freeDeviceBuffer(&nb->atdat->shiftVec);

                freeDeviceBuffer(&nb->atdat->q4);
                freeDeviceBuffer(&nb->atdat->atomTypes4);
                freeDeviceBuffer(&nb->nbparam->nbfp);

                delete nb->atdat;
                delete nb->nbparam;
                delete nb;

                checkOutput(&checker_, output);
            }
        }
    }
};

TEST_P(NonbondedFepGpuTest, testGpuKernel)
{
    testGpuKernel();
}

//! configurations to test
std::vector<ListInput> c_interaction = {
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::Cut, VanDerWaalsType::Cut, InteractionModifiers::None) },
    { ListInput(1e-6, 1e-8).setInteraction(CoulombInteractionType::Cut, VanDerWaalsType::Cut, InteractionModifiers::PotSwitch) }
};

//! test parameters
std::vector<real> c_fepLambdas                                  = { 0.0, 0.5, 1.0 };
std::vector<real> c_softcoreBeutlerAlphaOrGapsysLinpointScaling = { 0.0, 0.3 };
std::vector<bool> c_softcoreCoulomb                             = { true, false };
//  No Gapsys Softcore function support on GPU yet
std::vector<SoftcoreType> c_softcoreType = { SoftcoreType::Beutler };

//! Coordinates for testing
std::vector<PaddedVector<RVec>> c_coordinates = {
    { { 1.0, 1.0, 1.0 }, { 1.1, 1.15, 1.2 }, { 0.9, 0.85, 0.8 }, { 1.1, 1.15, 0.8 } }
};

INSTANTIATE_TEST_SUITE_P(NBInteraction,
                         NonbondedFepGpuTest,
                         ::testing::Combine(::testing::ValuesIn(c_softcoreType),
                                            ::testing::ValuesIn(c_interaction),
                                            ::testing::ValuesIn(c_coordinates),
                                            ::testing::ValuesIn(c_fepLambdas),
                                            ::testing::ValuesIn(c_softcoreBeutlerAlphaOrGapsysLinpointScaling),
                                            ::testing::ValuesIn(c_softcoreCoulomb)));


#endif

} // namespace

} // namespace test

} // namespace gmx

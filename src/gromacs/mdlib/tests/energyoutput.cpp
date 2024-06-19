/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief
 * Tests for energy output to log and .edr files.
 *
 * \todo Position and orientation restraints tests.
 * \todo Average and sum in edr file test.
 * \todo AWH output tests.
 * \todo The log and edr outputs are currently saved to the file on the disk and then read
 *       to compare with the reference data. This will be more elegant (and run faster) when we
 *       refactor the output routines to write to a stream interface, which can already be handled
 *       in-memory when running tests.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/mdlib/energyoutput.h"

#include <cassert>
#include <cstdio>

#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/enxio.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/ebin.h"
#include "gromacs/mdlib/makeconstraints.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/unique_cptr.h"

#include "testutils/refdata.h"
#include "testutils/setenv.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! Wraps fclose to discard the return value to use it as a deleter with gmx::unique_cptr.
void fcloseWrapper(FILE* fp)
{
    fclose(fp);
}

/*! \brief Test parameters space.
 *
 * The test will run on a set of combinations of this steucture parameters.
 */
struct EnergyOutputTestParameters
{
    //! Thermostat (enum)
    TemperatureCoupling temperatureCouplingScheme;
    //! Barostat (enum)
    PressureCoupling pressureCouplingScheme;
    //! Integrator
    IntegrationAlgorithm integrator;
    //! Number of saved energy frames (to test averages output).
    int numFrames;
    //! If output should be initialized as a rerun.
    bool isRerun;
    //! Is box triclinic (off-diagonal elements will be printed).
    bool isBoxTriclinic;
};

/*! \brief Sets of parameters on which to run the tests.
 *
 * Only several combinations of the parameters are used. Using all possible combinations will
 * require ~10 MB of test data and ~2 sec to run the tests.
 */
const EnergyOutputTestParameters parametersSets[] = {
    { TemperatureCoupling::No, PressureCoupling::No, IntegrationAlgorithm::MD, 1, false, false },
    { TemperatureCoupling::No, PressureCoupling::No, IntegrationAlgorithm::MD, 1, true, false },
    { TemperatureCoupling::No, PressureCoupling::No, IntegrationAlgorithm::MD, 1, false, true },
    { TemperatureCoupling::No, PressureCoupling::No, IntegrationAlgorithm::MD, 0, false, false },
    { TemperatureCoupling::No, PressureCoupling::No, IntegrationAlgorithm::MD, 10, false, false },
    { TemperatureCoupling::VRescale, PressureCoupling::No, IntegrationAlgorithm::MD, 1, false, false },
    { TemperatureCoupling::NoseHoover, PressureCoupling::No, IntegrationAlgorithm::MD, 1, false, false },
    { TemperatureCoupling::No, PressureCoupling::ParrinelloRahman, IntegrationAlgorithm::MD, 1, false, false },
    { TemperatureCoupling::No, PressureCoupling::Mttk, IntegrationAlgorithm::MD, 1, false, false },
    { TemperatureCoupling::No, PressureCoupling::No, IntegrationAlgorithm::VV, 1, false, false },
    { TemperatureCoupling::No, PressureCoupling::Mttk, IntegrationAlgorithm::VV, 1, false, false }
};

/*! \brief Test fixture to test energy output.
 *
 * The class is initialized to maximize amount of energy terms printed.
 * The test is run for different combinations of temperature and pressure control
 * schemes. Different number of printed steps is also tested.
 */
class EnergyOutputTest : public ::testing::TestWithParam<EnergyOutputTestParameters>
{
    static constexpr int                                      numTempCouplingGroups_ = 3;
    static constexpr std::array<real, numTempCouplingGroups_> tcgInit_{ 0.0_real, 0.0_real, 0.0_real };
    static constexpr real                                     cosAccel_ = 1.0;

public:
    //! File manager
    TestFileManager fileManager_;
    //! Energy (.edr) file
    ener_file_t energyFile_;
    //! Input data
    t_inputrec inputrec_;
    //! Topology
    gmx_mtop_t mtop_;
    //! Simulation time
    double time_;
    //! Total mass
    real tmass_;
    //! Potential energy data
    std::unique_ptr<gmx_enerdata_t> enerdata_;
    //! Kinetic energy data (for temperatures output)
    gmx_ekindata_t ekindata_;
    //! System state
    t_state state_;
    //! PBC box
    matrix box_;
    //! Total virial
    tensor totalVirial_;
    //! Pressure
    tensor pressure_;
    //! Names for the groups.
    std::vector<std::string> groupNameStrings_ = { "Protein", "Water", "Lipid" };
    //! Names for the groups as C strings.
    std::vector<std::vector<char>> groupNameCStrings_;
    //! Handles to the names as C strings in the way needed for SimulationGroups.
    std::vector<char*> groupNameHandles_;
    //! Total dipole momentum
    rvec muTotal_;
    //! Communication record
    t_commrec cr_;
    //! Constraints object (for constraints RMSD output in case of LINCS)
    std::unique_ptr<Constraints> constraints_;
    //! Temporary output filename
    std::string logFilename_;
    //! Temporary energy output filename
    std::string edrFilename_;
    //! Pointer to a temporary output file
    FILE* log_;
    //! Log file wrapper
    unique_cptr<FILE, fcloseWrapper> logFileGuard_;
    //! Reference data
    TestReferenceData refData_;
    //! Checker for reference data
    TestReferenceChecker checker_;

    EnergyOutputTest() :
        ekindata_(tcgInit_, EnsembleTemperatureSetting::NotAvailable, -1.0_real, false, cosAccel_, 1),
        logFilename_(fileManager_.getTemporaryFilePath(".log").string()),
        edrFilename_(fileManager_.getTemporaryFilePath(".edr").string()),
        log_(std::fopen(logFilename_.c_str(), "w")),
        logFileGuard_(log_),
        checker_(refData_.rootChecker())
    {
        // Input record
        inputrec_.delta_t = 0.001;

        // F_RF_EXCL will not be tested - group scheme is not supported any more
        inputrec_.cutoff_scheme = CutoffScheme::Verlet;
        // F_COUL_RECIP
        inputrec_.coulombtype = CoulombInteractionType::Pme;
        // F_LJ_RECIP
        inputrec_.vdwtype = VanDerWaalsType::Pme;

        // F_DVDL_COUL, F_DVDL_VDW, F_DVDL_BONDED, F_DVDL_RESTRAINT, F_DKDL and F_DVDL
        inputrec_.efep = FreeEnergyPerturbationType::Yes;
        inputrec_.fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Coul]      = true;
        inputrec_.fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Vdw]       = true;
        inputrec_.fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Bonded]    = true;
        inputrec_.fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Restraint] = true;
        inputrec_.fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Mass]      = true;
        inputrec_.fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Coul]      = true;
        inputrec_.fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Fep]       = true;

        // F_DISPCORR and F_PDISPCORR
        inputrec_.eDispCorr = DispersionCorrectionType::Ener;
        inputrec_.bRot      = true;

        // F_ECONSERVED
        inputrec_.pressureCouplingOptions.ref_p[YY][XX] = 0.0;
        inputrec_.pressureCouplingOptions.ref_p[ZZ][XX] = 0.0;
        inputrec_.pressureCouplingOptions.ref_p[ZZ][YY] = 0.0;

        // Dipole (mu)
        inputrec_.ewald_geometry = EwaldGeometry::ThreeDC;

        // To print constrain RMSD, constraints algorithm should be set to LINCS.
        inputrec_.eConstrAlg = ConstraintAlgorithm::Lincs;

        mtop_.bIntermolecularInteractions = false;

        // Constructing molecular topology
        gmx_moltype_t molType;

        molType.atoms.nr = 2;

        // F_CONSTR
        // This must be initialized so that Constraints object can be created below.
        InteractionList interactionListConstr;
        interactionListConstr.iatoms.resize(NRAL(F_CONSTR) + 1);
        interactionListConstr.iatoms[0] = 0;
        interactionListConstr.iatoms[1] = 0;
        interactionListConstr.iatoms[2] = 1;
        molType.ilist.at(F_CONSTR)      = interactionListConstr;

        InteractionList interactionListEmpty;
        interactionListEmpty.iatoms.resize(0);
        molType.ilist.at(F_CONSTRNC) = interactionListEmpty;
        molType.ilist.at(F_SETTLE)   = interactionListEmpty;

        // F_LJ14 and F_COUL14
        InteractionList interactionListLJ14;
        interactionListLJ14.iatoms.resize(NRAL(F_LJ14) + 1);
        molType.ilist.at(F_LJ14) = interactionListLJ14;

        // F_LJC14_Q
        InteractionList interactionListLJC14Q;
        interactionListLJC14Q.iatoms.resize(NRAL(F_LJC14_Q) + 1);
        molType.ilist.at(F_LJC14_Q) = interactionListLJC14Q;

        // TODO Do proper initialization for distance and orientation
        //      restraints and remove comments to enable their output
        // F_DISRES
        // InteractionList interactionListDISRES;
        // interactionListDISRES.iatoms.resize(NRAL(F_DISRES) + 1);
        // molType.ilist.at(F_DISRES)   = interactionListDISRES;
        //
        // F_ORIRES
        // InteractionList interactionListORIRES;
        // interactionListORIRES.iatoms.resize(NRAL(F_ORIRES) + 1);
        // molType.ilist.at(F_ORIRES)   = interactionListORIRES;

        mtop_.moltype.push_back(molType);

        gmx_molblock_t molBlock;
        molBlock.type = 0;
        molBlock.nmol = 1;
        mtop_.molblock.push_back(molBlock);

        // This is to keep constraints initialization happy
        mtop_.natoms = 2;
        mtop_.ffparams.iparams.resize(F_NRE);
        mtop_.ffparams.functype.resize(F_NRE);
        mtop_.ffparams.iparams.at(F_CONSTR).constr.dA   = 1.0;
        mtop_.ffparams.iparams.at(F_CONSTR).constr.dB   = 1.0;
        mtop_.ffparams.iparams.at(F_CONSTRNC).constr.dA = 1.0;
        mtop_.ffparams.iparams.at(F_CONSTRNC).constr.dB = 1.0;

        // Groups for energy output, temperature coupling and acceleration
        for (const auto& string : groupNameStrings_)
        {
            std::vector<char> cString(string.begin(), string.end());
            // Need to add null termination
            cString.push_back('\0');
            groupNameCStrings_.emplace_back(cString);
            groupNameHandles_.emplace_back(groupNameCStrings_.back().data());
        }
        for (auto& handle : groupNameHandles_)
        {
            mtop_.groups.groupNames.emplace_back(&handle);
        }

        mtop_.groups.groups[SimulationAtomGroupType::EnergyOutput].resize(numTempCouplingGroups_);
        mtop_.groups.groups[SimulationAtomGroupType::EnergyOutput][0] = 0;
        mtop_.groups.groups[SimulationAtomGroupType::EnergyOutput][1] = 1;
        mtop_.groups.groups[SimulationAtomGroupType::EnergyOutput][2] = 2;

        mtop_.groups.groups[SimulationAtomGroupType::TemperatureCoupling].resize(numTempCouplingGroups_);
        mtop_.groups.groups[SimulationAtomGroupType::TemperatureCoupling][0] = 0;
        mtop_.groups.groups[SimulationAtomGroupType::TemperatureCoupling][1] = 1;
        mtop_.groups.groups[SimulationAtomGroupType::TemperatureCoupling][2] = 2;

        // Nose-Hoover chains
        inputrec_.bPrintNHChains     = true;
        inputrec_.opts.nhchainlength = 2;
        state_.nosehoover_xi.resize(
                mtop_.groups.groups[SimulationAtomGroupType::TemperatureCoupling].size()
                * inputrec_.opts.nhchainlength);
        state_.nosehoover_vxi.resize(
                mtop_.groups.groups[SimulationAtomGroupType::TemperatureCoupling].size()
                * inputrec_.opts.nhchainlength);

        // This will be needed only with MTTK barostat
        state_.nhpres_xi.resize(1 * inputrec_.opts.nhchainlength);
        state_.nhpres_vxi.resize(1 * inputrec_.opts.nhchainlength);

        // Group pairs
        enerdata_ = std::make_unique<gmx_enerdata_t>(
                mtop_.groups.groups[SimulationAtomGroupType::EnergyOutput].size(), nullptr);

        // Kinetic energy and related data
        ekindata_.tcstat.resize(mtop_.groups.groups[SimulationAtomGroupType::TemperatureCoupling].size());

        // This is needed so that the ebin space will be allocated
        inputrec_.cos_accel = cosAccel_;

        // Group options for annealing output
        inputrec_.opts.ngtc = numTempCouplingGroups_;
        snew(inputrec_.opts.annealing, inputrec_.opts.ngtc);
        inputrec_.opts.annealing[0] = SimulatedAnnealing::No;
        inputrec_.opts.annealing[1] = SimulatedAnnealing::Single;
        inputrec_.opts.annealing[2] = SimulatedAnnealing::Periodic;

        // This is to keep done_inputrec happy (otherwise sfree() segfaults)
        snew(inputrec_.opts.anneal_time, inputrec_.opts.ngtc);
        snew(inputrec_.opts.anneal_temp, inputrec_.opts.ngtc);

        // Communication record (for Constraints constructor)
        cr_.nnodes = 1;
        cr_.dd     = nullptr;

        // Constraints object (to get constraints RMSD from)
        // TODO EnergyOutput should not take Constraints object
        // TODO This object will always return zero as RMSD value.
        //      It is more relevant to have non-zero value for testing.
        constraints_ = makeConstraints(
                mtop_, inputrec_, nullptr, false, false, nullptr, &cr_, false, nullptr, nullptr, nullptr, false, nullptr);
    }

    /*! \brief Helper function to generate synthetic data to output
     *
     * \param[in,out] testValue    Base value fr energy data.
     */
    void setStepData(double* testValue)
    {

        time_  = (*testValue += 0.1);
        tmass_ = (*testValue += 0.1);

        enerdata_->term[F_LJ]      = (*testValue += 0.1);
        enerdata_->term[F_COUL_SR] = (*testValue += 0.1);
        enerdata_->term[F_EPOT]    = (*testValue += 0.1);
        enerdata_->term[F_EKIN]    = (*testValue += 0.1);
        enerdata_->term[F_ETOT]    = (*testValue += 0.1);
        enerdata_->term[F_TEMP]    = (*testValue += 0.1);
        enerdata_->term[F_PRES]    = (*testValue += 0.1);

        enerdata_->term[F_BHAM]         = (*testValue += 0.1);
        enerdata_->term[F_EQM]          = (*testValue += 0.1);
        enerdata_->term[F_RF_EXCL]      = (*testValue += 0.1);
        enerdata_->term[F_COUL_RECIP]   = (*testValue += 0.1);
        enerdata_->term[F_LJ_RECIP]     = (*testValue += 0.1);
        enerdata_->term[F_LJ14]         = (*testValue += 0.1);
        enerdata_->term[F_COUL14]       = (*testValue += 0.1);
        enerdata_->term[F_LJC14_Q]      = (*testValue += 0.1);
        enerdata_->term[F_LJC_PAIRS_NB] = (*testValue += 0.1);

        enerdata_->term[F_DVDL_COUL]      = (*testValue += 0.1);
        enerdata_->term[F_DVDL_VDW]       = (*testValue += 0.1);
        enerdata_->term[F_DVDL_BONDED]    = (*testValue += 0.1);
        enerdata_->term[F_DVDL_RESTRAINT] = (*testValue += 0.1);
        enerdata_->term[F_DKDL]           = (*testValue += 0.1);
        enerdata_->term[F_DVDL]           = (*testValue += 0.1);

        enerdata_->term[F_DISPCORR]   = (*testValue += 0.1);
        enerdata_->term[F_PDISPCORR]  = (*testValue += 0.1);
        enerdata_->term[F_DISRESVIOL] = (*testValue += 0.1);
        enerdata_->term[F_ORIRESDEV]  = (*testValue += 0.1);
        enerdata_->term[F_COM_PULL]   = (*testValue += 0.1);
        enerdata_->term[F_ECONSERVED] = (*testValue += 0.1);

        // Group pairs
        for (int i = 0; i < enerdata_->grpp.nener; i++)
        {
            for (int k = 0; k < static_cast<int>(NonBondedEnergyTerms::Count); k++)
            {
                enerdata_->grpp.energyGroupPairTerms[k][i] = (*testValue += 0.1);
            }
        }

        // Kinetic energy and related data
        for (auto& tcstat : ekindata_.tcstat)
        {
            tcstat.T      = (*testValue += 0.1);
            tcstat.lambda = (*testValue += 0.1);
        }
        // Removing constant acceleration removed a total increment of 0.6
        // To avoid unnecessary changes in reference data, we keep the increment
        (*testValue += 0.6);

        // This conditional is to check whether the ebin was allocated.
        // Otherwise it will print cosacc data into the first bin.
        if (inputrec_.cos_accel != 0)
        {
            ekindata_.cosacc.cos_accel = (*testValue += 0.1);
            ekindata_.cosacc.vcos      = (*testValue += 0.1);
        }

        state_.box[XX][XX] = (*testValue += 0.1);
        state_.box[XX][YY] = (*testValue += 0.1);
        state_.box[XX][ZZ] = (*testValue += 0.1);
        state_.box[YY][XX] = (*testValue += 0.1);
        state_.box[YY][YY] = (*testValue += 0.1);
        state_.box[YY][ZZ] = (*testValue += 0.1);
        state_.box[ZZ][XX] = (*testValue += 0.1);
        state_.box[ZZ][YY] = (*testValue += 0.1);
        state_.box[ZZ][ZZ] = (*testValue += 0.1);

        box_[XX][XX] = (*testValue += 0.1);
        box_[XX][YY] = (*testValue += 0.1);
        box_[XX][ZZ] = (*testValue += 0.1);
        box_[YY][XX] = (*testValue += 0.1);
        box_[YY][YY] = (*testValue += 0.1);
        box_[YY][ZZ] = (*testValue += 0.1);
        box_[ZZ][XX] = (*testValue += 0.1);
        box_[ZZ][YY] = (*testValue += 0.1);
        box_[ZZ][ZZ] = (*testValue += 0.1);

        // Removing GMX_CONSTRVIR removed a total increment of 1.8
        // To avoid unnecessary changes in reference data, we keep the increment
        (*testValue += 1.8);

        totalVirial_[XX][XX] = (*testValue += 0.1);
        totalVirial_[XX][YY] = (*testValue += 0.1);
        totalVirial_[XX][ZZ] = (*testValue += 0.1);
        totalVirial_[YY][XX] = (*testValue += 0.1);
        totalVirial_[YY][YY] = (*testValue += 0.1);
        totalVirial_[YY][ZZ] = (*testValue += 0.1);
        totalVirial_[ZZ][XX] = (*testValue += 0.1);
        totalVirial_[ZZ][YY] = (*testValue += 0.1);
        totalVirial_[ZZ][ZZ] = (*testValue += 0.1);

        pressure_[XX][XX] = (*testValue += 0.1);
        pressure_[XX][YY] = (*testValue += 0.1);
        pressure_[XX][ZZ] = (*testValue += 0.1);
        pressure_[YY][XX] = (*testValue += 0.1);
        pressure_[YY][YY] = (*testValue += 0.1);
        pressure_[YY][ZZ] = (*testValue += 0.1);
        pressure_[ZZ][XX] = (*testValue += 0.1);
        pressure_[ZZ][YY] = (*testValue += 0.1);
        pressure_[ZZ][ZZ] = (*testValue += 0.1);

        muTotal_[XX] = (*testValue += 0.1);
        muTotal_[YY] = (*testValue += 0.1);
        muTotal_[ZZ] = (*testValue += 0.1);

        state_.boxv[XX][XX] = (*testValue += 0.1);
        state_.boxv[XX][YY] = (*testValue += 0.1);
        state_.boxv[XX][ZZ] = (*testValue += 0.1);
        state_.boxv[YY][XX] = (*testValue += 0.1);
        state_.boxv[YY][YY] = (*testValue += 0.1);
        state_.boxv[YY][ZZ] = (*testValue += 0.1);
        state_.boxv[ZZ][XX] = (*testValue += 0.1);
        state_.boxv[ZZ][YY] = (*testValue += 0.1);
        state_.boxv[ZZ][ZZ] = (*testValue += 0.1);

        for (int i = 0; i < inputrec_.opts.ngtc; i++)
        {
            *testValue += 0.1;
            ekindata_.setCurrentReferenceTemperature(i, *testValue);
        }

        for (Index k = 0; k < gmx::ssize(mtop_.groups.groups[SimulationAtomGroupType::TemperatureCoupling])
                                      * inputrec_.opts.nhchainlength;
             k++)
        {
            state_.nosehoover_xi[k]  = (*testValue += 0.1);
            state_.nosehoover_vxi[k] = (*testValue += 0.1);
        }
        for (int k = 0; k < inputrec_.opts.nhchainlength; k++)
        {
            state_.nhpres_xi[k]  = (*testValue += 0.1);
            state_.nhpres_vxi[k] = (*testValue += 0.1);
        }
    }

    /*! \brief Check if the contents of the .edr file correspond to the reference data.
     *
     * The code below is based on the 'gmx dump' tool.
     *
     * \param[in] fileName    Name of the file to check.
     * \param[in] frameCount  Number of frames to check.
     */
    void checkEdrFile(const char* fileName, int frameCount)
    {
        ener_file_t  edrFile;
        gmx_enxnm_t* energyTermsEdr = nullptr;
        int          numEnergyTermsEdr;

        edrFile = open_enx(fileName, "r");
        do_enxnms(edrFile, &numEnergyTermsEdr, &energyTermsEdr);
        assert(energyTermsEdr);

        // Check header
        TestReferenceChecker edrFileRef(checker_.checkCompound("File", "EnergyFile"));
        TestReferenceChecker energyTermsRef(
                edrFileRef.checkSequenceCompound("EnergyTerms", numEnergyTermsEdr));
        for (int i = 0; i < numEnergyTermsEdr; i++)
        {
            TestReferenceChecker energyTermRef(energyTermsRef.checkCompound("EnergyTerm", nullptr));
            energyTermRef.checkString(energyTermsEdr[i].name, "Name");
            energyTermRef.checkString(energyTermsEdr[i].unit, "Units");
        }

        // Check frames
        TestReferenceChecker framesRef(edrFileRef.checkSequenceCompound("Frames", frameCount));
        t_enxframe*          frameEdr;
        snew(frameEdr, 1);
        char buffer[22];
        for (int frameId = 0; frameId < frameCount; frameId++)
        {
            bool bCont = do_enx(edrFile, frameEdr);
            EXPECT_TRUE(bCont) << gmx::formatString("Cant read frame %d from .edr file.", frameId);

            TestReferenceChecker frameRef(framesRef.checkCompound("Frame", nullptr));
            frameRef.checkReal(frameEdr->t, "Time");
            frameRef.checkReal(frameEdr->dt, "Timestep");
            frameRef.checkString(gmx_step_str(frameEdr->step, buffer), "Step");
            frameRef.checkString(gmx_step_str(frameEdr->nsum, buffer), "NumSteps");

            EXPECT_EQ(frameEdr->nre, numEnergyTermsEdr)
                    << gmx::formatString("Wrong number of energy terms in frame %d.", frameId);
            TestReferenceChecker energyValuesRef(
                    frameRef.checkSequenceCompound("EnergyTerms", numEnergyTermsEdr));
            for (int i = 0; i < numEnergyTermsEdr; i++)
            {
                TestReferenceChecker energyValueRef(energyValuesRef.checkCompound("EnergyTerm", nullptr));
                energyValueRef.checkString(energyTermsEdr[i].name, "Name");
                energyValueRef.checkReal(frameEdr->ener[i].e, "Value");
            }
        }

        free_enxnms(numEnergyTermsEdr, energyTermsEdr);
        done_ener_file(edrFile);

        free_enxframe(frameEdr);
        sfree(frameEdr);
    }
};

TEST_P(EnergyOutputTest, CheckOutput)
{
    ASSERT_NE(log_, nullptr);
    // Binary output will be written to the temporary location
    energyFile_ = open_enx(edrFilename_.c_str(), "w");
    ASSERT_NE(energyFile_, nullptr);

    EnergyOutputTestParameters parameters = GetParam();
    inputrec_.etc                         = parameters.temperatureCouplingScheme;
    inputrec_.pressureCouplingOptions.epc = parameters.pressureCouplingScheme;
    inputrec_.eI                          = parameters.integrator;

    if (parameters.isBoxTriclinic)
    {
        inputrec_.pressureCouplingOptions.ref_p[YY][XX] = 1.0;
    }

    MDModulesNotifiers            mdModulesNotifiers;
    std::unique_ptr<EnergyOutput> energyOutput =
            std::make_unique<EnergyOutput>(energyFile_,
                                           mtop_,
                                           inputrec_,
                                           nullptr,
                                           nullptr,
                                           parameters.isRerun,
                                           StartingBehavior::NewSimulation,
                                           false,
                                           mdModulesNotifiers);

    // Add synthetic data for a single step
    double testValue = 10.0;
    for (int frame = 0; frame < parameters.numFrames; frame++)
    {
        setStepData(&testValue);
        energyOutput->addDataAtEnergyStep(false,
                                          true,
                                          time_,
                                          tmass_,
                                          enerdata_.get(),
                                          nullptr,
                                          box_,
                                          PTCouplingArrays({ state_.boxv,
                                                             state_.nosehoover_xi,
                                                             state_.nosehoover_vxi,
                                                             state_.nhpres_xi,
                                                             state_.nhpres_vxi }),
                                          state_.fep_state,
                                          totalVirial_,
                                          pressure_,
                                          &ekindata_,
                                          muTotal_,
                                          constraints_.get());

        energyOutput->printAnnealingTemperatures(log_, mtop_.groups, inputrec_.opts, ekindata_);
        energyOutput->printStepToEnergyFile(
                energyFile_, true, false, false, log_, 100 * frame, time_, nullptr, nullptr);
        time_ += 1.0;
    }

    energyOutput->printAnnealingTemperatures(log_, mtop_.groups, inputrec_.opts, ekindata_);
    energyOutput->printAverages(log_, &mtop_.groups);

    // We need to close the file before the contents are available.
    logFileGuard_.reset(nullptr);

    done_ener_file(energyFile_);

    // Set tolerance
    checker_.setDefaultTolerance(relativeToleranceAsFloatingPoint(testValue, 1.0e-5));

    if (parameters.numFrames > 0)
    {
        // Test binary output
        checkEdrFile(edrFilename_.c_str(), parameters.numFrames);
    }

    // Test printed values
    checker_.checkInteger(energyOutput->numEnergyTerms(), "Number of Energy Terms");
    checker_.checkString(TextReader::readFileToString(logFilename_), "log");
}

INSTANTIATE_TEST_SUITE_P(WithParameters, EnergyOutputTest, ::testing::ValuesIn(parametersSets));

} // namespace
} // namespace test
} // namespace gmx

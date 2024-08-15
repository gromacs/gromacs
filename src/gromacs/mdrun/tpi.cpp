/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
 *
 * \brief This file defines the integrator for test particle insertion
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include <cfenv>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/conformation_utilities.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdrunutility/printtime.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcebuffers.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/locality.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/taskassignment/include/gromacs/taskassignment/decidesimulationworkload.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "legacysimulator.h"

struct gmx_edsam;

namespace gmx
{

//! Global max algorithm
static void global_max(t_commrec* cr, int* n)
{
    std::vector<int> sum(cr->nnodes);
    sum[cr->nodeid] = *n;
    gmx_sumi(cr->nnodes, sum.data(), cr);
    *n = *std::max_element(sum.begin(), sum.end());
}

//! Computes and returns the RF exclusion energy for the last molecule starting at \p beginAtom
static real reactionFieldExclusionCorrection(gmx::ArrayRef<const gmx::RVec> x,
                                             const t_mdatoms&               mdatoms,
                                             const interaction_const_t&     ic,
                                             const int                      beginAtom)
{
    real energy = 0;

    for (int i = beginAtom; i < mdatoms.homenr; i++)
    {
        const real qi = mdatoms.chargeA[i];
        energy -= 0.5 * qi * qi * ic.reactionFieldShift;

        for (int j = i + 1; j < mdatoms.homenr; j++)
        {
            const real qj  = mdatoms.chargeA[j];
            const real rsq = distance2(x[i], x[j]);
            energy += qi * qj * (ic.reactionFieldCoefficient * rsq - ic.reactionFieldShift);
        }
    }

    return ic.epsfac * energy;
}

//! The limit in kT for the histogram of insertion energies
static constexpr real sc_bU_bin_limit = 50;
//! The limit in kT for the histogram of insertion energies including the log(volume) term
static constexpr real sc_bU_logV_bin_limit = sc_bU_bin_limit + 10;

//! Class for performing test particle insertions into trajectory frames
class TestParticleInsertion
{
public:
    TestParticleInsertion(const t_inputrec&         inputRec,
                          const gmx_mtop_t&         topGlobal,
                          const gmx_localtop_t&     top,
                          const t_mdatoms&          mdatoms,
                          const MDModulesNotifiers& mdModulesNotifiers,
                          t_forcerec*               forceRec,
                          gmx_enerdata_t*           enerd,
                          const Range<int>&         testAtomsRange,
                          ArrayRef<const RVec>      xMoleculeToInsert,
                          real                      beta,
                          real                      rfExclusionEnergy,
                          real                      referenceVolume,
                          int                       numTasks,
                          int                       taskIndex);

    //! Checks whether all inserted atoms belong to the same energy groups, prints a note when not
    void checkEnergyGroups(ArrayRef<const AtomInfoWithinMoleculeBlock> atomInfoForEachMoleculeBlock,
                           FILE*                                       fpLog) const;

    //! Opens the TPI output file, writes the legends, returns the file pointer
    FILE* openOutputFile(const char* fileName, const gmx_output_env_t* oenv) const;

private:
    /*! \brief Performs a single insertion into a trajecory frame
     *
     * Inserts randomly within a radius of \p rtpi around \p xInit.
     * Rotates the molecule when it consists of multiple atoms.
     *
     * Accumulates the energy term contributions into \p sum_UgembU_
     *
     * \returns the insertion energy and e^(-beta * energy)
     */
    std::pair<double, double> performSingleInsertion(double                 t,
                                                     int64_t                step,
                                                     bool                   hasStateChanged,
                                                     const RVec&            xInit,
                                                     t_state*               stateGlobal,
                                                     MdrunScheduleWorkload* runScheduleWork,
                                                     gmx_wallcycle*         wallCycleCounters,
                                                     t_nrnb*                nrnb);

public:
    //! Performs insertions into a trajecory frame, returns the sum of e^(-beta * energy)
    double insertIntoFrame(double                 t,
                           int64_t                step,
                           int64_t                rerunFrameStep,
                           ArrayRef<const RVec>   rerunX,
                           t_state*               stateGlobal,
                           MdrunScheduleWorkload* runScheduleWork,
                           gmx_wallcycle*         wallCycleCounters,
                           t_nrnb*                nrnb);

    //! Returns the sum over the insertions of separate energy terms time e^-beta*U
    ArrayRef<double> sum_UgembU() { return sum_UgembU_; }

    //! Returns the histogram of insertion energies
    std::vector<double>& bins() { return bins_; }

    //! The inverse bin width in 1/kT of the histogram
    real inverseBinWidth() const { return inverseBinWidth_; }

    //! Returns the number of insertions per block
    int stepBlockSize() const { return stepBlockSize_; }

    //! Returns the energy shift in kT with respect to the reference volume
    real referenceVolumeShift() const { return referenceVolumeShift_; }

private:
    //! The input record
    const t_inputrec& inputRec_;
    //! The global topology
    const gmx_mtop_t& topGlobal_;
    //! The local topology
    const gmx_localtop_t& top_;
    //! The MD-atoms data
    const t_mdatoms& mdatoms_;

    //! Notifiers for MDModules
    const MDModulesNotifiers& mdModulesNotifiers_;

    //! A non-parallel commrec for the energy calculations
    const t_commrec cr_;

    //! The force record
    t_forcerec& fr_;
    //! Energy output
    gmx_enerdata_t& enerd_;
    //! The force buffers
    gmx::ForceBuffers forceBuffers_;

    //! Random number generator with 16 bits internal counter => 2^16 * 2 = 131072 values per stream
    ThreeFry2x64<16> rng_;
    //! A uniform distribution
    UniformRealDistribution<real> dist_;

    //! The index range of the test atoms
    const Range<int> testAtomsRange_;
    //! Whether we insert in a cavity location and not in the whole volume
    bool insertIntoCavity_;
    //! The masses of the atoms that define the cavity
    ArrayRef<real> massesDefiningCavity_;
    //! The coordinates of the molecule to insert, centered at the origin
    std::vector<RVec> xMoleculeToInsert_;

    //! The number of energy groups
    int ngid_;
    //! The energy group index of the molecule to insert
    int gid_tp_;

    //! 1/(k_b T)
    const real beta_;

    //! Whether we use dispersion correction
    bool haveDispCorr_;
    //! Whether the molecule to insert has charges
    bool haveElectrostatics_;
    //! Whether we need to correct for exclusions for reaction-field
    bool haveRFExcl_;
    //! The reaction-field correction energy for exclusions in the inserted molecule
    real rfExclusionEnergy_;

    //! For NPT frames we need to correct for the volume, this is the shift in kT to the reference volume
    real referenceVolumeShift_;

    //! Sum over the insertions of separate energy terms times e^-beta*U
    std::vector<double> sum_UgembU_;

    //! The inverse bin width in 1/kT for the histogram of insertion energies
    const double inverseBinWidth_ = 10;
    //! The histogram of insertion energies
    std::vector<double> bins_;

    //! How many insertions we perform per block
    int stepBlockSize_;

    //! The number of insertion tasks (MPI ranks)
    const int numTasks_;
    //! The index of our task
    const int taskIndex_;

    //! Whether to dump PDB file for insertions with low energies
    bool dumpPdbs_;
    //! The threshold energy for dumping configurations to PDB file
    double dumpEnergyThreshold_;
};

namespace
{

//! Returns whether there are electrostatic contributions to the insertion energy
bool haveElectrostatics(const t_mdatoms& mdatoms, const Range<int>& testAtomsRange)
{
    return std::any_of(testAtomsRange.begin(), testAtomsRange.end(), [mdatoms](int i) {
        return mdatoms.chargeA[i] != 0 || (!mdatoms.chargeB.empty() && mdatoms.chargeB[i] != 0);
    });
}

} // namespace

TestParticleInsertion::TestParticleInsertion(const t_inputrec&         inputRec,
                                             const gmx_mtop_t&         topGlobal,
                                             const gmx_localtop_t&     top,
                                             const t_mdatoms&          mdatoms,
                                             const MDModulesNotifiers& mdModulesNotifiers,
                                             t_forcerec*               forceRec,
                                             gmx_enerdata_t*           enerd,
                                             const Range<int>&         testAtomsRange,
                                             ArrayRef<const RVec>      xMoleculeToInsert,
                                             const real                beta,
                                             const real                rfExclusionEnergy,
                                             const real                referenceVolume,
                                             const int                 numTasks,
                                             const int                 taskIndex) :
    inputRec_(inputRec),
    topGlobal_(topGlobal),
    top_(top),
    mdatoms_(mdatoms),
    mdModulesNotifiers_(mdModulesNotifiers),
    fr_(*forceRec),
    enerd_(*enerd),
    rng_(inputRec.ld_seed, RandomDomain::TestParticleInsertion),
    testAtomsRange_(testAtomsRange),
    insertIntoCavity_(inputRec.eI == IntegrationAlgorithm::TPIC),
    xMoleculeToInsert_(xMoleculeToInsert.begin(), xMoleculeToInsert.end()),
    ngid_(topGlobal.groups.groups[SimulationAtomGroupType::EnergyOutput].size()),
    beta_(beta),
    haveElectrostatics_(haveElectrostatics(mdatoms, testAtomsRange)),
    rfExclusionEnergy_(rfExclusionEnergy),
    referenceVolumeShift_(std::log(referenceVolume)),
    bins_(1),
    numTasks_(numTasks),
    taskIndex_(taskIndex)
{
    const auto& atomInfoInsertion = forceRec->atomInfoForEachMoleculeBlock.back();
    GMX_RELEASE_ASSERT(atomInfoInsertion.indexOfFirstAtomInMoleculeBlock == *testAtomsRange.begin(),
                       "The last molecule block should match the molecule to insert");

    forceBuffers_.resize(topGlobal_.natoms);

    gid_tp_ = atomInfoInsertion.atomInfo[0] & gmx::sc_atomInfo_EnergyGroupIdMask;

    haveDispCorr_ = (inputRec.eDispCorr != DispersionCorrectionType::No);
    haveRFExcl_   = (haveElectrostatics_ && usingRF(inputRec.coulombtype));

    int numEnergyTerms = 1 + ngid_;
    if (haveDispCorr_)
    {
        numEnergyTerms += 1;
    }
    if (haveElectrostatics_)
    {
        numEnergyTerms += ngid_;
        if (haveRFExcl_)
        {
            numEnergyTerms += 1;
        }
        if (usingFullElectrostatics(inputRec.coulombtype))
        {
            numEnergyTerms += 1;
        }
    }

    sum_UgembU_.resize(numEnergyTerms);

    switch (inputRec.eI)
    {
        case IntegrationAlgorithm::TPI: stepBlockSize_ = inputRec.nstlist; break;
        case IntegrationAlgorithm::TPIC: stepBlockSize_ = 1; break;
        default: GMX_RELEASE_ASSERT(false, "Unknown integrator");
    }

    /* The GMX_TPI_DUMP environment variable can be set to dump all configurations
     * to pdb with an insertion energy <= the value of GMX_TPI_DUMP.
     */
    const char* dumpPdbString = getenv("GMX_TPI_DUMP");
    dumpPdbs_                 = (dumpPdbString != nullptr);
    if (dumpPdbs_)
    {
        sscanf(dumpPdbString, "%20lf", &dumpEnergyThreshold_);
    }
}

void TestParticleInsertion::checkEnergyGroups(ArrayRef<const AtomInfoWithinMoleculeBlock> atomInfoForEachMoleculeBlock,
                                              FILE* fpLog) const
{
    const auto& atomInfoInsertion = atomInfoForEachMoleculeBlock.back();
    GMX_RELEASE_ASSERT(atomInfoInsertion.indexOfFirstAtomInMoleculeBlock == *testAtomsRange_.begin(),
                       "The last molecule block should match the molecule to insert");

    for (int a : testAtomsRange_)
    {
        if ((atomInfoInsertion.atomInfo[a - *testAtomsRange_.begin()] & gmx::sc_atomInfo_EnergyGroupIdMask)
            != gid_tp_)
        {
            fprintf(fpLog,
                    "NOTE: Atoms in the molecule to insert belong to different energy groups.\n"
                    "      Only contributions to the group of the first atom will be reported.\n");
            break;
        }
    }
}

FILE* TestParticleInsertion::openOutputFile(const char* fileName, const gmx_output_env_t* oenv) const
{
    FILE* fp = xvgropen(fileName, "TPI energies", "Time (ps)", "(kJ mol\\S-1\\N) / (nm\\S3\\N)", oenv);
    xvgr_subtitle(fp, "f. are averages over one frame", oenv);
    std::vector<std::string> leg;
    leg.emplace_back("-kT log(<Ve\\S-\\betaU\\N>/<V>)");
    leg.emplace_back("f. -kT log<e\\S-\\betaU\\N>");
    leg.emplace_back("f. <e\\S-\\betaU\\N>");
    leg.emplace_back("#f. V");
    leg.emplace_back("f. <Ue\\S-\\betaU\\N>");

    const SimulationGroups& groups = topGlobal_.groups;

    for (int i = 0; i < ngid_; i++)
    {
        leg.emplace_back(gmx::formatString(
                "f. <U\\sVdW %s\\Ne\\S-\\betaU\\N>",
                *(groups.groupNames[groups.groups[SimulationAtomGroupType::EnergyOutput][i]])));
    }
    if (haveDispCorr_)
    {
        leg.emplace_back("f. <U\\sdisp c\\Ne\\S-\\betaU\\N>");
    }
    if (haveElectrostatics_)
    {
        for (int i = 0; i < ngid_; i++)
        {
            leg.emplace_back(gmx::formatString(
                    "f. <U\\sCoul %s\\Ne\\S-\\betaU\\N>",
                    *(groups.groupNames[groups.groups[SimulationAtomGroupType::EnergyOutput][i]])));
        }
        if (haveRFExcl_)
        {
            leg.emplace_back("f. <U\\sRF excl\\Ne\\S-\\betaU\\N>");
        }
        if (usingFullElectrostatics(inputRec_.coulombtype))
        {
            leg.emplace_back("f. <U\\sCoul recip\\Ne\\S-\\betaU\\N>");
        }
    }
    xvgrLegend(fp, leg, oenv);

    return fp;
}

std::pair<double, double> TestParticleInsertion::performSingleInsertion(const double  t,
                                                                        const int64_t step,
                                                                        const bool  hasStateChanged,
                                                                        const RVec& xInit,
                                                                        t_state*    stateGlobal,
                                                                        MdrunScheduleWorkload* runScheduleWork,
                                                                        gmx_wallcycle* wallCycleCounters,
                                                                        t_nrnb*        nrnb)
{
    /* Add random displacement uniformly distributed in a sphere
     * of radius rtpi. We don't need to do this is we generate
     * a new center location every step.
     */
    const real drmax = inputRec_.rtpi;
    RVec       xTestParticle;
    if (insertIntoCavity_ || inputRec_.nstlist > 1)
    {
        /* Generate coordinates within |dx|=drmax of xInit */
        RVec dx;
        do
        {
            for (int d = 0; d < DIM; d++)
            {
                dx[d] = (2 * dist_(rng_) - 1) * drmax;
            }
        } while (norm2(dx) > drmax * drmax);
        xTestParticle = xInit + dx;
    }
    else
    {
        xTestParticle = xInit;
    }

    auto x = makeArrayRef(stateGlobal->x);

    if (testAtomsRange_.size() == 1)
    {
        /* Insert a single atom, just copy the insertion location */
        x[*testAtomsRange_.begin()] = xTestParticle;
    }
    else
    {
        /* Copy the coordinates from the top file */
        for (int i = *testAtomsRange_.begin(); i < *testAtomsRange_.end(); i++)
        {
            x[i] = xMoleculeToInsert_[i - *testAtomsRange_.begin()];
        }
        /* Rotate the molecule randomly */
        real angleX = 2 * M_PI * dist_(rng_);
        /* Draw uniform random number for sin(angleY) instead of angleY itself in order to
         * achieve the uniform distribution in the solid angles space. */
        real angleY = std::asin(2 * dist_(rng_) - 1);
        real angleZ = 2 * M_PI * dist_(rng_);
        rotate_conf(testAtomsRange_.size(),
                    stateGlobal->x.rvec_array() + *testAtomsRange_.begin(),
                    nullptr,
                    angleX,
                    angleY,
                    angleZ);
        /* Shift to the insertion location */
        for (int i : testAtomsRange_)
        {
            x[i] += xTestParticle;
        }
    }

    /* Note: NonLocal refers to the inserted molecule */
    fr_.nbv->convertCoordinates(AtomLocality::NonLocal, x);
    fr_.longRangeNonbondeds->updateAfterPartition(mdatoms_);

    // TPI might place a particle so close that the potential
    // is infinite. Since this is intended to happen, we
    // temporarily suppress any exceptions that the processor
    // might raise, then restore the old behaviour.
    std::fenv_t floatingPointEnvironment;
    std::feholdexcept(&floatingPointEnvironment);

    const int legacyForceFlags =
            GMX_FORCE_NONBONDED | GMX_FORCE_ENERGY | (hasStateChanged ? GMX_FORCE_STATECHANGED : 0);
    runScheduleWork->stepWork = setupStepWorkload(legacyForceFlags,
                                                  inputRec_.mtsLevels,
                                                  step,
                                                  runScheduleWork->domainWork,
                                                  runScheduleWork->simulationWork);

    tensor force_vir;
    clear_mat(force_vir);
    rvec mu_tot;
    do_force(nullptr,
             &cr_,
             nullptr,
             inputRec_,
             mdModulesNotifiers_,
             nullptr,
             nullptr,
             nullptr,
             nullptr,
             step,
             nrnb,
             wallCycleCounters,
             &top_,
             stateGlobal->box,
             stateGlobal->x.arrayRefWithPadding(),
             stateGlobal->v.arrayRefWithPadding().unpaddedArrayRef(),
             &stateGlobal->hist,
             &forceBuffers_.view(),
             force_vir,
             &mdatoms_,
             &enerd_,
             stateGlobal->lambda,
             &fr_,
             *runScheduleWork,
             nullptr,
             mu_tot,
             0.0,
             nullptr,
             fr_.longRangeNonbondeds.get(),
             DDBalanceRegionHandler(nullptr));
    std::feclearexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    std::feupdateenv(&floatingPointEnvironment);

    if (fr_.dispersionCorrection)
    {
        /* Calculate long range corrections to pressure and energy */
        const DispersionCorrection::Correction correction =
                fr_.dispersionCorrection->calculate(stateGlobal->box, 0);
        /* figure out how to rearrange the next 4 lines MRS 8/4/2009 */
        enerd_.term[F_DISPCORR] = correction.energy;
        enerd_.term[F_EPOT] += correction.energy;
        enerd_.term[F_PRES] += correction.pressure;
        enerd_.term[F_DVDL] += correction.dvdl;
    }
    else
    {
        enerd_.term[F_DISPCORR] = 0;
    }
    if (usingRF(fr_.ic->eeltype))
    {
        enerd_.term[F_EPOT] += rfExclusionEnergy_;
    }

    const double epot                = enerd_.term[F_EPOT];
    bool         isEnergyOutOfBounds = false;

    /* If the compiler doesn't optimize this check away
     * we catch the NAN energies.
     * The epot>GMX_REAL_MAX check catches inf values,
     * which should nicely result in embU=0 through the exp below,
     * but it does not hurt to check anyhow.
     */
    /* Non-bonded Interaction usually diverge at r=0.
     * With tabulated interaction functions the first few entries
     * should be capped in a consistent fashion between
     * repulsion, dispersion and Coulomb to avoid accidental
     * negative values in the total energy.
     * The table generation code in tables.c does this.
     * With user tbales the user should take care of this.
     */
    if (epot != epot || epot > GMX_REAL_MAX)
    {
        isEnergyOutOfBounds = TRUE;
    }
    double embU;
    if (isEnergyOutOfBounds)
    {
        if (debug)
        {
            fprintf(debug,
                    "\n  time %.3f, step %d: non-finite energy %f, using exp(-bU)=0\n",
                    t,
                    static_cast<int>(step),
                    epot);
        }
        embU = 0;
    }
    else
    {
        // Exponent argument is fine in SP range, but output can be in DP range
        embU = std::exp(static_cast<double>(-beta_ * epot));
        /* Determine the weighted energy contributions of each energy group */
        int e = 0;
        sum_UgembU_[e++] += epot * embU;
        if (fr_.haveBuckingham)
        {
            for (int i = 0; i < ngid_; i++)
            {
                sum_UgembU_[e++] +=
                        enerd_.grpp.energyGroupPairTerms[NonBondedEnergyTerms::BuckinghamSR][GID(i, gid_tp_, ngid_)]
                        * embU;
            }
        }
        else
        {
            for (int i = 0; i < ngid_; i++)
            {
                sum_UgembU_[e++] +=
                        enerd_.grpp.energyGroupPairTerms[NonBondedEnergyTerms::LJSR][GID(i, gid_tp_, ngid_)]
                        * embU;
            }
        }
        if (haveDispCorr_)
        {
            sum_UgembU_[e++] += enerd_.term[F_DISPCORR] * embU;
        }
        if (haveElectrostatics_)
        {
            for (int i = 0; i < ngid_; i++)
            {
                sum_UgembU_[e++] +=
                        enerd_.grpp.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR][GID(i, gid_tp_, ngid_)]
                        * embU;
            }
            if (haveRFExcl_)
            {
                sum_UgembU_[e++] += rfExclusionEnergy_ * embU;
            }
            if (usingFullElectrostatics(fr_.ic->eeltype))
            {
                sum_UgembU_[e++] += enerd_.term[F_COUL_RECIP] * embU;
            }
        }
    }

    if (debug)
    {
        fprintf(debug,
                "TPI %7d %12.5e %12.5f %12.5f %12.5f\n",
                static_cast<int>(step),
                epot,
                xTestParticle[XX],
                xTestParticle[YY],
                xTestParticle[ZZ]);
    }

    if (dumpPdbs_ && epot <= dumpEnergyThreshold_)
    {
        auto str  = formatString("t%g_step%d.pdb", t, static_cast<int>(step));
        auto str2 = formatString("t: %f step %d ener: %f", t, static_cast<int>(step), epot);
        write_sto_conf_mtop(str.c_str(),
                            str2.c_str(),
                            topGlobal_,
                            stateGlobal->x.rvec_array(),
                            stateGlobal->v.rvec_array(),
                            inputRec_.pbcType,
                            stateGlobal->box);
    }

    return { epot, embU };
}

double TestParticleInsertion::insertIntoFrame(const double           t,
                                              int64_t                step,
                                              const int64_t          rerunFrameStep,
                                              ArrayRef<const RVec>   rerunX,
                                              t_state*               stateGlobal,
                                              MdrunScheduleWorkload* runScheduleWork,
                                              gmx_wallcycle*         wallCycleCounters,
                                              t_nrnb*                nrnb)
{
    /* Copy the coordinates from the input trajectory */
    auto x = makeArrayRef(stateGlobal->x);
    for (Index i = 0; i < rerunX.ssize(); i++)
    {
        x[i] = rerunX[i];
    }
    const matrix& box = stateGlobal->box;

    const double V    = det(box);
    const double logV = std::log(V);

    bool hasStateChanged   = true;
    bool constructPairList = true;

    double sum_embU = 0;
    std::fill(sum_UgembU_.begin(), sum_UgembU_.end(), 0);

    // We initialize xInit to keep the static analyzer happy
    RVec xInit = { 0.0_real, 0.0_real, 0.0_real };
    while (step < inputRec_.nsteps)
    {
        /* Restart random engine using the frame and insertion step
         * as counters.
         * Note that we need to draw several random values per iteration,
         * but by using the internal subcounter functionality of ThreeFry2x64
         * we can draw 131072 unique 64-bit values before exhausting
         * the stream. This is a huge margin, and if something still goes
         * wrong you will get an exception when the stream is exhausted.
         */
        rng_.restart(rerunFrameStep, step);
        dist_.reset(); // erase any memory in the distribution

        if (!insertIntoCavity_)
        {
            /* Random insertion in the whole volume */
            constructPairList = (step % inputRec_.nstlist == 0);
            if (constructPairList)
            {
                /* Generate a random position in the box */
                for (int d = 0; d < DIM; d++)
                {
                    xInit[d] = dist_(rng_) * box[d][d];
                }
            }
        }
        else
        {
            /* Random insertion around a cavity location
             * given by the last coordinate of the trajectory.
             */
            if (step == 0)
            {
                if (massesDefiningCavity_.size() == 1)
                {
                    /* Copy the location of the cavity */
                    xInit = rerunX.back();
                }
                else
                {
                    /* Determine the center of mass of the last molecule */
                    const int numMasses = ssize(massesDefiningCavity_);
                    clear_rvec(xInit);
                    real totalMass = 0;
                    for (int i = 0; i < numMasses; i++)
                    {
                        for (int d = 0; d < DIM; d++)
                        {
                            xInit[d] += massesDefiningCavity_[i]
                                        * rerunX[rerunX.ssize() - numMasses + i][d];
                        }
                        totalMass += massesDefiningCavity_[i];
                    }
                    xInit /= totalMass;
                }
            }
        }

        if (constructPairList)
        {
            for (int a : testAtomsRange_)
            {
                x[a] = xInit;
            }

            /* Put the inserted molecule on it's own search grid */
            fr_.nbv->putAtomsOnGrid(
                    box, 1, xInit, xInit, nullptr, testAtomsRange_, *testAtomsRange_.end(), -1, fr_.atomInfo, x, nullptr);

            /* TODO: Avoid updating all atoms at every bNS step */
            fr_.nbv->setAtomProperties(mdatoms_.typeA, mdatoms_.chargeA, fr_.atomInfo);

            fr_.nbv->constructPairlist(InteractionLocality::Local, top_.excls, step, nrnb);

            constructPairList = FALSE;
        }

        const auto [epot, embU] = performSingleInsertion(
                t, step, hasStateChanged, xInit, stateGlobal, runScheduleWork, wallCycleCounters, nrnb);

        hasStateChanged = false;

        sum_embU += embU;

        if (embU == 0 || beta_ * epot > sc_bU_bin_limit)
        {
            bins_[0]++;
        }
        else
        {
            int i = roundToInt((sc_bU_logV_bin_limit - (beta_ * epot - logV + referenceVolumeShift_))
                               * inverseBinWidth_);
            if (i < 0)
            {
                i = 0;
            }
            bins_.resize(i + 1, 0);
            bins_[i]++;
        }

        step++;
        if ((step / stepBlockSize_) % numTasks_ != taskIndex_)
        {
            /* Skip all steps assigned to the other MPI ranks */
            step += (numTasks_ - 1) * stepBlockSize_;
        }
    }

    return sum_embU;
}

// TODO: Convert to use the nbnxm kernels by putting the system and the teset molecule on two separate search grids
void LegacySimulator::do_tpi()
{
    GMX_RELEASE_ASSERT(gmx_omp_nthreads_get(ModuleMultiThread::Default) == 1,
                       "TPI does not support OpenMP");

    GMX_UNUSED_VALUE(outputProvider_);

    if (usingLJPme(inputRec_->vdwtype))
    {
        gmx_fatal(FARGS, "Test particle insertion not implemented with LJ-PME");
    }
    if (haveEwaldSurfaceContribution(*inputRec_))
    {
        gmx_fatal(FARGS,
                  "TPI with PME currently only works in a 3D geometry with tin-foil "
                  "boundary conditions");
    }

    GMX_LOG(mdLog_.info)
            .asParagraph()
            .appendText(
                    "Note that it is planned to change the command gmx mdrun -tpi "
                    "(and -tpic) to make the functionality available in a different "
                    "form in a future version of GROMACS, e.g. gmx test-particle-insertion.");

    gmx_mtop_generate_local_top(topGlobal_, top_, inputRec_->efep != FreeEnergyPerturbationType::No);

    const bool        insertIntoCavity = (inputRec_->eI == IntegrationAlgorithm::TPIC);
    std::vector<real> massesDefiningCavity;
    if (insertIntoCavity)
    {
        char* ptr = getenv("GMX_TPIC_MASSES");
        if (ptr == nullptr)
        {
            // With a single atom the masses doesn't matter as long as it is !=0
            massesDefiningCavity.push_back(1);
        }
        else
        {
            /* Read (multiple) masses from env var GMX_TPIC_MASSES,
             * The center of mass of the last atoms is then used for TPIC.
             */
            double dbl;
            int    numCharactersRead;
            while (sscanf(ptr, "%20lf%n", &dbl, &numCharactersRead) > 0)
            {
                massesDefiningCavity.push_back(dbl);
                fprintf(fpLog_,
                        "mass[%d] = %f\n",
                        int(gmx::ssize(massesDefiningCavity)) + 1,
                        massesDefiningCavity.back());
                ptr += numCharactersRead;
            }
            if (massesDefiningCavity.empty())
            {
                gmx_fatal(FARGS, "Found zero masses in GMX_TPIC_MASSES");
            }
        }
    }

    /*
       init_em(fplog,TPI,inputrec,&lambda,nrnb,mu_tot,
       state_global->box,fr,mdatoms,top,cr,nfile,fnm,NULL,NULL);*/
    /* We never need full pbc for TPI */
    fr_->pbcType = PbcType::Xyz;
    /* Determine the temperature for the Boltzmann weighting */
    const real temp = constantEnsembleTemperature(*inputRec_);
    if (fpLog_)
    {
        for (int i = 1; i < inputRec_->opts.ngtc; i++)
        {
            if (inputRec_->opts.ref_t[i] != temp)
            {
                fprintf(fpLog_,
                        "\nWARNING: The temperatures of the different temperature coupling groups "
                        "are not identical\n\n");
                fprintf(stderr,
                        "\nWARNING: The temperatures of the different temperature coupling groups "
                        "are not identical\n\n");
            }
        }
        fprintf(fpLog_, "\n  The temperature for test particle insertion is %.3f K\n\n", temp);
    }
    const real beta = 1.0 / (gmx::c_boltz * temp);

    /* Number of insertions per frame */
    const int64_t nsteps = inputRec_->nsteps;

    /* Use the same neighborlist with more insertions points
     * in a sphere of radius drmax around the initial point
     */
    const real drmax = inputRec_->rtpi;

    auto* mdatoms = mdAtoms_->mdatoms();
    atoms2md(topGlobal_, *inputRec_, -1, {}, topGlobal_.natoms, mdAtoms_);
    const double initMassLambda =
            (inputRec_->efep == FreeEnergyPerturbationType::No
                     ? 0.0
                     : inputRec_->fepvals->initialLambda(FreeEnergyPerturbationCouplingType::Mass));
    update_mdatoms(mdatoms, initMassLambda);

    /* Print to log file  */
    walltime_accounting_start_time(wallTimeAccounting_);
    wallcycle_start(wallCycleCounters_, WallCycleCounter::Run);
    print_start(fpLog_, cr_, wallTimeAccounting_, "Test Particle Insertion");

    /* The last molecule is the molecule to be inserted */
    const t_atoms&   atomsToInsert  = topGlobal_.moltype[topGlobal_.molblock.back().type].atoms;
    const Range<int> testAtomsRange = { topGlobal_.natoms - atomsToInsert.nr, topGlobal_.natoms };
    if (debug)
    {
        fprintf(debug, "TPI atoms %d-%d\n", *testAtomsRange.begin(), *testAtomsRange.end());
    }

    auto x = makeArrayRef(stateGlobal_->x);

    if (usingPme(fr_->ic->eeltype))
    {
        gmx_pme_reinit_atoms(fr_->pmedata, *testAtomsRange.begin(), {}, {});
    }

    /* With reacion-field we have distance dependent potentials
     * between excluded atoms, we need to add these separately
     * for the inserted molecule.
     */
    real rfExclusionEnergy = 0;
    if (usingRF(fr_->ic->eeltype))
    {
        rfExclusionEnergy =
                reactionFieldExclusionCorrection(x, *mdatoms, *fr_->ic, *testAtomsRange.begin());
        if (debug)
        {
            fprintf(debug, "RF exclusion correction for inserted molecule: %f kJ/mol\n", rfExclusionEnergy);
        }
    }

    std::vector<RVec> x_mol(testAtomsRange.size());
    for (int a : testAtomsRange)
    {
        /* Copy the coordinates of the molecule to be insterted */
        x_mol[a - *testAtomsRange.begin()] = x[a];
    }

    // Calculate the center of geometry of the molecule to insert
    rvec cog = { 0, 0, 0 };
    for (int a : testAtomsRange)
    {
        rvec_inc(cog, x[a]);
    }
    svmul(1.0_real / testAtomsRange.size(), cog, cog);
    real molRadius = 0;
    for (int a : testAtomsRange)
    {
        molRadius = std::max(molRadius, distance2(x[a], cog));
    }
    molRadius = std::sqrt(molRadius);

    const real maxCutoff = std::max(inputRec_->rvdw, inputRec_->rcoulomb);
    if (insertIntoCavity)
    {
        if (norm(cog) > 0.5 * maxCutoff && fpLog_)
        {
            fprintf(fpLog_, "WARNING: Your TPI molecule is not centered at 0,0,0\n");
            fprintf(stderr, "WARNING: Your TPI molecule is not centered at 0,0,0\n");
        }
    }
    else
    {
        /* Center the molecule to be inserted at zero */
        for (int i = 0; i < testAtomsRange.size(); i++)
        {
            rvec_dec(x_mol[i], cog);
        }
    }

    if (fpLog_)
    {
        fprintf(fpLog_,
                "\nWill insert %d atoms %s partial charges\n",
                testAtomsRange.size(),
                haveElectrostatics(*mdatoms, testAtomsRange) ? "with" : "without");

        fprintf(fpLog_,
                "\nWill insert %" PRId64 " times in each frame of %s\n",
                nsteps,
                opt2fn("-rerun", nFile_, fnm_));
    }

    if (!insertIntoCavity)
    {
        if (inputRec_->nstlist > 1)
        {

            /* With the same pair list we insert in a sphere of radius rtpi  in different orientations */
            if (drmax == 0 && testAtomsRange.size() == 1)
            {
                gmx_fatal(FARGS,
                          "Re-using the neighborlist %d times for insertions of a single atom in a "
                          "sphere of radius %f does not make sense",
                          inputRec_->nstlist,
                          drmax);
            }
            if (fpLog_)
            {
                fprintf(fpLog_,
                        "Will use the same neighborlist for %d insertions in a sphere of radius "
                        "%f\n",
                        inputRec_->nstlist,
                        drmax);
            }
        }
    }
    else
    {
        if (fpLog_)
        {
            fprintf(fpLog_,
                    "Will insert randomly in a sphere of radius %f around the center of the "
                    "cavity\n",
                    drmax);
        }
    }

    /* With the same pair list we insert in a sphere of radius rtpi
     * in different orientations. We generate the pairlist with all
     * inserted atoms located in the center of the sphere, so we need
     * a buffer of size of the sphere and molecule radius.
     */
    fr_->rlist = maxCutoff + inputRec_->rtpi + molRadius;
    fr_->nbv->changePairlistRadii(fr_->rlist, fr_->rlist);

    // We don't compute listed forces, set up empty interaction lists
    const gmx_ffparams_t         emptyFFParams;
    const InteractionDefinitions emptyInteractionDefinitions(emptyFFParams);
    for (auto& listedForces : fr_->listedForces)
    {
        listedForces.setup(emptyInteractionDefinitions, 0, false);
    }

    double V_all     = 0;
    double VembU_all = 0;

    /* Avoid frame step numbers <= -1 */
    int64_t frame_step_prev = -1;

    t_trxstatus* status;
    t_trxframe   rerun_fr;
    bool         isNotLastFrame =
            read_first_frame(oenv_, &status, opt2fn("-rerun", nFile_, fnm_), &rerun_fr, TRX_NEED_X);
    int frame = 0;

    if (rerun_fr.natoms - (insertIntoCavity ? gmx::ssize(massesDefiningCavity) : 0)
        != mdatoms->nr - testAtomsRange.size())
    {
        gmx_fatal(FARGS,
                  "Number of atoms in trajectory (%d)%s "
                  "is not equal the number in the run input file (%d) "
                  "minus the number of atoms to insert (%d)\n",
                  rerun_fr.natoms,
                  insertIntoCavity ? " minus one" : "",
                  mdatoms->nr,
                  testAtomsRange.size());
    }

    TestParticleInsertion tpi(*inputRec_,
                              topGlobal_,
                              *top_,
                              *mdatoms,
                              mdModulesNotifiers_,
                              fr_,
                              enerd_,
                              testAtomsRange,
                              x_mol,
                              beta,
                              rfExclusionEnergy,
                              det(rerun_fr.box),
                              cr_->nnodes,
                              cr_->nodeid);

    tpi.checkEnergyGroups(fr_->atomInfoForEachMoleculeBlock, fpLog_);

    FILE* fp_tpi = nullptr;
    if (MAIN(cr_))
    {
        fp_tpi = tpi.openOutputFile(opt2fn("-tpi", nFile_, fnm_), oenv_);
    }

    while (isNotLastFrame)
    {
        int64_t frame_step = rerun_fr.step;
        if (frame_step <= frame_step_prev)
        {
            /* We don't have step number in the trajectory file,
             * or we have constant or decreasing step numbers.
             * Ensure we have increasing step numbers, since we use
             * the step numbers as a counter for random numbers.
             */
            frame_step = frame_step_prev + 1;
        }
        frame_step_prev = frame_step;

        /* Copy the coordinates from the input trajectory */
        auto x = makeArrayRef(stateGlobal_->x);
        for (int i = 0; i < rerun_fr.natoms; i++)
        {
            copy_rvec(rerun_fr.x[i], x[i]);
        }
        copy_mat(rerun_fr.box, stateGlobal_->box);
        const matrix& box    = stateGlobal_->box;
        const double  volume = det(box);

        put_atoms_in_box(fr_->pbcType, box, x);

        /* Put all atoms except for the inserted ones on the grid */
        rvec vzero       = { 0, 0, 0 };
        rvec boxDiagonal = { box[XX][XX], box[YY][YY], box[ZZ][ZZ] };
        fr_->nbv->putAtomsOnGrid(box,
                                 0,
                                 vzero,
                                 boxDiagonal,
                                 nullptr,
                                 { 0, *testAtomsRange.begin() },
                                 *testAtomsRange.begin(),
                                 -1,
                                 fr_->atomInfo,
                                 x,
                                 nullptr);

        gmx_edsam* const ed = nullptr;

        // TPI does not support DD so we only call this once, on the first step
        GMX_ASSERT(runScheduleWork_->simulationWork.havePpDomainDecomposition == false,
                   "We should not be using PP domain decomposition here");
        // TPI only computes non-bonded interaction energies, no other energies should be computed.

        // Note that lot of fr_ internal data (such as bondeds) is not fully set up, and will
        // not be set up later, because TPI runs only uses a narrow subset of functionality.
        // TPI also uses a different definition of local an non-local atoms from the rest of the code,
        // so care needs to be taken that members of domainWork get correctly initialized
        // for the TPI use-case.
        runScheduleWork_->domainWork = setupDomainLifetimeWorkload(
                *inputRec_, *fr_, pullWork_, ed, *mdatoms, runScheduleWork_->simulationWork);

        const int64_t step = cr_->nodeid * tpi.stepBlockSize();

        double sum_embU = tpi.insertIntoFrame(
                rerun_fr.time,
                step,
                rerun_fr.step,
                constArrayRefFromArray<RVec>(reinterpret_cast<RVec*>(rerun_fr.x), rerun_fr.natoms),
                stateGlobal_,
                runScheduleWork_,
                wallCycleCounters_,
                nrnb_);

        auto sum_UgembU = tpi.sum_UgembU();

        if (PAR(cr_))
        {
            /* When running in parallel sum the energies over the processes */
            gmx_sumd(1, &sum_embU, cr_);
            gmx_sumd(gmx::ssize(sum_UgembU), sum_UgembU.data(), cr_);
        }

        frame++;
        V_all += volume;
        VembU_all += volume * sum_embU / nsteps;

        if (fp_tpi)
        {
            if (mdrunOptions_.verbose || frame % 10 == 0 || frame < 10)
            {
                fprintf(stderr,
                        "mu %10.3e <mu> %10.3e\n",
                        -std::log(sum_embU / nsteps) / beta,
                        -std::log(VembU_all / V_all) / beta);
            }

            fprintf(fp_tpi,
                    "%10.3f %12.5e %12.5e %12.5e %12.5e",
                    rerun_fr.time,
                    VembU_all == 0 ? 20 / beta : -std::log(VembU_all / V_all) / beta,
                    sum_embU == 0 ? 20 / beta : -std::log(sum_embU / nsteps) / beta,
                    sum_embU / nsteps,
                    volume);
            for (double e : sum_UgembU)
            {
                fprintf(fp_tpi, " %12.5e", e / nsteps);
            }
            fprintf(fp_tpi, "\n");
            fflush(fp_tpi);
        }

        isNotLastFrame = read_next_frame(oenv_, status, &rerun_fr);
    } /* End of the loop  */
    walltime_accounting_end_time(wallTimeAccounting_);

    close_trx(status);

    if (fp_tpi != nullptr)
    {
        xvgrclose(fp_tpi);
    }

    if (fpLog_ != nullptr)
    {
        fprintf(fpLog_, "\n");
        fprintf(fpLog_, "  <V>  = %12.5e nm^3\n", V_all / frame);
        const double mu = -std::log(VembU_all / V_all) / beta;
        fprintf(fpLog_, "  <mu> = %12.5e kJ/mol\n", mu);

        if (!std::isfinite(mu))
        {
            fprintf(fpLog_,
                    "\nThe computed chemical potential is not finite - consider increasing the "
                    "number of steps and/or the number of frames to insert into.\n");
        }
    }

    auto& bins = tpi.bins();

    /* Write the Boltzmann factor histogram */
    if (PAR(cr_))
    {
        /* When running in parallel sum the bins over the processes */
        int i = gmx::ssize(bins);
        global_max(cr_, &i);
        bins.resize(i);
        gmx_sumd(gmx::ssize(bins), bins.data(), cr_);
    }
    if (MAIN(cr_))
    {
        fp_tpi   = xvgropen(opt2fn("-tpid", nFile_, fnm_),
                          "TPI energy distribution",
                          "\\betaU - log(V/<V>)",
                          "count",
                          oenv_);
        auto str = gmx::formatString("number \\betaU > %g: %9.3e", sc_bU_bin_limit, bins[0]);
        xvgr_subtitle(fp_tpi, str.c_str(), oenv_);
        std::array<std::string, 2> tpid_leg = { "direct", "reweighted" };
        xvgrLegend(fp_tpi, tpid_leg, oenv_);
        for (int i = gmx::ssize(bins) - 1; i > 0; i--)
        {
            double bUlogV = -i / tpi.inverseBinWidth() + sc_bU_logV_bin_limit
                            - tpi.referenceVolumeShift() + std::log(V_all / frame);
            fprintf(fp_tpi,
                    "%6.2f %10d %12.5e\n",
                    bUlogV,
                    roundToInt(bins[i]),
                    bins[i] * std::exp(-bUlogV) * V_all / VembU_all);
        }
        xvgrclose(fp_tpi);
    }

    walltime_accounting_set_nsteps_done(wallTimeAccounting_, frame * inputRec_->nsteps);
}

} // namespace gmx

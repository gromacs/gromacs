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
#include "gmxpre.h"

#include "config.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>

#include <array>
#include <optional>

#include "gromacs/applied_forces/awh/awh.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_coordinate_receiver_gpu.h"
#include "gromacs/ewald/pme_pp.h"
#include "gromacs/ewald/pme_pp_comm_gpu.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nonbonded/nb_free_energy.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/imd/imd.h"
#include "gromacs/listed_forces/disre.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/listed_forces/orires.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/calcmu.h"
#include "gromacs/mdlib/calcvir.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdlib/enerdata_utils.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdlib/wall.h"
#include "gromacs/mdlib/wholemoleculetransform.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcebuffers.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/state_propagator_data_gpu.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_rotation.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/wallcyclereporting.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/fixedcapacityvector.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

#include "gpuforcereduction.h"

using gmx::ArrayRef;
using gmx::AtomLocality;
using gmx::DomainLifetimeWorkload;
using gmx::ForceOutputs;
using gmx::ForceWithShiftForces;
using gmx::InteractionLocality;
using gmx::RVec;
using gmx::SimulationWorkload;
using gmx::StepWorkload;

// TODO: this environment variable allows us to verify before release
// that on less common architectures the total cost of polling is not larger than
// a blocking wait (so polling does not introduce overhead when the static
// PME-first ordering would suffice).
static const bool c_disableAlternatingWait = (getenv("GMX_DISABLE_ALTERNATING_GPU_WAIT") != nullptr);

static void sum_forces(ArrayRef<RVec> f, ArrayRef<const RVec> forceToAdd)
{
    GMX_ASSERT(f.size() >= forceToAdd.size(), "Accumulation buffer should be sufficiently large");
    const int end = forceToAdd.size();

    int gmx_unused nt = gmx_omp_nthreads_get(ModuleMultiThread::Default);
#pragma omp parallel for num_threads(nt) schedule(static)
    for (int i = 0; i < end; i++)
    {
        rvec_inc(f[i], forceToAdd[i]);
    }
}

static void calc_virial(int                              start,
                        int                              homenr,
                        const rvec                       x[],
                        const gmx::ForceWithShiftForces& forceWithShiftForces,
                        tensor                           vir_part,
                        const matrix                     box,
                        t_nrnb*                          nrnb,
                        const t_forcerec*                fr,
                        PbcType                          pbcType)
{
    /* The short-range virial from surrounding boxes */
    const rvec* fshift          = as_rvec_array(forceWithShiftForces.shiftForces().data());
    const rvec* shiftVecPointer = as_rvec_array(fr->shift_vec.data());
    calc_vir(gmx::c_numShiftVectors, shiftVecPointer, fshift, vir_part, pbcType == PbcType::Screw, box);
    inc_nrnb(nrnb, eNR_VIRIAL, gmx::c_numShiftVectors);

    /* Calculate partial virial, for local atoms only, based on short range.
     * Total virial is computed in global_stat, called from do_md
     */
    const rvec* f = as_rvec_array(forceWithShiftForces.force().data());
    f_calc_vir(start, start + homenr, x, f, vir_part, box);
    inc_nrnb(nrnb, eNR_VIRIAL, homenr);

    if (debug)
    {
        pr_rvecs(debug, 0, "vir_part", vir_part, DIM);
    }
}

static void pull_potential_wrapper(const t_commrec*               cr,
                                   const t_inputrec&              ir,
                                   const matrix                   box,
                                   gmx::ArrayRef<const gmx::RVec> x,
                                   const t_mdatoms*               mdatoms,
                                   gmx_enerdata_t*                enerd,
                                   pull_t*                        pull_work,
                                   const real*                    lambda,
                                   double                         t,
                                   gmx_wallcycle*                 wcycle)
{
    t_pbc pbc;
    real  dvdl;

    /* Calculate the center of mass forces, this requires communication,
     * which is why pull_potential is called close to other communication.
     */
    wallcycle_start(wcycle, WallCycleCounter::PullPot);
    set_pbc(&pbc, ir.pbcType, box);
    dvdl = 0;
    enerd->term[F_COM_PULL] +=
            pull_potential(pull_work,
                           mdatoms->massT,
                           pbc,
                           cr,
                           t,
                           lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Restraint)],
                           x,
                           &dvdl);
    enerd->dvdl_lin[FreeEnergyPerturbationCouplingType::Restraint] += dvdl;
    wallcycle_stop(wcycle, WallCycleCounter::PullPot);
}

static void pme_receive_force_ener(t_forcerec*           fr,
                                   const t_commrec*      cr,
                                   gmx::ForceWithVirial* forceWithVirial,
                                   gmx_enerdata_t*       enerd,
                                   bool                  useGpuPmePpComms,
                                   bool                  receivePmeForceToGpu,
                                   gmx_wallcycle*        wcycle)
{
    real  e_q, e_lj, dvdl_q, dvdl_lj;
    float cycles_ppdpme, cycles_seppme;

    cycles_ppdpme = wallcycle_stop(wcycle, WallCycleCounter::PpDuringPme);
    dd_cycles_add(cr->dd, cycles_ppdpme, ddCyclPPduringPME);

    /* In case of node-splitting, the PP nodes receive the long-range
     * forces, virial and energy from the PME nodes here.
     */
    wallcycle_start(wcycle, WallCycleCounter::PpPmeWaitRecvF);
    dvdl_q  = 0;
    dvdl_lj = 0;
    gmx_pme_receive_f(fr->pmePpCommGpu.get(),
                      cr,
                      forceWithVirial,
                      &e_q,
                      &e_lj,
                      &dvdl_q,
                      &dvdl_lj,
                      useGpuPmePpComms,
                      receivePmeForceToGpu,
                      &cycles_seppme);
    enerd->term[F_COUL_RECIP] += e_q;
    enerd->term[F_LJ_RECIP] += e_lj;
    enerd->dvdl_lin[FreeEnergyPerturbationCouplingType::Coul] += dvdl_q;
    enerd->dvdl_lin[FreeEnergyPerturbationCouplingType::Vdw] += dvdl_lj;

    if (wcycle)
    {
        dd_cycles_add(cr->dd, cycles_seppme, ddCyclPME);
    }
    wallcycle_stop(wcycle, WallCycleCounter::PpPmeWaitRecvF);
}

static void print_large_forces(FILE*                fp,
                               const t_mdatoms*     md,
                               const t_commrec*     cr,
                               int64_t              step,
                               real                 forceTolerance,
                               ArrayRef<const RVec> x,
                               ArrayRef<const RVec> f)
{
    real       force2Tolerance = gmx::square(forceTolerance);
    gmx::Index numNonFinite    = 0;
    for (int i = 0; i < md->homenr; i++)
    {
        real force2    = norm2(f[i]);
        bool nonFinite = !std::isfinite(force2);
        if (force2 >= force2Tolerance || nonFinite)
        {
            fprintf(fp,
                    "step %" PRId64 " atom %6d  x %8.3f %8.3f %8.3f  force %12.5e\n",
                    step,
                    ddglatnr(cr->dd, i),
                    x[i][XX],
                    x[i][YY],
                    x[i][ZZ],
                    std::sqrt(force2));
        }
        if (nonFinite)
        {
            numNonFinite++;
        }
    }
    if (numNonFinite > 0)
    {
        /* Note that with MPI this fatal call on one rank might interrupt
         * the printing on other ranks. But we can only avoid that with
         * an expensive MPI barrier that we would need at each step.
         */
        gmx_fatal(FARGS, "At step %" PRId64 " detected non-finite forces on %td atoms", step, numNonFinite);
    }
}

//! When necessary, spreads forces on vsites and computes the virial for \p forceOutputs->forceWithShiftForces()
static void postProcessForceWithShiftForces(t_nrnb*                   nrnb,
                                            gmx_wallcycle*            wcycle,
                                            const matrix              box,
                                            ArrayRef<const RVec>      x,
                                            ForceOutputs*             forceOutputs,
                                            tensor                    vir_force,
                                            const t_mdatoms&          mdatoms,
                                            const t_forcerec&         fr,
                                            gmx::VirtualSitesHandler* vsite,
                                            const StepWorkload&       stepWork)
{
    ForceWithShiftForces& forceWithShiftForces = forceOutputs->forceWithShiftForces();

    /* If we have NoVirSum forces, but we do not calculate the virial,
     * we later sum the forceWithShiftForces buffer together with
     * the noVirSum buffer and spread the combined vsite forces at once.
     */
    if (vsite && (!forceOutputs->haveForceWithVirial() || stepWork.computeVirial))
    {
        using VirialHandling = gmx::VirtualSitesHandler::VirialHandling;

        auto                 f      = forceWithShiftForces.force();
        auto                 fshift = forceWithShiftForces.shiftForces();
        const VirialHandling virialHandling =
                (stepWork.computeVirial ? VirialHandling::Pbc : VirialHandling::None);
        vsite->spreadForces(x, f, virialHandling, fshift, nullptr, nrnb, box, wcycle);
        forceWithShiftForces.haveSpreadVsiteForces() = true;
    }

    if (stepWork.computeVirial)
    {
        /* Calculation of the virial must be done after vsites! */
        calc_virial(
                0, mdatoms.homenr, as_rvec_array(x.data()), forceWithShiftForces, vir_force, box, nrnb, &fr, fr.pbcType);
    }
}

//! Spread, compute virial for and sum forces, when necessary
static void postProcessForces(const t_commrec*          cr,
                              int64_t                   step,
                              t_nrnb*                   nrnb,
                              gmx_wallcycle*            wcycle,
                              const matrix              box,
                              ArrayRef<const RVec>      x,
                              ForceOutputs*             forceOutputs,
                              tensor                    vir_force,
                              const t_mdatoms*          mdatoms,
                              const t_forcerec*         fr,
                              gmx::VirtualSitesHandler* vsite,
                              const StepWorkload&       stepWork)
{
    // Extract the final output force buffer, which is also the buffer for forces with shift forces
    ArrayRef<RVec> f = forceOutputs->forceWithShiftForces().force();

    if (forceOutputs->haveForceWithVirial())
    {
        auto& forceWithVirial = forceOutputs->forceWithVirial();

        if (vsite)
        {
            /* Spread the mesh force on virtual sites to the other particles...
             * This is parallellized. MPI communication is performed
             * if the constructing atoms aren't local.
             */
            GMX_ASSERT(!stepWork.computeVirial || f.data() != forceWithVirial.force_.data(),
                       "We need separate force buffers for shift and virial forces when "
                       "computing the virial");
            GMX_ASSERT(!stepWork.computeVirial
                               || forceOutputs->forceWithShiftForces().haveSpreadVsiteForces(),
                       "We should spread the force with shift forces separately when computing "
                       "the virial");
            const gmx::VirtualSitesHandler::VirialHandling virialHandling =
                    (stepWork.computeVirial ? gmx::VirtualSitesHandler::VirialHandling::NonLinear
                                            : gmx::VirtualSitesHandler::VirialHandling::None);
            matrix virial = { { 0 } };
            vsite->spreadForces(x, forceWithVirial.force_, virialHandling, {}, virial, nrnb, box, wcycle);
            forceWithVirial.addVirialContribution(virial);
        }

        if (stepWork.computeVirial)
        {
            /* Now add the forces, this is local */
            sum_forces(f, forceWithVirial.force_);

            /* Add the direct virial contributions */
            GMX_ASSERT(
                    forceWithVirial.computeVirial_,
                    "forceWithVirial should request virial computation when we request the virial");
            m_add(vir_force, forceWithVirial.getVirial(), vir_force);

            if (debug)
            {
                pr_rvecs(debug, 0, "vir_force", vir_force, DIM);
            }
        }
    }
    else
    {
        GMX_ASSERT(vsite == nullptr || forceOutputs->forceWithShiftForces().haveSpreadVsiteForces(),
                   "We should have spread the vsite forces (earlier)");
    }

    if (fr->print_force >= 0)
    {
        print_large_forces(stderr, mdatoms, cr, step, fr->print_force, x, f);
    }
}

static void do_nb_verlet(t_forcerec*                fr,
                         const interaction_const_t* ic,
                         gmx_enerdata_t*            enerd,
                         const StepWorkload&        stepWork,
                         const InteractionLocality  ilocality,
                         const int                  clearF,
                         const int64_t              step,
                         t_nrnb*                    nrnb,
                         gmx_wallcycle*             wcycle)
{
    if (!stepWork.computeNonbondedForces)
    {
        /* skip non-bonded calculation */
        return;
    }

    nonbonded_verlet_t* nbv = fr->nbv.get();

    /* GPU kernel launch overhead is already timed separately */
    if (!nbv->useGpu())
    {
        /* When dynamic pair-list  pruning is requested, we need to prune
         * at nstlistPrune steps.
         */
        if (nbv->isDynamicPruningStepCpu(step))
        {
            /* Prune the pair-list beyond fr->ic->rlistPrune using
             * the current coordinates of the atoms.
             */
            wallcycle_sub_start(wcycle, WallCycleSubCounter::NonbondedPruning);
            nbv->dispatchPruneKernelCpu(ilocality, fr->shift_vec);
            wallcycle_sub_stop(wcycle, WallCycleSubCounter::NonbondedPruning);
        }
    }

    nbv->dispatchNonbondedKernel(
            ilocality,
            *ic,
            stepWork,
            clearF,
            fr->shift_vec,
            enerd->grpp.energyGroupPairTerms[fr->haveBuckingham ? NonBondedEnergyTerms::BuckinghamSR
                                                                : NonBondedEnergyTerms::LJSR],
            enerd->grpp.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR],
            nrnb);
}

static inline void clearRVecs(ArrayRef<RVec> v, const bool useOpenmpThreading)
{
    int nth = gmx_omp_nthreads_get_simple_rvec_task(ModuleMultiThread::Default, v.ssize());

    /* Note that we would like to avoid this conditional by putting it
     * into the omp pragma instead, but then we still take the full
     * omp parallel for overhead (at least with gcc5).
     */
    if (!useOpenmpThreading || nth == 1)
    {
        for (RVec& elem : v)
        {
            clear_rvec(elem);
        }
    }
    else
    {
#pragma omp parallel for num_threads(nth) schedule(static)
        for (gmx::Index i = 0; i < v.ssize(); i++)
        {
            clear_rvec(v[i]);
        }
    }
}

/*! \brief Return an estimate of the average kinetic energy or 0 when unreliable
 *
 * \param groupOptions  Group options, containing T-coupling options
 */
static real averageKineticEnergyEstimate(const t_grpopts& groupOptions)
{
    real nrdfCoupled   = 0;
    real nrdfUncoupled = 0;
    real kineticEnergy = 0;
    for (int g = 0; g < groupOptions.ngtc; g++)
    {
        if (groupOptions.tau_t[g] >= 0)
        {
            nrdfCoupled += groupOptions.nrdf[g];
            kineticEnergy += groupOptions.nrdf[g] * 0.5 * groupOptions.ref_t[g] * gmx::c_boltz;
        }
        else
        {
            nrdfUncoupled += groupOptions.nrdf[g];
        }
    }

    /* This conditional with > also catches nrdf=0 */
    if (nrdfCoupled > nrdfUncoupled)
    {
        return kineticEnergy * (nrdfCoupled + nrdfUncoupled) / nrdfCoupled;
    }
    else
    {
        return 0;
    }
}

/*! \brief This routine checks that the potential energy is finite.
 *
 * Always checks that the potential energy is finite. If step equals
 * inputrec.init_step also checks that the magnitude of the potential energy
 * is reasonable. Terminates with a fatal error when a check fails.
 * Note that passing this check does not guarantee finite forces,
 * since those use slightly different arithmetics. But in most cases
 * there is just a narrow coordinate range where forces are not finite
 * and energies are finite.
 *
 * \param[in] step      The step number, used for checking and printing
 * \param[in] enerd     The energy data; the non-bonded group energies need to be added to
 *                      \c enerd.term[F_EPOT] before calling this routine
 * \param[in] inputrec  The input record
 */
static void checkPotentialEnergyValidity(int64_t step, const gmx_enerdata_t& enerd, const t_inputrec& inputrec)
{
    /* Threshold valid for comparing absolute potential energy against
     * the kinetic energy. Normally one should not consider absolute
     * potential energy values, but with a factor of one million
     * we should never get false positives.
     */
    constexpr real c_thresholdFactor = 1e6;

    bool energyIsNotFinite    = !std::isfinite(enerd.term[F_EPOT]);
    real averageKineticEnergy = 0;
    /* We only check for large potential energy at the initial step,
     * because that is by far the most likely step for this too occur
     * and because computing the average kinetic energy is not free.
     * Note: nstcalcenergy >> 1 often does not allow to catch large energies
     * before they become NaN.
     */
    if (step == inputrec.init_step && EI_DYNAMICS(inputrec.eI))
    {
        averageKineticEnergy = averageKineticEnergyEstimate(inputrec.opts);
    }

    if (energyIsNotFinite
        || (averageKineticEnergy > 0 && enerd.term[F_EPOT] > c_thresholdFactor * averageKineticEnergy))
    {
        GMX_THROW(gmx::InternalError(gmx::formatString(
                "Step %" PRId64
                ": The total potential energy is %g, which is %s. The LJ and electrostatic "
                "contributions to the energy are %g and %g, respectively. A %s potential energy "
                "can be caused by overlapping interactions in bonded interactions or very large%s "
                "coordinate values. Usually this is caused by a badly- or non-equilibrated initial "
                "configuration, incorrect interactions or parameters in the topology.",
                step,
                enerd.term[F_EPOT],
                energyIsNotFinite ? "not finite" : "extremely high",
                enerd.term[F_LJ],
                enerd.term[F_COUL_SR],
                energyIsNotFinite ? "non-finite" : "very high",
                energyIsNotFinite ? " or Nan" : "")));
    }
}

/*! \brief Compute forces and/or energies for special algorithms
 *
 * The intention is to collect all calls to algorithms that compute
 * forces on local atoms only and that do not contribute to the local
 * virial sum (but add their virial contribution separately).
 * Eventually these should likely all become ForceProviders.
 * Within this function the intention is to have algorithms that do
 * global communication at the end, so global barriers within the MD loop
 * are as close together as possible.
 *
 * \param[in]     fplog            The log file
 * \param[in]     cr               The communication record
 * \param[in]     inputrec         The input record
 * \param[in]     awh              The Awh module (nullptr if none in use).
 * \param[in]     enforcedRotation Enforced rotation module.
 * \param[in]     imdSession       The IMD session
 * \param[in]     pull_work        The pull work structure.
 * \param[in]     step             The current MD step
 * \param[in]     t                The current time
 * \param[in,out] wcycle           Wallcycle accounting struct
 * \param[in,out] forceProviders   Pointer to a list of force providers
 * \param[in]     box              The unit cell
 * \param[in]     x                The coordinates
 * \param[in]     mdatoms          Per atom properties
 * \param[in]     lambda           Array of free-energy lambda values
 * \param[in]     stepWork         Step schedule flags
 * \param[in,out] forceWithVirialMtsLevel0  Force and virial for MTS level0 forces
 * \param[in,out] forceWithVirialMtsLevel1  Force and virial for MTS level1 forces, can be nullptr
 * \param[in,out] enerd            Energy buffer
 * \param[in,out] ed               Essential dynamics pointer
 * \param[in]     didNeighborSearch Tells if we did neighbor searching this step, used for ED sampling
 *
 * \todo Remove didNeighborSearch, which is used incorrectly.
 * \todo Convert all other algorithms called here to ForceProviders.
 */
static void computeSpecialForces(FILE*                          fplog,
                                 const t_commrec*               cr,
                                 const t_inputrec&              inputrec,
                                 gmx::Awh*                      awh,
                                 gmx_enfrot*                    enforcedRotation,
                                 gmx::ImdSession*               imdSession,
                                 pull_t*                        pull_work,
                                 int64_t                        step,
                                 double                         t,
                                 gmx_wallcycle*                 wcycle,
                                 gmx::ForceProviders*           forceProviders,
                                 const matrix                   box,
                                 gmx::ArrayRef<const gmx::RVec> x,
                                 const t_mdatoms*               mdatoms,
                                 gmx::ArrayRef<const real>      lambda,
                                 const StepWorkload&            stepWork,
                                 gmx::ForceWithVirial*          forceWithVirialMtsLevel0,
                                 gmx::ForceWithVirial*          forceWithVirialMtsLevel1,
                                 gmx_enerdata_t*                enerd,
                                 gmx_edsam*                     ed,
                                 bool                           didNeighborSearch)
{
    /* NOTE: Currently all ForceProviders only provide forces.
     *       When they also provide energies, remove this conditional.
     */
    if (stepWork.computeForces)
    {
        gmx::ForceProviderInput forceProviderInput(
                x,
                mdatoms->homenr,
                gmx::makeArrayRef(mdatoms->chargeA).subArray(0, mdatoms->homenr),
                gmx::makeArrayRef(mdatoms->massT).subArray(0, mdatoms->homenr),
                t,
                step,
                box,
                *cr);
        gmx::ForceProviderOutput forceProviderOutput(forceWithVirialMtsLevel0, enerd);

        /* Collect forces from modules */
        forceProviders->calculateForces(forceProviderInput, &forceProviderOutput);
    }

    const int  pullMtsLevel = forceGroupMtsLevel(inputrec.mtsLevels, gmx::MtsForceGroups::Pull);
    const bool doPulling    = (inputrec.bPull && pull_have_potential(*pull_work)
                            && (pullMtsLevel == 0 || stepWork.computeSlowForces));

    /* pull_potential_wrapper(), awh->applyBiasForcesAndUpdateBias(), pull_apply_forces()
     * have to be called in this order
     */
    if (doPulling)
    {
        pull_potential_wrapper(cr, inputrec, box, x, mdatoms, enerd, pull_work, lambda.data(), t, wcycle);
    }
    if (awh && (pullMtsLevel == 0 || stepWork.computeSlowForces))
    {
        const bool          needForeignEnergyDifferences = awh->needForeignEnergyDifferences(step);
        std::vector<double> foreignLambdaDeltaH, foreignLambdaDhDl;
        if (needForeignEnergyDifferences)
        {
            enerd->foreignLambdaTerms.finalizePotentialContributions(
                    enerd->dvdl_lin, lambda, *inputrec.fepvals);
            std::tie(foreignLambdaDeltaH, foreignLambdaDhDl) = enerd->foreignLambdaTerms.getTerms(cr);
        }

        enerd->term[F_COM_PULL] += awh->applyBiasForcesAndUpdateBias(
                inputrec.pbcType, foreignLambdaDeltaH, foreignLambdaDhDl, box, t, step, wcycle, fplog);
    }
    if (doPulling)
    {
        wallcycle_start_nocount(wcycle, WallCycleCounter::PullPot);
        auto& forceWithVirial = (pullMtsLevel == 0) ? forceWithVirialMtsLevel0 : forceWithVirialMtsLevel1;
        pull_apply_forces(pull_work, mdatoms->massT, cr, forceWithVirial);
        wallcycle_stop(wcycle, WallCycleCounter::PullPot);
    }

    /* Add the forces from enforced rotation potentials (if any) */
    if (inputrec.bRot)
    {
        wallcycle_start(wcycle, WallCycleCounter::RotAdd);
        enerd->term[F_COM_PULL] +=
                add_rot_forces(enforcedRotation, forceWithVirialMtsLevel0->force_, cr, step, t);
        wallcycle_stop(wcycle, WallCycleCounter::RotAdd);
    }

    if (ed)
    {
        /* Note that since init_edsam() is called after the initialization
         * of forcerec, edsam doesn't request the noVirSum force buffer.
         * Thus if no other algorithm (e.g. PME) requires it, the forces
         * here will contribute to the virial.
         */
        do_flood(cr, inputrec, x, forceWithVirialMtsLevel0->force_, ed, box, step, didNeighborSearch);
    }

    /* Add forces from interactive molecular dynamics (IMD), if any */
    if (inputrec.bIMD && stepWork.computeForces)
    {
        imdSession->applyForces(forceWithVirialMtsLevel0->force_);
    }
}

/*! \brief Launch the prepare_step and spread stages of PME GPU.
 *
 * \param[in]  pmedata              The PME structure
 * \param[in]  box                  The box matrix
 * \param[in]  stepWork             Step schedule flags
 * \param[in]  xReadyOnDevice       Event synchronizer indicating that the coordinates are ready in the device memory.
 * \param[in]  lambdaQ              The Coulomb lambda of the current state.
 * \param[in]  useMdGpuGraph        Whether MD GPU Graph is in use.
 * \param[in]  wcycle               The wallcycle structure
 */
static inline void launchPmeGpuSpread(gmx_pme_t*            pmedata,
                                      const matrix          box,
                                      const StepWorkload&   stepWork,
                                      GpuEventSynchronizer* xReadyOnDevice,
                                      const real            lambdaQ,
                                      bool                  useMdGpuGraph,
                                      gmx_wallcycle*        wcycle)
{
    wallcycle_start(wcycle, WallCycleCounter::PmeGpuMesh);
    pme_gpu_prepare_computation(pmedata, box, wcycle, stepWork);
    bool                           useGpuDirectComm         = false;
    gmx::PmeCoordinateReceiverGpu* pmeCoordinateReceiverGpu = nullptr;
    pme_gpu_launch_spread(
            pmedata, xReadyOnDevice, wcycle, lambdaQ, useGpuDirectComm, pmeCoordinateReceiverGpu, useMdGpuGraph);
    wallcycle_stop(wcycle, WallCycleCounter::PmeGpuMesh);
}

/*! \brief Launch the FFT and gather stages of PME GPU
 *
 * This function only implements setting the output forces (no accumulation).
 *
 * \param[in]  pmedata        The PME structure
 * \param[in]  lambdaQ        The Coulomb lambda of the current system state.
 * \param[in]  wcycle         The wallcycle structure
 * \param[in]  stepWork       Step schedule flags
 */
static void launchPmeGpuFftAndGather(gmx_pme_t*               pmedata,
                                     const real               lambdaQ,
                                     gmx_wallcycle*           wcycle,
                                     const gmx::StepWorkload& stepWork)
{
    wallcycle_start_nocount(wcycle, WallCycleCounter::PmeGpuMesh);
    pme_gpu_launch_complex_transforms(pmedata, wcycle, stepWork);
    pme_gpu_launch_gather(pmedata, wcycle, lambdaQ, stepWork.computeVirial);
    wallcycle_stop(wcycle, WallCycleCounter::PmeGpuMesh);
}

/*! \brief
 * Blocks until PME GPU tasks are completed, and gets the output forces and virial/energy
 * (if they were to be computed).
 *
 * \param[in]  pme             The PME data structure.
 * \param[in]  stepWork        The required work for this simulation step
 * \param[in]  wcycle          The wallclock counter.
 * \param[out] forceWithVirial The output force and virial
 * \param[out] enerd           The output energies
 * \param[in]  lambdaQ         The Coulomb lambda to use when calculating the results.
 */
static void pmeGpuWaitAndReduce(gmx_pme_t*               pme,
                                const gmx::StepWorkload& stepWork,
                                gmx_wallcycle*           wcycle,
                                gmx::ForceWithVirial*    forceWithVirial,
                                gmx_enerdata_t*          enerd,
                                const real               lambdaQ)
{
    wallcycle_start_nocount(wcycle, WallCycleCounter::PmeGpuMesh);

    pme_gpu_wait_and_reduce(pme, stepWork, wcycle, forceWithVirial, enerd, lambdaQ);

    wallcycle_stop(wcycle, WallCycleCounter::PmeGpuMesh);
}

/*! \brief
 *  Polling wait for either of the PME or nonbonded GPU tasks.
 *
 * Instead of a static order in waiting for GPU tasks, this function
 * polls checking which of the two tasks completes first, and does the
 * associated force buffer reduction overlapped with the other task.
 * By doing that, unlike static scheduling order, it can always overlap
 * one of the reductions, regardless of the GPU task completion order.
 *
 * \param[in]     nbv              Nonbonded verlet structure
 * \param[in,out] pmedata          PME module data
 * \param[in,out] forceOutputsNonbonded  Force outputs for the non-bonded forces and shift forces
 * \param[in,out] forceOutputsPme  Force outputs for the PME forces and virial
 * \param[in,out] enerd            Energy data structure results are reduced into
 * \param[in]     lambdaQ          The Coulomb lambda of the current system state.
 * \param[in]     stepWork         Step schedule flags
 * \param[in]     wcycle           The wallcycle structure
 */
static void alternatePmeNbGpuWaitReduce(nonbonded_verlet_t* nbv,
                                        gmx_pme_t*          pmedata,
                                        gmx::ForceOutputs*  forceOutputsNonbonded,
                                        gmx::ForceOutputs*  forceOutputsPme,
                                        gmx_enerdata_t*     enerd,
                                        const real          lambdaQ,
                                        const StepWorkload& stepWork,
                                        gmx_wallcycle*      wcycle)
{
    bool isPmeGpuDone = false;
    bool isNbGpuDone  = false;

    gmx::ArrayRef<const gmx::RVec> pmeGpuForces;

    while (!isPmeGpuDone || !isNbGpuDone)
    {
        if (!isPmeGpuDone)
        {
            wallcycle_start_nocount(wcycle, WallCycleCounter::PmeGpuMesh);
            GpuTaskCompletion completionType =
                    (isNbGpuDone) ? GpuTaskCompletion::Wait : GpuTaskCompletion::Check;
            isPmeGpuDone = pme_gpu_try_finish_task(
                    pmedata, stepWork, wcycle, &forceOutputsPme->forceWithVirial(), enerd, lambdaQ, completionType);
            wallcycle_stop(wcycle, WallCycleCounter::PmeGpuMesh);
        }

        if (!isNbGpuDone)
        {
            auto&             forceBuffersNonbonded = forceOutputsNonbonded->forceWithShiftForces();
            GpuTaskCompletion completionType =
                    (isPmeGpuDone) ? GpuTaskCompletion::Wait : GpuTaskCompletion::Check;
            // To get the wcycle call count right, when in GpuTaskCompletion::Check mode,
            // we start without counting and only when the task finished we issue a
            // start/stop to increment.
            // GpuTaskCompletion::Wait mode the timing is expected to be done in the caller.
            wallcycle_start_nocount(wcycle, WallCycleCounter::WaitGpuNbL);
            isNbGpuDone = Nbnxm::gpu_try_finish_task(
                    nbv->gpuNbv(),
                    stepWork,
                    AtomLocality::Local,
                    enerd->grpp.energyGroupPairTerms[NonBondedEnergyTerms::LJSR].data(),
                    enerd->grpp.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR].data(),
                    forceBuffersNonbonded.shiftForces(),
                    completionType);
            wallcycle_stop(wcycle, WallCycleCounter::WaitGpuNbL);

            if (isNbGpuDone)
            {
                wallcycle_increment_event_count(wcycle, WallCycleCounter::WaitGpuNbL);
                nbv->atomdata_add_nbat_f_to_f(AtomLocality::Local, forceBuffersNonbonded.force());
            }
        }
    }
}

/*! \brief Set up the different force buffers; also does clearing.
 *
 * \param[in] forceHelperBuffers        Helper force buffers
 * \param[in] force                     force array
 * \param[in] domainWork                Domain lifetime workload flags
 * \param[in] stepWork                  Step schedule flags
 * \param[in] havePpDomainDecomposition Whether we have a PP domain decomposition
 * \param[out] wcycle                   wallcycle recording structure
 *
 * \returns                             Cleared force output structure
 */
static ForceOutputs setupForceOutputs(ForceHelperBuffers*                 forceHelperBuffers,
                                      gmx::ArrayRefWithPadding<gmx::RVec> force,
                                      const DomainLifetimeWorkload&       domainWork,
                                      const StepWorkload&                 stepWork,
                                      const bool                          havePpDomainDecomposition,
                                      gmx_wallcycle*                      wcycle)
{
    /* NOTE: We assume fr->shiftForces is all zeros here */
    gmx::ForceWithShiftForces forceWithShiftForces(
            force, stepWork.computeVirial, forceHelperBuffers->shiftForces());

    if (stepWork.computeForces
        && (domainWork.haveCpuLocalForceWork || !stepWork.useGpuFBufferOps
            || (havePpDomainDecomposition && !stepWork.useGpuFHalo)))
    {
        wallcycle_sub_start(wcycle, WallCycleSubCounter::ClearForceBuffer);
        /* Clear the short- and long-range forces */
        clearRVecs(forceWithShiftForces.force(), true);

        /* Clear the shift forces */
        clearRVecs(forceWithShiftForces.shiftForces(), false);
        wallcycle_sub_stop(wcycle, WallCycleSubCounter::ClearForceBuffer);
    }

    /* If we need to compute the virial, we might need a separate
     * force buffer for algorithms for which the virial is calculated
     * directly, such as PME. Otherwise, forceWithVirial uses the
     * the same force (f in legacy calls) buffer as other algorithms.
     */
    const bool useSeparateForceWithVirialBuffer =
            (stepWork.computeForces
             && (stepWork.computeVirial && forceHelperBuffers->haveDirectVirialContributions()));
    /* forceWithVirial uses the local atom range only */
    gmx::ForceWithVirial forceWithVirial(
            useSeparateForceWithVirialBuffer ? forceHelperBuffers->forceBufferForDirectVirialContributions()
                                             : force.unpaddedArrayRef(),
            stepWork.computeVirial);

    if (useSeparateForceWithVirialBuffer)
    {
        wallcycle_sub_start_nocount(wcycle, WallCycleSubCounter::ClearForceBuffer);
        /* TODO: update comment
         * We only compute forces on local atoms. Note that vsites can
         * spread to non-local atoms, but that part of the buffer is
         * cleared separately in the vsite spreading code.
         */
        clearRVecs(forceWithVirial.force_, true);
        wallcycle_sub_stop(wcycle, WallCycleSubCounter::ClearForceBuffer);
    }


    return ForceOutputs(
            forceWithShiftForces, forceHelperBuffers->haveDirectVirialContributions(), forceWithVirial);
}

/* \brief Launch end-of-step GPU tasks: buffer clearing and rolling pruning.
 *
 */
static void launchGpuEndOfStepTasks(nonbonded_verlet_t*               nbv,
                                    gmx::ListedForcesGpu*             listedForcesGpu,
                                    gmx_pme_t*                        pmedata,
                                    gmx_enerdata_t*                   enerd,
                                    const gmx::MdrunScheduleWorkload& runScheduleWork,
                                    int64_t                           step,
                                    gmx_wallcycle*                    wcycle)
{
    if (runScheduleWork.simulationWork.useGpuNonbonded && runScheduleWork.stepWork.computeNonbondedForces)
    {
        /* Launch pruning before buffer clearing because the API overhead of the
         * clear kernel launches can leave the GPU idle while it could be running
         * the prune kernel.
         */
        if (nbv->isDynamicPruningStepGpu(step))
        {
            nbv->dispatchPruneKernelGpu(step);
        }

        /* now clear the GPU outputs while we finish the step on the CPU */
        wallcycle_start_nocount(wcycle, WallCycleCounter::LaunchGpuPp);
        wallcycle_sub_start_nocount(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);
        Nbnxm::gpu_clear_outputs(nbv->gpuNbv(), runScheduleWork.stepWork.computeVirial);
        wallcycle_sub_stop(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);
        wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPp);
    }

    if (runScheduleWork.stepWork.haveGpuPmeOnThisRank)
    {
        wallcycle_start_nocount(wcycle, WallCycleCounter::PmeGpuMesh);
        bool gpuGraphWithSeparatePmeRank = false;
        pme_gpu_reinit_computation(pmedata, gpuGraphWithSeparatePmeRank, wcycle);
        wallcycle_stop(wcycle, WallCycleCounter::PmeGpuMesh);
    }

    if (runScheduleWork.domainWork.haveGpuBondedWork && runScheduleWork.stepWork.computeEnergy)
    {
        // in principle this should be included in the DD balancing region,
        // but generally it is infrequent so we'll omit it for the sake of
        // simpler code
        listedForcesGpu->waitAccumulateEnergyTerms(enerd);

        listedForcesGpu->clearEnergies();
    }
}

/*! \brief Compute the number of times the "local coordinates ready on device" GPU event will be used as a synchronization point.
 *
 * When some work is offloaded to GPU, force calculation should wait for the atom coordinates to
 * be ready on the device. The coordinates can come either from H2D copy at the beginning of the step,
 * or from the GPU integration at the end of the previous step.
 *
 * In GROMACS, we usually follow the "mark once - wait once" approach. But this event is "consumed"
 * (that is, waited upon either on host or on the device) multiple times, since many tasks
 * in different streams depend on the coordinates.
 *
 * This function return the number of times the event will be consumed based on this step's workload.
 *
 * \param simulationWork Simulation workload flags.
 * \param stepWork Step workload flags.
 * \param pmeSendCoordinatesFromGpu Whether peer-to-peer communication is used for PME coordinates.
 * \return
 */
static int getExpectedLocalXReadyOnDeviceConsumptionCount(gmx_used_in_debug const SimulationWorkload& simulationWork,
                                                          const StepWorkload& stepWork,
                                                          bool pmeSendCoordinatesFromGpu)
{
    int result = 0;
    if (stepWork.computeSlowForces)
    {
        if (pmeSendCoordinatesFromGpu)
        {
            GMX_ASSERT(simulationWork.haveSeparatePmeRank,
                       "GPU PME PP communications require having a separate PME rank");
            // Event is consumed by gmx_pme_send_coordinates for GPU PME PP Communications
            result++;
        }
        if (stepWork.haveGpuPmeOnThisRank)
        {
            // Event is consumed by launchPmeGpuSpread
            result++;
        }
        if (stepWork.computeNonbondedForces && stepWork.useGpuXBufferOps)
        {
            // Event is consumed by convertCoordinatesGpu
            result++;
        }
    }
    if (stepWork.useGpuXHalo)
    {
        // Event is consumed by communicateGpuHaloCoordinates
        result++;
        if (GMX_THREAD_MPI) // Issue #4262
        {
            result++;
        }
    }
    if (stepWork.clearGpuFBufferEarly && simulationWork.useGpuUpdate)
    {
        // Event is consumed by force clearing which waits for the update to complete
        result++;
    }
    return result;
}

/*! \brief Compute the number of times the "local forces ready on device" GPU event will be used as a synchronization point.
 *
 * In GROMACS, we usually follow the "mark once - wait once" approach. But this event is "consumed"
 * (that is, waited upon either on host or on the device) multiple times, since many tasks
 * in different streams depend on the local forces.
 *
 * \param simulationWork Simulation workload flags.
 * \param domainWork Domain workload flags.
 * \param stepWork Step workload flags.
 * \param useOrEmulateGpuNb Whether GPU non-bonded calculations are used or emulated.
 * \param alternateGpuWait Whether alternating wait/reduce scheme is used.
 * \return The number of times the event will be consumed based on this step's workload.
 */
static int getExpectedLocalFReadyOnDeviceConsumptionCount(const SimulationWorkload& simulationWork,
                                                          const DomainLifetimeWorkload& domainWork,
                                                          const StepWorkload&           stepWork,
                                                          bool useOrEmulateGpuNb,
                                                          bool alternateGpuWait)
{
    int  counter = 0;
    bool eventUsedInGpuForceReduction =
            (domainWork.haveCpuLocalForceWork
             || (simulationWork.havePpDomainDecomposition && !simulationWork.useGpuHaloExchange));
    bool gpuForceReductionUsed = useOrEmulateGpuNb && !alternateGpuWait && stepWork.useGpuFBufferOps
                                 && stepWork.computeNonbondedForces;
    if (gpuForceReductionUsed && eventUsedInGpuForceReduction)
    {
        counter++;
    }
    bool gpuForceHaloUsed = simulationWork.havePpDomainDecomposition && stepWork.computeForces
                            && stepWork.useGpuFHalo;
    if (gpuForceHaloUsed)
    {
        counter++;
    }
    return counter;
}

//! \brief Data structure to hold dipole-related data and staging arrays
struct DipoleData
{
    //! Dipole staging for fast summing over MPI
    gmx::DVec muStaging[2] = { { 0.0, 0.0, 0.0 } };
    //! Dipole staging for states A and B (index 0 and 1 resp.)
    gmx::RVec muStateAB[2] = { { 0.0_real, 0.0_real, 0.0_real } };
};


static void reduceAndUpdateMuTot(DipoleData*                   dipoleData,
                                 const t_commrec*              cr,
                                 const bool                    haveFreeEnergy,
                                 gmx::ArrayRef<const real>     lambda,
                                 rvec                          muTotal,
                                 const DDBalanceRegionHandler& ddBalanceRegionHandler)
{
    if (PAR(cr))
    {
        gmx_sumd(2 * DIM, dipoleData->muStaging[0], cr);
        ddBalanceRegionHandler.reopenRegionCpu();
    }
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            dipoleData->muStateAB[i][j] = dipoleData->muStaging[i][j];
        }
    }

    if (!haveFreeEnergy)
    {
        copy_rvec(dipoleData->muStateAB[0], muTotal);
    }
    else
    {
        for (int j = 0; j < DIM; j++)
        {
            muTotal[j] = (1.0 - lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)])
                                 * dipoleData->muStateAB[0][j]
                         + lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)]
                                   * dipoleData->muStateAB[1][j];
        }
    }
}

/*! \brief Combines MTS level0 and level1 force buffers into a full and MTS-combined force buffer.
 *
 * \param[in]     numAtoms        The number of atoms to combine forces for
 * \param[in,out] forceMtsLevel0  Input: F_level0, output: F_level0 + F_level1
 * \param[in,out] forceMts        Input: F_level1, output: F_level0 + mtsFactor * F_level1
 * \param[in]     mtsFactor       The factor between the level0 and level1 time step
 */
static void combineMtsForces(const int      numAtoms,
                             ArrayRef<RVec> forceMtsLevel0,
                             ArrayRef<RVec> forceMts,
                             const real     mtsFactor)
{
    const int gmx_unused numThreads = gmx_omp_nthreads_get(ModuleMultiThread::Default);
#pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int i = 0; i < numAtoms; i++)
    {
        const RVec forceMtsLevel0Tmp = forceMtsLevel0[i];
        forceMtsLevel0[i] += forceMts[i];
        forceMts[i] = forceMtsLevel0Tmp + mtsFactor * forceMts[i];
    }
}

/*! \brief Setup for the local GPU force reduction:
 * reinitialization plus the registration of forces and dependencies.
 *
 * \param [in] runScheduleWork     Schedule workload flag structure
 * \param [in] nbv                 Non-bonded Verlet object
 * \param [in] stateGpu            GPU state propagator object
 * \param [in] gpuForceReduction   GPU force reduction object
 * \param [in] pmePpCommGpu        PME-PP GPU communication object
 * \param [in] pmedata             PME data object
 * \param [in] dd                  Domain decomposition object
 */
static void setupLocalGpuForceReduction(const gmx::MdrunScheduleWorkload& runScheduleWork,
                                        nonbonded_verlet_t*               nbv,
                                        gmx::StatePropagatorDataGpu*      stateGpu,
                                        gmx::GpuForceReduction*           gpuForceReduction,
                                        gmx::PmePpCommGpu*                pmePpCommGpu,
                                        const gmx_pme_t*                  pmedata,
                                        const gmx_domdec_t*               dd)
{
    GMX_ASSERT(!runScheduleWork.simulationWork.useMts,
               "GPU force reduction is not compatible with MTS");

    // (re-)initialize local GPU force reduction
    const bool accumulate = runScheduleWork.domainWork.haveCpuLocalForceWork
                            || runScheduleWork.simulationWork.havePpDomainDecomposition;
    const int atomStart = 0;
    gpuForceReduction->reinit(stateGpu->getForces(),
                              nbv->getNumAtoms(AtomLocality::Local),
                              nbv->getGridIndices(),
                              atomStart,
                              accumulate,
                              stateGpu->fReducedOnDevice(AtomLocality::Local));

    // register forces and add dependencies
    gpuForceReduction->registerNbnxmForce(Nbnxm::gpu_get_f(nbv->gpuNbv()));

    DeviceBuffer<gmx::RVec> pmeForcePtr;
    GpuEventSynchronizer*   pmeSynchronizer     = nullptr;
    bool                    havePmeContribution = false;

    if (runScheduleWork.simulationWork.haveGpuPmeOnPpRank())
    {
        pmeForcePtr = pme_gpu_get_device_f(pmedata);
        if (pmeForcePtr)
        {
            pmeSynchronizer     = pme_gpu_get_f_ready_synchronizer(pmedata);
            havePmeContribution = true;
        }
    }
    else if (runScheduleWork.simulationWork.useGpuPmePpCommunication)
    {
        pmeForcePtr = pmePpCommGpu->getGpuForceStagingPtr();
        GMX_ASSERT(pmeForcePtr, "PME force for reduction has no data");
        if (GMX_THREAD_MPI)
        {
            pmeSynchronizer = pmePpCommGpu->getForcesReadySynchronizer();
        }
        havePmeContribution = true;
    }

    if (havePmeContribution)
    {
        gpuForceReduction->registerRvecForce(pmeForcePtr);
        if (runScheduleWork.simulationWork.useNvshmem)
        {
            DeviceBuffer<uint64_t> forcesReadyNvshmemFlags = pmePpCommGpu->getGpuForcesSyncObj();
            gpuForceReduction->registerForcesReadyNvshmemFlags(forcesReadyNvshmemFlags);
        }

        if (!runScheduleWork.simulationWork.useGpuPmePpCommunication || GMX_THREAD_MPI)
        {
            GMX_ASSERT(pmeSynchronizer != nullptr, "PME force ready cuda event should not be NULL");
            gpuForceReduction->addDependency(pmeSynchronizer);
        }
    }

    if (runScheduleWork.domainWork.haveCpuLocalForceWork
        || (runScheduleWork.simulationWork.havePpDomainDecomposition
            && !runScheduleWork.simulationWork.useGpuHaloExchange))
    {
        gpuForceReduction->addDependency(stateGpu->fReadyOnDevice(AtomLocality::Local));
    }

    if (runScheduleWork.simulationWork.useGpuHaloExchange)
    {
        gpuForceReduction->addDependency(dd->gpuHaloExchange[0][0]->getForcesReadyOnDeviceEvent());
    }
}

/*! \brief Setup for the non-local GPU force reduction:
 * reinitialization plus the registration of forces and dependencies.
 *
 * \param [in] runScheduleWork     Schedule workload flag structure
 * \param [in] nbv                 Non-bonded Verlet object
 * \param [in] stateGpu            GPU state propagator object
 * \param [in] gpuForceReduction   GPU force reduction object
 * \param [in] dd                  Domain decomposition object
 */
static void setupNonLocalGpuForceReduction(const gmx::MdrunScheduleWorkload& runScheduleWork,
                                           nonbonded_verlet_t*               nbv,
                                           gmx::StatePropagatorDataGpu*      stateGpu,
                                           gmx::GpuForceReduction*           gpuForceReduction,
                                           const gmx_domdec_t*               dd)
{
    // (re-)initialize non-local GPU force reduction
    const bool accumulate = runScheduleWork.domainWork.haveCpuNonLocalForceWork;
    const int  atomStart  = dd_numHomeAtoms(*dd);
    gpuForceReduction->reinit(stateGpu->getForces(),
                              nbv->getNumAtoms(AtomLocality::NonLocal),
                              nbv->getGridIndices(),
                              atomStart,
                              accumulate,
                              stateGpu->fReducedOnDevice(AtomLocality::NonLocal));

    // register forces and add dependencies
    gpuForceReduction->registerNbnxmForce(Nbnxm::gpu_get_f(nbv->gpuNbv()));

    if (runScheduleWork.domainWork.haveCpuNonLocalForceWork)
    {
        gpuForceReduction->addDependency(stateGpu->fReadyOnDevice(AtomLocality::NonLocal));
    }
}


/*! \brief Return the number of local atoms.
 */
static int getLocalAtomCount(const gmx_domdec_t* dd, const t_mdatoms& mdatoms, bool havePPDomainDecomposition)
{
    GMX_ASSERT(!(havePPDomainDecomposition && (dd == nullptr)),
               "Can't have PP decomposition with dd uninitialized!");
    return havePPDomainDecomposition ? dd_numAtomsZones(*dd) : mdatoms.homenr;
}

/*! \brief Does pair search and closely related activities required on search steps.
 */
static void doPairSearch(const t_commrec*                    cr,
                         const t_inputrec&                   inputrec,
                         const gmx::MDModulesNotifiers&      mdModulesNotifiers,
                         int64_t                             step,
                         t_nrnb*                             nrnb,
                         gmx_wallcycle*                      wcycle,
                         const gmx_localtop_t&               top,
                         const matrix                        box,
                         gmx::ArrayRefWithPadding<gmx::RVec> x,
                         gmx::ArrayRef<gmx::RVec>            v,
                         const t_mdatoms&                    mdatoms,
                         t_forcerec*                         fr,
                         const gmx::MdrunScheduleWorkload&   runScheduleWork)
{
    nonbonded_verlet_t* nbv = fr->nbv.get();

    gmx::StatePropagatorDataGpu* stateGpu = fr->stateGpu;

    const SimulationWorkload& simulationWork = runScheduleWork.simulationWork;
    const StepWorkload&       stepWork       = runScheduleWork.stepWork;

    if (gmx::needStateGpu(simulationWork))
    {
        // TODO refactor this to do_md, after partitioning.
        stateGpu->reinit(mdatoms.homenr,
                         getLocalAtomCount(cr->dd, mdatoms, simulationWork.havePpDomainDecomposition));
    }

    if (simulationWork.haveGpuPmeOnPpRank())
    {
        GMX_ASSERT(gmx::needStateGpu(simulationWork), "StatePropagatorDataGpu is needed");
        // TODO: This should be moved into PME setup function ( pme_gpu_prepare_computation(...) )
        pme_gpu_set_device_x(fr->pmedata, stateGpu->getCoordinates());
    }

    if (fr->pbcType != PbcType::No)
    {
        const bool calcCGCM = (stepWork.stateChanged && !haveDDAtomOrdering(*cr));
        if (calcCGCM)
        {
            put_atoms_in_box_omp(fr->pbcType,
                                 box,
                                 fr->haveBoxDeformation,
                                 inputrec.deform,
                                 x.unpaddedArrayRef().subArray(0, mdatoms.homenr),
                                 v.empty() ? ArrayRef<RVec>() : v.subArray(0, mdatoms.homenr),
                                 gmx_omp_nthreads_get(ModuleMultiThread::Default));
            inc_nrnb(nrnb, eNR_SHIFTX, mdatoms.homenr);
        }

        if (!haveDDAtomOrdering(*cr))
        {
            // Atoms might have changed periodic image, signal MDModules
            gmx::MDModulesAtomsRedistributedSignal mdModulesAtomsRedistributedSignal(
                    box, x.unpaddedArrayRef().subArray(0, mdatoms.homenr));
            mdModulesNotifiers.simulationSetupNotifier_.notify(mdModulesAtomsRedistributedSignal);
        }
    }

    if (fr->wholeMoleculeTransform && stepWork.stateChanged)
    {
        fr->wholeMoleculeTransform->updateForAtomPbcJumps(x.unpaddedArrayRef(), box);
    }

    wallcycle_start(wcycle, WallCycleCounter::NS);
    if (!haveDDAtomOrdering(*cr))
    {
        const rvec vzero       = { 0.0_real, 0.0_real, 0.0_real };
        const rvec boxDiagonal = { box[XX][XX], box[YY][YY], box[ZZ][ZZ] };
        wallcycle_sub_start(wcycle, WallCycleSubCounter::NBSGridLocal);
        nbv->putAtomsOnGrid(
                box, 0, vzero, boxDiagonal, nullptr, { 0, mdatoms.homenr }, -1, fr->atomInfo, x.unpaddedArrayRef(), 0, nullptr);
        wallcycle_sub_stop(wcycle, WallCycleSubCounter::NBSGridLocal);
    }
    else
    {
        wallcycle_sub_start(wcycle, WallCycleSubCounter::NBSGridNonLocal);
        nbnxn_put_on_grid_nonlocal(nbv, domdec_zones(cr->dd), fr->atomInfo, x.unpaddedArrayRef());
        wallcycle_sub_stop(wcycle, WallCycleSubCounter::NBSGridNonLocal);
    }

    nbv->setAtomProperties(mdatoms.typeA, mdatoms.chargeA, fr->atomInfo);

    wallcycle_stop(wcycle, WallCycleCounter::NS);

    /* initialize the GPU nbnxm atom data and bonded data structures */
    if (simulationWork.useGpuNonbonded)
    {
        // Note: cycle counting only nononbondeds, GPU listed forces counts internally
        wallcycle_start_nocount(wcycle, WallCycleCounter::LaunchGpuPp);
        wallcycle_sub_start_nocount(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);
        Nbnxm::gpu_init_atomdata(nbv->gpuNbv(), &nbv->nbat());
        wallcycle_sub_stop(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);
        wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPp);

        if (fr->listedForcesGpu)
        {
            /* Now we put all atoms on the grid, we can assign bonded
             * interactions to the GPU, where the grid order is
             * needed. Also the xq, f and fshift device buffers have
             * been reallocated if needed, so the bonded code can
             * learn about them. */
            // TODO the xq, f, and fshift buffers are now shared
            // resources, so they should be maintained by a
            // higher-level object than the nb module.
            fr->listedForcesGpu->updateInteractionListsAndDeviceBuffers(
                    nbv->getGridIndices(), top.idef, Nbnxm::gpuGetNBAtomData(nbv->gpuNbv()));
        }
    }

    wallcycle_start_nocount(wcycle, WallCycleCounter::NS);
    wallcycle_sub_start(wcycle, WallCycleSubCounter::NBSSearchLocal);
    /* Note that with a GPU the launch overhead of the list transfer is not timed separately */
    nbv->constructPairlist(InteractionLocality::Local, top.excls, step, nrnb);

    nbv->setupGpuShortRangeWork(fr->listedForcesGpu.get(), InteractionLocality::Local);

    wallcycle_sub_stop(wcycle, WallCycleSubCounter::NBSSearchLocal);
    wallcycle_stop(wcycle, WallCycleCounter::NS);

    if (simulationWork.useGpuXBufferOpsWhenAllowed)
    {
        nbv->atomdata_init_copy_x_to_nbat_x_gpu();
    }

    if (simulationWork.useGpuFBufferOpsWhenAllowed)
    {
        // with MPI, direct GPU communication, and separate PME ranks we need
        // gmx_pme_send_coordinates() to be called before we can set up force reduction
        bool delaySetupLocalGpuForceReduction = GMX_MPI && simulationWork.useGpuPmePpCommunication;
        if (!delaySetupLocalGpuForceReduction)
        {
            setupLocalGpuForceReduction(runScheduleWork,
                                        nbv,
                                        stateGpu,
                                        fr->gpuForceReduction[gmx::AtomLocality::Local].get(),
                                        fr->pmePpCommGpu.get(),
                                        fr->pmedata,
                                        cr->dd);
        }

        if (simulationWork.havePpDomainDecomposition)
        {
            setupNonLocalGpuForceReduction(runScheduleWork,
                                           nbv,
                                           stateGpu,
                                           fr->gpuForceReduction[gmx::AtomLocality::NonLocal].get(),
                                           cr->dd);
        }
    }

    /* do non-local pair search */
    if (simulationWork.havePpDomainDecomposition)
    {
        wallcycle_start_nocount(wcycle, WallCycleCounter::NS);
        wallcycle_sub_start(wcycle, WallCycleSubCounter::NBSSearchNonLocal);
        /* Note that with a GPU the launch overhead of the list transfer is not timed separately */
        nbv->constructPairlist(InteractionLocality::NonLocal, top.excls, step, nrnb);

        nbv->setupGpuShortRangeWork(fr->listedForcesGpu.get(), InteractionLocality::NonLocal);
        wallcycle_sub_stop(wcycle, WallCycleSubCounter::NBSSearchNonLocal);
        wallcycle_stop(wcycle, WallCycleCounter::NS);
        // TODO refactor this GPU halo exchange re-initialisation
        // to location in do_md where GPU halo exchange is
        // constructed at partitioning, after above stateGpu
        // re-initialization has similarly been refactored
        if (simulationWork.useGpuHaloExchange)
        {
            reinitGpuHaloExchange(*cr, stateGpu->getCoordinates(), stateGpu->getForces());
        }
    }

    // With FEP we set up the reduction over threads for local+non-local simultaneously,
    // so we need to do that here after the local and non-local pairlist construction.
    if (fr->efep != FreeEnergyPerturbationType::No)
    {
        wallcycle_sub_start(wcycle, WallCycleSubCounter::NonbondedFep);
        nbv->setupFepThreadedForceBuffer(fr->natoms_force_constr);
        wallcycle_sub_stop(wcycle, WallCycleSubCounter::NonbondedFep);
    }
}

void do_force(FILE*                               fplog,
              const t_commrec*                    cr,
              const gmx_multisim_t*               ms,
              const t_inputrec&                   inputrec,
              const gmx::MDModulesNotifiers&      mdModulesNotifiers,
              gmx::Awh*                           awh,
              gmx_enfrot*                         enforcedRotation,
              gmx::ImdSession*                    imdSession,
              pull_t*                             pull_work,
              int64_t                             step,
              t_nrnb*                             nrnb,
              gmx_wallcycle*                      wcycle,
              const gmx_localtop_t*               top,
              const matrix                        box,
              gmx::ArrayRefWithPadding<gmx::RVec> x,
              gmx::ArrayRef<gmx::RVec>            v,
              const history_t*                    hist,
              gmx::ForceBuffersView*              forceView,
              tensor                              vir_force,
              const t_mdatoms*                    mdatoms,
              gmx_enerdata_t*                     enerd,
              gmx::ArrayRef<const real>           lambda,
              t_forcerec*                         fr,
              const gmx::MdrunScheduleWorkload&   runScheduleWork,
              gmx::VirtualSitesHandler*           vsite,
              rvec                                muTotal,
              double                              t,
              gmx_edsam*                          ed,
              CpuPpLongRangeNonbondeds*           longRangeNonbondeds,
              const DDBalanceRegionHandler&       ddBalanceRegionHandler)
{
    auto force = forceView->forceWithPadding();
    GMX_ASSERT(force.unpaddedArrayRef().ssize() >= fr->natoms_force_constr,
               "The size of the force buffer should be at least the number of atoms to compute "
               "forces for");

    nonbonded_verlet_t*  nbv = fr->nbv.get();
    interaction_const_t* ic  = fr->ic.get();

    gmx::StatePropagatorDataGpu* stateGpu = fr->stateGpu;

    const SimulationWorkload& simulationWork = runScheduleWork.simulationWork;

    const gmx::DomainLifetimeWorkload& domainWork = runScheduleWork.domainWork;

    const StepWorkload& stepWork = runScheduleWork.stepWork;

    if (stepWork.doNeighborSearch)
    {
        doPairSearch(cr, inputrec, mdModulesNotifiers, step, nrnb, wcycle, *top, box, x, v, *mdatoms, fr, runScheduleWork);

        /* At a search step we need to start the first balancing region
         * somewhere early inside the step after communication during domain
         * decomposition (and not during the previous step as usual).
         */
        ddBalanceRegionHandler.openBeforeForceComputationCpu(DdAllowBalanceRegionReopen::yes);
    }

    const bool pmeSendCoordinatesFromGpu =
            simulationWork.useGpuPmePpCommunication && !stepWork.doNeighborSearch;
    auto* localXReadyOnDevice = (stepWork.haveGpuPmeOnThisRank || stepWork.useGpuXBufferOps
                                 || simulationWork.useGpuUpdate || pmeSendCoordinatesFromGpu)
                                        ? stateGpu->getCoordinatesReadyOnDeviceEvent(
                                                AtomLocality::Local, simulationWork, stepWork)
                                        : nullptr;

    if (stepWork.clearGpuFBufferEarly)
    {
        // GPU Force halo exchange will set a subset of local atoms with remote non-local data.
        // First clear local portion of force array, so that untouched atoms are zero.
        // The dependency for this is that forces from previous timestep have been consumed,
        // which is satisfied when localXReadyOnDevice has been marked for GPU update case.
        // For CPU update, the forces are consumed by the beginning of the step, so no extra sync needed.
        GpuEventSynchronizer* dependency = simulationWork.useGpuUpdate ? localXReadyOnDevice : nullptr;
        stateGpu->clearForcesOnGpu(AtomLocality::Local, dependency);
    }

    clear_mat(vir_force);

    if (fr->pbcType != PbcType::No)
    {
        /* Compute shift vectors every step,
         * because of pressure coupling or box deformation!
         */
        if (stepWork.haveDynamicBox && stepWork.stateChanged)
        {
            calc_shifts(box, fr->shift_vec);
        }
    }
    nbnxn_atomdata_copy_shiftvec(stepWork.haveDynamicBox, fr->shift_vec, &nbv->nbat());


    GMX_ASSERT(simulationWork.useGpuHaloExchange
                       == ((cr->dd != nullptr) && (!cr->dd->gpuHaloExchange[0].empty())),
               "The GPU halo exchange is active, but it has not been constructed.");

    bool gmx_used_in_debug haveCopiedXFromGpu = false;
    // Copy coordinate from the GPU if update is on the GPU and there
    // are forces to be computed on the CPU, or for the computation of
    // virial, or if host-side data will be transferred from this task
    // to a remote task for halo exchange or PME-PP communication. At
    // search steps the current coordinates are already on the host,
    // hence copy is not needed.
    if (simulationWork.useGpuUpdate && !stepWork.doNeighborSearch
        && (runScheduleWork.domainWork.haveCpuLocalForceWork || stepWork.computeVirial
            || simulationWork.useCpuPmePpCommunication || simulationWork.useCpuHaloExchange
            || simulationWork.computeMuTot))
    {
        stateGpu->copyCoordinatesFromGpu(x.unpaddedArrayRef(), AtomLocality::Local);
        haveCopiedXFromGpu = true;
    }

    // Coordinates on the device are needed if PME or BufferOps are offloaded.
    // The local coordinates can be copied right away.
    // NOTE: Consider moving this copy to right after they are updated and constrained,
    //       if the later is not offloaded.
    if (stepWork.haveGpuPmeOnThisRank || stepWork.useGpuXBufferOps || pmeSendCoordinatesFromGpu)
    {
        GMX_ASSERT(stateGpu != nullptr, "stateGpu should not be null");
        const int expectedLocalXReadyOnDeviceConsumptionCount =
                getExpectedLocalXReadyOnDeviceConsumptionCount(
                        simulationWork, stepWork, pmeSendCoordinatesFromGpu);

        // We need to copy coordinates when:
        // 1. Update is not offloaded
        // 2. The buffers were reinitialized on search step
        if (!simulationWork.useGpuUpdate || stepWork.doNeighborSearch)
        {
            stateGpu->copyCoordinatesToGpu(x.unpaddedArrayRef(),
                                           AtomLocality::Local,
                                           expectedLocalXReadyOnDeviceConsumptionCount);
        }
        else if (simulationWork.useGpuUpdate)
        {
            stateGpu->setXUpdatedOnDeviceEventExpectedConsumptionCount(
                    expectedLocalXReadyOnDeviceConsumptionCount);
        }
    }

    if (stepWork.computePmeOnSeparateRank)
    {
        /* Send particle coordinates to the pme nodes */
        if (!pmeSendCoordinatesFromGpu && !stepWork.doNeighborSearch && simulationWork.useGpuUpdate)
        {
            GMX_ASSERT(haveCopiedXFromGpu,
                       "a wait should only be triggered if copy has been scheduled");
            stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
        }

        const bool reinitGpuPmePpComms =
                simulationWork.useGpuPmePpCommunication && stepWork.doNeighborSearch;
        gmx_pme_send_coordinates(fr,
                                 cr,
                                 box,
                                 x.unpaddedArrayRef(),
                                 lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)],
                                 lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)],
                                 (stepWork.computeVirial || stepWork.computeEnergy),
                                 step,
                                 simulationWork.useGpuPmePpCommunication,
                                 reinitGpuPmePpComms,
                                 pmeSendCoordinatesFromGpu,
                                 stepWork.useGpuPmeFReduction,
                                 pmeSendCoordinatesFromGpu ? localXReadyOnDevice : nullptr,
                                 simulationWork.useMdGpuGraph,
                                 wcycle);
    }

    if (simulationWork.useGpuFBufferOpsWhenAllowed && stepWork.doNeighborSearch)
    {
        // with MPI, direct GPU communication, and separate PME ranks we need
        // gmx_pme_send_coordinates() to be called before we can set up force reduction
        bool doSetupLocalGpuForceReduction = GMX_MPI && simulationWork.useGpuPmePpCommunication;
        if (doSetupLocalGpuForceReduction)
        {
            setupLocalGpuForceReduction(runScheduleWork,
                                        fr->nbv.get(),
                                        stateGpu,
                                        fr->gpuForceReduction[gmx::AtomLocality::Local].get(),
                                        fr->pmePpCommGpu.get(),
                                        fr->pmedata,
                                        cr->dd);
        }
    }

    if (stepWork.haveGpuPmeOnThisRank)
    {
        launchPmeGpuSpread(fr->pmedata,
                           box,
                           stepWork,
                           localXReadyOnDevice,
                           lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)],
                           simulationWork.useMdGpuGraph,
                           wcycle);
    }


    if (!stepWork.doNeighborSearch && !EI_TPI(inputrec.eI) && stepWork.computeNonbondedForces)
    {
        if (stepWork.useGpuXBufferOps)
        {
            GMX_ASSERT(stateGpu, "stateGpu should be valid when buffer ops are offloaded");
            nbv->convertCoordinatesGpu(AtomLocality::Local, stateGpu->getCoordinates(), localXReadyOnDevice);
        }
        else
        {
            if (simulationWork.useGpuUpdate)
            {
                GMX_ASSERT(stateGpu, "need a valid stateGpu object");
                GMX_ASSERT(haveCopiedXFromGpu,
                           "a wait should only be triggered if copy has been scheduled");
                stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
            }
            nbv->convertCoordinates(AtomLocality::Local, x.unpaddedArrayRef());
        }
    }

    if (simulationWork.useGpuNonbonded && (stepWork.computeNonbondedForces || domainWork.haveGpuBondedWork))
    {
        ddBalanceRegionHandler.openBeforeForceComputationGpu();

        wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPp);
        wallcycle_sub_start(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);
        Nbnxm::gpu_upload_shiftvec(nbv->gpuNbv(), &nbv->nbat());
        if (!stepWork.useGpuXBufferOps)
        {
            Nbnxm::gpu_copy_xq_to_gpu(nbv->gpuNbv(), &nbv->nbat(), AtomLocality::Local);
        }
        wallcycle_sub_stop(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);
        wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPp);
        // with X buffer ops offloaded to the GPU on all but the search steps

        // bonded work not split into separate local and non-local, so with DD
        // we can only launch the kernel after non-local coordinates have been received.
        if (domainWork.haveGpuBondedWork && !simulationWork.havePpDomainDecomposition)
        {
            fr->listedForcesGpu->setPbcAndlaunchKernel(fr->pbcType, box, fr->bMolPBC, stepWork);
        }

        /* launch local nonbonded work on GPU */
        wallcycle_start_nocount(wcycle, WallCycleCounter::LaunchGpuPp);
        wallcycle_sub_start_nocount(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);
        do_nb_verlet(fr, ic, enerd, stepWork, InteractionLocality::Local, enbvClearFNo, step, nrnb, wcycle);
        wallcycle_sub_stop(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);
        wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPp);
    }

    if (stepWork.haveGpuPmeOnThisRank)
    {
        // In PME GPU and mixed mode we launch FFT / gather after the
        // X copy/transform to allow overlap as well as after the GPU NB
        // launch to avoid FFT launch overhead hijacking the CPU and delaying
        // the nonbonded kernel.
        launchPmeGpuFftAndGather(fr->pmedata,
                                 lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)],
                                 wcycle,
                                 stepWork);
    }

    /* Communicate coordinates and sum dipole if necessary */
    if (simulationWork.havePpDomainDecomposition)
    {
        if (!stepWork.doNeighborSearch)
        {
            GpuEventSynchronizer* gpuCoordinateHaloLaunched = nullptr;
            if (stepWork.useGpuXHalo)
            {
                // The following must be called after local setCoordinates (which records an event
                // when the coordinate data has been copied to the device).
                gpuCoordinateHaloLaunched = communicateGpuHaloCoordinates(*cr, box, localXReadyOnDevice);

                if (domainWork.haveCpuNonLocalForceWork)
                {
                    // non-local part of coordinate buffer must be copied back to host for CPU work
                    stateGpu->copyCoordinatesFromGpu(
                            x.unpaddedArrayRef(), AtomLocality::NonLocal, gpuCoordinateHaloLaunched);
                }
            }
            else
            {
                if (simulationWork.useGpuUpdate)
                {
                    GMX_ASSERT(haveCopiedXFromGpu,
                               "a wait should only be triggered if copy has been scheduled");
                    const bool haveAlreadyWaited =
                            (stepWork.computePmeOnSeparateRank && !pmeSendCoordinatesFromGpu);
                    if (!haveAlreadyWaited)
                    {
                        stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
                    }
                }
                dd_move_x(cr->dd, box, x.unpaddedArrayRef(), wcycle);
            }

            if (stepWork.useGpuXBufferOps)
            {
                if (!stepWork.useGpuXHalo)
                {
                    stateGpu->copyCoordinatesToGpu(x.unpaddedArrayRef(), AtomLocality::NonLocal);
                }
                GpuEventSynchronizer* xReadyOnDeviceEvent = stateGpu->getCoordinatesReadyOnDeviceEvent(
                        AtomLocality::NonLocal, simulationWork, stepWork, gpuCoordinateHaloLaunched);
                if (stepWork.useGpuXHalo && domainWork.haveCpuNonLocalForceWork)
                {
                    /* We already enqueued an event for Gpu Halo exchange completion into the
                     * NonLocal stream when D2H copying the coordinates. */
                    xReadyOnDeviceEvent = nullptr;
                }
                nbv->convertCoordinatesGpu(
                        AtomLocality::NonLocal, stateGpu->getCoordinates(), xReadyOnDeviceEvent);
            }
            else
            {
                nbv->convertCoordinates(AtomLocality::NonLocal, x.unpaddedArrayRef());
            }
        }

        if (simulationWork.useGpuNonbonded)
        {

            if (!stepWork.useGpuXBufferOps)
            {
                wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPp);
                wallcycle_sub_start(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);
                Nbnxm::gpu_copy_xq_to_gpu(nbv->gpuNbv(), &nbv->nbat(), AtomLocality::NonLocal);
                wallcycle_sub_stop(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);
                wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPp);
            }

            if (domainWork.haveGpuBondedWork)
            {
                fr->listedForcesGpu->setPbcAndlaunchKernel(fr->pbcType, box, fr->bMolPBC, stepWork);
            }

            /* launch non-local nonbonded tasks on GPU */
            wallcycle_start_nocount(wcycle, WallCycleCounter::LaunchGpuPp);
            wallcycle_sub_start(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);
            do_nb_verlet(fr, ic, enerd, stepWork, InteractionLocality::NonLocal, enbvClearFNo, step, nrnb, wcycle);
            wallcycle_sub_stop(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);
            wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPp);
        }
    }

    if (simulationWork.useGpuNonbonded && stepWork.computeNonbondedForces)
    {
        /* launch D2H copy-back F */
        wallcycle_start_nocount(wcycle, WallCycleCounter::LaunchGpuPp);
        wallcycle_sub_start_nocount(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);

        if (simulationWork.havePpDomainDecomposition)
        {
            Nbnxm::gpu_launch_cpyback(nbv->gpuNbv(), &nbv->nbat(), stepWork, AtomLocality::NonLocal);
        }
        Nbnxm::gpu_launch_cpyback(nbv->gpuNbv(), &nbv->nbat(), stepWork, AtomLocality::Local);
        wallcycle_sub_stop(wcycle, WallCycleSubCounter::LaunchGpuNonBonded);

        if (domainWork.haveGpuBondedWork && stepWork.computeEnergy)
        {
            fr->listedForcesGpu->launchEnergyTransfer();
        }
        wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPp);
    }

    gmx::ArrayRef<const gmx::RVec> xWholeMolecules;
    if (fr->wholeMoleculeTransform)
    {
        xWholeMolecules = fr->wholeMoleculeTransform->wholeMoleculeCoordinates(x.unpaddedArrayRef(), box);
    }

    // For the rest of the CPU tasks that depend on GPU-update produced coordinates,
    // this wait ensures that the D2H transfer is complete.
    if (simulationWork.useGpuUpdate && !stepWork.doNeighborSearch)
    {
        const bool needCoordsOnHost = (runScheduleWork.domainWork.haveCpuLocalForceWork
                                       || stepWork.computeVirial || simulationWork.computeMuTot);
        const bool haveAlreadyWaited =
                simulationWork.useCpuHaloExchange
                || (stepWork.computePmeOnSeparateRank && !pmeSendCoordinatesFromGpu);
        if (needCoordsOnHost && !haveAlreadyWaited)
        {
            GMX_ASSERT(haveCopiedXFromGpu,
                       "a wait should only be triggered if copy has been scheduled");
            stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
        }
    }

    DipoleData dipoleData;

    if (simulationWork.computeMuTot)
    {
        const int start = 0;

        /* Calculate total (local) dipole moment in a temporary common array.
         * This makes it possible to sum them over nodes faster.
         */
        gmx::ArrayRef<const gmx::RVec> xRef =
                (xWholeMolecules.empty() ? x.unpaddedArrayRef() : xWholeMolecules);
        calc_mu(start,
                mdatoms->homenr,
                xRef,
                mdatoms->chargeA,
                mdatoms->chargeB,
                mdatoms->nChargePerturbed != 0,
                dipoleData.muStaging[0],
                dipoleData.muStaging[1]);

        reduceAndUpdateMuTot(
                &dipoleData, cr, (fr->efep != FreeEnergyPerturbationType::No), lambda, muTotal, ddBalanceRegionHandler);
    }

    /* Reset energies */
    reset_enerdata(enerd);

    if (haveDDAtomOrdering(*cr) && simulationWork.haveSeparatePmeRank)
    {
        wallcycle_start(wcycle, WallCycleCounter::PpDuringPme);
        dd_force_flop_start(cr->dd, nrnb);
    }

    if (inputrec.bRot)
    {
        wallcycle_start(wcycle, WallCycleCounter::Rot);
        do_rotation(cr, enforcedRotation, box, x.unpaddedConstArrayRef(), t, step, stepWork.doNeighborSearch);
        wallcycle_stop(wcycle, WallCycleCounter::Rot);
    }

    /* Start the force cycle counter.
     * Note that a different counter is used for dynamic load balancing.
     */
    wallcycle_start(wcycle, WallCycleCounter::Force);

    /* Set up and clear force outputs:
     * forceOutMtsLevel0:  everything except what is in the other two outputs
     * forceOutMtsLevel1:  PME-mesh and listed-forces group 1
     * forceOutNonbonded: non-bonded forces
     * Without multiple time stepping all point to the same object.
     * With multiple time-stepping the use is different for MTS fast (level0 only) and slow steps.
     *
     * Note that CPU force buffer clearing needs to happen after the completion of the
     * previous step's CPU force H2D transfer (prior to force reduction).
     * In the current code this is ensured by the earlier waitCoordinatesReadyOnHost()
     * which is sufficient, but it is suboptimal as it prevents overlap of the force clearing
     * with independent GPU work (integration/constraints, x D2H copy).
     */
    ForceOutputs forceOutMtsLevel0 = setupForceOutputs(
            &fr->forceHelperBuffers[0], force, domainWork, stepWork, simulationWork.havePpDomainDecomposition, wcycle);

    // Force output for MTS combined forces, only set at level1 MTS steps
    std::optional<ForceOutputs> forceOutMts =
            (simulationWork.useMts && stepWork.computeSlowForces)
                    ? std::optional(setupForceOutputs(&fr->forceHelperBuffers[1],
                                                      forceView->forceMtsCombinedWithPadding(),
                                                      domainWork,
                                                      stepWork,
                                                      simulationWork.havePpDomainDecomposition,
                                                      wcycle))
                    : std::nullopt;

    ForceOutputs* forceOutMtsLevel1 =
            simulationWork.useMts ? (stepWork.computeSlowForces ? &forceOutMts.value() : nullptr)
                                  : &forceOutMtsLevel0;

    const bool nonbondedAtMtsLevel1 = runScheduleWork.simulationWork.computeNonbondedAtMtsLevel1;

    ForceOutputs* forceOutNonbonded = nonbondedAtMtsLevel1 ? forceOutMtsLevel1 : &forceOutMtsLevel0;

    if (inputrec.bPull && pull_have_constraint(*pull_work))
    {
        clear_pull_forces(pull_work);
    }

    wallcycle_stop(wcycle, WallCycleCounter::Force);

    /* We calculate the non-bonded forces, when done on the CPU, here.
     * We do this before calling do_force_lowlevel, because in that
     * function, the listed forces are calculated before PME, which
     * does communication.  With this order, non-bonded and listed
     * force calculation imbalance can be balanced out by the domain
     * decomposition load balancing.
     */

    const bool useOrEmulateGpuNb = simulationWork.useGpuNonbonded || fr->nbv->emulateGpu();

    if (!useOrEmulateGpuNb)
    {
        wallcycle_start_nocount(wcycle, WallCycleCounter::Force);
        do_nb_verlet(fr, ic, enerd, stepWork, InteractionLocality::Local, enbvClearFYes, step, nrnb, wcycle);
        wallcycle_stop(wcycle, WallCycleCounter::Force);
    }

    if (stepWork.useGpuXHalo && domainWork.haveCpuNonLocalForceWork)
    {
        /* Wait for non-local coordinate data to be copied from device */
        stateGpu->waitCoordinatesReadyOnHost(AtomLocality::NonLocal);
    }

    wallcycle_start_nocount(wcycle, WallCycleCounter::Force);
    if (fr->efep != FreeEnergyPerturbationType::No && stepWork.computeNonbondedForces)
    {
        /* Calculate the local and non-local free energy interactions here.
         * Happens here on the CPU both with and without GPU.
         */
        nbv->dispatchFreeEnergyKernels(x,
                                       &forceOutNonbonded->forceWithShiftForces(),
                                       fr->use_simd_kernels,
                                       fr->ntype,
                                       *fr->ic,
                                       fr->shift_vec,
                                       fr->nbfp,
                                       fr->ljpme_c6grid,
                                       mdatoms->chargeA,
                                       mdatoms->chargeB,
                                       mdatoms->typeA,
                                       mdatoms->typeB,
                                       lambda,
                                       enerd,
                                       stepWork,
                                       nrnb);
    }

    if (stepWork.computeNonbondedForces && !useOrEmulateGpuNb)
    {
        if (simulationWork.havePpDomainDecomposition)
        {
            do_nb_verlet(fr, ic, enerd, stepWork, InteractionLocality::NonLocal, enbvClearFNo, step, nrnb, wcycle);
        }

        if (stepWork.computeForces)
        {
            /* Add all the non-bonded force to the normal force array.
             * This can be split into a local and a non-local part when overlapping
             * communication with calculation with domain decomposition.
             */
            wallcycle_stop(wcycle, WallCycleCounter::Force);
            nbv->atomdata_add_nbat_f_to_f(AtomLocality::All,
                                          forceOutNonbonded->forceWithShiftForces().force());
            wallcycle_start_nocount(wcycle, WallCycleCounter::Force);
        }

        /* If there are multiple fshift output buffers we need to reduce them */
        if (stepWork.computeVirial)
        {
            /* This is not in a subcounter because it takes a
               negligible and constant-sized amount of time */
            nbnxn_atomdata_add_nbat_fshift_to_fshift(
                    nbv->nbat(), forceOutNonbonded->forceWithShiftForces().shiftForces());
        }
    }

    // Compute wall interactions, when present.
    // Note: should be moved to special forces.
    if (inputrec.nwall && stepWork.computeNonbondedForces)
    {
        /* foreign lambda component for walls */
        real dvdl_walls = do_walls(inputrec,
                                   *fr,
                                   box,
                                   mdatoms->typeA,
                                   mdatoms->typeB,
                                   mdatoms->cENER,
                                   mdatoms->homenr,
                                   mdatoms->nPerturbed,
                                   x.unpaddedConstArrayRef(),
                                   &forceOutMtsLevel0.forceWithVirial(),
                                   lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)],
                                   enerd->grpp.energyGroupPairTerms[NonBondedEnergyTerms::LJSR],
                                   nrnb);
        enerd->dvdl_lin[FreeEnergyPerturbationCouplingType::Vdw] += dvdl_walls;
    }

    if (stepWork.computeListedForces)
    {
        /* Check whether we need to take into account PBC in listed interactions */
        bool needMolPbc = false;
        for (const auto& listedForces : fr->listedForces)
        {
            if (listedForces.haveCpuListedForces(*fr->fcdata))
            {
                needMolPbc = fr->bMolPBC;
            }
        }

        t_pbc pbc;

        if (needMolPbc)
        {
            /* Since all atoms are in the rectangular or triclinic unit-cell,
             * only single box vector shifts (2 in x) are required.
             */
            set_pbc_dd(&pbc, fr->pbcType, haveDDAtomOrdering(*cr) ? cr->dd->numCells : nullptr, TRUE, box);
        }

        for (int mtsIndex = 0; mtsIndex < (simulationWork.useMts && stepWork.computeSlowForces ? 2 : 1);
             mtsIndex++)
        {
            ListedForces& listedForces = fr->listedForces[mtsIndex];
            ForceOutputs& forceOut     = (mtsIndex == 0 ? forceOutMtsLevel0 : *forceOutMtsLevel1);
            listedForces.calculate(wcycle,
                                   box,
                                   cr,
                                   ms,
                                   x,
                                   xWholeMolecules,
                                   fr->fcdata.get(),
                                   hist,
                                   &forceOut,
                                   fr,
                                   &pbc,
                                   enerd,
                                   nrnb,
                                   lambda,
                                   mdatoms->chargeA,
                                   mdatoms->chargeB,
                                   makeConstArrayRef(mdatoms->bPerturbed),
                                   mdatoms->cENER,
                                   mdatoms->nPerturbed,
                                   haveDDAtomOrdering(*cr) ? cr->dd->globalAtomIndices.data() : nullptr,
                                   stepWork);
        }
    }

    if (stepWork.computeSlowForces)
    {
        longRangeNonbondeds->calculate(fr->pmedata,
                                       cr,
                                       x.unpaddedConstArrayRef(),
                                       &forceOutMtsLevel1->forceWithVirial(),
                                       enerd,
                                       box,
                                       lambda,
                                       dipoleData.muStateAB,
                                       stepWork,
                                       ddBalanceRegionHandler);
    }

    wallcycle_stop(wcycle, WallCycleCounter::Force);

    // VdW dispersion correction, only computed on main rank to avoid double counting
    if ((stepWork.computeEnergy || stepWork.computeVirial) && fr->dispersionCorrection && MAIN(cr))
    {
        // Calculate long range corrections to pressure and energy
        const DispersionCorrection::Correction correction = fr->dispersionCorrection->calculate(
                box, lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Vdw)]);

        if (stepWork.computeEnergy)
        {
            enerd->term[F_DISPCORR] = correction.energy;
            enerd->term[F_DVDL_VDW] += correction.dvdl;
            enerd->dvdl_lin[FreeEnergyPerturbationCouplingType::Vdw] += correction.dvdl;
        }
        if (stepWork.computeVirial)
        {
            correction.correctVirial(vir_force);
            enerd->term[F_PDISPCORR] = correction.pressure;
        }
    }

    const bool needToReceivePmeResultsFromSeparateRank = (PAR(cr) && stepWork.computePmeOnSeparateRank);
    const bool needToReceivePmeResults =
            (stepWork.haveGpuPmeOnThisRank || needToReceivePmeResultsFromSeparateRank);

    /* When running free energy perturbations steered by AWH and doing PME calculations on the
     * GPU we must wait for the PME calculation (dhdl) results to finish before sampling the
     * FEP dimension with AWH. */
    const bool needEarlyPmeResults = (awh != nullptr && awh->hasFepLambdaDimension() && needToReceivePmeResults
                                      && stepWork.computeEnergy && stepWork.computeSlowForces);
    if (needEarlyPmeResults)
    {
        if (stepWork.haveGpuPmeOnThisRank)
        {
            pmeGpuWaitAndReduce(fr->pmedata,
                                stepWork,
                                wcycle,
                                &forceOutMtsLevel1->forceWithVirial(),
                                enerd,
                                lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)]);
        }
        else if (needToReceivePmeResultsFromSeparateRank)
        {
            /* In case of node-splitting, the PP nodes receive the long-range
             * forces, virial and energy from the PME nodes here.
             */
            pme_receive_force_ener(fr,
                                   cr,
                                   &forceOutMtsLevel1->forceWithVirial(),
                                   enerd,
                                   simulationWork.useGpuPmePpCommunication,
                                   stepWork.useGpuPmeFReduction,
                                   wcycle);
        }
    }

    computeSpecialForces(fplog,
                         cr,
                         inputrec,
                         awh,
                         enforcedRotation,
                         imdSession,
                         pull_work,
                         step,
                         t,
                         wcycle,
                         fr->forceProviders,
                         box,
                         x.unpaddedArrayRef(),
                         mdatoms,
                         lambda,
                         stepWork,
                         &forceOutMtsLevel0.forceWithVirial(),
                         forceOutMtsLevel1 ? &forceOutMtsLevel1->forceWithVirial() : nullptr,
                         enerd,
                         ed,
                         stepWork.doNeighborSearch);

    if (simulationWork.havePpDomainDecomposition && stepWork.computeForces && stepWork.useGpuFHalo
        && domainWork.haveCpuLocalForceWork)
    {
        stateGpu->copyForcesToGpu(forceOutMtsLevel0.forceWithShiftForces().force(), AtomLocality::Local);
    }

    GMX_ASSERT(!(nonbondedAtMtsLevel1 && stepWork.useGpuFBufferOps),
               "The schedule below does not allow for nonbonded MTS with GPU buffer ops");
    GMX_ASSERT(!(nonbondedAtMtsLevel1 && stepWork.useGpuFHalo),
               "The schedule below does not allow for nonbonded MTS with GPU halo exchange");
    // Will store the amount of cycles spent waiting for the GPU that
    // will be later used in the DLB accounting.
    float cycles_wait_gpu = 0;
    if (useOrEmulateGpuNb && stepWork.computeNonbondedForces)
    {
        auto& forceWithShiftForces = forceOutNonbonded->forceWithShiftForces();

        /* wait for non-local forces (or calculate in emulation mode) */
        if (simulationWork.havePpDomainDecomposition)
        {
            if (simulationWork.useGpuNonbonded)
            {
                cycles_wait_gpu += Nbnxm::gpu_wait_finish_task(
                        nbv->gpuNbv(),
                        stepWork,
                        AtomLocality::NonLocal,
                        enerd->grpp.energyGroupPairTerms[NonBondedEnergyTerms::LJSR].data(),
                        enerd->grpp.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR].data(),
                        forceWithShiftForces.shiftForces(),
                        wcycle);
            }
            else
            {
                wallcycle_start_nocount(wcycle, WallCycleCounter::Force);
                do_nb_verlet(
                        fr, ic, enerd, stepWork, InteractionLocality::NonLocal, enbvClearFYes, step, nrnb, wcycle);
                wallcycle_stop(wcycle, WallCycleCounter::Force);
            }

            if (stepWork.useGpuFBufferOps)
            {
                if (domainWork.haveCpuNonLocalForceWork)
                {
                    stateGpu->copyForcesToGpu(forceOutMtsLevel0.forceWithShiftForces().force(),
                                              AtomLocality::NonLocal);
                }


                fr->gpuForceReduction[gmx::AtomLocality::NonLocal]->execute();

                if (!stepWork.useGpuFHalo)
                {
                    /* We don't explicitly wait for the forces to be reduced on device,
                     * but wait for them to finish copying to CPU instead.
                     * So, we manually consume the event, see Issue #3988. */
                    stateGpu->consumeForcesReducedOnDeviceEvent(AtomLocality::NonLocal);
                    // copy from GPU input for dd_move_f()
                    stateGpu->copyForcesFromGpu(forceOutMtsLevel0.forceWithShiftForces().force(),
                                                AtomLocality::NonLocal);
                }
            }
            else
            {
                nbv->atomdata_add_nbat_f_to_f(AtomLocality::NonLocal, forceWithShiftForces.force());
            }

            if (fr->nbv->emulateGpu() && stepWork.computeVirial)
            {
                nbnxn_atomdata_add_nbat_fshift_to_fshift(nbv->nbat(), forceWithShiftForces.shiftForces());
            }
        }
    }

    /* Combining the forces for multiple time stepping before the halo exchange, when possible,
     * avoids an extra halo exchange (when DD is used) and post-processing step.
     */
    if (stepWork.combineMtsForcesBeforeHaloExchange)
    {
        wallcycle_start_nocount(wcycle, WallCycleCounter::Force);
        combineMtsForces(getLocalAtomCount(cr->dd, *mdatoms, simulationWork.havePpDomainDecomposition),
                         force.unpaddedArrayRef(),
                         forceView->forceMtsCombined(),
                         inputrec.mtsLevels[1].stepFactor);
        wallcycle_stop(wcycle, WallCycleCounter::Force);
    }

    // With both nonbonded and PME offloaded a GPU on the same rank, we use
    // an alternating wait/reduction scheme.
    // When running free energy perturbations steered by AWH and calculating PME on GPU,
    // i.e. if needEarlyPmeResults == true, the PME results have already been reduced above.
    const bool alternateGpuWait = (!c_disableAlternatingWait && stepWork.haveGpuPmeOnThisRank
                                   && simulationWork.useGpuNonbonded && !simulationWork.havePpDomainDecomposition
                                   && !stepWork.useGpuFBufferOps && !needEarlyPmeResults);


    const int expectedLocalFReadyOnDeviceConsumptionCount = getExpectedLocalFReadyOnDeviceConsumptionCount(
            simulationWork, domainWork, stepWork, useOrEmulateGpuNb, alternateGpuWait);
    // If expectedLocalFReadyOnDeviceConsumptionCount == 0, stateGpu can be uninitialized
    if (expectedLocalFReadyOnDeviceConsumptionCount > 0)
    {
        stateGpu->setFReadyOnDeviceEventExpectedConsumptionCount(
                AtomLocality::Local, expectedLocalFReadyOnDeviceConsumptionCount);
    }

    if (simulationWork.havePpDomainDecomposition)
    {
        /* We are done with the CPU compute.
         * We will now communicate the non-local forces.
         * If we use a GPU this will overlap with GPU work, so in that case
         * we do not close the DD force balancing region here.
         */
        ddBalanceRegionHandler.closeAfterForceComputationCpu();

        if (stepWork.computeForces)
        {

            if (stepWork.useGpuFHalo)
            {
                // If there exist CPU forces, data from halo exchange should accumulate into these
                bool accumulateForces = domainWork.haveCpuLocalForceWork;
                gmx::FixedCapacityVector<GpuEventSynchronizer*, 2> gpuForceHaloDependencies;
                // completion of both H2D copy and clearing is signaled by fReadyOnDevice
                if (domainWork.haveCpuLocalForceWork || stepWork.clearGpuFBufferEarly)
                {
                    gpuForceHaloDependencies.push_back(stateGpu->fReadyOnDevice(AtomLocality::Local));
                }
                gpuForceHaloDependencies.push_back(stateGpu->fReducedOnDevice(AtomLocality::NonLocal));

                communicateGpuHaloForces(*cr, accumulateForces, &gpuForceHaloDependencies);
            }
            else
            {
                if (stepWork.useGpuFBufferOps)
                {
                    stateGpu->waitForcesReadyOnHost(AtomLocality::NonLocal);
                }

                // Without MTS or with MTS at slow steps with uncombined forces we need to
                // communicate the fast forces
                if (!simulationWork.useMts || !stepWork.combineMtsForcesBeforeHaloExchange)
                {
                    dd_move_f(cr->dd, &forceOutMtsLevel0.forceWithShiftForces(), wcycle);
                }
                // With MTS we need to communicate the slow or combined (in forceOutMtsLevel1) forces
                if (simulationWork.useMts && stepWork.computeSlowForces)
                {
                    dd_move_f(cr->dd, &forceOutMtsLevel1->forceWithShiftForces(), wcycle);
                }
            }
        }
    }

    if (alternateGpuWait)
    {
        alternatePmeNbGpuWaitReduce(fr->nbv.get(),
                                    fr->pmedata,
                                    forceOutNonbonded,
                                    forceOutMtsLevel1,
                                    enerd,
                                    lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)],
                                    stepWork,
                                    wcycle);
    }

    if (!alternateGpuWait && stepWork.haveGpuPmeOnThisRank && !needEarlyPmeResults)
    {
        pmeGpuWaitAndReduce(fr->pmedata,
                            stepWork,
                            wcycle,
                            &forceOutMtsLevel1->forceWithVirial(),
                            enerd,
                            lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Coul)]);
    }

    /* Wait for local GPU NB outputs on the non-alternating wait path */
    if (!alternateGpuWait && stepWork.computeNonbondedForces && simulationWork.useGpuNonbonded)
    {
        /* Measured overhead on CUDA and OpenCL with(out) GPU sharing
         * is between 0.5 and 1.5 Mcycles. So 2 MCycles is an overestimate,
         * but even with a step of 0.1 ms the difference is less than 1%
         * of the step time.
         */
        const float gpuWaitApiOverheadMargin = 2e6F; /* cycles */
        const float waitCycles               = Nbnxm::gpu_wait_finish_task(
                nbv->gpuNbv(),
                stepWork,
                AtomLocality::Local,
                enerd->grpp.energyGroupPairTerms[NonBondedEnergyTerms::LJSR].data(),
                enerd->grpp.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR].data(),
                forceOutNonbonded->forceWithShiftForces().shiftForces(),
                wcycle);

        if (ddBalanceRegionHandler.useBalancingRegion())
        {
            DdBalanceRegionWaitedForGpu waitedForGpu = DdBalanceRegionWaitedForGpu::yes;
            if (stepWork.computeForces && waitCycles <= gpuWaitApiOverheadMargin)
            {
                /* We measured few cycles, it could be that the kernel
                 * and transfer finished earlier and there was no actual
                 * wait time, only API call overhead.
                 * Then the actual time could be anywhere between 0 and
                 * cycles_wait_est. We will use half of cycles_wait_est.
                 */
                waitedForGpu = DdBalanceRegionWaitedForGpu::no;
            }
            ddBalanceRegionHandler.closeAfterForceComputationGpu(cycles_wait_gpu, waitedForGpu);
        }
    }

    if (fr->nbv->emulateGpu())
    {
        // NOTE: emulation kernel is not included in the balancing region,
        // but emulation mode does not target performance anyway
        wallcycle_start_nocount(wcycle, WallCycleCounter::Force);
        do_nb_verlet(fr,
                     ic,
                     enerd,
                     stepWork,
                     InteractionLocality::Local,
                     haveDDAtomOrdering(*cr) ? enbvClearFNo : enbvClearFYes,
                     step,
                     nrnb,
                     wcycle);
        wallcycle_stop(wcycle, WallCycleCounter::Force);
    }

    // If on GPU PME-PP comms path, receive forces from PME before GPU buffer ops
    // TODO refactor this and unify with below default-path call to the same function
    // When running free energy perturbations steered by AWH and calculating PME on GPU,
    // i.e. if needEarlyPmeResults == true, the PME results have already been reduced above.
    if (needToReceivePmeResultsFromSeparateRank && simulationWork.useGpuPmePpCommunication && !needEarlyPmeResults)
    {
        /* In case of node-splitting, the PP nodes receive the long-range
         * forces, virial and energy from the PME nodes here.
         */
        pme_receive_force_ener(fr,
                               cr,
                               &forceOutMtsLevel1->forceWithVirial(),
                               enerd,
                               simulationWork.useGpuPmePpCommunication,
                               stepWork.useGpuPmeFReduction,
                               wcycle);
    }


    /* Do the nonbonded GPU (or emulation) force buffer reduction
     * on the non-alternating path. */
    GMX_ASSERT(!(nonbondedAtMtsLevel1 && stepWork.useGpuFBufferOps),
               "The schedule below does not allow for nonbonded MTS with GPU buffer ops");
    if (useOrEmulateGpuNb && !alternateGpuWait)
    {
        if (stepWork.useGpuFBufferOps)
        {
            ArrayRef<gmx::RVec> forceWithShift = forceOutNonbonded->forceWithShiftForces().force();

            // TODO: move these steps as early as possible:
            // - CPU f H2D should be as soon as all CPU-side forces are done
            // - wait for force reduction does not need to block host (at least not here, it's sufficient to wait
            //   before the next CPU task that consumes the forces: vsite spread or update)
            // - copy is not perfomed if GPU force halo exchange is active, because it would overwrite the result
            //   of the halo exchange. In that case the copy is instead performed above, before the exchange.
            //   These should be unified.
            if (domainWork.haveLocalForceContribInCpuBuffer && !stepWork.useGpuFHalo)
            {
                stateGpu->copyForcesToGpu(forceWithShift, AtomLocality::Local);
            }

            if (stepWork.computeNonbondedForces)
            {
                fr->gpuForceReduction[gmx::AtomLocality::Local]->execute();
            }

            // Copy forces to host if they are needed for update or if virtual sites are enabled.
            // If there are vsites, we need to copy forces every step to spread vsite forces on host.
            // TODO: When the output flags will be included in step workload, this copy can be combined with the
            //       copy call done in sim_utils(...) for the output.
            // NOTE: If there are virtual sites, the forces are modified on host after this D2H copy. Hence,
            //       they should not be copied in do_md(...) for the output.
            if (!simulationWork.useGpuUpdate
                || (simulationWork.useGpuUpdate && haveDDAtomOrdering(*cr) && simulationWork.useCpuPmePpCommunication)
                || vsite)
            {
                if (stepWork.computeNonbondedForces)
                {
                    /* We have previously issued force reduction on the GPU, but we will
                     * not use this event, instead relying on the stream being in-order.
                     * Issue #3988. */
                    stateGpu->consumeForcesReducedOnDeviceEvent(AtomLocality::Local);
                }
                stateGpu->copyForcesFromGpu(forceWithShift, AtomLocality::Local);
                stateGpu->waitForcesReadyOnHost(AtomLocality::Local);
            }
        }
        else if (stepWork.computeNonbondedForces)
        {
            ArrayRef<gmx::RVec> forceWithShift = forceOutNonbonded->forceWithShiftForces().force();
            nbv->atomdata_add_nbat_f_to_f(AtomLocality::Local, forceWithShift);
        }
    }

    if (expectedLocalFReadyOnDeviceConsumptionCount > 0)
    {
        /* The same fReadyOnDevice device synchronizer is later used to track buffer clearing,
         * so we reset the expected consumption value back to the default (1). */
        stateGpu->setFReadyOnDeviceEventExpectedConsumptionCount(AtomLocality::Local, 1);
    }

    launchGpuEndOfStepTasks(
            nbv, fr->listedForcesGpu.get(), fr->pmedata, enerd, runScheduleWork, step, wcycle);

    if (haveDDAtomOrdering(*cr))
    {
        dd_force_flop_stop(cr->dd, nrnb);
    }

    const bool haveCombinedMtsForces = (stepWork.computeForces && simulationWork.useMts && stepWork.computeSlowForces
                                        && stepWork.combineMtsForcesBeforeHaloExchange);
    if (stepWork.computeForces)
    {
        postProcessForceWithShiftForces(
                nrnb, wcycle, box, x.unpaddedArrayRef(), &forceOutMtsLevel0, vir_force, *mdatoms, *fr, vsite, stepWork);

        if (simulationWork.useMts && stepWork.computeSlowForces && !haveCombinedMtsForces)
        {
            postProcessForceWithShiftForces(
                    nrnb, wcycle, box, x.unpaddedArrayRef(), forceOutMtsLevel1, vir_force, *mdatoms, *fr, vsite, stepWork);
        }
    }

    // TODO refactor this and unify with above GPU PME-PP / GPU update path call to the same function
    // When running free energy perturbations steered by AWH and calculating PME on GPU,
    // i.e. if needEarlyPmeResults == true, the PME results have already been reduced above.
    if (needToReceivePmeResultsFromSeparateRank && simulationWork.useCpuPmePpCommunication && !needEarlyPmeResults)
    {
        /* In case of node-splitting, the PP nodes receive the long-range
         * forces, virial and energy from the PME nodes here.
         */
        pme_receive_force_ener(fr,
                               cr,
                               &forceOutMtsLevel1->forceWithVirial(),
                               enerd,
                               simulationWork.useGpuPmePpCommunication,
                               false,
                               wcycle);
    }

    if (stepWork.computeForces)
    {
        /* If we don't use MTS or if we already combined the MTS forces before, we only
         * need to post-process one ForceOutputs object here, called forceOutCombined,
         * otherwise we have to post-process two outputs and then combine them.
         */
        ForceOutputs& forceOutCombined = (haveCombinedMtsForces ? forceOutMts.value() : forceOutMtsLevel0);
        postProcessForces(
                cr, step, nrnb, wcycle, box, x.unpaddedArrayRef(), &forceOutCombined, vir_force, mdatoms, fr, vsite, stepWork);

        if (simulationWork.useMts && stepWork.computeSlowForces && !haveCombinedMtsForces)
        {
            postProcessForces(
                    cr, step, nrnb, wcycle, box, x.unpaddedArrayRef(), forceOutMtsLevel1, vir_force, mdatoms, fr, vsite, stepWork);

            combineMtsForces(mdatoms->homenr,
                             force.unpaddedArrayRef(),
                             forceView->forceMtsCombined(),
                             inputrec.mtsLevels[1].stepFactor);
        }
    }

    if (stepWork.computeEnergy)
    {
        /* Compute the final potential energy terms */
        accumulatePotentialEnergies(enerd, lambda, inputrec.fepvals.get());

        if (!EI_TPI(inputrec.eI))
        {
            checkPotentialEnergyValidity(step, *enerd, inputrec);
        }
    }

    /* In case we don't have constraints and are using GPUs, the next balancing
     * region starts here.
     * Some "special" work at the end of do_force_cuts?, such as vsite spread,
     * virial calculation and COM pulling, is not thus not included in
     * the balance timing, which is ok as most tasks do communication.
     */
    ddBalanceRegionHandler.openBeforeForceComputationCpu(DdAllowBalanceRegionReopen::no);
}

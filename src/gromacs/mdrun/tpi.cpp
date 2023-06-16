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
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <array>

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
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
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
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcebuffers.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "legacysimulator.h"

//! Global max algorithm
static void global_max(t_commrec* cr, int* n)
{
    int *sum, i;

    snew(sum, cr->nnodes);
    sum[cr->nodeid] = *n;
    gmx_sumi(cr->nnodes, sum, cr);
    for (i = 0; i < cr->nnodes; i++)
    {
        *n = std::max(*n, sum[i]);
    }

    sfree(sum);
}

//! Reallocate arrays.
static void realloc_bins(double** bin, int* nbin, int nbin_new)
{
    int i;

    if (nbin_new != *nbin)
    {
        srenew(*bin, nbin_new);
        for (i = *nbin; i < nbin_new; i++)
        {
            (*bin)[i] = 0;
        }
        *nbin = nbin_new;
    }
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

namespace gmx
{

// TODO: Convert to use the nbnxm kernels by putting the system and the teset molecule on two separate search grids
void LegacySimulator::do_tpi()
{
    GMX_RELEASE_ASSERT(gmx_omp_nthreads_get(ModuleMultiThread::Default) == 1,
                       "TPI does not support OpenMP");

    gmx::ForceBuffers          f;
    real                       lambda, t, temp, beta, drmax, epot;
    double                     embU, sum_embU, *sum_UgembU, V, V_all, VembU_all;
    t_trxstatus*               status;
    t_trxframe                 rerun_fr;
    gmx_bool                   bDispCorr, bCharge, bRFExcl, bNotLastFrame, bStateChanged, bNS;
    tensor                     force_vir, shake_vir, vir, pres;
    int                        a_tp0, a_tp1, ngid, gid_tp, nener, e;
    rvec*                      x_mol;
    rvec                       mu_tot, x_init, dx;
    int                        nnodes, frame;
    int64_t                    frame_step_prev, frame_step;
    int64_t                    nsteps, stepblocksize = 0, step;
    int64_t                    seed;
    int                        i;
    FILE*                      fp_tpi = nullptr;
    char *                     ptr, *dump_pdb;
    std::vector<std::string>   leg;
    double                     dbl, dump_ener;
    gmx_bool                   bCavity;
    int                        nat_cavity  = 0, d;
    real *                     mass_cavity = nullptr, mass_tot;
    int                        nbin;
    double                     invbinw, *bin, refvolshift, logV, bUlogV;
    gmx_bool                   bEnergyOutOfBounds;
    std::array<std::string, 2> tpid_leg = { "direct", "reweighted" };
    auto*                      mdatoms  = mdAtoms_->mdatoms();

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

    /* Since there is no upper limit to the insertion energies,
     * we need to set an upper limit for the distribution output.
     */
    real bU_bin_limit      = 50;
    real bU_logV_bin_limit = bU_bin_limit + 10;

    nnodes = cr_->nnodes;

    gmx_mtop_generate_local_top(topGlobal_, top_, inputRec_->efep != FreeEnergyPerturbationType::No);

    const SimulationGroups* groups = &topGlobal_.groups;

    bCavity = (inputRec_->eI == IntegrationAlgorithm::TPIC);
    if (bCavity)
    {
        ptr = getenv("GMX_TPIC_MASSES");
        if (ptr == nullptr)
        {
            nat_cavity = 1;
        }
        else
        {
            /* Read (multiple) masses from env var GMX_TPIC_MASSES,
             * The center of mass of the last atoms is then used for TPIC.
             */
            nat_cavity = 0;
            while (sscanf(ptr, "%20lf%n", &dbl, &i) > 0)
            {
                srenew(mass_cavity, nat_cavity + 1);
                mass_cavity[nat_cavity] = dbl;
                fprintf(fpLog_, "mass[%d] = %f\n", nat_cavity + 1, mass_cavity[nat_cavity]);
                nat_cavity++;
                ptr += i;
            }
            if (nat_cavity == 0)
            {
                gmx_fatal(FARGS, "Found %d masses in GMX_TPIC_MASSES", nat_cavity);
            }
        }
    }

    /*
       init_em(fplog,TPI,inputrec,&lambda,nrnb,mu_tot,
       state_global->box,fr,mdatoms,top,cr,nfile,fnm,NULL,NULL);*/
    /* We never need full pbc for TPI */
    fr_->pbcType = PbcType::Xyz;
    /* Determine the temperature for the Boltzmann weighting */
    temp = constantEnsembleTemperature(*inputRec_);
    if (fpLog_)
    {
        for (i = 1; (i < inputRec_->opts.ngtc); i++)
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
    beta = 1.0 / (gmx::c_boltz * temp);

    /* Number of insertions per frame */
    nsteps = inputRec_->nsteps;

    /* Use the same neighborlist with more insertions points
     * in a sphere of radius drmax around the initial point
     */
    /* This should be a proper mdp parameter */
    drmax = inputRec_->rtpi;

    /* An environment variable can be set to dump all configurations
     * to pdb with an insertion energy <= this value.
     */
    dump_pdb  = getenv("GMX_TPI_DUMP");
    dump_ener = 0;
    if (dump_pdb)
    {
        sscanf(dump_pdb, "%20lf", &dump_ener);
    }

    atoms2md(topGlobal_, *inputRec_, -1, {}, topGlobal_.natoms, mdAtoms_);
    update_mdatoms(mdatoms, inputRec_->fepvals->init_lambda);

    f.resize(topGlobal_.natoms);

    /* Print to log file  */
    walltime_accounting_start_time(wallTimeAccounting_);
    wallcycle_start(wallCycleCounters_, WallCycleCounter::Run);
    print_start(fpLog_, cr_, wallTimeAccounting_, "Test Particle Insertion");

    /* The last charge group is the group to be inserted */
    const t_atoms& atomsToInsert = topGlobal_.moltype[topGlobal_.molblock.back().type].atoms;
    a_tp0                        = topGlobal_.natoms - atomsToInsert.nr;
    a_tp1                        = topGlobal_.natoms;
    if (debug)
    {
        fprintf(debug, "TPI atoms %d-%d\n", a_tp0, a_tp1);
    }

    auto x = makeArrayRef(stateGlobal_->x);

    if (usingPme(fr_->ic->eeltype))
    {
        gmx_pme_reinit_atoms(fr_->pmedata, a_tp0, {}, {});
    }

    /* With reacion-field we have distance dependent potentials
     * between excluded atoms, we need to add these separately
     * for the inserted molecule.
     */
    real rfExclusionEnergy = 0;
    if (usingRF(fr_->ic->eeltype))
    {
        rfExclusionEnergy = reactionFieldExclusionCorrection(x, *mdatoms, *fr_->ic, a_tp0);
        if (debug)
        {
            fprintf(debug, "RF exclusion correction for inserted molecule: %f kJ/mol\n", rfExclusionEnergy);
        }
    }

    snew(x_mol, a_tp1 - a_tp0);

    bDispCorr = (inputRec_->eDispCorr != DispersionCorrectionType::No);
    bCharge   = FALSE;
    for (i = a_tp0; i < a_tp1; i++)
    {
        /* Copy the coordinates of the molecule to be insterted */
        copy_rvec(x[i], x_mol[i - a_tp0]);
        /* Check if we need to print electrostatic energies */
        bCharge |= (mdatoms->chargeA[i] != 0 || (!mdatoms->chargeB.empty() && mdatoms->chargeB[i] != 0));
    }
    bRFExcl = (bCharge && usingRF(fr_->ic->eeltype));

    // Calculate the center of geometry of the molecule to insert
    rvec cog = { 0, 0, 0 };
    for (int a = a_tp0; a < a_tp1; a++)
    {
        rvec_inc(cog, x[a]);
    }
    svmul(1.0_real / (a_tp1 - a_tp0), cog, cog);
    real molRadius = 0;
    for (int a = a_tp0; a < a_tp1; a++)
    {
        molRadius = std::max(molRadius, distance2(x[a], cog));
    }
    molRadius = std::sqrt(molRadius);

    const real maxCutoff = std::max(inputRec_->rvdw, inputRec_->rcoulomb);
    if (bCavity)
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
        for (i = 0; i < a_tp1 - a_tp0; i++)
        {
            rvec_dec(x_mol[i], cog);
        }
    }

    if (fpLog_)
    {
        fprintf(fpLog_, "\nWill insert %d atoms %s partial charges\n", a_tp1 - a_tp0, bCharge ? "with" : "without");

        fprintf(fpLog_,
                "\nWill insert %" PRId64 " times in each frame of %s\n",
                nsteps,
                opt2fn("-rerun", nFile_, fnm_));
    }

    if (!bCavity)
    {
        if (inputRec_->nstlist > 1)
        {

            /* With the same pair list we insert in a sphere of radius rtpi  in different orientations */
            if (drmax == 0 && a_tp1 - a_tp0 == 1)
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

    ngid   = groups->groups[SimulationAtomGroupType::EnergyOutput].size();
    gid_tp = fr_->atomInfo[a_tp0] & gmx::sc_atomInfo_EnergyGroupIdMask;
    for (int a = a_tp0 + 1; a < a_tp1; a++)
    {
        if ((fr_->atomInfo[a] & gmx::sc_atomInfo_EnergyGroupIdMask) != gid_tp)
        {
            fprintf(fpLog_,
                    "NOTE: Atoms in the molecule to insert belong to different energy groups.\n"
                    "      Only contributions to the group of the first atom will be reported.\n");
            break;
        }
    }
    nener = 1 + ngid;
    if (bDispCorr)
    {
        nener += 1;
    }
    if (bCharge)
    {
        nener += ngid;
        if (bRFExcl)
        {
            nener += 1;
        }
        if (usingFullElectrostatics(fr_->ic->eeltype))
        {
            nener += 1;
        }
    }
    snew(sum_UgembU, nener);

    /* Copy the random seed set by the user */
    seed = inputRec_->ld_seed;

    gmx::ThreeFry2x64<16> rng(
            seed, gmx::RandomDomain::TestParticleInsertion); // 16 bits internal counter => 2^16 * 2 = 131072 values per stream
    gmx::UniformRealDistribution<real> dist;

    if (MAIN(cr_))
    {
        fp_tpi = xvgropen(opt2fn("-tpi", nFile_, fnm_),
                          "TPI energies",
                          "Time (ps)",
                          "(kJ mol\\S-1\\N) / (nm\\S3\\N)",
                          oenv_);
        xvgr_subtitle(fp_tpi, "f. are averages over one frame", oenv_);
        leg.emplace_back("-kT log(<Ve\\S-\\betaU\\N>/<V>)");
        leg.emplace_back("f. -kT log<e\\S-\\betaU\\N>");
        leg.emplace_back("f. <e\\S-\\betaU\\N>");
        leg.emplace_back("#f. V");
        leg.emplace_back("f. <Ue\\S-\\betaU\\N>");
        for (i = 0; i < ngid; i++)
        {
            leg.emplace_back(gmx::formatString(
                    "f. <U\\sVdW %s\\Ne\\S-\\betaU\\N>",
                    *(groups->groupNames[groups->groups[SimulationAtomGroupType::EnergyOutput][i]])));
        }
        if (bDispCorr)
        {
            leg.emplace_back("f. <U\\sdisp c\\Ne\\S-\\betaU\\N>");
        }
        if (bCharge)
        {
            for (i = 0; i < ngid; i++)
            {
                leg.emplace_back(gmx::formatString(
                        "f. <U\\sCoul %s\\Ne\\S-\\betaU\\N>",
                        *(groups->groupNames[groups->groups[SimulationAtomGroupType::EnergyOutput][i]])));
            }
            if (bRFExcl)
            {
                leg.emplace_back("f. <U\\sRF excl\\Ne\\S-\\betaU\\N>");
            }
            if (usingFullElectrostatics(fr_->ic->eeltype))
            {
                leg.emplace_back("f. <U\\sCoul recip\\Ne\\S-\\betaU\\N>");
            }
        }
        xvgrLegend(fp_tpi, leg, oenv_);
    }
    clear_rvec(x_init);
    V_all     = 0;
    VembU_all = 0;

    invbinw = 10;
    nbin    = 10;
    snew(bin, nbin);

    /* Avoid frame step numbers <= -1 */
    frame_step_prev = -1;

    bNotLastFrame =
            read_first_frame(oenv_, &status, opt2fn("-rerun", nFile_, fnm_), &rerun_fr, TRX_NEED_X);
    frame = 0;

    if (rerun_fr.natoms - (bCavity ? nat_cavity : 0) != mdatoms->nr - (a_tp1 - a_tp0))
    {
        gmx_fatal(FARGS,
                  "Number of atoms in trajectory (%d)%s "
                  "is not equal the number in the run input file (%d) "
                  "minus the number of atoms to insert (%d)\n",
                  rerun_fr.natoms,
                  bCavity ? " minus one" : "",
                  mdatoms->nr,
                  a_tp1 - a_tp0);
    }

    refvolshift = log(det(rerun_fr.box));

    switch (inputRec_->eI)
    {
        case IntegrationAlgorithm::TPI: stepblocksize = inputRec_->nstlist; break;
        case IntegrationAlgorithm::TPIC: stepblocksize = 1; break;
        default: gmx_fatal(FARGS, "Unknown integrator %s", enumValueToString(inputRec_->eI));
    }

    while (bNotLastFrame)
    {
        frame_step = rerun_fr.step;
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

        lambda = rerun_fr.lambda;
        t      = rerun_fr.time;

        sum_embU = 0;
        for (e = 0; e < nener; e++)
        {
            sum_UgembU[e] = 0;
        }

        /* Copy the coordinates from the input trajectory */
        auto x = makeArrayRef(stateGlobal_->x);
        for (i = 0; i < rerun_fr.natoms; i++)
        {
            copy_rvec(rerun_fr.x[i], x[i]);
        }
        copy_mat(rerun_fr.box, stateGlobal_->box);
        const matrix& box = stateGlobal_->box;

        V    = det(box);
        logV = log(V);

        bStateChanged = TRUE;
        bNS           = TRUE;

        put_atoms_in_box(fr_->pbcType, box, x);

        /* Put all atoms except for the inserted ones on the grid */
        rvec vzero       = { 0, 0, 0 };
        rvec boxDiagonal = { box[XX][XX], box[YY][YY], box[ZZ][ZZ] };
        nbnxn_put_on_grid(
                fr_->nbv.get(), box, 0, vzero, boxDiagonal, nullptr, { 0, a_tp0 }, -1, fr_->atomInfo, x, 0, nullptr);

        step = cr_->nodeid * stepblocksize;
        while (step < nsteps)
        {
            /* Restart random engine using the frame and insertion step
             * as counters.
             * Note that we need to draw several random values per iteration,
             * but by using the internal subcounter functionality of ThreeFry2x64
             * we can draw 131072 unique 64-bit values before exhausting
             * the stream. This is a huge margin, and if something still goes
             * wrong you will get an exception when the stream is exhausted.
             */
            rng.restart(frame_step, step);
            dist.reset(); // erase any memory in the distribution

            if (!bCavity)
            {
                /* Random insertion in the whole volume */
                bNS = (step % inputRec_->nstlist == 0);
                if (bNS)
                {
                    /* Generate a random position in the box */
                    for (d = 0; d < DIM; d++)
                    {
                        x_init[d] = dist(rng) * stateGlobal_->box[d][d];
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
                    if (nat_cavity == 1)
                    {
                        /* Copy the location of the cavity */
                        copy_rvec(rerun_fr.x[rerun_fr.natoms - 1], x_init);
                    }
                    else
                    {
                        /* Determine the center of mass of the last molecule */
                        clear_rvec(x_init);
                        mass_tot = 0;
                        for (i = 0; i < nat_cavity; i++)
                        {
                            for (d = 0; d < DIM; d++)
                            {
                                x_init[d] += mass_cavity[i]
                                             * rerun_fr.x[rerun_fr.natoms - nat_cavity + i][d];
                            }
                            mass_tot += mass_cavity[i];
                        }
                        for (d = 0; d < DIM; d++)
                        {
                            x_init[d] /= mass_tot;
                        }
                    }
                }
            }

            if (bNS)
            {
                for (int a = a_tp0; a < a_tp1; a++)
                {
                    x[a] = x_init;
                }

                /* Put the inserted molecule on it's own search grid */
                nbnxn_put_on_grid(
                        fr_->nbv.get(), box, 1, x_init, x_init, nullptr, { a_tp0, a_tp1 }, -1, fr_->atomInfo, x, 0, nullptr);

                /* TODO: Avoid updating all atoms at every bNS step */
                fr_->nbv->setAtomProperties(mdatoms->typeA, mdatoms->chargeA, fr_->atomInfo);

                fr_->nbv->constructPairlist(InteractionLocality::Local, top_->excls, step, nrnb_);

                bNS = FALSE;
            }

            /* Add random displacement uniformly distributed in a sphere
             * of radius rtpi. We don't need to do this is we generate
             * a new center location every step.
             */
            rvec x_tp;
            if (bCavity || inputRec_->nstlist > 1)
            {
                /* Generate coordinates within |dx|=drmax of x_init */
                do
                {
                    for (d = 0; d < DIM; d++)
                    {
                        dx[d] = (2 * dist(rng) - 1) * drmax;
                    }
                } while (norm2(dx) > drmax * drmax);
                rvec_add(x_init, dx, x_tp);
            }
            else
            {
                copy_rvec(x_init, x_tp);
            }

            if (a_tp1 - a_tp0 == 1)
            {
                /* Insert a single atom, just copy the insertion location */
                copy_rvec(x_tp, x[a_tp0]);
            }
            else
            {
                /* Copy the coordinates from the top file */
                for (i = a_tp0; i < a_tp1; i++)
                {
                    copy_rvec(x_mol[i - a_tp0], x[i]);
                }
                /* Rotate the molecule randomly */
                real angleX = 2 * M_PI * dist(rng);
                /* Draw uniform random number for sin(angleY) instead of angleY itself in order to
                 * achieve the uniform distribution in the solid angles space. */
                real angleY = std::asin(2 * dist(rng) - 1);
                real angleZ = 2 * M_PI * dist(rng);
                rotate_conf(a_tp1 - a_tp0, stateGlobal_->x.rvec_array() + a_tp0, nullptr, angleX, angleY, angleZ);
                /* Shift to the insertion location */
                for (i = a_tp0; i < a_tp1; i++)
                {
                    rvec_inc(x[i], x_tp);
                }
            }

            /* Note: NonLocal refers to the inserted molecule */
            fr_->nbv->convertCoordinates(AtomLocality::NonLocal, x);
            fr_->longRangeNonbondeds->updateAfterPartition(*mdatoms);

            /* Clear some matrix variables  */
            clear_mat(force_vir);
            clear_mat(shake_vir);
            clear_mat(vir);
            clear_mat(pres);

            /* Calc energy (no forces) on new positions. */
            /* Make do_force do a single node force calculation */
            cr_->nnodes = 1;

            // TPI might place a particle so close that the potential
            // is infinite. Since this is intended to happen, we
            // temporarily suppress any exceptions that the processor
            // might raise, then restore the old behaviour.
            std::fenv_t floatingPointEnvironment;
            std::feholdexcept(&floatingPointEnvironment);
            do_force(fpLog_,
                     cr_,
                     ms_,
                     *inputRec_,
                     mdModulesNotifiers_,
                     nullptr,
                     nullptr,
                     imdSession_,
                     pullWork_,
                     step,
                     nrnb_,
                     wallCycleCounters_,
                     top_,
                     stateGlobal_->box,
                     stateGlobal_->x.arrayRefWithPadding(),
                     &stateGlobal_->hist,
                     &f.view(),
                     force_vir,
                     mdatoms,
                     enerd_,
                     stateGlobal_->lambda,
                     fr_,
                     runScheduleWork_,
                     nullptr,
                     mu_tot,
                     t,
                     nullptr,
                     fr_->longRangeNonbondeds.get(),
                     GMX_FORCE_NONBONDED | GMX_FORCE_ENERGY | (bStateChanged ? GMX_FORCE_STATECHANGED : 0),
                     DDBalanceRegionHandler(nullptr));
            std::feclearexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
            std::feupdateenv(&floatingPointEnvironment);

            cr_->nnodes   = nnodes;
            bStateChanged = FALSE;

            if (fr_->dispersionCorrection)
            {
                /* Calculate long range corrections to pressure and energy */
                const DispersionCorrection::Correction correction =
                        fr_->dispersionCorrection->calculate(stateGlobal_->box, lambda);
                /* figure out how to rearrange the next 4 lines MRS 8/4/2009 */
                enerd_->term[F_DISPCORR] = correction.energy;
                enerd_->term[F_EPOT] += correction.energy;
                enerd_->term[F_PRES] += correction.pressure;
                enerd_->term[F_DVDL] += correction.dvdl;
            }
            else
            {
                enerd_->term[F_DISPCORR] = 0;
            }
            if (usingRF(fr_->ic->eeltype))
            {
                enerd_->term[F_EPOT] += rfExclusionEnergy;
            }

            epot               = enerd_->term[F_EPOT];
            bEnergyOutOfBounds = FALSE;

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
                bEnergyOutOfBounds = TRUE;
            }
            if (bEnergyOutOfBounds)
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
                embU = exp(static_cast<double>(-beta * epot));
                sum_embU += embU;
                /* Determine the weighted energy contributions of each energy group */
                e = 0;
                sum_UgembU[e++] += epot * embU;
                if (fr_->haveBuckingham)
                {
                    for (i = 0; i < ngid; i++)
                    {
                        sum_UgembU[e++] +=
                                enerd_->grpp.energyGroupPairTerms[NonBondedEnergyTerms::BuckinghamSR]
                                                                 [GID(i, gid_tp, ngid)]
                                * embU;
                    }
                }
                else
                {
                    for (i = 0; i < ngid; i++)
                    {
                        sum_UgembU[e++] +=
                                enerd_->grpp.energyGroupPairTerms[NonBondedEnergyTerms::LJSR][GID(i, gid_tp, ngid)]
                                * embU;
                    }
                }
                if (bDispCorr)
                {
                    sum_UgembU[e++] += enerd_->term[F_DISPCORR] * embU;
                }
                if (bCharge)
                {
                    for (i = 0; i < ngid; i++)
                    {
                        sum_UgembU[e++] +=
                                enerd_->grpp.energyGroupPairTerms[NonBondedEnergyTerms::CoulombSR][GID(i, gid_tp, ngid)]
                                * embU;
                    }
                    if (bRFExcl)
                    {
                        sum_UgembU[e++] += rfExclusionEnergy * embU;
                    }
                    if (usingFullElectrostatics(fr_->ic->eeltype))
                    {
                        sum_UgembU[e++] += enerd_->term[F_COUL_RECIP] * embU;
                    }
                }
            }

            if (embU == 0 || beta * epot > bU_bin_limit)
            {
                bin[0]++;
            }
            else
            {
                i = gmx::roundToInt((bU_logV_bin_limit - (beta * epot - logV + refvolshift)) * invbinw);
                if (i < 0)
                {
                    i = 0;
                }
                if (i >= nbin)
                {
                    realloc_bins(&bin, &nbin, i + 10);
                }
                bin[i]++;
            }

            if (debug)
            {
                fprintf(debug,
                        "TPI %7d %12.5e %12.5f %12.5f %12.5f\n",
                        static_cast<int>(step),
                        epot,
                        x_tp[XX],
                        x_tp[YY],
                        x_tp[ZZ]);
            }

            if (dump_pdb && epot <= dump_ener)
            {
                auto str = gmx::formatString("t%g_step%d.pdb", t, static_cast<int>(step));
                auto str2 = gmx::formatString("t: %f step %d ener: %f", t, static_cast<int>(step), epot);
                write_sto_conf_mtop(str.c_str(),
                                    str2.c_str(),
                                    topGlobal_,
                                    stateGlobal_->x.rvec_array(),
                                    stateGlobal_->v.rvec_array(),
                                    inputRec_->pbcType,
                                    stateGlobal_->box);
            }

            step++;
            if ((step / stepblocksize) % cr_->nnodes != cr_->nodeid)
            {
                /* Skip all steps assigned to the other MPI ranks */
                step += (cr_->nnodes - 1) * stepblocksize;
            }
        }

        if (PAR(cr_))
        {
            /* When running in parallel sum the energies over the processes */
            gmx_sumd(1, &sum_embU, cr_);
            gmx_sumd(nener, sum_UgembU, cr_);
        }

        frame++;
        V_all += V;
        VembU_all += V * sum_embU / nsteps;

        if (fp_tpi)
        {
            if (mdrunOptions_.verbose || frame % 10 == 0 || frame < 10)
            {
                fprintf(stderr,
                        "mu %10.3e <mu> %10.3e\n",
                        -log(sum_embU / nsteps) / beta,
                        -log(VembU_all / V_all) / beta);
            }

            fprintf(fp_tpi,
                    "%10.3f %12.5e %12.5e %12.5e %12.5e",
                    t,
                    VembU_all == 0 ? 20 / beta : -log(VembU_all / V_all) / beta,
                    sum_embU == 0 ? 20 / beta : -log(sum_embU / nsteps) / beta,
                    sum_embU / nsteps,
                    V);
            for (e = 0; e < nener; e++)
            {
                fprintf(fp_tpi, " %12.5e", sum_UgembU[e] / nsteps);
            }
            fprintf(fp_tpi, "\n");
            fflush(fp_tpi);
        }

        bNotLastFrame = read_next_frame(oenv_, status, &rerun_fr);
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
        const double mu = -log(VembU_all / V_all) / beta;
        fprintf(fpLog_, "  <mu> = %12.5e kJ/mol\n", mu);

        if (!std::isfinite(mu))
        {
            fprintf(fpLog_,
                    "\nThe computed chemical potential is not finite - consider increasing the "
                    "number of steps and/or the number of frames to insert into.\n");
        }
    }

    /* Write the Boltzmann factor histogram */
    if (PAR(cr_))
    {
        /* When running in parallel sum the bins over the processes */
        i = nbin;
        global_max(cr_, &i);
        realloc_bins(&bin, &nbin, i);
        gmx_sumd(nbin, bin, cr_);
    }
    if (MAIN(cr_))
    {
        fp_tpi   = xvgropen(opt2fn("-tpid", nFile_, fnm_),
                          "TPI energy distribution",
                          "\\betaU - log(V/<V>)",
                          "count",
                          oenv_);
        auto str = gmx::formatString("number \\betaU > %g: %9.3e", bU_bin_limit, bin[0]);
        xvgr_subtitle(fp_tpi, str.c_str(), oenv_);
        xvgrLegend(fp_tpi, tpid_leg, oenv_);
        for (i = nbin - 1; i > 0; i--)
        {
            bUlogV = -i / invbinw + bU_logV_bin_limit - refvolshift + log(V_all / frame);
            fprintf(fp_tpi, "%6.2f %10d %12.5e\n", bUlogV, roundToInt(bin[i]), bin[i] * exp(-bUlogV) * V_all / VembU_all);
        }
        xvgrclose(fp_tpi);
    }
    sfree(bin);

    sfree(sum_UgembU);

    walltime_accounting_set_nsteps_done(wallTimeAccounting_, frame * inputRec_->nsteps);
}

} // namespace gmx

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
 * \brief This file defines the integrator for test particle insertion
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <algorithm>

#include <cfenv>

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
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"

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
        energy -= 0.5 * qi * qi * ic.c_rf;

        for (int j = i + 1; j < mdatoms.homenr; j++)
        {
            const real qj  = mdatoms.chargeA[j];
            const real rsq = distance2(x[i], x[j]);
            energy += qi * qj * (ic.k_rf * rsq - ic.c_rf);
        }
    }

    return ic.epsfac * energy;
}

namespace gmx
{

// TODO: Convert to use the nbnxm kernels by putting the system and the teset molecule on two separate search grids
void LegacySimulator::do_tpi()
{
    GMX_RELEASE_ASSERT(gmx_omp_nthreads_get(emntDefault) == 1, "TPI does not support OpenMP");

    gmx_localtop_t              top;
    PaddedHostVector<gmx::RVec> f{};
    real                        lambda, t, temp, beta, drmax, epot;
    double                      embU, sum_embU, *sum_UgembU, V, V_all, VembU_all;
    t_trxstatus*                status;
    t_trxframe                  rerun_fr;
    gmx_bool                    bDispCorr, bCharge, bRFExcl, bNotLastFrame, bStateChanged, bNS;
    tensor                      force_vir, shake_vir, vir, pres;
    int                         a_tp0, a_tp1, ngid, gid_tp, nener, e;
    rvec*                       x_mol;
    rvec                        mu_tot, x_init, dx;
    int                         nnodes, frame;
    int64_t                     frame_step_prev, frame_step;
    int64_t                     nsteps, stepblocksize = 0, step;
    int64_t                     seed;
    int                         i;
    FILE*                       fp_tpi = nullptr;
    char *                      ptr, *dump_pdb, **leg, str[STRLEN], str2[STRLEN];
    double                      dbl, dump_ener;
    gmx_bool                    bCavity;
    int                         nat_cavity  = 0, d;
    real *                      mass_cavity = nullptr, mass_tot;
    int                         nbin;
    double                      invbinw, *bin, refvolshift, logV, bUlogV;
    gmx_bool                    bEnergyOutOfBounds;
    const char*                 tpid_leg[2] = { "direct", "reweighted" };
    auto                        mdatoms     = mdAtoms->mdatoms();

    GMX_UNUSED_VALUE(outputProvider);

    GMX_LOG(mdlog.info)
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

    nnodes = cr->nnodes;

    gmx_mtop_generate_local_top(*top_global, &top, inputrec->efep != efepNO);

    SimulationGroups* groups = &top_global->groups;

    bCavity = (inputrec->eI == eiTPIC);
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
                fprintf(fplog, "mass[%d] = %f\n", nat_cavity + 1, mass_cavity[nat_cavity]);
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
    fr->ePBC = epbcXYZ;
    /* Determine the temperature for the Boltzmann weighting */
    temp = inputrec->opts.ref_t[0];
    if (fplog)
    {
        for (i = 1; (i < inputrec->opts.ngtc); i++)
        {
            if (inputrec->opts.ref_t[i] != temp)
            {
                fprintf(fplog,
                        "\nWARNING: The temperatures of the different temperature coupling groups "
                        "are not identical\n\n");
                fprintf(stderr,
                        "\nWARNING: The temperatures of the different temperature coupling groups "
                        "are not identical\n\n");
            }
        }
        fprintf(fplog, "\n  The temperature for test particle insertion is %.3f K\n\n", temp);
    }
    beta = 1.0 / (BOLTZ * temp);

    /* Number of insertions per frame */
    nsteps = inputrec->nsteps;

    /* Use the same neighborlist with more insertions points
     * in a sphere of radius drmax around the initial point
     */
    /* This should be a proper mdp parameter */
    drmax = inputrec->rtpi;

    /* An environment variable can be set to dump all configurations
     * to pdb with an insertion energy <= this value.
     */
    dump_pdb  = getenv("GMX_TPI_DUMP");
    dump_ener = 0;
    if (dump_pdb)
    {
        sscanf(dump_pdb, "%20lf", &dump_ener);
    }

    atoms2md(top_global, inputrec, -1, nullptr, top_global->natoms, mdAtoms);
    update_mdatoms(mdatoms, inputrec->fepvals->init_lambda);

    f.resizeWithPadding(top_global->natoms);

    /* Print to log file  */
    walltime_accounting_start_time(walltime_accounting);
    wallcycle_start(wcycle, ewcRUN);
    print_start(fplog, cr, walltime_accounting, "Test Particle Insertion");

    /* The last charge group is the group to be inserted */
    const t_atoms& atomsToInsert = top_global->moltype[top_global->molblock.back().type].atoms;
    a_tp0                        = top_global->natoms - atomsToInsert.nr;
    a_tp1                        = top_global->natoms;
    if (debug)
    {
        fprintf(debug, "TPI atoms %d-%d\n", a_tp0, a_tp1);
    }

    auto x = makeArrayRef(state_global->x);

    if (EEL_PME(fr->ic->eeltype))
    {
        gmx_pme_reinit_atoms(fr->pmedata, a_tp0, nullptr);
    }

    /* With reacion-field we have distance dependent potentials
     * between excluded atoms, we need to add these separately
     * for the inserted molecule.
     */
    real rfExclusionEnergy = 0;
    if (EEL_RF(fr->ic->eeltype))
    {
        rfExclusionEnergy = reactionFieldExclusionCorrection(x, *mdatoms, *fr->ic, a_tp0);
        if (debug)
        {
            fprintf(debug, "RF exclusion correction for inserted molecule: %f kJ/mol\n", rfExclusionEnergy);
        }
    }

    snew(x_mol, a_tp1 - a_tp0);

    bDispCorr = (inputrec->eDispCorr != edispcNO);
    bCharge   = FALSE;
    for (i = a_tp0; i < a_tp1; i++)
    {
        /* Copy the coordinates of the molecule to be insterted */
        copy_rvec(x[i], x_mol[i - a_tp0]);
        /* Check if we need to print electrostatic energies */
        bCharge |= (mdatoms->chargeA[i] != 0
                    || ((mdatoms->chargeB != nullptr) && mdatoms->chargeB[i] != 0));
    }
    bRFExcl = (bCharge && EEL_RF(fr->ic->eeltype));

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

    const real maxCutoff = std::max(inputrec->rvdw, inputrec->rcoulomb);
    if (bCavity)
    {
        if (norm(cog) > 0.5 * maxCutoff && fplog)
        {
            fprintf(fplog, "WARNING: Your TPI molecule is not centered at 0,0,0\n");
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

    if (fplog)
    {
        fprintf(fplog, "\nWill insert %d atoms %s partial charges\n", a_tp1 - a_tp0,
                bCharge ? "with" : "without");

        fprintf(fplog, "\nWill insert %" PRId64 " times in each frame of %s\n", nsteps,
                opt2fn("-rerun", nfile, fnm));
    }

    if (!bCavity)
    {
        if (inputrec->nstlist > 1)
        {

            /* With the same pair list we insert in a sphere of radius rtpi  in different orientations */
            if (drmax == 0 && a_tp1 - a_tp0 == 1)
            {
                gmx_fatal(FARGS,
                          "Re-using the neighborlist %d times for insertions of a single atom in a "
                          "sphere of radius %f does not make sense",
                          inputrec->nstlist, drmax);
            }
            if (fplog)
            {
                fprintf(fplog,
                        "Will use the same neighborlist for %d insertions in a sphere of radius "
                        "%f\n",
                        inputrec->nstlist, drmax);
            }
        }
    }
    else
    {
        if (fplog)
        {
            fprintf(fplog,
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
    inputrec->rlist = maxCutoff + 2 * inputrec->rtpi + 2 * molRadius;
    fr->rlist       = inputrec->rlist;
    fr->nbv->changePairlistRadii(inputrec->rlist, inputrec->rlist);

    ngid   = groups->groups[SimulationAtomGroupType::EnergyOutput].size();
    gid_tp = GET_CGINFO_GID(fr->cginfo[a_tp0]);
    for (int a = a_tp0 + 1; a < a_tp1; a++)
    {
        if (GET_CGINFO_GID(fr->cginfo[a]) != gid_tp)
        {
            fprintf(fplog,
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
        if (EEL_FULL(fr->ic->eeltype))
        {
            nener += 1;
        }
    }
    snew(sum_UgembU, nener);

    /* Copy the random seed set by the user */
    seed = inputrec->ld_seed;

    gmx::ThreeFry2x64<16> rng(
            seed, gmx::RandomDomain::TestParticleInsertion); // 16 bits internal counter => 2^16 * 2 = 131072 values per stream
    gmx::UniformRealDistribution<real> dist;

    if (MASTER(cr))
    {
        fp_tpi = xvgropen(opt2fn("-tpi", nfile, fnm), "TPI energies", "Time (ps)",
                          "(kJ mol\\S-1\\N) / (nm\\S3\\N)", oenv);
        xvgr_subtitle(fp_tpi, "f. are averages over one frame", oenv);
        snew(leg, 4 + nener);
        e = 0;
        sprintf(str, "-kT log(<Ve\\S-\\betaU\\N>/<V>)");
        leg[e++] = gmx_strdup(str);
        sprintf(str, "f. -kT log<e\\S-\\betaU\\N>");
        leg[e++] = gmx_strdup(str);
        sprintf(str, "f. <e\\S-\\betaU\\N>");
        leg[e++] = gmx_strdup(str);
        sprintf(str, "f. V");
        leg[e++] = gmx_strdup(str);
        sprintf(str, "f. <Ue\\S-\\betaU\\N>");
        leg[e++] = gmx_strdup(str);
        for (i = 0; i < ngid; i++)
        {
            sprintf(str, "f. <U\\sVdW %s\\Ne\\S-\\betaU\\N>",
                    *(groups->groupNames[groups->groups[SimulationAtomGroupType::EnergyOutput][i]]));
            leg[e++] = gmx_strdup(str);
        }
        if (bDispCorr)
        {
            sprintf(str, "f. <U\\sdisp c\\Ne\\S-\\betaU\\N>");
            leg[e++] = gmx_strdup(str);
        }
        if (bCharge)
        {
            for (i = 0; i < ngid; i++)
            {
                sprintf(str, "f. <U\\sCoul %s\\Ne\\S-\\betaU\\N>",
                        *(groups->groupNames[groups->groups[SimulationAtomGroupType::EnergyOutput][i]]));
                leg[e++] = gmx_strdup(str);
            }
            if (bRFExcl)
            {
                sprintf(str, "f. <U\\sRF excl\\Ne\\S-\\betaU\\N>");
                leg[e++] = gmx_strdup(str);
            }
            if (EEL_FULL(fr->ic->eeltype))
            {
                sprintf(str, "f. <U\\sCoul recip\\Ne\\S-\\betaU\\N>");
                leg[e++] = gmx_strdup(str);
            }
        }
        xvgr_legend(fp_tpi, 4 + nener, leg, oenv);
        for (i = 0; i < 4 + nener; i++)
        {
            sfree(leg[i]);
        }
        sfree(leg);
    }
    clear_rvec(x_init);
    V_all     = 0;
    VembU_all = 0;

    invbinw = 10;
    nbin    = 10;
    snew(bin, nbin);

    /* Avoid frame step numbers <= -1 */
    frame_step_prev = -1;

    bNotLastFrame = read_first_frame(oenv, &status, opt2fn("-rerun", nfile, fnm), &rerun_fr, TRX_NEED_X);
    frame         = 0;

    if (rerun_fr.natoms - (bCavity ? nat_cavity : 0) != mdatoms->nr - (a_tp1 - a_tp0))
    {
        gmx_fatal(FARGS,
                  "Number of atoms in trajectory (%d)%s "
                  "is not equal the number in the run input file (%d) "
                  "minus the number of atoms to insert (%d)\n",
                  rerun_fr.natoms, bCavity ? " minus one" : "", mdatoms->nr, a_tp1 - a_tp0);
    }

    refvolshift = log(det(rerun_fr.box));

    switch (inputrec->eI)
    {
        case eiTPI: stepblocksize = inputrec->nstlist; break;
        case eiTPIC: stepblocksize = 1; break;
        default: gmx_fatal(FARGS, "Unknown integrator %s", ei_names[inputrec->eI]);
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
        auto x = makeArrayRef(state_global->x);
        for (i = 0; i < rerun_fr.natoms; i++)
        {
            copy_rvec(rerun_fr.x[i], x[i]);
        }
        copy_mat(rerun_fr.box, state_global->box);
        const matrix& box = state_global->box;

        V    = det(box);
        logV = log(V);

        bStateChanged = TRUE;
        bNS           = TRUE;

        put_atoms_in_box(fr->ePBC, box, x);

        /* Put all atoms except for the inserted ones on the grid */
        rvec vzero       = { 0, 0, 0 };
        rvec boxDiagonal = { box[XX][XX], box[YY][YY], box[ZZ][ZZ] };
        nbnxn_put_on_grid(fr->nbv.get(), box, 0, vzero, boxDiagonal, nullptr, { 0, a_tp0 }, -1,
                          fr->cginfo, x, 0, nullptr);

        step = cr->nodeid * stepblocksize;
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
                bNS = (step % inputrec->nstlist == 0);
                if (bNS)
                {
                    /* Generate a random position in the box */
                    for (d = 0; d < DIM; d++)
                    {
                        x_init[d] = dist(rng) * state_global->box[d][d];
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
                nbnxn_put_on_grid(fr->nbv.get(), box, 1, x_init, x_init, nullptr, { a_tp0, a_tp1 },
                                  -1, fr->cginfo, x, 0, nullptr);

                /* TODO: Avoid updating all atoms at every bNS step */
                fr->nbv->setAtomProperties(*mdatoms, fr->cginfo);

                fr->nbv->constructPairlist(InteractionLocality::Local, &top.excls, step, nrnb);

                bNS = FALSE;
            }

            /* Add random displacement uniformly distributed in a sphere
             * of radius rtpi. We don't need to do this is we generate
             * a new center location every step.
             */
            rvec x_tp;
            if (bCavity || inputrec->nstlist > 1)
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
                real angleY = 2 * M_PI * dist(rng);
                real angleZ = 2 * M_PI * dist(rng);
                rotate_conf(a_tp1 - a_tp0, state_global->x.rvec_array() + a_tp0, nullptr, angleX,
                            angleY, angleZ);
                /* Shift to the insertion location */
                for (i = a_tp0; i < a_tp1; i++)
                {
                    rvec_inc(x[i], x_tp);
                }
            }

            /* Note: NonLocal refers to the inserted molecule */
            fr->nbv->convertCoordinates(AtomLocality::NonLocal, false, x);

            /* Clear some matrix variables  */
            clear_mat(force_vir);
            clear_mat(shake_vir);
            clear_mat(vir);
            clear_mat(pres);

            /* Calc energy (no forces) on new positions.
             * Since we only need the intermolecular energy
             * and the RF exclusion terms of the inserted molecule occur
             * within a single charge group we can pass NULL for the graph.
             * This also avoids shifts that would move charge groups
             * out of the box. */
            /* Make do_force do a single node force calculation */
            cr->nnodes = 1;

            // TPI might place a particle so close that the potential
            // is infinite. Since this is intended to happen, we
            // temporarily suppress any exceptions that the processor
            // might raise, then restore the old behaviour.
            std::fenv_t floatingPointEnvironment;
            std::feholdexcept(&floatingPointEnvironment);
            do_force(fplog, cr, ms, inputrec, nullptr, nullptr, imdSession, pull_work, step, nrnb,
                     wcycle, &top, state_global->box, state_global->x.arrayRefWithPadding(),
                     &state_global->hist, f.arrayRefWithPadding(), force_vir, mdatoms, enerd, fcd,
                     state_global->lambda, nullptr, fr, runScheduleWork, nullptr, mu_tot, t, nullptr,
                     GMX_FORCE_NONBONDED | GMX_FORCE_ENERGY | (bStateChanged ? GMX_FORCE_STATECHANGED : 0),
                     DDBalanceRegionHandler(nullptr));
            std::feclearexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
            std::feupdateenv(&floatingPointEnvironment);

            cr->nnodes    = nnodes;
            bStateChanged = FALSE;

            if (fr->dispersionCorrection)
            {
                /* Calculate long range corrections to pressure and energy */
                const DispersionCorrection::Correction correction =
                        fr->dispersionCorrection->calculate(state_global->box, lambda);
                /* figure out how to rearrange the next 4 lines MRS 8/4/2009 */
                enerd->term[F_DISPCORR] = correction.energy;
                enerd->term[F_EPOT] += correction.energy;
                enerd->term[F_PRES] += correction.pressure;
                enerd->term[F_DVDL] += correction.dvdl;
            }
            else
            {
                enerd->term[F_DISPCORR] = 0;
            }
            if (EEL_RF(fr->ic->eeltype))
            {
                enerd->term[F_EPOT] += rfExclusionEnergy;
            }

            epot               = enerd->term[F_EPOT];
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
                            "\n  time %.3f, step %d: non-finite energy %f, using exp(-bU)=0\n", t,
                            static_cast<int>(step), epot);
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
                if (fr->bBHAM)
                {
                    for (i = 0; i < ngid; i++)
                    {
                        sum_UgembU[e++] += enerd->grpp.ener[egBHAMSR][GID(i, gid_tp, ngid)] * embU;
                    }
                }
                else
                {
                    for (i = 0; i < ngid; i++)
                    {
                        sum_UgembU[e++] += enerd->grpp.ener[egLJSR][GID(i, gid_tp, ngid)] * embU;
                    }
                }
                if (bDispCorr)
                {
                    sum_UgembU[e++] += enerd->term[F_DISPCORR] * embU;
                }
                if (bCharge)
                {
                    for (i = 0; i < ngid; i++)
                    {
                        sum_UgembU[e++] += enerd->grpp.ener[egCOULSR][GID(i, gid_tp, ngid)] * embU;
                    }
                    if (bRFExcl)
                    {
                        sum_UgembU[e++] += rfExclusionEnergy * embU;
                    }
                    if (EEL_FULL(fr->ic->eeltype))
                    {
                        sum_UgembU[e++] += enerd->term[F_COUL_RECIP] * embU;
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
                fprintf(debug, "TPI %7d %12.5e %12.5f %12.5f %12.5f\n", static_cast<int>(step),
                        epot, x_tp[XX], x_tp[YY], x_tp[ZZ]);
            }

            if (dump_pdb && epot <= dump_ener)
            {
                sprintf(str, "t%g_step%d.pdb", t, static_cast<int>(step));
                sprintf(str2, "t: %f step %d ener: %f", t, static_cast<int>(step), epot);
                write_sto_conf_mtop(str, str2, top_global, state_global->x.rvec_array(),
                                    state_global->v.rvec_array(), inputrec->ePBC, state_global->box);
            }

            step++;
            if ((step / stepblocksize) % cr->nnodes != cr->nodeid)
            {
                /* Skip all steps assigned to the other MPI ranks */
                step += (cr->nnodes - 1) * stepblocksize;
            }
        }

        if (PAR(cr))
        {
            /* When running in parallel sum the energies over the processes */
            gmx_sumd(1, &sum_embU, cr);
            gmx_sumd(nener, sum_UgembU, cr);
        }

        frame++;
        V_all += V;
        VembU_all += V * sum_embU / nsteps;

        if (fp_tpi)
        {
            if (mdrunOptions.verbose || frame % 10 == 0 || frame < 10)
            {
                fprintf(stderr, "mu %10.3e <mu> %10.3e\n", -log(sum_embU / nsteps) / beta,
                        -log(VembU_all / V_all) / beta);
            }

            fprintf(fp_tpi, "%10.3f %12.5e %12.5e %12.5e %12.5e", t,
                    VembU_all == 0 ? 20 / beta : -log(VembU_all / V_all) / beta,
                    sum_embU == 0 ? 20 / beta : -log(sum_embU / nsteps) / beta, sum_embU / nsteps, V);
            for (e = 0; e < nener; e++)
            {
                fprintf(fp_tpi, " %12.5e", sum_UgembU[e] / nsteps);
            }
            fprintf(fp_tpi, "\n");
            fflush(fp_tpi);
        }

        bNotLastFrame = read_next_frame(oenv, status, &rerun_fr);
    } /* End of the loop  */
    walltime_accounting_end_time(walltime_accounting);

    close_trx(status);

    if (fp_tpi != nullptr)
    {
        xvgrclose(fp_tpi);
    }

    if (fplog != nullptr)
    {
        fprintf(fplog, "\n");
        fprintf(fplog, "  <V>  = %12.5e nm^3\n", V_all / frame);
        const double mu = -log(VembU_all / V_all) / beta;
        fprintf(fplog, "  <mu> = %12.5e kJ/mol\n", mu);

        if (!std::isfinite(mu))
        {
            fprintf(fplog,
                    "\nThe computed chemical potential is not finite - consider increasing the "
                    "number of steps and/or the number of frames to insert into.\n");
        }
    }

    /* Write the Boltzmann factor histogram */
    if (PAR(cr))
    {
        /* When running in parallel sum the bins over the processes */
        i = nbin;
        global_max(cr, &i);
        realloc_bins(&bin, &nbin, i);
        gmx_sumd(nbin, bin, cr);
    }
    if (MASTER(cr))
    {
        fp_tpi = xvgropen(opt2fn("-tpid", nfile, fnm), "TPI energy distribution",
                          "\\betaU - log(V/<V>)", "count", oenv);
        sprintf(str, "number \\betaU > %g: %9.3e", bU_bin_limit, bin[0]);
        xvgr_subtitle(fp_tpi, str, oenv);
        xvgr_legend(fp_tpi, 2, tpid_leg, oenv);
        for (i = nbin - 1; i > 0; i--)
        {
            bUlogV = -i / invbinw + bU_logV_bin_limit - refvolshift + log(V_all / frame);
            fprintf(fp_tpi, "%6.2f %10d %12.5e\n", bUlogV, roundToInt(bin[i]),
                    bin[i] * exp(-bUlogV) * V_all / VembU_all);
        }
        xvgrclose(fp_tpi);
    }
    sfree(bin);

    sfree(sum_UgembU);

    walltime_accounting_set_nsteps_done(walltime_accounting, frame * inputrec->nsteps);
}

} // namespace gmx

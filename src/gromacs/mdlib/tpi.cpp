/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/conformation-utilities.h"
#include "gromacs/legacyheaders/chargegroup.h"
#include "gromacs/legacyheaders/constr.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/mdatoms.h"
#include "gromacs/legacyheaders/mdebin.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/ns.h"
#include "gromacs/legacyheaders/tgroup.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/update.h"
#include "gromacs/legacyheaders/vsite.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/random/random.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static void global_max(t_commrec *cr, int *n)
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

static void realloc_bins(double **bin, int *nbin, int nbin_new)
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

double do_tpi(FILE *fplog, t_commrec *cr,
              int nfile, const t_filenm fnm[],
              const output_env_t oenv, gmx_bool bVerbose, gmx_bool gmx_unused bCompact,
              int gmx_unused nstglobalcomm,
              gmx_vsite_t gmx_unused *vsite, gmx_constr_t gmx_unused constr,
              int gmx_unused stepout,
              t_inputrec *inputrec,
              gmx_mtop_t *top_global, t_fcdata *fcd,
              t_state *state,
              t_mdatoms *mdatoms,
              t_nrnb *nrnb, gmx_wallcycle_t wcycle,
              gmx_edsam_t gmx_unused ed,
              t_forcerec *fr,
              int gmx_unused repl_ex_nst, int gmx_unused repl_ex_nex, int gmx_unused repl_ex_seed,
              gmx_membed_t gmx_unused membed,
              real gmx_unused cpt_period, real gmx_unused max_hours,
              int gmx_unused imdport,
              unsigned long gmx_unused Flags,
              gmx_walltime_accounting_t walltime_accounting)
{
    gmx_localtop_t *top;
    gmx_groups_t   *groups;
    gmx_enerdata_t *enerd;
    rvec           *f;
    real            lambda, t, temp, beta, drmax, epot;
    double          embU, sum_embU, *sum_UgembU, V, V_all, VembU_all;
    t_trxstatus    *status;
    t_trxframe      rerun_fr;
    gmx_bool        bDispCorr, bCharge, bRFExcl, bNotLastFrame, bStateChanged, bNS;
    tensor          force_vir, shake_vir, vir, pres;
    int             cg_tp, a_tp0, a_tp1, ngid, gid_tp, nener, e;
    rvec           *x_mol;
    rvec            mu_tot, x_init, dx, x_tp;
    int             nnodes, frame;
    gmx_int64_t     frame_step_prev, frame_step;
    gmx_int64_t     nsteps, stepblocksize = 0, step;
    gmx_int64_t     rnd_count_stride, rnd_count;
    gmx_int64_t     seed;
    double          rnd[4];
    int             i;
    FILE           *fp_tpi = NULL;
    char           *ptr, *dump_pdb, **leg, str[STRLEN], str2[STRLEN];
    double          dbl, dump_ener;
    gmx_bool        bCavity;
    int             nat_cavity  = 0, d;
    real           *mass_cavity = NULL, mass_tot;
    int             nbin;
    double          invbinw, *bin, refvolshift, logV, bUlogV;
    real            prescorr, enercorr, dvdlcorr;
    gmx_bool        bEnergyOutOfBounds;
    const char     *tpid_leg[2] = {"direct", "reweighted"};

    /* Since there is no upper limit to the insertion energies,
     * we need to set an upper limit for the distribution output.
     */
    real bU_bin_limit      = 50;
    real bU_logV_bin_limit = bU_bin_limit + 10;

    if (inputrec->cutoff_scheme == ecutsVERLET)
    {
        gmx_fatal(FARGS, "TPI does not work (yet) with the Verlet cut-off scheme");
    }

    nnodes = cr->nnodes;

    top = gmx_mtop_generate_local_top(top_global, inputrec);

    groups = &top_global->groups;

    bCavity = (inputrec->eI == eiTPIC);
    if (bCavity)
    {
        ptr = getenv("GMX_TPIC_MASSES");
        if (ptr == NULL)
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
                srenew(mass_cavity, nat_cavity+1);
                mass_cavity[nat_cavity] = dbl;
                fprintf(fplog, "mass[%d] = %f\n",
                        nat_cavity+1, mass_cavity[nat_cavity]);
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
       state->box,fr,mdatoms,top,cr,nfile,fnm,NULL,NULL);*/
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
                fprintf(fplog, "\nWARNING: The temperatures of the different temperature coupling groups are not identical\n\n");
                fprintf(stderr, "\nWARNING: The temperatures of the different temperature coupling groups are not identical\n\n");
            }
        }
        fprintf(fplog,
                "\n  The temperature for test particle insertion is %.3f K\n\n",
                temp);
    }
    beta = 1.0/(BOLTZ*temp);

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

    atoms2md(top_global, inputrec, 0, NULL, top_global->natoms, mdatoms);
    update_mdatoms(mdatoms, inputrec->fepvals->init_lambda);

    snew(enerd, 1);
    init_enerdata(groups->grps[egcENER].nr, inputrec->fepvals->n_lambda, enerd);
    snew(f, top_global->natoms);

    /* Print to log file  */
    walltime_accounting_start(walltime_accounting);
    wallcycle_start(wcycle, ewcRUN);
    print_start(fplog, cr, walltime_accounting, "Test Particle Insertion");

    /* The last charge group is the group to be inserted */
    cg_tp = top->cgs.nr - 1;
    a_tp0 = top->cgs.index[cg_tp];
    a_tp1 = top->cgs.index[cg_tp+1];
    if (debug)
    {
        fprintf(debug, "TPI cg %d, atoms %d-%d\n", cg_tp, a_tp0, a_tp1);
    }
    if (a_tp1 - a_tp0 > 1 &&
        (inputrec->rlist < inputrec->rcoulomb ||
         inputrec->rlist < inputrec->rvdw))
    {
        gmx_fatal(FARGS, "Can not do TPI for multi-atom molecule with a twin-range cut-off");
    }
    snew(x_mol, a_tp1-a_tp0);

    bDispCorr = (inputrec->eDispCorr != edispcNO);
    bCharge   = FALSE;
    for (i = a_tp0; i < a_tp1; i++)
    {
        /* Copy the coordinates of the molecule to be insterted */
        copy_rvec(state->x[i], x_mol[i-a_tp0]);
        /* Check if we need to print electrostatic energies */
        bCharge |= (mdatoms->chargeA[i] != 0 ||
                    (mdatoms->chargeB && mdatoms->chargeB[i] != 0));
    }
    bRFExcl = (bCharge && EEL_RF(fr->eeltype) && fr->eeltype != eelRF_NEC);

    calc_cgcm(fplog, cg_tp, cg_tp+1, &(top->cgs), state->x, fr->cg_cm);
    if (bCavity)
    {
        if (norm(fr->cg_cm[cg_tp]) > 0.5*inputrec->rlist && fplog)
        {
            fprintf(fplog, "WARNING: Your TPI molecule is not centered at 0,0,0\n");
            fprintf(stderr, "WARNING: Your TPI molecule is not centered at 0,0,0\n");
        }
    }
    else
    {
        /* Center the molecule to be inserted at zero */
        for (i = 0; i < a_tp1-a_tp0; i++)
        {
            rvec_dec(x_mol[i], fr->cg_cm[cg_tp]);
        }
    }

    if (fplog)
    {
        fprintf(fplog, "\nWill insert %d atoms %s partial charges\n",
                a_tp1-a_tp0, bCharge ? "with" : "without");

        fprintf(fplog, "\nWill insert %d times in each frame of %s\n",
                (int)nsteps, opt2fn("-rerun", nfile, fnm));
    }

    if (!bCavity)
    {
        if (inputrec->nstlist > 1)
        {
            if (drmax == 0 && a_tp1-a_tp0 == 1)
            {
                gmx_fatal(FARGS, "Re-using the neighborlist %d times for insertions of a single atom in a sphere of radius %f does not make sense", inputrec->nstlist, drmax);
            }
            if (fplog)
            {
                fprintf(fplog, "Will use the same neighborlist for %d insertions in a sphere of radius %f\n", inputrec->nstlist, drmax);
            }
        }
    }
    else
    {
        if (fplog)
        {
            fprintf(fplog, "Will insert randomly in a sphere of radius %f around the center of the cavity\n", drmax);
        }
    }

    ngid   = groups->grps[egcENER].nr;
    gid_tp = GET_CGINFO_GID(fr->cginfo[cg_tp]);
    nener  = 1 + ngid;
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
        if (EEL_FULL(fr->eeltype))
        {
            nener += 1;
        }
    }
    snew(sum_UgembU, nener);

    /* Copy the random seed set by the user */
    seed = inputrec->ld_seed;
    /* We use the frame step number as one random counter.
     * The second counter use the insertion (step) count. But we
     * need multiple random numbers per insertion. This number is
     * not fixed, since we generate random locations in a sphere
     * by putting locations in a cube and some of these fail.
     * A count of 20 is already extremely unlikely, so 10000 is
     * a safe margin for random numbers per insertion.
     */
    rnd_count_stride = 10000;

    if (MASTER(cr))
    {
        fp_tpi = xvgropen(opt2fn("-tpi", nfile, fnm),
                          "TPI energies", "Time (ps)",
                          "(kJ mol\\S-1\\N) / (nm\\S3\\N)", oenv);
        xvgr_subtitle(fp_tpi, "f. are averages over one frame", oenv);
        snew(leg, 4+nener);
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
                    *(groups->grpname[groups->grps[egcENER].nm_ind[i]]));
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
                        *(groups->grpname[groups->grps[egcENER].nm_ind[i]]));
                leg[e++] = gmx_strdup(str);
            }
            if (bRFExcl)
            {
                sprintf(str, "f. <U\\sRF excl\\Ne\\S-\\betaU\\N>");
                leg[e++] = gmx_strdup(str);
            }
            if (EEL_FULL(fr->eeltype))
            {
                sprintf(str, "f. <U\\sCoul recip\\Ne\\S-\\betaU\\N>");
                leg[e++] = gmx_strdup(str);
            }
        }
        xvgr_legend(fp_tpi, 4+nener, (const char**)leg, oenv);
        for (i = 0; i < 4+nener; i++)
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

    bNotLastFrame = read_first_frame(oenv, &status, opt2fn("-rerun", nfile, fnm),
                                     &rerun_fr, TRX_NEED_X);
    frame = 0;

    if (rerun_fr.natoms - (bCavity ? nat_cavity : 0) !=
        mdatoms->nr - (a_tp1 - a_tp0))
    {
        gmx_fatal(FARGS, "Number of atoms in trajectory (%d)%s "
                  "is not equal the number in the run input file (%d) "
                  "minus the number of atoms to insert (%d)\n",
                  rerun_fr.natoms, bCavity ? " minus one" : "",
                  mdatoms->nr, a_tp1-a_tp0);
    }

    refvolshift = log(det(rerun_fr.box));

    switch (inputrec->eI)
    {
        case eiTPI:
            stepblocksize = inputrec->nstlist;
            break;
        case eiTPIC:
            stepblocksize = 1;
            break;
        default:
            gmx_fatal(FARGS, "Unknown integrator %s", ei_names[inputrec->eI]);
    }

    while (bNotLastFrame)
    {
        frame_step      = rerun_fr.step;
        if (frame_step <= frame_step_prev)
        {
            /* We don't have step number in the trajectory file,
             * or we have constant or decreasing step numbers.
             * Ensure we have increasing step numbers, since we use
             * the step numbers as a counter for random numbers.
             */
            frame_step  = frame_step_prev + 1;
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
        for (i = 0; i < rerun_fr.natoms; i++)
        {
            copy_rvec(rerun_fr.x[i], state->x[i]);
        }
        copy_mat(rerun_fr.box, state->box);

        V    = det(state->box);
        logV = log(V);

        bStateChanged = TRUE;
        bNS           = TRUE;

        step = cr->nodeid*stepblocksize;
        while (step < nsteps)
        {
            /* Initialize the second counter for random numbers using
             * the insertion step index. This ensures that we get
             * the same random numbers independently of how many
             * MPI ranks we use. Also for the same seed, we get
             * the same initial random sequence for different nsteps.
             */
            rnd_count = step*rnd_count_stride;

            if (!bCavity)
            {
                /* Random insertion in the whole volume */
                bNS = (step % inputrec->nstlist == 0);
                if (bNS)
                {
                    /* Generate a random position in the box */
                    gmx_rng_cycle_2uniform(frame_step, rnd_count++, seed, RND_SEED_TPI, rnd);
                    gmx_rng_cycle_2uniform(frame_step, rnd_count++, seed, RND_SEED_TPI, rnd+2);
                    for (d = 0; d < DIM; d++)
                    {
                        x_init[d] = rnd[d]*state->box[d][d];
                    }
                }
                if (inputrec->nstlist == 1)
                {
                    copy_rvec(x_init, x_tp);
                }
                else
                {
                    /* Generate coordinates within |dx|=drmax of x_init */
                    do
                    {
                        gmx_rng_cycle_2uniform(frame_step, rnd_count++, seed, RND_SEED_TPI, rnd);
                        gmx_rng_cycle_2uniform(frame_step, rnd_count++, seed, RND_SEED_TPI, rnd+2);
                        for (d = 0; d < DIM; d++)
                        {
                            dx[d] = (2*rnd[d] - 1)*drmax;
                        }
                    }
                    while (norm2(dx) > drmax*drmax);
                    rvec_add(x_init, dx, x_tp);
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
                        copy_rvec(rerun_fr.x[rerun_fr.natoms-1], x_init);
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
                                x_init[d] +=
                                    mass_cavity[i]*rerun_fr.x[rerun_fr.natoms-nat_cavity+i][d];
                            }
                            mass_tot += mass_cavity[i];
                        }
                        for (d = 0; d < DIM; d++)
                        {
                            x_init[d] /= mass_tot;
                        }
                    }
                }
                /* Generate coordinates within |dx|=drmax of x_init */
                do
                {
                    gmx_rng_cycle_2uniform(frame_step, rnd_count++, seed, RND_SEED_TPI, rnd);
                    gmx_rng_cycle_2uniform(frame_step, rnd_count++, seed, RND_SEED_TPI, rnd+2);
                    for (d = 0; d < DIM; d++)
                    {
                        dx[d] = (2*rnd[d] - 1)*drmax;
                    }
                }
                while (norm2(dx) > drmax*drmax);
                rvec_add(x_init, dx, x_tp);
            }

            if (a_tp1 - a_tp0 == 1)
            {
                /* Insert a single atom, just copy the insertion location */
                copy_rvec(x_tp, state->x[a_tp0]);
            }
            else
            {
                /* Copy the coordinates from the top file */
                for (i = a_tp0; i < a_tp1; i++)
                {
                    copy_rvec(x_mol[i-a_tp0], state->x[i]);
                }
                /* Rotate the molecule randomly */
                gmx_rng_cycle_2uniform(frame_step, rnd_count++, seed, RND_SEED_TPI, rnd);
                gmx_rng_cycle_2uniform(frame_step, rnd_count++, seed, RND_SEED_TPI, rnd+2);
                rotate_conf(a_tp1-a_tp0, state->x+a_tp0, NULL,
                            2*M_PI*rnd[0],
                            2*M_PI*rnd[1],
                            2*M_PI*rnd[2]);
                /* Shift to the insertion location */
                for (i = a_tp0; i < a_tp1; i++)
                {
                    rvec_inc(state->x[i], x_tp);
                }
            }

            /* Clear some matrix variables  */
            clear_mat(force_vir);
            clear_mat(shake_vir);
            clear_mat(vir);
            clear_mat(pres);

            /* Set the charge group center of mass of the test particle */
            copy_rvec(x_init, fr->cg_cm[top->cgs.nr-1]);

            /* Calc energy (no forces) on new positions.
             * Since we only need the intermolecular energy
             * and the RF exclusion terms of the inserted molecule occur
             * within a single charge group we can pass NULL for the graph.
             * This also avoids shifts that would move charge groups
             * out of the box.
             *
             * Some checks above ensure than we can not have
             * twin-range interactions together with nstlist > 1,
             * therefore we do not need to remember the LR energies.
             */
            /* Make do_force do a single node force calculation */
            cr->nnodes = 1;
            do_force(fplog, cr, inputrec,
                     step, nrnb, wcycle, top, &top_global->groups,
                     state->box, state->x, &state->hist,
                     f, force_vir, mdatoms, enerd, fcd,
                     state->lambda,
                     NULL, fr, NULL, mu_tot, t, NULL, NULL, FALSE,
                     GMX_FORCE_NONBONDED | GMX_FORCE_ENERGY |
                     (bNS ? GMX_FORCE_DYNAMICBOX | GMX_FORCE_NS | GMX_FORCE_DO_LR : 0) |
                     (bStateChanged ? GMX_FORCE_STATECHANGED : 0));
            cr->nnodes    = nnodes;
            bStateChanged = FALSE;
            bNS           = FALSE;

            /* Calculate long range corrections to pressure and energy */
            calc_dispcorr(inputrec, fr, top_global->natoms, state->box,
                          lambda, pres, vir, &prescorr, &enercorr, &dvdlcorr);
            /* figure out how to rearrange the next 4 lines MRS 8/4/2009 */
            enerd->term[F_DISPCORR]  = enercorr;
            enerd->term[F_EPOT]     += enercorr;
            enerd->term[F_PRES]     += prescorr;
            enerd->term[F_DVDL_VDW] += dvdlcorr;

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
                    fprintf(debug, "\n  time %.3f, step %d: non-finite energy %f, using exp(-bU)=0\n", t, (int)step, epot);
                }
                embU = 0;
            }
            else
            {
                embU      = exp(-beta*epot);
                sum_embU += embU;
                /* Determine the weighted energy contributions of each energy group */
                e                = 0;
                sum_UgembU[e++] += epot*embU;
                if (fr->bBHAM)
                {
                    for (i = 0; i < ngid; i++)
                    {
                        sum_UgembU[e++] +=
                            (enerd->grpp.ener[egBHAMSR][GID(i, gid_tp, ngid)] +
                             enerd->grpp.ener[egBHAMLR][GID(i, gid_tp, ngid)])*embU;
                    }
                }
                else
                {
                    for (i = 0; i < ngid; i++)
                    {
                        sum_UgembU[e++] +=
                            (enerd->grpp.ener[egLJSR][GID(i, gid_tp, ngid)] +
                             enerd->grpp.ener[egLJLR][GID(i, gid_tp, ngid)])*embU;
                    }
                }
                if (bDispCorr)
                {
                    sum_UgembU[e++] += enerd->term[F_DISPCORR]*embU;
                }
                if (bCharge)
                {
                    for (i = 0; i < ngid; i++)
                    {
                        sum_UgembU[e++] +=
                            (enerd->grpp.ener[egCOULSR][GID(i, gid_tp, ngid)] +
                             enerd->grpp.ener[egCOULLR][GID(i, gid_tp, ngid)])*embU;
                    }
                    if (bRFExcl)
                    {
                        sum_UgembU[e++] += enerd->term[F_RF_EXCL]*embU;
                    }
                    if (EEL_FULL(fr->eeltype))
                    {
                        sum_UgembU[e++] += enerd->term[F_COUL_RECIP]*embU;
                    }
                }
            }

            if (embU == 0 || beta*epot > bU_bin_limit)
            {
                bin[0]++;
            }
            else
            {
                i = (int)((bU_logV_bin_limit
                           - (beta*epot - logV + refvolshift))*invbinw
                          + 0.5);
                if (i < 0)
                {
                    i = 0;
                }
                if (i >= nbin)
                {
                    realloc_bins(&bin, &nbin, i+10);
                }
                bin[i]++;
            }

            if (debug)
            {
                fprintf(debug, "TPI %7d %12.5e %12.5f %12.5f %12.5f\n",
                        (int)step, epot, x_tp[XX], x_tp[YY], x_tp[ZZ]);
            }

            if (dump_pdb && epot <= dump_ener)
            {
                sprintf(str, "t%g_step%d.pdb", t, (int)step);
                sprintf(str2, "t: %f step %d ener: %f", t, (int)step, epot);
                write_sto_conf_mtop(str, str2, top_global, state->x, state->v,
                                    inputrec->ePBC, state->box);
            }

            step++;
            if ((step/stepblocksize) % cr->nnodes != cr->nodeid)
            {
                /* Skip all steps assigned to the other MPI ranks */
                step += (cr->nnodes - 1)*stepblocksize;
            }
        }

        if (PAR(cr))
        {
            /* When running in parallel sum the energies over the processes */
            gmx_sumd(1,    &sum_embU, cr);
            gmx_sumd(nener, sum_UgembU, cr);
        }

        frame++;
        V_all     += V;
        VembU_all += V*sum_embU/nsteps;

        if (fp_tpi)
        {
            if (bVerbose || frame%10 == 0 || frame < 10)
            {
                fprintf(stderr, "mu %10.3e <mu> %10.3e\n",
                        -log(sum_embU/nsteps)/beta, -log(VembU_all/V_all)/beta);
            }

            fprintf(fp_tpi, "%10.3f %12.5e %12.5e %12.5e %12.5e",
                    t,
                    VembU_all == 0 ? 20/beta : -log(VembU_all/V_all)/beta,
                    sum_embU == 0  ? 20/beta : -log(sum_embU/nsteps)/beta,
                    sum_embU/nsteps, V);
            for (e = 0; e < nener; e++)
            {
                fprintf(fp_tpi, " %12.5e", sum_UgembU[e]/nsteps);
            }
            fprintf(fp_tpi, "\n");
            fflush(fp_tpi);
        }

        bNotLastFrame = read_next_frame(oenv, status, &rerun_fr);
    } /* End of the loop  */
    walltime_accounting_end(walltime_accounting);

    close_trj(status);

    if (fp_tpi != NULL)
    {
        xvgrclose(fp_tpi);
    }

    if (fplog != NULL)
    {
        fprintf(fplog, "\n");
        fprintf(fplog, "  <V>  = %12.5e nm^3\n", V_all/frame);
        fprintf(fplog, "  <mu> = %12.5e kJ/mol\n", -log(VembU_all/V_all)/beta);
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
        fp_tpi = xvgropen(opt2fn("-tpid", nfile, fnm),
                          "TPI energy distribution",
                          "\\betaU - log(V/<V>)", "count", oenv);
        sprintf(str, "number \\betaU > %g: %9.3e", bU_bin_limit, bin[0]);
        xvgr_subtitle(fp_tpi, str, oenv);
        xvgr_legend(fp_tpi, 2, (const char **)tpid_leg, oenv);
        for (i = nbin-1; i > 0; i--)
        {
            bUlogV = -i/invbinw + bU_logV_bin_limit - refvolshift + log(V_all/frame);
            fprintf(fp_tpi, "%6.2f %10d %12.5e\n",
                    bUlogV,
                    (int)(bin[i]+0.5),
                    bin[i]*exp(-bUlogV)*V_all/VembU_all);
        }
        xvgrclose(fp_tpi);
    }
    sfree(bin);

    sfree(sum_UgembU);

    walltime_accounting_set_nsteps_done(walltime_accounting, frame*inputrec->nsteps);

    return 0;
}

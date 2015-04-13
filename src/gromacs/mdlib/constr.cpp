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

#include "gromacs/legacyheaders/constr.h"

#include <assert.h>
#include <stdlib.h>

#include <cmath>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/splitter.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/invblock.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

typedef struct gmx_constr {
    int                ncon_tot;       /* The total number of constraints    */
    int                nflexcon;       /* The number of flexible constraints */
    int                n_at2con_mt;    /* The size of at2con = #moltypes     */
    t_blocka          *at2con_mt;      /* A list of atoms to constraints     */
    int                n_at2settle_mt; /* The size of at2settle = #moltypes  */
    int              **at2settle_mt;   /* A list of atoms to settles         */
    gmx_bool           bInterCGsettles;
    gmx_lincsdata_t    lincsd;         /* LINCS data                         */
    gmx_shakedata_t    shaked;         /* SHAKE data                         */
    gmx_settledata_t   settled;        /* SETTLE data                        */
    int                nblocks;        /* The number of SHAKE blocks         */
    int               *sblock;         /* The SHAKE blocks                   */
    int                sblock_nalloc;  /* The allocation size of sblock      */
    real              *lagr;           /* -2 times the Lagrange multipliers for SHAKE */
    int                lagr_nalloc;    /* The allocation size of lagr        */
    int                maxwarn;        /* The maximum number of warnings     */
    int                warncount_lincs;
    int                warncount_settle;
    gmx_edsam_t        ed;            /* The essential dynamics data        */

    tensor            *vir_r_m_dr_th; /* Thread local working data          */
    int               *settle_error;  /* Thread local working data          */

    gmx_mtop_t        *warn_mtop;     /* Only used for printing warnings    */
} t_gmx_constr;

typedef struct {
    atom_id iatom[3];
    atom_id blocknr;
} t_sortblock;

static int pcomp(const void *p1, const void *p2)
{
    int          db;
    atom_id      min1, min2, max1, max2;
    t_sortblock *a1 = (t_sortblock *)p1;
    t_sortblock *a2 = (t_sortblock *)p2;

    db = a1->blocknr-a2->blocknr;

    if (db != 0)
    {
        return db;
    }

    min1 = std::min(a1->iatom[1], a1->iatom[2]);
    max1 = std::max(a1->iatom[1], a1->iatom[2]);
    min2 = std::min(a2->iatom[1], a2->iatom[2]);
    max2 = std::max(a2->iatom[1], a2->iatom[2]);

    if (min1 == min2)
    {
        return max1-max2;
    }
    else
    {
        return min1-min2;
    }
}

int n_flexible_constraints(struct gmx_constr *constr)
{
    int nflexcon;

    if (constr)
    {
        nflexcon = constr->nflexcon;
    }
    else
    {
        nflexcon = 0;
    }

    return nflexcon;
}

static void clear_constraint_quantity_nonlocal(gmx_domdec_t *dd, rvec *q)
{
    int nonlocal_at_start, nonlocal_at_end, at;

    dd_get_constraint_range(dd, &nonlocal_at_start, &nonlocal_at_end);

    for (at = nonlocal_at_start; at < nonlocal_at_end; at++)
    {
        clear_rvec(q[at]);
    }
}

void too_many_constraint_warnings(int eConstrAlg, int warncount)
{
    gmx_fatal(FARGS,
              "Too many %s warnings (%d)\n"
              "If you know what you are doing you can %s"
              "set the environment variable GMX_MAXCONSTRWARN to -1,\n"
              "but normally it is better to fix the problem",
              (eConstrAlg == econtLINCS) ? "LINCS" : "SETTLE", warncount,
              (eConstrAlg == econtLINCS) ?
              "adjust the lincs warning threshold in your mdp file\nor " : "\n");
}

static void write_constr_pdb(const char *fn, const char *title,
                             gmx_mtop_t *mtop,
                             int start, int homenr, t_commrec *cr,
                             rvec x[], matrix box)
{
    char          fname[STRLEN];
    FILE         *out;
    int           dd_ac0 = 0, dd_ac1 = 0, i, ii, resnr;
    gmx_domdec_t *dd;
    char         *anm, *resnm;

    dd = NULL;
    if (DOMAINDECOMP(cr))
    {
        dd = cr->dd;
        dd_get_constraint_range(dd, &dd_ac0, &dd_ac1);
        start  = 0;
        homenr = dd_ac1;
    }

    if (PAR(cr))
    {
        sprintf(fname, "%s_n%d.pdb", fn, cr->sim_nodeid);
    }
    else
    {
        sprintf(fname, "%s.pdb", fn);
    }

    out = gmx_fio_fopen(fname, "w");

    fprintf(out, "TITLE     %s\n", title);
    gmx_write_pdb_box(out, -1, box);
    for (i = start; i < start+homenr; i++)
    {
        if (dd != NULL)
        {
            if (i >= dd->nat_home && i < dd_ac0)
            {
                continue;
            }
            ii = dd->gatindex[i];
        }
        else
        {
            ii = i;
        }
        gmx_mtop_atominfo_global(mtop, ii, &anm, &resnr, &resnm);
        gmx_fprintf_pdb_atomline(out, epdbATOM, ii+1, anm, ' ', resnm, ' ', resnr, ' ',
                                 10*x[i][XX], 10*x[i][YY], 10*x[i][ZZ], 1.0, 0.0, "");
    }
    fprintf(out, "TER\n");

    gmx_fio_fclose(out);
}

static void dump_confs(FILE *fplog, gmx_int64_t step, gmx_mtop_t *mtop,
                       int start, int homenr, t_commrec *cr,
                       rvec x[], rvec xprime[], matrix box)
{
    char  buf[256], buf2[22];

    char *env = getenv("GMX_SUPPRESS_DUMP");
    if (env)
    {
        return;
    }

    sprintf(buf, "step%sb", gmx_step_str(step, buf2));
    write_constr_pdb(buf, "initial coordinates",
                     mtop, start, homenr, cr, x, box);
    sprintf(buf, "step%sc", gmx_step_str(step, buf2));
    write_constr_pdb(buf, "coordinates after constraining",
                     mtop, start, homenr, cr, xprime, box);
    if (fplog)
    {
        fprintf(fplog, "Wrote pdb files with previous and current coordinates\n");
    }
    fprintf(stderr, "Wrote pdb files with previous and current coordinates\n");
}

static void pr_sortblock(FILE *fp, const char *title, int nsb, t_sortblock sb[])
{
    int i;

    fprintf(fp, "%s\n", title);
    for (i = 0; (i < nsb); i++)
    {
        fprintf(fp, "i: %5d, iatom: (%5d %5d %5d), blocknr: %5d\n",
                i, sb[i].iatom[0], sb[i].iatom[1], sb[i].iatom[2],
                sb[i].blocknr);
    }
}

gmx_bool constrain(FILE *fplog, gmx_bool bLog, gmx_bool bEner,
                   struct gmx_constr *constr,
                   t_idef *idef, t_inputrec *ir,
                   t_commrec *cr,
                   gmx_int64_t step, int delta_step,
                   real step_scaling,
                   t_mdatoms *md,
                   rvec *x, rvec *xprime, rvec *min_proj,
                   gmx_bool bMolPBC, matrix box,
                   real lambda, real *dvdlambda,
                   rvec *v, tensor *vir,
                   t_nrnb *nrnb, int econq)
{
    gmx_bool    bOK, bDump;
    int         start, homenr;
    int         i, j;
    int         settle_error;
    tensor      vir_r_m_dr;
    real        scaled_delta_t;
    real        invdt, vir_fac = 0, t;
    t_ilist    *settle;
    int         nsettle;
    t_pbc       pbc, *pbc_null;
    char        buf[22];
    int         nth, th;

    if (econq == econqForceDispl && !EI_ENERGY_MINIMIZATION(ir->eI))
    {
        gmx_incons("constrain called for forces displacements while not doing energy minimization, can not do this while the LINCS and SETTLE constraint connection matrices are mass weighted");
    }

    bOK   = TRUE;
    bDump = FALSE;

    start  = 0;
    homenr = md->homenr;

    scaled_delta_t = step_scaling * ir->delta_t;

    /* Prepare time step for use in constraint implementations, and
       avoid generating inf when ir->delta_t = 0. */
    if (ir->delta_t == 0)
    {
        invdt = 0.0;
    }
    else
    {
        invdt = 1.0/scaled_delta_t;
    }

    if (ir->efep != efepNO && EI_DYNAMICS(ir->eI))
    {
        /* Set the constraint lengths for the step at which this configuration
         * is meant to be. The invmasses should not be changed.
         */
        lambda += delta_step*ir->fepvals->delta_lambda;
    }

    if (vir != NULL)
    {
        clear_mat(vir_r_m_dr);
    }

    where();

    settle  = &idef->il[F_SETTLE];
    nsettle = settle->nr/(1+NRAL(F_SETTLE));

    if (nsettle > 0)
    {
        nth = gmx_omp_nthreads_get(emntSETTLE);
    }
    else
    {
        nth = 1;
    }

    if (nth > 1 && constr->vir_r_m_dr_th == NULL)
    {
        snew(constr->vir_r_m_dr_th, nth);
        snew(constr->settle_error, nth);
    }

    settle_error = -1;

    /* We do not need full pbc when constraints do not cross charge groups,
     * i.e. when dd->constraint_comm==NULL.
     * Note that PBC for constraints is different from PBC for bondeds.
     * For constraints there is both forward and backward communication.
     */
    if (ir->ePBC != epbcNONE &&
        (cr->dd || bMolPBC) && !(cr->dd && cr->dd->constraint_comm == NULL))
    {
        /* With pbc=screw the screw has been changed to a shift
         * by the constraint coordinate communication routine,
         * so that here we can use normal pbc.
         */
        pbc_null = set_pbc_dd(&pbc, ir->ePBC, cr->dd, FALSE, box);
    }
    else
    {
        pbc_null = NULL;
    }

    /* Communicate the coordinates required for the non-local constraints
     * for LINCS and/or SETTLE.
     */
    if (cr->dd)
    {
        dd_move_x_constraints(cr->dd, box, x, xprime, econq == econqCoord);

        if (v != NULL)
        {
            /* We need to initialize the non-local components of v.
             * We never actually use these values, but we do increment them,
             * so we should avoid uninitialized variables and overflows.
             */
            clear_constraint_quantity_nonlocal(cr->dd, v);
        }
    }

    if (constr->lincsd != NULL)
    {
        bOK = constrain_lincs(fplog, bLog, bEner, ir, step, constr->lincsd, md, cr,
                              x, xprime, min_proj,
                              box, pbc_null, lambda, dvdlambda,
                              invdt, v, vir != NULL, vir_r_m_dr,
                              econq, nrnb,
                              constr->maxwarn, &constr->warncount_lincs);
        if (!bOK && constr->maxwarn >= 0)
        {
            if (fplog != NULL)
            {
                fprintf(fplog, "Constraint error in algorithm %s at step %s\n",
                        econstr_names[econtLINCS], gmx_step_str(step, buf));
            }
            bDump = TRUE;
        }
    }

    if (constr->nblocks > 0)
    {
        switch (econq)
        {
            case (econqCoord):
                bOK = bshakef(fplog, constr->shaked,
                              md->invmass, constr->nblocks, constr->sblock,
                              idef, ir, x, xprime, nrnb,
                              constr->lagr, lambda, dvdlambda,
                              invdt, v, vir != NULL, vir_r_m_dr,
                              constr->maxwarn >= 0, econq);
                break;
            case (econqVeloc):
                bOK = bshakef(fplog, constr->shaked,
                              md->invmass, constr->nblocks, constr->sblock,
                              idef, ir, x, min_proj, nrnb,
                              constr->lagr, lambda, dvdlambda,
                              invdt, NULL, vir != NULL, vir_r_m_dr,
                              constr->maxwarn >= 0, econq);
                break;
            default:
                gmx_fatal(FARGS, "Internal error, SHAKE called for constraining something else than coordinates");
                break;
        }

        if (!bOK && constr->maxwarn >= 0)
        {
            if (fplog != NULL)
            {
                fprintf(fplog, "Constraint error in algorithm %s at step %s\n",
                        econstr_names[econtSHAKE], gmx_step_str(step, buf));
            }
            bDump = TRUE;
        }
    }

    if (nsettle > 0)
    {
        int calcvir_atom_end;

        if (vir == NULL)
        {
            calcvir_atom_end = 0;
        }
        else
        {
            calcvir_atom_end = md->homenr;
        }

        switch (econq)
        {
            case econqCoord:
#pragma omp parallel for num_threads(nth) schedule(static)
                for (th = 0; th < nth; th++)
                {
                    int start_th, end_th;

                    if (th > 0)
                    {
                        clear_mat(constr->vir_r_m_dr_th[th]);
                    }

                    start_th = (nsettle* th   )/nth;
                    end_th   = (nsettle*(th+1))/nth;
                    if (start_th >= 0 && end_th - start_th > 0)
                    {
                        csettle(constr->settled,
                                end_th-start_th,
                                settle->iatoms+start_th*(1+NRAL(F_SETTLE)),
                                pbc_null,
                                x[0], xprime[0],
                                invdt, v ? v[0] : NULL, calcvir_atom_end,
                                th == 0 ? vir_r_m_dr : constr->vir_r_m_dr_th[th],
                                th == 0 ? &settle_error : &constr->settle_error[th]);
                    }
                }
                inc_nrnb(nrnb, eNR_SETTLE, nsettle);
                if (v != NULL)
                {
                    inc_nrnb(nrnb, eNR_CONSTR_V, nsettle*3);
                }
                if (vir != NULL)
                {
                    inc_nrnb(nrnb, eNR_CONSTR_VIR, nsettle*3);
                }
                break;
            case econqVeloc:
            case econqDeriv:
            case econqForce:
            case econqForceDispl:
#pragma omp parallel for num_threads(nth) schedule(static)
                for (th = 0; th < nth; th++)
                {
                    int start_th, end_th;

                    if (th > 0)
                    {
                        clear_mat(constr->vir_r_m_dr_th[th]);
                    }

                    start_th = (nsettle* th   )/nth;
                    end_th   = (nsettle*(th+1))/nth;

                    if (start_th >= 0 && end_th - start_th > 0)
                    {
                        settle_proj(constr->settled, econq,
                                    end_th-start_th,
                                    settle->iatoms+start_th*(1+NRAL(F_SETTLE)),
                                    pbc_null,
                                    x,
                                    xprime, min_proj, calcvir_atom_end,
                                    th == 0 ? vir_r_m_dr : constr->vir_r_m_dr_th[th]);
                    }
                }
                /* This is an overestimate */
                inc_nrnb(nrnb, eNR_SETTLE, nsettle);
                break;
            case econqDeriv_FlexCon:
                /* Nothing to do, since the are no flexible constraints in settles */
                break;
            default:
                gmx_incons("Unknown constraint quantity for settle");
        }
    }

    if (settle->nr > 0)
    {
        /* Combine virial and error info of the other threads */
        for (i = 1; i < nth; i++)
        {
            settle_error = constr->settle_error[i];
        }
        if (vir != NULL)
        {
            for (i = 1; i < nth; i++)
            {
                m_add(vir_r_m_dr, constr->vir_r_m_dr_th[i], vir_r_m_dr);
            }
        }

        if (econq == econqCoord && settle_error >= 0)
        {
            bOK = FALSE;
            if (constr->maxwarn >= 0)
            {
                char buf[256];
                sprintf(buf,
                        "\nstep " "%" GMX_PRId64 ": Water molecule starting at atom %d can not be "
                        "settled.\nCheck for bad contacts and/or reduce the timestep if appropriate.\n",
                        step, ddglatnr(cr->dd, settle->iatoms[settle_error*(1+NRAL(F_SETTLE))+1]));
                if (fplog)
                {
                    fprintf(fplog, "%s", buf);
                }
                fprintf(stderr, "%s", buf);
                constr->warncount_settle++;
                if (constr->warncount_settle > constr->maxwarn)
                {
                    too_many_constraint_warnings(-1, constr->warncount_settle);
                }
                bDump = TRUE;
            }
        }
    }

    if (vir != NULL)
    {
        /* The normal uses of constrain() pass step_scaling = 1.0.
         * The call to constrain() for SD1 that passes step_scaling =
         * 0.5 also passes vir = NULL, so cannot reach this
         * assertion. This assertion should remain until someone knows
         * that this path works for their intended purpose, and then
         * they can use scaled_delta_t instead of ir->delta_t
         * below. */
        assert(gmx_within_tol(step_scaling, 1.0, GMX_REAL_EPS));
        switch (econq)
        {
            case econqCoord:
                vir_fac = 0.5/(ir->delta_t*ir->delta_t);
                break;
            case econqVeloc:
                vir_fac = 0.5/ir->delta_t;
                break;
            case econqForce:
            case econqForceDispl:
                vir_fac = 0.5;
                break;
            default:
                gmx_incons("Unsupported constraint quantity for virial");
        }

        if (EI_VV(ir->eI))
        {
            vir_fac *= 2;  /* only constraining over half the distance here */
        }
        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                (*vir)[i][j] = vir_fac*vir_r_m_dr[i][j];
            }
        }
    }

    if (bDump)
    {
        dump_confs(fplog, step, constr->warn_mtop, start, homenr, cr, x, xprime, box);
    }

    if (econq == econqCoord)
    {
        if (ir->bPull && pull_have_constraint(ir->pull_work))
        {
            if (EI_DYNAMICS(ir->eI))
            {
                t = ir->init_t + (step + delta_step)*ir->delta_t;
            }
            else
            {
                t = ir->init_t;
            }
            set_pbc(&pbc, ir->ePBC, box);
            pull_constraint(ir->pull_work, md, &pbc, cr, ir->delta_t, t, x, xprime, v, *vir);
        }
        if (constr->ed && delta_step > 0)
        {
            /* apply the essential dynamcs constraints here */
            do_edsam(ir, step, cr, xprime, v, box, constr->ed);
        }
    }

    return bOK;
}

real *constr_rmsd_data(struct gmx_constr *constr)
{
    if (constr->lincsd)
    {
        return lincs_rmsd_data(constr->lincsd);
    }
    else
    {
        return NULL;
    }
}

real constr_rmsd(struct gmx_constr *constr, gmx_bool bSD2)
{
    if (constr->lincsd)
    {
        return lincs_rmsd(constr->lincsd, bSD2);
    }
    else
    {
        return 0;
    }
}

static void make_shake_sblock_serial(struct gmx_constr *constr,
                                     t_idef *idef, t_mdatoms *md)
{
    int          i, j, m, ncons;
    int          bstart, bnr;
    t_blocka     sblocks;
    t_sortblock *sb;
    t_iatom     *iatom;
    atom_id     *inv_sblock;

    /* Since we are processing the local topology,
     * the F_CONSTRNC ilist has been concatenated to the F_CONSTR ilist.
     */
    ncons = idef->il[F_CONSTR].nr/3;

    init_blocka(&sblocks);
    gen_sblocks(NULL, 0, md->homenr, idef, &sblocks, FALSE);

    /*
       bstart=(idef->nodeid > 0) ? blocks->multinr[idef->nodeid-1] : 0;
       nblocks=blocks->multinr[idef->nodeid] - bstart;
     */
    bstart          = 0;
    constr->nblocks = sblocks.nr;
    if (debug)
    {
        fprintf(debug, "ncons: %d, bstart: %d, nblocks: %d\n",
                ncons, bstart, constr->nblocks);
    }

    /* Calculate block number for each atom */
    inv_sblock = make_invblocka(&sblocks, md->nr);

    done_blocka(&sblocks);

    /* Store the block number in temp array and
     * sort the constraints in order of the sblock number
     * and the atom numbers, really sorting a segment of the array!
     */
#ifdef DEBUGIDEF
    pr_idef(fplog, 0, "Before Sort", idef);
#endif
    iatom = idef->il[F_CONSTR].iatoms;
    snew(sb, ncons);
    for (i = 0; (i < ncons); i++, iatom += 3)
    {
        for (m = 0; (m < 3); m++)
        {
            sb[i].iatom[m] = iatom[m];
        }
        sb[i].blocknr = inv_sblock[iatom[1]];
    }

    /* Now sort the blocks */
    if (debug)
    {
        pr_sortblock(debug, "Before sorting", ncons, sb);
        fprintf(debug, "Going to sort constraints\n");
    }

    qsort(sb, ncons, (size_t)sizeof(*sb), pcomp);

    if (debug)
    {
        pr_sortblock(debug, "After sorting", ncons, sb);
    }

    iatom = idef->il[F_CONSTR].iatoms;
    for (i = 0; (i < ncons); i++, iatom += 3)
    {
        for (m = 0; (m < 3); m++)
        {
            iatom[m] = sb[i].iatom[m];
        }
    }
#ifdef DEBUGIDEF
    pr_idef(fplog, 0, "After Sort", idef);
#endif

    j = 0;
    snew(constr->sblock, constr->nblocks+1);
    bnr = -2;
    for (i = 0; (i < ncons); i++)
    {
        if (sb[i].blocknr != bnr)
        {
            bnr                 = sb[i].blocknr;
            constr->sblock[j++] = 3*i;
        }
    }
    /* Last block... */
    constr->sblock[j++] = 3*ncons;

    if (j != (constr->nblocks+1))
    {
        fprintf(stderr, "bstart: %d\n", bstart);
        fprintf(stderr, "j: %d, nblocks: %d, ncons: %d\n",
                j, constr->nblocks, ncons);
        for (i = 0; (i < ncons); i++)
        {
            fprintf(stderr, "i: %5d  sb[i].blocknr: %5d\n", i, sb[i].blocknr);
        }
        for (j = 0; (j <= constr->nblocks); j++)
        {
            fprintf(stderr, "sblock[%3d]=%5d\n", j, (int)constr->sblock[j]);
        }
        gmx_fatal(FARGS, "DEATH HORROR: "
                  "sblocks does not match idef->il[F_CONSTR]");
    }
    sfree(sb);
    sfree(inv_sblock);
}

static void make_shake_sblock_dd(struct gmx_constr *constr,
                                 t_ilist *ilcon, t_block *cgs,
                                 gmx_domdec_t *dd)
{
    int      ncons, c, cg;
    t_iatom *iatom;

    if (dd->ncg_home+1 > constr->sblock_nalloc)
    {
        constr->sblock_nalloc = over_alloc_dd(dd->ncg_home+1);
        srenew(constr->sblock, constr->sblock_nalloc);
    }

    ncons           = ilcon->nr/3;
    iatom           = ilcon->iatoms;
    constr->nblocks = 0;
    cg              = 0;
    for (c = 0; c < ncons; c++)
    {
        if (c == 0 || iatom[1] >= cgs->index[cg+1])
        {
            constr->sblock[constr->nblocks++] = 3*c;
            while (iatom[1] >= cgs->index[cg+1])
            {
                cg++;
            }
        }
        iatom += 3;
    }
    constr->sblock[constr->nblocks] = 3*ncons;
}

t_blocka make_at2con(int start, int natoms,
                     t_ilist *ilist, t_iparams *iparams,
                     gmx_bool bDynamics, int *nflexiblecons)
{
    int      *count, ncon, con, con_tot, nflexcon, ftype, i, a;
    t_iatom  *ia;
    t_blocka  at2con;
    gmx_bool  bFlexCon;

    snew(count, natoms);
    nflexcon = 0;
    for (ftype = F_CONSTR; ftype <= F_CONSTRNC; ftype++)
    {
        ncon = ilist[ftype].nr/3;
        ia   = ilist[ftype].iatoms;
        for (con = 0; con < ncon; con++)
        {
            bFlexCon = (iparams[ia[0]].constr.dA == 0 &&
                        iparams[ia[0]].constr.dB == 0);
            if (bFlexCon)
            {
                nflexcon++;
            }
            if (bDynamics || !bFlexCon)
            {
                for (i = 1; i < 3; i++)
                {
                    a = ia[i] - start;
                    count[a]++;
                }
            }
            ia += 3;
        }
    }
    *nflexiblecons = nflexcon;

    at2con.nr           = natoms;
    at2con.nalloc_index = at2con.nr+1;
    snew(at2con.index, at2con.nalloc_index);
    at2con.index[0] = 0;
    for (a = 0; a < natoms; a++)
    {
        at2con.index[a+1] = at2con.index[a] + count[a];
        count[a]          = 0;
    }
    at2con.nra      = at2con.index[natoms];
    at2con.nalloc_a = at2con.nra;
    snew(at2con.a, at2con.nalloc_a);

    /* The F_CONSTRNC constraints have constraint numbers
     * that continue after the last F_CONSTR constraint.
     */
    con_tot = 0;
    for (ftype = F_CONSTR; ftype <= F_CONSTRNC; ftype++)
    {
        ncon = ilist[ftype].nr/3;
        ia   = ilist[ftype].iatoms;
        for (con = 0; con < ncon; con++)
        {
            bFlexCon = (iparams[ia[0]].constr.dA == 0 &&
                        iparams[ia[0]].constr.dB == 0);
            if (bDynamics || !bFlexCon)
            {
                for (i = 1; i < 3; i++)
                {
                    a = ia[i] - start;
                    at2con.a[at2con.index[a]+count[a]++] = con_tot;
                }
            }
            con_tot++;
            ia += 3;
        }
    }

    sfree(count);

    return at2con;
}

static int *make_at2settle(int natoms, const t_ilist *ilist)
{
    int *at2s;
    int  a, stride, s;

    snew(at2s, natoms);
    /* Set all to no settle */
    for (a = 0; a < natoms; a++)
    {
        at2s[a] = -1;
    }

    stride = 1 + NRAL(F_SETTLE);

    for (s = 0; s < ilist->nr; s += stride)
    {
        at2s[ilist->iatoms[s+1]] = s/stride;
        at2s[ilist->iatoms[s+2]] = s/stride;
        at2s[ilist->iatoms[s+3]] = s/stride;
    }

    return at2s;
}

void set_constraints(struct gmx_constr *constr,
                     gmx_localtop_t *top, t_inputrec *ir,
                     t_mdatoms *md, t_commrec *cr)
{
    t_idef  *idef;
    int      ncons;
    t_ilist *settle;
    int      iO, iH;

    idef = &top->idef;

    if (constr->ncon_tot > 0)
    {
        /* We are using the local topology,
         * so there are only F_CONSTR constraints.
         */
        ncons = idef->il[F_CONSTR].nr/3;

        /* With DD we might also need to call LINCS with ncons=0 for
         * communicating coordinates to other nodes that do have constraints.
         */
        if (ir->eConstrAlg == econtLINCS)
        {
            set_lincs(idef, md, EI_DYNAMICS(ir->eI), cr, constr->lincsd);
        }
        if (ir->eConstrAlg == econtSHAKE)
        {
            if (cr->dd)
            {
                make_shake_sblock_dd(constr, &idef->il[F_CONSTR], &top->cgs, cr->dd);
            }
            else
            {
                make_shake_sblock_serial(constr, idef, md);
            }
            if (ncons > constr->lagr_nalloc)
            {
                constr->lagr_nalloc = over_alloc_dd(ncons);
                srenew(constr->lagr, constr->lagr_nalloc);
            }
        }
    }

    if (idef->il[F_SETTLE].nr > 0 && constr->settled == NULL)
    {
        settle          = &idef->il[F_SETTLE];
        iO              = settle->iatoms[1];
        iH              = settle->iatoms[2];
        constr->settled =
            settle_init(md->massT[iO], md->massT[iH],
                        md->invmass[iO], md->invmass[iH],
                        idef->iparams[settle->iatoms[0]].settle.doh,
                        idef->iparams[settle->iatoms[0]].settle.dhh);
    }

    /* Make a selection of the local atoms for essential dynamics */
    if (constr->ed && cr->dd)
    {
        dd_make_local_ed_indices(cr->dd, constr->ed);
    }
}

static void constr_recur(t_blocka *at2con,
                         t_ilist *ilist, t_iparams *iparams, gmx_bool bTopB,
                         int at, int depth, int nc, int *path,
                         real r0, real r1, real *r2max,
                         int *count)
{
    int      ncon1;
    t_iatom *ia1, *ia2;
    int      c, con, a1;
    gmx_bool bUse;
    t_iatom *ia;
    real     len, rn0, rn1;

    (*count)++;

    ncon1 = ilist[F_CONSTR].nr/3;
    ia1   = ilist[F_CONSTR].iatoms;
    ia2   = ilist[F_CONSTRNC].iatoms;

    /* Loop over all constraints connected to this atom */
    for (c = at2con->index[at]; c < at2con->index[at+1]; c++)
    {
        con = at2con->a[c];
        /* Do not walk over already used constraints */
        bUse = TRUE;
        for (a1 = 0; a1 < depth; a1++)
        {
            if (con == path[a1])
            {
                bUse = FALSE;
            }
        }
        if (bUse)
        {
            ia = constr_iatomptr(ncon1, ia1, ia2, con);
            /* Flexible constraints currently have length 0, which is incorrect */
            if (!bTopB)
            {
                len = iparams[ia[0]].constr.dA;
            }
            else
            {
                len = iparams[ia[0]].constr.dB;
            }
            /* In the worst case the bond directions alternate */
            if (nc % 2 == 0)
            {
                rn0 = r0 + len;
                rn1 = r1;
            }
            else
            {
                rn0 = r0;
                rn1 = r1 + len;
            }
            /* Assume angles of 120 degrees between all bonds */
            if (rn0*rn0 + rn1*rn1 + rn0*rn1 > *r2max)
            {
                *r2max = rn0*rn0 + rn1*rn1 + r0*rn1;
                if (debug)
                {
                    fprintf(debug, "Found longer constraint distance: r0 %5.3f r1 %5.3f rmax %5.3f\n", rn0, rn1, sqrt(*r2max));
                    for (a1 = 0; a1 < depth; a1++)
                    {
                        fprintf(debug, " %d %5.3f",
                                path[a1],
                                iparams[constr_iatomptr(ncon1, ia1, ia2, con)[0]].constr.dA);
                    }
                    fprintf(debug, " %d %5.3f\n", con, len);
                }
            }
            /* Limit the number of recursions to 1000*nc,
             * so a call does not take more than a second,
             * even for highly connected systems.
             */
            if (depth + 1 < nc && *count < 1000*nc)
            {
                if (ia[1] == at)
                {
                    a1 = ia[2];
                }
                else
                {
                    a1 = ia[1];
                }
                /* Recursion */
                path[depth] = con;
                constr_recur(at2con, ilist, iparams,
                             bTopB, a1, depth+1, nc, path, rn0, rn1, r2max, count);
                path[depth] = -1;
            }
        }
    }
}

static real constr_r_max_moltype(gmx_moltype_t *molt, t_iparams *iparams,
                                 t_inputrec *ir)
{
    int      natoms, nflexcon, *path, at, count;

    t_blocka at2con;
    real     r0, r1, r2maxA, r2maxB, rmax, lam0, lam1;

    if (molt->ilist[F_CONSTR].nr   == 0 &&
        molt->ilist[F_CONSTRNC].nr == 0)
    {
        return 0;
    }

    natoms = molt->atoms.nr;

    at2con = make_at2con(0, natoms, molt->ilist, iparams,
                         EI_DYNAMICS(ir->eI), &nflexcon);
    snew(path, 1+ir->nProjOrder);
    for (at = 0; at < 1+ir->nProjOrder; at++)
    {
        path[at] = -1;
    }

    r2maxA = 0;
    for (at = 0; at < natoms; at++)
    {
        r0 = 0;
        r1 = 0;

        count = 0;
        constr_recur(&at2con, molt->ilist, iparams,
                     FALSE, at, 0, 1+ir->nProjOrder, path, r0, r1, &r2maxA, &count);
    }
    if (ir->efep == efepNO)
    {
        rmax = sqrt(r2maxA);
    }
    else
    {
        r2maxB = 0;
        for (at = 0; at < natoms; at++)
        {
            r0    = 0;
            r1    = 0;
            count = 0;
            constr_recur(&at2con, molt->ilist, iparams,
                         TRUE, at, 0, 1+ir->nProjOrder, path, r0, r1, &r2maxB, &count);
        }
        lam0 = ir->fepvals->init_lambda;
        if (EI_DYNAMICS(ir->eI))
        {
            lam0 += ir->init_step*ir->fepvals->delta_lambda;
        }
        rmax = (1 - lam0)*sqrt(r2maxA) + lam0*sqrt(r2maxB);
        if (EI_DYNAMICS(ir->eI))
        {
            lam1 = ir->fepvals->init_lambda + (ir->init_step + ir->nsteps)*ir->fepvals->delta_lambda;
            rmax = std::max(rmax, (1 - lam1)*std::sqrt(r2maxA) + lam1*std::sqrt(r2maxB));
        }
    }

    done_blocka(&at2con);
    sfree(path);

    return rmax;
}

real constr_r_max(FILE *fplog, gmx_mtop_t *mtop, t_inputrec *ir)
{
    int  mt;
    real rmax;

    rmax = 0;
    for (mt = 0; mt < mtop->nmoltype; mt++)
    {
        rmax = std::max(rmax,
                        constr_r_max_moltype(&mtop->moltype[mt],
                                             mtop->ffparams.iparams, ir));
    }

    if (fplog)
    {
        fprintf(fplog, "Maximum distance for %d constraints, at 120 deg. angles, all-trans: %.3f nm\n", 1+ir->nProjOrder, rmax);
    }

    return rmax;
}

gmx_constr_t init_constraints(FILE *fplog,
                              gmx_mtop_t *mtop, t_inputrec *ir,
                              gmx_edsam_t ed, t_state *state,
                              t_commrec *cr)
{
    int                  ncon, nset, nmol, settle_type, i, mt, nflexcon;
    struct gmx_constr   *constr;
    char                *env;
    t_ilist             *ilist;
    gmx_mtop_ilistloop_t iloop;

    ncon =
        gmx_mtop_ftype_count(mtop, F_CONSTR) +
        gmx_mtop_ftype_count(mtop, F_CONSTRNC);
    nset = gmx_mtop_ftype_count(mtop, F_SETTLE);

    if (ncon+nset == 0 &&
        !(ir->bPull && pull_have_constraint(ir->pull_work)) &&
        ed == NULL)
    {
        return NULL;
    }

    snew(constr, 1);

    constr->ncon_tot = ncon;
    constr->nflexcon = 0;
    if (ncon > 0)
    {
        constr->n_at2con_mt = mtop->nmoltype;
        snew(constr->at2con_mt, constr->n_at2con_mt);
        for (mt = 0; mt < mtop->nmoltype; mt++)
        {
            constr->at2con_mt[mt] = make_at2con(0, mtop->moltype[mt].atoms.nr,
                                                mtop->moltype[mt].ilist,
                                                mtop->ffparams.iparams,
                                                EI_DYNAMICS(ir->eI), &nflexcon);
            for (i = 0; i < mtop->nmolblock; i++)
            {
                if (mtop->molblock[i].type == mt)
                {
                    constr->nflexcon += mtop->molblock[i].nmol*nflexcon;
                }
            }
        }

        if (constr->nflexcon > 0)
        {
            if (fplog)
            {
                fprintf(fplog, "There are %d flexible constraints\n",
                        constr->nflexcon);
                if (ir->fc_stepsize == 0)
                {
                    fprintf(fplog, "\n"
                            "WARNING: step size for flexible constraining = 0\n"
                            "         All flexible constraints will be rigid.\n"
                            "         Will try to keep all flexible constraints at their original length,\n"
                            "         but the lengths may exhibit some drift.\n\n");
                    constr->nflexcon = 0;
                }
            }
            if (constr->nflexcon > 0)
            {
                please_cite(fplog, "Hess2002");
            }
        }

        if (ir->eConstrAlg == econtLINCS)
        {
            constr->lincsd = init_lincs(fplog, mtop,
                                        constr->nflexcon, constr->at2con_mt,
                                        DOMAINDECOMP(cr) && cr->dd->bInterCGcons,
                                        ir->nLincsIter, ir->nProjOrder);
        }

        if (ir->eConstrAlg == econtSHAKE)
        {
            if (DOMAINDECOMP(cr) && cr->dd->bInterCGcons)
            {
                gmx_fatal(FARGS, "SHAKE is not supported with domain decomposition and constraint that cross charge group boundaries, use LINCS");
            }
            if (constr->nflexcon)
            {
                gmx_fatal(FARGS, "For this system also velocities and/or forces need to be constrained, this can not be done with SHAKE, you should select LINCS");
            }
            please_cite(fplog, "Ryckaert77a");
            if (ir->bShakeSOR)
            {
                please_cite(fplog, "Barth95a");
            }

            constr->shaked = shake_init();
        }
    }

    if (nset > 0)
    {
        please_cite(fplog, "Miyamoto92a");

        constr->bInterCGsettles = inter_charge_group_settles(mtop);

        /* Check that we have only one settle type */
        settle_type = -1;
        iloop       = gmx_mtop_ilistloop_init(mtop);
        while (gmx_mtop_ilistloop_next(iloop, &ilist, &nmol))
        {
            for (i = 0; i < ilist[F_SETTLE].nr; i += 4)
            {
                if (settle_type == -1)
                {
                    settle_type = ilist[F_SETTLE].iatoms[i];
                }
                else if (ilist[F_SETTLE].iatoms[i] != settle_type)
                {
                    gmx_fatal(FARGS,
                              "The [molecules] section of your topology specifies more than one block of\n"
                              "a [moleculetype] with a [settles] block. Only one such is allowed. If you\n"
                              "are trying to partition your solvent into different *groups* (e.g. for\n"
                              "freezing, T-coupling, etc.) then you are using the wrong approach. Index\n"
                              "files specify groups. Otherwise, you may wish to change the least-used\n"
                              "block of molecules with SETTLE constraints into 3 normal constraints.");
                }
            }
        }

        constr->n_at2settle_mt = mtop->nmoltype;
        snew(constr->at2settle_mt, constr->n_at2settle_mt);
        for (mt = 0; mt < mtop->nmoltype; mt++)
        {
            constr->at2settle_mt[mt] =
                make_at2settle(mtop->moltype[mt].atoms.nr,
                               &mtop->moltype[mt].ilist[F_SETTLE]);
        }
    }

    if ((ncon + nset) > 0 && ir->epc == epcMTTK)
    {
        gmx_fatal(FARGS, "Constraints are not implemented with MTTK pressure control.");
    }

    constr->maxwarn = 999;
    env             = getenv("GMX_MAXCONSTRWARN");
    if (env)
    {
        constr->maxwarn = 0;
        sscanf(env, "%8d", &constr->maxwarn);
        if (fplog)
        {
            fprintf(fplog,
                    "Setting the maximum number of constraint warnings to %d\n",
                    constr->maxwarn);
        }
        if (MASTER(cr))
        {
            fprintf(stderr,
                    "Setting the maximum number of constraint warnings to %d\n",
                    constr->maxwarn);
        }
    }
    if (constr->maxwarn < 0 && fplog)
    {
        fprintf(fplog, "maxwarn < 0, will not stop on constraint errors\n");
    }
    constr->warncount_lincs  = 0;
    constr->warncount_settle = 0;

    /* Initialize the essential dynamics sampling.
     * Put the pointer to the ED struct in constr */
    constr->ed = ed;
    if (ed != NULL || state->edsamstate.nED > 0)
    {
        init_edsam(mtop, ir, cr, ed, state->x, state->box, &state->edsamstate);
    }

    constr->warn_mtop = mtop;

    return constr;
}

const t_blocka *atom2constraints_moltype(gmx_constr_t constr)
{
    return constr->at2con_mt;
}

const int **atom2settle_moltype(gmx_constr_t constr)
{
    return (const int **)constr->at2settle_mt;
}


gmx_bool inter_charge_group_constraints(const gmx_mtop_t *mtop)
{
    const gmx_moltype_t *molt;
    const t_block       *cgs;
    const t_ilist       *il;
    int                  mb;
    int                 *at2cg, cg, a, ftype, i;
    gmx_bool             bInterCG;

    bInterCG = FALSE;
    for (mb = 0; mb < mtop->nmolblock && !bInterCG; mb++)
    {
        molt = &mtop->moltype[mtop->molblock[mb].type];

        if (molt->ilist[F_CONSTR].nr   > 0 ||
            molt->ilist[F_CONSTRNC].nr > 0 ||
            molt->ilist[F_SETTLE].nr > 0)
        {
            cgs  = &molt->cgs;
            snew(at2cg, molt->atoms.nr);
            for (cg = 0; cg < cgs->nr; cg++)
            {
                for (a = cgs->index[cg]; a < cgs->index[cg+1]; a++)
                {
                    at2cg[a] = cg;
                }
            }

            for (ftype = F_CONSTR; ftype <= F_CONSTRNC; ftype++)
            {
                il = &molt->ilist[ftype];
                for (i = 0; i < il->nr && !bInterCG; i += 1+NRAL(ftype))
                {
                    if (at2cg[il->iatoms[i+1]] != at2cg[il->iatoms[i+2]])
                    {
                        bInterCG = TRUE;
                    }
                }
            }

            sfree(at2cg);
        }
    }

    return bInterCG;
}

gmx_bool inter_charge_group_settles(const gmx_mtop_t *mtop)
{
    const gmx_moltype_t *molt;
    const t_block       *cgs;
    const t_ilist       *il;
    int                  mb;
    int                 *at2cg, cg, a, ftype, i;
    gmx_bool             bInterCG;

    bInterCG = FALSE;
    for (mb = 0; mb < mtop->nmolblock && !bInterCG; mb++)
    {
        molt = &mtop->moltype[mtop->molblock[mb].type];

        if (molt->ilist[F_SETTLE].nr > 0)
        {
            cgs  = &molt->cgs;
            snew(at2cg, molt->atoms.nr);
            for (cg = 0; cg < cgs->nr; cg++)
            {
                for (a = cgs->index[cg]; a < cgs->index[cg+1]; a++)
                {
                    at2cg[a] = cg;
                }
            }

            for (ftype = F_SETTLE; ftype <= F_SETTLE; ftype++)
            {
                il = &molt->ilist[ftype];
                for (i = 0; i < il->nr && !bInterCG; i += 1+NRAL(F_SETTLE))
                {
                    if (at2cg[il->iatoms[i+1]] != at2cg[il->iatoms[i+2]] ||
                        at2cg[il->iatoms[i+1]] != at2cg[il->iatoms[i+3]])
                    {
                        bInterCG = TRUE;
                    }
                }
            }

            sfree(at2cg);
        }
    }

    return bInterCG;
}

/* helper functions for andersen temperature control, because the
 * gmx_constr construct is only defined in constr.c. Return the list
 * of blocks (get_sblock) and the number of blocks (get_nblocks).  */

extern int *get_sblock(struct gmx_constr *constr)
{
    return constr->sblock;
}

extern int get_nblocks(struct gmx_constr *constr)
{
    return constr->nblocks;
}

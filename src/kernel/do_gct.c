/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "xmdrun.h"
#include "futil.h"
#include "xvgr.h"
#include "macros.h"
#include "physics.h"
#include "network.h"
#include "smalloc.h"
#include "mdrun.h"
#include "string2.h"
#include "readinp.h"
#include "filenm.h"
#include "update.h"
#include "vec.h"
#include "main.h"
#include "txtdump.h"

/*#define DEBUGGCT*/
t_coupl_rec *init_coupling(FILE *log, int nfile, const t_filenm fnm[],
                           t_commrec *cr, t_forcerec *fr,
                           t_mdatoms *md, t_idef *idef)
{
    int          i, nc, index, j;
    int          ati, atj;
    t_coupl_rec *tcr;

    snew(tcr, 1);
    if (MASTER(cr))
    {
        read_gct (opt2fn("-j", nfile, fnm), tcr);
        write_gct(opt2fn("-jo", nfile, fnm), tcr, idef);
    }
    /* Update all processors with coupling info */
    if (PAR(cr))
    {
        comm_tcr(log, cr, &tcr);
    }

    /* Copy information from the coupling to the force field stuff */
    copy_ff(tcr, fr, md, idef);

    return tcr;
}

static real Ecouple(t_coupl_rec *tcr, real ener[])
{
    if (tcr->bInter)
    {
        return ener[F_COUL_SR]+ener[F_LJ]+ener[F_COUL_LR]+ener[F_LJ_LR]+
               ener[F_RF_EXCL]+ener[F_COUL_RECIP];
    }
    else
    {
        return ener[F_EPOT];
    }
}

static char *mk_gct_nm(const char *fn, int ftp, int ati, int atj)
{
    static char buf[256];

    strcpy(buf, fn);
    if (atj == -1)
    {
        sprintf(buf+strlen(fn)-4, "%d.%s", ati, ftp2ext(ftp));
    }
    else
    {
        sprintf(buf+strlen(fn)-4, "%d_%d.%s", ati, atj, ftp2ext(ftp));
    }

    return buf;
}

static void pr_ff(t_coupl_rec *tcr, real time, t_idef *idef,
                  t_commrec *cr, int nfile, const t_filenm fnm[],
                  const output_env_t oenv)
{
    static FILE  *prop;
    static FILE **out = NULL;
    static FILE **qq  = NULL;
    static FILE **ip  = NULL;
    t_coupl_LJ   *tclj;
    t_coupl_BU   *tcbu;
    char          buf[256];
    const char   *leg[]  =  { "C12", "C6" };
    const char   *eleg[] =  { "Epsilon", "Sigma" };
    const char   *bleg[] = { "A", "B", "C" };
    char        **raleg;
    int           i, j, index;

    if ((prop == NULL) && (out == NULL) && (qq == NULL) && (ip == NULL))
    {
        prop = xvgropen(opt2fn("-runav", nfile, fnm),
                        "Properties and Running Averages", "Time (ps)", "", oenv);
        snew(raleg, 2*eoObsNR);
        for (i = j = 0; (i < eoObsNR); i++)
        {
            if (tcr->bObsUsed[i])
            {
                raleg[j++] = strdup(eoNames[i]);
                sprintf(buf, "RA-%s", eoNames[i]);
                raleg[j++] = strdup(buf);
            }
        }
        xvgr_legend(prop, j, (const char**)raleg, oenv);
        for (i = 0; (i < j); i++)
        {
            sfree(raleg[i]);
        }
        sfree(raleg);

        if (tcr->nLJ)
        {
            snew(out, tcr->nLJ);
            for (i = 0; (i < tcr->nLJ); i++)
            {
                if (tcr->tcLJ[i].bPrint)
                {
                    tclj   = &(tcr->tcLJ[i]);
                    out[i] =
                        xvgropen(mk_gct_nm(opt2fn("-ffout", nfile, fnm),
                                           efXVG, tclj->at_i, tclj->at_j),
                                 "General Coupling Lennard Jones", "Time (ps)",
                                 "Force constant (units)", oenv);
                    fprintf(out[i], "@ subtitle \"Interaction between types %d and %d\"\n",
                            tclj->at_i, tclj->at_j);
                    if (tcr->combrule == 1)
                    {
                        xvgr_legend(out[i], asize(leg), leg, oenv);
                    }
                    else
                    {
                        xvgr_legend(out[i], asize(eleg), eleg, oenv);
                    }
                    fflush(out[i]);
                }
            }
        }
        else if (tcr->nBU)
        {
            snew(out, tcr->nBU);
            for (i = 0; (i < tcr->nBU); i++)
            {
                if (tcr->tcBU[i].bPrint)
                {
                    tcbu   = &(tcr->tcBU[i]);
                    out[i] =
                        xvgropen(mk_gct_nm(opt2fn("-ffout", nfile, fnm), efXVG,
                                           tcbu->at_i, tcbu->at_j),
                                 "General Coupling Buckingham", "Time (ps)",
                                 "Force constant (units)", oenv);
                    fprintf(out[i], "@ subtitle \"Interaction between types %d and %d\"\n",
                            tcbu->at_i, tcbu->at_j);
                    xvgr_legend(out[i], asize(bleg), bleg, oenv);
                    fflush(out[i]);
                }
            }
        }
        snew(qq, tcr->nQ);
        for (i = 0; (i < tcr->nQ); i++)
        {
            if (tcr->tcQ[i].bPrint)
            {
                qq[i] = xvgropen(mk_gct_nm(opt2fn("-ffout", nfile, fnm), efXVG,
                                           tcr->tcQ[i].at_i, -1),
                                 "General Coupling Charge", "Time (ps)", "Charge (e)",
                                 oenv);
                fprintf(qq[i], "@ subtitle \"Type %d\"\n", tcr->tcQ[i].at_i);
                fflush(qq[i]);
            }
        }
        snew(ip, tcr->nIP);
        for (i = 0; (i < tcr->nIP); i++)
        {
            sprintf(buf, "gctIP%d", tcr->tIP[i].type);
            ip[i] = xvgropen(mk_gct_nm(opt2fn("-ffout", nfile, fnm), efXVG, 0, -1),
                             "General Coupling iparams", "Time (ps)", "ip ()", oenv);
            index = tcr->tIP[i].type;
            fprintf(ip[i], "@ subtitle \"Coupling to %s\"\n",
                    interaction_function[idef->functype[index]].longname);
            fflush(ip[i]);
        }
    }
    /* Write properties to file */
    fprintf(prop, "%10.3f", time);
    for (i = 0; (i < eoObsNR); i++)
    {
        if (tcr->bObsUsed[i])
        {
            fprintf(prop, "  %10.3e  %10.3e", tcr->act_value[i], tcr->av_value[i]);
        }
    }
    fprintf(prop, "\n");
    fflush(prop);

    for (i = 0; (i < tcr->nLJ); i++)
    {
        tclj = &(tcr->tcLJ[i]);
        if (tclj->bPrint)
        {
            if (tcr->combrule == 1)
            {
                fprintf(out[i], "%14.7e  %14.7e  %14.7e\n",
                        time, tclj->c12, tclj->c6);
            }
            else
            {
                double sigma   = pow(tclj->c12/tclj->c6, 1.0/6.0);
                double epsilon = 0.25*sqr(tclj->c6)/tclj->c12;
                fprintf(out[i], "%14.7e  %14.7e  %14.7e\n",
                        time, epsilon, sigma);
            }
            fflush(out[i]);
        }
    }
    for (i = 0; (i < tcr->nBU); i++)
    {
        tcbu = &(tcr->tcBU[i]);
        if (tcbu->bPrint)
        {
            fprintf(out[i], "%14.7e  %14.7e  %14.7e  %14.7e\n",
                    time, tcbu->a, tcbu->b, tcbu->c);
            fflush(out[i]);
        }
    }
    for (i = 0; (i < tcr->nQ); i++)
    {
        if (tcr->tcQ[i].bPrint)
        {
            fprintf(qq[i], "%14.7e  %14.7e\n", time, tcr->tcQ[i].Q);
            fflush(qq[i]);
        }
    }
    for (i = 0; (i < tcr->nIP); i++)
    {
        fprintf(ip[i], "%10g  ", time);
        index = tcr->tIP[i].type;
        switch (idef->functype[index])
        {
            case F_BONDS:
                fprintf(ip[i], "%10g  %10g\n", tcr->tIP[i].iprint.harmonic.krA,
                        tcr->tIP[i].iprint.harmonic.rA);
                break;
            default:
                break;
        }
        fflush(ip[i]);
    }
}

static void pr_dev(t_coupl_rec *tcr,
                   real t, real dev[eoObsNR], t_commrec *cr, int nfile,
                   const t_filenm fnm[], const output_env_t oenv)
{
    static FILE *fp = NULL;
    char       **ptr;
    int          i, j;

    if (!fp)
    {
        fp = xvgropen(opt2fn("-devout", nfile, fnm),
                      "Deviations from target value", "Time (ps)", "", oenv);
        snew(ptr, eoObsNR);
        for (i = j = 0; (i < eoObsNR); i++)
        {
            if (tcr->bObsUsed[i])
            {
                ptr[j++] = strdup(eoNames[i]);
            }
        }
        xvgr_legend(fp, j, (const char**)ptr, oenv);
        for (i = 0; i < j; i++)
        {
            sfree(ptr[i]);
        }
        sfree(ptr);
    }
    fprintf(fp, "%10.3f", t);
    for (i = 0; (i < eoObsNR); i++)
    {
        if (tcr->bObsUsed[i])
        {
            fprintf(fp, "  %10.3e", dev[i]);
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
}

static void upd_nbfplj(FILE *log, real *nbfp, int atnr, real f6[], real f12[],
                       int combrule)
{
    double *sigma, *epsilon, c6, c12, eps, sig, sig6;
    int     n, m, k;

    /* Update the nonbonded force parameters */
    switch (combrule)
    {
        case 1:
            for (k = n = 0; (n < atnr); n++)
            {
                for (m = 0; (m < atnr); m++, k++)
                {
                    C6 (nbfp, atnr, n, m) *= f6[k];
                    C12(nbfp, atnr, n, m) *= f12[k];
                }
            }
            break;
        case 2:
        case 3:
            /* Convert to sigma and epsilon */
            snew(sigma, atnr);
            snew(epsilon, atnr);
            for (n = 0; (n < atnr); n++)
            {
                k   = n*(atnr+1);
                c6  = C6 (nbfp, atnr, n, n) * f6[k];
                c12 = C12(nbfp, atnr, n, n) * f12[k];
                if ((c6 == 0) || (c12 == 0))
                {
                    gmx_fatal(FARGS, "You can not use combination rule %d with zero C6 (%f) or C12 (%f)", combrule, c6, c12);
                }
                sigma[n]   = pow(c12/c6, 1.0/6.0);
                epsilon[n] = 0.25*(c6*c6/c12);
            }
            for (k = n = 0; (n < atnr); n++)
            {
                for (m = 0; (m < atnr); m++, k++)
                {
                    eps  = sqrt(epsilon[n]*epsilon[m]);
                    if (combrule == 2)
                    {
                        sig  = 0.5*(sigma[n]+sigma[m]);
                    }
                    else
                    {
                        sig  = sqrt(sigma[n]*sigma[m]);
                    }
                    sig6 = pow(sig, 6.0);
                    /* nbfp now includes the 6.0/12.0 derivative prefactors */
                    C6 (nbfp, atnr, n, m) = 4*eps*sig6/6.0;
                    C12(nbfp, atnr, n, m) = 4*eps*sig6*sig6/12.0;
                }
            }
            sfree(sigma);
            sfree(epsilon);
            break;
        default:
            gmx_fatal(FARGS, "Combination rule should be 1,2 or 3 instead of %d",
                      combrule);
    }
}

static void upd_nbfpbu(FILE *log, real *nbfp, int atnr,
                       real fa[], real fb[], real fc[])
{
    int n, m, k;

    /* Update the nonbonded force parameters */
    for (k = n = 0; (n < atnr); n++)
    {
        for (m = 0; (m < atnr); m++, k++)
        {
            BHAMA(nbfp, atnr, n, m) *= fa[k];
            BHAMB(nbfp, atnr, n, m) *= fb[k];
            BHAMC(nbfp, atnr, n, m) *= fc[k];
        }
    }
}

void gprod(t_commrec *cr, int n, real f[])
{
    /* Compute the global product of all elements in an array
     * such that after gprod f[i] = PROD_j=1,nprocs f[i][j]
     */
    static int   nbuf = 0;
    static real *buf  = NULL;
    int          i;

    if (n > nbuf)
    {
        nbuf = n;
        srenew(buf, nbuf);
    }

#ifdef GMX_MPI
#ifdef GMX_DOUBLE
    MPI_Allreduce(f, buf, n, MPI_DOUBLE, MPI_PROD, cr->mpi_comm_mygroup);
#else
    MPI_Allreduce(f, buf, n, MPI_FLOAT, MPI_PROD, cr->mpi_comm_mygroup);
#endif
    for (i = 0; (i < n); i++)
    {
        f[i] = buf[i];
    }
#endif
}

static void set_factor_matrix(int ntypes, real f[], real fmult, int ati, int atj)
{
#define FMIN 0.95
#define FMAX 1.05
    int i;

    fmult = min(FMAX, max(FMIN, fmult));
    if (atj != -1)
    {
        f[ntypes*ati+atj] *= fmult;
        f[ntypes*atj+ati] *= fmult;
    }
    else
    {
        for (i = 0; (i < ntypes); i++)
        {
            f[ntypes*ati+i] *= fmult;
            f[ntypes*i+ati] *= fmult;
        }
    }
#undef FMIN
#undef FMAX
}

static real calc_deviation(real xav, real xt, real x0)
{
    /* This may prevent overshooting in GCT coupling... */

    /* real dev;

       if (xav > x0) {
       if (xt > x0)
        dev = min(xav-x0,xt-x0);
       else
        dev = 0;
       }
       else {
       if (xt < x0)
        dev = max(xav-x0,xt-x0);
       else
        dev = 0;
       }
     */
    return x0-xav;
}

static real calc_dist(FILE *log, rvec x[])
{
    static gmx_bool bFirst = TRUE;
    static gmx_bool bDist;
    static int      i1, i2;
    char           *buf;
    rvec            dx;

    if (bFirst)
    {
        if ((buf = getenv("DISTGCT")) == NULL)
        {
            bDist = FALSE;
        }
        else
        {
            bDist  = (sscanf(buf, "%d%d", &i1, &i2) == 2);
            if (bDist)
            {
                fprintf(log, "GCT: Will couple to distance between %d and %d\n", i1, i2);
            }
            else
            {
                fprintf(log, "GCT: Will not couple to distances\n");
            }
        }
        bFirst = FALSE;
    }
    if (bDist)
    {
        rvec_sub(x[i1], x[i2], dx);
        return norm(dx);
    }
    else
    {
        return 0.0;
    }
}

real run_aver(real old, real cur, int step, int nmem)
{
    nmem   = max(1, nmem);

    return ((nmem-1)*old+cur)/nmem;
}

static void set_act_value(t_coupl_rec *tcr, int index, real val, int step)
{
    tcr->act_value[index] = val;
    tcr->av_value[index]  = run_aver(tcr->av_value[index], val, step, tcr->nmemory);
}

static void upd_f_value(FILE *log, int atnr, real xi, real dt, real factor,
                        real ff[], int ati, int atj)
{
    real fff;

    if (xi != 0)
    {
        fff = 1 + (dt/xi)  * factor;
        if (fff > 0)
        {
            set_factor_matrix(atnr, ff, sqrt(fff), ati, atj);
        }
    }
}

static void dump_fm(FILE *fp, int n, real f[], char *s)
{
    int i, j;

    fprintf(fp, "Factor matrix (all numbers -1) %s\n", s);
    for (i = 0; (i < n); i++)
    {
        for (j = 0; (j < n); j++)
        {
            fprintf(fp, "  %10.3e", f[n*i+j]-1.0);
        }
        fprintf(fp, "\n");
    }
}

void do_coupling(FILE *log, const output_env_t oenv, int nfile,
                 const t_filenm fnm[], t_coupl_rec *tcr, real t,
                 int step, real ener[], t_forcerec *fr, t_inputrec *ir,
                 gmx_bool bMaster, t_mdatoms *md, t_idef *idef, real mu_aver, int nmols,
                 t_commrec *cr, matrix box, tensor virial,
                 tensor pres, rvec mu_tot,
                 rvec x[], rvec f[], gmx_bool bDoIt)
{
#define enm2Debye 48.0321
#define d2e(x) (x)/enm2Debye
#define enm2kjmol(x) (x)*0.0143952 /* = 2.0*4.0*M_PI*EPSILON0 */

    static real     *f6, *f12, *fa, *fb, *fc, *fq;
    static gmx_bool  bFirst = TRUE;

    int              i, j, ati, atj, atnr2, type, ftype;
    real             deviation[eoObsNR], prdev[eoObsNR], epot0, dist, rmsf;
    real             ff6, ff12, ffa, ffb, ffc, ffq, factor, dt, mu_ind;
    real             Epol, Eintern, Virial, muabs, xiH = -1, xiS = -1, xi6, xi12;
    rvec             fmol[2];
    gmx_bool         bTest, bPrint;
    t_coupl_LJ      *tclj;
    t_coupl_BU      *tcbu;
    t_coupl_Q       *tcq;
    t_coupl_iparams *tip;

    atnr2 = idef->atnr * idef->atnr;
    if (bFirst)
    {
        if (PAR(cr))
        {
            fprintf(log, "GCT: this is parallel\n");
        }
        else
        {
            fprintf(log, "GCT: this is not parallel\n");
        }
        fflush(log);
        snew(f6, atnr2);
        snew(f12, atnr2);
        snew(fa, atnr2);
        snew(fb, atnr2);
        snew(fc, atnr2);
        snew(fq, idef->atnr);

        if (tcr->bVirial)
        {
            int  nrdf = 0;
            real TTT  = 0;
            real Vol  = det(box);

            for (i = 0; (i < ir->opts.ngtc); i++)
            {
                nrdf += ir->opts.nrdf[i];
                TTT  += ir->opts.nrdf[i]*ir->opts.ref_t[i];
            }
            TTT /= nrdf;

            /* Calculate reference virial from reference temperature and pressure */
            tcr->ref_value[eoVir] = 0.5*BOLTZ*nrdf*TTT - (3.0/2.0)*
                Vol*tcr->ref_value[eoPres];

            fprintf(log, "GCT: TTT = %g, nrdf = %d, vir0 = %g,  Vol = %g\n",
                    TTT, nrdf, tcr->ref_value[eoVir], Vol);
            fflush(log);
        }
        bFirst = FALSE;
    }

    bPrint = MASTER(cr) && do_per_step(step, ir->nstlog);
    dt     = ir->delta_t;

    /* Initiate coupling to the reference pressure and temperature to start
     * coupling slowly.
     */
    if (step == 0)
    {
        for (i = 0; (i < eoObsNR); i++)
        {
            tcr->av_value[i] = tcr->ref_value[i];
        }
        if ((tcr->ref_value[eoDipole]) != 0.0)
        {
            mu_ind                 = mu_aver - d2e(tcr->ref_value[eoDipole]); /* in e nm */
            Epol                   = mu_ind*mu_ind/(enm2kjmol(tcr->ref_value[eoPolarizability]));
            tcr->av_value[eoEpot] -= Epol;
            fprintf(log, "GCT: mu_aver = %g(D), mu_ind = %g(D), Epol = %g (kJ/mol)\n",
                    mu_aver*enm2Debye, mu_ind*enm2Debye, Epol);
        }
    }

    /* We want to optimize the LJ params, usually to the Vaporization energy
     * therefore we only count intermolecular degrees of freedom.
     * Note that this is now optional. switch UseEinter to yes in your gct file
     * if you want this.
     */
    dist      = calc_dist(log, x);
    muabs     = norm(mu_tot);
    Eintern   = Ecouple(tcr, ener)/nmols;
    Virial    = virial[XX][XX]+virial[YY][YY]+virial[ZZ][ZZ];

    /*calc_force(md->nr,f,fmol);*/
    clear_rvec(fmol[0]);

    /* Use a memory of tcr->nmemory steps, so we actually couple to the
     * average observable over the last tcr->nmemory steps. This may help
     * in avoiding local minima in parameter space.
     */
    set_act_value(tcr, eoPres, ener[F_PRES], step);
    set_act_value(tcr, eoEpot, Eintern,     step);
    set_act_value(tcr, eoVir,  Virial,      step);
    set_act_value(tcr, eoDist, dist,        step);
    set_act_value(tcr, eoMu,   muabs,       step);
    set_act_value(tcr, eoFx,   fmol[0][XX], step);
    set_act_value(tcr, eoFy,   fmol[0][YY], step);
    set_act_value(tcr, eoFz,   fmol[0][ZZ], step);
    set_act_value(tcr, eoPx,   pres[XX][XX], step);
    set_act_value(tcr, eoPy,   pres[YY][YY], step);
    set_act_value(tcr, eoPz,   pres[ZZ][ZZ], step);

    epot0 = tcr->ref_value[eoEpot];
    /* If dipole != 0.0 assume we want to use polarization corrected coupling */
    if ((tcr->ref_value[eoDipole]) != 0.0)
    {
        mu_ind = mu_aver - d2e(tcr->ref_value[eoDipole]); /* in e nm */

        Epol = mu_ind*mu_ind/(enm2kjmol(tcr->ref_value[eoPolarizability]));

        epot0 -= Epol;
        if (debug)
        {
            fprintf(debug, "mu_ind: %g (%g D) mu_aver: %g (%g D)\n",
                    mu_ind, mu_ind*enm2Debye, mu_aver, mu_aver*enm2Debye);
            fprintf(debug, "Eref %g Epol %g Erunav %g Eact %g\n",
                    tcr->ref_value[eoEpot], Epol, tcr->av_value[eoEpot],
                    tcr->act_value[eoEpot]);
        }
    }

    if (bPrint)
    {
        pr_ff(tcr, t, idef, cr, nfile, fnm, oenv);
    }
    /* Calculate the deviation of average value from the target value */
    for (i = 0; (i < eoObsNR); i++)
    {
        deviation[i] = calc_deviation(tcr->av_value[i], tcr->act_value[i],
                                      tcr->ref_value[i]);
        prdev[i]     = tcr->ref_value[i] - tcr->act_value[i];
    }
    deviation[eoEpot] = calc_deviation(tcr->av_value[eoEpot], tcr->act_value[eoEpot],
                                       epot0);
    prdev[eoEpot]     = epot0 - tcr->act_value[eoEpot];

    if (bPrint)
    {
        pr_dev(tcr, t, prdev, cr, nfile, fnm, oenv);
    }

    /* First set all factors to 1 */
    for (i = 0; (i < atnr2); i++)
    {
        f6[i] = f12[i] = fa[i] = fb[i] = fc[i] = 1.0;
    }
    for (i = 0; (i < idef->atnr); i++)
    {
        fq[i] = 1.0;
    }

    /* Now compute the actual coupling compononents */
    if (!fr->bBHAM)
    {
        if (bDoIt)
        {
            for (i = 0; (i < tcr->nLJ); i++)
            {
                tclj = &(tcr->tcLJ[i]);

                ati = tclj->at_i;
                atj = tclj->at_j;

                ff6 = ff12 = 1.0;

                if (tclj->eObs == eoForce)
                {
                    gmx_fatal(FARGS, "Hack code for this to work again ");
                    if (debug)
                    {
                        fprintf(debug, "Have computed derivatives: xiH = %g, xiS = %g\n", xiH, xiS);
                    }
                    if (ati == 1)
                    {
                        /* Hydrogen */
                        ff12 += xiH;
                    }
                    else if (ati == 2)
                    {
                        /* Shell */
                        ff12 += xiS;
                    }
                    else
                    {
                        gmx_fatal(FARGS, "No H, no Shell, edit code at %s, line %d\n",
                                  __FILE__, __LINE__);
                    }
                    if (ff6 > 0)
                    {
                        set_factor_matrix(idef->atnr, f6, sqrt(ff6), ati, atj);
                    }
                    if (ff12 > 0)
                    {
                        set_factor_matrix(idef->atnr, f12, sqrt(ff12), ati, atj);
                    }
                }
                else
                {
                    if (debug)
                    {
                        fprintf(debug, "tcr->tcLJ[%d].xi_6 = %g, xi_12 = %g deviation = %g\n", i,
                                tclj->xi_6, tclj->xi_12, deviation[tclj->eObs]);
                    }
                    factor = deviation[tclj->eObs];

                    upd_f_value(log, idef->atnr, tclj->xi_6, dt, factor, f6, ati, atj);
                    upd_f_value(log, idef->atnr, tclj->xi_12, dt, factor, f12, ati, atj);
                }
            }
        }
        if (PAR(cr))
        {
            gprod(cr, atnr2, f6);
            gprod(cr, atnr2, f12);
#ifdef DEBUGGCT
            dump_fm(log, idef->atnr, f6, "f6");
            dump_fm(log, idef->atnr, f12, "f12");
#endif
        }
        upd_nbfplj(log, fr->nbfp, idef->atnr, f6, f12, tcr->combrule);

        /* Copy for printing */
        for (i = 0; (i < tcr->nLJ); i++)
        {
            tclj = &(tcr->tcLJ[i]);
            ati  = tclj->at_i;
            atj  = tclj->at_j;
            if (atj == -1)
            {
                atj = ati;
            }
            /* nbfp now includes the 6.0/12.0 derivative prefactors */
            tclj->c6  =  C6(fr->nbfp, fr->ntype, ati, atj)/6.0;
            tclj->c12 = C12(fr->nbfp, fr->ntype, ati, atj)/12.0;
        }
    }
    else
    {
        if (bDoIt)
        {
            for (i = 0; (i < tcr->nBU); i++)
            {
                tcbu   = &(tcr->tcBU[i]);
                factor = deviation[tcbu->eObs];
                ati    = tcbu->at_i;
                atj    = tcbu->at_j;

                upd_f_value(log, idef->atnr, tcbu->xi_a, dt, factor, fa, ati, atj);
                upd_f_value(log, idef->atnr, tcbu->xi_b, dt, factor, fb, ati, atj);
                upd_f_value(log, idef->atnr, tcbu->xi_c, dt, factor, fc, ati, atj);
            }
        }
        if (PAR(cr))
        {
            gprod(cr, atnr2, fa);
            gprod(cr, atnr2, fb);
            gprod(cr, atnr2, fc);
        }
        upd_nbfpbu(log, fr->nbfp, idef->atnr, fa, fb, fc);
        /* Copy for printing */
        for (i = 0; (i < tcr->nBU); i++)
        {
            tcbu = &(tcr->tcBU[i]);
            ati  = tcbu->at_i;
            atj  = tcbu->at_j;
            if (atj == -1)
            {
                atj = ati;
            }
            /* nbfp now includes the 6.0 derivative prefactors */
            tcbu->a = BHAMA(fr->nbfp, fr->ntype, ati, atj);
            tcbu->b = BHAMB(fr->nbfp, fr->ntype, ati, atj);
            tcbu->c = BHAMC(fr->nbfp, fr->ntype, ati, atj)/6.0;
            if (debug)
            {
                fprintf(debug, "buck (type=%d) = %e, %e, %e\n",
                        tcbu->at_i, tcbu->a, tcbu->b, tcbu->c);
            }
        }
    }
    if (bDoIt)
    {
        for (i = 0; (i < tcr->nQ); i++)
        {
            tcq = &(tcr->tcQ[i]);
            if (tcq->xi_Q)
            {
                ffq = 1.0 + (dt/tcq->xi_Q) * deviation[tcq->eObs];
            }
            else
            {
                ffq = 1.0;
            }
            fq[tcq->at_i] *= ffq;
        }
    }
    if (PAR(cr))
    {
        gprod(cr, idef->atnr, fq);
    }

    for (j = 0; (j < md->nr); j++)
    {
        md->chargeA[j] *= fq[md->typeA[j]];
    }
    for (i = 0; (i < tcr->nQ); i++)
    {
        tcq = &(tcr->tcQ[i]);
        for (j = 0; (j < md->nr); j++)
        {
            if (md->typeA[j] == tcq->at_i)
            {
                tcq->Q = md->chargeA[j];
                break;
            }
        }
        if (j == md->nr)
        {
            gmx_fatal(FARGS, "Coupling type %d not found", tcq->at_i);
        }
    }
    for (i = 0; (i < tcr->nIP); i++)
    {
        tip    = &(tcr->tIP[i]);
        type   = tip->type;
        ftype  = idef->functype[type];
        factor = dt*deviation[tip->eObs];

        switch (ftype)
        {
            case F_BONDS:
                if (tip->xi.harmonic.krA)
                {
                    idef->iparams[type].harmonic.krA *= (1+factor/tip->xi.harmonic.krA);
                }
                if (tip->xi.harmonic.rA)
                {
                    idef->iparams[type].harmonic.rA *= (1+factor/tip->xi.harmonic.rA);
                }
                break;
            default:
                break;
        }
        tip->iprint = idef->iparams[type];
    }
}

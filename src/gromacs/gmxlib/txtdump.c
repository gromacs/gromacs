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

/* This file is completely threadsafe - please keep it that way! */

#include "gromacs/legacyheaders/txtdump.h"

#include <stdio.h>
#include <stdlib.h>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

int pr_indent(FILE *fp, int n)
{
    int i;

    for (i = 0; i < n; i++)
    {
        (void) fprintf(fp, " ");
    }
    return n;
}

int available(FILE *fp, void *p, int indent, const char *title)
{
    if (!p)
    {
        if (indent > 0)
        {
            pr_indent(fp, indent);
        }
        (void) fprintf(fp, "%s: not available\n", title);
    }
    return (p != NULL);
}

int pr_title(FILE *fp, int indent, const char *title)
{
    (void) pr_indent(fp, indent);
    (void) fprintf(fp, "%s:\n", title);
    return (indent+INDENT);
}

int pr_title_n(FILE *fp, int indent, const char *title, int n)
{
    (void) pr_indent(fp, indent);
    (void) fprintf(fp, "%s (%d):\n", title, n);
    return (indent+INDENT);
}

int pr_title_nxn(FILE *fp, int indent, const char *title, int n1, int n2)
{
    (void) pr_indent(fp, indent);
    (void) fprintf(fp, "%s (%dx%d):\n", title, n1, n2);
    return (indent+INDENT);
}

void pr_ivec(FILE *fp, int indent, const char *title, int vec[], int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            (void) pr_indent(fp, indent);
            (void) fprintf(fp, "%s[%d]=%d\n", title, bShowNumbers ? i : -1, vec[i]);
        }
    }
}

void pr_ivec_block(FILE *fp, int indent, const char *title, int vec[], int n, gmx_bool bShowNumbers)
{
    int i, j;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        i      = 0;
        while (i < n)
        {
            j = i+1;
            while (j < n && vec[j] == vec[j-1]+1)
            {
                j++;
            }
            /* Print consecutive groups of 3 or more as blocks */
            if (j - i < 3)
            {
                while (i < j)
                {
                    (void) pr_indent(fp, indent);
                    (void) fprintf(fp, "%s[%d]=%d\n",
                                   title, bShowNumbers ? i : -1, vec[i]);
                    i++;
                }
            }
            else
            {
                (void) pr_indent(fp, indent);
                (void) fprintf(fp, "%s[%d,...,%d] = {%d,...,%d}\n",
                               title,
                               bShowNumbers ? i : -1,
                               bShowNumbers ? j-1 : -1,
                               vec[i], vec[j-1]);
                i = j;
            }
        }
    }
}

void pr_bvec(FILE *fp, int indent, const char *title, gmx_bool vec[], int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            (void) pr_indent(fp, indent);
            (void) fprintf(fp, "%s[%d]=%s\n", title, bShowNumbers ? i : -1,
                           EBOOL(vec[i]));
        }
    }
}

void pr_ivecs(FILE *fp, int indent, const char *title, ivec vec[], int n, gmx_bool bShowNumbers)
{
    int i, j;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, DIM);
        for (i = 0; i < n; i++)
        {
            (void) pr_indent(fp, indent);
            (void) fprintf(fp, "%s[%d]={", title, bShowNumbers ? i : -1);
            for (j = 0; j < DIM; j++)
            {
                if (j != 0)
                {
                    (void) fprintf(fp, ", ");
                }
                fprintf(fp, "%d", vec[i][j]);
            }
            (void) fprintf(fp, "}\n");
        }
    }
}

void pr_rvec(FILE *fp, int indent, const char *title, real vec[], int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]=%12.5e\n", title, bShowNumbers ? i : -1, vec[i]);
        }
    }
}

void pr_dvec(FILE *fp, int indent, const char *title, double vec[], int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]=%12.5e\n", title, bShowNumbers ? i : -1, vec[i]);
        }
    }
}


/*
   void pr_mat(FILE *fp,int indent,char *title,matrix m)
   {
   int i,j;

   if (available(fp,m,indent,title)) {
    indent=pr_title_n(fp,indent,title,n);
    for(i=0; i<n; i++) {
      pr_indent(fp,indent);
      fprintf(fp,"%s[%d]=%12.5e %12.5e %12.5e\n",
          title,bShowNumbers?i:-1,m[i][XX],m[i][YY],m[i][ZZ]);
    }
   }
   }
 */

void pr_rvecs_len(FILE *fp, int indent, const char *title, rvec vec[], int n)
{
    int i, j;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, DIM);
        for (i = 0; i < n; i++)
        {
            (void) pr_indent(fp, indent);
            (void) fprintf(fp, "%s[%5d]={", title, i);
            for (j = 0; j < DIM; j++)
            {
                if (j != 0)
                {
                    (void) fprintf(fp, ", ");
                }
                (void) fprintf(fp, "%12.5e", vec[i][j]);
            }
            (void) fprintf(fp, "} len=%12.5e\n", norm(vec[i]));
        }
    }
}

void pr_rvecs(FILE *fp, int indent, const char *title, rvec vec[], int n)
{
    const char *fshort = "%12.5e";
    const char *flong  = "%15.8e";
    const char *format;
    int         i, j;

    if (getenv("GMX_PRINT_LONGFORMAT") != NULL)
    {
        format = flong;
    }
    else
    {
        format = fshort;
    }

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, DIM);
        for (i = 0; i < n; i++)
        {
            (void) pr_indent(fp, indent);
            (void) fprintf(fp, "%s[%5d]={", title, i);
            for (j = 0; j < DIM; j++)
            {
                if (j != 0)
                {
                    (void) fprintf(fp, ", ");
                }
                (void) fprintf(fp, format, vec[i][j]);
            }
            (void) fprintf(fp, "}\n");
        }
    }
}


void pr_rvecs_of_dim(FILE *fp, int indent, const char *title, rvec vec[], int n, int dim)
{
    const char *fshort = "%12.5e";
    const char *flong  = "%15.8e";
    const char *format;
    int         i, j;

    if (getenv("GMX_PRINT_LONGFORMAT") != NULL)
    {
        format = flong;
    }
    else
    {
        format = fshort;
    }

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, dim);
        for (i = 0; i < n; i++)
        {
            (void) pr_indent(fp, indent);
            (void) fprintf(fp, "%s[%5d]={", title, i);
            for (j = 0; j < dim; j++)
            {
                if (j != 0)
                {
                    (void) fprintf(fp, ", ");
                }
                (void) fprintf(fp, format, vec[i][j]);
            }
            (void) fprintf(fp, "}\n");
        }
    }
}

void pr_reals(FILE *fp, int indent, const char *title, real *vec, int n)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "%s:\t", title);
        for (i = 0; i < n; i++)
        {
            fprintf(fp, "  %10g", vec[i]);
        }
        (void) fprintf(fp, "\n");
    }
}

void pr_doubles(FILE *fp, int indent, const char *title, double *vec, int n)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "%s:\t", title);
        for (i = 0; i < n; i++)
        {
            fprintf(fp, "  %10g", vec[i]);
        }
        (void) fprintf(fp, "\n");
    }
}

void pr_reals_of_dim(FILE *fp, int indent, const char *title, real *vec, int n, int dim)
{
    int         i, j;
    const char *fshort = "%12.5e";
    const char *flong  = "%15.8e";
    const char *format;

    if (getenv("GMX_PRINT_LONGFORMAT") != NULL)
    {
        format = flong;
    }
    else
    {
        format = fshort;
    }

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, dim);
        for (i = 0; i < n; i++)
        {
            (void) pr_indent(fp, indent);
            (void) fprintf(fp, "%s[%5d]={", title, i);
            for (j = 0; j < dim; j++)
            {
                if (j != 0)
                {
                    (void) fprintf(fp, ", ");
                }
                (void) fprintf(fp, format, vec[i * dim  + j]);
            }
            (void) fprintf(fp, "}\n");
        }
    }
}

static void pr_int(FILE *fp, int indent, const char *title, int i)
{
    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %d\n", title, i);
}

static void pr_int64(FILE *fp, int indent, const char *title, gmx_int64_t i)
{
    char buf[STEPSTRSIZE];

    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %s\n", title, gmx_step_str(i, buf));
}

static void pr_real(FILE *fp, int indent, const char *title, real r)
{
    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %g\n", title, r);
}

static void pr_double(FILE *fp, int indent, const char *title, double d)
{
    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %g\n", title, d);
}

static void pr_str(FILE *fp, int indent, const char *title, const char *s)
{
    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %s\n", title, s);
}

void pr_qm_opts(FILE *fp, int indent, const char *title, t_grpopts *opts)
{
    int i, m, j;

    fprintf(fp, "%s:\n", title);

    pr_int(fp, indent, "ngQM", opts->ngQM);
    if (opts->ngQM > 0)
    {
        pr_ivec(fp, indent, "QMmethod", opts->QMmethod, opts->ngQM, FALSE);
        pr_ivec(fp, indent, "QMbasis", opts->QMbasis, opts->ngQM, FALSE);
        pr_ivec(fp, indent, "QMcharge", opts->QMcharge, opts->ngQM, FALSE);
        pr_ivec(fp, indent, "QMmult", opts->QMmult, opts->ngQM, FALSE);
        pr_bvec(fp, indent, "SH", opts->bSH, opts->ngQM, FALSE);
        pr_ivec(fp, indent, "CASorbitals", opts->CASorbitals, opts->ngQM, FALSE);
        pr_ivec(fp, indent, "CASelectrons", opts->CASelectrons, opts->ngQM, FALSE);
        pr_rvec(fp, indent, "SAon", opts->SAon, opts->ngQM, FALSE);
        pr_rvec(fp, indent, "SAoff", opts->SAoff, opts->ngQM, FALSE);
        pr_ivec(fp, indent, "SAsteps", opts->SAsteps, opts->ngQM, FALSE);
        pr_bvec(fp, indent, "bOPT", opts->bOPT, opts->ngQM, FALSE);
        pr_bvec(fp, indent, "bTS", opts->bTS, opts->ngQM, FALSE);
    }
}

static void pr_grp_opts(FILE *out, int indent, const char *title, t_grpopts *opts,
                        gmx_bool bMDPformat)
{
    int i, m, j;

    if (!bMDPformat)
    {
        fprintf(out, "%s:\n", title);
    }

    pr_indent(out, indent);
    fprintf(out, "nrdf%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10g", opts->nrdf[i]);
    }
    fprintf(out, "\n");

    pr_indent(out, indent);
    fprintf(out, "ref-t%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10g", opts->ref_t[i]);
    }
    fprintf(out, "\n");

    pr_indent(out, indent);
    fprintf(out, "tau-t%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10g", opts->tau_t[i]);
    }
    fprintf(out, "\n");

    /* Pretty-print the simulated annealing info */
    fprintf(out, "annealing%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10s", EANNEAL(opts->annealing[i]));
    }
    fprintf(out, "\n");

    fprintf(out, "annealing-npoints%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10d", opts->anneal_npoints[i]);
    }
    fprintf(out, "\n");

    for (i = 0; (i < opts->ngtc); i++)
    {
        if (opts->anneal_npoints[i] > 0)
        {
            fprintf(out, "annealing-time [%d]:\t", i);
            for (j = 0; (j < opts->anneal_npoints[i]); j++)
            {
                fprintf(out, "  %10.1f", opts->anneal_time[i][j]);
            }
            fprintf(out, "\n");
            fprintf(out, "annealing-temp [%d]:\t", i);
            for (j = 0; (j < opts->anneal_npoints[i]); j++)
            {
                fprintf(out, "  %10.1f", opts->anneal_temp[i][j]);
            }
            fprintf(out, "\n");
        }
    }

    pr_indent(out, indent);
    fprintf(out, "acc:\t");
    for (i = 0; (i < opts->ngacc); i++)
    {
        for (m = 0; (m < DIM); m++)
        {
            fprintf(out, "  %10g", opts->acc[i][m]);
        }
    }
    fprintf(out, "\n");

    pr_indent(out, indent);
    fprintf(out, "nfreeze:");
    for (i = 0; (i < opts->ngfrz); i++)
    {
        for (m = 0; (m < DIM); m++)
        {
            fprintf(out, "  %10s", opts->nFreeze[i][m] ? "Y" : "N");
        }
    }
    fprintf(out, "\n");


    for (i = 0; (i < opts->ngener); i++)
    {
        pr_indent(out, indent);
        fprintf(out, "energygrp-flags[%3d]:", i);
        for (m = 0; (m < opts->ngener); m++)
        {
            fprintf(out, " %d", opts->egp_flags[opts->ngener*i+m]);
        }
        fprintf(out, "\n");
    }

    fflush(out);
}

static void pr_matrix(FILE *fp, int indent, const char *title, rvec *m,
                      gmx_bool bMDPformat)
{
    if (bMDPformat)
    {
        fprintf(fp, "%-10s    = %g %g %g %g %g %g\n", title,
                m[XX][XX], m[YY][YY], m[ZZ][ZZ], m[XX][YY], m[XX][ZZ], m[YY][ZZ]);
    }
    else
    {
        pr_rvecs(fp, indent, title, m, DIM);
    }
}

static void pr_cosine(FILE *fp, int indent, const char *title, t_cosines *cos,
                      gmx_bool bMDPformat)
{
    int j;

    if (bMDPformat)
    {
        fprintf(fp, "%s = %d\n", title, cos->n);
    }
    else
    {
        indent = pr_title(fp, indent, title);
        (void) pr_indent(fp, indent);
        fprintf(fp, "n = %d\n", cos->n);
        if (cos->n > 0)
        {
            (void) pr_indent(fp, indent+2);
            fprintf(fp, "a =");
            for (j = 0; (j < cos->n); j++)
            {
                fprintf(fp, " %e", cos->a[j]);
            }
            fprintf(fp, "\n");
            (void) pr_indent(fp, indent+2);
            fprintf(fp, "phi =");
            for (j = 0; (j < cos->n); j++)
            {
                fprintf(fp, " %e", cos->phi[j]);
            }
            fprintf(fp, "\n");
        }
    }
}

#define PS(t, s) pr_str(fp, indent, t, s)
#define PI(t, s) pr_int(fp, indent, t, s)
#define PSTEP(t, s) pr_int64(fp, indent, t, s)
#define PR(t, s) pr_real(fp, indent, t, s)
#define PD(t, s) pr_double(fp, indent, t, s)

static void pr_pull_group(FILE *fp, int indent, int g, t_pull_group *pgrp)
{
    pr_indent(fp, indent);
    fprintf(fp, "pull-group %d:\n", g);
    indent += 2;
    pr_ivec_block(fp, indent, "atom", pgrp->ind, pgrp->nat, TRUE);
    pr_rvec(fp, indent, "weight", pgrp->weight, pgrp->nweight, TRUE);
    PI("pbcatom", pgrp->pbcatom);
}

static void pr_pull_coord(FILE *fp, int indent, int c, t_pull_coord *pcrd)
{
    pr_indent(fp, indent);
    fprintf(fp, "pull-coord %d:\n", c);
    PI("group[0]", pcrd->group[0]);
    PI("group[1]", pcrd->group[1]);
    if (pcrd->eGeom == epullgDIRRELATIVE)
    {
        PI("group[2]", pcrd->group[2]);
        PI("group[3]", pcrd->group[3]);
    }
    PS("type", EPULLTYPE(pcrd->eType));
    PS("geometry", EPULLGEOM(pcrd->eGeom));
    pr_ivec(fp, indent, "dim", pcrd->dim, DIM, TRUE);
    pr_rvec(fp, indent, "origin", pcrd->origin, DIM, TRUE);
    pr_rvec(fp, indent, "vec", pcrd->vec, DIM, TRUE);
    PS("start", EBOOL(pcrd->bStart));
    PR("init", pcrd->init);
    PR("rate", pcrd->rate);
    PR("k", pcrd->k);
    PR("kB", pcrd->kB);
}

static void pr_simtempvals(FILE *fp, int indent, t_simtemp *simtemp, int n_lambda)
{
    PS("simulated-tempering-scaling", ESIMTEMP(simtemp->eSimTempScale));
    PR("sim-temp-low", simtemp->simtemp_low);
    PR("sim-temp-high", simtemp->simtemp_high);
    pr_rvec(fp, indent, "simulated tempering temperatures", simtemp->temperatures, n_lambda, TRUE);
}

static void pr_expandedvals(FILE *fp, int indent, t_expanded *expand, int n_lambda)
{

    PI("nstexpanded", expand->nstexpanded);
    PS("lmc-stats", elamstats_names[expand->elamstats]);
    PS("lmc-move", elmcmove_names[expand->elmcmove]);
    PS("lmc-weights-equil", elmceq_names[expand->elmceq]);
    if (expand->elmceq == elmceqNUMATLAM)
    {
        PI("weight-equil-number-all-lambda", expand->equil_n_at_lam);
    }
    if (expand->elmceq == elmceqSAMPLES)
    {
        PI("weight-equil-number-samples", expand->equil_samples);
    }
    if (expand->elmceq == elmceqSTEPS)
    {
        PI("weight-equil-number-steps", expand->equil_steps);
    }
    if (expand->elmceq == elmceqWLDELTA)
    {
        PR("weight-equil-wl-delta", expand->equil_wl_delta);
    }
    if (expand->elmceq == elmceqRATIO)
    {
        PR("weight-equil-count-ratio", expand->equil_ratio);
    }
    PI("lmc-seed", expand->lmc_seed);
    PR("mc-temperature", expand->mc_temp);
    PI("lmc-repeats", expand->lmc_repeats);
    PI("lmc-gibbsdelta", expand->gibbsdeltalam);
    PI("lmc-forced-nstart", expand->lmc_forced_nstart);
    PS("symmetrized-transition-matrix", EBOOL(expand->bSymmetrizedTMatrix));
    PI("nst-transition-matrix", expand->nstTij);
    PI("mininum-var-min", expand->minvarmin); /*default is reasonable */
    PI("weight-c-range", expand->c_range);    /* default is just C=0 */
    PR("wl-scale", expand->wl_scale);
    PR("wl-ratio", expand->wl_ratio);
    PR("init-wl-delta", expand->init_wl_delta);
    PS("wl-oneovert", EBOOL(expand->bWLoneovert));

    pr_indent(fp, indent);
    pr_rvec(fp, indent, "init-lambda-weights", expand->init_lambda_weights, n_lambda, TRUE);
    PS("init-weights", EBOOL(expand->bInit_weights));
}

static void pr_fepvals(FILE *fp, int indent, t_lambda *fep, gmx_bool bMDPformat)
{
    int i, j;

    PR("init-lambda", fep->init_lambda);
    PI("init-lambda-state", fep->init_fep_state);
    PR("delta-lambda", fep->delta_lambda);
    PI("nstdhdl", fep->nstdhdl);

    if (!bMDPformat)
    {
        PI("n-lambdas", fep->n_lambda);
    }
    if (fep->n_lambda > 0)
    {
        pr_indent(fp, indent);
        fprintf(fp, "separate-dvdl%s\n", bMDPformat ? " = " : ":");
        for (i = 0; i < efptNR; i++)
        {
            fprintf(fp, "%18s = ", efpt_names[i]);
            if (fep->separate_dvdl[i])
            {
                fprintf(fp, "  TRUE");
            }
            else
            {
                fprintf(fp, "  FALSE");
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "all-lambdas%s\n", bMDPformat ? " = " : ":");
        for (i = 0; i < efptNR; i++)
        {
            fprintf(fp, "%18s = ", efpt_names[i]);
            for (j = 0; j < fep->n_lambda; j++)
            {
                fprintf(fp, "  %10g", fep->all_lambda[i][j]);
            }
            fprintf(fp, "\n");
        }
    }
    PI("calc-lambda-neighbors", fep->lambda_neighbors);
    PS("dhdl-print-energy", edHdLPrintEnergy_names[fep->edHdLPrintEnergy]);
    PR("sc-alpha", fep->sc_alpha);
    PI("sc-power", fep->sc_power);
    PR("sc-r-power", fep->sc_r_power);
    PR("sc-sigma", fep->sc_sigma);
    PR("sc-sigma-min", fep->sc_sigma_min);
    PS("sc-coul", EBOOL(fep->bScCoul));
    PI("dh-hist-size", fep->dh_hist_size);
    PD("dh-hist-spacing", fep->dh_hist_spacing);
    PS("separate-dhdl-file", SEPDHDLFILETYPE(fep->separate_dhdl_file));
    PS("dhdl-derivatives", DHDLDERIVATIVESTYPE(fep->dhdl_derivatives));
};

static void pr_pull(FILE *fp, int indent, pull_params_t *pull)
{
    int g;

    PR("pull-cylinder-r", pull->cylinder_r);
    PR("pull-constr-tol", pull->constr_tol);
    PS("pull-print-COM1", EBOOL(pull->bPrintCOM1));
    PS("pull-print-COM2", EBOOL(pull->bPrintCOM2));
    PS("pull-print-ref-value", EBOOL(pull->bPrintRefValue));
    PS("pull-print-components", EBOOL(pull->bPrintComp));
    PI("pull-nstxout", pull->nstxout);
    PI("pull-nstfout", pull->nstfout);
    PI("pull-ngroups", pull->ngroup);
    for (g = 0; g < pull->ngroup; g++)
    {
        pr_pull_group(fp, indent, g, &pull->group[g]);
    }
    PI("pull-ncoords", pull->ncoord);
    for (g = 0; g < pull->ncoord; g++)
    {
        pr_pull_coord(fp, indent, g, &pull->coord[g]);
    }
}

static void pr_rotgrp(FILE *fp, int indent, int g, t_rotgrp *rotg)
{
    pr_indent(fp, indent);
    fprintf(fp, "rot-group %d:\n", g);
    indent += 2;
    PS("rot-type", EROTGEOM(rotg->eType));
    PS("rot-massw", EBOOL(rotg->bMassW));
    pr_ivec_block(fp, indent, "atom", rotg->ind, rotg->nat, TRUE);
    pr_rvecs(fp, indent, "x-ref", rotg->x_ref, rotg->nat);
    pr_rvec(fp, indent, "rot-vec", rotg->vec, DIM, TRUE);
    pr_rvec(fp, indent, "rot-pivot", rotg->pivot, DIM, TRUE);
    PR("rot-rate", rotg->rate);
    PR("rot-k", rotg->k);
    PR("rot-slab-dist", rotg->slab_dist);
    PR("rot-min-gauss", rotg->min_gaussian);
    PR("rot-eps", rotg->eps);
    PS("rot-fit-method", EROTFIT(rotg->eFittype));
    PI("rot-potfit-nstep", rotg->PotAngle_nstep);
    PR("rot-potfit-step", rotg->PotAngle_step);
}

static void pr_rot(FILE *fp, int indent, t_rot *rot)
{
    int g;

    PI("rot-nstrout", rot->nstrout);
    PI("rot-nstsout", rot->nstsout);
    PI("rot-ngroups", rot->ngrp);
    for (g = 0; g < rot->ngrp; g++)
    {
        pr_rotgrp(fp, indent, g, &rot->grp[g]);
    }
}


static void pr_swap(FILE *fp, int indent, t_swapcoords *swap)
{
    int  i, j;
    char str[STRLEN];


    PI("swap-frequency", swap->nstswap);
    for (j = 0; j < 2; j++)
    {
        sprintf(str, "massw_split%d", j);
        PS(str, EBOOL(swap->massw_split[j]));
        sprintf(str, "split atoms group %d", j);
        pr_ivec_block(fp, indent, str, swap->ind_split[j], swap->nat_split[j], TRUE);
    }
    pr_ivec_block(fp, indent, "swap atoms", swap->ind, swap->nat, TRUE);
    pr_ivec_block(fp, indent, "solvent atoms", swap->ind_sol, swap->nat_sol, TRUE);
    PR("cyl0-r", swap->cyl0r);
    PR("cyl0-up", swap->cyl0u);
    PR("cyl0-down", swap->cyl0l);
    PR("cyl1-r", swap->cyl1r);
    PR("cyl1-up", swap->cyl1u);
    PR("cyl1-down", swap->cyl1l);
    PI("coupl-steps", swap->nAverage);
    for (j = 0; j < 2; j++)
    {
        sprintf(str, "anions%c", j+'A');
        PI(str, swap->nanions[j]);
        sprintf(str, "cations%c", j+'A');
        PI(str, swap->ncations[j]);
    }
    PR("threshold", swap->threshold);
}


static void pr_imd(FILE *fp, int indent, t_IMD *imd)
{
    PI("IMD-atoms", imd->nat);
    pr_ivec_block(fp, indent, "atom", imd->ind, imd->nat, TRUE);
}


void pr_inputrec(FILE *fp, int indent, const char *title, t_inputrec *ir,
                 gmx_bool bMDPformat)
{
    const char *infbuf = "inf";
    int         i;

    if (available(fp, ir, indent, title))
    {
        if (!bMDPformat)
        {
            indent = pr_title(fp, indent, title);
        }
        /* Try to make this list appear in the same order as the
         * options are written in the default mdout.mdp, and with
         * the same user-exposed names to facilitate debugging.
         */
        PS("integrator", EI(ir->eI));
        PR("tinit", ir->init_t);
        PR("dt", ir->delta_t);
        PSTEP("nsteps", ir->nsteps);
        PSTEP("init-step", ir->init_step);
        PI("simulation-part", ir->simulation_part);
        PS("comm-mode", ECOM(ir->comm_mode));
        PI("nstcomm", ir->nstcomm);

        /* Langevin dynamics */
        PR("bd-fric", ir->bd_fric);
        PSTEP("ld-seed", ir->ld_seed);

        /* Energy minimization */
        PR("emtol", ir->em_tol);
        PR("emstep", ir->em_stepsize);
        PI("niter", ir->niter);
        PR("fcstep", ir->fc_stepsize);
        PI("nstcgsteep", ir->nstcgsteep);
        PI("nbfgscorr", ir->nbfgscorr);

        /* Test particle insertion */
        PR("rtpi", ir->rtpi);

        /* Output control */
        PI("nstxout", ir->nstxout);
        PI("nstvout", ir->nstvout);
        PI("nstfout", ir->nstfout);
        PI("nstlog", ir->nstlog);
        PI("nstcalcenergy", ir->nstcalcenergy);
        PI("nstenergy", ir->nstenergy);
        PI("nstxout-compressed", ir->nstxout_compressed);
        PR("compressed-x-precision", ir->x_compression_precision);

        /* Neighborsearching parameters */
        PS("cutoff-scheme", ECUTSCHEME(ir->cutoff_scheme));
        PI("nstlist", ir->nstlist);
        PS("ns-type", ENS(ir->ns_type));
        PS("pbc", EPBC(ir->ePBC));
        PS("periodic-molecules", EBOOL(ir->bPeriodicMols));
        PR("verlet-buffer-tolerance", ir->verletbuf_tol);
        PR("rlist", ir->rlist);
        PR("rlistlong", ir->rlistlong);
        PR("nstcalclr", ir->nstcalclr);

        /* Options for electrostatics and VdW */
        PS("coulombtype", EELTYPE(ir->coulombtype));
        PS("coulomb-modifier", INTMODIFIER(ir->coulomb_modifier));
        PR("rcoulomb-switch", ir->rcoulomb_switch);
        PR("rcoulomb", ir->rcoulomb);
        if (ir->epsilon_r != 0)
        {
            PR("epsilon-r", ir->epsilon_r);
        }
        else
        {
            PS("epsilon-r", infbuf);
        }
        if (ir->epsilon_rf != 0)
        {
            PR("epsilon-rf", ir->epsilon_rf);
        }
        else
        {
            PS("epsilon-rf", infbuf);
        }
        PS("vdw-type", EVDWTYPE(ir->vdwtype));
        PS("vdw-modifier", INTMODIFIER(ir->vdw_modifier));
        PR("rvdw-switch", ir->rvdw_switch);
        PR("rvdw", ir->rvdw);
        PS("DispCorr", EDISPCORR(ir->eDispCorr));
        PR("table-extension", ir->tabext);

        PR("fourierspacing", ir->fourier_spacing);
        PI("fourier-nx", ir->nkx);
        PI("fourier-ny", ir->nky);
        PI("fourier-nz", ir->nkz);
        PI("pme-order", ir->pme_order);
        PR("ewald-rtol", ir->ewald_rtol);
        PR("ewald-rtol-lj", ir->ewald_rtol_lj);
        PS("lj-pme-comb-rule", ELJPMECOMBNAMES(ir->ljpme_combination_rule));
        PR("ewald-geometry", ir->ewald_geometry);
        PR("epsilon-surface", ir->epsilon_surface);

        /* Implicit solvent */
        PS("implicit-solvent", EIMPLICITSOL(ir->implicit_solvent));

        /* Generalized born electrostatics */
        PS("gb-algorithm", EGBALGORITHM(ir->gb_algorithm));
        PI("nstgbradii", ir->nstgbradii);
        PR("rgbradii", ir->rgbradii);
        PR("gb-epsilon-solvent", ir->gb_epsilon_solvent);
        PR("gb-saltconc", ir->gb_saltconc);
        PR("gb-obc-alpha", ir->gb_obc_alpha);
        PR("gb-obc-beta", ir->gb_obc_beta);
        PR("gb-obc-gamma", ir->gb_obc_gamma);
        PR("gb-dielectric-offset", ir->gb_dielectric_offset);
        PS("sa-algorithm", ESAALGORITHM(ir->sa_algorithm));
        PR("sa-surface-tension", ir->sa_surface_tension);

        /* Options for weak coupling algorithms */
        PS("tcoupl", ETCOUPLTYPE(ir->etc));
        PI("nsttcouple", ir->nsttcouple);
        PI("nh-chain-length", ir->opts.nhchainlength);
        PS("print-nose-hoover-chain-variables", EBOOL(ir->bPrintNHChains));

        PS("pcoupl", EPCOUPLTYPE(ir->epc));
        PS("pcoupltype", EPCOUPLTYPETYPE(ir->epct));
        PI("nstpcouple", ir->nstpcouple);
        PR("tau-p", ir->tau_p);
        pr_matrix(fp, indent, "compressibility", ir->compress, bMDPformat);
        pr_matrix(fp, indent, "ref-p", ir->ref_p, bMDPformat);
        PS("refcoord-scaling", EREFSCALINGTYPE(ir->refcoord_scaling));

        if (bMDPformat)
        {
            fprintf(fp, "posres-com  = %g %g %g\n", ir->posres_com[XX],
                    ir->posres_com[YY], ir->posres_com[ZZ]);
            fprintf(fp, "posres-comB = %g %g %g\n", ir->posres_comB[XX],
                    ir->posres_comB[YY], ir->posres_comB[ZZ]);
        }
        else
        {
            pr_rvec(fp, indent, "posres-com", ir->posres_com, DIM, TRUE);
            pr_rvec(fp, indent, "posres-comB", ir->posres_comB, DIM, TRUE);
        }

        /* QMMM */
        PS("QMMM", EBOOL(ir->bQMMM));
        PI("QMconstraints", ir->QMconstraints);
        PI("QMMMscheme", ir->QMMMscheme);
        PR("MMChargeScaleFactor", ir->scalefactor);
        pr_qm_opts(fp, indent, "qm-opts", &(ir->opts));

        /* CONSTRAINT OPTIONS */
        PS("constraint-algorithm", ECONSTRTYPE(ir->eConstrAlg));
        PS("continuation", EBOOL(ir->bContinuation));

        PS("Shake-SOR", EBOOL(ir->bShakeSOR));
        PR("shake-tol", ir->shake_tol);
        PI("lincs-order", ir->nProjOrder);
        PI("lincs-iter", ir->nLincsIter);
        PR("lincs-warnangle", ir->LincsWarnAngle);

        /* Walls */
        PI("nwall", ir->nwall);
        PS("wall-type", EWALLTYPE(ir->wall_type));
        PR("wall-r-linpot", ir->wall_r_linpot);
        /* wall-atomtype */
        PI("wall-atomtype[0]", ir->wall_atomtype[0]);
        PI("wall-atomtype[1]", ir->wall_atomtype[1]);
        /* wall-density */
        PR("wall-density[0]", ir->wall_density[0]);
        PR("wall-density[1]", ir->wall_density[1]);
        PR("wall-ewald-zfac", ir->wall_ewald_zfac);

        /* COM PULLING */
        PS("pull", EBOOL(ir->bPull));
        if (ir->bPull)
        {
            pr_pull(fp, indent, ir->pull);
        }

        /* ENFORCED ROTATION */
        PS("rotation", EBOOL(ir->bRot));
        if (ir->bRot)
        {
            pr_rot(fp, indent, ir->rot);
        }

        /* INTERACTIVE MD */
        PS("interactiveMD", EBOOL(ir->bIMD));
        if (ir->bIMD)
        {
            pr_imd(fp, indent, ir->imd);
        }

        /* NMR refinement stuff */
        PS("disre", EDISRETYPE(ir->eDisre));
        PS("disre-weighting", EDISREWEIGHTING(ir->eDisreWeighting));
        PS("disre-mixed", EBOOL(ir->bDisreMixed));
        PR("dr-fc", ir->dr_fc);
        PR("dr-tau", ir->dr_tau);
        PR("nstdisreout", ir->nstdisreout);

        PR("orire-fc", ir->orires_fc);
        PR("orire-tau", ir->orires_tau);
        PR("nstorireout", ir->nstorireout);

        /* FREE ENERGY VARIABLES */
        PS("free-energy", EFEPTYPE(ir->efep));
        if (ir->efep != efepNO || ir->bSimTemp)
        {
            pr_fepvals(fp, indent, ir->fepvals, bMDPformat);
        }
        if (ir->bExpanded)
        {
            pr_expandedvals(fp, indent, ir->expandedvals, ir->fepvals->n_lambda);
        }

        /* NON-equilibrium MD stuff */
        PR("cos-acceleration", ir->cos_accel);
        pr_matrix(fp, indent, "deform", ir->deform, bMDPformat);

        /* SIMULATED TEMPERING */
        PS("simulated-tempering", EBOOL(ir->bSimTemp));
        if (ir->bSimTemp)
        {
            pr_simtempvals(fp, indent, ir->simtempvals, ir->fepvals->n_lambda);
        }

        /* ELECTRIC FIELDS */
        pr_cosine(fp, indent, "E-x", &(ir->ex[XX]), bMDPformat);
        pr_cosine(fp, indent, "E-xt", &(ir->et[XX]), bMDPformat);
        pr_cosine(fp, indent, "E-y", &(ir->ex[YY]), bMDPformat);
        pr_cosine(fp, indent, "E-yt", &(ir->et[YY]), bMDPformat);
        pr_cosine(fp, indent, "E-z", &(ir->ex[ZZ]), bMDPformat);
        pr_cosine(fp, indent, "E-zt", &(ir->et[ZZ]), bMDPformat);

        /* ION/WATER SWAPPING FOR COMPUTATIONAL ELECTROPHYSIOLOGY */
        PS("swapcoords", ESWAPTYPE(ir->eSwapCoords));
        if (ir->eSwapCoords != eswapNO)
        {
            pr_swap(fp, indent, ir->swap);
        }

        /* AdResS PARAMETERS */
        PS("adress", EBOOL(ir->bAdress));
        if (ir->bAdress)
        {
            PS("adress-type", EADRESSTYPE(ir->adress->type));
            PR("adress-const-wf", ir->adress->const_wf);
            PR("adress-ex-width", ir->adress->ex_width);
            PR("adress-hy-width", ir->adress->hy_width);
            PR("adress-ex-forcecap", ir->adress->ex_forcecap);
            PS("adress-interface-correction", EADRESSICTYPE(ir->adress->icor));
            PS("adress-site", EADRESSSITETYPE(ir->adress->site));
            pr_rvec(fp, indent, "adress-reference-coords", ir->adress->refs, DIM, TRUE);
            PS("adress-do-hybridpairs", EBOOL(ir->adress->do_hybridpairs));
        }

        /* USER-DEFINED THINGIES */
        PI("userint1", ir->userint1);
        PI("userint2", ir->userint2);
        PI("userint3", ir->userint3);
        PI("userint4", ir->userint4);
        PR("userreal1", ir->userreal1);
        PR("userreal2", ir->userreal2);
        PR("userreal3", ir->userreal3);
        PR("userreal4", ir->userreal4);

        pr_grp_opts(fp, indent, "grpopts", &(ir->opts), bMDPformat);
    }
}
#undef PS
#undef PR
#undef PI

static void pr_harm(FILE *fp, t_iparams *iparams, const char *r, const char *kr)
{
    fprintf(fp, "%sA=%12.5e, %sA=%12.5e, %sB=%12.5e, %sB=%12.5e\n",
            r, iparams->harmonic.rA, kr, iparams->harmonic.krA,
            r, iparams->harmonic.rB, kr, iparams->harmonic.krB);
}

void pr_iparams(FILE *fp, t_functype ftype, t_iparams *iparams)
{
    int  i;
    real VA[4], VB[4], *rbcA, *rbcB;

    switch (ftype)
    {
        case F_ANGLES:
        case F_G96ANGLES:
            pr_harm(fp, iparams, "th", "ct");
            break;
        case F_CROSS_BOND_BONDS:
            fprintf(fp, "r1e=%15.8e, r2e=%15.8e, krr=%15.8e\n",
                    iparams->cross_bb.r1e, iparams->cross_bb.r2e,
                    iparams->cross_bb.krr);
            break;
        case F_CROSS_BOND_ANGLES:
            fprintf(fp, "r1e=%15.8e, r1e=%15.8e, r3e=%15.8e, krt=%15.8e\n",
                    iparams->cross_ba.r1e, iparams->cross_ba.r2e,
                    iparams->cross_ba.r3e, iparams->cross_ba.krt);
            break;
        case F_LINEAR_ANGLES:
            fprintf(fp, "klinA=%15.8e, aA=%15.8e, klinB=%15.8e, aB=%15.8e\n",
                    iparams->linangle.klinA, iparams->linangle.aA,
                    iparams->linangle.klinB, iparams->linangle.aB);
            break;
        case F_UREY_BRADLEY:
            fprintf(fp, "thetaA=%15.8e, kthetaA=%15.8e, r13A=%15.8e, kUBA=%15.8e, thetaB=%15.8e, kthetaB=%15.8e, r13B=%15.8e, kUBB=%15.8e\n", iparams->u_b.thetaA, iparams->u_b.kthetaA, iparams->u_b.r13A, iparams->u_b.kUBA, iparams->u_b.thetaB, iparams->u_b.kthetaB, iparams->u_b.r13B, iparams->u_b.kUBB);
            break;
        case F_QUARTIC_ANGLES:
            fprintf(fp, "theta=%15.8e", iparams->qangle.theta);
            for (i = 0; i < 5; i++)
            {
                fprintf(fp, ", c%c=%15.8e", '0'+i, iparams->qangle.c[i]);
            }
            fprintf(fp, "\n");
            break;
        case F_BHAM:
            fprintf(fp, "a=%15.8e, b=%15.8e, c=%15.8e\n",
                    iparams->bham.a, iparams->bham.b, iparams->bham.c);
            break;
        case F_BONDS:
        case F_G96BONDS:
        case F_HARMONIC:
            pr_harm(fp, iparams, "b0", "cb");
            break;
        case F_IDIHS:
            pr_harm(fp, iparams, "xi", "cx");
            break;
        case F_MORSE:
            fprintf(fp, "b0A=%15.8e, cbA=%15.8e, betaA=%15.8e, b0B=%15.8e, cbB=%15.8e, betaB=%15.8e\n",
                    iparams->morse.b0A, iparams->morse.cbA, iparams->morse.betaA,
                    iparams->morse.b0B, iparams->morse.cbB, iparams->morse.betaB);
            break;
        case F_CUBICBONDS:
            fprintf(fp, "b0=%15.8e, kb=%15.8e, kcub=%15.8e\n",
                    iparams->cubic.b0, iparams->cubic.kb, iparams->cubic.kcub);
            break;
        case F_CONNBONDS:
            fprintf(fp, "\n");
            break;
        case F_FENEBONDS:
            fprintf(fp, "bm=%15.8e, kb=%15.8e\n", iparams->fene.bm, iparams->fene.kb);
            break;
        case F_RESTRBONDS:
            fprintf(fp, "lowA=%15.8e, up1A=%15.8e, up2A=%15.8e, kA=%15.8e, lowB=%15.8e, up1B=%15.8e, up2B=%15.8e, kB=%15.8e,\n",
                    iparams->restraint.lowA, iparams->restraint.up1A,
                    iparams->restraint.up2A, iparams->restraint.kA,
                    iparams->restraint.lowB, iparams->restraint.up1B,
                    iparams->restraint.up2B, iparams->restraint.kB);
            break;
        case F_TABBONDS:
        case F_TABBONDSNC:
        case F_TABANGLES:
        case F_TABDIHS:
            fprintf(fp, "tab=%d, kA=%15.8e, kB=%15.8e\n",
                    iparams->tab.table, iparams->tab.kA, iparams->tab.kB);
            break;
        case F_POLARIZATION:
            fprintf(fp, "alpha=%15.8e\n", iparams->polarize.alpha);
            break;
        case F_ANHARM_POL:
            fprintf(fp, "alpha=%15.8e drcut=%15.8e khyp=%15.8e\n",
                    iparams->anharm_polarize.alpha,
                    iparams->anharm_polarize.drcut,
                    iparams->anharm_polarize.khyp);
            break;
        case F_THOLE_POL:
            fprintf(fp, "a=%15.8e, alpha1=%15.8e, alpha2=%15.8e, rfac=%15.8e\n",
                    iparams->thole.a, iparams->thole.alpha1, iparams->thole.alpha2,
                    iparams->thole.rfac);
            break;
        case F_WATER_POL:
            fprintf(fp, "al_x=%15.8e, al_y=%15.8e, al_z=%15.8e, rOH=%9.6f, rHH=%9.6f, rOD=%9.6f\n",
                    iparams->wpol.al_x, iparams->wpol.al_y, iparams->wpol.al_z,
                    iparams->wpol.rOH, iparams->wpol.rHH, iparams->wpol.rOD);
            break;
        case F_LJ:
            fprintf(fp, "c6=%15.8e, c12=%15.8e\n", iparams->lj.c6, iparams->lj.c12);
            break;
        case F_LJ14:
            fprintf(fp, "c6A=%15.8e, c12A=%15.8e, c6B=%15.8e, c12B=%15.8e\n",
                    iparams->lj14.c6A, iparams->lj14.c12A,
                    iparams->lj14.c6B, iparams->lj14.c12B);
            break;
        case F_LJC14_Q:
            fprintf(fp, "fqq=%15.8e, qi=%15.8e, qj=%15.8e, c6=%15.8e, c12=%15.8e\n",
                    iparams->ljc14.fqq,
                    iparams->ljc14.qi, iparams->ljc14.qj,
                    iparams->ljc14.c6, iparams->ljc14.c12);
            break;
        case F_LJC_PAIRS_NB:
            fprintf(fp, "qi=%15.8e, qj=%15.8e, c6=%15.8e, c12=%15.8e\n",
                    iparams->ljcnb.qi, iparams->ljcnb.qj,
                    iparams->ljcnb.c6, iparams->ljcnb.c12);
            break;
        case F_PDIHS:
        case F_PIDIHS:
        case F_ANGRES:
        case F_ANGRESZ:
            fprintf(fp, "phiA=%15.8e, cpA=%15.8e, phiB=%15.8e, cpB=%15.8e, mult=%d\n",
                    iparams->pdihs.phiA, iparams->pdihs.cpA,
                    iparams->pdihs.phiB, iparams->pdihs.cpB,
                    iparams->pdihs.mult);
            break;
        case F_DISRES:
            fprintf(fp, "label=%4d, type=%1d, low=%15.8e, up1=%15.8e, up2=%15.8e, fac=%15.8e)\n",
                    iparams->disres.label, iparams->disres.type,
                    iparams->disres.low, iparams->disres.up1,
                    iparams->disres.up2, iparams->disres.kfac);
            break;
        case F_ORIRES:
            fprintf(fp, "ex=%4d, label=%d, power=%4d, c=%15.8e, obs=%15.8e, kfac=%15.8e)\n",
                    iparams->orires.ex, iparams->orires.label, iparams->orires.power,
                    iparams->orires.c, iparams->orires.obs, iparams->orires.kfac);
            break;
        case F_DIHRES:
            fprintf(fp, "phiA=%15.8e, dphiA=%15.8e, kfacA=%15.8e, phiB=%15.8e, dphiB=%15.8e, kfacB=%15.8e\n",
                    iparams->dihres.phiA, iparams->dihres.dphiA, iparams->dihres.kfacA,
                    iparams->dihres.phiB, iparams->dihres.dphiB, iparams->dihres.kfacB);
            break;
        case F_POSRES:
            fprintf(fp, "pos0A=(%15.8e,%15.8e,%15.8e), fcA=(%15.8e,%15.8e,%15.8e), pos0B=(%15.8e,%15.8e,%15.8e), fcB=(%15.8e,%15.8e,%15.8e)\n",
                    iparams->posres.pos0A[XX], iparams->posres.pos0A[YY],
                    iparams->posres.pos0A[ZZ], iparams->posres.fcA[XX],
                    iparams->posres.fcA[YY], iparams->posres.fcA[ZZ],
                    iparams->posres.pos0B[XX], iparams->posres.pos0B[YY],
                    iparams->posres.pos0B[ZZ], iparams->posres.fcB[XX],
                    iparams->posres.fcB[YY], iparams->posres.fcB[ZZ]);
            break;
        case F_FBPOSRES:
            fprintf(fp, "pos0=(%15.8e,%15.8e,%15.8e), geometry=%d, r=%15.8e, k=%15.8e\n",
                    iparams->fbposres.pos0[XX], iparams->fbposres.pos0[YY],
                    iparams->fbposres.pos0[ZZ], iparams->fbposres.geom,
                    iparams->fbposres.r,        iparams->fbposres.k);
            break;
        case F_RBDIHS:
            for (i = 0; i < NR_RBDIHS; i++)
            {
                fprintf(fp, "%srbcA[%d]=%15.8e", i == 0 ? "" : ", ", i, iparams->rbdihs.rbcA[i]);
            }
            fprintf(fp, "\n");
            for (i = 0; i < NR_RBDIHS; i++)
            {
                fprintf(fp, "%srbcB[%d]=%15.8e", i == 0 ? "" : ", ", i, iparams->rbdihs.rbcB[i]);
            }
            fprintf(fp, "\n");
            break;
        case F_FOURDIHS:
            /* Use the OPLS -> Ryckaert-Bellemans formula backwards to get the
             * OPLS potential constants back.
             */
            rbcA = iparams->rbdihs.rbcA;
            rbcB = iparams->rbdihs.rbcB;

            VA[3] = -0.25*rbcA[4];
            VA[2] = -0.5*rbcA[3];
            VA[1] = 4.0*VA[3]-rbcA[2];
            VA[0] = 3.0*VA[2]-2.0*rbcA[1];

            VB[3] = -0.25*rbcB[4];
            VB[2] = -0.5*rbcB[3];
            VB[1] = 4.0*VB[3]-rbcB[2];
            VB[0] = 3.0*VB[2]-2.0*rbcB[1];

            for (i = 0; i < NR_FOURDIHS; i++)
            {
                fprintf(fp, "%sFourA[%d]=%15.8e", i == 0 ? "" : ", ", i, VA[i]);
            }
            fprintf(fp, "\n");
            for (i = 0; i < NR_FOURDIHS; i++)
            {
                fprintf(fp, "%sFourB[%d]=%15.8e", i == 0 ? "" : ", ", i, VB[i]);
            }
            fprintf(fp, "\n");
            break;

        case F_CONSTR:
        case F_CONSTRNC:
            fprintf(fp, "dA=%15.8e, dB=%15.8e\n", iparams->constr.dA, iparams->constr.dB);
            break;
        case F_SETTLE:
            fprintf(fp, "doh=%15.8e, dhh=%15.8e\n", iparams->settle.doh,
                    iparams->settle.dhh);
            break;
        case F_VSITE2:
            fprintf(fp, "a=%15.8e\n", iparams->vsite.a);
            break;
        case F_VSITE3:
        case F_VSITE3FD:
        case F_VSITE3FAD:
            fprintf(fp, "a=%15.8e, b=%15.8e\n", iparams->vsite.a, iparams->vsite.b);
            break;
        case F_VSITE3OUT:
        case F_VSITE4FD:
        case F_VSITE4FDN:
            fprintf(fp, "a=%15.8e, b=%15.8e, c=%15.8e\n",
                    iparams->vsite.a, iparams->vsite.b, iparams->vsite.c);
            break;
        case F_VSITEN:
            fprintf(fp, "n=%2d, a=%15.8e\n", iparams->vsiten.n, iparams->vsiten.a);
            break;
        case F_GB12:
        case F_GB13:
        case F_GB14:
            fprintf(fp, "sar=%15.8e, st=%15.8e, pi=%15.8e, gbr=%15.8e, bmlt=%15.8e\n", iparams->gb.sar, iparams->gb.st, iparams->gb.pi, iparams->gb.gbr, iparams->gb.bmlt);
            break;
        case F_CMAP:
            fprintf(fp, "cmapA=%1d, cmapB=%1d\n", iparams->cmap.cmapA, iparams->cmap.cmapB);
            break;
        case  F_RESTRANGLES:
            pr_harm(fp, iparams, "ktheta", "costheta0");
            break;
        case  F_RESTRDIHS:
            fprintf(fp, "phiA=%15.8e, cpA=%15.8e",
                    iparams->pdihs.phiA, iparams->pdihs.cpA);
            break;
        case  F_CBTDIHS:
            fprintf(fp, "kphi=%15.8e", iparams->cbtdihs.cbtcA[0]);
            for (i = 1; i < NR_CBTDIHS; i++)
            {
                fprintf(fp, ", cbtcA[%d]=%15.8e", i-1, iparams->cbtdihs.cbtcA[i]);
            }
            fprintf(fp, "\n");
            break;
        default:
            gmx_fatal(FARGS, "unknown function type %d (%s) in %s line %d",
                      ftype, interaction_function[ftype].name, __FILE__, __LINE__);
    }
}

void pr_ilist(FILE *fp, int indent, const char *title,
              t_functype *functype, t_ilist *ilist, gmx_bool bShowNumbers)
{
    int      i, j, k, type, ftype;
    t_iatom *iatoms;

    if (available(fp, ilist, indent, title) && ilist->nr > 0)
    {
        indent = pr_title(fp, indent, title);
        (void) pr_indent(fp, indent);
        fprintf(fp, "nr: %d\n", ilist->nr);
        if (ilist->nr > 0)
        {
            (void) pr_indent(fp, indent);
            fprintf(fp, "iatoms:\n");
            iatoms = ilist->iatoms;
            for (i = j = 0; i < ilist->nr; )
            {
#ifndef DEBUG
                (void) pr_indent(fp, indent+INDENT);
                type  = *(iatoms++);
                ftype = functype[type];
                (void) fprintf(fp, "%d type=%d (%s)",
                               bShowNumbers ? j : -1, bShowNumbers ? type : -1,
                               interaction_function[ftype].name);
                j++;
                for (k = 0; k < interaction_function[ftype].nratoms; k++)
                {
                    (void) fprintf(fp, " %d", *(iatoms++));
                }
                (void) fprintf(fp, "\n");
                i += 1+interaction_function[ftype].nratoms;
#else
                fprintf(fp, "%5d%5d\n", i, iatoms[i]);
                i++;
#endif
            }
        }
    }
}

static void pr_cmap(FILE *fp, int indent, const char *title,
                    gmx_cmap_t *cmap_grid, gmx_bool bShowNumbers)
{
    int  i, j, nelem;
    real dx, idx;

    dx    = 360.0 / cmap_grid->grid_spacing;
    nelem = cmap_grid->grid_spacing*cmap_grid->grid_spacing;

    if (available(fp, cmap_grid, indent, title))
    {
        fprintf(fp, "%s\n", title);

        for (i = 0; i < cmap_grid->ngrid; i++)
        {
            idx = -180.0;
            fprintf(fp, "%8s %8s %8s %8s\n", "V", "dVdx", "dVdy", "d2dV");

            fprintf(fp, "grid[%3d]={\n", bShowNumbers ? i : -1);

            for (j = 0; j < nelem; j++)
            {
                if ( (j%cmap_grid->grid_spacing) == 0)
                {
                    fprintf(fp, "%8.1f\n", idx);
                    idx += dx;
                }

                fprintf(fp, "%8.3f ", cmap_grid->cmapdata[i].cmap[j*4]);
                fprintf(fp, "%8.3f ", cmap_grid->cmapdata[i].cmap[j*4+1]);
                fprintf(fp, "%8.3f ", cmap_grid->cmapdata[i].cmap[j*4+2]);
                fprintf(fp, "%8.3f\n", cmap_grid->cmapdata[i].cmap[j*4+3]);
            }
            fprintf(fp, "\n");
        }
    }

}

void pr_ffparams(FILE *fp, int indent, const char *title,
                 gmx_ffparams_t *ffparams,
                 gmx_bool bShowNumbers)
{
    int i, j;

    indent = pr_title(fp, indent, title);
    (void) pr_indent(fp, indent);
    (void) fprintf(fp, "atnr=%d\n", ffparams->atnr);
    (void) pr_indent(fp, indent);
    (void) fprintf(fp, "ntypes=%d\n", ffparams->ntypes);
    for (i = 0; i < ffparams->ntypes; i++)
    {
        (void) pr_indent(fp, indent+INDENT);
        (void) fprintf(fp, "functype[%d]=%s, ",
                       bShowNumbers ? i : -1,
                       interaction_function[ffparams->functype[i]].name);
        pr_iparams(fp, ffparams->functype[i], &ffparams->iparams[i]);
    }
    (void) pr_double(fp, indent, "reppow", ffparams->reppow);
    (void) pr_real(fp, indent, "fudgeQQ", ffparams->fudgeQQ);
    pr_cmap(fp, indent, "cmap", &ffparams->cmap_grid, bShowNumbers);
}

void pr_idef(FILE *fp, int indent, const char *title, t_idef *idef, gmx_bool bShowNumbers)
{
    int i, j;

    if (available(fp, idef, indent, title))
    {
        indent = pr_title(fp, indent, title);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "atnr=%d\n", idef->atnr);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "ntypes=%d\n", idef->ntypes);
        for (i = 0; i < idef->ntypes; i++)
        {
            (void) pr_indent(fp, indent+INDENT);
            (void) fprintf(fp, "functype[%d]=%s, ",
                           bShowNumbers ? i : -1,
                           interaction_function[idef->functype[i]].name);
            pr_iparams(fp, idef->functype[i], &idef->iparams[i]);
        }
        (void) pr_real(fp, indent, "fudgeQQ", idef->fudgeQQ);

        for (j = 0; (j < F_NRE); j++)
        {
            pr_ilist(fp, indent, interaction_function[j].longname,
                     idef->functype, &idef->il[j], bShowNumbers);
        }
    }
}

static int pr_block_title(FILE *fp, int indent, const char *title, t_block *block)
{
    int i;

    if (available(fp, block, indent, title))
    {
        indent = pr_title(fp, indent, title);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "nr=%d\n", block->nr);
    }
    return indent;
}

static int pr_blocka_title(FILE *fp, int indent, const char *title, t_blocka *block)
{
    int i;

    if (available(fp, block, indent, title))
    {
        indent = pr_title(fp, indent, title);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "nr=%d\n", block->nr);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "nra=%d\n", block->nra);
    }
    return indent;
}

static void low_pr_blocka(FILE *fp, int indent, const char *title, t_blocka *block, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, block, indent, title))
    {
        indent = pr_blocka_title(fp, indent, title, block);
        for (i = 0; i <= block->nr; i++)
        {
            (void) pr_indent(fp, indent+INDENT);
            (void) fprintf(fp, "%s->index[%d]=%d\n",
                           title, bShowNumbers ? i : -1, block->index[i]);
        }
        for (i = 0; i < block->nra; i++)
        {
            (void) pr_indent(fp, indent+INDENT);
            (void) fprintf(fp, "%s->a[%d]=%d\n",
                           title, bShowNumbers ? i : -1, block->a[i]);
        }
    }
}

void pr_block(FILE *fp, int indent, const char *title, t_block *block, gmx_bool bShowNumbers)
{
    int i, j, ok, size, start, end;

    if (available(fp, block, indent, title))
    {
        indent = pr_block_title(fp, indent, title, block);
        start  = 0;
        end    = start;
        if ((ok = (block->index[start] == 0)) == 0)
        {
            (void) fprintf(fp, "block->index[%d] should be 0\n", start);
        }
        else
        {
            for (i = 0; i < block->nr; i++)
            {
                end  = block->index[i+1];
                size = pr_indent(fp, indent);
                if (end <= start)
                {
                    size += fprintf(fp, "%s[%d]={}\n", title, i);
                }
                else
                {
                    size += fprintf(fp, "%s[%d]={%d..%d}\n",
                                    title, bShowNumbers ? i : -1,
                                    bShowNumbers ? start : -1, bShowNumbers ? end-1 : -1);
                }
                start = end;
            }
        }
    }
}

void pr_blocka(FILE *fp, int indent, const char *title, t_blocka *block, gmx_bool bShowNumbers)
{
    int i, j, ok, size, start, end;

    if (available(fp, block, indent, title))
    {
        indent = pr_blocka_title(fp, indent, title, block);
        start  = 0;
        end    = start;
        if ((ok = (block->index[start] == 0)) == 0)
        {
            (void) fprintf(fp, "block->index[%d] should be 0\n", start);
        }
        else
        {
            for (i = 0; i < block->nr; i++)
            {
                end  = block->index[i+1];
                size = pr_indent(fp, indent);
                if (end <= start)
                {
                    size += fprintf(fp, "%s[%d]={", title, i);
                }
                else
                {
                    size += fprintf(fp, "%s[%d][%d..%d]={",
                                    title, bShowNumbers ? i : -1,
                                    bShowNumbers ? start : -1, bShowNumbers ? end-1 : -1);
                }
                for (j = start; j < end; j++)
                {
                    if (j > start)
                    {
                        size += fprintf(fp, ", ");
                    }
                    if ((size) > (USE_WIDTH))
                    {
                        (void) fprintf(fp, "\n");
                        size = pr_indent(fp, indent+INDENT);
                    }
                    size += fprintf(fp, "%d", block->a[j]);
                }
                (void) fprintf(fp, "}\n");
                start = end;
            }
        }
        if ((end != block->nra) || (!ok))
        {
            (void) pr_indent(fp, indent);
            (void) fprintf(fp, "tables inconsistent, dumping complete tables:\n");
            low_pr_blocka(fp, indent, title, block, bShowNumbers);
        }
    }
}

static void pr_strings(FILE *fp, int indent, const char *title, char ***nm, int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, nm, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            (void) pr_indent(fp, indent);
            (void) fprintf(fp, "%s[%d]={name=\"%s\"}\n",
                           title, bShowNumbers ? i : -1, *(nm[i]));
        }
    }
}

static void pr_strings2(FILE *fp, int indent, const char *title,
                        char ***nm, char ***nmB, int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, nm, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            (void) pr_indent(fp, indent);
            (void) fprintf(fp, "%s[%d]={name=\"%s\",nameB=\"%s\"}\n",
                           title, bShowNumbers ? i : -1, *(nm[i]), *(nmB[i]));
        }
    }
}

static void pr_resinfo(FILE *fp, int indent, const char *title, t_resinfo *resinfo, int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, resinfo, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            (void) pr_indent(fp, indent);
            (void) fprintf(fp, "%s[%d]={name=\"%s\", nr=%d, ic='%c'}\n",
                           title, bShowNumbers ? i : -1,
                           *(resinfo[i].name), resinfo[i].nr,
                           (resinfo[i].ic == '\0') ? ' ' : resinfo[i].ic);
        }
    }
}

static void pr_atom(FILE *fp, int indent, const char *title, t_atom *atom, int n)
{
    int i, j;

    if (available(fp, atom, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            (void) pr_indent(fp, indent);
            fprintf(fp, "%s[%6d]={type=%3d, typeB=%3d, ptype=%8s, m=%12.5e, "
                    "q=%12.5e, mB=%12.5e, qB=%12.5e, resind=%5d, atomnumber=%3d}\n",
                    title, i, atom[i].type, atom[i].typeB, ptype_str[atom[i].ptype],
                    atom[i].m, atom[i].q, atom[i].mB, atom[i].qB,
                    atom[i].resind, atom[i].atomnumber);
        }
    }
}

static void pr_grps(FILE *fp, const char *title, t_grps grps[], char **grpname[])
{
    int i, j;

    for (i = 0; (i < egcNR); i++)
    {
        fprintf(fp, "%s[%-12s] nr=%d, name=[", title, gtypes[i], grps[i].nr);
        for (j = 0; (j < grps[i].nr); j++)
        {
            fprintf(fp, " %s", *(grpname[grps[i].nm_ind[j]]));
        }
        fprintf(fp, "]\n");
    }
}

static void pr_groups(FILE *fp, int indent,
                      gmx_groups_t *groups,
                      gmx_bool bShowNumbers)
{
    int grpnr[egcNR];
    int nat_max, i, g;

    pr_grps(fp, "grp", groups->grps, groups->grpname);
    pr_strings(fp, indent, "grpname", groups->grpname, groups->ngrpname, bShowNumbers);

    (void) pr_indent(fp, indent);
    fprintf(fp, "groups          ");
    for (g = 0; g < egcNR; g++)
    {
        printf(" %5.5s", gtypes[g]);
    }
    printf("\n");

    (void) pr_indent(fp, indent);
    fprintf(fp, "allocated       ");
    nat_max = 0;
    for (g = 0; g < egcNR; g++)
    {
        printf(" %5d", groups->ngrpnr[g]);
        nat_max = max(nat_max, groups->ngrpnr[g]);
    }
    printf("\n");

    if (nat_max == 0)
    {
        (void) pr_indent(fp, indent);
        fprintf(fp, "groupnr[%5s] =", "*");
        for (g = 0; g < egcNR; g++)
        {
            fprintf(fp, "  %3d ", 0);
        }
        fprintf(fp, "\n");
    }
    else
    {
        for (i = 0; i < nat_max; i++)
        {
            (void) pr_indent(fp, indent);
            fprintf(fp, "groupnr[%5d] =", i);
            for (g = 0; g < egcNR; g++)
            {
                fprintf(fp, "  %3d ",
                        groups->grpnr[g] ? groups->grpnr[g][i] : 0);
            }
            fprintf(fp, "\n");
        }
    }
}

void pr_atoms(FILE *fp, int indent, const char *title, t_atoms *atoms,
              gmx_bool bShownumbers)
{
    if (available(fp, atoms, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_atom(fp, indent, "atom", atoms->atom, atoms->nr);
        pr_strings(fp, indent, "atom", atoms->atomname, atoms->nr, bShownumbers);
        pr_strings2(fp, indent, "type", atoms->atomtype, atoms->atomtypeB, atoms->nr, bShownumbers);
        pr_resinfo(fp, indent, "residue", atoms->resinfo, atoms->nres, bShownumbers);
    }
}


void pr_atomtypes(FILE *fp, int indent, const char *title, t_atomtypes *atomtypes,
                  gmx_bool bShowNumbers)
{
    int i;
    if (available(fp, atomtypes, indent, title))
    {
        indent = pr_title(fp, indent, title);
        for (i = 0; i < atomtypes->nr; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp,
                    "atomtype[%3d]={radius=%12.5e, volume=%12.5e, gb_radius=%12.5e, surftens=%12.5e, atomnumber=%4d, S_hct=%12.5e)}\n",
                    bShowNumbers ? i : -1, atomtypes->radius[i], atomtypes->vol[i],
                    atomtypes->gb_radius[i],
                    atomtypes->surftens[i], atomtypes->atomnumber[i], atomtypes->S_hct[i]);
        }
    }
}

static void pr_moltype(FILE *fp, int indent, const char *title,
                       gmx_moltype_t *molt, int n,
                       gmx_ffparams_t *ffparams,
                       gmx_bool bShowNumbers)
{
    int j;

    indent = pr_title_n(fp, indent, title, n);
    (void) pr_indent(fp, indent);
    (void) fprintf(fp, "name=\"%s\"\n", *(molt->name));
    pr_atoms(fp, indent, "atoms", &(molt->atoms), bShowNumbers);
    pr_block(fp, indent, "cgs", &molt->cgs, bShowNumbers);
    pr_blocka(fp, indent, "excls", &molt->excls, bShowNumbers);
    for (j = 0; (j < F_NRE); j++)
    {
        pr_ilist(fp, indent, interaction_function[j].longname,
                 ffparams->functype, &molt->ilist[j], bShowNumbers);
    }
}

static void pr_molblock(FILE *fp, int indent, const char *title,
                        gmx_molblock_t *molb, int n,
                        gmx_moltype_t *molt)
{
    indent = pr_title_n(fp, indent, title, n);
    (void) pr_indent(fp, indent);
    (void) fprintf(fp, "%-20s = %d \"%s\"\n",
                   "moltype", molb->type, *(molt[molb->type].name));
    pr_int(fp, indent, "#molecules", molb->nmol);
    pr_int(fp, indent, "#atoms_mol", molb->natoms_mol);
    pr_int(fp, indent, "#posres_xA", molb->nposres_xA);
    if (molb->nposres_xA > 0)
    {
        pr_rvecs(fp, indent, "posres_xA", molb->posres_xA, molb->nposres_xA);
    }
    pr_int(fp, indent, "#posres_xB", molb->nposres_xB);
    if (molb->nposres_xB > 0)
    {
        pr_rvecs(fp, indent, "posres_xB", molb->posres_xB, molb->nposres_xB);
    }
}

void pr_mtop(FILE *fp, int indent, const char *title, gmx_mtop_t *mtop,
             gmx_bool bShowNumbers)
{
    int mt, mb, j;

    if (available(fp, mtop, indent, title))
    {
        indent = pr_title(fp, indent, title);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "name=\"%s\"\n", *(mtop->name));
        pr_int(fp, indent, "#atoms", mtop->natoms);
        pr_int(fp, indent, "#molblock", mtop->nmolblock);
        for (mb = 0; mb < mtop->nmolblock; mb++)
        {
            pr_molblock(fp, indent, "molblock", &mtop->molblock[mb], mb, mtop->moltype);
        }
        pr_str(fp, indent, "bIntermolecularInteractions", EBOOL(mtop->bIntermolecularInteractions));
        if (mtop->bIntermolecularInteractions)
        {
            for (j = 0; (j < F_NRE); j++)
            {
                pr_ilist(fp, indent, interaction_function[j].longname,
                         mtop->ffparams.functype,
                         &mtop->intermolecular_ilist[j], bShowNumbers);
            }
        }
        pr_ffparams(fp, indent, "ffparams", &(mtop->ffparams), bShowNumbers);
        pr_atomtypes(fp, indent, "atomtypes", &(mtop->atomtypes), bShowNumbers);
        for (mt = 0; mt < mtop->nmoltype; mt++)
        {
            pr_moltype(fp, indent, "moltype", &mtop->moltype[mt], mt,
                       &mtop->ffparams, bShowNumbers);
        }
        pr_groups(fp, indent, &mtop->groups, bShowNumbers);
    }
}

void pr_top(FILE *fp, int indent, const char *title, t_topology *top, gmx_bool bShowNumbers)
{
    if (available(fp, top, indent, title))
    {
        indent = pr_title(fp, indent, title);
        (void) pr_indent(fp, indent);
        (void) fprintf(fp, "name=\"%s\"\n", *(top->name));
        pr_atoms(fp, indent, "atoms", &(top->atoms), bShowNumbers);
        pr_atomtypes(fp, indent, "atomtypes", &(top->atomtypes), bShowNumbers);
        pr_block(fp, indent, "cgs", &top->cgs, bShowNumbers);
        pr_block(fp, indent, "mols", &top->mols, bShowNumbers);
        pr_str(fp, indent, "bIntermolecularInteractions", EBOOL(top->bIntermolecularInteractions));
        pr_blocka(fp, indent, "excls", &top->excls, bShowNumbers);
        pr_idef(fp, indent, "idef", &top->idef, bShowNumbers);
    }
}

void pr_header(FILE *fp, int indent, const char *title, t_tpxheader *sh)
{
    char buf[22];

    if (available(fp, sh, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_indent(fp, indent);
        fprintf(fp, "bIr    = %spresent\n", sh->bIr ? "" : "not ");
        pr_indent(fp, indent);
        fprintf(fp, "bBox   = %spresent\n", sh->bBox ? "" : "not ");
        pr_indent(fp, indent);
        fprintf(fp, "bTop   = %spresent\n", sh->bTop ? "" : "not ");
        pr_indent(fp, indent);
        fprintf(fp, "bX     = %spresent\n", sh->bX ? "" : "not ");
        pr_indent(fp, indent);
        fprintf(fp, "bV     = %spresent\n", sh->bV ? "" : "not ");
        pr_indent(fp, indent);
        fprintf(fp, "bF     = %spresent\n", sh->bF ? "" : "not ");

        pr_indent(fp, indent);
        fprintf(fp, "natoms = %d\n", sh->natoms);
        pr_indent(fp, indent);
        fprintf(fp, "lambda = %e\n", sh->lambda);
    }
}

void pr_commrec(FILE *fp, int indent, t_commrec *cr)
{
    pr_indent(fp, indent);
    fprintf(fp, "commrec:\n");
    indent += 2;
    pr_indent(fp, indent);
    fprintf(fp, "rank      = %d\n", cr->nodeid);
    pr_indent(fp, indent);
    fprintf(fp, "number of ranks = %d\n", cr->nnodes);
    pr_indent(fp, indent);
    fprintf(fp, "PME-only ranks = %d\n", cr->npmenodes);
    /*
       pr_indent(fp,indent);
       fprintf(fp,"threadid  = %d\n",cr->threadid);
       pr_indent(fp,indent);
       fprintf(fp,"nthreads  = %d\n",cr->nthreads);
     */
}

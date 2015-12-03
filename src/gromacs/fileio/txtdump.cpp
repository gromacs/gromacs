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

#include "txtdump.h"

#include <cstdio>
#include <cstdlib>

#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/stringutil.h"

//! Macro to select a bool name
#define EBOOL(e)       gmx::boolToString(e)

int pr_indent(FILE *fp, int n)
{
    int i;

    for (i = 0; i < n; i++)
    {
        fprintf(fp, " ");
    }
    return n;
}

int available(FILE *fp, const void *p, int indent, const char *title)
{
    if (!p)
    {
        if (indent > 0)
        {
            pr_indent(fp, indent);
        }
        fprintf(fp, "%s: not available\n", title);
    }
    return (p != NULL);
}

int pr_title(FILE *fp, int indent, const char *title)
{
    pr_indent(fp, indent);
    fprintf(fp, "%s:\n", title);
    return (indent+INDENT);
}

int pr_title_n(FILE *fp, int indent, const char *title, int n)
{
    pr_indent(fp, indent);
    fprintf(fp, "%s (%d):\n", title, n);
    return (indent+INDENT);
}

int pr_title_nxn(FILE *fp, int indent, const char *title, int n1, int n2)
{
    pr_indent(fp, indent);
    fprintf(fp, "%s (%dx%d):\n", title, n1, n2);
    return (indent+INDENT);
}

void pr_ivec(FILE *fp, int indent, const char *title, const int vec[], int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]=%d\n", title, bShowNumbers ? i : -1, vec[i]);
        }
    }
}

void pr_ivec_block(FILE *fp, int indent, const char *title, const int vec[], int n, gmx_bool bShowNumbers)
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
                    pr_indent(fp, indent);
                    fprintf(fp, "%s[%d]=%d\n",
                            title, bShowNumbers ? i : -1, vec[i]);
                    i++;
                }
            }
            else
            {
                pr_indent(fp, indent);
                fprintf(fp, "%s[%d,...,%d] = {%d,...,%d}\n",
                        title,
                        bShowNumbers ? i : -1,
                        bShowNumbers ? j-1 : -1,
                        vec[i], vec[j-1]);
                i = j;
            }
        }
    }
}

void pr_bvec(FILE *fp, int indent, const char *title, const gmx_bool vec[], int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]=%s\n", title, bShowNumbers ? i : -1,
                    EBOOL(vec[i]));
        }
    }
}

void pr_ivecs(FILE *fp, int indent, const char *title, const ivec vec[], int n, gmx_bool bShowNumbers)
{
    int i, j;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, DIM);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]={", title, bShowNumbers ? i : -1);
            for (j = 0; j < DIM; j++)
            {
                if (j != 0)
                {
                    fprintf(fp, ", ");
                }
                fprintf(fp, "%d", vec[i][j]);
            }
            fprintf(fp, "}\n");
        }
    }
}

void pr_rvec(FILE *fp, int indent, const char *title, const real vec[], int n, gmx_bool bShowNumbers)
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

void pr_dvec(FILE *fp, int indent, const char *title, const double vec[], int n, gmx_bool bShowNumbers)
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


void pr_rvecs_len(FILE *fp, int indent, const char *title, const rvec vec[], int n)
{
    int i, j;

    if (available(fp, vec, indent, title))
    {
        indent = pr_title_nxn(fp, indent, title, n, DIM);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%5d]={", title, i);
            for (j = 0; j < DIM; j++)
            {
                if (j != 0)
                {
                    fprintf(fp, ", ");
                }
                fprintf(fp, "%12.5e", vec[i][j]);
            }
            fprintf(fp, "} len=%12.5e\n", norm(vec[i]));
        }
    }
}

void pr_rvecs(FILE *fp, int indent, const char *title, const rvec vec[], int n)
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
            pr_indent(fp, indent);
            fprintf(fp, "%s[%5d]={", title, i);
            for (j = 0; j < DIM; j++)
            {
                if (j != 0)
                {
                    fprintf(fp, ", ");
                }
                fprintf(fp, format, vec[i][j]);
            }
            fprintf(fp, "}\n");
        }
    }
}


void pr_rvecs_of_dim(FILE *fp, int indent, const char *title, const rvec vec[], int n, int dim)
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
            pr_indent(fp, indent);
            fprintf(fp, "%s[%5d]={", title, i);
            for (j = 0; j < dim; j++)
            {
                if (j != 0)
                {
                    fprintf(fp, ", ");
                }
                fprintf(fp, format, vec[i][j]);
            }
            fprintf(fp, "}\n");
        }
    }
}

void pr_reals(FILE *fp, int indent, const char *title, const real *vec, int n)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        pr_indent(fp, indent);
        fprintf(fp, "%s:\t", title);
        for (i = 0; i < n; i++)
        {
            fprintf(fp, "  %10g", vec[i]);
        }
        fprintf(fp, "\n");
    }
}

void pr_doubles(FILE *fp, int indent, const char *title, const double *vec, int n)
{
    int i;

    if (available(fp, vec, indent, title))
    {
        pr_indent(fp, indent);
        fprintf(fp, "%s:\t", title);
        for (i = 0; i < n; i++)
        {
            fprintf(fp, "  %10g", vec[i]);
        }
        fprintf(fp, "\n");
    }
}

void pr_reals_of_dim(FILE *fp, int indent, const char *title, const real *vec, int n, int dim)
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
            pr_indent(fp, indent);
            fprintf(fp, "%s[%5d]={", title, i);
            for (j = 0; j < dim; j++)
            {
                if (j != 0)
                {
                    fprintf(fp, ", ");
                }
                fprintf(fp, format, vec[i * dim  + j]);
            }
            fprintf(fp, "}\n");
        }
    }
}

void pr_int(FILE *fp, int indent, const char *title, int i)
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

void pr_real(FILE *fp, int indent, const char *title, real r)
{
    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %g\n", title, r);
}

void pr_double(FILE *fp, int indent, const char *title, double d)
{
    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %g\n", title, d);
}

void pr_str(FILE *fp, int indent, const char *title, const char *s)
{
    pr_indent(fp, indent);
    fprintf(fp, "%-30s = %s\n", title, s);
}

void pr_qm_opts(FILE *fp, int indent, const char *title, const t_grpopts *opts)
{
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

static void pr_grp_opts(FILE *out, int indent, const char *title, const t_grpopts *opts,
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

static void pr_matrix(FILE *fp, int indent, const char *title, const rvec *m,
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

static void pr_cosine(FILE *fp, int indent, const char *title, const t_cosines *cos,
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
        pr_indent(fp, indent);
        fprintf(fp, "n = %d\n", cos->n);
        if (cos->n > 0)
        {
            pr_indent(fp, indent+2);
            fprintf(fp, "a =");
            for (j = 0; (j < cos->n); j++)
            {
                fprintf(fp, " %e", cos->a[j]);
            }
            fprintf(fp, "\n");
            pr_indent(fp, indent+2);
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

static void pr_pull_group(FILE *fp, int indent, int g, const t_pull_group *pgrp)
{
    pr_indent(fp, indent);
    fprintf(fp, "pull-group %d:\n", g);
    indent += 2;
    pr_ivec_block(fp, indent, "atom", pgrp->ind, pgrp->nat, TRUE);
    pr_rvec(fp, indent, "weight", pgrp->weight, pgrp->nweight, TRUE);
    PI("pbcatom", pgrp->pbcatom);
}

static void pr_pull_coord(FILE *fp, int indent, int c, const t_pull_coord *pcrd)
{
    int g;

    pr_indent(fp, indent);
    fprintf(fp, "pull-coord %d:\n", c);
    PS("type", EPULLTYPE(pcrd->eType));
    PS("geometry", EPULLGEOM(pcrd->eGeom));
    for (g = 0; g < pcrd->ngroup; g++)
    {
        char buf[10];

        sprintf(buf, "group[%d]", g);
        PI(buf, pcrd->group[g]);
    }
    pr_ivec(fp, indent, "dim", pcrd->dim, DIM, TRUE);
    pr_rvec(fp, indent, "origin", pcrd->origin, DIM, TRUE);
    pr_rvec(fp, indent, "vec", pcrd->vec, DIM, TRUE);
    PS("start", EBOOL(pcrd->bStart));
    PR("init", pcrd->init);
    PR("rate", pcrd->rate);
    PR("k", pcrd->k);
    PR("kB", pcrd->kB);
}

static void pr_simtempvals(FILE *fp, int indent, const t_simtemp *simtemp, int n_lambda)
{
    PS("simulated-tempering-scaling", ESIMTEMP(simtemp->eSimTempScale));
    PR("sim-temp-low", simtemp->simtemp_low);
    PR("sim-temp-high", simtemp->simtemp_high);
    pr_rvec(fp, indent, "simulated tempering temperatures", simtemp->temperatures, n_lambda, TRUE);
}

static void pr_expandedvals(FILE *fp, int indent, const t_expanded *expand, int n_lambda)
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

static void pr_fepvals(FILE *fp, int indent, const t_lambda *fep, gmx_bool bMDPformat)
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

static void pr_pull(FILE *fp, int indent, const pull_params_t *pull)
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

static void pr_rotgrp(FILE *fp, int indent, int g, const t_rotgrp *rotg)
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

static void pr_rot(FILE *fp, int indent, const t_rot *rot)
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


static void pr_swap(FILE *fp, int indent, const t_swapcoords *swap)
{
    char str[STRLEN];

    /* Enums for better readability of the code */
    enum {
        eCompA = 0, eCompB
    };


    PI("swap-frequency", swap->nstswap);

    /* The split groups that define the compartments */
    for (int j = 0; j < 2; j++)
    {
        snprintf(str, STRLEN, "massw_split%d", j);
        PS(str, EBOOL(swap->massw_split[j]));
        snprintf(str, STRLEN, "split atoms group %d", j);
        pr_ivec_block(fp, indent, str, swap->grp[j].ind, swap->grp[j].nat, TRUE);
    }

    /* The solvent group */
    snprintf(str, STRLEN, "solvent group %s", swap->grp[eGrpSolvent].molname);
    pr_ivec_block(fp, indent, str, swap->grp[eGrpSolvent].ind, swap->grp[eGrpSolvent].nat, TRUE);

    /* Now print the indices for all the ion groups: */
    for (int ig = eSwapFixedGrpNR; ig < swap->ngrp; ig++)
    {
        snprintf(str, STRLEN, "ion group %s", swap->grp[ig].molname);
        pr_ivec_block(fp, indent, str, swap->grp[ig].ind, swap->grp[ig].nat, TRUE);
    }

    PR("cyl0-r", swap->cyl0r);
    PR("cyl0-up", swap->cyl0u);
    PR("cyl0-down", swap->cyl0l);
    PR("cyl1-r", swap->cyl1r);
    PR("cyl1-up", swap->cyl1u);
    PR("cyl1-down", swap->cyl1l);
    PI("coupl-steps", swap->nAverage);

    /* Print the requested ion counts for both compartments */
    for (int ic = eCompA; ic <= eCompB; ic++)
    {
        for (int ig = eSwapFixedGrpNR; ig < swap->ngrp; ig++)
        {
            snprintf(str, STRLEN, "%s-in-%c", swap->grp[ig].molname, 'A'+ic);
            PI(str, swap->grp[ig].nmolReq[ic]);
        }
    }

    PR("threshold", swap->threshold);
    PR("bulk-offsetA", swap->bulkOffset[eCompA]);
    PR("bulk-offsetB", swap->bulkOffset[eCompB]);
}


static void pr_imd(FILE *fp, int indent, const t_IMD *imd)
{
    PI("IMD-atoms", imd->nat);
    pr_ivec_block(fp, indent, "atom", imd->ind, imd->nat, TRUE);
}


void pr_inputrec(FILE *fp, int indent, const char *title, const t_inputrec *ir,
                 gmx_bool bMDPformat)
{
    const char *infbuf = "inf";

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
        PS("pbc", epbc_names[ir->ePBC]);
        PS("periodic-molecules", EBOOL(ir->bPeriodicMols));
        PR("verlet-buffer-tolerance", ir->verletbuf_tol);
        PR("rlist", ir->rlist);

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

void pr_strings(FILE *fp, int indent, const char *title, char ***nm, int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, nm, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]={name=\"%s\"}\n",
                    title, bShowNumbers ? i : -1, *(nm[i]));
        }
    }
}

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

#include "gromacs/commandline/pargs.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/mdebin.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static real       minthird = -1.0/3.0, minsixth = -1.0/6.0;

typedef struct {
    real sum;
    real sum2;
} exactsum_t;

typedef struct {
    real       *ener;
    exactsum_t *es;
    gmx_bool    bExactStat;
    double      av;
    double      rmsd;
    double      ee;
    double      slope;
} enerdat_t;

typedef struct {
    gmx_int64_t      nsteps;
    gmx_int64_t      npoints;
    int              nframes;
    int             *step;
    int             *steps;
    int             *points;
    enerdat_t       *s;
    gmx_bool         bHaveSums;
} enerdata_t;

static double mypow(double x, double y)
{
    if (x > 0)
    {
        return pow(x, y);
    }
    else
    {
        return 0.0;
    }
}

static int *select_it(int nre, char *nm[], int *nset)
{
    gmx_bool *bE;
    int       n, k, j, i;
    int      *set;
    gmx_bool  bVerbose = TRUE;

    if ((getenv("GMX_ENER_VERBOSE")) != NULL)
    {
        bVerbose = FALSE;
    }

    fprintf(stderr, "Select the terms you want from the following list\n");
    fprintf(stderr, "End your selection with 0\n");

    if (bVerbose)
    {
        for (k = 0; (k < nre); )
        {
            for (j = 0; (j < 4) && (k < nre); j++, k++)
            {
                fprintf(stderr, " %3d=%14s", k+1, nm[k]);
            }
            fprintf(stderr, "\n");
        }
    }

    snew(bE, nre);
    do
    {
        if (1 != scanf("%d", &n))
        {
            gmx_fatal(FARGS, "Error reading user input");
        }
        if ((n > 0) && (n <= nre))
        {
            bE[n-1] = TRUE;
        }
    }
    while (n != 0);

    snew(set, nre);
    for (i = (*nset) = 0; (i < nre); i++)
    {
        if (bE[i])
        {
            set[(*nset)++] = i;
        }
    }

    sfree(bE);

    return set;
}

static void chomp(char *buf)
{
    int len = strlen(buf);

    while ((len > 0) && (buf[len-1] == '\n'))
    {
        buf[len-1] = '\0';
        len--;
    }
}

static int *select_by_name(int nre, gmx_enxnm_t *nm, int *nset)
{
    gmx_bool   *bE;
    int         n, k, kk, j, i, nmatch, nind, nss;
    int        *set;
    gmx_bool    bEOF, bVerbose = TRUE, bLong = FALSE;
    char       *ptr, buf[STRLEN];
    const char *fm4   = "%3d  %-14s";
    const char *fm2   = "%3d  %-34s";
    char      **newnm = NULL;

    if ((getenv("GMX_ENER_VERBOSE")) != NULL)
    {
        bVerbose = FALSE;
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "Select the terms you want from the following list by\n");
    fprintf(stderr, "selecting either (part of) the name or the number or a combination.\n");
    fprintf(stderr, "End your selection with an empty line or a zero.\n");
    fprintf(stderr, "-------------------------------------------------------------------\n");

    snew(newnm, nre);
    j = 0;
    for (k = 0; k < nre; k++)
    {
        newnm[k] = gmx_strdup(nm[k].name);
        /* Insert dashes in all the names */
        while ((ptr = strchr(newnm[k], ' ')) != NULL)
        {
            *ptr = '-';
        }
        if (bVerbose)
        {
            if (j == 0)
            {
                if (k > 0)
                {
                    fprintf(stderr, "\n");
                }
                bLong = FALSE;
                for (kk = k; kk < k+4; kk++)
                {
                    if (kk < nre && strlen(nm[kk].name) > 14)
                    {
                        bLong = TRUE;
                    }
                }
            }
            else
            {
                fprintf(stderr, " ");
            }
            if (!bLong)
            {
                fprintf(stderr, fm4, k+1, newnm[k]);
                j++;
                if (j == 4)
                {
                    j = 0;
                }
            }
            else
            {
                fprintf(stderr, fm2, k+1, newnm[k]);
                j++;
                if (j == 2)
                {
                    j = 0;
                }
            }
        }
    }
    if (bVerbose)
    {
        fprintf(stderr, "\n\n");
    }

    snew(bE, nre);

    bEOF = FALSE;
    while (!bEOF && (fgets2(buf, STRLEN-1, stdin)))
    {
        /* Remove newlines */
        chomp(buf);

        /* Remove spaces */
        trim(buf);

        /* Empty line means end of input */
        bEOF = (strlen(buf) == 0);
        if (!bEOF)
        {
            ptr = buf;
            do
            {
                if (!bEOF)
                {
                    /* First try to read an integer */
                    nss   = sscanf(ptr, "%d", &nind);
                    if (nss == 1)
                    {
                        /* Zero means end of input */
                        if (nind == 0)
                        {
                            bEOF = TRUE;
                        }
                        else if ((1 <= nind) && (nind <= nre))
                        {
                            bE[nind-1] = TRUE;
                        }
                        else
                        {
                            fprintf(stderr, "number %d is out of range\n", nind);
                        }
                    }
                    else
                    {
                        /* Now try to read a string */
                        i      = strlen(ptr);
                        nmatch = 0;
                        for (nind = 0; nind < nre; nind++)
                        {
                            if (gmx_strcasecmp(newnm[nind], ptr) == 0)
                            {
                                bE[nind] = TRUE;
                                nmatch++;
                            }
                        }
                        if (nmatch == 0)
                        {
                            i      = strlen(ptr);
                            nmatch = 0;
                            for (nind = 0; nind < nre; nind++)
                            {
                                if (gmx_strncasecmp(newnm[nind], ptr, i) == 0)
                                {
                                    bE[nind] = TRUE;
                                    nmatch++;
                                }
                            }
                            if (nmatch == 0)
                            {
                                fprintf(stderr, "String '%s' does not match anything\n", ptr);
                            }
                        }
                    }
                }
                /* Look for the first space, and remove spaces from there */
                if ((ptr = strchr(ptr, ' ')) != NULL)
                {
                    trim(ptr);
                }
            }
            while (!bEOF && (ptr && (strlen(ptr) > 0)));
        }
    }

    snew(set, nre);
    for (i = (*nset) = 0; (i < nre); i++)
    {
        if (bE[i])
        {
            set[(*nset)++] = i;
        }
    }

    sfree(bE);

    if (*nset == 0)
    {
        gmx_fatal(FARGS, "No energy terms selected");
    }

    for (i = 0; (i < nre); i++)
    {
        sfree(newnm[i]);
    }
    sfree(newnm);

    return set;
}

static void get_dhdl_parms(const char *topnm, t_inputrec *ir)
{
    gmx_mtop_t  mtop;
    int         natoms;
    t_iatom    *iatom;
    matrix      box;

    /* all we need is the ir to be able to write the label */
    read_tpx(topnm, ir, box, &natoms, NULL, NULL, NULL, &mtop);
}

static void get_orires_parms(const char *topnm,
                             int *nor, int *nex, int **label, real **obs)
{
    gmx_mtop_t      mtop;
    gmx_localtop_t *top;
    t_inputrec      ir;
    t_iparams      *ip;
    int             natoms, i;
    t_iatom        *iatom;
    int             nb;
    matrix          box;

    read_tpx(topnm, &ir, box, &natoms, NULL, NULL, NULL, &mtop);
    top = gmx_mtop_generate_local_top(&mtop, &ir);

    ip       = top->idef.iparams;
    iatom    = top->idef.il[F_ORIRES].iatoms;

    /* Count how many distance restraint there are... */
    nb = top->idef.il[F_ORIRES].nr;
    if (nb == 0)
    {
        gmx_fatal(FARGS, "No orientation restraints in topology!\n");
    }

    *nor = nb/3;
    *nex = 0;
    snew(*label, *nor);
    snew(*obs, *nor);
    for (i = 0; i < nb; i += 3)
    {
        (*label)[i/3] = ip[iatom[i]].orires.label;
        (*obs)[i/3]   = ip[iatom[i]].orires.obs;
        if (ip[iatom[i]].orires.ex >= *nex)
        {
            *nex = ip[iatom[i]].orires.ex+1;
        }
    }
    fprintf(stderr, "Found %d orientation restraints with %d experiments",
            *nor, *nex);
}

static int get_bounds(const char *topnm,
                      real **bounds, int **index, int **dr_pair, int *npairs,
                      gmx_mtop_t *mtop, gmx_localtop_t **ltop, t_inputrec *ir)
{
    gmx_localtop_t *top;
    t_functype     *functype;
    t_iparams      *ip;
    int             natoms, i, j, k, type, ftype, natom;
    t_ilist        *disres;
    t_iatom        *iatom;
    real           *b;
    int            *ind, *pair;
    int             nb, label1;
    matrix          box;

    read_tpx(topnm, ir, box, &natoms, NULL, NULL, NULL, mtop);
    snew(*ltop, 1);
    top   = gmx_mtop_generate_local_top(mtop, ir);
    *ltop = top;

    functype = top->idef.functype;
    ip       = top->idef.iparams;

    /* Count how many distance restraint there are... */
    nb = top->idef.il[F_DISRES].nr;
    if (nb == 0)
    {
        gmx_fatal(FARGS, "No distance restraints in topology!\n");
    }

    /* Allocate memory */
    snew(b, nb);
    snew(ind, nb);
    snew(pair, nb+1);

    /* Fill the bound array */
    nb = 0;
    for (i = 0; (i < top->idef.ntypes); i++)
    {
        ftype = functype[i];
        if (ftype == F_DISRES)
        {

            label1  = ip[i].disres.label;
            b[nb]   = ip[i].disres.up1;
            ind[nb] = label1;
            nb++;
        }
    }
    *bounds = b;

    /* Fill the index array */
    label1  = -1;
    disres  = &(top->idef.il[F_DISRES]);
    iatom   = disres->iatoms;
    for (i = j = k = 0; (i < disres->nr); )
    {
        type  = iatom[i];
        ftype = top->idef.functype[type];
        natom = interaction_function[ftype].nratoms+1;
        if (label1 != top->idef.iparams[type].disres.label)
        {
            pair[j] = k;
            label1  = top->idef.iparams[type].disres.label;
            j++;
        }
        k++;
        i += natom;
    }
    pair[j]  = k;
    *npairs  = k;
    if (j != nb)
    {
        gmx_incons("get_bounds for distance restraints");
    }

    *index   = ind;
    *dr_pair = pair;

    return nb;
}

static void calc_violations(real rt[], real rav3[], int nb, int index[],
                            real bounds[], real *viol, double *st, double *sa)
{
    const   real sixth = 1.0/6.0;
    int          i, j;
    double       rsum, rav, sumaver, sumt;

    sumaver = 0;
    sumt    = 0;
    for (i = 0; (i < nb); i++)
    {
        rsum = 0.0;
        rav  = 0.0;
        for (j = index[i]; (j < index[i+1]); j++)
        {
            if (viol)
            {
                viol[j] += mypow(rt[j], -3.0);
            }
            rav     += sqr(rav3[j]);
            rsum    += mypow(rt[j], -6);
        }
        rsum    = max(0.0, mypow(rsum, -sixth)-bounds[i]);
        rav     = max(0.0, mypow(rav, -sixth)-bounds[i]);

        sumt    += rsum;
        sumaver += rav;
    }
    *st = sumt;
    *sa = sumaver;
}

static void analyse_disre(const char *voutfn,    int nframes,
                          real violaver[], real bounds[], int index[],
                          int pair[],      int nbounds,
                          const output_env_t oenv)
{
    FILE   *vout;
    double  sum, sumt, sumaver;
    int     i, j;

    /* Subtract bounds from distances, to calculate violations */
    calc_violations(violaver, violaver,
                    nbounds, pair, bounds, NULL, &sumt, &sumaver);

#ifdef DEBUG
    fprintf(stdout, "\nSum of violations averaged over simulation: %g nm\n",
            sumaver);
    fprintf(stdout, "Largest violation averaged over simulation: %g nm\n\n",
            sumt);
#endif
    vout = xvgropen(voutfn, "r\\S-3\\N average violations", "DR Index", "nm",
                    oenv);
    sum  = 0.0;
    sumt = 0.0;
    for (i = 0; (i < nbounds); i++)
    {
        /* Do ensemble averaging */
        sumaver = 0;
        for (j = pair[i]; (j < pair[i+1]); j++)
        {
            sumaver += sqr(violaver[j]/nframes);
        }
        sumaver = max(0.0, mypow(sumaver, minsixth)-bounds[i]);

        sumt   += sumaver;
        sum     = max(sum, sumaver);
        fprintf(vout, "%10d  %10.5e\n", index[i], sumaver);
    }
#ifdef DEBUG
    for (j = 0; (j < dr.ndr); j++)
    {
        fprintf(vout, "%10d  %10.5e\n", j, mypow(violaver[j]/nframes, minthird));
    }
#endif
    xvgrclose(vout);

    fprintf(stdout, "\nSum of violations averaged over simulation: %g nm\n",
            sumt);
    fprintf(stdout, "Largest violation averaged over simulation: %g nm\n\n", sum);

    do_view(oenv, voutfn, "-graphtype bar");
}

static void einstein_visco(const char *fn, const char *fni, int nsets,
                           int nint, real **eneint,
                           real V, real T, double dt,
                           const output_env_t oenv)
{
    FILE  *fp0, *fp1;
    double av[4], avold[4];
    double fac, di;
    int    i, j, m, nf4;

    nf4 = nint/4 + 1;

    for (i = 0; i <= nsets; i++)
    {
        avold[i] = 0;
    }
    fp0 = xvgropen(fni, "Shear viscosity integral",
                   "Time (ps)", "(kg m\\S-1\\N s\\S-1\\N ps)", oenv);
    fp1 = xvgropen(fn, "Shear viscosity using Einstein relation",
                   "Time (ps)", "(kg m\\S-1\\N s\\S-1\\N)", oenv);
    for (i = 0; i < nf4; i++)
    {
        for (m = 0; m <= nsets; m++)
        {
            av[m] = 0;
        }
        for (j = 0; j < nint - i; j++)
        {
            for (m = 0; m < nsets; m++)
            {
                di   = sqr(eneint[m][j+i] - eneint[m][j]);

                av[m]     += di;
                av[nsets] += di/nsets;
            }
        }
        /* Convert to SI for the viscosity */
        fac = (V*NANO*NANO*NANO*PICO*1e10)/(2*BOLTZMANN*T)/(nint - i);
        fprintf(fp0, "%10g", i*dt);
        for (m = 0; (m <= nsets); m++)
        {
            av[m] = fac*av[m];
            fprintf(fp0, "  %10g", av[m]);
        }
        fprintf(fp0, "\n");
        fprintf(fp1, "%10g", (i + 0.5)*dt);
        for (m = 0; (m <= nsets); m++)
        {
            fprintf(fp1, "  %10g", (av[m]-avold[m])/dt);
            avold[m] = av[m];
        }
        fprintf(fp1, "\n");
    }
    xvgrclose(fp0);
    xvgrclose(fp1);
}

typedef struct {
    gmx_int64_t     np;
    double          sum;
    double          sav;
    double          sav2;
} ee_sum_t;

typedef struct {
    int             b;
    ee_sum_t        sum;
    gmx_int64_t     nst;
    gmx_int64_t     nst_min;
} ener_ee_t;

static void clear_ee_sum(ee_sum_t *ees)
{
    ees->sav  = 0;
    ees->sav2 = 0;
    ees->np   = 0;
    ees->sum  = 0;
}

static void add_ee_sum(ee_sum_t *ees, double sum, int np)
{
    ees->np  += np;
    ees->sum += sum;
}

static void add_ee_av(ee_sum_t *ees)
{
    double av;

    av         = ees->sum/ees->np;
    ees->sav  += av;
    ees->sav2 += av*av;
    ees->np    = 0;
    ees->sum   = 0;
}

static double calc_ee2(int nb, ee_sum_t *ees)
{
    return (ees->sav2/nb - dsqr(ees->sav/nb))/(nb - 1);
}

static void set_ee_av(ener_ee_t *eee)
{
    if (debug)
    {
        char buf[STEPSTRSIZE];
        fprintf(debug, "Storing average for err.est.: %s steps\n",
                gmx_step_str(eee->nst, buf));
    }
    add_ee_av(&eee->sum);
    eee->b++;
    if (eee->b == 1 || eee->nst < eee->nst_min)
    {
        eee->nst_min = eee->nst;
    }
    eee->nst = 0;
}

static void calc_averages(int nset, enerdata_t *edat, int nbmin, int nbmax)
{
    int             nb, i, f, nee;
    double          sum, sum2, sump, see2;
    gmx_int64_t     steps, np, p, bound_nb;
    enerdat_t      *ed;
    exactsum_t     *es;
    gmx_bool        bAllZero;
    double          x, sx, sy, sxx, sxy;
    ener_ee_t      *eee;

    /* Check if we have exact statistics over all points */
    for (i = 0; i < nset; i++)
    {
        ed             = &edat->s[i];
        ed->bExactStat = FALSE;
        if (edat->bHaveSums)
        {
            /* All energy file sum entries 0 signals no exact sums.
             * But if all energy values are 0, we still have exact sums.
             */
            bAllZero = TRUE;
            for (f = 0; f < edat->nframes && !ed->bExactStat; f++)
            {
                if (ed->ener[i] != 0)
                {
                    bAllZero = FALSE;
                }
                ed->bExactStat = (ed->es[f].sum != 0);
            }
            if (bAllZero)
            {
                ed->bExactStat = TRUE;
            }
        }
    }

    snew(eee, nbmax+1);
    for (i = 0; i < nset; i++)
    {
        ed = &edat->s[i];

        sum  = 0;
        sum2 = 0;
        np   = 0;
        sx   = 0;
        sy   = 0;
        sxx  = 0;
        sxy  = 0;
        for (nb = nbmin; nb <= nbmax; nb++)
        {
            eee[nb].b     = 0;
            clear_ee_sum(&eee[nb].sum);
            eee[nb].nst     = 0;
            eee[nb].nst_min = 0;
        }
        for (f = 0; f < edat->nframes; f++)
        {
            es = &ed->es[f];

            if (ed->bExactStat)
            {
                /* Add the sum and the sum of variances to the totals. */
                p     = edat->points[f];
                sump  = es->sum;
                sum2 += es->sum2;
                if (np > 0)
                {
                    sum2 += dsqr(sum/np - (sum + es->sum)/(np + p))
                        *np*(np + p)/p;
                }
            }
            else
            {
                /* Add a single value to the sum and sum of squares. */
                p     = 1;
                sump  = ed->ener[f];
                sum2 += dsqr(sump);
            }

            /* sum has to be increased after sum2 */
            np  += p;
            sum += sump;

            /* For the linear regression use variance 1/p.
             * Note that sump is the sum, not the average, so we don't need p*.
             */
            x    = edat->step[f] - 0.5*(edat->steps[f] - 1);
            sx  += p*x;
            sy  += sump;
            sxx += p*x*x;
            sxy += x*sump;

            for (nb = nbmin; nb <= nbmax; nb++)
            {
                /* Check if the current end step is closer to the desired
                 * block boundary than the next end step.
                 */
                bound_nb = (edat->step[0]-1)*nb + edat->nsteps*(eee[nb].b+1);
                if (eee[nb].nst > 0 &&
                    bound_nb - edat->step[f-1]*nb < edat->step[f]*nb - bound_nb)
                {
                    set_ee_av(&eee[nb]);
                }
                if (f == 0)
                {
                    eee[nb].nst = 1;
                }
                else
                {
                    eee[nb].nst += edat->step[f] - edat->step[f-1];
                }
                if (ed->bExactStat)
                {
                    add_ee_sum(&eee[nb].sum, es->sum, edat->points[f]);
                }
                else
                {
                    add_ee_sum(&eee[nb].sum, edat->s[i].ener[f], 1);
                }
                bound_nb = (edat->step[0]-1)*nb + edat->nsteps*(eee[nb].b+1);
                if (edat->step[f]*nb >= bound_nb)
                {
                    set_ee_av(&eee[nb]);
                }
            }
        }

        edat->s[i].av = sum/np;
        if (ed->bExactStat)
        {
            edat->s[i].rmsd = sqrt(sum2/np);
        }
        else
        {
            edat->s[i].rmsd = sqrt(sum2/np - dsqr(edat->s[i].av));
        }

        if (edat->nframes > 1)
        {
            edat->s[i].slope = (np*sxy - sx*sy)/(np*sxx - sx*sx);
        }
        else
        {
            edat->s[i].slope = 0;
        }

        nee  = 0;
        see2 = 0;
        for (nb = nbmin; nb <= nbmax; nb++)
        {
            /* Check if we actually got nb blocks and if the smallest
             * block is not shorter than 80% of the average.
             */
            if (debug)
            {
                char buf1[STEPSTRSIZE], buf2[STEPSTRSIZE];
                fprintf(debug, "Requested %d blocks, we have %d blocks, min %s nsteps %s\n",
                        nb, eee[nb].b,
                        gmx_step_str(eee[nb].nst_min, buf1),
                        gmx_step_str(edat->nsteps, buf2));
            }
            if (eee[nb].b == nb && 5*nb*eee[nb].nst_min >= 4*edat->nsteps)
            {
                see2 += calc_ee2(nb, &eee[nb].sum);
                nee++;
            }
        }
        if (nee > 0)
        {
            edat->s[i].ee = sqrt(see2/nee);
        }
        else
        {
            edat->s[i].ee = -1;
        }
    }
    sfree(eee);
}

static enerdata_t *calc_sum(int nset, enerdata_t *edat, int nbmin, int nbmax)
{
    enerdata_t *esum;
    enerdat_t  *s;
    int         f, i;
    double      sum;

    snew(esum, 1);
    *esum = *edat;
    snew(esum->s, 1);
    s = &esum->s[0];
    snew(s->ener, esum->nframes);
    snew(s->es, esum->nframes);

    s->bExactStat = TRUE;
    s->slope      = 0;
    for (i = 0; i < nset; i++)
    {
        if (!edat->s[i].bExactStat)
        {
            s->bExactStat = FALSE;
        }
        s->slope += edat->s[i].slope;
    }

    for (f = 0; f < edat->nframes; f++)
    {
        sum = 0;
        for (i = 0; i < nset; i++)
        {
            sum += edat->s[i].ener[f];
        }
        s->ener[f] = sum;
        sum        = 0;
        for (i = 0; i < nset; i++)
        {
            sum += edat->s[i].es[f].sum;
        }
        s->es[f].sum  = sum;
        s->es[f].sum2 = 0;
    }

    calc_averages(1, esum, nbmin, nbmax);

    return esum;
}

static char *ee_pr(double ee, char *buf)
{
    char   tmp[100];
    double rnd;

    if (ee < 0)
    {
        sprintf(buf, "%s", "--");
    }
    else
    {
        /* Round to two decimals by printing. */
        sprintf(tmp, "%.1e", ee);
        sscanf(tmp, "%lf", &rnd);
        sprintf(buf, "%g", rnd);
    }

    return buf;
}

static void remove_drift(int nset, int nbmin, int nbmax, real dt, enerdata_t *edat)
{
/* Remove the drift by performing a fit to y = ax+b.
   Uses 5 iterations. */
    int    i, j, k;
    double delta, d, sd, sd2;

    edat->npoints = edat->nframes;
    edat->nsteps  = edat->nframes;

    for (k = 0; (k < 5); k++)
    {
        for (i = 0; (i < nset); i++)
        {
            delta = edat->s[i].slope*dt;

            if (NULL != debug)
            {
                fprintf(debug, "slope for set %d is %g\n", i, edat->s[i].slope);
            }

            for (j = 0; (j < edat->nframes); j++)
            {
                edat->s[i].ener[j]   -= j*delta;
                edat->s[i].es[j].sum  = 0;
                edat->s[i].es[j].sum2 = 0;
            }
        }
        calc_averages(nset, edat, nbmin, nbmax);
    }
}

static void calc_fluctuation_props(FILE *fp,
                                   gmx_bool bDriftCorr, real dt,
                                   int nset, int nmol,
                                   char **leg, enerdata_t *edat,
                                   int nbmin, int nbmax)
{
    int    i, j;
    double vv, v, h, varv, hh, varh, tt, cv, cp, alpha, kappa, dcp, et, varet;
    double NANO3;
    enum {
        eVol, eEnth, eTemp, eEtot, eNR
    };
    const char *my_ener[] = { "Volume", "Enthalpy", "Temperature", "Total Energy" };
    int         ii[eNR];

    NANO3 = NANO*NANO*NANO;
    if (!bDriftCorr)
    {
        fprintf(fp, "\nYou may want to use the -driftcorr flag in order to correct\n"
                "for spurious drift in the graphs. Note that this is not\n"
                "a substitute for proper equilibration and sampling!\n");
    }
    else
    {
        remove_drift(nset, nbmin, nbmax, dt, edat);
    }
    for (i = 0; (i < eNR); i++)
    {
        for (ii[i] = 0; (ii[i] < nset &&
                         (gmx_strcasecmp(leg[ii[i]], my_ener[i]) != 0)); ii[i]++)
        {
            ;
        }
/*        if (ii[i] < nset)
            fprintf(fp,"Found %s data.\n",my_ener[i]);
 */ }
    /* Compute it all! */
    alpha = kappa = cp = dcp = cv = NOTSET;

    /* Temperature */
    tt = NOTSET;
    if (ii[eTemp] < nset)
    {
        tt    = edat->s[ii[eTemp]].av;
    }
    /* Volume */
    vv = varv = NOTSET;
    if ((ii[eVol] < nset) && (ii[eTemp] < nset))
    {
        vv    = edat->s[ii[eVol]].av*NANO3;
        varv  = dsqr(edat->s[ii[eVol]].rmsd*NANO3);
        kappa = (varv/vv)/(BOLTZMANN*tt);
    }
    /* Enthalpy */
    hh = varh = NOTSET;
    if ((ii[eEnth] < nset) && (ii[eTemp] < nset))
    {
        hh    = KILO*edat->s[ii[eEnth]].av/AVOGADRO;
        varh  = dsqr(KILO*edat->s[ii[eEnth]].rmsd/AVOGADRO);
        cp    = AVOGADRO*((varh/nmol)/(BOLTZMANN*tt*tt));
    }
    /* Total energy */
    et = varet = NOTSET;
    if ((ii[eEtot] < nset) && (hh == NOTSET) && (tt != NOTSET))
    {
        /* Only compute cv in constant volume runs, which we can test
           by checking whether the enthalpy was computed.
         */
        et    = edat->s[ii[eEtot]].av;
        varet = sqr(edat->s[ii[eEtot]].rmsd);
        cv    = KILO*((varet/nmol)/(BOLTZ*tt*tt));
    }
    /* Alpha, dcp */
    if ((ii[eVol] < nset) && (ii[eEnth] < nset) && (ii[eTemp] < nset))
    {
        double v_sum, h_sum, vh_sum, v_aver, h_aver, vh_aver;
        vh_sum = v_sum = h_sum = 0;
        for (j = 0; (j < edat->nframes); j++)
        {
            v       = edat->s[ii[eVol]].ener[j]*NANO3;
            h       = KILO*edat->s[ii[eEnth]].ener[j]/AVOGADRO;
            v_sum  += v;
            h_sum  += h;
            vh_sum += (v*h);
        }
        vh_aver = vh_sum / edat->nframes;
        v_aver  = v_sum  / edat->nframes;
        h_aver  = h_sum  / edat->nframes;
        alpha   = (vh_aver-v_aver*h_aver)/(v_aver*BOLTZMANN*tt*tt);
        dcp     = (v_aver*AVOGADRO/nmol)*tt*sqr(alpha)/(kappa);
    }

    if (tt != NOTSET)
    {
        if (nmol < 2)
        {
            fprintf(fp, "\nWARNING: nmol = %d, this may not be what you want.\n",
                    nmol);
        }
        fprintf(fp, "\nTemperature dependent fluctuation properties at T = %g.\n", tt);
        fprintf(fp, "\nHeat capacities obtained from fluctuations do *not* include\n");
        fprintf(fp, "quantum corrections. If you want to get a more accurate estimate\n");
        fprintf(fp, "please use the g_dos program.\n\n");
        fprintf(fp, "WARNING: Please verify that your simulations are converged and perform\n"
                "a block-averaging error analysis (not implemented in g_energy yet)\n");

        if (debug != NULL)
        {
            if (varv != NOTSET)
            {
                fprintf(fp, "varv  =  %10g (m^6)\n", varv*AVOGADRO/nmol);
            }
        }
        if (vv != NOTSET)
        {
            fprintf(fp, "Volume                                   = %10g m^3/mol\n",
                    vv*AVOGADRO/nmol);
        }
        if (varh != NOTSET)
        {
            fprintf(fp, "Enthalpy                                 = %10g kJ/mol\n",
                    hh*AVOGADRO/(KILO*nmol));
        }
        if (alpha != NOTSET)
        {
            fprintf(fp, "Coefficient of Thermal Expansion Alpha_P = %10g (1/K)\n",
                    alpha);
        }
        if (kappa != NOTSET)
        {
            fprintf(fp, "Isothermal Compressibility Kappa         = %10g (J/m^3)\n",
                    kappa);
            fprintf(fp, "Adiabatic bulk modulus                   = %10g (m^3/J)\n",
                    1.0/kappa);
        }
        if (cp != NOTSET)
        {
            fprintf(fp, "Heat capacity at constant pressure Cp    = %10g J/mol K\n",
                    cp);
        }
        if (cv != NOTSET)
        {
            fprintf(fp, "Heat capacity at constant volume Cv      = %10g J/mol K\n",
                    cv);
        }
        if (dcp != NOTSET)
        {
            fprintf(fp, "Cp-Cv                                    =  %10g J/mol K\n",
                    dcp);
        }
        please_cite(fp, "Allen1987a");
    }
    else
    {
        fprintf(fp, "You should select the temperature in order to obtain fluctuation properties.\n");
    }
}

static void analyse_ener(gmx_bool bCorr, const char *corrfn,
                         gmx_bool bFee, gmx_bool bSum, gmx_bool bFluct,
                         gmx_bool bVisco, const char *visfn, int nmol,
                         gmx_int64_t start_step, double start_t,
                         gmx_int64_t step, double t,
                         real reftemp,
                         enerdata_t *edat,
                         int nset, int set[], gmx_bool *bIsEner,
                         char **leg, gmx_enxnm_t *enm,
                         real Vaver, real ezero,
                         int nbmin, int nbmax,
                         const output_env_t oenv)
{
    FILE           *fp;
    /* Check out the printed manual for equations! */
    double          Dt, aver, stddev, errest, delta_t, totaldrift;
    enerdata_t     *esum = NULL;
    real            xxx, integral, intBulk, Temp = 0, Pres = 0;
    real            sfrac, oldfrac, diffsum, diffav, fstep, pr_aver, pr_stddev, pr_errest;
    double          beta = 0, expE, expEtot, *fee = NULL;
    gmx_int64_t     nsteps;
    int             nexact, nnotexact;
    double          x1m, x1mk;
    int             i, j, k, nout;
    real            chi2;
    char            buf[256], eebuf[100];

    nsteps  = step - start_step + 1;
    if (nsteps < 1)
    {
        fprintf(stdout, "Not enough steps (%s) for statistics\n",
                gmx_step_str(nsteps, buf));
    }
    else
    {
        /* Calculate the time difference */
        delta_t = t - start_t;

        fprintf(stdout, "\nStatistics over %s steps [ %.4f through %.4f ps ], %d data sets\n",
                gmx_step_str(nsteps, buf), start_t, t, nset);

        calc_averages(nset, edat, nbmin, nbmax);

        if (bSum)
        {
            esum = calc_sum(nset, edat, nbmin, nbmax);
        }

        if (!edat->bHaveSums)
        {
            nexact    = 0;
            nnotexact = nset;
        }
        else
        {
            nexact    = 0;
            nnotexact = 0;
            for (i = 0; (i < nset); i++)
            {
                if (edat->s[i].bExactStat)
                {
                    nexact++;
                }
                else
                {
                    nnotexact++;
                }
            }
        }

        if (nnotexact == 0)
        {
            fprintf(stdout, "All statistics are over %s points\n",
                    gmx_step_str(edat->npoints, buf));
        }
        else if (nexact == 0 || edat->npoints == edat->nframes)
        {
            fprintf(stdout, "All statistics are over %d points (frames)\n",
                    edat->nframes);
        }
        else
        {
            fprintf(stdout, "The term%s", nnotexact == 1 ? "" : "s");
            for (i = 0; (i < nset); i++)
            {
                if (!edat->s[i].bExactStat)
                {
                    fprintf(stdout, " '%s'", leg[i]);
                }
            }
            fprintf(stdout, " %s has statistics over %d points (frames)\n",
                    nnotexact == 1 ? "is" : "are", edat->nframes);
            fprintf(stdout, "All other statistics are over %s points\n",
                    gmx_step_str(edat->npoints, buf));
        }
        fprintf(stdout, "\n");

        fprintf(stdout, "%-24s %10s %10s %10s %10s",
                "Energy", "Average", "Err.Est.", "RMSD", "Tot-Drift");
        if (bFee)
        {
            fprintf(stdout, "  %10s\n", "-kT ln<e^(E/kT)>");
        }
        else
        {
            fprintf(stdout, "\n");
        }
        fprintf(stdout, "-------------------------------------------------------------------------------\n");

        /* Initiate locals, only used with -sum */
        expEtot = 0;
        if (bFee)
        {
            beta = 1.0/(BOLTZ*reftemp);
            snew(fee, nset);
        }
        for (i = 0; (i < nset); i++)
        {
            aver   = edat->s[i].av;
            stddev = edat->s[i].rmsd;
            errest = edat->s[i].ee;

            if (bFee)
            {
                expE = 0;
                for (j = 0; (j < edat->nframes); j++)
                {
                    expE += exp(beta*(edat->s[i].ener[j] - aver)/nmol);
                }
                if (bSum)
                {
                    expEtot += expE/edat->nframes;
                }

                fee[i] = log(expE/edat->nframes)/beta + aver/nmol;
            }
            if (strstr(leg[i], "empera") != NULL)
            {
                Temp = aver;
            }
            else if (strstr(leg[i], "olum") != NULL)
            {
                Vaver = aver;
            }
            else if (strstr(leg[i], "essure") != NULL)
            {
                Pres = aver;
            }
            if (bIsEner[i])
            {
                pr_aver   = aver/nmol-ezero;
                pr_stddev = stddev/nmol;
                pr_errest = errest/nmol;
            }
            else
            {
                pr_aver   = aver;
                pr_stddev = stddev;
                pr_errest = errest;
            }

            /* Multiply the slope in steps with the number of steps taken */
            totaldrift = (edat->nsteps - 1)*edat->s[i].slope;
            if (bIsEner[i])
            {
                totaldrift /= nmol;
            }

            fprintf(stdout, "%-24s %10g %10s %10g %10g",
                    leg[i], pr_aver, ee_pr(pr_errest, eebuf), pr_stddev, totaldrift);
            if (bFee)
            {
                fprintf(stdout, "  %10g", fee[i]);
            }

            fprintf(stdout, "  (%s)\n", enm[set[i]].unit);

            if (bFluct)
            {
                for (j = 0; (j < edat->nframes); j++)
                {
                    edat->s[i].ener[j] -= aver;
                }
            }
        }
        if (bSum)
        {
            totaldrift = (edat->nsteps - 1)*esum->s[0].slope;
            fprintf(stdout, "%-24s %10g %10s %10s %10g  (%s)",
                    "Total", esum->s[0].av/nmol, ee_pr(esum->s[0].ee/nmol, eebuf),
                    "--", totaldrift/nmol, enm[set[0]].unit);
            /* pr_aver,pr_stddev,a,totaldrift */
            if (bFee)
            {
                fprintf(stdout, "  %10g  %10g\n",
                        log(expEtot)/beta + esum->s[0].av/nmol, log(expEtot)/beta);
            }
            else
            {
                fprintf(stdout, "\n");
            }
        }

        /* Do correlation function */
        if (edat->nframes > 1)
        {
            Dt = delta_t/(edat->nframes - 1);
        }
        else
        {
            Dt = 0;
        }
        if (bVisco)
        {
            const char* leg[] = { "Shear", "Bulk" };
            real        factor;
            real      **eneset;
            real      **eneint;

            /* Assume pressure tensor is in Pxx Pxy Pxz Pyx Pyy Pyz Pzx Pzy Pzz */

            /* Symmetrise tensor! (and store in first three elements)
             * And subtract average pressure!
             */
            snew(eneset, 12);
            for (i = 0; i < 12; i++)
            {
                snew(eneset[i], edat->nframes);
            }
            for (i = 0; (i < edat->nframes); i++)
            {
                eneset[0][i] = 0.5*(edat->s[1].ener[i]+edat->s[3].ener[i]);
                eneset[1][i] = 0.5*(edat->s[2].ener[i]+edat->s[6].ener[i]);
                eneset[2][i] = 0.5*(edat->s[5].ener[i]+edat->s[7].ener[i]);
                for (j = 3; j <= 11; j++)
                {
                    eneset[j][i] = edat->s[j].ener[i];
                }
                eneset[11][i] -= Pres;
            }

            /* Determine integrals of the off-diagonal pressure elements */
            snew(eneint, 3);
            for (i = 0; i < 3; i++)
            {
                snew(eneint[i], edat->nframes + 1);
            }
            eneint[0][0] = 0;
            eneint[1][0] = 0;
            eneint[2][0] = 0;
            for (i = 0; i < edat->nframes; i++)
            {
                eneint[0][i+1] = eneint[0][i] + 0.5*(edat->s[1].es[i].sum + edat->s[3].es[i].sum)*Dt/edat->points[i];
                eneint[1][i+1] = eneint[1][i] + 0.5*(edat->s[2].es[i].sum + edat->s[6].es[i].sum)*Dt/edat->points[i];
                eneint[2][i+1] = eneint[2][i] + 0.5*(edat->s[5].es[i].sum + edat->s[7].es[i].sum)*Dt/edat->points[i];
            }

            einstein_visco("evisco.xvg", "eviscoi.xvg",
                           3, edat->nframes+1, eneint, Vaver, Temp, Dt, oenv);

            for (i = 0; i < 3; i++)
            {
                sfree(eneint[i]);
            }
            sfree(eneint);

            /*do_autocorr(corrfn,buf,nenergy,3,eneset,Dt,eacNormal,TRUE);*/
            /* Do it for shear viscosity */
            strcpy(buf, "Shear Viscosity");
            low_do_autocorr(corrfn, oenv, buf, edat->nframes, 3,
                            (edat->nframes+1)/2, eneset, Dt,
                            eacNormal, 1, TRUE, FALSE, FALSE, 0.0, 0.0, 0);

            /* Now for bulk viscosity */
            strcpy(buf, "Bulk Viscosity");
            low_do_autocorr(corrfn, oenv, buf, edat->nframes, 1,
                            (edat->nframes+1)/2, &(eneset[11]), Dt,
                            eacNormal, 1, TRUE, FALSE, FALSE, 0.0, 0.0, 0);

            factor = (Vaver*1e-26/(BOLTZMANN*Temp))*Dt;
            fp     = xvgropen(visfn, buf, "Time (ps)", "\\8h\\4 (cp)", oenv);
            xvgr_legend(fp, asize(leg), leg, oenv);

            /* Use trapezium rule for integration */
            integral = 0;
            intBulk  = 0;
            nout     = get_acfnout();
            if ((nout < 2) || (nout >= edat->nframes/2))
            {
                nout = edat->nframes/2;
            }
            for (i = 1; (i < nout); i++)
            {
                integral += 0.5*(eneset[0][i-1]  + eneset[0][i])*factor;
                intBulk  += 0.5*(eneset[11][i-1] + eneset[11][i])*factor;
                fprintf(fp, "%10g  %10g  %10g\n", (i*Dt), integral, intBulk);
            }
            xvgrclose(fp);

            for (i = 0; i < 12; i++)
            {
                sfree(eneset[i]);
            }
            sfree(eneset);
        }
        else if (bCorr)
        {
            if (bFluct)
            {
                strcpy(buf, "Autocorrelation of Energy Fluctuations");
            }
            else
            {
                strcpy(buf, "Energy Autocorrelation");
            }
#if 0
            do_autocorr(corrfn, oenv, buf, edat->nframes,
                        bSum ? 1                 : nset,
                        bSum ? &edat->s[nset-1].ener : eneset,
                        (delta_t/edat->nframes), eacNormal, FALSE);
#endif
        }
    }
}

static void print_time(FILE *fp, double t)
{
    fprintf(fp, "%12.6f", t);
}

static void print1(FILE *fp, gmx_bool bDp, real e)
{
    if (bDp)
    {
        fprintf(fp, "  %16.12f", e);
    }
    else
    {
        fprintf(fp, "  %10.6f", e);
    }
}

static void fec(const char *ene2fn, const char *runavgfn,
                real reftemp, int nset, int set[], char *leg[],
                enerdata_t *edat, double time[],
                const output_env_t oenv)
{
    const char * ravgleg[] = {
        "\\8D\\4E = E\\sB\\N-E\\sA\\N",
        "<e\\S-\\8D\\4E/kT\\N>\\s0..t\\N"
    };
    FILE        *fp;
    ener_file_t  enx;
    int          nre, timecheck, step, nenergy, nenergy2, maxenergy;
    int          i, j;
    gmx_bool     bCont;
    real         aver, beta;
    real       **eneset2;
    double       dE, sum;
    gmx_enxnm_t *enm = NULL;
    t_enxframe  *fr;
    char         buf[22];

    /* read second energy file */
    snew(fr, 1);
    enm = NULL;
    enx = open_enx(ene2fn, "r");
    do_enxnms(enx, &(fr->nre), &enm);

    snew(eneset2, nset+1);
    nenergy2  = 0;
    maxenergy = 0;
    timecheck = 0;
    do
    {
        /* This loop searches for the first frame (when -b option is given),
         * or when this has been found it reads just one energy frame
         */
        do
        {
            bCont = do_enx(enx, fr);

            if (bCont)
            {
                timecheck = check_times(fr->t);
            }

        }
        while (bCont && (timecheck < 0));

        /* Store energies for analysis afterwards... */
        if ((timecheck == 0) && bCont)
        {
            if (fr->nre > 0)
            {
                if (nenergy2 >= maxenergy)
                {
                    maxenergy += 1000;
                    for (i = 0; i <= nset; i++)
                    {
                        srenew(eneset2[i], maxenergy);
                    }
                }
                if (fr->t != time[nenergy2])
                {
                    fprintf(stderr, "\nWARNING time mismatch %g!=%g at frame %s\n",
                            fr->t, time[nenergy2], gmx_step_str(fr->step, buf));
                }
                for (i = 0; i < nset; i++)
                {
                    eneset2[i][nenergy2] = fr->ener[set[i]].e;
                }
                nenergy2++;
            }
        }
    }
    while (bCont && (timecheck == 0));

    /* check */
    if (edat->nframes != nenergy2)
    {
        fprintf(stderr, "\nWARNING file length mismatch %d!=%d\n",
                edat->nframes, nenergy2);
    }
    nenergy = min(edat->nframes, nenergy2);

    /* calculate fe difference dF = -kT ln < exp(-(E_B-E_A)/kT) >_A */
    fp = NULL;
    if (runavgfn)
    {
        fp = xvgropen(runavgfn, "Running average free energy difference",
                      "Time (" unit_time ")", "\\8D\\4E (" unit_energy ")", oenv);
        xvgr_legend(fp, asize(ravgleg), ravgleg, oenv);
    }
    fprintf(stdout, "\n%-24s %10s\n",
            "Energy", "dF = -kT ln < exp(-(EB-EA)/kT) >A");
    sum  = 0;
    beta = 1.0/(BOLTZ*reftemp);
    for (i = 0; i < nset; i++)
    {
        if (gmx_strcasecmp(leg[i], enm[set[i]].name) != 0)
        {
            fprintf(stderr, "\nWARNING energy set name mismatch %s!=%s\n",
                    leg[i], enm[set[i]].name);
        }
        for (j = 0; j < nenergy; j++)
        {
            dE   = eneset2[i][j] - edat->s[i].ener[j];
            sum += exp(-dE*beta);
            if (fp)
            {
                fprintf(fp, "%10g %10g %10g\n",
                        time[j], dE, -BOLTZ*reftemp*log(sum/(j+1)) );
            }
        }
        aver = -BOLTZ*reftemp*log(sum/nenergy);
        fprintf(stdout, "%-24s %10g\n", leg[i], aver);
    }
    if (fp)
    {
        xvgrclose(fp);
    }
    sfree(fr);
}


static void do_dhdl(t_enxframe *fr, t_inputrec *ir, FILE **fp_dhdl,
                    const char *filename, gmx_bool bDp,
                    int *blocks, int *hists, int *samples, int *nlambdas,
                    const output_env_t oenv)
{
    const char  *dhdl = "dH/d\\lambda", *deltag = "\\DeltaH", *lambda = "\\lambda";
    char         title[STRLEN], label_x[STRLEN], label_y[STRLEN], legend[STRLEN];
    char         buf[STRLEN];
    gmx_bool     first       = FALSE;
    int          nblock_hist = 0, nblock_dh = 0, nblock_dhcoll = 0;
    int          i, j, k;
    /* coll data */
    double       temp              = 0, start_time = 0, delta_time = 0, start_lambda = 0, delta_lambda = 0;
    static int   setnr             = 0;
    double      *native_lambda_vec = NULL;
    const char **lambda_components = NULL;
    int          n_lambda_vec      = 0;
    gmx_bool     changing_lambda   = FALSE;
    int          lambda_fep_state;

    /* now count the blocks & handle the global dh data */
    for (i = 0; i < fr->nblock; i++)
    {
        if (fr->block[i].id == enxDHHIST)
        {
            nblock_hist++;
        }
        else if (fr->block[i].id == enxDH)
        {
            nblock_dh++;
        }
        else if (fr->block[i].id == enxDHCOLL)
        {
            nblock_dhcoll++;
            if ( (fr->block[i].nsub < 1) ||
                 (fr->block[i].sub[0].type != xdr_datatype_double) ||
                 (fr->block[i].sub[0].nr < 5))
            {
                gmx_fatal(FARGS, "Unexpected block data");
            }

            /* read the data from the DHCOLL block */
            temp            =         fr->block[i].sub[0].dval[0];
            start_time      =   fr->block[i].sub[0].dval[1];
            delta_time      =   fr->block[i].sub[0].dval[2];
            start_lambda    = fr->block[i].sub[0].dval[3];
            delta_lambda    = fr->block[i].sub[0].dval[4];
            changing_lambda = (delta_lambda != 0);
            if (fr->block[i].nsub > 1)
            {
                lambda_fep_state = fr->block[i].sub[1].ival[0];
                if (n_lambda_vec == 0)
                {
                    n_lambda_vec = fr->block[i].sub[1].ival[1];
                }
                else
                {
                    if (n_lambda_vec != fr->block[i].sub[1].ival[1])
                    {
                        gmx_fatal(FARGS,
                                  "Unexpected change of basis set in lambda");
                    }
                }
                if (lambda_components == NULL)
                {
                    snew(lambda_components, n_lambda_vec);
                }
                if (native_lambda_vec == NULL)
                {
                    snew(native_lambda_vec, n_lambda_vec);
                }
                for (j = 0; j < n_lambda_vec; j++)
                {
                    native_lambda_vec[j] = fr->block[i].sub[0].dval[5+j];
                    lambda_components[j] =
                        efpt_singular_names[fr->block[i].sub[1].ival[2+j]];
                }
            }
        }
    }

    if (nblock_hist == 0 && nblock_dh == 0)
    {
        /* don't do anything */
        return;
    }
    if (nblock_hist > 0 && nblock_dh > 0)
    {
        gmx_fatal(FARGS, "This energy file contains both histogram dhdl data and non-histogram dhdl data. Don't know what to do.");
    }
    if (!*fp_dhdl)
    {
        if (nblock_dh > 0)
        {
            /* we have standard, non-histogram data --
               call open_dhdl to open the file */
            /* TODO this is an ugly hack that needs to be fixed: this will only
               work if the order of data is always the same and if we're
               only using the g_energy compiled with the mdrun that produced
               the ener.edr. */
            *fp_dhdl = open_dhdl(filename, ir, oenv);
        }
        else
        {
            sprintf(title, "N(%s)", deltag);
            sprintf(label_x, "%s (%s)", deltag, unit_energy);
            sprintf(label_y, "Samples");
            *fp_dhdl = xvgropen_type(filename, title, label_x, label_y, exvggtXNY, oenv);
            sprintf(buf, "T = %g (K), %s = %g", temp, lambda, start_lambda);
            xvgr_subtitle(*fp_dhdl, buf, oenv);
        }
    }

    (*hists)   += nblock_hist;
    (*blocks)  += nblock_dh;
    (*nlambdas) = nblock_hist+nblock_dh;

    /* write the data */
    if (nblock_hist > 0)
    {
        gmx_int64_t sum = 0;
        /* histograms */
        for (i = 0; i < fr->nblock; i++)
        {
            t_enxblock *blk = &(fr->block[i]);
            if (blk->id == enxDHHIST)
            {
                double          foreign_lambda, dx;
                gmx_int64_t     x0;
                int             nhist, derivative;

                /* check the block types etc. */
                if ( (blk->nsub < 2) ||
                     (blk->sub[0].type != xdr_datatype_double) ||
                     (blk->sub[1].type != xdr_datatype_int64) ||
                     (blk->sub[0].nr < 2)  ||
                     (blk->sub[1].nr < 2) )
                {
                    gmx_fatal(FARGS, "Unexpected block data in file");
                }
                foreign_lambda = blk->sub[0].dval[0];
                dx             = blk->sub[0].dval[1];
                nhist          = blk->sub[1].lval[0];
                derivative     = blk->sub[1].lval[1];
                for (j = 0; j < nhist; j++)
                {
                    const char *lg[1];
                    x0 = blk->sub[1].lval[2+j];

                    if (!derivative)
                    {
                        sprintf(legend, "N(%s(%s=%g) | %s=%g)",
                                deltag, lambda, foreign_lambda,
                                lambda, start_lambda);
                    }
                    else
                    {
                        sprintf(legend, "N(%s | %s=%g)",
                                dhdl, lambda, start_lambda);
                    }

                    lg[0] = legend;
                    xvgr_new_dataset(*fp_dhdl, setnr, 1, lg, oenv);
                    setnr++;
                    for (k = 0; k < blk->sub[j+2].nr; k++)
                    {
                        int    hist;
                        double xmin, xmax;

                        hist = blk->sub[j+2].ival[k];
                        xmin = (x0+k)*dx;
                        xmax = (x0+k+1)*dx;
                        fprintf(*fp_dhdl, "%g %d\n%g %d\n", xmin, hist,
                                xmax, hist);
                        sum += hist;
                    }
                    /* multiple histogram data blocks in one histogram
                       mean that the second one is the reverse of the first one:
                       for dhdl derivatives, it's important to know both the
                       maximum and minimum values */
                    dx = -dx;
                }
            }
        }
        (*samples) += (int)(sum/nblock_hist);
    }
    else
    {
        /* raw dh */
        int    len      = 0;
        char **setnames = NULL;
        int    nnames   = nblock_dh;

        for (i = 0; i < fr->nblock; i++)
        {
            t_enxblock *blk = &(fr->block[i]);
            if (blk->id == enxDH)
            {
                if (len == 0)
                {
                    len = blk->sub[2].nr;
                }
                else
                {
                    if (len != blk->sub[2].nr)
                    {
                        gmx_fatal(FARGS, "Length inconsistency in dhdl data");
                    }
                }
            }
        }
        (*samples) += len;

        for (i = 0; i < len; i++)
        {
            double time = start_time + delta_time*i;

            fprintf(*fp_dhdl, "%.4f ", time);

            for (j = 0; j < fr->nblock; j++)
            {
                t_enxblock *blk = &(fr->block[j]);
                if (blk->id == enxDH)
                {
                    double value;
                    if (blk->sub[2].type == xdr_datatype_float)
                    {
                        value = blk->sub[2].fval[i];
                    }
                    else
                    {
                        value = blk->sub[2].dval[i];
                    }
                    /* we need to decide which data type it is based on the count*/

                    if (j == 1 && ir->bExpanded)
                    {
                        fprintf(*fp_dhdl, "%4d", (int)value);   /* if expanded ensembles and zero, this is a state value, it's an integer. We need a cleaner conditional than if j==1! */
                    }
                    else
                    {
                        if (bDp)
                        {
                            fprintf(*fp_dhdl, " %#.12g", value);   /* print normal precision */
                        }
                        else
                        {
                            fprintf(*fp_dhdl, " %#.8g", value);   /* print normal precision */
                        }
                    }
                }
            }
            fprintf(*fp_dhdl, "\n");
        }
    }
}


int gmx_energy(int argc, char *argv[])
{
    const char        *desc[] = {
        "[THISMODULE] extracts energy components or distance restraint",
        "data from an energy file. The user is prompted to interactively",
        "select the desired energy terms.[PAR]",

        "Average, RMSD, and drift are calculated with full precision from the",
        "simulation (see printed manual). Drift is calculated by performing",
        "a least-squares fit of the data to a straight line. The reported total drift",
        "is the difference of the fit at the first and last point.",
        "An error estimate of the average is given based on a block averages",
        "over 5 blocks using the full-precision averages. The error estimate",
        "can be performed over multiple block lengths with the options",
        "[TT]-nbmin[tt] and [TT]-nbmax[tt].",
        "[BB]Note[bb] that in most cases the energy files contains averages over all",
        "MD steps, or over many more points than the number of frames in",
        "energy file. This makes the [THISMODULE] statistics output more accurate",
        "than the [REF].xvg[ref] output. When exact averages are not present in the energy",
        "file, the statistics mentioned above are simply over the single, per-frame",
        "energy values.[PAR]",

        "The term fluctuation gives the RMSD around the least-squares fit.[PAR]",

        "Some fluctuation-dependent properties can be calculated provided",
        "the correct energy terms are selected, and that the command line option",
        "[TT]-fluct_props[tt] is given. The following properties",
        "will be computed:",
        "",
        "===============================  ===================",
        "Property                         Energy terms needed",
        "===============================  ===================",
        "Heat capacity C[SUB]p[sub] (NPT sims):    Enthalpy, Temp",
        "Heat capacity C[SUB]v[sub] (NVT sims):    Etot, Temp",
        "Thermal expansion coeff. (NPT):  Enthalpy, Vol, Temp",
        "Isothermal compressibility:      Vol, Temp",
        "Adiabatic bulk modulus:          Vol, Temp",
        "===============================  ===================",
        "",
        "You always need to set the number of molecules [TT]-nmol[tt].",
        "The C[SUB]p[sub]/C[SUB]v[sub] computations do [BB]not[bb] include any corrections",
        "for quantum effects. Use the [gmx-dos] program if you need that (and you do).[PAR]"
        "When the [TT]-viol[tt] option is set, the time averaged",
        "violations are plotted and the running time-averaged and",
        "instantaneous sum of violations are recalculated. Additionally",
        "running time-averaged and instantaneous distances between",
        "selected pairs can be plotted with the [TT]-pairs[tt] option.[PAR]",

        "Options [TT]-ora[tt], [TT]-ort[tt], [TT]-oda[tt], [TT]-odr[tt] and",
        "[TT]-odt[tt] are used for analyzing orientation restraint data.",
        "The first two options plot the orientation, the last three the",
        "deviations of the orientations from the experimental values.",
        "The options that end on an 'a' plot the average over time",
        "as a function of restraint. The options that end on a 't'",
        "prompt the user for restraint label numbers and plot the data",
        "as a function of time. Option [TT]-odr[tt] plots the RMS",
        "deviation as a function of restraint.",
        "When the run used time or ensemble averaged orientation restraints,",
        "option [TT]-orinst[tt] can be used to analyse the instantaneous,",
        "not ensemble-averaged orientations and deviations instead of",
        "the time and ensemble averages.[PAR]",

        "Option [TT]-oten[tt] plots the eigenvalues of the molecular order",
        "tensor for each orientation restraint experiment. With option",
        "[TT]-ovec[tt] also the eigenvectors are plotted.[PAR]",

        "Option [TT]-odh[tt] extracts and plots the free energy data",
        "(Hamiltoian differences and/or the Hamiltonian derivative dhdl)",
        "from the [TT]ener.edr[tt] file.[PAR]",

        "With [TT]-fee[tt] an estimate is calculated for the free-energy",
        "difference with an ideal gas state::",
        "",
        "  [GRK]Delta[grk] A = A(N,V,T) - A[SUB]idealgas[sub](N,V,T) = kT [LN][CHEVRON][EXP]U[SUB]pot[sub]/kT[exp][chevron][ln]",
        "  [GRK]Delta[grk] G = G(N,p,T) - G[SUB]idealgas[sub](N,p,T) = kT [LN][CHEVRON][EXP]U[SUB]pot[sub]/kT[exp][chevron][ln]",
        "",
        "where k is Boltzmann's constant, T is set by [TT]-fetemp[tt] and",
        "the average is over the ensemble (or time in a trajectory).",
        "Note that this is in principle",
        "only correct when averaging over the whole (Boltzmann) ensemble",
        "and using the potential energy. This also allows for an entropy",
        "estimate using::",
        "",
        "  [GRK]Delta[grk] S(N,V,T) = S(N,V,T) - S[SUB]idealgas[sub](N,V,T) = ([CHEVRON]U[SUB]pot[sub][chevron] - [GRK]Delta[grk] A)/T",
        "  [GRK]Delta[grk] S(N,p,T) = S(N,p,T) - S[SUB]idealgas[sub](N,p,T) = ([CHEVRON]U[SUB]pot[sub][chevron] + pV - [GRK]Delta[grk] G)/T",
        "",

        "When a second energy file is specified ([TT]-f2[tt]), a free energy",
        "difference is calculated::",
        "",
        "  dF = -kT [LN][CHEVRON][EXP]-(E[SUB]B[sub]-E[SUB]A[sub])/kT[exp][chevron][SUB]A[sub][ln] ,",
        "",
        "where E[SUB]A[sub] and E[SUB]B[sub] are the energies from the first and second energy",
        "files, and the average is over the ensemble A. The running average",
        "of the free energy difference is printed to a file specified by [TT]-ravg[tt].",
        "[BB]Note[bb] that the energies must both be calculated from the same trajectory."

    };
    static gmx_bool    bSum    = FALSE, bFee = FALSE, bPrAll = FALSE, bFluct = FALSE, bDriftCorr = FALSE;
    static gmx_bool    bDp     = FALSE, bMutot = FALSE, bOrinst = FALSE, bOvec = FALSE, bFluctProps = FALSE;
    static int         skip    = 0, nmol = 1, nbmin = 5, nbmax = 5;
    static real        reftemp = 300.0, ezero = 0;
    t_pargs            pa[]    = {
        { "-fee",   FALSE, etBOOL,  {&bFee},
          "Do a free energy estimate" },
        { "-fetemp", FALSE, etREAL, {&reftemp},
          "Reference temperature for free energy calculation" },
        { "-zero", FALSE, etREAL, {&ezero},
          "Subtract a zero-point energy" },
        { "-sum",  FALSE, etBOOL, {&bSum},
          "Sum the energy terms selected rather than display them all" },
        { "-dp",   FALSE, etBOOL, {&bDp},
          "Print energies in high precision" },
        { "-nbmin", FALSE, etINT, {&nbmin},
          "Minimum number of blocks for error estimate" },
        { "-nbmax", FALSE, etINT, {&nbmax},
          "Maximum number of blocks for error estimate" },
        { "-mutot", FALSE, etBOOL, {&bMutot},
          "Compute the total dipole moment from the components" },
        { "-skip", FALSE, etINT,  {&skip},
          "Skip number of frames between data points" },
        { "-aver", FALSE, etBOOL, {&bPrAll},
          "Also print the exact average and rmsd stored in the energy frames (only when 1 term is requested)" },
        { "-nmol", FALSE, etINT,  {&nmol},
          "Number of molecules in your sample: the energies are divided by this number" },
        { "-fluct_props", FALSE, etBOOL, {&bFluctProps},
          "Compute properties based on energy fluctuations, like heat capacity" },
        { "-driftcorr", FALSE, etBOOL, {&bDriftCorr},
          "Useful only for calculations of fluctuation properties. The drift in the observables will be subtracted before computing the fluctuation properties."},
        { "-fluc", FALSE, etBOOL, {&bFluct},
          "Calculate autocorrelation of energy fluctuations rather than energy itself" },
        { "-orinst", FALSE, etBOOL, {&bOrinst},
          "Analyse instantaneous orientation data" },
        { "-ovec", FALSE, etBOOL, {&bOvec},
          "Also plot the eigenvectors with [TT]-oten[tt]" }
    };
    const char       * drleg[] = {
        "Running average",
        "Instantaneous"
    };
    static const char *setnm[] = {
        "Pres-XX", "Pres-XY", "Pres-XZ", "Pres-YX", "Pres-YY",
        "Pres-YZ", "Pres-ZX", "Pres-ZY", "Pres-ZZ", "Temperature",
        "Volume",  "Pressure"
    };

    FILE              *out     = NULL, *fp_pairs = NULL, *fort = NULL, *fodt = NULL, *foten = NULL;
    FILE              *fp_dhdl = NULL;
    FILE             **drout;
    ener_file_t        fp;
    int                timecheck = 0;
    gmx_mtop_t         mtop;
    gmx_localtop_t    *top = NULL;
    t_inputrec         ir;
    t_energy         **ee;
    enerdata_t         edat;
    gmx_enxnm_t       *enm = NULL;
    t_enxframe        *frame, *fr = NULL;
    int                cur = 0;
#define NEXT (1-cur)
    int                nre, teller, teller_disre, nfr;
    gmx_int64_t        start_step;
    int                nor = 0, nex = 0, norfr = 0, enx_i = 0;
    real               start_t;
    real              *bounds  = NULL, *violaver = NULL, *oobs = NULL, *orient = NULL, *odrms = NULL;
    int               *index   = NULL, *pair = NULL, norsel = 0, *orsel = NULL, *or_label = NULL;
    int                nbounds = 0, npairs;
    gmx_bool           bDisRe, bDRAll, bORA, bORT, bODA, bODR, bODT, bORIRE, bOTEN, bDHDL;
    gmx_bool           bFoundStart, bCont, bEDR, bVisco;
    double             sum, sumaver, sumt, ener, dbl;
    double            *time = NULL;
    real               Vaver;
    int               *set     = NULL, i, j, k, nset, sss;
    gmx_bool          *bIsEner = NULL;
    char             **pairleg, **odtleg, **otenleg;
    char             **leg = NULL;
    char             **nms;
    char              *anm_j, *anm_k, *resnm_j, *resnm_k;
    int                resnr_j, resnr_k;
    const char        *orinst_sub = "@ subtitle \"instantaneous\"\n";
    char               buf[256];
    output_env_t       oenv;
    t_enxblock        *blk       = NULL;
    t_enxblock        *blk_disre = NULL;
    int                ndisre    = 0;
    int                dh_blocks = 0, dh_hists = 0, dh_samples = 0, dh_lambdas = 0;

    t_filenm           fnm[] = {
        { efEDR, "-f",    NULL,      ffREAD  },
        { efEDR, "-f2",   NULL,      ffOPTRD },
        { efTPR, "-s",    NULL,      ffOPTRD },
        { efXVG, "-o",    "energy",  ffWRITE },
        { efXVG, "-viol", "violaver", ffOPTWR },
        { efXVG, "-pairs", "pairs",   ffOPTWR },
        { efXVG, "-ora",  "orienta", ffOPTWR },
        { efXVG, "-ort",  "orientt", ffOPTWR },
        { efXVG, "-oda",  "orideva", ffOPTWR },
        { efXVG, "-odr",  "oridevr", ffOPTWR },
        { efXVG, "-odt",  "oridevt", ffOPTWR },
        { efXVG, "-oten", "oriten",  ffOPTWR },
        { efXVG, "-corr", "enecorr", ffOPTWR },
        { efXVG, "-vis",  "visco",   ffOPTWR },
        { efXVG, "-ravg", "runavgdf", ffOPTWR },
        { efXVG, "-odh",  "dhdl", ffOPTWR }
    };
#define NFILE asize(fnm)
    int                npargs;
    t_pargs           *ppa;

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);
    if (!parse_common_args(&argc, argv,
                           PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END,
                           NFILE, fnm, npargs, ppa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    bDRAll = opt2bSet("-pairs", NFILE, fnm);
    bDisRe = opt2bSet("-viol", NFILE, fnm) || bDRAll;
    bORA   = opt2bSet("-ora", NFILE, fnm);
    bORT   = opt2bSet("-ort", NFILE, fnm);
    bODA   = opt2bSet("-oda", NFILE, fnm);
    bODR   = opt2bSet("-odr", NFILE, fnm);
    bODT   = opt2bSet("-odt", NFILE, fnm);
    bORIRE = bORA || bORT || bODA || bODR || bODT;
    bOTEN  = opt2bSet("-oten", NFILE, fnm);
    bDHDL  = opt2bSet("-odh", NFILE, fnm);

    nset = 0;

    snew(frame, 2);
    fp = open_enx(ftp2fn(efEDR, NFILE, fnm), "r");
    do_enxnms(fp, &nre, &enm);

    Vaver = -1;

    bVisco = opt2bSet("-vis", NFILE, fnm);

    if ((!bDisRe) && (!bDHDL))
    {
        if (bVisco)
        {
            nset = asize(setnm);
            snew(set, nset);
            /* This is nasty code... To extract Pres tensor, Volume and Temperature */
            for (j = 0; j < nset; j++)
            {
                for (i = 0; i < nre; i++)
                {
                    if (strstr(enm[i].name, setnm[j]))
                    {
                        set[j] = i;
                        break;
                    }
                }
                if (i == nre)
                {
                    if (gmx_strcasecmp(setnm[j], "Volume") == 0)
                    {
                        printf("Enter the box volume (" unit_volume "): ");
                        if (1 != scanf("%lf", &dbl))
                        {
                            gmx_fatal(FARGS, "Error reading user input");
                        }
                        Vaver = dbl;
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Could not find term %s for viscosity calculation",
                                  setnm[j]);
                    }
                }
            }
        }
        else
        {
            set = select_by_name(nre, enm, &nset);
        }
        /* Print all the different units once */
        sprintf(buf, "(%s)", enm[set[0]].unit);
        for (i = 1; i < nset; i++)
        {
            for (j = 0; j < i; j++)
            {
                if (strcmp(enm[set[i]].unit, enm[set[j]].unit) == 0)
                {
                    break;
                }
            }
            if (j == i)
            {
                strcat(buf, ", (");
                strcat(buf, enm[set[i]].unit);
                strcat(buf, ")");
            }
        }
        out = xvgropen(opt2fn("-o", NFILE, fnm), "GROMACS Energies", "Time (ps)", buf,
                       oenv);

        snew(leg, nset+1);
        for (i = 0; (i < nset); i++)
        {
            leg[i] = enm[set[i]].name;
        }
        if (bSum)
        {
            leg[nset] = gmx_strdup("Sum");
            xvgr_legend(out, nset+1, (const char**)leg, oenv);
        }
        else
        {
            xvgr_legend(out, nset, (const char**)leg, oenv);
        }

        snew(bIsEner, nset);
        for (i = 0; (i < nset); i++)
        {
            bIsEner[i] = FALSE;
            for (j = 0; (j <= F_ETOT); j++)
            {
                bIsEner[i] = bIsEner[i] ||
                    (gmx_strcasecmp(interaction_function[j].longname, leg[i]) == 0);
            }
        }

        if (bPrAll && nset > 1)
        {
            gmx_fatal(FARGS, "Printing averages can only be done when a single set is selected");
        }

        time = NULL;

        if (bORIRE || bOTEN)
        {
            get_orires_parms(ftp2fn(efTPR, NFILE, fnm), &nor, &nex, &or_label, &oobs);
        }

        if (bORIRE)
        {
            if (bOrinst)
            {
                enx_i = enxORI;
            }
            else
            {
                enx_i = enxOR;
            }

            if (bORA || bODA)
            {
                snew(orient, nor);
            }
            if (bODR)
            {
                snew(odrms, nor);
            }
            if (bORT || bODT)
            {
                fprintf(stderr, "Select the orientation restraint labels you want (-1 is all)\n");
                fprintf(stderr, "End your selection with 0\n");
                j     = -1;
                orsel = NULL;
                do
                {
                    j++;
                    srenew(orsel, j+1);
                    if (1 != scanf("%d", &(orsel[j])))
                    {
                        gmx_fatal(FARGS, "Error reading user input");
                    }
                }
                while (orsel[j] > 0);
                if (orsel[0] == -1)
                {
                    fprintf(stderr, "Selecting all %d orientation restraints\n", nor);
                    norsel = nor;
                    srenew(orsel, nor);
                    for (i = 0; i < nor; i++)
                    {
                        orsel[i] = i;
                    }
                }
                else
                {
                    /* Build the selection */
                    norsel = 0;
                    for (i = 0; i < j; i++)
                    {
                        for (k = 0; k < nor; k++)
                        {
                            if (or_label[k] == orsel[i])
                            {
                                orsel[norsel] = k;
                                norsel++;
                                break;
                            }
                        }
                        if (k == nor)
                        {
                            fprintf(stderr, "Orientation restraint label %d not found\n",
                                    orsel[i]);
                        }
                    }
                }
                snew(odtleg, norsel);
                for (i = 0; i < norsel; i++)
                {
                    snew(odtleg[i], 256);
                    sprintf(odtleg[i], "%d", or_label[orsel[i]]);
                }
                if (bORT)
                {
                    fort = xvgropen(opt2fn("-ort", NFILE, fnm), "Calculated orientations",
                                    "Time (ps)", "", oenv);
                    if (bOrinst && output_env_get_print_xvgr_codes(oenv))
                    {
                        fprintf(fort, "%s", orinst_sub);
                    }
                    xvgr_legend(fort, norsel, (const char**)odtleg, oenv);
                }
                if (bODT)
                {
                    fodt = xvgropen(opt2fn("-odt", NFILE, fnm),
                                    "Orientation restraint deviation",
                                    "Time (ps)", "", oenv);
                    if (bOrinst && output_env_get_print_xvgr_codes(oenv))
                    {
                        fprintf(fodt, "%s", orinst_sub);
                    }
                    xvgr_legend(fodt, norsel, (const char**)odtleg, oenv);
                }
            }
        }
        if (bOTEN)
        {
            foten = xvgropen(opt2fn("-oten", NFILE, fnm),
                             "Order tensor", "Time (ps)", "", oenv);
            snew(otenleg, bOvec ? nex*12 : nex*3);
            for (i = 0; i < nex; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    sprintf(buf, "eig%d", j+1);
                    otenleg[(bOvec ? 12 : 3)*i+j] = gmx_strdup(buf);
                }
                if (bOvec)
                {
                    for (j = 0; j < 9; j++)
                    {
                        sprintf(buf, "vec%d%s", j/3+1, j%3 == 0 ? "x" : (j%3 == 1 ? "y" : "z"));
                        otenleg[12*i+3+j] = gmx_strdup(buf);
                    }
                }
            }
            xvgr_legend(foten, bOvec ? nex*12 : nex*3, (const char**)otenleg, oenv);
        }
    }
    else if (bDisRe)
    {
        nbounds = get_bounds(ftp2fn(efTPR, NFILE, fnm), &bounds, &index, &pair, &npairs,
                             &mtop, &top, &ir);
        snew(violaver, npairs);
        out = xvgropen(opt2fn("-o", NFILE, fnm), "Sum of Violations",
                       "Time (ps)", "nm", oenv);
        xvgr_legend(out, 2, drleg, oenv);
        if (bDRAll)
        {
            fp_pairs = xvgropen(opt2fn("-pairs", NFILE, fnm), "Pair Distances",
                                "Time (ps)", "Distance (nm)", oenv);
            if (output_env_get_print_xvgr_codes(oenv))
            {
                fprintf(fp_pairs, "@ subtitle \"averaged (tau=%g) and instantaneous\"\n",
                        ir.dr_tau);
            }
        }
    }
    else if (bDHDL)
    {
        get_dhdl_parms(ftp2fn(efTPR, NFILE, fnm), &ir);
    }

    /* Initiate energies and set them to zero */
    edat.nsteps    = 0;
    edat.npoints   = 0;
    edat.nframes   = 0;
    edat.step      = NULL;
    edat.steps     = NULL;
    edat.points    = NULL;
    edat.bHaveSums = TRUE;
    snew(edat.s, nset);

    /* Initiate counters */
    teller       = 0;
    teller_disre = 0;
    bFoundStart  = FALSE;
    start_step   = 0;
    start_t      = 0;
    do
    {
        /* This loop searches for the first frame (when -b option is given),
         * or when this has been found it reads just one energy frame
         */
        do
        {
            bCont = do_enx(fp, &(frame[NEXT]));
            if (bCont)
            {
                timecheck = check_times(frame[NEXT].t);
            }
        }
        while (bCont && (timecheck < 0));

        if ((timecheck == 0) && bCont)
        {
            /* We read a valid frame, so we can use it */
            fr = &(frame[NEXT]);

            if (fr->nre > 0)
            {
                /* The frame contains energies, so update cur */
                cur  = NEXT;

                if (edat.nframes % 1000 == 0)
                {
                    srenew(edat.step, edat.nframes+1000);
                    memset(&(edat.step[edat.nframes]), 0, 1000*sizeof(edat.step[0]));
                    srenew(edat.steps, edat.nframes+1000);
                    memset(&(edat.steps[edat.nframes]), 0, 1000*sizeof(edat.steps[0]));
                    srenew(edat.points, edat.nframes+1000);
                    memset(&(edat.points[edat.nframes]), 0, 1000*sizeof(edat.points[0]));

                    for (i = 0; i < nset; i++)
                    {
                        srenew(edat.s[i].ener, edat.nframes+1000);
                        memset(&(edat.s[i].ener[edat.nframes]), 0,
                               1000*sizeof(edat.s[i].ener[0]));
                        srenew(edat.s[i].es, edat.nframes+1000);
                        memset(&(edat.s[i].es[edat.nframes]), 0,
                               1000*sizeof(edat.s[i].es[0]));
                    }
                }

                nfr            = edat.nframes;
                edat.step[nfr] = fr->step;

                if (!bFoundStart)
                {
                    bFoundStart = TRUE;
                    /* Initiate the previous step data */
                    start_step = fr->step;
                    start_t    = fr->t;
                    /* Initiate the energy sums */
                    edat.steps[nfr]  = 1;
                    edat.points[nfr] = 1;
                    for (i = 0; i < nset; i++)
                    {
                        sss                    = set[i];
                        edat.s[i].es[nfr].sum  = fr->ener[sss].e;
                        edat.s[i].es[nfr].sum2 = 0;
                    }
                    edat.nsteps  = 1;
                    edat.npoints = 1;
                }
                else
                {
                    edat.steps[nfr] = fr->nsteps;

                    if (fr->nsum <= 1)
                    {
                        /* mdrun only calculated the energy at energy output
                         * steps. We don't need to check step intervals.
                         */
                        edat.points[nfr] = 1;
                        for (i = 0; i < nset; i++)
                        {
                            sss                    = set[i];
                            edat.s[i].es[nfr].sum  = fr->ener[sss].e;
                            edat.s[i].es[nfr].sum2 = 0;
                        }
                        edat.npoints  += 1;
                        edat.bHaveSums = FALSE;
                    }
                    else if (fr->step - start_step + 1 == edat.nsteps + fr->nsteps)
                    {
                        /* We have statistics  to the previous frame */
                        edat.points[nfr] = fr->nsum;
                        for (i = 0; i < nset; i++)
                        {
                            sss                    = set[i];
                            edat.s[i].es[nfr].sum  = fr->ener[sss].esum;
                            edat.s[i].es[nfr].sum2 = fr->ener[sss].eav;
                        }
                        edat.npoints += fr->nsum;
                    }
                    else
                    {
                        /* The interval does not match fr->nsteps:
                         * can not do exact averages.
                         */
                        edat.bHaveSums = FALSE;
                    }

                    edat.nsteps = fr->step - start_step + 1;
                }
                for (i = 0; i < nset; i++)
                {
                    edat.s[i].ener[nfr] = fr->ener[set[i]].e;
                }
            }
            /*
             * Define distance restraint legends. Can only be done after
             * the first frame has been read... (Then we know how many there are)
             */
            blk_disre = find_block_id_enxframe(fr, enxDISRE, NULL);
            if (bDisRe && bDRAll && !leg && blk_disre)
            {
                t_iatom   *fa;
                t_iparams *ip;

                fa = top->idef.il[F_DISRES].iatoms;
                ip = top->idef.iparams;
                if (blk_disre->nsub != 2 ||
                    (blk_disre->sub[0].nr != blk_disre->sub[1].nr) )
                {
                    gmx_incons("Number of disre sub-blocks not equal to 2");
                }

                ndisre = blk_disre->sub[0].nr;
                if (ndisre != top->idef.il[F_DISRES].nr/3)
                {
                    gmx_fatal(FARGS, "Number of disre pairs in the energy file (%d) does not match the number in the run input file (%d)\n",
                              ndisre, top->idef.il[F_DISRES].nr/3);
                }
                snew(pairleg, ndisre);
                for (i = 0; i < ndisre; i++)
                {
                    snew(pairleg[i], 30);
                    j = fa[3*i+1];
                    k = fa[3*i+2];
                    gmx_mtop_atominfo_global(&mtop, j, &anm_j, &resnr_j, &resnm_j);
                    gmx_mtop_atominfo_global(&mtop, k, &anm_k, &resnr_k, &resnm_k);
                    sprintf(pairleg[i], "%d %s %d %s (%d)",
                            resnr_j, anm_j, resnr_k, anm_k,
                            ip[fa[3*i]].disres.label);
                }
                set = select_it(ndisre, pairleg, &nset);
                snew(leg, 2*nset);
                for (i = 0; (i < nset); i++)
                {
                    snew(leg[2*i], 32);
                    sprintf(leg[2*i],  "a %s", pairleg[set[i]]);
                    snew(leg[2*i+1], 32);
                    sprintf(leg[2*i+1], "i %s", pairleg[set[i]]);
                }
                xvgr_legend(fp_pairs, 2*nset, (const char**)leg, oenv);
            }

            /*
             * Store energies for analysis afterwards...
             */
            if (!bDisRe && !bDHDL && (fr->nre > 0))
            {
                if (edat.nframes % 1000 == 0)
                {
                    srenew(time, edat.nframes+1000);
                }
                time[edat.nframes] = fr->t;
                edat.nframes++;
            }
            /*
             * Printing time, only when we do not want to skip frames
             */
            if (!skip || teller % skip == 0)
            {
                if (bDisRe)
                {
                    /*******************************************
                     * D I S T A N C E   R E S T R A I N T S
                     *******************************************/
                    if (ndisre > 0)
                    {
 #ifndef GMX_DOUBLE
                        float  *disre_rt     =     blk_disre->sub[0].fval;
                        float  *disre_rm3tav = blk_disre->sub[1].fval;
 #else
                        double *disre_rt     =     blk_disre->sub[0].dval;
                        double *disre_rm3tav = blk_disre->sub[1].dval;
 #endif

                        print_time(out, fr->t);
                        if (violaver == NULL)
                        {
                            snew(violaver, ndisre);
                        }

                        /* Subtract bounds from distances, to calculate violations */
                        calc_violations(disre_rt, disre_rm3tav,
                                        nbounds, pair, bounds, violaver, &sumt, &sumaver);

                        fprintf(out, "  %8.4f  %8.4f\n", sumaver, sumt);
                        if (bDRAll)
                        {
                            print_time(fp_pairs, fr->t);
                            for (i = 0; (i < nset); i++)
                            {
                                sss = set[i];
                                fprintf(fp_pairs, "  %8.4f", mypow(disre_rm3tav[sss], minthird));
                                fprintf(fp_pairs, "  %8.4f", disre_rt[sss]);
                            }
                            fprintf(fp_pairs, "\n");
                        }
                        teller_disre++;
                    }
                }
                else if (bDHDL)
                {
                    do_dhdl(fr, &ir, &fp_dhdl, opt2fn("-odh", NFILE, fnm), bDp, &dh_blocks, &dh_hists, &dh_samples, &dh_lambdas, oenv);
                }

                /*******************************************
                 * E N E R G I E S
                 *******************************************/
                else
                {
                    if (fr->nre > 0)
                    {
                        if (bPrAll)
                        {
                            /* We skip frames with single points (usually only the first frame),
                             * since they would result in an average plot with outliers.
                             */
                            if (fr->nsum > 1)
                            {
                                print_time(out, fr->t);
                                print1(out, bDp, fr->ener[set[0]].e);
                                print1(out, bDp, fr->ener[set[0]].esum/fr->nsum);
                                print1(out, bDp, sqrt(fr->ener[set[0]].eav/fr->nsum));
                                fprintf(out, "\n");
                            }
                        }
                        else
                        {
                            print_time(out, fr->t);
                            if (bSum)
                            {
                                sum = 0;
                                for (i = 0; i < nset; i++)
                                {
                                    sum += fr->ener[set[i]].e;
                                }
                                print1(out, bDp, sum/nmol-ezero);
                            }
                            else
                            {
                                for (i = 0; (i < nset); i++)
                                {
                                    if (bIsEner[i])
                                    {
                                        print1(out, bDp, (fr->ener[set[i]].e)/nmol-ezero);
                                    }
                                    else
                                    {
                                        print1(out, bDp, fr->ener[set[i]].e);
                                    }
                                }
                            }
                            fprintf(out, "\n");
                        }
                    }
                    blk = find_block_id_enxframe(fr, enx_i, NULL);
                    if (bORIRE && blk)
                    {
#ifndef GMX_DOUBLE
                        xdr_datatype dt = xdr_datatype_float;
#else
                        xdr_datatype dt = xdr_datatype_double;
#endif
                        real        *vals;

                        if ( (blk->nsub != 1) || (blk->sub[0].type != dt) )
                        {
                            gmx_fatal(FARGS, "Orientational restraints read in incorrectly");
                        }
#ifndef GMX_DOUBLE
                        vals = blk->sub[0].fval;
#else
                        vals = blk->sub[0].dval;
#endif

                        if (blk->sub[0].nr != (size_t)nor)
                        {
                            gmx_fatal(FARGS, "Number of orientation restraints in energy file (%d) does not match with the topology (%d)", blk->sub[0].nr);
                        }
                        if (bORA || bODA)
                        {
                            for (i = 0; i < nor; i++)
                            {
                                orient[i] += vals[i];
                            }
                        }
                        if (bODR)
                        {
                            for (i = 0; i < nor; i++)
                            {
                                odrms[i] += sqr(vals[i]-oobs[i]);
                            }
                        }
                        if (bORT)
                        {
                            fprintf(fort, "  %10f", fr->t);
                            for (i = 0; i < norsel; i++)
                            {
                                fprintf(fort, " %g", vals[orsel[i]]);
                            }
                            fprintf(fort, "\n");
                        }
                        if (bODT)
                        {
                            fprintf(fodt, "  %10f", fr->t);
                            for (i = 0; i < norsel; i++)
                            {
                                fprintf(fodt, " %g", vals[orsel[i]]-oobs[orsel[i]]);
                            }
                            fprintf(fodt, "\n");
                        }
                        norfr++;
                    }
                    blk = find_block_id_enxframe(fr, enxORT, NULL);
                    if (bOTEN && blk)
                    {
#ifndef GMX_DOUBLE
                        xdr_datatype dt = xdr_datatype_float;
#else
                        xdr_datatype dt = xdr_datatype_double;
#endif
                        real        *vals;

                        if ( (blk->nsub != 1) || (blk->sub[0].type != dt) )
                        {
                            gmx_fatal(FARGS, "Orientational restraints read in incorrectly");
                        }
#ifndef GMX_DOUBLE
                        vals = blk->sub[0].fval;
#else
                        vals = blk->sub[0].dval;
#endif

                        if (blk->sub[0].nr != (size_t)(nex*12))
                        {
                            gmx_fatal(FARGS, "Number of orientation experiments in energy file (%g) does not match with the topology (%d)",
                                      blk->sub[0].nr/12, nex);
                        }
                        fprintf(foten, "  %10f", fr->t);
                        for (i = 0; i < nex; i++)
                        {
                            for (j = 0; j < (bOvec ? 12 : 3); j++)
                            {
                                fprintf(foten, " %g", vals[i*12+j]);
                            }
                        }
                        fprintf(foten, "\n");
                    }
                }
            }
            teller++;
        }
    }
    while (bCont && (timecheck == 0));

    fprintf(stderr, "\n");
    close_enx(fp);
    if (out)
    {
        xvgrclose(out);
    }

    if (bDRAll)
    {
        xvgrclose(fp_pairs);
    }

    if (bORT)
    {
        xvgrclose(fort);
    }
    if (bODT)
    {
        xvgrclose(fodt);
    }
    if (bORA)
    {
        out = xvgropen(opt2fn("-ora", NFILE, fnm),
                       "Average calculated orientations",
                       "Restraint label", "", oenv);
        if (bOrinst && output_env_get_print_xvgr_codes(oenv))
        {
            fprintf(out, "%s", orinst_sub);
        }
        for (i = 0; i < nor; i++)
        {
            fprintf(out, "%5d  %g\n", or_label[i], orient[i]/norfr);
        }
        xvgrclose(out);
    }
    if (bODA)
    {
        out = xvgropen(opt2fn("-oda", NFILE, fnm),
                       "Average restraint deviation",
                       "Restraint label", "", oenv);
        if (bOrinst && output_env_get_print_xvgr_codes(oenv))
        {
            fprintf(out, "%s", orinst_sub);
        }
        for (i = 0; i < nor; i++)
        {
            fprintf(out, "%5d  %g\n", or_label[i], orient[i]/norfr-oobs[i]);
        }
        xvgrclose(out);
    }
    if (bODR)
    {
        out = xvgropen(opt2fn("-odr", NFILE, fnm),
                       "RMS orientation restraint deviations",
                       "Restraint label", "", oenv);
        if (bOrinst && output_env_get_print_xvgr_codes(oenv))
        {
            fprintf(out, "%s", orinst_sub);
        }
        for (i = 0; i < nor; i++)
        {
            fprintf(out, "%5d  %g\n", or_label[i], sqrt(odrms[i]/norfr));
        }
        xvgrclose(out);
    }
    if (bOTEN)
    {
        xvgrclose(foten);
    }

    if (bDisRe)
    {
        analyse_disre(opt2fn("-viol", NFILE, fnm),
                      teller_disre, violaver, bounds, index, pair, nbounds, oenv);
    }
    else if (bDHDL)
    {
        if (fp_dhdl)
        {
            gmx_fio_fclose(fp_dhdl);
            printf("\n\nWrote %d lambda values with %d samples as ",
                   dh_lambdas, dh_samples);
            if (dh_hists > 0)
            {
                printf("%d dH histograms ", dh_hists);
            }
            if (dh_blocks > 0)
            {
                printf("%d dH data blocks ", dh_blocks);
            }
            printf("to %s\n", opt2fn("-odh", NFILE, fnm));

        }
        else
        {
            gmx_fatal(FARGS, "No dH data in %s\n", opt2fn("-f", NFILE, fnm));
        }

    }
    else
    {
        double dt = (frame[cur].t-start_t)/(edat.nframes-1);
        analyse_ener(opt2bSet("-corr", NFILE, fnm), opt2fn("-corr", NFILE, fnm),
                     bFee, bSum, opt2parg_bSet("-nmol", npargs, ppa),
                     bVisco, opt2fn("-vis", NFILE, fnm),
                     nmol,
                     start_step, start_t, frame[cur].step, frame[cur].t,
                     reftemp, &edat,
                     nset, set, bIsEner, leg, enm, Vaver, ezero, nbmin, nbmax,
                     oenv);
        if (bFluctProps)
        {
            calc_fluctuation_props(stdout, bDriftCorr, dt, nset, nmol, leg, &edat,
                                   nbmin, nbmax);
        }
    }
    if (opt2bSet("-f2", NFILE, fnm))
    {
        fec(opt2fn("-f2", NFILE, fnm), opt2fn("-ravg", NFILE, fnm),
            reftemp, nset, set, leg, &edat, time, oenv);
    }

    {
        const char *nxy = "-nxy";

        do_view(oenv, opt2fn("-o", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-ravg", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-ora", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-ort", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-oda", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-odr", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-odt", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-oten", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-odh", NFILE, fnm), nxy);
    }

    return 0;
}

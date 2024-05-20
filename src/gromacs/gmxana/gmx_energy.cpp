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

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

static const int NOTSET = -23451;

typedef struct
{
    real sum;
    real sum2;
} exactsum_t;

typedef struct
{
    real*       ener;
    exactsum_t* es;
    gmx_bool    bExactStat;
    double      av;
    double      rmsd;
    double      ee;
    double      slope;
} enerdat_t;

typedef struct
{
    int64_t    nsteps;
    int64_t    npoints;
    int        nframes;
    int*       step;
    int*       steps;
    int*       points;
    enerdat_t* s;
    gmx_bool   bHaveSums;
} enerdata_t;

static void done_enerdata_t(int nset, enerdata_t* edat)
{
    sfree(edat->step);
    sfree(edat->steps);
    sfree(edat->points);
    for (int i = 0; i < nset; i++)
    {
        sfree(edat->s[i].ener);
        sfree(edat->s[i].es);
    }
    sfree(edat->s);
}

static void chomp(char* buf)
{
    int len = std::strlen(buf);

    while ((len > 0) && (buf[len - 1] == '\n'))
    {
        buf[len - 1] = '\0';
        len--;
    }
}

static int* select_by_name(int nre, gmx_enxnm_t* nm, int* nset)
{
    gmx_bool*   bE;
    int         k, kk, j, i, nmatch, nind, nss;
    int*        set;
    gmx_bool    bEOF, bVerbose = TRUE, bLong = FALSE;
    char *      ptr, buf[STRLEN];
    const char* fm4   = "%3d  %-14s";
    const char* fm2   = "%3d  %-34s";
    char**      newnm = nullptr;

    if ((getenv("GMX_ENER_VERBOSE")) != nullptr)
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
        while ((ptr = std::strchr(newnm[k], ' ')) != nullptr)
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
                for (kk = k; kk < k + 4; kk++)
                {
                    if (kk < nre && std::strlen(nm[kk].name) > 14)
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
                fprintf(stderr, fm4, k + 1, newnm[k]);
                j++;
                if (j == 4)
                {
                    j = 0;
                }
            }
            else
            {
                fprintf(stderr, fm2, k + 1, newnm[k]);
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
    while (!bEOF && (fgets2(buf, STRLEN - 1, stdin)))
    {
        /* Remove newlines */
        chomp(buf);

        /* Remove spaces */
        trim(buf);

        /* Empty line means end of input */
        bEOF = (std::strlen(buf) == 0);
        if (!bEOF)
        {
            ptr = buf;
            do
            {
                if (!bEOF)
                {
                    /* First, try to match an exact field name */
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
                        /* Second, try to read an integer */
                        nss = sscanf(ptr, "%d", &nind);
                        if (nss == 1)
                        {
                            /* Zero means end of input */
                            if (nind == 0)
                            {
                                bEOF = TRUE;
                            }
                            else if ((1 <= nind) && (nind <= nre))
                            {
                                bE[nind - 1] = TRUE;
                            }
                            else
                            {
                                fprintf(stderr, "number %d is out of range\n", nind);
                            }
                        }
                        else
                        {
                            /* Finally, match on part of the field name */
                            i      = std::strlen(ptr);
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
                if ((ptr = std::strchr(ptr, ' ')) != nullptr)
                {
                    trim(ptr);
                }
            } while (!bEOF && ((ptr != nullptr) && (std::strlen(ptr) > 0)));
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

static void get_dhdl_parms(const char* topnm, t_inputrec* ir)
{
    gmx_mtop_t mtop;
    int        natoms;
    matrix     box;

    /* all we need is the ir to be able to write the label */
    read_tpx(topnm, ir, box, &natoms, nullptr, nullptr, &mtop);
}

// Computes and writes the shear viscosity using the Einstein relation
static void einstein_visco(const char*             fn,
                           const char*             fni,
                           int                     nsets,
                           const enerdata_t&       edat,
                           const real              volume,
                           const real              temperature,
                           const int               numRestarts,
                           double                  dt,
                           const gmx_output_env_t* oenv)
{
    constexpr int c_numSets = 3;

    GMX_RELEASE_ASSERT(nsets == c_numSets, "Only nsets=3 is currently supported");

    /* Determine integrals of the off-diagonal pressure elements */
    const int                          nint = edat.nframes + 1;
    std::array<std::vector<double>, 3> eneint;
    for (int i = 0; i < c_numSets; i++)
    {
        eneint[i].resize(nint, 0.0);
    }
    for (int i = 0; i < edat.nframes; i++)
    {
        const double fac = dt / edat.points[i];
        eneint[0][i + 1] = eneint[0][i] + 0.5 * (edat.s[1].es[i].sum + edat.s[3].es[i].sum) * fac;
        eneint[1][i + 1] = eneint[1][i] + 0.5 * (edat.s[2].es[i].sum + edat.s[6].es[i].sum) * fac;
        eneint[2][i + 1] = eneint[2][i] + 0.5 * (edat.s[5].es[i].sum + edat.s[7].es[i].sum) * fac;
    }

    const int nf4 = nint / 4 + 1;

    if (numRestarts <= 0)
    {
        GMX_THROW(
                gmx::InvalidInputError("The number of restarts for computing the viscosity using "
                                       "Einstein should be positive"));
    }

    const int stepSize = std::max(nf4 / numRestarts, 1);

    printf("\n");
    printf("Computing shear viscosity using the Einstein relation with %d start points separated "
           "by %g ps\n",
           (nf4 + stepSize - 1) / stepSize,
           stepSize * dt);


    std::array<double, c_numSets + 1> avold = { 0.0 };

    FILE* fp0 = xvgropen(fni, "Shear viscosity integral", "Time (ps)", "(kg m\\S-1\\N s\\S-1\\N ps)", oenv);
    FILE* fp1 = xvgropen(
            fn, "Shear viscosity using Einstein relation", "Time (ps)", "(kg m\\S-1\\N s\\S-1\\N)", oenv);
    for (int i = 0; i < nf4; i += stepSize)
    {
        std::array<double, c_numSets + 1> av = { 0.0 };

        for (int m = 0; m < c_numSets; m++)
        {
            for (int j = 0; j < nint - i; j++)
            {
                double di = gmx::square(eneint[m][j + i] - eneint[m][j]);

                av[m] += di;
                av[c_numSets] += di / c_numSets;
            }
        }
        /* Convert to SI for the viscosity */
        const double fac = (volume * gmx::c_nano * gmx::c_nano * gmx::c_nano * gmx::c_pico * 1e10)
                           / (2 * gmx::c_boltzmann * temperature) / (nint - i);
        fprintf(fp0, "%10g", i * dt);
        for (int m = 0; m <= c_numSets; m++)
        {
            av[m] = fac * av[m];
            fprintf(fp0, "  %10g", av[m]);
        }
        fprintf(fp0, "\n");
        fprintf(fp1, "%10g", (i + 0.5) * dt);
        for (int m = 0; m <= c_numSets; m++)
        {
            fprintf(fp1, "  %10g", (av[m] - avold[m]) / (stepSize * dt));
            avold[m] = av[m];
        }
        fprintf(fp1, "\n");
    }
    xvgrclose(fp0);
    xvgrclose(fp1);
}

typedef struct
{
    int64_t np;
    double  sum;
    double  sav;
    double  sav2;
} ee_sum_t;

typedef struct
{
    int      b;
    ee_sum_t sum;
    int64_t  nst;
    int64_t  nst_min;
} ener_ee_t;

static void clear_ee_sum(ee_sum_t* ees)
{
    ees->sav  = 0;
    ees->sav2 = 0;
    ees->np   = 0;
    ees->sum  = 0;
}

static void add_ee_sum(ee_sum_t* ees, double sum, int np)
{
    ees->np += np;
    ees->sum += sum;
}

static void add_ee_av(ee_sum_t* ees)
{
    double av;

    av = ees->sum / ees->np;
    ees->sav += av;
    ees->sav2 += av * av;
    ees->np  = 0;
    ees->sum = 0;
}

static double calc_ee2(int nb, ee_sum_t* ees)
{
    return (ees->sav2 / nb - gmx::square(ees->sav / nb)) / (nb - 1);
}

static void set_ee_av(ener_ee_t* eee)
{
    if (debug)
    {
        char buf[STEPSTRSIZE];
        fprintf(debug, "Storing average for err.est.: %s steps\n", gmx_step_str(eee->nst, buf));
    }
    add_ee_av(&eee->sum);
    eee->b++;
    if (eee->b == 1 || eee->nst < eee->nst_min)
    {
        eee->nst_min = eee->nst;
    }
    eee->nst = 0;
}

static void calc_averages(int nset, enerdata_t* edat, int nbmin, int nbmax)
{
    int         nb, i, f, nee;
    double      sum, sum2, sump, see2;
    int64_t     np, p, bound_nb;
    enerdat_t*  ed;
    exactsum_t* es;
    gmx_bool    bAllZero;
    double      x, sx, sy, sxx, sxy;
    ener_ee_t*  eee;

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

    snew(eee, nbmax + 1);
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
            eee[nb].b = 0;
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
                p    = edat->points[f];
                sump = es->sum;
                sum2 += es->sum2;
                if (np > 0)
                {
                    sum2 += gmx::square(sum / np - (sum + es->sum) / (np + p)) * np * (np + p) / p;
                }
            }
            else
            {
                /* Add a single value to the sum and sum of squares. */
                p    = 1;
                sump = ed->ener[f];
                sum2 += gmx::square(sump);
            }

            /* sum has to be increased after sum2 */
            np += p;
            sum += sump;

            /* For the linear regression use variance 1/p.
             * Note that sump is the sum, not the average, so we don't need p*.
             */
            x = edat->step[f] - 0.5 * (edat->steps[f] - 1);
            sx += p * x;
            sy += sump;
            sxx += p * x * x;
            sxy += x * sump;

            for (nb = nbmin; nb <= nbmax; nb++)
            {
                /* Check if the current end step is closer to the desired
                 * block boundary than the next end step.
                 */
                bound_nb = (edat->step[0] - 1) * nb + edat->nsteps * (eee[nb].b + 1);
                if (eee[nb].nst > 0 && bound_nb - edat->step[f - 1] * nb < edat->step[f] * nb - bound_nb)
                {
                    set_ee_av(&eee[nb]);
                }
                if (f == 0)
                {
                    eee[nb].nst = 1;
                }
                else
                {
                    eee[nb].nst += edat->step[f] - edat->step[f - 1];
                }
                if (ed->bExactStat)
                {
                    add_ee_sum(&eee[nb].sum, es->sum, edat->points[f]);
                }
                else
                {
                    add_ee_sum(&eee[nb].sum, edat->s[i].ener[f], 1);
                }
                bound_nb = (edat->step[0] - 1) * nb + edat->nsteps * (eee[nb].b + 1);
                if (edat->step[f] * nb >= bound_nb)
                {
                    set_ee_av(&eee[nb]);
                }
            }
        }

        edat->s[i].av = sum / np;
        if (ed->bExactStat)
        {
            edat->s[i].rmsd = std::sqrt(sum2 / np);
        }
        else
        {
            edat->s[i].rmsd = std::sqrt(std::max(sum2 / np - gmx::square(edat->s[i].av), 0.0));
        }

        if (edat->nframes > 1)
        {
            edat->s[i].slope = (np * sxy - sx * sy) / (np * sxx - sx * sx);
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
                fprintf(debug,
                        "Requested %d blocks, we have %d blocks, min %s nsteps %s\n",
                        nb,
                        eee[nb].b,
                        gmx_step_str(eee[nb].nst_min, buf1),
                        gmx_step_str(edat->nsteps, buf2));
            }
            if (eee[nb].b == nb && 5 * nb * eee[nb].nst_min >= 4 * edat->nsteps)
            {
                see2 += calc_ee2(nb, &eee[nb].sum);
                nee++;
            }
        }
        if (nee > 0)
        {
            edat->s[i].ee = std::sqrt(see2 / nee);
        }
        else
        {
            edat->s[i].ee = -1;
        }
    }
    sfree(eee);
}

static enerdata_t* calc_sum(int nset, enerdata_t* edat, int nbmin, int nbmax)
{
    enerdata_t* esum;
    enerdat_t*  s;
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

static void ee_pr(double ee, int buflen, char* buf)
{
    snprintf(buf, buflen, "%s", "--");
    if (ee >= 0)
    {
        /* Round to two decimals by printing. */
        char tmp[100];
        snprintf(tmp, sizeof(tmp), "%.1e", ee);
        double rnd = gmx::doubleFromString(tmp);
        snprintf(buf, buflen, "%g", rnd);
    }
}

static void remove_drift(int nset, int nbmin, int nbmax, real dt, enerdata_t* edat)
{
    /* Remove the drift by performing a fit to y = ax+b.
       Uses 5 iterations. */
    int    i, j, k;
    double delta;

    edat->npoints = edat->nframes;
    edat->nsteps  = edat->nframes;

    for (k = 0; (k < 5); k++)
    {
        for (i = 0; (i < nset); i++)
        {
            delta = edat->s[i].slope * dt;

            if (nullptr != debug)
            {
                fprintf(debug, "slope for set %d is %g\n", i, edat->s[i].slope);
            }

            for (j = 0; (j < edat->nframes); j++)
            {
                edat->s[i].ener[j] -= j * delta;
                edat->s[i].es[j].sum  = 0;
                edat->s[i].es[j].sum2 = 0;
            }
        }
        calc_averages(nset, edat, nbmin, nbmax);
    }
}

static void calc_fluctuation_props(FILE*                            fp,
                                   gmx_bool                         bDriftCorr,
                                   real                             dt,
                                   int                              nset,
                                   int                              nmol,
                                   gmx::ArrayRef<const std::string> leg,
                                   enerdata_t*                      edat,
                                   int                              nbmin,
                                   int                              nbmax)
{
    int    i, j;
    double vv, v, h, varv, hh, varh, tt, cv, cp, alpha, kappa, dcp, varet;
    double NANO3;
    enum
    {
        eVol,
        eEnth,
        eTemp,
        eEtot,
        eNR
    };
    const char* my_ener[] = { "Volume", "Enthalpy", "Temperature", "Total Energy" };
    int         ii[eNR];

    NANO3 = gmx::c_nano * gmx::c_nano * gmx::c_nano;
    if (!bDriftCorr)
    {
        fprintf(fp,
                "\nYou may want to use the -driftcorr flag in order to correct\n"
                "for spurious drift in the graphs. Note that this is not\n"
                "a substitute for proper equilibration and sampling!\n");
    }
    else
    {
        remove_drift(nset, nbmin, nbmax, dt, edat);
    }
    for (i = 0; (i < eNR); i++)
    {
        for (ii[i] = 0; (ii[i] < nset && (gmx_strcasecmp(leg[ii[i]].c_str(), my_ener[i]) != 0)); ii[i]++)
        {
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
    tt = edat->s[ii[eTemp]].av;
}
/* Volume */
vv = varv = NOTSET;
if ((ii[eVol] < nset) && (ii[eTemp] < nset))
{
    vv    = edat->s[ii[eVol]].av * NANO3;
    varv  = gmx::square(edat->s[ii[eVol]].rmsd * NANO3);
    kappa = (varv / vv) / (gmx::c_boltzmann * tt);
}
/* Enthalpy */
hh = varh = NOTSET;
if ((ii[eEnth] < nset) && (ii[eTemp] < nset))
{
    hh   = gmx::c_kilo * edat->s[ii[eEnth]].av / gmx::c_avogadro;
    varh = gmx::square(gmx::c_kilo * edat->s[ii[eEnth]].rmsd / gmx::c_avogadro);
    cp   = gmx::c_avogadro * ((varh / nmol) / (gmx::c_boltzmann * tt * tt));
}
/* Total energy */
if ((ii[eEtot] < nset) && (hh == NOTSET) && (tt != NOTSET))
{
    /* Only compute cv in constant volume runs, which we can test
       by checking whether the enthalpy was computed.
     */
    varet = gmx::square(edat->s[ii[eEtot]].rmsd);
    cv    = gmx::c_kilo * ((varet / nmol) / (gmx::c_boltz * tt * tt));
}
/* Alpha, dcp */
if ((ii[eVol] < nset) && (ii[eEnth] < nset) && (ii[eTemp] < nset))
{
    double v_sum, h_sum, vh_sum, v_aver, h_aver, vh_aver;
    vh_sum = v_sum = h_sum = 0;
    for (j = 0; (j < edat->nframes); j++)
    {
        v = edat->s[ii[eVol]].ener[j] * NANO3;
        h = gmx::c_kilo * edat->s[ii[eEnth]].ener[j] / gmx::c_avogadro;
        v_sum += v;
        h_sum += h;
        vh_sum += (v * h);
    }
    vh_aver = vh_sum / edat->nframes;
    v_aver  = v_sum / edat->nframes;
    h_aver  = h_sum / edat->nframes;
    alpha   = (vh_aver - v_aver * h_aver) / (v_aver * gmx::c_boltzmann * tt * tt);
    dcp     = (v_aver * gmx::c_avogadro / nmol) * tt * gmx::square(alpha) / (kappa);
}

if (tt != NOTSET)
{
    if (nmol < 2)
    {
        fprintf(fp, "\nWARNING: nmol = %d, this may not be what you want.\n", nmol);
    }
    fprintf(fp, "\nTemperature dependent fluctuation properties at T = %g.\n", tt);
    fprintf(fp, "\nHeat capacities obtained from fluctuations do *not* include\n");
    fprintf(fp, "quantum corrections. If you want to get a more accurate estimate\n");
    fprintf(fp, "please use the gmx dos program.\n\n");
    fprintf(fp,
            "WARNING: Please verify that your simulations are converged and perform\n"
            "a block-averaging error analysis (not implemented in gmx energy yet)\n");

    if (debug != nullptr)
    {
        if (varv != NOTSET)
        {
            fprintf(fp, "varv  =  %10g (m^6)\n", varv * gmx::c_avogadro / nmol);
        }
    }
    if (vv != NOTSET)
    {
        fprintf(fp, "Volume                                   = %10g m^3/mol\n", vv * gmx::c_avogadro / nmol);
    }
    if (varh != NOTSET)
    {
        fprintf(fp,
                "Enthalpy                                 = %10g kJ/mol\n",
                hh * gmx::c_avogadro / (gmx::c_kilo * nmol));
    }
    if (alpha != NOTSET)
    {
        fprintf(fp, "Coefficient of Thermal Expansion Alpha_P = %10g (1/K)\n", alpha);
    }
    if (kappa != NOTSET)
    {
        fprintf(fp, "Isothermal Compressibility Kappa         = %10g (m^3/J)\n", kappa);
        fprintf(fp, "Adiabatic bulk modulus                   = %10g (J/m^3)\n", 1.0 / kappa);
    }
    if (cp != NOTSET)
    {
        fprintf(fp, "Heat capacity at constant pressure Cp    = %10g J/(mol K)\n", cp);
    }
    if (cv != NOTSET)
    {
        fprintf(fp, "Heat capacity at constant volume Cv      = %10g J/(mol K)\n", cv);
    }
    if (dcp != NOTSET)
    {
        fprintf(fp, "Cp-Cv                                    =  %10g J/(mol K)\n", dcp);
    }
    please_cite(fp, "Allen1987a");
}
else
{
    fprintf(fp, "You should select the temperature in order to obtain fluctuation properties.\n");
}
}

static void analyse_ener(gmx_bool                         bCorr,
                         const char*                      corrfn,
                         const char*                      eviscofn,
                         const char*                      eviscoifn,
                         gmx_bool                         bFee,
                         gmx_bool                         bSum,
                         gmx_bool                         bFluct,
                         const bool                       computeACViscosity,
                         const bool                       computeEinsteinViscosity,
                         const int                        einsteinRestarts,
                         const char*                      visfn,
                         int                              nmol,
                         int64_t                          start_step,
                         double                           start_t,
                         int64_t                          step,
                         double                           t,
                         real                             reftemp,
                         enerdata_t*                      edat,
                         int                              nset,
                         const int                        set[],
                         const gmx_bool*                  bIsEner,
                         gmx::ArrayRef<const std::string> leg,
                         gmx_enxnm_t*                     enm,
                         real                             Vaver,
                         real                             ezero,
                         int                              nbmin,
                         int                              nbmax,
                         const gmx_output_env_t*          oenv)
{
    FILE* fp;
    /* Check out the printed manual for equations! */
    double      Dt, aver, stddev, errest, delta_t, totaldrift;
    enerdata_t* esum = nullptr;
    real        integral, intBulk, Temp = 0, Pres = 0;
    real        pr_aver, pr_stddev, pr_errest;
    double      beta = 0, expE, expEtot, *fee = nullptr;
    int64_t     nsteps;
    int         nexact, nnotexact;
    int         i, j, nout;
    char        buf[256], eebuf[100];

    nsteps = step - start_step + 1;
    if (nsteps < 1)
    {
        fprintf(stdout, "Not enough steps (%s) for statistics\n", gmx_step_str(nsteps, buf));
    }
    else
    {
        /* Calculate the time difference */
        delta_t = t - start_t;

        fprintf(stdout,
                "\nStatistics over %s steps [ %.4f through %.4f ps ], %d data sets\n",
                gmx_step_str(nsteps, buf),
                start_t,
                t,
                nset);

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
            fprintf(stdout, "All statistics are over %s points\n", gmx_step_str(edat->npoints, buf));
        }
        else if (nexact == 0 || edat->npoints == edat->nframes)
        {
            fprintf(stdout, "All statistics are over %d points (frames)\n", edat->nframes);
        }
        else
        {
            fprintf(stdout, "The term%s", nnotexact == 1 ? "" : "s");
            for (i = 0; (i < nset); i++)
            {
                if (!edat->s[i].bExactStat)
                {
                    fprintf(stdout, " '%s'", leg[i].c_str());
                }
            }
            fprintf(stdout,
                    " %s has statistics over %d points (frames)\n",
                    nnotexact == 1 ? "is" : "are",
                    edat->nframes);
            fprintf(stdout, "All other statistics are over %s points\n", gmx_step_str(edat->npoints, buf));
        }
        fprintf(stdout, "\n");

        fprintf(stdout,
                "%-24s %10s %10s %10s %10s",
                "Energy",
                "Average",
                "Err.Est.",
                "RMSD",
                "Tot-Drift");
        if (bFee)
        {
            fprintf(stdout, "  %10s\n", "-kT ln<e^(E/kT)>");
        }
        else
        {
            fprintf(stdout, "\n");
        }
        fprintf(stdout,
                "-------------------------------------------------------------------------------"
                "\n");

        /* Initiate locals, only used with -sum */
        expEtot = 0;
        if (bFee)
        {
            beta = 1.0 / (gmx::c_boltz * reftemp);
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
                    expE += std::exp(beta * (edat->s[i].ener[j] - aver) / nmol);
                }
                if (bSum)
                {
                    expEtot += expE / edat->nframes;
                }

                fee[i] = std::log(expE / edat->nframes) / beta + aver / nmol;
            }
            if (std::strstr(leg[i].c_str(), "empera") != nullptr)
            {
                Temp = aver;
            }
            else if (std::strstr(leg[i].c_str(), "olum") != nullptr)
            {
                Vaver = aver;
            }
            else if (std::strstr(leg[i].c_str(), "essure") != nullptr)
            {
                Pres = aver;
            }
            if (bIsEner[i])
            {
                pr_aver   = aver / nmol - ezero;
                pr_stddev = stddev / nmol;
                pr_errest = errest / nmol;
            }
            else
            {
                pr_aver   = aver;
                pr_stddev = stddev;
                pr_errest = errest;
            }

            /* Multiply the slope in steps with the number of steps taken */
            totaldrift = (edat->nsteps - 1) * edat->s[i].slope;
            if (bIsEner[i])
            {
                totaldrift /= nmol;
            }

            ee_pr(pr_errest, sizeof(eebuf), eebuf);
            fprintf(stdout, "%-24s %10g %10s %10g %10g", leg[i].c_str(), pr_aver, eebuf, pr_stddev, totaldrift);
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
            totaldrift = (edat->nsteps - 1) * esum->s[0].slope;
            ee_pr(esum->s[0].ee / nmol, sizeof(eebuf), eebuf);
            fprintf(stdout,
                    "%-24s %10g %10s %10s %10g  (%s)",
                    "Total",
                    esum->s[0].av / nmol,
                    eebuf,
                    "--",
                    totaldrift / nmol,
                    enm[set[0]].unit);
            /* pr_aver,pr_stddev,a,totaldrift */
            if (bFee)
            {
                fprintf(stdout,
                        "  %10g  %10g\n",
                        std::log(expEtot) / beta + esum->s[0].av / nmol,
                        std::log(expEtot) / beta);
            }
            else
            {
                fprintf(stdout, "\n");
            }
        }

        /* Do correlation function */
        if (edat->nframes > 1)
        {
            Dt = delta_t / (edat->nframes - 1);
        }
        else
        {
            Dt = 0;
        }
        if (computeACViscosity || computeEinsteinViscosity)
        {
            std::array<std::string, 2> localLeg = { "Shear", "Bulk" };
            real                       factor;
            real**                     eneset;

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
                eneset[0][i] = 0.5 * (edat->s[1].ener[i] + edat->s[3].ener[i]);
                eneset[1][i] = 0.5 * (edat->s[2].ener[i] + edat->s[6].ener[i]);
                eneset[2][i] = 0.5 * (edat->s[5].ener[i] + edat->s[7].ener[i]);
                for (j = 3; j <= 11; j++)
                {
                    eneset[j][i] = edat->s[j].ener[i];
                }
                eneset[11][i] -= Pres;
            }

            if (computeEinsteinViscosity)
            {
                einstein_visco(eviscofn, eviscoifn, 3, *edat, Vaver, Temp, einsteinRestarts, Dt, oenv);
            }

            if (computeACViscosity)
            {
                /*do_autocorr(corrfn,buf,nenergy,3,eneset,Dt,eacNormal,TRUE);*/
                /* Do it for shear viscosity */
                std::strcpy(buf, "Shear Viscosity");
                low_do_autocorr(corrfn,
                                oenv,
                                buf,
                                edat->nframes,
                                3,
                                (edat->nframes + 1) / 2,
                                eneset,
                                Dt,
                                eacNormal,
                                1,
                                TRUE,
                                FALSE,
                                FALSE,
                                0.0,
                                0.0,
                                0);

                /* Now for bulk viscosity */
                std::strcpy(buf, "Bulk Viscosity");
                low_do_autocorr(corrfn,
                                oenv,
                                buf,
                                edat->nframes,
                                1,
                                (edat->nframes + 1) / 2,
                                &(eneset[11]),
                                Dt,
                                eacNormal,
                                1,
                                TRUE,
                                FALSE,
                                FALSE,
                                0.0,
                                0.0,
                                0);

                factor = (Vaver * 1e-26 / (gmx::c_boltzmann * Temp)) * Dt;
                fp     = xvgropen(visfn, buf, "Time (ps)", "\\8h\\4 (cp)", oenv);
                xvgrLegend(fp, localLeg, oenv);

                /* Use trapezium rule for integration */
                integral = 0;
                intBulk  = 0;
                nout     = get_acfnout();
                if ((nout < 2) || (nout >= edat->nframes / 2))
                {
                    nout = edat->nframes / 2;
                }
                for (i = 1; (i < nout); i++)
                {
                    integral += 0.5 * (eneset[0][i - 1] + eneset[0][i]) * factor;
                    intBulk += 0.5 * (eneset[11][i - 1] + eneset[11][i]) * factor;
                    fprintf(fp, "%10g  %10g  %10g\n", (i * Dt), integral, intBulk);
                }
                xvgrclose(fp);
            }

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
                std::strcpy(buf, "Autocorrelation of Energy Fluctuations");
            }
            else
            {
                std::strcpy(buf, "Energy Autocorrelation");
            }
#if 0
            do_autocorr(corrfn, oenv, buf, edat->nframes,
                        bSum ? 1                 : nset,
                        bSum ? &edat->s[nset-1].energyGroupPairTerms : eneset,
                        (delta_t/edat->nframes), eacNormal, FALSE);
#endif
        }
    }
}

static void print_time(FILE* fp, double t)
{
    fprintf(fp, "%12.6f", t);
}

static void print1(FILE* fp, gmx_bool bDp, real e)
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

static void fec(const char*                      ene2fn,
                const char*                      runavgfn,
                real                             reftemp,
                int                              nset,
                const int                        set[],
                gmx::ArrayRef<const std::string> leg,
                enerdata_t*                      edat,
                double                           time[],
                const gmx_output_env_t*          oenv)
{
    std::array<std::string, 2> ravgleg = { "\\8D\\4E = E\\sB\\N-E\\sA\\N",
                                           "<e\\S-\\8D\\4E/kT\\N>\\s0..t\\N" };
    FILE*                      fp;
    ener_file_t                enx;
    int                        timecheck, nenergy, nenergy2, maxenergy;
    int                        i, j;
    gmx_bool                   bCont;
    real                       aver, beta;
    real**                     eneset2;
    double                     dE, sum;
    gmx_enxnm_t*               enm = nullptr;
    t_enxframe*                fr;
    char                       buf[22];

    /* read second energy file */
    snew(fr, 1);
    enm = nullptr;
    enx = open_enx(ene2fn, "r");
    do_enxnms(enx, &(fr->nre), &enm);

    snew(eneset2, nset + 1);
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

        } while (bCont && (timecheck < 0));

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
                GMX_RELEASE_ASSERT(time != nullptr, "trying to dereference NULL time pointer");

                if (fr->t != time[nenergy2])
                {
                    fprintf(stderr,
                            "\nWARNING time mismatch %g!=%g at frame %s\n",
                            fr->t,
                            time[nenergy2],
                            gmx_step_str(fr->step, buf));
                }
                for (i = 0; i < nset; i++)
                {
                    eneset2[i][nenergy2] = fr->ener[set[i]].e;
                }
                nenergy2++;
            }
        }
    } while (bCont && (timecheck == 0));

    /* check */
    if (edat->nframes != nenergy2)
    {
        fprintf(stderr, "\nWARNING file length mismatch %d!=%d\n", edat->nframes, nenergy2);
    }
    nenergy = std::min(edat->nframes, nenergy2);

    /* calculate fe difference dF = -kT ln < exp(-(E_B-E_A)/kT) >_A */
    fp = nullptr;
    if (runavgfn)
    {
        fp = xvgropen(runavgfn,
                      "Running average free energy difference",
                      "Time (" unit_time ")",
                      "\\8D\\4E (" unit_energy ")",
                      oenv);
        xvgrLegend(fp, ravgleg, oenv);
    }
    fprintf(stdout, "\n%-24s %10s\n", "Energy", "dF = -kT ln < exp(-(EB-EA)/kT) >A");
    sum  = 0;
    beta = 1.0 / (gmx::c_boltz * reftemp);
    for (i = 0; i < nset; i++)
    {
        if (gmx_strcasecmp(leg[i].c_str(), enm[set[i]].name) != 0)
        {
            fprintf(stderr,
                    "\nWARNING energy set name mismatch %s!=%s\n",
                    leg[i].c_str(),
                    enm[set[i]].name);
        }
        for (j = 0; j < nenergy; j++)
        {
            dE = eneset2[i][j] - edat->s[i].ener[j];
            sum += std::exp(-dE * beta);
            if (fp)
            {
                fprintf(fp, "%10g %10g %10g\n", time[j], dE, -gmx::c_boltz * reftemp * std::log(sum / (j + 1)));
            }
        }
        aver = -gmx::c_boltz * reftemp * std::log(sum / nenergy);
        fprintf(stdout, "%-24s %10g\n", leg[i].c_str(), aver);
    }
    if (fp)
    {
        xvgrclose(fp);
    }
    sfree(fr);
}


static void do_dhdl(t_enxframe*             fr,
                    const t_inputrec*       ir,
                    FILE**                  fp_dhdl,
                    const char*             filename,
                    gmx_bool                bDp,
                    int*                    blocks,
                    int*                    hists,
                    int*                    samples,
                    int*                    nlambdas,
                    const gmx_output_env_t* oenv)
{
    const char *dhdl = "dH/d\\lambda", *deltag = "\\DeltaH", *lambda = "\\lambda";
    char        title[STRLEN], label_x[STRLEN], label_y[STRLEN];
    char        buf[STRLEN];
    int         nblock_hist = 0, nblock_dh = 0;
    int         i, j, k;
    /* coll data */
    double       temp = 0, start_time = 0, delta_time = 0, start_lambda = 0;
    static int   setnr             = 0;
    double*      native_lambda_vec = nullptr;
    const char** lambda_components = nullptr;
    int          n_lambda_vec      = 0;
    bool         firstPass         = true;

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
            if ((fr->block[i].nsub < 1) || (fr->block[i].sub[0].type != XdrDataType::Double)
                || (fr->block[i].sub[0].nr < 5))
            {
                gmx_fatal(FARGS, "Unexpected block data");
            }

            /* read the data from the DHCOLL block */
            temp         = fr->block[i].sub[0].dval[0];
            start_time   = fr->block[i].sub[0].dval[1];
            delta_time   = fr->block[i].sub[0].dval[2];
            start_lambda = fr->block[i].sub[0].dval[3];
            if (fr->block[i].nsub > 1)
            {
                if (firstPass)
                {
                    n_lambda_vec = fr->block[i].sub[1].ival[1];
                    snew(lambda_components, n_lambda_vec);
                    snew(native_lambda_vec, n_lambda_vec);
                    firstPass = false;
                }
                else
                {
                    if (n_lambda_vec != fr->block[i].sub[1].ival[1])
                    {
                        gmx_fatal(FARGS, "Unexpected change of basis set in lambda");
                    }
                }
                for (j = 0; j < n_lambda_vec; j++)
                {
                    native_lambda_vec[j] = fr->block[i].sub[0].dval[5 + j];
                    lambda_components[j] =
                            enumValueToStringSingular(static_cast<FreeEnergyPerturbationCouplingType>(
                                    fr->block[i].sub[1].ival[2 + j]));
                }
            }
        }
    }
    // Clean up!
    sfree(native_lambda_vec);
    sfree(lambda_components);
    if (nblock_hist == 0 && nblock_dh == 0)
    {
        /* don't do anything */
        return;
    }
    if (nblock_hist > 0 && nblock_dh > 0)
    {
        gmx_fatal(FARGS,
                  "This energy file contains both histogram dhdl data and non-histogram dhdl data. "
                  "Don't know what to do.");
    }
    if (!*fp_dhdl)
    {
        if (nblock_dh > 0)
        {
            /* we have standard, non-histogram data --
               call open_dhdl to open the file */
            /* TODO this is an ugly hack that needs to be fixed: this will only
               work if the order of data is always the same and if we're
               only using the gmx energy compiled with the mdrun that produced
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

    (*hists) += nblock_hist;
    (*blocks) += nblock_dh;
    (*nlambdas) = nblock_hist + nblock_dh;

    /* write the data */
    if (nblock_hist > 0)
    {
        int64_t sum = 0;
        /* histograms */
        for (i = 0; i < fr->nblock; i++)
        {
            t_enxblock* blk = &(fr->block[i]);
            if (blk->id == enxDHHIST)
            {
                double  foreign_lambda, dx;
                int64_t x0;
                int     nhist, derivative;

                /* check the block types etc. */
                if ((blk->nsub < 2) || (blk->sub[0].type != XdrDataType::Double)
                    || (blk->sub[1].type != XdrDataType::Int64) || (blk->sub[0].nr < 2)
                    || (blk->sub[1].nr < 2))
                {
                    gmx_fatal(FARGS, "Unexpected block data in file");
                }
                foreign_lambda = blk->sub[0].dval[0];
                dx             = blk->sub[0].dval[1];
                nhist          = blk->sub[1].lval[0];
                derivative     = blk->sub[1].lval[1];
                for (j = 0; j < nhist; j++)
                {
                    const std::string legend =
                            derivative ? gmx::formatString("N(%s | %s=%g)", dhdl, lambda, start_lambda)
                                       : gmx::formatString("N(%s(%s=%g) | %s=%g)",
                                                           deltag,
                                                           lambda,
                                                           foreign_lambda,
                                                           lambda,
                                                           start_lambda);
                    x0 = blk->sub[1].lval[2 + j];

                    xvgrNewDataset(*fp_dhdl, setnr, gmx::arrayRefFromArray(&legend, 1), oenv);
                    setnr++;
                    for (k = 0; k < blk->sub[j + 2].nr; k++)
                    {
                        int    hist;
                        double xmin, xmax;

                        hist = blk->sub[j + 2].ival[k];
                        xmin = (x0 + k) * dx;
                        xmax = (x0 + k + 1) * dx;
                        fprintf(*fp_dhdl, "%g %d\n%g %d\n", xmin, hist, xmax, hist);
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
        (*samples) += static_cast<int>(sum / nblock_hist);
    }
    else
    {
        /* raw dh */
        int len = 0;

        for (i = 0; i < fr->nblock; i++)
        {
            t_enxblock* blk = &(fr->block[i]);
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
            double time = start_time + delta_time * i;

            fprintf(*fp_dhdl, "%.4f ", time);

            for (j = 0; j < fr->nblock; j++)
            {
                t_enxblock* blk = &(fr->block[j]);
                if (blk->id == enxDH)
                {
                    double value;
                    if (blk->sub[2].type == XdrDataType::Float)
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
                        fprintf(*fp_dhdl, "%4d", static_cast<int>(value)); /* if expanded ensembles and zero, this is a state value, it's an integer. We need a cleaner conditional than if j==1! */
                    }
                    else
                    {
                        if (bDp)
                        {
                            fprintf(*fp_dhdl, " %#.12g", value); /* print normal precision */
                        }
                        else
                        {
                            fprintf(*fp_dhdl, " %#.8g", value); /* print normal precision */
                        }
                    }
                }
            }
            fprintf(*fp_dhdl, "\n");
        }
    }
}


int gmx_energy(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] extracts energy components",
        "from an energy file. The user is prompted to interactively",
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
        "for quantum effects. Use the [gmx-dos] program if you need that (and you do).[PAR]",

        "Option [TT]-odh[tt] extracts and plots the free energy data",
        "(Hamiltoian differences and/or the Hamiltonian derivative dhdl)",
        "from the [TT]ener.edr[tt] file.[PAR]",

        "With [TT]-fee[tt] an estimate is calculated for the free-energy",
        "difference with an ideal gas state::",
        "",
        "  [GRK]Delta[grk] A = A(N,V,T) - A[SUB]idealgas[sub](N,V,T) = kT ",
        "  [LN][CHEVRON][EXP]U[SUB]pot[sub]/kT[exp][chevron][ln]",
        "  [GRK]Delta[grk] G = G(N,p,T) - G[SUB]idealgas[sub](N,p,T) = kT ",
        "  [LN][CHEVRON][EXP]U[SUB]pot[sub]/kT[exp][chevron][ln]",
        "",
        "where k is Boltzmann's constant, T is set by [TT]-fetemp[tt] and",
        "the average is over the ensemble (or time in a trajectory).",
        "Note that this is in principle",
        "only correct when averaging over the whole (Boltzmann) ensemble",
        "and using the potential energy. This also allows for an entropy",
        "estimate using::",
        "",
        "  [GRK]Delta[grk] S(N,V,T) = S(N,V,T) - S[SUB]idealgas[sub](N,V,T) = ",
        "  ([CHEVRON]U[SUB]pot[sub][chevron] - [GRK]Delta[grk] A)/T",
        "  [GRK]Delta[grk] S(N,p,T) = S(N,p,T) - S[SUB]idealgas[sub](N,p,T) = ",
        "  ([CHEVRON]U[SUB]pot[sub][chevron] + pV - [GRK]Delta[grk] G)/T",
        "",

        "When a second energy file is specified ([TT]-f2[tt]), a free energy",
        "difference is calculated::",
        "",
        "  dF = -kT ",
        "  [LN][CHEVRON][EXP]-(E[SUB]B[sub]-E[SUB]A[sub]) / ",
        "  kT[exp][chevron][SUB]A[sub][ln],",
        "",
        "where E[SUB]A[sub] and E[SUB]B[sub] are the energies from the first and second energy",
        "files, and the average is over the ensemble A. The running average",
        "of the free energy difference is printed to a file specified by [TT]-ravg[tt].",
        "[BB]Note[bb] that the energies must both be calculated from the same trajectory.[PAR]",

        "For liquids, viscosities can be calculated by integrating the auto-correlation function ",
        "of, or by using the Einstein formula for, the off-diagonal pressure elements. ",
        "The option [TT]-vis[tt] turns calculation of the shear and bulk viscosity through ",
        "integration of the auto-correlation function. For accurate results, this requires ",
        "extremely frequent computation and output of the pressure tensor. ",
        "The Einstein formula does not require frequent output and is therefore more convenient. ",
        "Note that frequent pressure calculation (nstcalcenergy mdp parameter) is still needed. ",
        "Option [TT]-evicso[tt] gives this shear viscosity estimate and option [TT]-eviscoi[tt] ",
        "the integral. Using one of these two options also triggers the other. ",
        "The viscosity is computed from integrals averaged over [TT]-einstein_restarts[tt] ",
        "starting points uniformly distributed over the first quarter of the trajectory."
    };
    static gmx_bool bSum = FALSE, bFee = FALSE, bPrAll = FALSE, bFluct = FALSE, bDriftCorr = FALSE;
    static gmx_bool bDp = FALSE, bMutot = FALSE, bOrinst = FALSE, bOvec = FALSE, bFluctProps = FALSE;
    static int      nmol = 1, nbmin = 5, nbmax = 5;
    static real     reftemp = 300.0, ezero = 0;
    static int      einsteinRestarts = 100;
    t_pargs         pa[]             = {
        { "-fee", FALSE, etBOOL, { &bFee }, "Do a free energy estimate" },
        { "-fetemp",
          FALSE,
          etREAL,
          { &reftemp },
          "Reference temperature for free energy calculation" },
        { "-zero", FALSE, etREAL, { &ezero }, "Subtract a zero-point energy" },
        { "-sum",
          FALSE,
          etBOOL,
          { &bSum },
          "Sum the energy terms selected rather than display them all" },
        { "-dp", FALSE, etBOOL, { &bDp }, "Print energies in high precision" },
        { "-nbmin", FALSE, etINT, { &nbmin }, "Minimum number of blocks for error estimate" },
        { "-nbmax", FALSE, etINT, { &nbmax }, "Maximum number of blocks for error estimate" },
        { "-mutot",
          FALSE,
          etBOOL,
          { &bMutot },
          "Compute the total dipole moment from the components" },
        { "-aver",
          FALSE,
          etBOOL,
          { &bPrAll },
          "Also print the exact average and rmsd stored in the energy frames (only when 1 term is "
          "requested)" },
        { "-nmol",
          FALSE,
          etINT,
          { &nmol },
          "Number of molecules in your sample: the energies are divided by this number" },
        { "-fluct_props",
          FALSE,
          etBOOL,
          { &bFluctProps },
          "Compute properties based on energy fluctuations, like heat capacity" },
        { "-driftcorr",
          FALSE,
          etBOOL,
          { &bDriftCorr },
          "Useful only for calculations of fluctuation properties. The drift in the observables "
          "will be subtracted before computing the fluctuation properties." },
        { "-fluc",
          FALSE,
          etBOOL,
          { &bFluct },
          "Calculate autocorrelation of energy fluctuations rather than energy itself" },
        { "-orinst", FALSE, etBOOL, { &bOrinst }, "Analyse instantaneous orientation data" },
        { "-ovec", FALSE, etBOOL, { &bOvec }, "Also plot the eigenvectors with [TT]-oten[tt]" },
        { "-einstein_restarts",
          FALSE,
          etINT,
          { &einsteinRestarts },
          "Number of restarts for computing the viscosity using the Einstein relation" }
    };
    static const char* setnm[] = { "Pres-XX", "Pres-XY",     "Pres-XZ", "Pres-YX",
                                   "Pres-YY", "Pres-YZ",     "Pres-ZX", "Pres-ZY",
                                   "Pres-ZZ", "Temperature", "Volume",  "Pressure" };

    FILE*        out     = nullptr;
    FILE*        fp_dhdl = nullptr;
    ener_file_t  fp;
    int          timecheck = 0;
    enerdata_t   edat;
    gmx_enxnm_t* enm = nullptr;
    t_enxframe * frame, *fr = nullptr;
    int          cur = 0;
#define NEXT (1 - cur)
    int                      nre, nfr;
    int64_t                  start_step;
    real                     start_t;
    gmx_bool                 bDHDL;
    gmx_bool                 bFoundStart, bCont;
    double                   sum, dbl;
    double*                  time = nullptr;
    real                     Vaver;
    int *                    set     = nullptr, i, j, nset, sss;
    gmx_bool*                bIsEner = nullptr;
    std::vector<std::string> leg;
    char                     buf[256];
    gmx_output_env_t*        oenv;
    int                      dh_blocks = 0, dh_hists = 0, dh_samples = 0, dh_lambdas = 0;

    t_filenm fnm[] = {
        { efEDR, "-f", nullptr, ffREAD },        { efEDR, "-f2", nullptr, ffOPTRD },
        { efTPR, "-s", nullptr, ffOPTRD },       { efXVG, "-o", "energy", ffWRITE },
        { efXVG, "-viol", "violaver", ffOPTWR }, { efXVG, "-pairs", "pairs", ffOPTWR },
        { efXVG, "-corr", "enecorr", ffOPTWR },  { efXVG, "-vis", "visco", ffOPTWR },
        { efXVG, "-evisco", "evisco", ffOPTWR }, { efXVG, "-eviscoi", "eviscoi", ffOPTWR },
        { efXVG, "-ravg", "runavgdf", ffOPTWR }, { efXVG, "-odh", "dhdl", ffOPTWR }
    };
#define NFILE asize(fnm)
    int      npargs;
    t_pargs* ppa;

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);
    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END, NFILE, fnm, npargs, ppa, asize(desc), desc, 0, nullptr, &oenv))
    {
        sfree(ppa);
        return 0;
    }

    bDHDL = opt2bSet("-odh", NFILE, fnm);

    nset = 0;

    snew(frame, 2);
    fp = open_enx(ftp2fn(efEDR, NFILE, fnm), "r");
    do_enxnms(fp, &nre, &enm);

    Vaver = -1;

    const bool computeACViscosity = opt2bSet("-vis", NFILE, fnm);
    // Assign 'true' if either -evisco or -eviscoi flags are specified
    const bool computeEinsteinViscosity =
            opt2bSet("-evisco", NFILE, fnm) || opt2bSet("-eviscoi", NFILE, fnm);

    t_inputrec  irInstance;
    t_inputrec* ir = &irInstance;

    if (!bDHDL)
    {
        if (computeACViscosity || computeEinsteinViscosity)
        {
            nset = asize(setnm);
            snew(set, nset);
            /* This is nasty code... To extract Pres tensor, Volume and Temperature */
            for (j = 0; j < nset; j++)
            {
                for (i = 0; i < nre; i++)
                {
                    if (std::strstr(enm[i].name, setnm[j]))
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
                        gmx_fatal(FARGS, "Could not find term %s for viscosity calculation", setnm[j]);
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
                if (std::strcmp(enm[set[i]].unit, enm[set[j]].unit) == 0)
                {
                    break;
                }
            }
            if (j == i)
            {
                std::strcat(buf, ", (");
                std::strcat(buf, enm[set[i]].unit);
                std::strcat(buf, ")");
            }
        }
        out = xvgropen(opt2fn("-o", NFILE, fnm), "GROMACS Energies", "Time (ps)", buf, oenv);

        for (i = 0; (i < nset); i++)
        {
            leg.emplace_back(enm[set[i]].name);
        }
        if (bSum)
        {
            leg.emplace_back("Sum");
        }
        xvgrLegend(out, leg, oenv);

        snew(bIsEner, nset);
        for (i = 0; (i < nset); i++)
        {
            bIsEner[i] = FALSE;
            for (j = 0; (j <= F_ETOT); j++)
            {
                bIsEner[i] = bIsEner[i]
                             || (gmx_strcasecmp(interaction_function[j].longname, leg[i].c_str()) == 0);
            }
            bIsEner[i] = bIsEner[i] || gmx::equalCaseInsensitive(pvEnergyFieldName, leg[i]);
            bIsEner[i] = bIsEner[i] || gmx::equalCaseInsensitive(enthalpyEnergyFieldName, leg[i]);
            for (const char* name : virialEnergyFieldNames)
            {
                bIsEner[i] = bIsEner[i] || gmx::equalCaseInsensitive(name, leg[i]);
            }
        }
        if (bPrAll && nset > 1)
        {
            gmx_fatal(FARGS, "Printing averages can only be done when a single set is selected");
        }
    }
    else if (bDHDL)
    {
        get_dhdl_parms(ftp2fn(efTPR, NFILE, fnm), ir);
    }

    /* Initiate energies and set them to zero */
    edat.nsteps    = 0;
    edat.npoints   = 0;
    edat.nframes   = 0;
    edat.step      = nullptr;
    edat.steps     = nullptr;
    edat.points    = nullptr;
    edat.bHaveSums = TRUE;
    snew(edat.s, nset);

    /* Initiate counters */
    bFoundStart = FALSE;
    start_step  = 0;
    start_t     = 0;
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
        } while (bCont && (timecheck < 0));

        if ((timecheck == 0) && bCont)
        {
            /* We read a valid frame, so we can use it */
            fr = &(frame[NEXT]);

            if (fr->nre > 0)
            {
                /* The frame contains energies, so update cur */
                cur = NEXT;

                if (edat.nframes % 1000 == 0)
                {
                    srenew(edat.step, edat.nframes + 1000);
                    std::memset(&(edat.step[edat.nframes]), 0, 1000 * sizeof(edat.step[0]));
                    srenew(edat.steps, edat.nframes + 1000);
                    std::memset(&(edat.steps[edat.nframes]), 0, 1000 * sizeof(edat.steps[0]));
                    srenew(edat.points, edat.nframes + 1000);
                    std::memset(&(edat.points[edat.nframes]), 0, 1000 * sizeof(edat.points[0]));

                    for (i = 0; i < nset; i++)
                    {
                        srenew(edat.s[i].ener, edat.nframes + 1000);
                        std::memset(&(edat.s[i].ener[edat.nframes]), 0, 1000 * sizeof(edat.s[i].ener[0]));
                        srenew(edat.s[i].es, edat.nframes + 1000);
                        std::memset(&(edat.s[i].es[edat.nframes]), 0, 1000 * sizeof(edat.s[i].es[0]));
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
                        edat.npoints += 1;
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
             * Store energies for analysis afterwards...
             */
            if (!bDHDL && (fr->nre > 0))
            {
                if (edat.nframes % 1000 == 0)
                {
                    srenew(time, edat.nframes + 1000);
                }
                time[edat.nframes] = fr->t;
                edat.nframes++;
            }
            if (bDHDL)
            {
                do_dhdl(fr, ir, &fp_dhdl, opt2fn("-odh", NFILE, fnm), bDp, &dh_blocks, &dh_hists, &dh_samples, &dh_lambdas, oenv);
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
                            print1(out, bDp, fr->ener[set[0]].esum / fr->nsum);
                            print1(out, bDp, std::sqrt(fr->ener[set[0]].eav / fr->nsum));
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
                            print1(out, bDp, sum / nmol - ezero);
                        }
                        else
                        {
                            for (i = 0; (i < nset); i++)
                            {
                                if (bIsEner[i])
                                {
                                    print1(out, bDp, (fr->ener[set[i]].e) / nmol - ezero);
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
            }
        }
    } while (bCont && (timecheck == 0));

    fprintf(stderr, "\n");
    done_ener_file(fp);
    if (out)
    {
        xvgrclose(out);
    }

    if (bDHDL)
    {
        if (fp_dhdl)
        {
            gmx_fio_fclose(fp_dhdl);
            printf("\n\nWrote %d lambda values with %d samples as ", dh_lambdas, dh_samples);
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
        double dt = (frame[cur].t - start_t) / (edat.nframes - 1);
        analyse_ener(opt2bSet("-corr", NFILE, fnm),
                     opt2fn("-corr", NFILE, fnm),
                     opt2fn("-evisco", NFILE, fnm),
                     opt2fn("-eviscoi", NFILE, fnm),
                     bFee,
                     bSum,
                     bFluct,
                     computeACViscosity,
                     computeEinsteinViscosity,
                     einsteinRestarts,
                     opt2fn("-vis", NFILE, fnm),
                     nmol,
                     start_step,
                     start_t,
                     frame[cur].step,
                     frame[cur].t,
                     reftemp,
                     &edat,
                     nset,
                     set,
                     bIsEner,
                     leg,
                     enm,
                     Vaver,
                     ezero,
                     nbmin,
                     nbmax,
                     oenv);
        if (bFluctProps)
        {
            calc_fluctuation_props(stdout, bDriftCorr, dt, nset, nmol, leg, &edat, nbmin, nbmax);
        }
    }
    if (opt2bSet("-f2", NFILE, fnm))
    {
        fec(opt2fn("-f2", NFILE, fnm), opt2fn("-ravg", NFILE, fnm), reftemp, nset, set, leg, &edat, time, oenv);
    }
    // Clean up!
    done_enerdata_t(nset, &edat);
    sfree(time);
    free_enxframe(&frame[0]);
    free_enxframe(&frame[1]);
    sfree(frame);
    free_enxnms(nre, enm);
    sfree(ppa);
    sfree(set);
    sfree(bIsEner);
    {
        const char* nxy = "-nxy";

        do_view(oenv, opt2fn("-o", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-ravg", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-odh", NFILE, fnm), nxy);
    }
    output_env_done(oenv);

    return 0;
}

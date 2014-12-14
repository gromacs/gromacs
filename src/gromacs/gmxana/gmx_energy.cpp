/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>

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
#include "gromacs/mdlib/mdebin.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"
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
        return std::pow(x, y);
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

    if ((getenv("GMX_ENER_VERBOSE")) != nullptr)
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
    int len = std::strlen(buf);

    while ((len > 0) && (buf[len-1] == '\n'))
    {
        buf[len-1] = '\0';
        len--;
    }
}

static int *select_by_name(int nre, gmx_enxnm_t *nm, int *nset)
{
    gmx_bool   *bE;
    int         k, kk, j, i, nmatch, nind, nss;
    int        *set;
    gmx_bool    bEOF, bVerbose = TRUE, bLong = FALSE;
    char       *ptr, buf[STRLEN];
    const char *fm4   = "%3d  %-14s";
    const char *fm2   = "%3d  %-34s";
    char      **newnm = nullptr;

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
                for (kk = k; kk < k+4; kk++)
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
        bEOF = (std::strlen(buf) == 0);
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
            }
            while (!bEOF && (ptr && (std::strlen(ptr) > 0)));
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

static void get_orires_parms(const char *topnm, t_inputrec *ir,
                             int *nor, int *nex, int **label, real **obs)
{
    gmx_mtop_t      mtop;
    gmx_localtop_t *top;
    t_iparams      *ip;
    int             natoms, i;
    t_iatom        *iatom;
    int             nb;
    matrix          box;

    read_tpx(topnm, ir, box, &natoms, nullptr, nullptr, &mtop);
    top = gmx_mtop_generate_local_top(&mtop, ir->efep != efepNO);

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

    read_tpx(topnm, ir, box, &natoms, nullptr, nullptr, mtop);
    snew(*ltop, 1);
    top   = gmx_mtop_generate_local_top(mtop, ir->efep != efepNO);
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
            rav     += gmx::square(rav3[j]);
            rsum    += mypow(rt[j], -6);
        }
        rsum    = std::max(0.0, mypow(rsum, -sixth)-bounds[i]);
        rav     = std::max(0.0, mypow(rav, -sixth)-bounds[i]);

        sumt    += rsum;
        sumaver += rav;
    }
    *st = sumt;
    *sa = sumaver;
}

static void analyse_disre(const char *voutfn,    int nframes,
                          real violaver[], real bounds[], int index[],
                          int pair[],      int nbounds,
                          const gmx_output_env_t *oenv)
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
            sumaver += gmx::square(violaver[j]/nframes);
        }
        sumaver = std::max(0.0, mypow(sumaver, minsixth)-bounds[i]);

        sumt   += sumaver;
        sum     = std::max(sum, sumaver);
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
    return (ees->sav2/nb - gmx::square(ees->sav/nb))/(nb - 1);
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

int gmx_nmr(int argc, char *argv[])
{
    const char        *desc[] = {
        "[THISMODULE] extracts distance restraint",
        "data from an energy file. The user is prompted to interactively",
        "select the desired terms.[PAR]",

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


    };
    static gmx_bool    bSum    = FALSE, bPrAll = FALSE;
    static gmx_bool    bDp     = FALSE, bOrinst = FALSE, bOvec = FALSE;
    static int         skip    = 0;
    static real        ezero   = 0;
    t_pargs            pa[]    = {
        { "-dp",   FALSE, etBOOL, {&bDp},
          "Print energies in high precision" },
        { "-skip", FALSE, etINT,  {&skip},
          "Skip number of frames between data points" },
        { "-aver", FALSE, etBOOL, {&bPrAll},
          "Also print the exact average and rmsd stored in the energy frames (only when 1 term is requested)" },
        { "-orinst", FALSE, etBOOL, {&bOrinst},
          "Analyse instantaneous orientation data" },
        { "-ovec", FALSE, etBOOL, {&bOvec},
          "Also plot the eigenvectors with [TT]-oten[tt]" }
    };
    const char       * drleg[] = {
        "Running average",
        "Instantaneous"
    };

    FILE              *out     = NULL, *fp_pairs = NULL, *fort = NULL, *fodt = NULL, *foten = NULL;
    ener_file_t        fp;
    int                timecheck = 0;
    gmx_mtop_t         mtop;
    gmx_localtop_t    *top = nullptr;
    enerdata_t         edat;
    gmx_enxnm_t       *enm = nullptr;
    t_enxframe        *frame, *fr = nullptr;
    int                cur = 0;
#define NEXT (1-cur)
    int                nre, teller, teller_disre, nfr;
    gmx_int64_t        start_step;
    int                nor     = 0, nex = 0, norfr = 0, enx_i = 0;
    real              *bounds  = NULL, *violaver = NULL, *oobs = NULL, *orient = NULL, *odrms = NULL;
    int               *index   = NULL, *pair = NULL, norsel = 0, *orsel = NULL, *or_label = NULL;
    int                nbounds = 0, npairs;
    gmx_bool           bDisRe, bDRAll, bORA, bORT, bODA, bODR, bODT, bORIRE, bOTEN;
    gmx_bool           bFoundStart, bCont;
    double             sum, sumaver, sumt;
    double            *time    = NULL;
    int               *set     = NULL, i, j, k, nset, sss;
    gmx_bool          *bIsEner = NULL;
    char             **pairleg, **odtleg, **otenleg;
    char             **leg = nullptr;
    const char        *anm_j, *anm_k, *resnm_j, *resnm_k;
    int                resnr_j, resnr_k;
    const char        *orinst_sub = "@ subtitle \"instantaneous\"\n";
    char               buf[256];
    gmx_output_env_t  *oenv;
    t_enxblock        *blk       = nullptr;
    t_enxblock        *blk_disre = nullptr;
    int                ndisre    = 0;

    t_filenm           fnm[] = {
        { efEDR, "-f",    nullptr,      ffREAD  },
        { efEDR, "-f2",   nullptr,      ffOPTRD },
        { efTPR, "-s",    nullptr,      ffOPTRD },
        { efXVG, "-o",    "energy",  ffWRITE },
        { efXVG, "-viol", "violaver", ffOPTWR },
        { efXVG, "-pairs", "pairs",   ffOPTWR },
        { efXVG, "-ora",  "orienta", ffOPTWR },
        { efXVG, "-ort",  "orientt", ffOPTWR },
        { efXVG, "-oda",  "orideva", ffOPTWR },
        { efXVG, "-odr",  "oridevr", ffOPTWR },
        { efXVG, "-odt",  "oridevt", ffOPTWR },
        { efXVG, "-oten", "oriten",  ffOPTWR }
    };
#define NFILE asize(fnm)
    int                npargs;
    t_pargs           *ppa;

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);
    if (!parse_common_args(&argc, argv,
                           PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END,
                           NFILE, fnm, npargs, ppa, asize(desc), desc, 0, nullptr, &oenv))
    {
        sfree(ppa);
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

    nset = 0;

    snew(frame, 2);
    fp = open_enx(ftp2fn(efEDR, NFILE, fnm), "r");
    do_enxnms(fp, &nre, &enm);

    gmx::MDModules  mdModules;
    t_inputrec     *ir = mdModules.inputrec();

    if (!bDisRe)
    {
        set = select_by_name(nre, enm, &nset);
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
        out = xvgropen(opt2fn("-o", NFILE, fnm), "GROMACS Energies", "Time (ps)", buf,
                       oenv);

        snew(leg, nset+1);
        for (i = 0; (i < nset); i++)
        {
            leg[i] = enm[set[i]].name;
        }
        xvgr_legend(out, nset, (const char**)leg, oenv);

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

        time = nullptr;

        if (bORIRE || bOTEN)
        {
            get_orires_parms(ftp2fn(efTPR, NFILE, fnm), ir,
                             &nor, &nex, &or_label, &oobs);
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
                orsel = nullptr;
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
    else
    {
        nbounds = get_bounds(ftp2fn(efTPR, NFILE, fnm), &bounds, &index, &pair, &npairs,
                             &mtop, &top, ir);
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
                        ir->dr_tau);
            }
        }
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
    teller       = 0;
    teller_disre = 0;
    bFoundStart  = FALSE;
    start_step   = 0;
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
                    std::memset(&(edat.step[edat.nframes]), 0, 1000*sizeof(edat.step[0]));
                    srenew(edat.steps, edat.nframes+1000);
                    std::memset(&(edat.steps[edat.nframes]), 0, 1000*sizeof(edat.steps[0]));
                    srenew(edat.points, edat.nframes+1000);
                    std::memset(&(edat.points[edat.nframes]), 0, 1000*sizeof(edat.points[0]));

                    for (i = 0; i < nset; i++)
                    {
                        srenew(edat.s[i].ener, edat.nframes+1000);
                        std::memset(&(edat.s[i].ener[edat.nframes]), 0,
                                    1000*sizeof(edat.s[i].ener[0]));
                        srenew(edat.s[i].es, edat.nframes+1000);
                        std::memset(&(edat.s[i].es[edat.nframes]), 0,
                                    1000*sizeof(edat.s[i].es[0]));
                    }
                }

                nfr            = edat.nframes;
                edat.step[nfr] = fr->step;

                if (!bFoundStart)
                {
                    bFoundStart = TRUE;
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
            blk_disre = find_block_id_enxframe(fr, enxDISRE, nullptr);
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
                int molb = 0;
                for (i = 0; i < ndisre; i++)
                {
                    snew(pairleg[i], 30);
                    j = fa[3*i+1];
                    k = fa[3*i+2];
                    mtopGetAtomAndResidueName(&mtop, j, &molb, &anm_j, &resnr_j, &resnm_j, nullptr);
                    mtopGetAtomAndResidueName(&mtop, k, &molb, &anm_k, &resnr_k, &resnm_k, nullptr);
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
            if (!bDisRe && (fr->nre > 0))
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
                        GMX_RELEASE_ASSERT(blk_disre != nullptr, "Trying to dereference NULL blk_disre pointer");
 #if !GMX_DOUBLE
                        float  *disre_rt     = blk_disre->sub[0].fval;
                        float  *disre_rm3tav = blk_disre->sub[1].fval;
 #else
                        double *disre_rt     = blk_disre->sub[0].dval;
                        double *disre_rm3tav = blk_disre->sub[1].dval;
 #endif

                        print_time(out, fr->t);
                        if (violaver == nullptr)
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
                                print1(out, bDp, std::sqrt(fr->ener[set[0]].eav/fr->nsum));
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
                                print1(out, bDp, sum-ezero);
                            }
                            else
                            {
                                for (i = 0; (i < nset); i++)
                                {
                                    if (bIsEner[i])
                                    {
                                        print1(out, bDp, (fr->ener[set[i]].e)-ezero);
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
                    blk = find_block_id_enxframe(fr, enx_i, nullptr);
                    if (bORIRE && blk)
                    {
#if !GMX_DOUBLE
                        xdr_datatype dt = xdr_datatype_float;
#else
                        xdr_datatype dt = xdr_datatype_double;
#endif
                        real        *vals;

                        if ( (blk->nsub != 1) || (blk->sub[0].type != dt) )
                        {
                            gmx_fatal(FARGS, "Orientational restraints read in incorrectly");
                        }
#if !GMX_DOUBLE
                        vals = blk->sub[0].fval;
#else
                        vals = blk->sub[0].dval;
#endif

                        if (blk->sub[0].nr != nor)
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
                                odrms[i] += gmx::square(vals[i]-oobs[i]);
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
                    blk = find_block_id_enxframe(fr, enxORT, nullptr);
                    if (bOTEN && blk)
                    {
#if !GMX_DOUBLE
                        xdr_datatype dt = xdr_datatype_float;
#else
                        xdr_datatype dt = xdr_datatype_double;
#endif
                        real        *vals;

                        if ( (blk->nsub != 1) || (blk->sub[0].type != dt) )
                        {
                            gmx_fatal(FARGS, "Orientational restraints read in incorrectly");
                        }
#if !GMX_DOUBLE
                        vals = blk->sub[0].fval;
#else
                        vals = blk->sub[0].dval;
#endif

                        if (blk->sub[0].nr != nex*12)
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
            fprintf(out, "%5d  %g\n", or_label[i], std::sqrt(odrms[i]/norfr));
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

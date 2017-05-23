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
#include "gromacs/utility/strconvert.h"

static real       minthird = -1.0/3.0, minsixth = -1.0/6.0;

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

static real blk_value(t_enxblock *blk, int sub, int index)
{
    range_check(index, 0, blk->sub[sub].nr);
    if (blk->sub[sub].type == xdr_datatype_float)
    {
        return blk->sub[sub].fval[index];
    }
    else if (blk->sub[sub].type == xdr_datatype_double)
    {
        return blk->sub[sub].dval[index];
    }
    else
    {
        gmx_incons("Unknown datatype in t_enxblock");
    }
    return 0.0;
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

static void get_orires_parms(const char *topnm, t_inputrec *ir,
                             int *nor, int *nex, int **label, real **obs)
{
    gmx_mtop_t      mtop;
    t_topology      top;
    t_iparams      *ip;
    int             natoms, i;
    t_iatom        *iatom;
    int             nb;
    matrix          box;

    read_tpx(topnm, ir, box, &natoms, nullptr, nullptr, &mtop);
    top = gmx_mtop_t_to_t_topology(&mtop, FALSE);

    ip       = top.idef.iparams;
    iatom    = top.idef.il[F_ORIRES].iatoms;

    /* Count how many distance restraint there are... */
    nb = top.idef.il[F_ORIRES].nr;
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
    done_top_mtop(&top, &mtop);
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
                    nbounds, pair, bounds, nullptr, &sumt, &sumaver);

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

static void print_time(FILE *fp, double t)
{
    fprintf(fp, "%12.6f", t);
}

int gmx_nmr(int argc, char *argv[])
{
    const char        *desc[] = {
        "[THISMODULE] extracts distance or orientation restraint",
        "data from an energy file. The user is prompted to interactively",
        "select the desired terms.[PAR]",

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


    };
    static gmx_bool    bPrAll  = FALSE;
    static gmx_bool    bDp     = FALSE, bOrinst = FALSE, bOvec = FALSE;
    static int         skip    = 0;
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

    FILE              /* *out     = NULL,*/ *out_disre = NULL, *fp_pairs = NULL, *fort = NULL, *fodt = NULL, *foten = NULL;
    ener_file_t        fp;
    int                timecheck = 0;
    gmx_mtop_t         mtop;
    gmx_localtop_t    *top = nullptr;
    gmx_enxnm_t       *enm = nullptr;
    t_enxframe         fr;
    int                nre, teller, teller_disre;
    int                nor     = 0, nex = 0, norfr = 0, enx_i = 0;
    real              *bounds  = NULL, *violaver = NULL, *oobs = NULL, *orient = NULL, *odrms = NULL;
    int               *index   = NULL, *pair = NULL, norsel = 0, *orsel = NULL, *or_label = NULL;
    int                nbounds = 0, npairs;
    gmx_bool           bDisRe, bDRAll, bORA, bORT, bODA, bODR, bODT, bORIRE, bOTEN;
    gmx_bool           bCont;
    double             sumaver, sumt;
    int               *set     = NULL, i, j, k, nset, sss;
    char             **pairleg, **odtleg, **otenleg;
    char             **leg = nullptr;
    const char        *anm_j, *anm_k, *resnm_j, *resnm_k;
    int                resnr_j, resnr_k;
    const char        *orinst_sub = "@ subtitle \"instantaneous\"\n";
    char               buf[256];
    gmx_output_env_t  *oenv;
    t_enxblock        *blk_disre = nullptr;
    int                ndisre    = 0;

    t_filenm           fnm[] = {
        { efEDR, "-f",    nullptr,      ffREAD  },
        { efEDR, "-f2",   nullptr,      ffOPTRD },
        { efTPR, "-s",    nullptr,      ffOPTRD },
        // { efXVG, "-o",    "energy",  ffWRITE },
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

    npargs = asize(pa);
    if (!parse_common_args(&argc, argv,
                           PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END,
                           NFILE, fnm, npargs, pa, asize(desc), desc, 0, nullptr, &oenv))
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
    if (!(bDRAll || bDisRe || bORA || bORT || bODA || bODR || bODT || bORIRE || bOTEN))
    {
        printf("No output selected. Run with -h to see options. Terminating program.\n");
        return 0;
    }
    nset = 0;

    fp = open_enx(ftp2fn(efEDR, NFILE, fnm), "r");
    do_enxnms(fp, &nre, &enm);
    free_enxnms(nre, enm);

    t_inputrec  irInstance;
    t_inputrec *ir = &irInstance;
    init_enxframe(&fr);

    if (!bDisRe)
    {
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
                for (i = 0; i < norsel; i++)
                {
                    sfree(odtleg[i]);
                }
                sfree(odtleg);
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
            for (j = 0; j < 3; j++)
            {
                sfree(otenleg[j]);
            }
            sfree(otenleg);
        }
    }
    else
    {
        nbounds = get_bounds(ftp2fn(efTPR, NFILE, fnm), &bounds, &index, &pair, &npairs,
                             &mtop, &top, ir);
        snew(violaver, npairs);
        out_disre = xvgropen(opt2fn("-o", NFILE, fnm), "Sum of Violations",
                             "Time (ps)", "nm", oenv);
        xvgr_legend(out_disre, 2, drleg, oenv);
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

    /* Initiate counters */
    teller       = 0;
    teller_disre = 0;
    do
    {
        /* This loop searches for the first frame (when -b option is given),
         * or when this has been found it reads just one energy frame
         */
        do
        {
            bCont = do_enx(fp, &fr);
            if (bCont)
            {
                timecheck = check_times(fr.t);
            }
        }
        while (bCont && (timecheck < 0));

        if ((timecheck == 0) && bCont)
        {
            /*
             * Define distance restraint legends. Can only be done after
             * the first frame has been read... (Then we know how many there are)
             */
            blk_disre = find_block_id_enxframe(&fr, enxDISRE, nullptr);
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

                        print_time(out_disre, fr.t);
                        if (violaver == nullptr)
                        {
                            snew(violaver, ndisre);
                        }

                        /* Subtract bounds from distances, to calculate violations */
                        calc_violations(disre_rt, disre_rm3tav,
                                        nbounds, pair, bounds, violaver, &sumt, &sumaver);

                        fprintf(out_disre, "  %8.4f  %8.4f\n", sumaver, sumt);
                        if (bDRAll)
                        {
                            print_time(fp_pairs, fr.t);
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
                 * O R I R E S
                 *******************************************/
                else
                {
                    t_enxblock *blk = find_block_id_enxframe(&fr, enx_i, nullptr);
                    if (bORIRE && blk)
                    {
                        if (blk->nsub != 1)
                        {
                            gmx_fatal(FARGS, "Orientational restraints read in incorrectly.");
                        }

                        if (blk->sub[0].nr != nor)
                        {
                            gmx_fatal(FARGS, "Number of orientation restraints in energy file (%d) does not match with the topology (%d)", blk->sub[0].nr);
                        }
                        if (bORA || bODA)
                        {
                            for (i = 0; i < nor; i++)
                            {
                                orient[i] += blk_value(blk, 0, i);
                            }
                        }
                        if (bODR)
                        {
                            for (i = 0; i < nor; i++)
                            {
                                real v = blk_value(blk, 0, i);
                                odrms[i] += gmx::square(v - oobs[i]);
                            }
                        }
                        if (bORT)
                        {
                            fprintf(fort, "  %10f", fr.t);
                            for (i = 0; i < norsel; i++)
                            {
                                fprintf(fort, " %g", blk_value(blk, 0, orsel[i]));
                            }
                            fprintf(fort, "\n");
                        }
                        if (bODT)
                        {
                            fprintf(fodt, "  %10f", fr.t);
                            for (i = 0; i < norsel; i++)
                            {
                                fprintf(fodt, " %g", blk_value(blk, 0, orsel[i])-oobs[orsel[i]]);
                            }
                            fprintf(fodt, "\n");
                        }
                        norfr++;
                    }
                    blk = find_block_id_enxframe(&fr, enxORT, nullptr);
                    if (bOTEN && blk)
                    {
                        if (blk->nsub != 1)
                        {
                            gmx_fatal(FARGS, "Orientational restraints read in incorrectly");
                        }

                        if (blk->sub[0].nr != nex*12)
                        {
                            gmx_fatal(FARGS, "Number of orientation experiments in energy file (%g) does not match with the topology (%d)",
                                      blk->sub[0].nr/12, nex);
                        }
                        fprintf(foten, "  %10f", fr.t);
                        for (i = 0; i < nex; i++)
                        {
                            for (j = 0; j < (bOvec ? 12 : 3); j++)
                            {
                                fprintf(foten, " %g", blk_value(blk, 0, i*12+j));
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
    free_enxframe(&fr);

    fprintf(stderr, "\n");
    done_ener_file(fp);
    if (out_disre)
    {
        xvgrclose(out_disre);
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
        FILE *out = xvgropen(opt2fn("-ora", NFILE, fnm),
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
        FILE *out = xvgropen(opt2fn("-oda", NFILE, fnm),
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
        FILE *out = xvgropen(opt2fn("-odr", NFILE, fnm),
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
    // Clean up orires variables.
    sfree(or_label);
    sfree(oobs);
    sfree(orient);
    sfree(odrms);
    sfree(orsel);
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

        do_view(oenv, opt2fn_null("-ora", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-ort", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-oda", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-odr", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-odt", NFILE, fnm), nxy);
        do_view(oenv, opt2fn_null("-oten", NFILE, fnm), nxy);
    }
    output_env_done(oenv);
    done_filenms(NFILE, fnm);

    return 0;
}

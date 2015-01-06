/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

int gmx_dyecoupl(int argc, char *argv[])
{
    const char *desc[] =
    {
        "[THISMODULE] extracts dye dynamics from trajectory files.",
        "Currently, R and kappa^2 between dyes is extracted for (F)RET",
        "simulations with assumed dipolar coupling as in the Foerster equation.",
        "It further allows the calculation of R(t) and kappa^2(t), R and",
        "kappa^2 histograms and averages, as well as the instantaneous FRET",
        "efficiency E(t) for a specified Foerster radius R_0 (switch [TT]-R0[tt]).",
        "The input dyes have to be whole (see res and mol pbc options",
        "in [TT]trjconv[tt]).",
        "The dye transition dipole moment has to be defined by at least",
        "a single atom pair, however multiple atom pairs can be provided ",
        "in the index file. The distance R is calculated on the basis of",
        "the COMs of the given atom pairs.",
        "The [TT]-pbcdist[tt] option calculates distances to the nearest periodic",
        "image instead to the distance in the box. This works however only,"
        "for periodic boundaries in all 3 dimensions.",
        "The [TT]-norm[tt] option (area-) normalizes the histograms."
    };

    static gmx_bool bPBCdist = FALSE, bNormHist = FALSE;
    int             histbins = 50;
    output_env_t    oenv;
    real            R0 = -1;

    t_pargs         pa[] =
    {
        { "-pbcdist", FALSE, etBOOL, { &bPBCdist }, "Distance R based on PBC" },
        { "-norm", FALSE, etBOOL, { &bNormHist }, "Normalize histograms" },
        { "-bins", FALSE, etINT, {&histbins}, "# of histogram bins" },
        { "-R0", FALSE, etREAL, {&R0}, "Foerster radius including kappa^2=2/3 in nm" }
    };
#define NPA asize(pa)

    t_filenm fnm[] =
    {
        { efTRX, "-f", NULL, ffREAD },
        { efNDX, NULL, NULL, ffREAD },
        { efXVG, "-ot", "rkappa", ffOPTWR },
        { efXVG, "-oe", "insteff", ffOPTWR },
        { efDAT, "-o", "rkappa", ffOPTWR },
        { efXVG, "-rhist", "rhist", ffOPTWR },
        { efXVG, "-khist", "khist", ffOPTWR }
    };
#define NFILE asize(fnm)


    const char  *in_trajfile, *in_ndxfile, *out_xvgrkfile = NULL, *out_xvginstefffile = NULL, *out_xvgrhistfile = NULL, *out_xvgkhistfile = NULL, *out_datfile = NULL;
    gmx_bool     bHaveFirstFrame, bHaveNextFrame, indexOK = TRUE;
    int          ndon, nacc;
    atom_id     *donindex, *accindex;
    char        *grpnm;
    t_trxstatus *status;
    t_trxframe   fr;

    int          flags;
    int          allocblock = 1000;
    real         histexpand = 1e-6;
    rvec         donvec, accvec, donpos, accpos, dist, distnorm;
    int          natoms;

    /*we rely on PBC autodetection (...currently)*/
    int         ePBC = -1;

    real       *rvalues = NULL, *kappa2values = NULL, *rhist = NULL, *khist = NULL;
    t_pbc      *pbc     = NULL;
    int         i, bin;
    FILE       *rkfp = NULL, *rhfp = NULL, *khfp = NULL, *datfp = NULL, *iefp = NULL;
    gmx_bool    bRKout, bRhistout, bKhistout, bDatout, bInstEffout, grident;

    const char *rkleg[2] = { "R", "\\f{Symbol}k\\f{}\\S2\\N" };
    const char *rhleg[1] = { "p(R)" };
    const char *khleg[1] = { "p(\\f{Symbol}k\\f{}\\S2\\N)" };
    const char *ieleg[1] = { "E\\sRET\\N(t)" };

    real        R, kappa2, insteff, Rs = 0., kappa2s = 0., insteffs = 0., rmax, rmin, kmin = 0., kmax = 4.,
                rrange, krange, rincr, kincr, Rfrac;
    int         rkcount = 0, rblocksallocated = 0, kblocksallocated = 0;

    if (!parse_common_args(&argc, argv, PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_VIEW | PCA_TIME_UNIT,
                           NFILE, fnm, NPA, pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }


    /* Check command line options for filenames and set bool flags when switch used*/
    in_trajfile        = opt2fn("-f", NFILE, fnm);
    in_ndxfile         = opt2fn("-n", NFILE, fnm);
    out_xvgrkfile      = opt2fn("-ot", NFILE, fnm);
    out_xvgrhistfile   = opt2fn("-rhist", NFILE, fnm);
    out_xvgkhistfile   = opt2fn("-khist", NFILE, fnm);
    out_xvginstefffile = opt2fn("-oe", NFILE, fnm);
    out_datfile        = opt2fn("-o", NFILE, fnm);

    bRKout      = opt2bSet("-ot", NFILE, fnm);
    bRhistout   = opt2bSet("-rhist", NFILE, fnm);
    bKhistout   = opt2bSet("-khist", NFILE, fnm);
    bDatout     = opt2bSet("-o", NFILE, fnm);
    bInstEffout = opt2bSet("-oe", NFILE, fnm);


    /* PBC warning. */
    if (bPBCdist)
    {
        printf("Calculating distances to periodic image.\n");
        printf("Be careful! This produces only valid results for PBC in all three dimensions\n");
    }


    if (bInstEffout && R0 <= 0.)
    {
        gmx_fatal(FARGS, "You have to specify R0 and R0 has to be larger than 0 nm.\n\n");
    }

    printf("Select group with donor atom pairs defining the transition moment\n");
    get_index(NULL, ftp2fn_null(efNDX, NFILE, fnm), 1, &ndon, &donindex, &grpnm);

    printf("Select group with acceptor atom pairs defining the transition moment\n");
    get_index(NULL, ftp2fn_null(efNDX, NFILE, fnm), 1, &nacc, &accindex, &grpnm);

    /*check if groups are identical*/
    grident = TRUE;

    if (ndon == nacc)
    {
        for (i = 0; i < nacc; i++)
        {
            if (accindex[i] != donindex[i])
            {
                grident = FALSE;
                break;
            }
        }
    }

    if (grident)
    {
        gmx_fatal(FARGS, "Donor and acceptor group are identical. This makes no sense.");
    }

    printf("Reading first frame\n");
    /* open trx file for reading */
    flags           = 0;
    flags           = flags | TRX_READ_X;
    bHaveFirstFrame = read_first_frame(oenv, &status, in_trajfile, &fr, flags);

    if (bHaveFirstFrame)
    {
        printf("First frame is OK\n");
        natoms = fr.natoms;
        if ((ndon % 2 != 0) || (nacc % 2 != 0))
        {
            indexOK = FALSE;
        }
        else
        {
            for (i = 0; i < ndon; i++)
            {
                if (donindex[i] >= natoms)
                {
                    indexOK = FALSE;
                }
            }
            for (i = 0; i < nacc; i++)
            {
                if (accindex[i] >= natoms)
                {
                    indexOK = FALSE;
                }
            }
        }

        if (indexOK)
        {

            if (bDatout)
            {
                datfp = fopen(out_datfile, "w");
            }

            if (bRKout)
            {
                rkfp = xvgropen(out_xvgrkfile,
                                "Distance and \\f{Symbol}k\\f{}\\S2\\N trajectory",
                                "Time (ps)", "Distance (nm) / \\f{Symbol}k\\f{}\\S2\\N",
                                oenv);
                xvgr_legend(rkfp, 2, rkleg, oenv);
            }

            if (bInstEffout)
            {
                iefp = xvgropen(out_xvginstefffile,
                                "Instantaneous RET Efficiency",
                                "Time (ps)", "RET Efficiency",
                                oenv);
                xvgr_legend(iefp, 1, ieleg, oenv);
            }


            if (bRhistout)
            {
                snew(rvalues, allocblock);
                rblocksallocated += 1;
                snew(rhist, histbins);
            }

            if (bKhistout)
            {
                snew(kappa2values, allocblock);
                kblocksallocated += 1;
                snew(khist, histbins);
            }

            do
            {
                clear_rvec(donvec);
                clear_rvec(accvec);
                clear_rvec(donpos);
                clear_rvec(accpos);
                for (i = 0; i < ndon / 2; i++)
                {
                    rvec_sub(donvec, fr.x[donindex[2 * i]], donvec);
                    rvec_add(donvec, fr.x[donindex[2 * i + 1]], donvec);
                    rvec_add(donpos, fr.x[donindex[2 * i]], donpos);
                    rvec_add(donpos, fr.x[donindex[2 * i + 1]], donpos);
                }

                for (i = 0; i < nacc / 2; i++)
                {
                    rvec_sub(accvec, fr.x[accindex[2 * i]], accvec);
                    rvec_add(accvec, fr.x[accindex[2 * i + 1]], accvec);
                    rvec_add(accpos, fr.x[accindex[2 * i]], accpos);
                    rvec_add(accpos, fr.x[accindex[2 * i + 1]], accpos);
                }

                unitv(donvec, donvec);
                unitv(accvec, accvec);

                svmul((real) 1. / ndon, donpos, donpos);
                svmul((real) 1. / nacc, accpos, accpos);

                if (bPBCdist)
                {
                    set_pbc(pbc, ePBC, fr.box);
                    pbc_dx(pbc, donpos, accpos, dist);
                }
                else
                {
                    rvec_sub(donpos, accpos, dist);
                }

                unitv(dist, distnorm);
                R       = norm(dist);
                kappa2  = iprod(donvec, accvec)- 3.* (iprod(donvec, distnorm) * iprod(distnorm, accvec));
                kappa2 *= kappa2;
                if (R0 > 0)
                {
                    Rfrac     = R/R0;
                    insteff   = 1/(1+(Rfrac*Rfrac*Rfrac*Rfrac*Rfrac*Rfrac)*2/3/kappa2);
                    insteffs += insteff;

                    if (bInstEffout)
                    {
                        fprintf(iefp, "%12.7f %12.7f\n", fr.time, insteff);
                    }
                }


                Rs      += R;
                kappa2s += kappa2;
                rkcount++;

                if (bRKout)
                {
                    fprintf(rkfp, "%12.7f %12.7f %12.7f\n", fr.time, R, kappa2);
                }

                if (bDatout)
                {
                    fprintf(datfp, "%12.7f %12.7f %12.7f\n", fr.time, R, kappa2);
                }

                if (bRhistout)
                {
                    rvalues[rkcount-1] = R;
                    if (rkcount % allocblock == 0)
                    {
                        srenew(rvalues, allocblock*(rblocksallocated+1));
                        rblocksallocated += 1;
                    }
                }

                if (bKhistout)
                {
                    kappa2values[rkcount-1] = kappa2;
                    if (rkcount % allocblock == 0)
                    {
                        srenew(kappa2values, allocblock*(kblocksallocated+1));
                        kblocksallocated += 1;
                    }
                }

                bHaveNextFrame = read_next_frame(oenv, status, &fr);
            }
            while (bHaveNextFrame);

            if (bRKout)
            {
                xvgrclose(rkfp);
            }

            if (bDatout)
            {
                gmx_ffclose(datfp);
            }

            if (bInstEffout)
            {
                xvgrclose(iefp);
            }


            if (bRhistout)
            {
                printf("Writing R-Histogram\n");
                rmin = rvalues[0];
                rmax = rvalues[0];
                for (i = 1; i < rkcount; i++)
                {
                    if (rvalues[i] < rmin)
                    {
                        rmin = rvalues[i];
                    }
                    else if (rvalues[i] > rmax)
                    {
                        rmax = rvalues[i];
                    }
                }
                rmin -= histexpand;
                rmax += histexpand;

                rrange = rmax - rmin;
                rincr  = rrange / histbins;

                for (i = 1; i < rkcount; i++)
                {
                    bin         = (int) ((rvalues[i] - rmin) / rincr);
                    rhist[bin] += 1;
                }
                if (bNormHist)
                {
                    for (i = 0; i < histbins; i++)
                    {
                        rhist[i] /= rkcount * rrange/histbins;
                    }
                    rhfp = xvgropen(out_xvgrhistfile, "Distance Distribution",
                                    "R (nm)", "Normalized Probability", oenv);
                }
                else
                {
                    rhfp = xvgropen(out_xvgrhistfile, "Distance Distribution",
                                    "R (nm)", "Probability", oenv);
                }
                xvgr_legend(rhfp, 1, rhleg, oenv);
                for (i = 0; i < histbins; i++)
                {
                    fprintf(rhfp, "%12.7f %12.7f\n", (i + 0.5) * rincr + rmin,
                            rhist[i]);
                }
                xvgrclose(rhfp);
            }

            if (bKhistout)
            {
                printf("Writing kappa^2-Histogram\n");
                krange = kmax - kmin;
                kincr  = krange / histbins;

                for (i = 1; i < rkcount; i++)
                {
                    bin         = (int) ((kappa2values[i] - kmin) / kincr);
                    khist[bin] += 1;
                }
                if (bNormHist)
                {
                    for (i = 0; i < histbins; i++)
                    {
                        khist[i] /= rkcount * krange/histbins;
                    }
                    khfp = xvgropen(out_xvgkhistfile,
                                    "\\f{Symbol}k\\f{}\\S2\\N Distribution",
                                    "\\f{Symbol}k\\f{}\\S2\\N",
                                    "Normalized Probability", oenv);
                }
                else
                {
                    khfp = xvgropen(out_xvgkhistfile,
                                    "\\f{Symbol}k\\f{}\\S2\\N Distribution",
                                    "\\f{Symbol}k\\f{}\\S2\\N", "Probability", oenv);
                }
                xvgr_legend(khfp, 1, khleg, oenv);
                for (i = 0; i < histbins; i++)
                {
                    fprintf(khfp, "%12.7f %12.7f\n", (i + 0.5) * kincr + kmin,
                            khist[i]);
                }
                xvgrclose(khfp);
            }

            printf("\nAverages:\n");
            printf("R_avg   = %8.4f nm\nKappa^2 = %8.4f\n", Rs / rkcount,
                   kappa2s / rkcount);
            if (R0 > 0)
            {
                printf("E_RETavg   = %8.4f\n", insteffs / rkcount);
            }
            please_cite(stdout, "Hoefling2011");
        }
        else
        {
            gmx_fatal(FARGS, "Index file invalid, check your index file for correct pairs.\n");
        }
    }
    else
    {
        gmx_fatal(FARGS, "Could not read first frame of the trajectory.\n");
    }

    return 0;
}

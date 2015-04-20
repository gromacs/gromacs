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

#include "config.h"

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/gmxana/nsfactor.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

int gmx_sans(int argc, char *argv[])
{
    const char          *desc[] = {
        "[THISMODULE] computes SANS spectra using Debye formula.",
        "It currently uses topology file (since it need to assigne element for each atom).",
        "[PAR]",
        "Parameters:[PAR]"
        "[TT]-pr[tt] Computes normalized g(r) function averaged over trajectory[PAR]",
        "[TT]-prframe[tt] Computes normalized g(r) function for each frame[PAR]",
        "[TT]-sq[tt] Computes SANS intensity curve averaged over trajectory[PAR]",
        "[TT]-sqframe[tt] Computes SANS intensity curve for each frame[PAR]",
        "[TT]-startq[tt] Starting q value in nm[PAR]",
        "[TT]-endq[tt] Ending q value in nm[PAR]",
        "[TT]-qstep[tt] Stepping in q space[PAR]",
        "Note: When using Debye direct method computational cost increases as",
        "1/2 * N * (N - 1) where N is atom number in group of interest.",
        "[PAR]",
        "WARNING: If sq or pr specified this tool can produce large number of files! Up to two times larger than number of frames!"
    };
    static gmx_bool      bPBC     = TRUE;
    static gmx_bool      bNORM    = FALSE;
    static real          binwidth = 0.2, grid = 0.05; /* bins shouldnt be smaller then smallest bond (~0.1nm) length */
    static real          start_q  = 0.0, end_q = 2.0, q_step = 0.01;
    static real          mcover   = -1;
    static unsigned int  seed     = 0;
    static int           nthreads = -1;

    static const char   *emode[]   = { NULL, "direct", "mc", NULL };
    static const char   *emethod[] = { NULL, "debye", "fft", NULL };

    gmx_neutron_atomic_structurefactors_t    *gnsf;
    gmx_sans_t                               *gsans;

#define NPA asize(pa)

    t_pargs                               pa[] = {
        { "-bin", FALSE, etREAL, {&binwidth},
          "[HIDDEN]Binwidth (nm)" },
        { "-mode", FALSE, etENUM, {emode},
          "Mode for sans spectra calculation" },
        { "-mcover", FALSE, etREAL, {&mcover},
          "Monte-Carlo coverage should be -1(default) or (0,1]"},
        { "-method", FALSE, etENUM, {emethod},
          "[HIDDEN]Method for sans spectra calculation" },
        { "-pbc", FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances" },
        { "-grid", FALSE, etREAL, {&grid},
          "[HIDDEN]Grid spacing (in nm) for FFTs" },
        {"-startq", FALSE, etREAL, {&start_q},
         "Starting q (1/nm) "},
        {"-endq", FALSE, etREAL, {&end_q},
         "Ending q (1/nm)"},
        { "-qstep", FALSE, etREAL, {&q_step},
          "Stepping in q (1/nm)"},
        { "-seed",     FALSE, etINT,  {&seed},
          "Random seed for Monte-Carlo"},
#ifdef GMX_OPENMP
        { "-nt",  FALSE, etINT, {&nthreads},
          "Number of threads to start"},
#endif
    };
    FILE                                 *fp;
    const char                           *fnTPX, *fnNDX, *fnTRX, *fnDAT = NULL;
    t_trxstatus                          *status;
    t_topology                           *top  = NULL;
    t_atom                               *atom = NULL;
    gmx_rmpbc_t                           gpbc = NULL;
    gmx_bool                              bTPX;
    gmx_bool                              bFFT = FALSE, bDEBYE = FALSE;
    gmx_bool                              bMC  = FALSE;
    int                                   ePBC = -1;
    matrix                                box;
    char                                  title[STRLEN];
    rvec                                 *x;
    int                                   natoms;
    real                                  t;
    char                                **grpname = NULL;
    atom_id                              *index   = NULL;
    int                                   isize;
    int                                   i, j;
    char                                 *hdr            = NULL;
    char                                 *suffix         = NULL;
    t_filenm                             *fnmdup         = NULL;
    gmx_radial_distribution_histogram_t  *prframecurrent = NULL, *pr = NULL;
    gmx_static_structurefactor_t         *sqframecurrent = NULL, *sq = NULL;
    output_env_t                          oenv;

#define NFILE asize(fnm)

    t_filenm   fnm[] = {
        { efTPR,  "-s",       NULL,       ffREAD },
        { efTRX,  "-f",       NULL,       ffREAD },
        { efNDX,  NULL,       NULL,       ffOPTRD },
        { efDAT,  "-d",       "nsfactor", ffOPTRD },
        { efXVG,  "-pr",      "pr",       ffWRITE },
        { efXVG,  "-sq",       "sq",      ffWRITE },
        { efXVG,  "-prframe", "prframe",  ffOPTWR },
        { efXVG,  "-sqframe", "sqframe",  ffOPTWR }
    };

    nthreads = gmx_omp_get_max_threads();

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_TIME_UNIT,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    /* check that binwidth not smaller than smallers distance */
    check_binwidth(binwidth);
    check_mcover(mcover);

    /* setting number of omp threads globaly */
    gmx_omp_set_num_threads(nthreads);

    /* Now try to parse opts for modes */
    switch (emethod[0][0])
    {
        case 'd':
            bDEBYE = TRUE;
            switch (emode[0][0])
            {
                case 'd':
                    bMC = FALSE;
                    break;
                case 'm':
                    bMC = TRUE;
                    break;
                default:
                    break;
            }
            break;
        case 'f':
            bFFT = TRUE;
            break;
        default:
            break;
    }

    if (bDEBYE)
    {
        if (bMC)
        {
            fprintf(stderr, "Using Monte Carlo Debye method to calculate spectrum\n");
        }
        else
        {
            fprintf(stderr, "Using direct Debye method to calculate spectrum\n");
        }
    }
    else if (bFFT)
    {
        gmx_fatal(FARGS, "FFT method not implemented!");
    }
    else
    {
        gmx_fatal(FARGS, "Unknown combination for mode and method!");
    }

    /* Try to read files */
    fnDAT = ftp2fn(efDAT, NFILE, fnm);
    fnTPX = ftp2fn(efTPR, NFILE, fnm);
    fnTRX = ftp2fn(efTRX, NFILE, fnm);

    gnsf = gmx_neutronstructurefactors_init(fnDAT);
    fprintf(stderr, "Read %d atom names from %s with neutron scattering parameters\n\n", gnsf->nratoms, fnDAT);

    snew(top, 1);
    snew(grpname, 1);
    snew(index, 1);

    bTPX = read_tps_conf(fnTPX, title, top, &ePBC, &x, NULL, box, TRUE);

    printf("\nPlease select group for SANS spectra calculation:\n");
    get_index(&(top->atoms), ftp2fn_null(efNDX, NFILE, fnm), 1, &isize, &index, grpname);

    gsans = gmx_sans_init(top, gnsf);

    /* Prepare reference frame */
    if (bPBC)
    {
        gpbc = gmx_rmpbc_init(&top->idef, ePBC, top->atoms.nr);
        gmx_rmpbc(gpbc, top->atoms.nr, box, x);
    }

    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);
    if (natoms != top->atoms.nr)
    {
        fprintf(stderr, "\nWARNING: number of atoms in tpx (%d) and trajectory (%d) do not match\n", natoms, top->atoms.nr);
    }

    do
    {
        if (bPBC)
        {
            gmx_rmpbc(gpbc, top->atoms.nr, box, x);
        }
        /* allocate memory for pr */
        if (pr == NULL)
        {
            /* in case its first frame to read */
            snew(pr, 1);
        }
        /*  realy calc p(r) */
        prframecurrent = calc_radial_distribution_histogram(gsans, x, box, index, isize, binwidth, bMC, bNORM, mcover, seed);
        /* copy prframecurrent -> pr and summ up pr->gr[i] */
        /* allocate and/or resize memory for pr->gr[i] and pr->r[i] */
        if (pr->gr == NULL)
        {
            /* check if we use pr->gr first time */
            snew(pr->gr, prframecurrent->grn);
            snew(pr->r, prframecurrent->grn);
        }
        else
        {
            /* resize pr->gr and pr->r if needed to preven overruns */
            if (prframecurrent->grn > pr->grn)
            {
                srenew(pr->gr, prframecurrent->grn);
                srenew(pr->r, prframecurrent->grn);
            }
        }
        pr->grn      = prframecurrent->grn;
        pr->binwidth = prframecurrent->binwidth;
        /* summ up gr and fill r */
        for (i = 0; i < prframecurrent->grn; i++)
        {
            pr->gr[i] += prframecurrent->gr[i];
            pr->r[i]   = prframecurrent->r[i];
        }
        /* normalize histo */
        normalize_probability(prframecurrent->grn, prframecurrent->gr);
        /* convert p(r) to sq */
        sqframecurrent = convert_histogram_to_intensity_curve(prframecurrent, start_q, end_q, q_step);
        /* print frame data if needed */
        if (opt2fn_null("-prframe", NFILE, fnm))
        {
            snew(hdr, 25);
            snew(suffix, GMX_PATH_MAX);
            /* prepare header */
            sprintf(hdr, "g(r), t = %f", t);
            /* prepare output filename */
            fnmdup = dup_tfn(NFILE, fnm);
            sprintf(suffix, "-t%.2f", t);
            add_suffix_to_output_names(fnmdup, NFILE, suffix);
            fp = xvgropen(opt2fn_null("-prframe", NFILE, fnmdup), hdr, "Distance (nm)", "Probability", oenv);
            for (i = 0; i < prframecurrent->grn; i++)
            {
                fprintf(fp, "%10.6f%10.6f\n", prframecurrent->r[i], prframecurrent->gr[i]);
            }
            done_filenms(NFILE, fnmdup);
            xvgrclose(fp);
            sfree(hdr);
            sfree(suffix);
            sfree(fnmdup);
        }
        if (opt2fn_null("-sqframe", NFILE, fnm))
        {
            snew(hdr, 25);
            snew(suffix, GMX_PATH_MAX);
            /* prepare header */
            sprintf(hdr, "I(q), t = %f", t);
            /* prepare output filename */
            fnmdup = dup_tfn(NFILE, fnm);
            sprintf(suffix, "-t%.2f", t);
            add_suffix_to_output_names(fnmdup, NFILE, suffix);
            fp = xvgropen(opt2fn_null("-sqframe", NFILE, fnmdup), hdr, "q (nm^-1)", "s(q)/s(0)", oenv);
            for (i = 0; i < sqframecurrent->qn; i++)
            {
                fprintf(fp, "%10.6f%10.6f\n", sqframecurrent->q[i], sqframecurrent->s[i]);
            }
            done_filenms(NFILE, fnmdup);
            xvgrclose(fp);
            sfree(hdr);
            sfree(suffix);
            sfree(fnmdup);
        }
        /* free pr structure */
        sfree(prframecurrent->gr);
        sfree(prframecurrent->r);
        sfree(prframecurrent);
        /* free sq structure */
        sfree(sqframecurrent->q);
        sfree(sqframecurrent->s);
        sfree(sqframecurrent);
    }
    while (read_next_x(oenv, status, &t, x, box));
    close_trj(status);

    /* normalize histo */
    normalize_probability(pr->grn, pr->gr);
    sq = convert_histogram_to_intensity_curve(pr, start_q, end_q, q_step);
    /* prepare pr.xvg */
    fp = xvgropen(opt2fn_null("-pr", NFILE, fnm), "G(r)", "Distance (nm)", "Probability", oenv);
    for (i = 0; i < pr->grn; i++)
    {
        fprintf(fp, "%10.6f%10.6f\n", pr->r[i], pr->gr[i]);
    }
    xvgrclose(fp);

    /* prepare sq.xvg */
    fp = xvgropen(opt2fn_null("-sq", NFILE, fnm), "I(q)", "q (nm^-1)", "s(q)/s(0)", oenv);
    for (i = 0; i < sq->qn; i++)
    {
        fprintf(fp, "%10.6f%10.6f\n", sq->q[i], sq->s[i]);
    }
    xvgrclose(fp);
    /*
     * Clean up memory
     */
    sfree(pr->gr);
    sfree(pr->r);
    sfree(pr);
    sfree(sq->q);
    sfree(sq->s);
    sfree(sq);

    please_cite(stdout, "Garmay2012");

    return 0;
}

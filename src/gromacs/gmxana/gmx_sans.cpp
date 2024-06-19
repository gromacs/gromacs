/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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

#include "config.h"

#include <cstdio>

#include <array>
#include <filesystem>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/nsfactor.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

struct gmx_output_env_t;

int gmx_sans(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] computes SANS spectra using Debye formula.",
        "It currently uses topology file (since it need to assign element for each atom).",
        "[PAR]",
        "Parameters:[PAR]",
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
        "WARNING: If sq or pr specified this tool can produce large number of files! Up to ",
        "two times larger than number of frames!"
    };
    static gmx_bool bPBC     = TRUE;
    static gmx_bool bNORM    = FALSE;
    static real     binwidth = 0.2,
                grid = 0.05; /* bins shouldn't be smaller then smallest bond (~0.1nm) length */
    static real         start_q = 0.0, end_q = 2.0, q_step = 0.01;
    static real         mcover   = -1;
    static unsigned int seed     = 0;
    static int          nthreads = -1;

    static const char* emode[]   = { nullptr, "direct", "mc", nullptr };
    static const char* emethod[] = { nullptr, "debye", "fft", nullptr };

    gmx_neutron_atomic_structurefactors_t* gnsf;
    gmx_sans_t*                            gsans;

    t_pargs pa[] = {
        { "-bin", FALSE, etREAL, { &binwidth }, "HIDDENBinwidth (nm)" },
        { "-mode", FALSE, etENUM, { emode }, "Mode for sans spectra calculation" },
        { "-mcover",
          FALSE,
          etREAL,
          { &mcover },
          "Monte-Carlo coverage should be -1(default) or (0,1]" },
        { "-method", FALSE, etENUM, { emethod }, "HIDDENMethod for sans spectra calculation" },
        { "-pbc",
          FALSE,
          etBOOL,
          { &bPBC },
          "Use periodic boundary conditions for computing distances" },
        { "-grid", FALSE, etREAL, { &grid }, "HIDDENGrid spacing (in nm) for FFTs" },
        { "-startq", FALSE, etREAL, { &start_q }, "Starting q (1/nm) " },
        { "-endq", FALSE, etREAL, { &end_q }, "Ending q (1/nm)" },
        { "-qstep", FALSE, etREAL, { &q_step }, "Stepping in q (1/nm)" },
        { "-seed", FALSE, etINT, { &seed }, "Random seed for Monte-Carlo" },
#if GMX_OPENMP
        { "-nt", FALSE, etINT, { &nthreads }, "Number of threads to start" },
#endif
    };
    FILE*                                fp;
    const char *                         fnTPX, *fnTRX, *fnDAT = nullptr;
    t_trxstatus*                         status;
    t_topology*                          top  = nullptr;
    gmx_rmpbc_t                          gpbc = nullptr;
    gmx_bool                             bFFT = FALSE, bDEBYE = FALSE;
    gmx_bool                             bMC     = FALSE;
    PbcType                              pbcType = PbcType::Unset;
    matrix                               box;
    rvec*                                x;
    int                                  natoms;
    real                                 t;
    char**                               grpname = nullptr;
    int*                                 index   = nullptr;
    int                                  isize;
    int                                  i;
    gmx_radial_distribution_histogram_t *prframecurrent = nullptr, *pr = nullptr;
    gmx_static_structurefactor_t *       sqframecurrent = nullptr, *sq = nullptr;
    gmx_output_env_t*                    oenv;

    std::array<t_filenm, 8> filenames = { { { efTPR, "-s", nullptr, ffREAD },
                                            { efTRX, "-f", nullptr, ffREAD },
                                            { efNDX, nullptr, nullptr, ffOPTRD },
                                            { efDAT, "-d", "nsfactor", ffOPTRD },
                                            { efXVG, "-pr", "pr", ffWRITE },
                                            { efXVG, "-sq", "sq", ffWRITE },
                                            { efXVG, "-prframe", "prframe", ffOPTWR },
                                            { efXVG, "-sqframe", "sqframe", ffOPTWR } } };
    t_filenm*               fnm       = filenames.data();

    const auto NFILE = filenames.size();

    nthreads = gmx_omp_get_max_threads();

    if (!parse_common_args(
                &argc, argv, PCA_CAN_TIME | PCA_TIME_UNIT, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    std::fprintf(stdout,
                 "You are going to use a deprecated gmx tool. Please migrate to the new one, gmx "
                 "scattering");

    /* check that binwidth not smaller than smallers distance */
    check_binwidth(binwidth);
    check_mcover(mcover);

    /* setting number of omp threads globaly */
    gmx_omp_set_num_threads(nthreads);

    /* Now try to parse opts for modes */
    GMX_RELEASE_ASSERT(emethod[0] != nullptr, "Options inconsistency; emethod[0] is NULL");
    switch (emethod[0][0])
    {
        case 'd':
            bDEBYE = TRUE;
            switch (emode[0][0])
            {
                case 'd': bMC = FALSE; break;
                case 'm': bMC = TRUE; break;
                default: break;
            }
            break;
        case 'f': bFFT = TRUE; break;
        default: break;
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

    read_tps_conf(fnTPX, top, &pbcType, &x, nullptr, box, TRUE);

    printf("\nPlease select group for SANS spectra calculation:\n");
    get_index(&(top->atoms), ftp2fn_null(efNDX, NFILE, fnm), 1, &isize, &index, grpname);

    gsans = gmx_sans_init(top, gnsf);

    /* Prepare reference frame */
    if (bPBC)
    {
        gpbc = gmx_rmpbc_init(&top->idef, pbcType, top->atoms.nr);
        gmx_rmpbc_apply(gpbc, top->atoms.nr, box, x);
    }

    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);
    if (natoms != top->atoms.nr)
    {
        fprintf(stderr,
                "\nWARNING: number of atoms in tpx (%d) and trajectory (%d) do not match\n",
                natoms,
                top->atoms.nr);
    }

    do
    {
        if (bPBC)
        {
            gmx_rmpbc_apply(gpbc, top->atoms.nr, box, x);
        }
        /* allocate memory for pr */
        if (pr == nullptr)
        {
            /* in case its first frame to read */
            snew(pr, 1);
        }
        /*  realy calc p(r) */
        prframecurrent = calc_radial_distribution_histogram(
                gsans, x, box, index, isize, binwidth, bMC, bNORM, mcover, seed);
        /* copy prframecurrent -> pr and summ up pr->gr[i] */
        /* allocate and/or resize memory for pr->gr[i] and pr->r[i] */
        if (pr->gr == nullptr)
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
            pr->r[i] = prframecurrent->r[i];
        }
        /* normalize histo */
        normalize_probability(prframecurrent->grn, prframecurrent->gr);
        /* convert p(r) to sq */
        sqframecurrent = convert_histogram_to_intensity_curve(prframecurrent, start_q, end_q, q_step);
        /* print frame data if needed */
        if (opt2fn_null("-prframe", NFILE, fnm))
        {
            /* prepare header */
            auto hdr    = gmx::formatString("g(r), t = %f", t);
            auto fnmdup = filenames;
            auto suffix = gmx::formatString("-t%.2f", t);
            add_suffix_to_output_names(fnmdup.data(), NFILE, suffix.c_str());
            fp = xvgropen(opt2fn_null("-prframe", NFILE, fnmdup.data()),
                          hdr.c_str(),
                          "Distance (nm)",
                          "Probability",
                          oenv);
            for (i = 0; i < prframecurrent->grn; i++)
            {
                fprintf(fp, "%10.6f%10.6f\n", prframecurrent->r[i], prframecurrent->gr[i]);
            }
            xvgrclose(fp);
        }
        if (opt2fn_null("-sqframe", NFILE, fnm))
        {
            /* prepare header */
            auto hdr = gmx::formatString("I(q), t = %f", t);
            /* prepare output filename */
            auto fnmdup = filenames;
            auto suffix = gmx::formatString("-t%.2f", t);
            add_suffix_to_output_names(fnmdup.data(), NFILE, suffix.c_str());
            fp = xvgropen(opt2fn_null("-sqframe", NFILE, fnmdup.data()),
                          hdr.c_str(),
                          "q (nm^-1)",
                          "s(q)/s(0)",
                          oenv);
            for (i = 0; i < sqframecurrent->qn; i++)
            {
                fprintf(fp, "%10.6f%10.6f\n", sqframecurrent->q[i], sqframecurrent->s[i]);
            }
            xvgrclose(fp);
        }
        /* free pr structure */
        sfree(prframecurrent->gr);
        sfree(prframecurrent->r);
        sfree(prframecurrent);
        /* free sq structure */
        sfree(sqframecurrent->q);
        sfree(sqframecurrent->s);
        sfree(sqframecurrent);
    } while (read_next_x(oenv, status, &t, x, box));
    close_trx(status);

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

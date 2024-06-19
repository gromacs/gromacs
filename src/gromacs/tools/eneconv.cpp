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

#include "eneconv.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/listed_forces/disre.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/int64_to_int.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"

struct gmx_output_env_t;

#define TIME_EXPLICIT 0
#define TIME_CONTINUE 1
#define TIME_LAST 2
#ifndef FLT_MAX
#    define FLT_MAX 1e36
#endif

static int* select_it(int nre, gmx_enxnm_t* nm, int* nset)
{
    gmx_bool* bE;
    int       n, k, j, i;
    int*      set;
    gmx_bool  bVerbose = TRUE;

    if ((getenv("GMX_ENER_VERBOSE")) != nullptr)
    {
        bVerbose = FALSE;
    }

    fprintf(stderr, "Select the terms you want to scale from the following list\n");
    fprintf(stderr, "End your selection with 0\n");

    if (bVerbose)
    {
        for (k = 0; (k < nre);)
        {
            for (j = 0; (j < 4) && (k < nre); j++, k++)
            {
                fprintf(stderr, " %3d=%14s", k + 1, nm[k].name);
            }
            fprintf(stderr, "\n");
        }
    }

    snew(bE, nre);
    do
    {
        if (1 != scanf("%d", &n))
        {
            gmx_fatal(FARGS, "Cannot read energy term");
        }
        if ((n > 0) && (n <= nre))
        {
            bE[n - 1] = TRUE;
        }
    } while (n != 0);

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

static void sort_files(gmx::ArrayRef<std::string> files, real* settime)
{
    for (gmx::Index i = 0; i < files.ssize(); i++)
    {
        gmx::Index minidx = i;
        for (gmx::Index j = i + 1; j < files.ssize(); j++)
        {
            if (settime[j] < settime[minidx])
            {
                minidx = j;
            }
        }
        if (minidx != i)
        {
            real timeswap   = settime[i];
            settime[i]      = settime[minidx];
            settime[minidx] = timeswap;
            std::string tmp = files[i];
            files[i]        = files[minidx];
            files[minidx]   = tmp;
        }
    }
}


static int scan_ene_files(const std::vector<std::string>& files, real* readtime, real* timestep, int* nremax)
{
    /* Check number of energy terms and start time of all files */
    int          nre, nremin = 0, nresav = 0;
    ener_file_t  in;
    real         t1, t2;
    char         inputstring[STRLEN];
    gmx_enxnm_t* enm;
    t_enxframe*  fr;

    snew(fr, 1);

    for (size_t f = 0; f < files.size(); f++)
    {
        in  = open_enx(files[f].c_str(), "r");
        enm = nullptr;
        do_enxnms(in, &nre, &enm);

        if (f == 0)
        {
            nresav  = nre;
            nremin  = nre;
            *nremax = nre;
            do_enx(in, fr);
            t1 = fr->t;
            do_enx(in, fr);
            t2          = fr->t;
            *timestep   = t2 - t1;
            readtime[f] = t1;
            close_enx(in);
        }
        else
        {
            nremin  = std::min(nremin, fr->nre);
            *nremax = std::max(*nremax, fr->nre);
            if (nre != nresav)
            {
                fprintf(stderr,
                        "Energy files don't match, different number of energies:\n"
                        " %s: %d\n %s: %d\n",
                        files[f - 1].c_str(),
                        nresav,
                        files[f].c_str(),
                        fr->nre);
                fprintf(stderr,
                        "\nContinue conversion using only the first %d terms (n/y)?\n"
                        "(you should be sure that the energy terms match)\n",
                        nremin);
                if (nullptr == fgets(inputstring, STRLEN - 1, stdin))
                {
                    gmx_fatal(FARGS, "Error reading user input");
                }
                if (inputstring[0] != 'y' && inputstring[0] != 'Y')
                {
                    fprintf(stderr, "Will not convert\n");
                    exit(0);
                }
                nresav = fr->nre;
            }
            do_enx(in, fr);
            readtime[f] = fr->t;
            close_enx(in);
        }
        fprintf(stderr, "\n");
        free_enxnms(nre, enm);
    }

    free_enxframe(fr);
    sfree(fr);

    return nremin;
}


static void edit_files(gmx::ArrayRef<std::string> files,
                       real*                      readtime,
                       real*                      settime,
                       int*                       cont_type,
                       gmx_bool                   bSetTime,
                       gmx_bool                   bSort)
{
    gmx_bool ok;
    char     inputstring[STRLEN], *chptr;

    if (bSetTime)
    {
        if (files.size() == 1)
        {
            fprintf(stderr, "\n\nEnter the new start time:\n\n");
        }
        else
        {
            fprintf(stderr,
                    "\n\nEnter the new start time for each file.\n"
                    "There are two special options, both disables sorting:\n\n"
                    "c (continue) - The start time is taken from the end\n"
                    "of the previous file. Use it when your continuation run\n"
                    "restarts with t=0 and there is no overlap.\n\n"
                    "l (last) - The time in this file will be changed the\n"
                    "same amount as in the previous. Use it when the time in the\n"
                    "new run continues from the end of the previous one,\n"
                    "since this takes possible overlap into account.\n\n");
        }

        fprintf(stderr,
                "          File             Current start       New start\n"
                "---------------------------------------------------------\n");

        for (gmx::Index i = 0; i < files.ssize(); i++)
        {
            fprintf(stderr, "%25s   %10.3f             ", files[i].c_str(), readtime[i]);
            ok = FALSE;
            do
            {
                if (nullptr == fgets(inputstring, STRLEN - 1, stdin))
                {
                    gmx_fatal(FARGS, "Error reading user input");
                }
                inputstring[std::strlen(inputstring) - 1] = 0;

                if (inputstring[0] == 'c' || inputstring[0] == 'C')
                {
                    cont_type[i] = TIME_CONTINUE;
                    bSort        = FALSE;
                    ok           = TRUE;
                    settime[i]   = FLT_MAX;
                }
                else if (inputstring[0] == 'l' || inputstring[0] == 'L')
                {
                    cont_type[i] = TIME_LAST;
                    bSort        = FALSE;
                    ok           = TRUE;
                    settime[i]   = FLT_MAX;
                }
                else
                {
                    settime[i] = strtod(inputstring, &chptr);
                    if (chptr == inputstring)
                    {
                        fprintf(stderr, "Try that again: ");
                    }
                    else
                    {
                        cont_type[i] = TIME_EXPLICIT;
                        ok           = TRUE;
                    }
                }
            } while (!ok);
        }
        if (cont_type[0] != TIME_EXPLICIT)
        {
            cont_type[0] = TIME_EXPLICIT;
            settime[0]   = 0;
        }
    }
    else
    {
        for (gmx::Index i = 0; i < files.ssize(); i++)
        {
            settime[i] = readtime[i];
        }
    }

    if (bSort && files.size() > 1)
    {
        sort_files(files, settime);
    }
    else
    {
        fprintf(stderr, "Sorting disabled.\n");
    }


    /* Write out the new order and start times */
    fprintf(stderr,
            "\nSummary of files and start times used:\n\n"
            "          File                Start time\n"
            "-----------------------------------------\n");
    for (gmx::Index i = 0; i < files.ssize(); i++)
    {
        switch (cont_type[i])
        {
            case TIME_EXPLICIT:
                fprintf(stderr, "%25s   %10.3f\n", files[i].c_str(), settime[i]);
                break;
            case TIME_CONTINUE:
                fprintf(stderr, "%25s        Continue from end of last file\n", files[i].c_str());
                break;
            case TIME_LAST:
                fprintf(stderr, "%25s        Change by same amount as last file\n", files[i].c_str());
                break;
        }
    }
    fprintf(stderr, "\n");

    settime[files.size()]   = FLT_MAX;
    cont_type[files.size()] = TIME_EXPLICIT;
    readtime[files.size()]  = FLT_MAX;
}


static void update_ee_sum(int         nre,
                          int64_t*    ee_sum_step,
                          int64_t*    ee_sum_nsteps,
                          int64_t*    ee_sum_nsum,
                          t_energy*   ee_sum,
                          t_enxframe* fr,
                          int         out_step)
{
    int64_t nsteps, nsum, fr_nsum;
    int     i;

    nsteps = *ee_sum_nsteps;
    nsum   = *ee_sum_nsum;

    fr_nsum = fr->nsum;
    if (fr_nsum == 0)
    {
        fr_nsum = 1;
    }

    if (nsteps == 0)
    {
        if (fr_nsum == 1)
        {
            for (i = 0; i < nre; i++)
            {
                ee_sum[i].esum = fr->ener[i].e;
                ee_sum[i].eav  = 0;
            }
        }
        else
        {
            for (i = 0; i < nre; i++)
            {
                ee_sum[i].esum = fr->ener[i].esum;
                ee_sum[i].eav  = fr->ener[i].eav;
            }
        }
        nsteps = fr->nsteps;
        nsum   = fr_nsum;
    }
    else if (out_step + *ee_sum_nsum - *ee_sum_step == nsteps + fr->nsteps)
    {
        if (fr_nsum == 1)
        {
            for (i = 0; i < nre; i++)
            {
                ee_sum[i].eav += gmx::square(ee_sum[i].esum / nsum
                                             - (ee_sum[i].esum + fr->ener[i].e) / (nsum + 1))
                                 * nsum * (nsum + 1);
                ee_sum[i].esum += fr->ener[i].e;
            }
        }
        else
        {
            for (i = 0; i < fr->nre; i++)
            {
                ee_sum[i].eav += fr->ener[i].eav
                                 + gmx::square(ee_sum[i].esum / nsum
                                               - (ee_sum[i].esum + fr->ener[i].esum) / (nsum + fr->nsum))
                                           * nsum * (nsum + fr->nsum) / static_cast<double>(fr->nsum);
                ee_sum[i].esum += fr->ener[i].esum;
            }
        }
        nsteps += fr->nsteps;
        nsum += fr_nsum;
    }
    else
    {
        if (fr->nsum != 0)
        {
            fprintf(stderr, "\nWARNING: missing energy sums at time %f\n", fr->t);
        }
        nsteps = 0;
        nsum   = 0;
    }

    *ee_sum_step   = out_step;
    *ee_sum_nsteps = nsteps;
    *ee_sum_nsum   = nsum;
}

int gmx_eneconv(int argc, char* argv[])
{
    const char* desc[] = {
        "With [IT]multiple files[it] specified for the [TT]-f[tt] option:[PAR]",
        "Concatenates several energy files in sorted order.",
        "In the case of double time frames, the one",
        "in the later file is used. By specifying [TT]-settime[tt] you will be",
        "asked for the start time of each file. The input files are taken",
        "from the command line,",
        "such that the command [TT]gmx eneconv -f *.edr -o fixed.edr[tt] should do",
        "the trick. [PAR]",
        "With [IT]one file[it] specified for [TT]-f[tt]:[PAR]",
        "Reads one energy file and writes another, applying the [TT]-dt[tt],",
        "[TT]-offset[tt], [TT]-t0[tt] and [TT]-settime[tt] options and",
        "converting to a different format if necessary (indicated by file",
        "extensions).[PAR]",
        "[TT]-settime[tt] is applied first, then [TT]-dt[tt]/[TT]-offset[tt]",
        "followed by [TT]-b[tt] and [TT]-e[tt] to select which frames to write."
    };
    const char* bugs[] = {
        "When combining trajectories the sigma and E^2 (necessary for statistics) are not "
        "updated correctly. Only the actual energy is correct. One thus has to compute "
        "statistics in another way."
    };
    ener_file_t  in = nullptr, out = nullptr;
    gmx_enxnm_t* enm = nullptr;
#if 0
    ener_file_t       in, out = NULL;
    gmx_enxnm_t      *enm = NULL;
#endif
    t_enxframe *      fr, *fro;
    int64_t           ee_sum_step = 0, ee_sum_nsteps, ee_sum_nsum;
    t_energy*         ee_sum;
    int64_t           lastfilestep, laststep, startstep_file = 0;
    int               noutfr;
    int               nre, nremax, this_nre, i, kkk, nset, *set = nullptr;
    double            last_t;
    real *            readtime, *settime, timestep, tadjust;
    char              buf[22], buf2[22];
    int*              cont_type;
    gmx_bool          bNewFile, bFirst, bNewOutput;
    gmx_output_env_t* oenv;
    gmx_bool          warned_about_dh = FALSE;
    t_enxblock*       blocks          = nullptr;
    int               nblocks         = 0;
    int               nblocks_alloc   = 0;

    t_filenm fnm[] = {
        { efEDR, "-f", nullptr, ffRDMULT },
        { efEDR, "-o", "fixed", ffWRITE },
    };

#define NFILE asize(fnm)
    gmx_bool        bWrite;
    static real     delta_t = 0.0, toffset = 0, scalefac = 1;
    static gmx_bool bSetTime = FALSE;
    static gmx_bool bSort = TRUE, bError = TRUE;
    static real     begin     = -1;
    static real     end       = -1;
    gmx_bool        remove_dh = FALSE;

    t_pargs pa[] = {
        { "-b", FALSE, etREAL, { &begin }, "First time to use" },
        { "-e", FALSE, etREAL, { &end }, "Last time to use" },
        { "-dt", FALSE, etREAL, { &delta_t }, "Only write out frame when t MOD dt = offset" },
        { "-offset", FALSE, etREAL, { &toffset }, "Time offset for [TT]-dt[tt] option" },
        { "-settime", FALSE, etBOOL, { &bSetTime }, "Change starting time interactively" },
        { "-sort", FALSE, etBOOL, { &bSort }, "Sort energy files (not frames)" },
        { "-rmdh", FALSE, etBOOL, { &remove_dh }, "Remove free energy block data" },
        { "-scalefac", FALSE, etREAL, { &scalefac }, "Multiply energy component by this factor" },
        { "-error", FALSE, etBOOL, { &bError }, "Stop on errors in the file" }
    };

    if (!parse_common_args(
                &argc, argv, 0, NFILE, fnm, asize(pa), pa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }
    fprintf(stdout,
            "Note that major changes are planned in future for "
            "eneconv, to improve usability and utility.");

    tadjust      = 0;
    nremax       = 0;
    nset         = 0;
    timestep     = 0.0;
    lastfilestep = 0;
    laststep     = 0;

    auto files = gmx::copyOf(opt2fns("-f", NFILE, fnm));

    if (files.empty())
    {
        gmx_fatal(FARGS, "No input files!");
    }

    snew(settime, files.size() + 1);
    snew(readtime, files.size() + 1);
    snew(cont_type, files.size() + 1);

    nre = scan_ene_files(files, readtime, &timestep, &nremax);
    edit_files(files, readtime, settime, cont_type, bSetTime, bSort);

    ee_sum_nsteps = 0;
    ee_sum_nsum   = 0;
    snew(ee_sum, nremax);

    snew(fr, 1);
    snew(fro, 1);
    fro->t   = -1e20;
    fro->nre = nre;
    snew(fro->ener, nremax);

    noutfr = 0;
    bFirst = TRUE;

    last_t = fro->t;
    for (size_t f = 0; f < files.size(); f++)
    {
        bNewFile   = TRUE;
        bNewOutput = TRUE;
        in         = open_enx(files[f].c_str(), "r");
        enm        = nullptr;
        do_enxnms(in, &this_nre, &enm);
        if (f == 0)
        {
            if (scalefac != 1)
            {
                set = select_it(nre, enm, &nset);
            }

            /* write names to the output file */
            out = open_enx(opt2fn("-o", NFILE, fnm), "w");
            do_enxnms(out, &nre, &enm);
        }

        /* start reading from the next file */
        while ((fro->t <= (settime[f + 1] + GMX_REAL_EPS)) && do_enx(in, fr))
        {
            if (bNewFile)
            {
                startstep_file = fr->step;
                tadjust        = settime[f] - fr->t;
                if (cont_type[f + 1] == TIME_LAST)
                {
                    settime[f + 1]   = readtime[f + 1] - readtime[f] + settime[f];
                    cont_type[f + 1] = TIME_EXPLICIT;
                }
                bNewFile = FALSE;
            }

            if (tadjust + fr->t <= last_t)
            {
                /* Skip this frame, since we already have it / past it */
                if (debug)
                {
                    fprintf(debug, "fr->step %s, fr->t %.4f\n", gmx_step_str(fr->step, buf), fr->t);
                    fprintf(debug, "tadjust %12.6e + fr->t %12.6e <= t %12.6e\n", tadjust, fr->t, last_t);
                }
                continue;
            }

            fro->step = lastfilestep + fr->step - startstep_file;
            fro->t    = tadjust + fr->t;

            bWrite = ((begin < 0 || (fro->t >= begin - GMX_REAL_EPS))
                      && (end < 0 || (fro->t <= end + GMX_REAL_EPS))
                      && (fro->t <= settime[f + 1] + 0.5 * timestep));

            if (debug)
            {
                fprintf(debug,
                        "fr->step %s, fr->t %.4f, fro->step %s fro->t %.4f, w %s\n",
                        gmx_step_str(fr->step, buf),
                        fr->t,
                        gmx_step_str(fro->step, buf2),
                        fro->t,
                        gmx::boolToString(bWrite));
            }

            if (bError)
            {
                if ((end > 0) && (fro->t > end + GMX_REAL_EPS))
                {
                    f = files.size();
                    break;
                }
            }

            if (fro->t >= begin - GMX_REAL_EPS)
            {
                if (bFirst)
                {
                    bFirst = FALSE;
                }
                if (bWrite)
                {
                    update_ee_sum(nre, &ee_sum_step, &ee_sum_nsteps, &ee_sum_nsum, ee_sum, fr, fro->step);
                }
            }

            /* determine if we should write it */
            if (bWrite && (delta_t == 0 || bRmod(fro->t, toffset, delta_t)))
            {
                laststep = fro->step;
                last_t   = fro->t;
                if (bNewOutput)
                {
                    bNewOutput = FALSE;
                    fprintf(stderr,
                            "\nContinue writing frames from t=%g, step=%s\n",
                            fro->t,
                            gmx_step_str(fro->step, buf));
                }

                /* Copy the energies */
                for (i = 0; i < nre; i++)
                {
                    fro->ener[i].e = fr->ener[i].e;
                }

                fro->nsteps = ee_sum_nsteps;
                fro->dt     = fr->dt;

                if (ee_sum_nsum <= 1)
                {
                    fro->nsum = 0;
                }
                else
                {
                    fro->nsum = int64_to_int(ee_sum_nsum, "energy average summation");
                    /* Copy the energy sums */
                    for (i = 0; i < nre; i++)
                    {
                        fro->ener[i].esum = ee_sum[i].esum;
                        fro->ener[i].eav  = ee_sum[i].eav;
                    }
                }
                /* We wrote the energies, so reset the counts */
                ee_sum_nsteps = 0;
                ee_sum_nsum   = 0;

                if (scalefac != 1)
                {
                    for (kkk = 0; kkk < nset; kkk++)
                    {
                        fro->ener[set[kkk]].e *= scalefac;
                        if (fro->nsum > 0)
                        {
                            fro->ener[set[kkk]].eav *= scalefac * scalefac;
                            fro->ener[set[kkk]].esum *= scalefac;
                        }
                    }
                }
                /* Copy restraint stuff */
                /*fro->ndisre       = fr->ndisre;
                   fro->disre_rm3tav = fr->disre_rm3tav;
                   fro->disre_rt     = fr->disre_rt;*/
                fro->nblock = fr->nblock;
                /*fro->nr           = fr->nr;*/
                fro->block = fr->block;

                /* check if we have blocks with delta_h data and are throwing
                   away data */
                if (fro->nblock > 0)
                {
                    if (remove_dh)
                    {
                        int i;
                        if (!blocks || nblocks_alloc < fr->nblock)
                        {
                            /* we pre-allocate the blocks */
                            nblocks_alloc = fr->nblock;
                            snew(blocks, nblocks_alloc);
                        }
                        nblocks = 0; /* number of blocks so far */

                        for (i = 0; i < fr->nblock; i++)
                        {
                            if ((fr->block[i].id != enxDHCOLL) && (fr->block[i].id != enxDH)
                                && (fr->block[i].id != enxDHHIST))
                            {
                                /* copy everything verbatim */
                                blocks[nblocks] = fr->block[i];
                                nblocks++;
                            }
                        }
                        /* now set the block pointer to the new blocks */
                        fro->nblock = nblocks;
                        fro->block  = blocks;
                    }
                    else if (delta_t > 0)
                    {
                        if (!warned_about_dh)
                        {
                            for (i = 0; i < fr->nblock; i++)
                            {
                                if (fr->block[i].id == enxDH || fr->block[i].id == enxDHHIST)
                                {
                                    int size;
                                    if (fr->block[i].id == enxDH)
                                    {
                                        size = fr->block[i].sub[2].nr;
                                    }
                                    else
                                    {
                                        size = fr->nsteps;
                                    }
                                    if (size > 0)
                                    {
                                        printf("\nWARNING: %s contains delta H blocks or "
                                               "histograms for which\n"
                                               "         some data is thrown away on a "
                                               "block-by-block basis, where each block\n"
                                               "         contains up to %d samples.\n"
                                               "         This is almost certainly not what you "
                                               "want.\n"
                                               "         Use the -rmdh option to throw all delta H "
                                               "samples away.\n"
                                               "         Use gmx energy -odh option to extract "
                                               "these samples.\n",
                                               files[f].c_str(),
                                               size);
                                        warned_about_dh = TRUE;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }

                do_enx(out, fro);
                if (noutfr % 1000 == 0)
                {
                    fprintf(stderr, "Writing frame time %g    ", fro->t);
                }
                noutfr++;
            }
        }
        if (f == files.size())
        {
            f--;
        }
        printf("\nLast step written from %s: t %g, step %s\n",
               files[f].c_str(),
               last_t,
               gmx_step_str(laststep, buf));
        lastfilestep = laststep;

        /* set the next time from the last in previous file */
        if (cont_type[f + 1] == TIME_CONTINUE)
        {
            settime[f + 1] = fro->t;
            /* in this case we have already written the last frame of
             * previous file, so update begin to avoid doubling it
             * with the start of the next file
             */
            begin = fro->t + 0.5 * timestep;
            /* cont_type[f+1]==TIME_EXPLICIT; */
        }

        if ((fro->t < end) && (f < files.size() - 1) && (fro->t < settime[f + 1] - 1.5 * timestep))
        {
            fprintf(stderr, "\nWARNING: There might be a gap around t=%g\n", fro->t);
        }

        /* move energies to lastee */
        close_enx(in);
        free_enxnms(this_nre, enm);

        fprintf(stderr, "\n");
    }
    if (noutfr == 0)
    {
        fprintf(stderr, "No frames written.\n");
    }
    else
    {
        fprintf(stderr, "Last frame written was at step %s, time %f\n", gmx_step_str(fro->step, buf), fro->t);
        fprintf(stderr, "Wrote %d frames\n", noutfr);
    }

    return 0;
}

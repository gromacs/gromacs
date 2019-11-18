/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include "trjcat.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <string>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/tngio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/index.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#define TIME_EXPLICIT 0
#define TIME_CONTINUE 1
#define TIME_LAST 2
#ifndef FLT_MAX
#    define FLT_MAX 1e36
#endif
#define FLAGS (TRX_READ_X | TRX_READ_V | TRX_READ_F)

static void scan_trj_files(gmx::ArrayRef<const std::string> files,
                           real*                            readtime,
                           real*                            timestep,
                           int                              imax,
                           const gmx_output_env_t*          oenv)
{
    /* Check start time of all files */
    int          natoms = 0;
    t_trxstatus* status;
    t_trxframe   fr;
    bool         ok;

    for (gmx::index i = 0; i < files.ssize(); i++)
    {
        ok = read_first_frame(oenv, &status, files[i].c_str(), &fr, FLAGS);

        if (!ok)
        {
            gmx_fatal(FARGS, "\nCouldn't read frame from file.");
        }
        if (fr.bTime)
        {
            readtime[i] = fr.time;
        }
        else
        {
            readtime[i] = 0;
            fprintf(stderr, "\nWARNING: Couldn't find a time in the frame.\n");
        }

        if (i == 0)
        {
            natoms = fr.natoms;
        }
        else
        {
            if (imax == -1)
            {
                if (natoms != fr.natoms)
                {
                    gmx_fatal(FARGS, "\nDifferent numbers of atoms (%d/%d) in files", natoms, fr.natoms);
                }
            }
            else
            {
                if (fr.natoms <= imax)
                {
                    gmx_fatal(FARGS, "\nNot enough atoms (%d) for index group (%d)", fr.natoms, imax);
                }
            }
        }
        ok = read_next_frame(oenv, status, &fr);
        if (ok && fr.bTime)
        {
            timestep[i] = fr.time - readtime[i];
        }
        else
        {
            timestep[i] = 0;
        }

        close_trx(status);
        if (fr.bX)
        {
            sfree(fr.x);
        }
        if (fr.bV)
        {
            sfree(fr.v);
        }
        if (fr.bF)
        {
            sfree(fr.f);
        }
    }
    fprintf(stderr, "\n");
}

static void sort_files(gmx::ArrayRef<std::string> files, real* settime)
{
    for (gmx::index i = 0; i < files.ssize(); i++)
    {
        gmx::index minidx = i;
        for (gmx::index j = i + 1; j < files.ssize(); j++)
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
            std::swap(files[i], files[minidx]);
        }
    }
}

static void edit_files(gmx::ArrayRef<std::string> files,
                       real*                      readtime,
                       real*                      timestep,
                       real*                      settime,
                       int*                       cont_type,
                       gmx_bool                   bSetTime,
                       gmx_bool                   bSort,
                       const gmx_output_env_t*    oenv)
{
    gmx_bool ok;
    char     inputstring[STRLEN], *chptr;

    auto timeUnit = output_env_get_time_unit(oenv);
    if (bSetTime)
    {
        fprintf(stderr,
                "\n\nEnter the new start time (%s) for each file.\n"
                "There are two special options, both disable sorting:\n\n"
                "c (continue) - The start time is taken from the end\n"
                "of the previous file. Use it when your continuation run\n"
                "restarts with t=0.\n\n"
                "l (last) - The time in this file will be changed the\n"
                "same amount as in the previous. Use it when the time in the\n"
                "new run continues from the end of the previous one,\n"
                "since this takes possible overlap into account.\n\n",
                timeUnit.c_str());

        fprintf(stderr,
                "          File             Current start (%s)  New start (%s)\n"
                "---------------------------------------------------------\n",
                timeUnit.c_str(), timeUnit.c_str());

        for (gmx::index i = 0; i < files.ssize(); i++)
        {
            fprintf(stderr, "%25s   %10.3f %s          ", files[i].c_str(),
                    output_env_conv_time(oenv, readtime[i]), timeUnit.c_str());
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
                    settime[i] = strtod(inputstring, &chptr) * output_env_get_time_invfactor(oenv);
                    if (chptr == inputstring)
                    {
                        fprintf(stderr,
                                "'%s' not recognized as a floating point number, 'c' or 'l'. "
                                "Try again: ",
                                inputstring);
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
        for (gmx::index i = 0; i < files.ssize(); i++)
        {
            settime[i] = readtime[i];
        }
    }
    if (!bSort)
    {
        fprintf(stderr, "Sorting disabled.\n");
    }
    else
    {
        sort_files(files, settime);
    }
    /* Write out the new order and start times */
    fprintf(stderr,
            "\nSummary of files and start times used:\n\n"
            "          File                Start time       Time step\n"
            "---------------------------------------------------------\n");
    for (gmx::index i = 0; i < files.ssize(); i++)
    {
        switch (cont_type[i])
        {
            case TIME_EXPLICIT:
                fprintf(stderr, "%25s   %10.3f %s   %10.3f %s", files[i].c_str(),
                        output_env_conv_time(oenv, settime[i]), timeUnit.c_str(),
                        output_env_conv_time(oenv, timestep[i]), timeUnit.c_str());
                if (i > 0 && cont_type[i - 1] == TIME_EXPLICIT && settime[i] == settime[i - 1])
                {
                    fprintf(stderr, " WARNING: same Start time as previous");
                }
                fprintf(stderr, "\n");
                break;
            case TIME_CONTINUE:
                fprintf(stderr, "%25s        Continue from last file\n", files[i].c_str());
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

static void do_demux(gmx::ArrayRef<const std::string> inFiles,
                     gmx::ArrayRef<const std::string> outFiles,
                     int                              nval,
                     real**                           value,
                     real*                            time,
                     real                             dt_remd,
                     int                              isize,
                     int                              index[],
                     real                             dt,
                     const gmx_output_env_t*          oenv)
{
    int           k, natoms;
    t_trxstatus **fp_in, **fp_out;
    gmx_bool      bCont, *bSet;
    real          t, first_time = 0;
    t_trxframe*   trx;

    snew(fp_in, inFiles.size());
    snew(trx, inFiles.size());
    snew(bSet, inFiles.size());
    natoms = -1;
    t      = -1;
    for (gmx::index i = 0; i < inFiles.ssize(); i++)
    {
        read_first_frame(oenv, &(fp_in[i]), inFiles[i].c_str(), &(trx[i]), TRX_NEED_X);
        if (natoms == -1)
        {
            natoms     = trx[i].natoms;
            first_time = trx[i].time;
        }
        else if (natoms != trx[i].natoms)
        {
            gmx_fatal(FARGS, "Trajectory file %s has %d atoms while previous trajs had %d atoms",
                      inFiles[i].c_str(), trx[i].natoms, natoms);
        }
        if (t == -1)
        {
            t = trx[i].time;
        }
        else if (t != trx[i].time)
        {
            gmx_fatal(FARGS, "Trajectory file %s has time %f while previous trajs had time %f",
                      inFiles[i].c_str(), trx[i].time, t);
        }
    }

    snew(fp_out, inFiles.size());
    for (gmx::index i = 0; i < inFiles.ssize(); i++)
    {
        fp_out[i] = open_trx(outFiles[i].c_str(), "w");
    }
    k = 0;
    if (std::round(time[k] - t) != 0)
    {
        gmx_fatal(FARGS, "First time in demuxing table does not match trajectories");
    }
    do
    {
        while ((k + 1 < nval) && ((trx[0].time - time[k + 1]) > dt_remd * 0.1))
        {
            k++;
        }
        if (debug)
        {
            fprintf(debug, "trx[0].time = %g, time[k] = %g\n", trx[0].time, time[k]);
        }
        for (gmx::index i = 0; i < inFiles.ssize(); i++)
        {
            bSet[i] = FALSE;
        }
        for (gmx::index i = 0; i < inFiles.ssize(); i++)
        {
            int j = gmx::roundToInt(value[i][k]);
            range_check(j, 0, inFiles.size());
            if (bSet[j])
            {
                gmx_fatal(FARGS, "Demuxing the same replica %d twice at time %f", j, trx[0].time);
            }
            bSet[j] = TRUE;

            if (dt == 0 || bRmod(trx[i].time, first_time, dt))
            {
                if (index)
                {
                    write_trxframe_indexed(fp_out[j], &trx[i], isize, index, nullptr);
                }
                else
                {
                    write_trxframe(fp_out[j], &trx[i], nullptr);
                }
            }
        }

        bCont = (k < nval);
        for (gmx::index i = 0; i < inFiles.ssize(); i++)
        {
            bCont = bCont && read_next_frame(oenv, fp_in[i], &trx[i]);
        }
    } while (bCont);

    for (gmx::index i = 0; i < inFiles.ssize(); i++)
    {
        close_trx(fp_in[i]);
        close_trx(fp_out[i]);
    }
}

int gmx_trjcat(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] concatenates several input trajectory files in sorted order. ",
        "In case of double time frames the one in the later file is used. ",
        "By specifying [TT]-settime[tt] you will be asked for the start time ",
        "of each file. The input files are taken from the command line, ",
        "such that a command like [TT]gmx trjcat -f *.trr -o fixed.trr[tt] should do ",
        "the trick. Using [TT]-cat[tt], you can simply paste several files ",
        "together without removal of frames with identical time stamps.[PAR]",
        "One important option is inferred when the output file is amongst the",
        "input files. In that case that particular file will be appended to",
        "which implies you do not need to store double the amount of data.",
        "Obviously the file to append to has to be the one with lowest starting",
        "time since one can only append at the end of a file.[PAR]",
        "If the [TT]-demux[tt] option is given, the N trajectories that are",
        "read, are written in another order as specified in the [REF].xvg[ref] file.",
        "The [REF].xvg[ref] file should contain something like::",
        "",
        "    0  0  1  2  3  4  5",
        "    2  1  0  2  3  5  4",
        "",
        "The first number is the time, and subsequent numbers point to",
        "trajectory indices.",
        "The frames corresponding to the numbers present at the first line",
        "are collected into the output trajectory. If the number of frames in",
        "the trajectory does not match that in the [REF].xvg[ref] file then the program",
        "tries to be smart. Beware."
    };
    static gmx_bool bCat            = FALSE;
    static gmx_bool bSort           = TRUE;
    static gmx_bool bKeepLast       = FALSE;
    static gmx_bool bKeepLastAppend = FALSE;
    static gmx_bool bOverwrite      = FALSE;
    static gmx_bool bSetTime        = FALSE;
    static gmx_bool bDeMux;
    static real     begin = -1;
    static real     end   = -1;
    static real     dt    = 0;

    t_pargs pa[] = {
        { "-b", FALSE, etTIME, { &begin }, "First time to use (%t)" },
        { "-e", FALSE, etTIME, { &end }, "Last time to use (%t)" },
        { "-dt", FALSE, etTIME, { &dt }, "Only write frame when t MOD dt = first time (%t)" },
        { "-settime", FALSE, etBOOL, { &bSetTime }, "Change starting time interactively" },
        { "-sort", FALSE, etBOOL, { &bSort }, "Sort trajectory files (not frames)" },
        { "-keeplast",
          FALSE,
          etBOOL,
          { &bKeepLast },
          "Keep overlapping frames at end of trajectory" },
        { "-overwrite",
          FALSE,
          etBOOL,
          { &bOverwrite },
          "Overwrite overlapping frames during appending" },
        { "-cat", FALSE, etBOOL, { &bCat }, "Do not discard double time frames" }
    };
#define npargs asize(pa)
    int               ftpin, i, frame, frame_out;
    t_trxstatus *     status, *trxout = nullptr;
    real              t_corr;
    t_trxframe        fr, frout;
    int               n_append;
    gmx_bool          bNewFile, bIndex, bWrite;
    int*              cont_type;
    real *            readtime, *timest, *settime;
    real              first_time = 0, lasttime = 0, last_ok_t = -1, timestep;
    gmx_bool          lastTimeSet = FALSE;
    real              last_frame_time, searchtime;
    int               isize = 0, j;
    int *             index = nullptr, imax;
    char*             grpname;
    real **           val = nullptr, *t = nullptr, dt_remd;
    int               n, nset, ftpout = -1, prevEndStep = 0, filetype;
    gmx_off_t         fpos;
    gmx_output_env_t* oenv;
    t_filenm          fnm[] = { { efTRX, "-f", nullptr, ffRDMULT },
                       { efTRO, "-o", nullptr, ffWRMULT },
                       { efNDX, "-n", "index", ffOPTRD },
                       { efXVG, "-demux", "remd", ffOPTRD } };

#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, PCA_TIME_UNIT, NFILE, fnm, asize(pa), pa, asize(desc), desc,
                           0, nullptr, &oenv))
    {
        return 0;
    }
    fprintf(stdout,
            "Note that major changes are planned in future for "
            "trjcat, to improve usability and utility.");

    auto timeUnit = output_env_get_time_unit(oenv);

    bIndex = ftp2bSet(efNDX, NFILE, fnm);
    bDeMux = ftp2bSet(efXVG, NFILE, fnm);
    bSort  = bSort && !bDeMux;

    imax = -1;
    if (bIndex)
    {
        printf("Select group for output\n");
        rd_index(ftp2fn(efNDX, NFILE, fnm), 1, &isize, &index, &grpname);
        /* scan index */
        imax = index[0];
        for (i = 1; i < isize; i++)
        {
            imax = std::max(imax, index[i]);
        }
    }
    if (bDeMux)
    {
        nset    = 0;
        dt_remd = 0;
        val     = read_xvg_time(opt2fn("-demux", NFILE, fnm), TRUE, opt2parg_bSet("-b", npargs, pa),
                            begin, opt2parg_bSet("-e", npargs, pa), end, 1, &nset, &n, &dt_remd, &t);
        printf("Read %d sets of %d points, dt = %g\n\n", nset, n, dt_remd);
        if (debug)
        {
            fprintf(debug, "Dump of replica_index.xvg\n");
            for (i = 0; (i < n); i++)
            {
                fprintf(debug, "%10g", t[i]);
                for (j = 0; (j < nset); j++)
                {
                    fprintf(debug, "  %3d", static_cast<int>(std::round(val[j][i])));
                }
                fprintf(debug, "\n");
            }
        }
    }

    gmx::ArrayRef<const std::string> inFiles = opt2fns("-f", NFILE, fnm);
    if (inFiles.empty())
    {
        gmx_fatal(FARGS, "No input files!");
    }

    if (bDeMux && ssize(inFiles) != nset)
    {
        gmx_fatal(FARGS, "You have specified %td files and %d entries in the demux table",
                  inFiles.ssize(), nset);
    }

    ftpin = fn2ftp(inFiles[0].c_str());

    if (ftpin != efTRR && ftpin != efXTC && ftpin != efTNG)
    {
        gmx_fatal(FARGS, "gmx trjcat can only handle binary trajectory formats (trr, xtc, tng)");
    }

    for (const std::string& inFile : inFiles)
    {
        if (ftpin != fn2ftp(inFile.c_str()))
        {
            gmx_fatal(FARGS, "All input files must be of the same (trr, xtc or tng) format");
        }
    }

    gmx::ArrayRef<const std::string> outFiles = opt2fns("-o", NFILE, fnm);
    if (outFiles.empty())
    {
        gmx_fatal(FARGS, "No output files!");
    }
    if ((outFiles.size() > 1) && !bDeMux)
    {
        gmx_fatal(FARGS,
                  "Don't know what to do with more than 1 output file if  not demultiplexing");
    }
    else if (bDeMux && ssize(outFiles) != nset && outFiles.size() != 1)
    {
        gmx_fatal(FARGS, "Number of output files should be 1 or %d (#input files), not %td", nset,
                  outFiles.ssize());
    }
    if (bDeMux)
    {
        auto outFilesDemux = gmx::copyOf(outFiles);
        if (gmx::ssize(outFilesDemux) != nset)
        {
            std::string name = outFilesDemux[0];
            outFilesDemux.resize(nset);
            for (i = 0; (i < nset); i++)
            {
                outFilesDemux[0] = gmx::formatString("%d_%s", i, name.c_str());
            }
        }
        do_demux(inFiles, outFilesDemux, n, val, t, dt_remd, isize, index, dt, oenv);
    }
    else
    {
        snew(readtime, inFiles.size() + 1);
        snew(timest, inFiles.size() + 1);
        scan_trj_files(inFiles, readtime, timest, imax, oenv);

        snew(settime, inFiles.size() + 1);
        snew(cont_type, inFiles.size() + 1);
        auto inFilesEdited = gmx::copyOf(inFiles);
        edit_files(inFilesEdited, readtime, timest, settime, cont_type, bSetTime, bSort, oenv);

        /* Check whether the output file is amongst the input files
         * This has to be done after sorting etc.
         */
        const char* out_file = outFiles[0].c_str();
        ftpout               = fn2ftp(out_file);
        n_append             = -1;
        for (size_t i = 0; i < inFilesEdited.size() && n_append == -1; i++)
        {
            if (std::strcmp(inFilesEdited[i].c_str(), out_file) == 0)
            {
                n_append = i;
            }
        }
        if (n_append == 0)
        {
            fprintf(stderr, "Will append to %s rather than creating a new file\n", out_file);
        }
        else if (n_append != -1)
        {
            gmx_fatal(FARGS, "Can only append to the first file which is %s (not %s)",
                      inFilesEdited[0].c_str(), out_file);
        }

        /* Not checking input format, could be dangerous :-) */
        /* Not checking output format, equally dangerous :-) */

        frame     = -1;
        frame_out = -1;
        /* the default is not to change the time at all,
         * but this is overridden by the edit_files routine
         */
        t_corr = 0;

        if (n_append == -1)
        {
            if (ftpout == efTNG)
            {
                if (ftpout != ftpin)
                {
                    gmx_fatal(FARGS, "When writing TNG the input file format must also be TNG");
                }
                if (bIndex)
                {
                    trxout = trjtools_gmx_prepare_tng_writing(
                            out_file, 'w', nullptr, inFilesEdited[0].c_str(), isize, nullptr,
                            gmx::arrayRefFromArray(index, isize), grpname);
                }
                else
                {
                    trxout = trjtools_gmx_prepare_tng_writing(
                            out_file, 'w', nullptr, inFilesEdited[0].c_str(), -1, nullptr, {}, nullptr);
                }
            }
            else
            {
                trxout = open_trx(out_file, "w");
            }
            std::memset(&frout, 0, sizeof(frout));
        }
        else
        {
            t_fileio* stfio;

            if (!read_first_frame(oenv, &status, out_file, &fr, FLAGS))
            {
                gmx_fatal(FARGS, "Reading first frame from %s", out_file);
            }

            stfio = trx_get_fileio(status);
            if (!bKeepLast && !bOverwrite)
            {
                fprintf(stderr,
                        "\n\nWARNING: Appending without -overwrite implies -keeplast "
                        "between the first two files. \n"
                        "If the trajectories have an overlap and have not been written binary \n"
                        "reproducible this will produce an incorrect trajectory!\n\n");

                filetype = gmx_fio_getftp(stfio);
                /* Fails if last frame is incomplete
                 * We can't do anything about it without overwriting
                 * */
                if (filetype == efXTC || filetype == efTNG)
                {
                    lasttime = trx_get_time_of_final_frame(status);
                    fr.time  = lasttime;
                }
                else
                {
                    while (read_next_frame(oenv, status, &fr)) {}
                    lasttime = fr.time;
                }
                lastTimeSet     = TRUE;
                bKeepLastAppend = TRUE;
                close_trx(status);
                trxout = open_trx(out_file, "a");
            }
            else if (bOverwrite)
            {
                if (gmx_fio_getftp(stfio) != efXTC)
                {
                    gmx_fatal(FARGS, "Overwrite only supported for XTC.");
                }
                last_frame_time = trx_get_time_of_final_frame(status);

                /* xtc_seek_time broken for trajectories containing only 1 or 2 frames
                 *     or when seek time = 0 */
                if (inFilesEdited.size() > 1 && settime[1] < last_frame_time + timest[0] * 0.5)
                {
                    /* Jump to one time-frame before the start of next
                     *  trajectory file */
                    searchtime = settime[1] - timest[0] * 1.25;
                }
                else
                {
                    searchtime = last_frame_time;
                }
                if (xtc_seek_time(stfio, searchtime, fr.natoms, TRUE))
                {
                    gmx_fatal(FARGS, "Error seeking to append position.");
                }
                read_next_frame(oenv, status, &fr);
                if (std::abs(searchtime - fr.time) > timest[0] * 0.5)
                {
                    gmx_fatal(FARGS, "Error seeking: attempted to seek to %f but got %f.",
                              searchtime, fr.time);
                }
                lasttime    = fr.time;
                lastTimeSet = TRUE;
                fpos        = gmx_fio_ftell(stfio);
                close_trx(status);
                trxout = open_trx(out_file, "r+");
                if (gmx_fio_seek(trx_get_fileio(trxout), fpos))
                {
                    gmx_fatal(FARGS, "Error seeking to append position.");
                }
            }
            if (lastTimeSet)
            {
                printf("\n Will append after %f \n", lasttime);
            }
            frout = fr;
        }
        /* Lets stitch up some files */
        timestep = timest[0];
        for (size_t i = n_append + 1; i < inFilesEdited.size(); i++)
        {
            /* Open next file */

            /* set the next time from the last frame in previous file */
            if (i > 0)
            {
                /* When writing TNG the step determine which frame to write. Use an
                 * offset to be able to increase steps properly when changing files. */
                if (ftpout == efTNG)
                {
                    prevEndStep = frout.step;
                }

                if (frame_out >= 0)
                {
                    if (cont_type[i] == TIME_CONTINUE)
                    {
                        begin = frout.time;
                        begin += 0.5 * timestep;
                        settime[i]   = frout.time;
                        cont_type[i] = TIME_EXPLICIT;
                    }
                    else if (cont_type[i] == TIME_LAST)
                    {
                        begin = frout.time;
                        begin += 0.5 * timestep;
                    }
                    /* Or, if the time in the next part should be changed by the
                     * same amount, start at half a timestep from the last time
                     * so we dont repeat frames.
                     */
                    /* I don't understand the comment above, but for all the cases
                     * I tried the code seems to work properly. B. Hess 2008-4-2.
                     */
                }
                /* Or, if time is set explicitly, we check for overlap/gap */
                if (cont_type[i] == TIME_EXPLICIT)
                {
                    if (i < inFilesEdited.size() && frout.time < settime[i] - 1.5 * timestep)
                    {
                        fprintf(stderr,
                                "WARNING: Frames around t=%f %s have a different "
                                "spacing than the rest,\n"
                                "might be a gap or overlap that couldn't be corrected "
                                "automatically.\n",
                                output_env_conv_time(oenv, frout.time), timeUnit.c_str());
                    }
                }
            }

            /* if we don't have a timestep in the current file, use the old one */
            if (timest[i] != 0)
            {
                timestep = timest[i];
            }
            read_first_frame(oenv, &status, inFilesEdited[i].c_str(), &fr, FLAGS);
            if (!fr.bTime)
            {
                fr.time = 0;
                fprintf(stderr, "\nWARNING: Couldn't find a time in the frame.\n");
            }

            if (cont_type[i] == TIME_EXPLICIT)
            {
                t_corr = settime[i] - fr.time;
            }
            /* t_corr is the amount we want to change the time.
             * If the user has chosen not to change the time for
             * this part of the trajectory t_corr remains at
             * the value it had in the last part, changing this
             * by the same amount.
             * If no value was given for the first trajectory part
             * we let the time start at zero, see the edit_files routine.
             */

            bNewFile = TRUE;

            if (!lastTimeSet)
            {
                lasttime    = 0;
                lastTimeSet = true;
            }
            printf("\n");
            printf("lasttime %g\n", lasttime);

            do
            {
                /* copy the input frame to the output frame */
                frout = fr;
                /* set the new time by adding the correct calculated above */
                frout.time += t_corr;
                if (ftpout == efTNG)
                {
                    frout.step += prevEndStep;
                }
                /* quit if we have reached the end of what should be written */
                if ((end > 0) && (frout.time > end + GMX_REAL_EPS))
                {
                    i = inFilesEdited.size();
                    break;
                }

                /* determine if we should write this frame (dt is handled elsewhere) */
                if (bCat) /* write all frames of all files */
                {
                    bWrite = TRUE;
                }
                else if (bKeepLast || (bKeepLastAppend && i == 1))
                /* write till last frame of this traj
                   and skip first frame(s) of next traj */
                {
                    bWrite = (frout.time > lasttime + 0.5 * timestep);
                }
                else /* write till first frame of next traj */
                {
                    bWrite = (frout.time < settime[i + 1] - 0.5 * timestep);
                }

                if (bWrite && (frout.time >= begin))
                {
                    frame++;
                    if (frame_out == -1)
                    {
                        first_time = frout.time;
                    }
                    lasttime    = frout.time;
                    lastTimeSet = TRUE;
                    if (dt == 0 || bRmod(frout.time, first_time, dt))
                    {
                        frame_out++;
                        last_ok_t = frout.time;
                        if (bNewFile)
                        {
                            fprintf(stderr,
                                    "\nContinue writing frames from %s t=%g %s, "
                                    "frame=%d      \n",
                                    inFilesEdited[i].c_str(),
                                    output_env_conv_time(oenv, frout.time), timeUnit.c_str(), frame);
                            bNewFile = FALSE;
                        }

                        if (bIndex)
                        {
                            write_trxframe_indexed(trxout, &frout, isize, index, nullptr);
                        }
                        else
                        {
                            write_trxframe(trxout, &frout, nullptr);
                        }
                        if (((frame % 10) == 0) || (frame < 10))
                        {
                            fprintf(stderr, " ->  frame %6d time %8.3f %s     \r", frame_out,
                                    output_env_conv_time(oenv, frout.time), timeUnit.c_str());
                            fflush(stderr);
                        }
                    }
                }
            } while (read_next_frame(oenv, status, &fr));

            close_trx(status);
        }
        if (trxout)
        {
            close_trx(trxout);
        }
        fprintf(stderr, "\nLast frame written was %d, time %f %s\n", frame,
                output_env_conv_time(oenv, last_ok_t), timeUnit.c_str());
    }

    return 0;
}

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/g96io.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/groio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/tngio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

enum {
    euSel, euRect, euTric, euCompact, euNR
};


static void mk_filenm(char *base, const char *ext, int ndigit, int file_nr,
                      char out_file[])
{
    char nbuf[128];
    int  nd = 0, fnr;

    std::strcpy(out_file, base);
    fnr = file_nr;
    do
    {
        fnr /= 10;
        nd++;
    }
    while (fnr > 0);

    if (nd < ndigit)
    {
        std::strncat(out_file, "00000000000", ndigit-nd);
    }
    sprintf(nbuf, "%d.", file_nr);
    std::strcat(out_file, nbuf);
    std::strcat(out_file, ext);
}

static void check_trr(const char *fn)
{
    if (fn2ftp(fn) != efTRR)
    {
        gmx_fatal(FARGS, "%s is not a trajectory file, exiting\n", fn);
    }
}

static void do_trunc(const char *fn, real t0)
{
    t_fileio        *in;
    FILE            *fp;
    gmx_bool         bStop, bOK;
    gmx_trr_header_t sh;
    gmx_off_t        fpos;
    char             yesno[256];
    int              j;
    real             t = 0;

    if (t0 == -1)
    {
        gmx_fatal(FARGS, "You forgot to set the truncation time");
    }

    /* Check whether this is a .trr file */
    check_trr(fn);

    in   = gmx_trr_open(fn, "r");
    fp   = gmx_fio_getfp(in);
    if (fp == nullptr)
    {
        fprintf(stderr, "Sorry, can not trunc %s, truncation of this filetype is not supported\n", fn);
        gmx_trr_close(in);
    }
    else
    {
        j     = 0;
        fpos  = gmx_fio_ftell(in);
        bStop = FALSE;
        while (!bStop && gmx_trr_read_frame_header(in, &sh, &bOK))
        {
            gmx_trr_read_frame_data(in, &sh, nullptr, nullptr, nullptr, nullptr);
            fpos = gmx_ftell(fp);
            t    = sh.t;
            if (t >= t0)
            {
                gmx_fseek(fp, fpos, SEEK_SET);
                bStop = TRUE;
            }
        }
        if (bStop)
        {
            fprintf(stderr, "Do you REALLY want to truncate this trajectory (%s) at:\n"
                    "frame %d, time %g, bytes %ld ??? (type YES if so)\n",
                    fn, j, t, (long int)fpos);
            if (1 != scanf("%s", yesno))
            {
                gmx_fatal(FARGS, "Error reading user input");
            }
            if (std::strcmp(yesno, "YES") == 0)
            {
                fprintf(stderr, "Once again, I'm gonna DO this...\n");
                gmx_trr_close(in);
                if (0 != gmx_truncate(fn, fpos))
                {
                    gmx_fatal(FARGS, "Error truncating file %s", fn);
                }
            }
            else
            {
                fprintf(stderr, "Ok, I'll forget about it\n");
            }
        }
        else
        {
            fprintf(stderr, "Already at end of file (t=%g)...\n", t);
            gmx_trr_close(in);
        }
    }
}

/*! \brief Read a full molecular topology if useful and available.
 *
 * If the input trajectory file is not in TNG format, and the output
 * file is in TNG format, then we want to try to read a full topology
 * (if available), so that we can write molecule information to the
 * output file. The full topology provides better molecule information
 * than is available from the normal t_topology data used by GROMACS
 * tools.
 *
 * Also, the t_topology is only read under (different) particular
 * conditions. If both apply, then a .tpr file might be read
 * twice. Trying to fix this redundancy while trjconv is still an
 * all-purpose tool does not seem worthwhile.
 *
 * Because of the way gmx_prepare_tng_writing is implemented, the case
 * where the input TNG file has no molecule information will never
 * lead to an output TNG file having molecule information. Since
 * molecule information will generally be present if the input TNG
 * file was written by a GROMACS tool, this seems like reasonable
 * behaviour. */
static gmx_mtop_t *read_mtop_for_tng(const char *tps_file,
                                     const char *input_file,
                                     const char *output_file)
{
    gmx_mtop_t *mtop = nullptr;

    if (fn2bTPX(tps_file) &&
        efTNG != fn2ftp(input_file) &&
        efTNG == fn2ftp(output_file))
    {
        int temp_natoms = -1;
        snew(mtop, 1);
        read_tpx(tps_file, nullptr, nullptr, &temp_natoms,
                 nullptr, nullptr, mtop);
    }

    return mtop;
}

int gmx_trjconv(int argc, char *argv[])
{
    const char *desc[] = {
        "[THISMODULE] can convert trajectory files between different formats, and supports",
        "some additonal options:",
	"",
        "* reduce the number of frames",
        "* change the timestamps of the frames ([TT]-t0[tt] and [TT]-timestep[tt])",
        "* cut the trajectory in small subtrajectories according",
        "  to information in an index file. This allows subsequent analysis of",
        "  the subtrajectories that could, for example, be the result of a",
        "  cluster analysis. Use option [TT]-sub[tt].",
        "  This assumes that the entries in the index file are frame numbers and",
        "  dumps each group in the index file to a separate trajectory file.",
        "* select frames within a certain range of a quantity given",
        "  in an [REF].xvg[ref] file.",
        "",
        "[gmx-trjcat] is better suited for concatenating multiple trajectory files.",
        "[PAR]",

        "The following formats are supported for input and output:",
        "[REF].xtc[ref], [REF].trr[ref], [REF].gro[ref], [TT].g96[tt]",
        "and [REF].pdb[ref].",
        "The file formats are detected from the file extension.",
        "The precision of [REF].xtc[ref] and [REF].gro[ref] output is taken from the",
        "input file for [REF].xtc[ref], [REF].gro[ref] and [REF].pdb[ref],",
        "and from the [TT]-ndec[tt] option for other input formats. The precision",
        "is always taken from [TT]-ndec[tt], when this option is set.",
        "All other formats have fixed precision. [REF].trr[ref]",
        "output can be single or double precision, depending on the precision",
        "of the [THISMODULE] binary.",
        "Note that velocities are only supported in",
        "[REF].trr[ref], [REF].gro[ref] and [TT].g96[tt] files.[PAR]",

        "Option [TT]-sep[tt] can be used to write every frame to a separate",
        "[TT].gro, .g96[tt] or [REF].pdb[ref] file. By default, all frames all written to one file.",
        "[REF].pdb[ref] files with all frames concatenated can be viewed with",
        "[TT]rasmol -nmrpdb[tt].[PAR]",

        "[BB]ALWAYS[bb] put the original trajectory on tape!",
        "We recommend to use the portable [REF].xtc[ref] format for your analysis",
        "to save disk space and to have portable files.[PAR]",

        "With [TT]-dt[tt], it is possible to reduce the number of ",
        "frames in the output. This option relies on the accuracy of the times",
        "in your input trajectory, so if these are inaccurate use the",
        "[TT]-timestep[tt] option to modify the time (this can be done",
        "simultaneously). For making smooth movies, the program [gmx-filter]",
        "can reduce the number of frames while using low-pass frequency",
        "filtering, this reduces aliasing of high frequency motions.[PAR]",

        "Using [TT]-trunc[tt] [THISMODULE] can truncate [REF].trr[ref] in place, i.e.",
        "without copying the file. This is useful when a run has crashed",
        "during disk I/O (i.e. full disk), or when two contiguous",
        "trajectories must be concatenated without having double frames.[PAR]",

        "Option [TT]-dump[tt] can be used to extract a frame at or near",
        "one specific time from your trajectory, but only works reliably",
        "if the time interval between frames is uniform.[PAR]",

        "Option [TT]-drop[tt] reads an [REF].xvg[ref] file with times and values.",
        "When options [TT]-dropunder[tt] and/or [TT]-dropover[tt] are set,",
        "frames with a value below and above the value of the respective options",
        "will not be written."
    };

    static gmx_bool  bSeparate     = FALSE, bVels = TRUE, bForce = FALSE, bCONECT = FALSE;
    static int       skip_nr       = 1, ndec = 3, nzero = 0;
    static real      tzero         = 0, delta_t = 0, timestep = 0, ttrunc = -1, tdump = -1, split_t = 0;
    static char     *exec_command  = nullptr;
    static real      dropunder     = 0, dropover = 0;
    static gmx_bool  bRound        = FALSE;

    t_pargs
        pa[] =
    {
        { "-skip", FALSE, etINT,
          { &skip_nr }, "Only write every nr-th frame" },
        { "-dt", FALSE, etTIME,
          { &delta_t },
          "Only write frame when t MOD dt = first time (%t)" },
        { "-round", FALSE, etBOOL,
          { &bRound }, "Round measurements to nearest picosecond"},
        { "-dump", FALSE, etTIME,
          { &tdump }, "Dump frame nearest specified time (%t)" },
        { "-t0", FALSE, etTIME,
          { &tzero },
          "Starting time (%t) (default: don't change)" },
        { "-timestep", FALSE, etTIME,
          { &timestep },
          "Change time step between input frames (%t)" },
        { "-ndec", FALSE, etINT,
          { &ndec },
          "Precision for .xtc and .gro writing in number of "
          "decimal places" },
        { "-vel", FALSE, etBOOL,
          { &bVels }, "Read and write velocities if possible" },
        { "-force", FALSE, etBOOL,
          { &bForce }, "Read and write forces if possible" },
        { "-trunc", FALSE, etTIME,
          { &ttrunc },
          "Truncate input trajectory file after this time (%t)" },
        { "-exec", FALSE, etSTR,
          { &exec_command },
          "Execute command for every output frame with the "
          "frame number as argument" },
        { "-split", FALSE, etTIME,
          { &split_t },
          "Start writing new file when t MOD split = first "
          "time (%t)" },
        { "-sep", FALSE, etBOOL,
          { &bSeparate },
          "Write each frame to a separate .gro, .g96 or .pdb "
          "file" },
        { "-nzero", FALSE, etINT,
          { &nzero },
          "If the -sep flag is set, use these many digits "
          "for the file numbers and prepend zeros as needed" },
        { "-dropunder", FALSE, etREAL,
          { &dropunder }, "Drop all frames below this value" },
        { "-dropover", FALSE, etREAL,
          { &dropover }, "Drop all frames above this value" },
        { "-conect", FALSE, etBOOL,
          { &bCONECT },
          "Add conect records when writing [REF].pdb[ref] files. Useful "
          "for visualization of non-standard molecules, e.g. "
          "coarse grained ones" }
    };
#define NPA asize(pa)

    FILE             *out    = nullptr;
    t_trxstatus      *trxout = nullptr;
    t_trxstatus      *trxin;
    int               file_nr;
    t_trxframe        fr, frout;
    int               flags;
    rvec             *xmem  = nullptr, *vmem = nullptr, *fmem = nullptr;
    rvec             *xp    = nullptr, x_shift;
    int               i, frame, outframe, natoms, nout, newstep = 0, model_nr;
#define SKIP 10
    t_topology        top;
    gmx_mtop_t       *mtop  = nullptr;
    gmx_conect        gc    = nullptr;
    int               ePBC  = -1;
    t_atoms          *atoms = nullptr, useatoms;
    matrix            top_box;
    int              *index = nullptr, *cindex = nullptr;
    char             *grpnm = nullptr;
    int              *frindex, nrfri;
    char             *frname;
    int               my_clust = -1;
    t_cluster_ndx    *clust           = nullptr;
    t_trxstatus     **clust_status    = nullptr;
    int              *clust_status_id = nullptr;
    int               ntrxopen        = 0;
    int              *nfwritten       = nullptr;
    int               ndrop           = 0, ncol, drop0 = 0, drop1 = 0, dropuse = 0;
    double          **dropval;
    real              tshift = 0, dt = -1, prec;
    gmx_bool          bCopy, bDoIt, bIndex, bTDump, bSetTime, bTPS = FALSE, bDTset = FALSE;
    gmx_bool          bExec, bTimeStep = FALSE, bDumpFrame = FALSE, bSetPrec, bNeedPrec;
    gmx_bool          bHaveFirstFrame, bHaveNextFrame, bSplit = FALSE;
    gmx_bool          bSubTraj = FALSE, bDropUnder = FALSE, bDropOver = FALSE;
    gmx_bool          bWriteFrame, bSplitHere;
    const char       *top_file, *in_file, *out_file = nullptr;
    char              out_file2[256], *charpt;
    char             *outf_base = nullptr;
    const char       *outf_ext  = nullptr;
    char              top_title[256], title[256], timestr[32], stepstr[32], filemode[5];
    gmx_output_env_t *oenv;

    t_filenm          fnm[] = {
        { efTRX, "-f",   nullptr,      ffREAD  },
        { efTRO, "-o",   nullptr,      ffWRITE },
        { efTPS, nullptr,   nullptr,      ffOPTRD },
        { efNDX, nullptr,   nullptr,      ffOPTRD },
        { efNDX, "-fr",  "frames",  ffOPTRD },
        { efNDX, "-sub", "cluster", ffOPTRD },
        { efXVG, "-drop", "drop",    ffOPTRD }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv,
                           PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_VIEW |
                           PCA_TIME_UNIT,
                           NFILE, fnm, NPA, pa, asize(desc), desc,
                           0, nullptr, &oenv))
    {
        return 0;
    }

    top_file = ftp2fn(efTPS, NFILE, fnm);

    /* Check command line */
    in_file = opt2fn("-f", NFILE, fnm);

    if (ttrunc != -1)
    {
        do_trunc(in_file, ttrunc);
    }
    else
    {
        /* mark active cmdline options */

        bSetTime   = opt2parg_bSet("-t0", NPA, pa);
        bSetPrec   = opt2parg_bSet("-ndec", NPA, pa);

        bExec      = opt2parg_bSet("-exec", NPA, pa);
        bTimeStep  = opt2parg_bSet("-timestep", NPA, pa);
        bTDump     = opt2parg_bSet("-dump", NPA, pa);
        bDropUnder = opt2parg_bSet("-dropunder", NPA, pa);
        bDropOver  = opt2parg_bSet("-dropover", NPA, pa);

        bSplit     = (split_t != 0);

        /* ndec is in nr of decimal places, prec is a multiplication factor: */
        prec = 1;
        for (i = 0; i < ndec; i++)
        {
            prec *= 10;
        }

        bIndex = ftp2bSet(efNDX, NFILE, fnm);


        /* Determine output type */
        out_file = opt2fn("-o", NFILE, fnm);
        int ftp  = fn2ftp(out_file);
        fprintf(stderr, "Will write %s: %s\n", ftp2ext(ftp), ftp2desc(ftp));
        bNeedPrec = (ftp == efXTC || ftp == efGRO);
        int ftpin = fn2ftp(in_file);
        if (bVels)
        {
            /* check if velocities are possible in input and output files */
            bVels = (ftp == efTRR || ftp == efGRO ||
                     ftp == efG96 || ftp == efTNG)
                && (ftpin == efTRR || ftpin == efGRO ||
                    ftpin == efG96 || ftpin == efTNG || ftpin == efCPT);
        }
        if (bSeparate || bSplit)
        {
            outf_ext = std::strrchr(out_file, '.');
            if (outf_ext == nullptr)
            {
                gmx_fatal(FARGS, "Output file name '%s' does not contain a '.'", out_file);
            }
            outf_base = gmx_strdup(out_file);
            outf_base[outf_ext - out_file] = '\0';
        }
        bSubTraj = opt2bSet("-sub", NFILE, fnm);
        if (bSubTraj)
        {
            if ((ftp != efXTC) && (ftp != efTRR))
            {
                /* It seems likely that other trajectory file types
                 * could work here. */
                gmx_fatal(FARGS, "Can only use the sub option with output file types "
                          "xtc and trr");
            }
            clust = cluster_index(nullptr, opt2fn("-sub", NFILE, fnm));

            /* Check for number of files disabled, as FOPEN_MAX is not the correct
             * number to check for. In my linux box it is only 16.
             */
            if (0 && (clust->clust->nr > FOPEN_MAX-4))
            {
                gmx_fatal(FARGS, "Can not open enough (%d) files to write all the"
                          " trajectories.\ntry splitting the index file in %d parts.\n"
                          "FOPEN_MAX = %d",
                          clust->clust->nr, 1+clust->clust->nr/FOPEN_MAX, FOPEN_MAX);
            }
            gmx_warning("The -sub option could require as many open output files as there are\n"
                        "index groups in the file (%d). If you get I/O errors opening new files,\n"
                        "try reducing the number of index groups in the file, and perhaps\n"
                        "using trjconv -sub several times on different chunks of your index file.\n",
                        clust->clust->nr);

            snew(clust_status, clust->clust->nr);
            snew(clust_status_id, clust->clust->nr);
            snew(nfwritten, clust->clust->nr);
            for (i = 0; (i < clust->clust->nr); i++)
            {
                clust_status[i]    = nullptr;
                clust_status_id[i] = -1;
            }
            bSeparate = bSplit = FALSE;
        }



        /* skipping */
// what is this used for????
        if (skip_nr <= 0)
        {
        }

        mtop = read_mtop_for_tng(top_file, in_file, out_file);

        /* Determine whether to read a topology */
        bTPS = (ftp2bSet(efTPS, NFILE, fnm) ||
                (ftp == efGRO) || (ftp == efPDB) || bCONECT);

        /* Determine if when can read index groups */
        bIndex = (bIndex || bTPS);

        if (bTPS)
        {
            read_tps_conf(top_file, &top, &ePBC, &xp, nullptr, top_box, false);
            std::strncpy(top_title, *top.name, 255);
            top_title[255] = '\0';
            atoms          = &top.atoms;

            /* top_title is only used for gro and pdb,
             * the header in such a file is top_title, followed by
             * t= ... and/or step= ...
             * to prevent double t= or step=, remove it from top_title.
             * From GROMACS-2018 we only write t/step when the frame actually
             * has a valid time/step, so we need to check for both separately.
             */
            if ((charpt = std::strstr(top_title, " t= ")))
            {
                charpt[0] = '\0';
            }
            if ((charpt = std::strstr(top_title, " step= ")))
            {
                charpt[0] = '\0';
            }

            if (bCONECT)
            {
                gc = gmx_conect_generate(&top);
            }
        }

        /* get frame number index */
        frindex = nullptr;
        if (opt2bSet("-fr", NFILE, fnm))
        {
            printf("Select groups of frame number indices:\n");
            rd_index(opt2fn("-fr", NFILE, fnm), 1, &nrfri, (int **)&frindex, &frname);
            if (debug)
            {
                for (i = 0; i < nrfri; i++)
                {
                    fprintf(debug, "frindex[%4d]=%4d\n", i, frindex[i]);
                }
            }
        }
        if (bIndex)
        {
            printf("Select group for output\n");
            get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm),
                      1, &nout, &index, &grpnm);
        }
        else
        {
            /* no index file, so read natoms from TRX */
            if (!read_first_frame(oenv, &trxin, in_file, &fr, TRX_DONT_SKIP))
            {
                gmx_fatal(FARGS, "Could not read a frame from %s", in_file);
            }
            natoms = fr.natoms;
            close_trx(trxin);
            sfree(fr.x);
            snew(index, natoms);
            for (i = 0; i < natoms; i++)
            {
                index[i] = i;
            }
            nout = natoms;
       }

        {
            clear_rvec(x_shift);
        }

        if (bDropUnder || bDropOver)
        {
            /* Read the .xvg file with the drop values */
            fprintf(stderr, "\nReading drop file ...");
            ndrop = read_xvg(opt2fn("-drop", NFILE, fnm), &dropval, &ncol);
            fprintf(stderr, " %d time points\n", ndrop);
            if (ndrop == 0 || ncol < 2)
            {
                gmx_fatal(FARGS, "Found no data points in %s",
                          opt2fn("-drop", NFILE, fnm));
            }
            drop0 = 0;
            drop1 = 0;
        }

        /* Make atoms struct for output in GRO or PDB files */
        if ((ftp == efGRO) || ((ftp == efG96) && bTPS) || (ftp == efPDB))
        {
            /* get memory for stuff to go in .pdb file, and initialize
             * the pdbinfo structure part if the input has it.
             */
            init_t_atoms(&useatoms, atoms->nr, atoms->havePdbInfo);
            sfree(useatoms.resinfo);
            useatoms.resinfo = atoms->resinfo;
            for (i = 0; (i < nout); i++)
            {
                useatoms.atomname[i] = atoms->atomname[index[i]];
                useatoms.atom[i]     = atoms->atom[index[i]];
                if (atoms->havePdbInfo)
                {
                    useatoms.pdbinfo[i]  = atoms->pdbinfo[index[i]];
                }
                useatoms.nres        = std::max(useatoms.nres, useatoms.atom[i].resind+1);
            }
            useatoms.nr = nout;
        }
        /* select what to read */
        if (ftp == efTRR)
        {
            flags = TRX_READ_X;
        }
        else
        {
            flags = TRX_NEED_X;
        }
        if (bVels)
        {
            flags = flags | TRX_READ_V;
        }
        if (bForce)
        {
            flags = flags | TRX_READ_F;
        }

        /* open trx file for reading */
        bHaveFirstFrame = read_first_frame(oenv, &trxin, in_file, &fr, flags);
        if (fr.bPrec)
        {
            fprintf(stderr, "\nPrecision of %s is %g (nm)\n", in_file, 1/fr.prec);
        }
        if (bNeedPrec)
        {
            if (bSetPrec || !fr.bPrec)
            {
                fprintf(stderr, "\nSetting output precision to %g (nm)\n", 1/prec);
            }
            else
            {
                fprintf(stderr, "Using output precision of %g (nm)\n", 1/prec);
            }
        }

        if (bHaveFirstFrame)
        {
            if (bTDump)
            {
                // Determine timestep (assuming constant spacing for now) if we
                // need to dump frames based on time. This is required so we do not
                // skip the first frame if that was the one that should have been dumped
                double firstFrameTime = fr.time;
                if (read_next_frame(oenv, trxin, &fr))
                {
                    dt     = fr.time - firstFrameTime;
                    bDTset = TRUE;
                    if (dt <= 0)
                    {
                        fprintf(stderr, "Warning: Frame times are not incrementing - will dump first frame.\n");
                    }
                }
                // Now close and reopen so we are at first frame again
                close_trx(trxin);
                done_frame(&fr);
                // Reopen at first frame (We already know it exists if we got here)
                read_first_frame(oenv, &trxin, in_file, &fr, flags);
            }

            set_trxframe_ePBC(&fr, ePBC);
            natoms = fr.natoms;

            if (bSetTime)
            {
                tshift = tzero-fr.time;
            }
            else
            {
                tzero = fr.time;
            }

            bCopy = FALSE;
            if (bIndex)
            {
                /* check if index is meaningful */
                for (i = 0; i < nout; i++)
                {
                    if (index[i] >= natoms)
                    {
                        gmx_fatal(FARGS,
                                  "Index[%d] %d is larger than the number of atoms in the\n"
                                  "trajectory file (%d). There is a mismatch in the contents\n"
                                  "of your -f, -s and/or -n files.", i, index[i]+1, natoms);
                    }
                    bCopy = bCopy || (i != index[i]);
                }
            }

            /* open output for writing */
            std::strcpy(filemode, "w");
            switch (ftp)
            {
                case efTNG:
                    trjtools_gmx_prepare_tng_writing(out_file,
                                                     filemode[0],
                                                     trxin,
                                                     &trxout,
                                                     nullptr,
                                                     nout,
                                                     mtop,
                                                     index,
                                                     grpnm);
                    break;
                case efXTC:
                case efTRR:
                    out = nullptr;
                    if (!bSplit && !bSubTraj)
                    {
                        trxout = open_trx(out_file, filemode);
                    }
                    break;
                case efGRO:
                case efG96:
                case efPDB:
                    if (( !bSeparate && !bSplit ) && !bSubTraj)
                    {
                        out = gmx_ffopen(out_file, filemode);
                    }
                    break;
                default:
                    gmx_incons("Illegal output file format");
            }

            if (bCopy)
            {
                snew(xmem, nout);
                if (bVels)
                {
                    snew(vmem, nout);
                }
                if (bForce)
                {
                    snew(fmem, nout);
                }
            }

            /* Start the big loop over frames */
            file_nr  =  0;
            frame    =  0;
            outframe =  0;
            model_nr =  0;

            /* Main loop over frames */
            do
            {
                if (!fr.bStep)
                {
                    /* set the step */
                    fr.step = newstep;
                    newstep++;
                }
                if (bSubTraj)
                {
                    /*if (frame >= clust->clust->nra)
                       gmx_fatal(FARGS,"There are more frames in the trajectory than in the cluster index file\n");*/
                    if (frame > clust->maxframe)
                    {
                        my_clust = -1;
                    }
                    else
                    {
                        my_clust = clust->inv_clust[frame];
                    }
                    if ((my_clust < 0) || (my_clust >= clust->clust->nr) ||
                        (my_clust == -1))
                    {
                        my_clust = -1;
                    }
                }

                 if (bTDump)
                {
                    // If we could not read two frames or times are not incrementing
                    // we have almost no idea what to do,
                    // but dump the first frame so output is not broken.
                    if (dt <= 0 || !bDTset)
                    {
                        bDumpFrame = true;
                    }
                    else
                    {
                        // Dump the frame if we are less than half a frame time
                        // below it. This will also ensure we at least dump a
                        // somewhat reasonable frame if the spacing is unequal
                        // and we have overrun the frame time. Once we dump one
                        // frame based on time we quit, so it does not matter
                        // that this might be true for all subsequent frames too.
                        bDumpFrame = (fr.time > tdump-0.5*dt);
                    }
                }
                else
                {
                    bDumpFrame = FALSE;
                }


                if (frindex)
                {
                    /* see if we have a frame from the frame index group */
                    for (i = 0; i < nrfri && !bDumpFrame; i++)
                    {
                        bDumpFrame = frame == frindex[i];
                    }
                }
                if (debug && bDumpFrame)
                {
                    fprintf(debug, "dumping %d\n", frame);
                }

                bWriteFrame =
                    ( ( !bTDump && !frindex && frame % skip_nr == 0 ) || bDumpFrame );

                if (bWriteFrame && (bDropUnder || bDropOver))
                {
                    while (dropval[0][drop1] < fr.time && drop1+1 < ndrop)
                    {
                        drop0 = drop1;
                        drop1++;
                    }
                    if (std::abs(dropval[0][drop0] - fr.time)
                        < std::abs(dropval[0][drop1] - fr.time))
                    {
                        dropuse = drop0;
                    }
                    else
                    {
                        dropuse = drop1;
                    }
                    if ((bDropUnder && dropval[1][dropuse] < dropunder) ||
                        (bDropOver && dropval[1][dropuse] > dropover))
                    {
                        bWriteFrame = FALSE;
                    }
                }

                if (bWriteFrame)
                {
                    /* We should avoid modifying the input frame,
                     * but since here we don't have the output frame yet,
                     * we introduce a temporary output frame time variable.
                     */
                    real frout_time;

                    frout_time = fr.time;

                    /* calc new time */
                    if (bTimeStep)
                    {
                        frout_time = tzero + frame*timestep;
                    }
                    else
                    if (bSetTime)
                    {
                        frout_time += tshift;
                    }

                    if (bTDump)
                    {
                        fprintf(stderr, "\nDumping frame at t= %g %s\n",
                                output_env_conv_time(oenv, frout_time), output_env_get_time_unit(oenv).c_str());
                    }

                    /* check for writing at each delta_t */
                    bDoIt = (delta_t == 0);
                    if (!bDoIt)
                    {
                        if (!bRound)
                        {
                            bDoIt = bRmod(frout_time, tzero, delta_t);
                        }
                        else
                        {
                            /* round() is not C89 compatible, so we do this:  */
                            bDoIt = bRmod(std::floor(frout_time+0.5), std::floor(tzero+0.5),
                                          std::floor(delta_t+0.5));
                        }
                    }

                    if (bDoIt || bTDump)
                    {
                        /* print sometimes */
                        if ( ((outframe % SKIP) == 0) || (outframe < SKIP) )
                        {
                            fprintf(stderr, " ->  frame %6d time %8.3f      \r",
                                    outframe, output_env_conv_time(oenv, frout_time));
                            fflush(stderr);
                        }

                        /* Copy the input trxframe struct to the output trxframe struct */
                        frout        = fr;
                        frout.time   = frout_time;
                        frout.bV     = (frout.bV && bVels);
                        frout.bF     = (frout.bF && bForce);
                        frout.natoms = nout;
                        if (bNeedPrec && (bSetPrec || !fr.bPrec))
                        {
                            frout.bPrec = TRUE;
                            frout.prec  = prec;
                        }
                        if (bCopy)
                        {
                            frout.x = xmem;
                            if (frout.bV)
                            {
                                frout.v = vmem;
                            }
                            if (frout.bF)
                            {
                                frout.f = fmem;
                            }
                            for (i = 0; i < nout; i++)
                            {
                                copy_rvec(fr.x[index[i]], frout.x[i]);
                                if (frout.bV)
                                {
                                    copy_rvec(fr.v[index[i]], frout.v[i]);
                                }
                                if (frout.bF)
                                {
                                    copy_rvec(fr.f[index[i]], frout.f[i]);
                                }
                            }
                        }


                        if (!bRound)
                        {
                            bSplitHere = bSplit && bRmod(frout.time, tzero, split_t);
                        }
                        else
                        {
                            /* round() is not C89 compatible, so we do this: */
                            bSplitHere = bSplit && bRmod(std::floor(frout.time+0.5),
                                                         std::floor(tzero+0.5),
                                                         std::floor(split_t+0.5));
                        }
                        if (bSeparate || bSplitHere)
                        {
                            mk_filenm(outf_base, ftp2ext(ftp), nzero, file_nr, out_file2);
                        }

                        switch (ftp)
                        {
                            case efTNG:
                                write_tng_frame(trxout, &frout);
                                // TODO when trjconv behaves better: work how to read and write lambda
                                break;
                            case efTRR:
                            case efXTC:
                                if (bSplitHere)
                                {
                                    if (trxout)
                                    {
                                        close_trx(trxout);
                                    }
                                    trxout = open_trx(out_file2, filemode);
                                }
                                if (bSubTraj)
                                {
                                    if (my_clust != -1)
                                    {
                                        char buf[STRLEN];
                                        if (clust_status_id[my_clust] == -1)
                                        {
                                            sprintf(buf, "%s.%s", clust->grpname[my_clust], ftp2ext(ftp));
                                            clust_status[my_clust]    = open_trx(buf, "w");
                                            clust_status_id[my_clust] = 1;
                                            ntrxopen++;
                                        }
                                        else if (clust_status_id[my_clust] == -2)
                                        {
                                            gmx_fatal(FARGS, "File %s.xtc should still be open (%d open .xtc files)\n" "in order to write frame %d. my_clust = %d",
                                                      clust->grpname[my_clust], ntrxopen, frame,
                                                      my_clust);
                                        }
                                        write_trxframe(clust_status[my_clust], &frout, gc);
                                        nfwritten[my_clust]++;
                                        if (nfwritten[my_clust] ==
                                            (clust->clust->index[my_clust+1]-
                                             clust->clust->index[my_clust]))
                                        {
                                            close_trx(clust_status[my_clust]);
                                            clust_status[my_clust]    = nullptr;
                                            clust_status_id[my_clust] = -2;
                                            ntrxopen--;
                                            if (ntrxopen < 0)
                                            {
                                                gmx_fatal(FARGS, "Less than zero open .xtc files!");
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    write_trxframe(trxout, &frout, gc);
                                }
                                break;
                            case efGRO:
                            case efG96:
                            case efPDB:
                                // Only add a generator statement if title is empty,
                                // to avoid multiple generated-by statements from various programs
                                if (std::strlen(top_title) == 0)
                                {
                                    sprintf(top_title, "Generated by trjconv");
                                }
                                if (frout.bTime)
                                {
                                    sprintf(timestr, " t= %9.5f", frout.time);
                                }
                                else
                                {
                                    std::strcpy(timestr, "");
                                }
                                if (frout.bStep)
                                {
                                    sprintf(stepstr, " step= %" GMX_PRId64, frout.step);
                                }
                                else
                                {
                                    std::strcpy(stepstr, "");
                                }
                                snprintf(title, 256, "%s%s%s", top_title, timestr, stepstr);
                                if (bSeparate || bSplitHere)
                                {
                                    out = gmx_ffopen(out_file2, "w");
                                }
                                // switch is own the same file type as the first output file type
                                switch (ftp)
                                {
                                    case efGRO:
                                        write_hconf_p(out, title, &useatoms,
                                                      frout.x, frout.bV ? frout.v : nullptr, frout.box);
                                        break;
                                    case efPDB:
                                        fprintf(out, "REMARK    GENERATED BY TRJCONV\n");
                                        /* if reading from pdb, we want to keep the original
                                           model numbering else we write the output frame
                                           number plus one, because model 0 is not allowed in pdb */
                                        if (ftpin == efPDB && fr.bStep && fr.step > model_nr)
                                        {
                                            model_nr = fr.step;
                                        }
                                        else
                                        {
                                            model_nr++;
                                        }
                                        write_pdbfile(out, title, &useatoms, frout.x,
                                                      frout.ePBC, frout.box, ' ', model_nr, gc, TRUE);
                                        break;
                                    case efG96:
                                        const char *outputTitle = "";
                                        if (bSeparate || bTDump)
                                        {
                                            outputTitle = title;
                                            if (bTPS)
                                            {
                                                frout.bAtoms = TRUE;
                                            }
                                            frout.atoms  = &useatoms;
                                            frout.bStep  = FALSE;
                                            frout.bTime  = FALSE;
                                        }
                                        else
                                        {
                                            if (outframe == 0)
                                            {
                                                outputTitle = title;
                                            }
                                            frout.bAtoms = FALSE;
                                            frout.bStep  = TRUE;
                                            frout.bTime  = TRUE;
                                        }
                                        write_g96_conf(out, outputTitle, &frout, -1, nullptr);
                                }
                                if (bSeparate || bSplitHere)
                                {
                                    gmx_ffclose(out);
                                    out = nullptr;
                                }
                                break;
                            default:
                                gmx_fatal(FARGS, "DHE, ftp=%d\n", ftp);
                        }
                        if (bSeparate || bSplitHere)
                        {
                            file_nr++;
                        }

                        /* execute command */
                        if (bExec)
                        {
                            char c[255];
                            sprintf(c, "%s  %d", exec_command, file_nr-1);
                            /*fprintf(stderr,"Executing '%s'\n",c);*/
                            if (0 != system(c))
                            {
                                gmx_fatal(FARGS, "Error executing command: %s", c);
                            }
                        }
                        outframe++;
                    }
                }
                frame++;
                bHaveNextFrame = read_next_frame(oenv, trxin, &fr);
            }
            while (!(bTDump && bDumpFrame) && bHaveNextFrame);
        }

        if (!bHaveFirstFrame || (bTDump && !bDumpFrame))
        {
            fprintf(stderr, "\nWARNING no output, "
                    "last frame read at t=%g\n", fr.time);
        }
        fprintf(stderr, "\n");

        close_trx(trxin);
        sfree(outf_base);

        if (trxout)
        {
            close_trx(trxout);
        }
        else if (out != nullptr)
        {
            gmx_ffclose(out);
        }
        if (bSubTraj)
        {
            for (i = 0; (i < clust->clust->nr); i++)
            {
                if (clust_status_id[i] >= 0)
                {
                    close_trx(clust_status[i]);
                }
            }
        }
    }

    sfree(mtop);
    done_top(&top);
    sfree(xp);
    sfree(xmem);
    sfree(vmem);
    sfree(fmem);
    sfree(grpnm);
    sfree(index);
    sfree(cindex);
    done_filenms(NFILE, fnm);
    done_frame(&fr);

    do_view(oenv, out_file, nullptr);

    output_env_done(oenv);
    return 0;
}

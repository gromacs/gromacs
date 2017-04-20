/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2013,2014,2015,2017, by the GROMACS development team, led by
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

#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/dir_separator.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

static int strip_dssp(char *dsspfile, int nres,
                      gmx_bool bPhobres[], real t,
                      real *acc, FILE *fTArea,
                      t_matrix *mat, int average_area[],
                      const gmx_output_env_t *oenv)
{
    static gmx_bool bFirst = TRUE;
    static char    *ssbuf;
    FILE           *tapeout;
    static int      xsize, frame;
    char            buf[STRLEN+1];
    char            SSTP;
    int             i, nr, iacc, nresidues;
    int             naccf, naccb; /* Count hydrophobic and hydrophilic residues */
    real            iaccf, iaccb;
    t_xpmelmt       c;

    tapeout = gmx_ffopen(dsspfile, "r");

    /* Skip header */
    do
    {
        fgets2(buf, STRLEN, tapeout);
    }
    while (std::strstr(buf, "KAPPA") == nullptr);
    if (bFirst)
    {
        /* Since we also have empty lines in the dssp output (temp) file,
         * and some line content is saved to the ssbuf variable,
         * we need more memory than just nres elements. To be shure,
         * we allocate 2*nres-1, since for each chain there is a
         * separating line in the temp file. (At most each residue
         * could have been defined as a separate chain.) */
        snew(ssbuf, 2*nres-1);
    }

    iaccb     = iaccf = 0;
    nresidues = 0;
    naccf     = 0;
    naccb     = 0;
    for (nr = 0; (fgets2(buf, STRLEN, tapeout) != nullptr); nr++)
    {
        if (buf[13] == '!') /* Chain separator line has '!' at pos. 13 */
        {
            SSTP = '=';     /* Chain separator sign '=' */
        }
        else
        {
            SSTP = buf[16] == ' ' ? '~' : buf[16];
        }
        ssbuf[nr] = SSTP;

        buf[39] = '\0';

        /* Only calculate solvent accessible area if needed */
        if ((nullptr != acc) && (buf[13] != '!'))
        {
            sscanf(&(buf[34]), "%d", &iacc);
            acc[nr] = iacc;
            /* average_area and bPhobres are counted from 0...nres-1 */
            average_area[nresidues] += iacc;
            if (bPhobres[nresidues])
            {
                naccb++;
                iaccb += iacc;
            }
            else
            {
                naccf++;
                iaccf += iacc;
            }
            /* Keep track of the residue number (do not count chain separator lines '!') */
            nresidues++;
        }
    }
    ssbuf[nr] = '\0';

    if (bFirst)
    {
        if (nullptr != acc)
        {
            fprintf(stderr, "%d residues were classified as hydrophobic and %d as hydrophilic.\n", naccb, naccf);
        }

        sprintf(mat->title, "Secondary structure");
        mat->legend[0] = 0;
        sprintf(mat->label_x, "%s", output_env_get_time_label(oenv).c_str());
        sprintf(mat->label_y, "Residue");
        mat->bDiscrete = TRUE;
        mat->ny        = nr;
        snew(mat->axis_y, nr);
        for (i = 0; i < nr; i++)
        {
            mat->axis_y[i] = i+1;
        }
        mat->axis_x = nullptr;
        mat->matrix = nullptr;
        xsize       = 0;
        frame       = 0;
        bFirst      = FALSE;
    }
    if (frame >= xsize)
    {
        xsize += 10;
        srenew(mat->axis_x, xsize);
        srenew(mat->matrix, xsize);
    }
    mat->axis_x[frame] = t;
    snew(mat->matrix[frame], nr);
    c.c2 = 0;
    for (i = 0; i < nr; i++)
    {
        c.c1                  = ssbuf[i];
        mat->matrix[frame][i] = std::max(static_cast<t_matelmt>(0), searchcmap(mat->nmap, mat->map, c));
    }
    frame++;
    mat->nx = frame;

    if (fTArea)
    {
        fprintf(fTArea, "%10g  %10g  %10g\n", t, 0.01*iaccb, 0.01*iaccf);
    }
    gmx_ffclose(tapeout);

    /* Return the number of lines found in the dssp file (i.e. number
     * of redidues plus chain separator lines).
     * This is the number of y elements needed for the area xpm file */
    return nr;
}

static gmx_bool *bPhobics(t_atoms *atoms)
{
    int         i, nb;
    char      **cb;
    gmx_bool   *bb;


    nb = get_lines("phbres.dat", &cb);
    snew(bb, atoms->nres);

    for (i = 0; (i < atoms->nres); i++)
    {
        if (-1 != search_str(nb, cb, *atoms->resinfo[i].name) )
        {
            bb[i] = TRUE;
        }
    }
    return bb;
}

static void check_oo(t_atoms *atoms)
{
    char *OOO;

    int   i;

    OOO = gmx_strdup("O");

    for (i = 0; (i < atoms->nr); i++)
    {
        if (std::strcmp(*(atoms->atomname[i]), "OXT") == 0)
        {
            *atoms->atomname[i] = OOO;
        }
        else if (std::strcmp(*(atoms->atomname[i]), "O1") == 0)
        {
            *atoms->atomname[i] = OOO;
        }
        else if (std::strcmp(*(atoms->atomname[i]), "OC1") == 0)
        {
            *atoms->atomname[i] = OOO;
        }
    }
}

static void norm_acc(t_atoms *atoms, int nres,
                     real av_area[], real norm_av_area[])
{
    int     i, n, n_surf;

    char    surffn[] = "surface.dat";
    char  **surf_res, **surf_lines;
    double *surf;

    n_surf = get_lines(surffn, &surf_lines);
    snew(surf, n_surf);
    snew(surf_res, n_surf);
    for (i = 0; (i < n_surf); i++)
    {
        snew(surf_res[i], 5);
        sscanf(surf_lines[i], "%s %lf", surf_res[i], &surf[i]);
    }

    for (i = 0; (i < nres); i++)
    {
        n = search_str(n_surf, surf_res, *atoms->resinfo[i].name);
        if (n != -1)
        {
            norm_av_area[i] = av_area[i] / surf[n];
        }
        else
        {
            fprintf(stderr, "Residue %s not found in surface database (%s)\n",
                    *atoms->resinfo[i].name, surffn);
        }
    }
}

void prune_ss_legend(t_matrix *mat)
{
    gmx_bool  *present;
    int       *newnum;
    int        i, r, f, newnmap;
    t_mapping *newmap;

    snew(present, mat->nmap);
    snew(newnum, mat->nmap);

    for (f = 0; f < mat->nx; f++)
    {
        for (r = 0; r < mat->ny; r++)
        {
            present[mat->matrix[f][r]] = TRUE;
        }
    }

    newnmap = 0;
    newmap  = nullptr;
    for (i = 0; i < mat->nmap; i++)
    {
        newnum[i] = -1;
        if (present[i])
        {
            newnum[i] = newnmap;
            newnmap++;
            srenew(newmap, newnmap);
            newmap[newnmap-1] = mat->map[i];
        }
    }
    if (newnmap != mat->nmap)
    {
        mat->nmap = newnmap;
        mat->map  = newmap;
        for (f = 0; f < mat->nx; f++)
        {
            for (r = 0; r < mat->ny; r++)
            {
                mat->matrix[f][r] = newnum[mat->matrix[f][r]];
            }
        }
    }
}

void write_sas_mat(const char *fn, real **accr, int nframe, int nres, t_matrix *mat)
{
    real  lo, hi;
    int   i, j, nlev;
    t_rgb rlo = {1, 1, 1}, rhi = {0, 0, 0};
    FILE *fp;

    if (fn)
    {
        hi = lo = accr[0][0];
        for (i = 0; i < nframe; i++)
        {
            for (j = 0; j < nres; j++)
            {
                lo = std::min(lo, accr[i][j]);
                hi = std::max(hi, accr[i][j]);
            }
        }
        fp   = gmx_ffopen(fn, "w");
        nlev = static_cast<int>(hi-lo+1);
        write_xpm(fp, 0, "Solvent Accessible Surface", "Surface (A^2)",
                  "Time", "Residue Index", nframe, nres,
                  mat->axis_x, mat->axis_y, accr, lo, hi, rlo, rhi, &nlev);
        gmx_ffclose(fp);
    }
}

void analyse_ss(const char *outfile, t_matrix *mat, const char *ss_string,
                const gmx_output_env_t *oenv)
{
    FILE        *fp;
    t_mapping   *map;
    int          f, r, *count, *total, ss_count, total_count;
    size_t       s;
    const char** leg;

    map = mat->map;
    snew(count, mat->nmap);
    snew(total, mat->nmap);
    snew(leg, mat->nmap+1);
    leg[0] = "Structure";
    for (s = 0; s < (size_t)mat->nmap; s++)
    {
        leg[s+1] = gmx_strdup(map[s].desc);
    }

    fp = xvgropen(outfile, "Secondary Structure",
                  output_env_get_xvgr_tlabel(oenv), "Number of Residues", oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(fp, "@ subtitle \"Structure = ");
    }
    for (s = 0; s < std::strlen(ss_string); s++)
    {
        if (s > 0)
        {
            fprintf(fp, " + ");
        }
        for (f = 0; f < mat->nmap; f++)
        {
            if (ss_string[s] == map[f].code.c1)
            {
                fprintf(fp, "%s", map[f].desc);
            }
        }
    }
    fprintf(fp, "\"\n");
    xvgr_legend(fp, mat->nmap+1, leg, oenv);

    total_count = 0;
    for (s = 0; s < (size_t)mat->nmap; s++)
    {
        total[s] = 0;
    }
    for (f = 0; f < mat->nx; f++)
    {
        ss_count = 0;
        for (s = 0; s < static_cast<size_t>(mat->nmap); s++)
        {
            count[s] = 0;
        }
        for (r = 0; r < mat->ny; r++)
        {
            count[mat->matrix[f][r]]++;
            total[mat->matrix[f][r]]++;
        }
        for (s = 0; s < static_cast<size_t>(mat->nmap); s++)
        {
            if (std::strchr(ss_string, map[s].code.c1))
            {
                ss_count    += count[s];
                total_count += count[s];
            }
        }
        fprintf(fp, "%8g %5d", mat->axis_x[f], ss_count);
        for (s = 0; s < static_cast<size_t>(mat->nmap); s++)
        {
            fprintf(fp, " %5d", count[s]);
        }
        fprintf(fp, "\n");
    }
    /* now print column totals */
    fprintf(fp, "%-8s %5d", "# Totals", total_count);
    for (s = 0; s < static_cast<size_t>(mat->nmap); s++)
    {
        fprintf(fp, " %5d", total[s]);
    }
    fprintf(fp, "\n");

    /* now print probabilities */
    fprintf(fp, "%-8s %5.2f", "# SS pr.", total_count / static_cast<real>(mat->nx * mat->ny));
    for (s = 0; s < static_cast<size_t>(mat->nmap); s++)
    {
        fprintf(fp, " %5.2f", total[s] / static_cast<real>(mat->nx * mat->ny));
    }
    fprintf(fp, "\n");

    xvgrclose(fp);
    sfree(leg);
    sfree(count);
}

int gmx_do_dssp(int argc, char *argv[])
{
    const char        *desc[] = {
        "[THISMODULE] ",
        "reads a trajectory file and computes the secondary structure for",
        "each time frame ",
        "calling the dssp program. If you do not have the dssp program,",
        "get it from http://swift.cmbi.ru.nl/gv/dssp. [THISMODULE] assumes ",
        "that the dssp executable is located in ",
        "[TT]/usr/local/bin/dssp[tt]. If this is not the case, then you should",
        "set an environment variable [TT]DSSP[tt] pointing to the dssp",
        "executable, e.g.: [PAR]",
        "[TT]setenv DSSP /opt/dssp/bin/dssp[tt][PAR]",
        "Since version 2.0.0, dssp is invoked with a syntax that differs",
        "from earlier versions. If you have an older version of dssp,",
        "use the [TT]-ver[tt] option to direct do_dssp to use the older syntax.",
        "By default, do_dssp uses the syntax introduced with version 2.0.0.",
        "Even newer versions (which at the time of writing are not yet released)",
        "are assumed to have the same syntax as 2.0.0.[PAR]",
        "The structure assignment for each residue and time is written to an",
        "[REF].xpm[ref] matrix file. This file can be visualized with for instance",
        "[TT]xv[tt] and can be converted to postscript with [TT]xpm2ps[tt].",
        "Individual chains are separated by light grey lines in the [REF].xpm[ref] and",
        "postscript files.",
        "The number of residues with each secondary structure type and the",
        "total secondary structure ([TT]-sss[tt]) count as a function of",
        "time are also written to file ([TT]-sc[tt]).[PAR]",
        "Solvent accessible surface (SAS) per residue can be calculated, both in",
        "absolute values (A^2) and in fractions of the maximal accessible",
        "surface of a residue. The maximal accessible surface is defined as",
        "the accessible surface of a residue in a chain of glycines.",
        "[BB]Note[bb] that the program [gmx-sas] can also compute SAS",
        "and that is more efficient.[PAR]",
        "Finally, this program can dump the secondary structure in a special file",
        "[TT]ssdump.dat[tt] for usage in the program [gmx-chi]. Together",
        "these two programs can be used to analyze dihedral properties as a",
        "function of secondary structure type."
    };
    static gmx_bool    bVerbose;
    static const char *ss_string   = "HEBT";
    static int         dsspVersion = 2;
    t_pargs            pa[]        = {
        { "-v",  FALSE, etBOOL, {&bVerbose},
          "HIDDENGenerate miles of useless information" },
        { "-sss", FALSE, etSTR, {&ss_string},
          "Secondary structures for structure count"},
        { "-ver", FALSE, etINT, {&dsspVersion},
          "DSSP major version. Syntax changed with version 2"}
    };

    t_trxstatus       *status;
    FILE              *tapein;
    FILE              *ss, *acc, *fTArea, *tmpf;
    const char        *fnSCount, *fnArea, *fnTArea, *fnAArea;
    const char        *leg[] = { "Phobic", "Phylic" };
    t_topology         top;
    int                ePBC;
    t_atoms           *atoms;
    t_matrix           mat;
    int                nres, nr0, naccr, nres_plus_separators;
    gmx_bool          *bPhbres, bDoAccSurf;
    real               t;
    int                i, j, natoms, nframe = 0;
    matrix             box = {{0}};
    int                gnx;
    char              *grpnm, *ss_str;
    int               *index;
    rvec              *xp, *x;
    int               *average_area;
    real             **accr, *accr_ptr = nullptr, *av_area, *norm_av_area;
    char               pdbfile[32], tmpfile[32];
    char               dssp[256];
    const char        *dptr;
    gmx_output_env_t  *oenv;
    gmx_rmpbc_t        gpbc = nullptr;

    t_filenm           fnm[] = {
        { efTRX, "-f",   nullptr,      ffREAD },
        { efTPS, nullptr,   nullptr,      ffREAD },
        { efNDX, nullptr,   nullptr,      ffOPTRD },
        { efDAT, "-ssdump", "ssdump", ffOPTWR },
        { efMAP, "-map", "ss",      ffLIBRD },
        { efXPM, "-o",   "ss",      ffWRITE },
        { efXVG, "-sc",  "scount",  ffWRITE },
        { efXPM, "-a",   "area",    ffOPTWR },
        { efXVG, "-ta",  "totarea", ffOPTWR },
        { efXVG, "-aa",  "averarea", ffOPTWR }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv,
                           PCA_CAN_TIME | PCA_CAN_VIEW | PCA_TIME_UNIT,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }
    fnSCount   = opt2fn("-sc", NFILE, fnm);
    fnArea     = opt2fn_null("-a", NFILE, fnm);
    fnTArea    = opt2fn_null("-ta", NFILE, fnm);
    fnAArea    = opt2fn_null("-aa", NFILE, fnm);
    bDoAccSurf = (fnArea || fnTArea || fnAArea);

    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &ePBC, &xp, nullptr, box, FALSE);
    atoms = &(top.atoms);
    check_oo(atoms);
    bPhbres = bPhobics(atoms);

    get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &gnx, &index, &grpnm);
    nres = 0;
    nr0  = -1;
    for (i = 0; (i < gnx); i++)
    {
        if (atoms->atom[index[i]].resind != nr0)
        {
            nr0 = atoms->atom[index[i]].resind;
            nres++;
        }
    }
    fprintf(stderr, "There are %d residues in your selected group\n", nres);

    std::strcpy(pdbfile, "ddXXXXXX");
    gmx_tmpnam(pdbfile);
    if ((tmpf = fopen(pdbfile, "w")) == nullptr)
    {
        sprintf(pdbfile, "%ctmp%cfilterXXXXXX", DIR_SEPARATOR, DIR_SEPARATOR);
        gmx_tmpnam(pdbfile);
        if ((tmpf = fopen(pdbfile, "w")) == nullptr)
        {
            gmx_fatal(FARGS, "Can not open tmp file %s", pdbfile);
        }
    }
    else
    {
        fclose(tmpf);
    }

    std::strcpy(tmpfile, "ddXXXXXX");
    gmx_tmpnam(tmpfile);
    if ((tmpf = fopen(tmpfile, "w")) == nullptr)
    {
        sprintf(tmpfile, "%ctmp%cfilterXXXXXX", DIR_SEPARATOR, DIR_SEPARATOR);
        gmx_tmpnam(tmpfile);
        if ((tmpf = fopen(tmpfile, "w")) == nullptr)
        {
            gmx_fatal(FARGS, "Can not open tmp file %s", tmpfile);
        }
    }
    else
    {
        fclose(tmpf);
    }

    if ((dptr = getenv("DSSP")) == nullptr)
    {
        dptr = "/usr/local/bin/dssp";
    }
    if (!gmx_fexist(dptr))
    {
        gmx_fatal(FARGS, "DSSP executable (%s) does not exist (use setenv DSSP)",
                  dptr);
    }
    if (dsspVersion >= 2)
    {
        if (dsspVersion > 2)
        {
            printf("\nWARNING: You use DSSP version %d, which is not explicitly\nsupported by do_dssp. Assuming version 2 syntax.\n\n", dsspVersion);
        }

        sprintf(dssp, "%s -i %s -o %s > /dev/null %s",
                dptr, pdbfile, tmpfile, bVerbose ? "" : "2> /dev/null");
    }
    else
    {
        sprintf(dssp, "%s %s %s %s > /dev/null %s",
                dptr, bDoAccSurf ? "" : "-na", pdbfile, tmpfile, bVerbose ? "" : "2> /dev/null");

    }
    fprintf(stderr, "dssp cmd='%s'\n", dssp);

    if (fnTArea)
    {
        fTArea = xvgropen(fnTArea, "Solvent Accessible Surface Area",
                          output_env_get_xvgr_tlabel(oenv), "Area (nm\\S2\\N)", oenv);
        xvgr_legend(fTArea, 2, leg, oenv);
    }
    else
    {
        fTArea = nullptr;
    }

    mat.map  = nullptr;
    mat.nmap = readcmap(opt2fn("-map", NFILE, fnm), &(mat.map));

    natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);
    if (natoms > atoms->nr)
    {
        gmx_fatal(FARGS, "\nTrajectory does not match topology!");
    }
    if (gnx > natoms)
    {
        gmx_fatal(FARGS, "\nTrajectory does not match selected group!");
    }

    snew(average_area, atoms->nres);
    snew(av_area, atoms->nres);
    snew(norm_av_area, atoms->nres);
    accr  = nullptr;
    naccr = 0;

    gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms);
    do
    {
        t = output_env_conv_time(oenv, t);
        if (bDoAccSurf && nframe >= naccr)
        {
            naccr += 10;
            srenew(accr, naccr);
            for (i = naccr-10; i < naccr; i++)
            {
                snew(accr[i], 2*atoms->nres-1);
            }
        }
        gmx_rmpbc(gpbc, natoms, box, x);
        tapein = gmx_ffopen(pdbfile, "w");
        write_pdbfile_indexed(tapein, nullptr, atoms, x, ePBC, box, ' ', -1, gnx, index, nullptr, TRUE);
        gmx_ffclose(tapein);

        if (0 != system(dssp))
        {
            gmx_fatal(FARGS, "Failed to execute command: %s\n",
                      "Try specifying your dssp version with the -ver option.", dssp);
        }

        /* strip_dssp returns the number of lines found in the dssp file, i.e.
         * the number of residues plus the separator lines */

        if (bDoAccSurf)
        {
            accr_ptr = accr[nframe];
        }

        nres_plus_separators = strip_dssp(tmpfile, nres, bPhbres, t,
                                          accr_ptr, fTArea, &mat, average_area, oenv);
        remove(tmpfile);
        remove(pdbfile);
        nframe++;
    }
    while (read_next_x(oenv, status, &t, x, box));
    fprintf(stderr, "\n");
    close_trj(status);
    if (fTArea)
    {
        xvgrclose(fTArea);
    }
    gmx_rmpbc_done(gpbc);

    prune_ss_legend(&mat);

    ss        = opt2FILE("-o", NFILE, fnm, "w");
    mat.flags = 0;
    write_xpm_m(ss, mat);
    gmx_ffclose(ss);

    if (opt2bSet("-ssdump", NFILE, fnm))
    {
        ss = opt2FILE("-ssdump", NFILE, fnm, "w");
        snew(ss_str, nres+1);
        fprintf(ss, "%d\n", nres);
        for (j = 0; j < mat.nx; j++)
        {
            for (i = 0; (i < mat.ny); i++)
            {
                ss_str[i] = mat.map[mat.matrix[j][i]].code.c1;
            }
            ss_str[i] = '\0';
            fprintf(ss, "%s\n", ss_str);
        }
        gmx_ffclose(ss);
        sfree(ss_str);
    }
    analyse_ss(fnSCount, &mat, ss_string, oenv);

    if (bDoAccSurf)
    {
        write_sas_mat(fnArea, accr, nframe, nres_plus_separators, &mat);

        for (i = 0; i < atoms->nres; i++)
        {
            av_area[i] = (average_area[i] / static_cast<real>(nframe));
        }

        norm_acc(atoms, nres, av_area, norm_av_area);

        if (fnAArea)
        {
            acc = xvgropen(fnAArea, "Average Accessible Area",
                           "Residue", "A\\S2", oenv);
            for (i = 0; (i < nres); i++)
            {
                fprintf(acc, "%5d  %10g %10g\n", i+1, av_area[i], norm_av_area[i]);
            }
            xvgrclose(acc);
        }
    }

    view_all(oenv, NFILE, fnm);

    return 0;
}

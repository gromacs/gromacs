/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2013,2014,2015,2017,2018, by the GROMACS development team, led by
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
#include <sstream>
#include <istream>
#include <iostream>


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
#include "gromacs/gmxana/dssp/structure.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/dir_separator.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/residuetypes.h"

#include <string>

namespace
{

int strip_dssp(const std::string  &dsspString,
               const gmx_bool bPhobres[], real t,
               real *acc, FILE *fTArea,
               t_matrix *mat, int average_area[],
               const gmx_output_env_t *oenv)
{
    std::istringstream is {
        dsspString
    };
    static gmx_bool bFirst = TRUE;

    std::string     ssbuf = {};
    static int      xsize;
    static int      frame;
    int             iacc = {};
    int             nresidues = {};
    int             naccf = {};
    int             naccb = {}; /* Count hydrophobic and hydrophilic residues */
    real            iaccf = {};
    real            iaccb = {};
    t_xpmelmt       c = {};

    for (int nr = 0; !is.eof(); nr++)
    {
        std::string line;
        std::getline(is, line);
        if (line[13] == '!') /* Chain separator line has '!' at pos. 13 */
        {
            ssbuf += "=";    /* Chain separator sign '=' */
        }
        else
        {
            if (line[16] == ' ')
            {
                ssbuf += "~";
            }
            else
            {
                ssbuf += line[16];
            }
        }

        line[39] = '\0';

        /* Only calculate solvent accessible area if needed */
        if ((nullptr != acc) && (line[13] != '!'))
        {
            sscanf(&(line.c_str()[34]), "%d", &iacc);
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
        mat->ny        = ssbuf.size();
        snew(mat->axis_y, ssbuf.size());
        for (size_t i = 0; i < ssbuf.size(); i++)
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
    snew(mat->matrix[frame], ssbuf.size());
    c.c2 = 0;
    for (size_t i = 0; i < ssbuf.size(); i++)
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

    /* Return the number of lines found in the dssp file (i.e. number
     * of redidues plus chain separator lines).
     * This is the number of y elements needed for the area xpm file */
    return ssbuf.size();
}
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
                     const real av_area[], real norm_av_area[])
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

static void prune_ss_legend(t_matrix *mat)
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

static void write_sas_mat(const char *fn, real **accr, int nframe, int nres, t_matrix *mat)
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

static void analyse_ss(const char *outfile, t_matrix *mat, const char *ss_string,
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
    for (s = 0; s < static_cast<size_t>(mat->nmap); s++)
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
    for (s = 0; s < static_cast<size_t>(mat->nmap); s++)
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

namespace
{
static void xlate_atomname_gmx2pdb(char *name)
{
    int  i, length;
    char temp;

    length = std::strlen(name);
    if (length > 3 && std::isdigit(name[length-1]))
    {
        temp = name[length-1];
        for (i = length-1; i > 0; --i)
        {
            name[i] = name[i-1];
        }
        name[0] = temp;
    }
}
static const char *pdbtp[epdbNR] = {
    "ATOM  ", "HETATM", "ANISOU", "CRYST1",
    "COMPND", "MODEL", "ENDMDL", "TER", "HEADER", "TITLE", "REMARK",
    "CONECT"
};

std::string pdbAtomLine(enum PDB_record   record,
                        int               atom_seq_number,
                        const char *      atom_name,
                        char              alternate_location,
                        const char *      res_name,
                        char              chain_id,
                        int               res_seq_number,
                        char              res_insertion_code,
                        real              x,
                        real              y,
                        real              z,
                        real              occupancy,
                        real              b_factor,
                        const char *      element)
{
    char     tmp_atomname[6], tmp_resname[6];
    gmx_bool start_name_in_col13;

    if (record != epdbATOM && record != epdbHETATM)
    {
        gmx_fatal(FARGS, "Can only print PDB atom lines as ATOM or HETATM records");
    }

    /* Format atom name */
    if (atom_name != nullptr)
    {
        /* If the atom name is an element name with two chars, it should start already in column 13.
         * Otherwise it should start in column 14, unless the name length is 4 chars.
         */
        if ( (element != nullptr) && (std::strlen(element) >= 2) && (gmx_strncasecmp(atom_name, element, 2) == 0) )
        {
            start_name_in_col13 = TRUE;
        }
        else
        {
            start_name_in_col13 = (std::strlen(atom_name) >= 4);
        }
        snprintf(tmp_atomname, sizeof(tmp_atomname), start_name_in_col13 ? "" : " ");
        std::strncat(tmp_atomname, atom_name, 4);
        tmp_atomname[5] = '\0';
    }
    else
    {
        tmp_atomname[0] = '\0';
    }

    /* Format residue name */
    std::strncpy(tmp_resname, (res_name != nullptr) ? res_name : "", 4);
    /* Make sure the string is terminated if strlen was > 4 */
    tmp_resname[4] = '\0';
    /* String is properly terminated, so now we can use strcat. By adding a
     * space we can write it right-justified, and if the original name was
     * three characters or less there will be a space added on the right side.
     */
    std::strcat(tmp_resname, " ");

    /* Truncate integers so they fit */
    atom_seq_number = atom_seq_number % 100000;
    res_seq_number  = res_seq_number % 10000;

    char line[STRLEN];
    sprintf(line, "%-6s%5d %-4.4s%c%4.4s%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
            pdbtp[record],
            atom_seq_number,
            tmp_atomname,
            alternate_location,
            tmp_resname,
            chain_id,
            res_seq_number,
            res_insertion_code,
            x, y, z,
            occupancy,
            b_factor,
            (element != nullptr) ? element : "");

    return line;
}

std::string AtomsAsPdbString(const t_atoms *atoms, const rvec x[], int nindex, const int index[] )
{
    std::string       pdbString;
    char              resnm[6], nm[6];
    int               i, ii;
    int               resind, resnr;
    enum PDB_record   type;
    unsigned char     resic, ch;
    char              altloc;
    real              occup, bfac;
    gmx_bool          bOccup;
    int               chainnum, lastchainnum;
    gmx_residuetype_t*rt;
    const char       *p_restype;
    const char       *p_lastrestype;

    gmx_residuetype_init(&rt);
    pdbString += "MODEL\n";
    if (atoms->havePdbInfo)
    {
        /* Check whether any occupancies are set, in that case leave it as is,
         * otherwise set them all to one
         */
        bOccup = TRUE;
        for (ii = 0; (ii < nindex) && bOccup; ii++)
        {
            i      = index[ii];
            bOccup = bOccup && (atoms->pdbinfo[i].occup == 0.0);
        }
    }
    else
    {
        bOccup = FALSE;
    }

    lastchainnum      = -1;
    p_restype         = nullptr;

    for (ii = 0; ii < nindex; ii++)
    {
        i             = index[ii];
        resind        = atoms->atom[i].resind;
        chainnum      = atoms->resinfo[resind].chainnum;
        p_lastrestype = p_restype;
        gmx_residuetype_get_type(rt, *atoms->resinfo[resind].name, &p_restype);

        /* Add a TER record if we changed chain, and if either the previous or this chain is protein/DNA/RNA. */
        if (ii > 0 && chainnum != lastchainnum)
        {
            /* Only add TER if the previous chain contained protein/DNA/RNA. */
            if (gmx_residuetype_is_protein(rt, p_lastrestype) || gmx_residuetype_is_dna(rt, p_lastrestype) || gmx_residuetype_is_rna(rt, p_lastrestype))
            {
                pdbString += "TER\n";
            }
            lastchainnum    = chainnum;
        }

        strncpy(resnm, *atoms->resinfo[resind].name, sizeof(resnm)-1);
        resnm[sizeof(resnm)-1] = 0;
        strncpy(nm, *atoms->atomname[i], sizeof(nm)-1);
        nm[sizeof(nm)-1] = 0;

        /* rename HG12 to 2HG1, etc. */
        xlate_atomname_gmx2pdb(nm);
        resnr = atoms->resinfo[resind].nr;
        resic = atoms->resinfo[resind].ic;
        ch    = atoms->resinfo[resind].chainid;

        if (ch == 0)
        {
            ch = ' ';
        }
        if (resnr >= 10000)
        {
            resnr = resnr % 10000;
        }
        t_pdbinfo pdbinfo;
        if (atoms->pdbinfo != nullptr)
        {
            pdbinfo = atoms->pdbinfo[i];
        }
        else
        {
            gmx_pdbinfo_init_default(&pdbinfo);
        }
        type   = static_cast<enum PDB_record>(pdbinfo.type);
        altloc = pdbinfo.altloc;
        if (!isalnum(altloc))
        {
            altloc = ' ';
        }
        occup      = bOccup ? 1.0 : pdbinfo.occup;
        bfac       = pdbinfo.bfac;
        pdbString += pdbAtomLine(type, i+1, nm, altloc, resnm, ch, resnr, resic,
                                 10*x[i][XX], 10*x[i][YY], 10*x[i][ZZ], occup, bfac, atoms->atom[i].elem);
    }

    pdbString += "TER\n";
    pdbString += "ENDMDL\n";

    gmx_residuetype_destroy(rt);
    return pdbString;
}
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
    FILE              *ss, *acc, *fTArea;
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
    bDoAccSurf = ((fnArea != nullptr) || (fnTArea != nullptr) || (fnAArea != nullptr));

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

        std::string pdbString = AtomsAsPdbString(atoms, x, gnx, index);
        MProtein    protein;
        protein.ReadPDB(pdbString);
        protein.CalculateSecondaryStructure();
        std::string dsspString = WriteDSSP(protein);
        if (bDoAccSurf)
        {
            accr_ptr = accr[nframe];
        }
        /* strip_dssp returns the number of lines found in the dssp file, i.e.
         * the number of residues plus the separator lines */
        nres_plus_separators = strip_dssp(dsspString, bPhbres, t,
                                          accr_ptr, fTArea, &mat, average_area, oenv);
        nframe++;
    }
    while (read_next_x(oenv, status, &t, x, box));
    fprintf(stderr, "\n");
    close_trx(status);
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

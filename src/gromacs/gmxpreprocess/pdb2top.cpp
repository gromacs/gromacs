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

#include "pdb2top.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxpreprocess/add_par.h"
#include "gromacs/gmxpreprocess/fflibutil.h"
#include "gromacs/gmxpreprocess/gen_ad.h"
#include "gromacs/gmxpreprocess/gen_vsite.h"
#include "gromacs/gmxpreprocess/gpp_nextnb.h"
#include "gromacs/gmxpreprocess/h_db.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/pgutil.h"
#include "gromacs/gmxpreprocess/resall.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/topio.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/binaryinformation.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/dir_separator.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/niceheader.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

/* this must correspond to enum in pdb2top.h */
const char *hh[ehisNR]   = { "HISD", "HISE", "HISH", "HIS1" };

static int missing_atoms(t_restp *rp, int resind, t_atoms *at, int i0, int i)
{
    int      j, k, nmiss;
    char    *name;
    bool     bFound;

    nmiss = 0;
    for (j = 0; j < rp->natom; j++)
    {
        name   = *(rp->atomname[j]);
        bFound = FALSE;
        for (k = i0; k < i; k++)
        {
            bFound = (bFound || (gmx_strcasecmp(*(at->atomname[k]), name) == 0));
        }
        if (!bFound)
        {
            nmiss++;
            fprintf(stderr, "\nWARNING: "
                    "atom %s is missing in residue %s %d in the pdb file\n",
                    name, *(at->resinfo[resind].name), at->resinfo[resind].nr);
            if (name[0] == 'H' || name[0] == 'h')
            {
                fprintf(stderr, "         You might need to add atom %s to the hydrogen database of building block %s\n"
                        "         in the file %s.hdb (see the manual)\n",
                        name, *(at->resinfo[resind].rtp), rp->filebase);
            }
            fprintf(stderr, "\n");
        }
    }

    return nmiss;
}

bool is_int(double x)
{
    const double tol = 1e-4;
    int          ix;

    if (x < 0)
    {
        x = -x;
    }
    ix = gmx::roundToInt(x);

    return (fabs(x-ix) < tol);
}

static void
choose_ff_impl(const char *ffsel,
               char *forcefield, int ff_maxlen,
               char *ffdir, int ffdir_maxlen)
{
    std::vector<gmx::DataFileInfo> ffdirs = fflib_enumerate_forcefields();
    const int nff = static_cast<int>(ffdirs.size());

    /* Replace with unix path separators */
#if DIR_SEPARATOR != '/'
    for (int i = 0; i < nff; ++i)
    {
        std::replace(ffdirs[i].dir.begin(), ffdirs[i].dir.end(), DIR_SEPARATOR, '/');
    }
#endif

    /* Store the force field names in ffs */
    std::vector<std::string> ffs;
    ffs.reserve(ffdirs.size());
    for (int i = 0; i < nff; ++i)
    {
        ffs.push_back(gmx::stripSuffixIfPresent(ffdirs[i].name,
                                                fflib_forcefield_dir_ext()));
    }

    int sel;
    if (ffsel != nullptr)
    {
        sel         = -1;
        int cwdsel  = -1;
        int nfound  = 0;
        for (int i = 0; i < nff; ++i)
        {
            if (ffs[i] == ffsel)
            {
                /* Matching ff name */
                sel = i;
                nfound++;

                if (ffdirs[i].dir == ".")
                {
                    cwdsel = i;
                }
            }
        }

        if (cwdsel != -1)
        {
            sel = cwdsel;
        }

        if (nfound > 1)
        {
            if (cwdsel != -1)
            {
                fprintf(stderr,
                        "Force field '%s' occurs in %d places. pdb2gmx is using the one in the\n"
                        "current directory. Use interactive selection (not the -ff option) if\n"
                        "you would prefer a different one.\n", ffsel, nfound);
            }
            else
            {
                std::string message = gmx::formatString(
                            "Force field '%s' occurs in %d places, but not in "
                            "the current directory.\n"
                            "Run without the -ff switch and select the force "
                            "field interactively.", ffsel, nfound);
                GMX_THROW(gmx::InconsistentInputError(message));
            }
        }
        else if (nfound == 0)
        {
            std::string message = gmx::formatString(
                        "Could not find force field '%s' in current directory, "
                        "install tree or GMXLIB path.", ffsel);
            GMX_THROW(gmx::InconsistentInputError(message));
        }
    }
    else if (nff > 1)
    {
        std::vector<std::string> desc;
        desc.reserve(ffdirs.size());
        for (int i = 0; i < nff; ++i)
        {
            std::string docFileName(
                    gmx::Path::join(ffdirs[i].dir, ffdirs[i].name,
                                    fflib_forcefield_doc()));
            // TODO: Just try to open the file with a method that does not
            // throw/bail out with a fatal error instead of multiple checks.
            if (gmx::File::exists(docFileName, gmx::File::returnFalseOnError))
            {
                // TODO: Use a C++ API without such an intermediate/fixed-length buffer.
                char  buf[STRLEN];
                /* We don't use fflib_open, because we don't want printf's */
                FILE *fp = gmx_ffopen(docFileName, "r");
                get_a_line(fp, buf, STRLEN);
                gmx_ffclose(fp);
                desc.emplace_back(buf);
            }
            else
            {
                desc.push_back(ffs[i]);
            }
        }
        /* Order force fields from the same dir alphabetically
         * and put deprecated force fields at the end.
         */
        for (int i = 0; i < nff; ++i)
        {
            for (int j = i + 1; j < nff; ++j)
            {
                if (ffdirs[i].dir == ffdirs[j].dir &&
                    ((desc[i][0] == '[' && desc[j][0] != '[') ||
                     ((desc[i][0] == '[' || desc[j][0] != '[') &&
                      gmx_strcasecmp(desc[i].c_str(), desc[j].c_str()) > 0)))
                {
                    std::swap(ffdirs[i].name, ffdirs[j].name);
                    std::swap(ffs[i], ffs[j]);
                    std::swap(desc[i], desc[j]);
                }
            }
        }

        printf("\nSelect the Force Field:\n");
        for (int i = 0; i < nff; ++i)
        {
            if (i == 0 || ffdirs[i-1].dir != ffdirs[i].dir)
            {
                if (ffdirs[i].dir == ".")
                {
                    printf("From current directory:\n");
                }
                else
                {
                    printf("From '%s':\n", ffdirs[i].dir.c_str());
                }
            }
            printf("%2d: %s\n", i+1, desc[i].c_str());
        }

        sel = -1;
        // TODO: Add a C++ API for this.
        char   buf[STRLEN];
        char  *pret;
        do
        {
            pret = fgets(buf, STRLEN, stdin);

            if (pret != nullptr)
            {
                sel = strtol(buf, nullptr, 10);
                sel--;
            }
        }
        while (pret == nullptr || (sel < 0) || (sel >= nff));

        /* Check for a current limitation of the fflib code.
         * It will always read from the first ff directory in the list.
         * This check assumes that the order of ffs matches the order
         * in which fflib_open searches ff library files.
         */
        for (int i = 0; i < sel; i++)
        {
            if (ffs[i] == ffs[sel])
            {
                std::string message = gmx::formatString(
                            "Can only select the first of multiple force "
                            "field entries with directory name '%s%s' in "
                            "the list. If you want to use the next entry, "
                            "run pdb2gmx in a different directory, set GMXLIB "
                            "to point to the desired force field first, and/or "
                            "rename or move the force field directory present "
                            "in the current working directory.",
                            ffs[sel].c_str(), fflib_forcefield_dir_ext());
                GMX_THROW(gmx::NotImplementedError(message));
            }
        }
    }
    else
    {
        sel = 0;
    }

    if (ffs[sel].length() >= static_cast<size_t>(ff_maxlen))
    {
        std::string message = gmx::formatString(
                    "Length of force field name (%d) >= maxlen (%d)",
                    static_cast<int>(ffs[sel].length()), ff_maxlen);
        GMX_THROW(gmx::InvalidInputError(message));
    }
    strcpy(forcefield, ffs[sel].c_str());

    std::string ffpath;
    if (ffdirs[sel].bFromDefaultDir)
    {
        ffpath = ffdirs[sel].name;
    }
    else
    {
        ffpath = gmx::Path::join(ffdirs[sel].dir, ffdirs[sel].name);
    }
    if (ffpath.length() >= static_cast<size_t>(ffdir_maxlen))
    {
        std::string message = gmx::formatString(
                    "Length of force field dir (%d) >= maxlen (%d)",
                    static_cast<int>(ffpath.length()), ffdir_maxlen);
        GMX_THROW(gmx::InvalidInputError(message));
    }
    strcpy(ffdir, ffpath.c_str());
}

void
choose_ff(const char *ffsel,
          char *forcefield, int ff_maxlen,
          char *ffdir, int ffdir_maxlen)
{
    try
    {
        choose_ff_impl(ffsel, forcefield, ff_maxlen, ffdir, ffdir_maxlen);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

void choose_watermodel(const char *wmsel, const char *ffdir,
                       char **watermodel)
{
    const char *fn_watermodels = "watermodels.dat";
    FILE       *fp;
    char        buf[STRLEN];
    int         nwm, sel, i;
    char      **model;
    char       *pret;

    if (strcmp(wmsel, "none") == 0)
    {
        *watermodel = nullptr;

        return;
    }
    else if (strcmp(wmsel, "select") != 0)
    {
        *watermodel = gmx_strdup(wmsel);

        return;
    }

    std::string filename = gmx::Path::join(ffdir, fn_watermodels);
    if (!fflib_fexist(filename))
    {
        fprintf(stderr, "No file '%s' found, will not include a water model\n",
                fn_watermodels);
        *watermodel = nullptr;

        return;
    }

    fp = fflib_open(filename);
    printf("\nSelect the Water Model:\n");
    nwm   = 0;
    model = nullptr;
    while (get_a_line(fp, buf, STRLEN))
    {
        srenew(model, nwm+1);
        snew(model[nwm], STRLEN);
        sscanf(buf, "%s%n", model[nwm], &i);
        if (i > 0)
        {
            ltrim(buf+i);
            fprintf(stderr, "%2d: %s\n", nwm+1, buf+i);
            nwm++;
        }
        else
        {
            sfree(model[nwm]);
        }
    }
    gmx_ffclose(fp);
    fprintf(stderr, "%2d: %s\n", nwm+1, "None");

    sel = -1;
    do
    {
        pret = fgets(buf, STRLEN, stdin);

        if (pret != nullptr)
        {
            sel = strtol(buf, nullptr, 10);
            sel--;
        }
    }
    while (pret == nullptr || sel < 0 || sel > nwm);

    if (sel == nwm)
    {
        *watermodel = nullptr;
    }
    else
    {
        *watermodel = gmx_strdup(model[sel]);
    }

    for (i = 0; i < nwm; i++)
    {
        sfree(model[i]);
    }
    sfree(model);
}

static int name2type(t_atoms *at, int **cgnr,
                     t_restp restp[], gmx_residuetype_t *rt)
{
    int         i, j, prevresind, resind, i0, prevcg, cg, curcg;
    char       *name;
    bool        bNterm;
    double      qt;
    int         nmissat;

    nmissat = 0;

    resind = -1;
    bNterm = FALSE;
    i0     = 0;
    snew(*cgnr, at->nr);
    qt     = 0;
    curcg  = 0;
    cg     = -1;

    for (i = 0; (i < at->nr); i++)
    {
        prevresind = resind;
        if (at->atom[i].resind != resind)
        {
            bool bProt;
            resind = at->atom[i].resind;
            bProt  = gmx_residuetype_is_protein(rt, *(at->resinfo[resind].name));
            bNterm = bProt && (resind == 0);
            if (resind > 0)
            {
                nmissat += missing_atoms(&restp[prevresind], prevresind, at, i0, i);
            }
            i0 = i;
        }
        if (at->atom[i].m == 0)
        {
            qt               = 0;
            prevcg           = cg;
            name             = *(at->atomname[i]);
            j                = search_jtype(&restp[resind], name, bNterm);
            at->atom[i].type = restp[resind].atom[j].type;
            at->atom[i].q    = restp[resind].atom[j].q;
            at->atom[i].m    = restp[resind].atom[j].m;
            cg               = restp[resind].cgnr[j];
            /* A charge group number -1 signals a separate charge group
             * for this atom.
             */
            if ( (cg == -1) || (cg != prevcg) || (resind != prevresind) )
            {
                curcg++;
            }
        }
        else
        {
            cg = -1;
            if (is_int(qt))
            {
                qt = 0;
                curcg++;
            }
            qt += at->atom[i].q;
        }
        (*cgnr)[i]        = curcg;
        at->atom[i].typeB = at->atom[i].type;
        at->atom[i].qB    = at->atom[i].q;
        at->atom[i].mB    = at->atom[i].m;
    }
    nmissat += missing_atoms(&restp[resind], resind, at, i0, i);

    return nmissat;
}

static void print_top_heavy_H(FILE *out, real mHmult)
{
    if (mHmult == 2.0)
    {
        fprintf(out, "; Using deuterium instead of hydrogen\n\n");
    }
    else if (mHmult == 4.0)
    {
        fprintf(out, "#define HEAVY_H\n\n");
    }
    else if (mHmult != 1.0)
    {
        fprintf(stderr, "WARNING: unsupported proton mass multiplier (%g) "
                "in pdb2top\n", mHmult);
    }
}

void print_top_comment(FILE       *out,
                       const char *filename,
                       const char *ffdir,
                       bool        bITP)
{
    char  ffdir_parent[STRLEN];
    char *p;

    try
    {
        gmx::TextWriter writer(out);
        gmx::niceHeader(&writer, filename, ';');
        writer.writeLine(gmx::formatString(";\tThis is a %s topology file", bITP ? "include" : "standalone"));
        writer.writeLine(";");

        gmx::BinaryInformationSettings settings;
        settings.generatedByHeader(true);
        settings.linePrefix(";\t");
        gmx::printBinaryInformation(&writer, gmx::getProgramContext(), settings);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    if (strchr(ffdir, '/') == nullptr)
    {
        fprintf(out, ";\tForce field was read from the standard GROMACS share directory.\n;\n\n");
    }
    else if (ffdir[0] == '.')
    {
        fprintf(out, ";\tForce field was read from current directory or a relative path - path added.\n;\n\n");
    }
    else
    {
        strncpy(ffdir_parent, ffdir, STRLEN-1);
        ffdir_parent[STRLEN-1] = '\0'; /*make sure it is 0-terminated even for long string*/
        p = strrchr(ffdir_parent, '/');

        *p = '\0';

        fprintf(out,
                ";\tForce field data was read from:\n"
                ";\t%s\n"
                ";\n"
                ";\tNote:\n"
                ";\tThis might be a non-standard force field location. When you use this topology, the\n"
                ";\tforce field must either be present in the current directory, or the location\n"
                ";\tspecified in the GMXLIB path variable or with the 'include' mdp file option.\n;\n\n",
                ffdir_parent);
    }
}

void print_top_header(FILE *out, const char *filename,
                      bool bITP, const char *ffdir, real mHmult)
{
    const char *p;

    print_top_comment(out, filename, ffdir, bITP);

    print_top_heavy_H(out, mHmult);
    fprintf(out, "; Include forcefield parameters\n");

    p = strrchr(ffdir, '/');
    p = (ffdir[0] == '.' || p == nullptr) ? ffdir : p+1;

    fprintf(out, "#include \"%s/%s\"\n\n", p, fflib_forcefield_itp());
}

static void print_top_posre(FILE *out, const char *pr)
{
    fprintf(out, "; Include Position restraint file\n");
    fprintf(out, "#ifdef POSRES\n");
    fprintf(out, "#include \"%s\"\n", pr);
    fprintf(out, "#endif\n\n");
}

static void print_top_water(FILE *out, const char *ffdir, const char *water)
{
    const char *p;
    char        buf[STRLEN];

    fprintf(out, "; Include water topology\n");

    p = strrchr(ffdir, '/');
    p = (ffdir[0] == '.' || p == nullptr) ? ffdir : p+1;
    fprintf(out, "#include \"%s/%s.itp\"\n", p, water);

    fprintf(out, "\n");
    fprintf(out, "#ifdef POSRES_WATER\n");
    fprintf(out, "; Position restraint for each water oxygen\n");
    fprintf(out, "[ position_restraints ]\n");
    fprintf(out, ";%3s %5s %9s %10s %10s\n", "i", "funct", "fcx", "fcy", "fcz");
    fprintf(out, "%4d %4d %10g %10g %10g\n", 1, 1, 1000.0, 1000.0, 1000.0);
    fprintf(out, "#endif\n");
    fprintf(out, "\n");

    sprintf(buf, "%s/ions.itp", p);

    if (fflib_fexist(buf))
    {
        fprintf(out, "; Include topology for ions\n");
        fprintf(out, "#include \"%s\"\n", buf);
        fprintf(out, "\n");
    }
}

static void print_top_system(FILE *out, const char *title)
{
    fprintf(out, "[ %s ]\n", dir2str(d_system));
    fprintf(out, "; Name\n");
    fprintf(out, "%s\n\n", title[0] ? title : "Protein");
}

void print_top_mols(FILE *out,
                    const char *title, const char *ffdir, const char *water,
                    int nincl, char **incls, int nmol, t_mols *mols)
{

    if (nincl > 0)
    {
        fprintf(out, "; Include chain topologies\n");
        for (int i = 0; i < nincl; i++)
        {
            fprintf(out, "#include \"%s\"\n", gmx::Path::getFilename(incls[i]).c_str());
        }
        fprintf(out, "\n");
    }

    if (water)
    {
        print_top_water(out, ffdir, water);
    }
    print_top_system(out, title);

    if (nmol)
    {
        fprintf(out, "[ %s ]\n", dir2str(d_molecules));
        fprintf(out, "; %-15s %5s\n", "Compound", "#mols");
        for (int i = 0; i < nmol; i++)
        {
            fprintf(out, "%-15s %5d\n", mols[i].name, mols[i].nr);
        }
    }
}

void write_top(FILE *out, const char *pr, const char *molname,
               t_atoms *at, bool bRTPresname,
               int bts[], t_params plist[], t_excls excls[],
               gpp_atomtype_t atype, int *cgnr, int nrexcl)
/* NOTE: nrexcl is not the size of *excl! */
{
    if (at && atype && cgnr)
    {
        fprintf(out, "[ %s ]\n", dir2str(d_moleculetype));
        fprintf(out, "; %-15s %5s\n", "Name", "nrexcl");
        fprintf(out, "%-15s %5d\n\n", molname ? molname : "Protein", nrexcl);

        print_atoms(out, atype, at, cgnr, bRTPresname);
        print_bondeds(out, at->nr, d_bonds,      F_BONDS,    bts[ebtsBONDS], plist);
        print_bondeds(out, at->nr, d_constraints, F_CONSTR,   0,              plist);
        print_bondeds(out, at->nr, d_constraints, F_CONSTRNC, 0,              plist);
        print_bondeds(out, at->nr, d_pairs,      F_LJ14,     0,              plist);
        print_excl(out, at->nr, excls);
        print_bondeds(out, at->nr, d_angles,     F_ANGLES,   bts[ebtsANGLES], plist);
        print_bondeds(out, at->nr, d_dihedrals,  F_PDIHS,    bts[ebtsPDIHS], plist);
        print_bondeds(out, at->nr, d_dihedrals,  F_IDIHS,    bts[ebtsIDIHS], plist);
        print_bondeds(out, at->nr, d_cmap,       F_CMAP,     bts[ebtsCMAP],  plist);
        print_bondeds(out, at->nr, d_polarization, F_POLARIZATION,   0,       plist);
        print_bondeds(out, at->nr, d_thole_polarization, F_THOLE_POL, 0,       plist);
        print_bondeds(out, at->nr, d_vsites2,    F_VSITE2,   0,              plist);
        print_bondeds(out, at->nr, d_vsites3,    F_VSITE3,   0,              plist);
        print_bondeds(out, at->nr, d_vsites3,    F_VSITE3FD, 0,              plist);
        print_bondeds(out, at->nr, d_vsites3,    F_VSITE3FAD, 0,              plist);
        print_bondeds(out, at->nr, d_vsites3,    F_VSITE3OUT, 0,              plist);
        print_bondeds(out, at->nr, d_vsites4,    F_VSITE4FD, 0,              plist);
        print_bondeds(out, at->nr, d_vsites4,    F_VSITE4FDN, 0,             plist);

        if (pr)
        {
            print_top_posre(out, pr);
        }
    }
}



static void do_ssbonds(t_params *ps, t_atoms *atoms,
                       int nssbonds, t_ssbond *ssbonds, bool bAllowMissing)
{
    int     i, ri, rj;
    int     ai, aj;

    for (i = 0; (i < nssbonds); i++)
    {
        ri = ssbonds[i].res1;
        rj = ssbonds[i].res2;
        ai = search_res_atom(ssbonds[i].a1, ri, atoms,
                             "special bond", bAllowMissing);
        aj = search_res_atom(ssbonds[i].a2, rj, atoms,
                             "special bond", bAllowMissing);
        if ((ai == -1) || (aj == -1))
        {
            gmx_fatal(FARGS, "Trying to make impossible special bond (%s-%s)!",
                      ssbonds[i].a1, ssbonds[i].a2);
        }
        add_param(ps, ai, aj, nullptr, nullptr);
    }
}

static void at2bonds(t_params *psb, t_hackblock *hb,
                     t_atoms *atoms,
                     rvec x[],
                     real long_bond_dist, real short_bond_dist)
{
    int         resind, i, j, k;
    int         ai, aj;
    real        dist2, long_bond_dist2, short_bond_dist2;
    const char *ptr;

    long_bond_dist2  = gmx::square(long_bond_dist);
    short_bond_dist2 = gmx::square(short_bond_dist);

    if (debug)
    {
        ptr = "bond";
    }
    else
    {
        ptr = "check";
    }

    fprintf(stderr, "Making bonds...\n");
    i = 0;
    for (resind = 0; (resind < atoms->nres) && (i < atoms->nr); resind++)
    {
        /* add bonds from list of bonded interactions */
        for (j = 0; j < hb[resind].rb[ebtsBONDS].nb; j++)
        {
            /* Unfortunately we can not issue errors or warnings
             * for missing atoms in bonds, as the hydrogens and terminal atoms
             * have not been added yet.
             */
            ai = search_atom(hb[resind].rb[ebtsBONDS].b[j].a[0], i, atoms,
                             ptr, TRUE);
            aj = search_atom(hb[resind].rb[ebtsBONDS].b[j].a[1], i, atoms,
                             ptr, TRUE);
            if (ai != -1 && aj != -1)
            {
                dist2 = distance2(x[ai], x[aj]);
                if (dist2 > long_bond_dist2)

                {
                    fprintf(stderr, "Warning: Long Bond (%d-%d = %g nm)\n",
                            ai+1, aj+1, std::sqrt(dist2));
                }
                else if (dist2 < short_bond_dist2)
                {
                    fprintf(stderr, "Warning: Short Bond (%d-%d = %g nm)\n",
                            ai+1, aj+1, std::sqrt(dist2));
                }
                add_param(psb, ai, aj, nullptr, hb[resind].rb[ebtsBONDS].b[j].s);
            }
        }
        /* add bonds from list of hacks (each added atom gets a bond) */
        while ( (i < atoms->nr) && (atoms->atom[i].resind == resind) )
        {
            for (j = 0; j < hb[resind].nhack; j++)
            {
                if ( ( hb[resind].hack[j].tp > 0 ||
                       hb[resind].hack[j].oname == nullptr ) &&
                     strcmp(hb[resind].hack[j].a[0], *(atoms->atomname[i])) == 0)
                {
                    switch (hb[resind].hack[j].tp)
                    {
                        case 9:                                         /* COOH terminus */
                            add_param(psb, i, i+1, nullptr, nullptr);   /* C-O  */
                            add_param(psb, i, i+2, nullptr, nullptr);   /* C-OA */
                            add_param(psb, i+2, i+3, nullptr, nullptr); /* OA-H */
                            break;
                        default:
                            for (k = 0; (k < hb[resind].hack[j].nr); k++)
                            {
                                add_param(psb, i, i+k+1, nullptr, nullptr);
                            }
                    }
                }
            }
            i++;
        }
        /* we're now at the start of the next residue */
    }
}

static int pcompar(const void *a, const void *b)
{
    const t_param *pa, *pb;
    int            d;
    pa = static_cast<const t_param *>(a);
    pb = static_cast<const t_param *>(b);

    d = pa->a[0] - pb->a[0];
    if (d == 0)
    {
        d = pa->a[1] - pb->a[1];
    }
    if (d == 0)
    {
        return strlen(pb->s) - strlen(pa->s);
    }
    else
    {
        return d;
    }
}

static void clean_bonds(t_params *ps)
{
    int     i, j;
    int     a;

    if (ps->nr > 0)
    {
        /* swap atomnumbers in bond if first larger than second: */
        for (i = 0; (i < ps->nr); i++)
        {
            if (ps->param[i].a[1] < ps->param[i].a[0])
            {
                a                 = ps->param[i].a[0];
                ps->param[i].a[0] = ps->param[i].a[1];
                ps->param[i].a[1] = a;
            }
        }

        /* Sort bonds */
        qsort(ps->param, ps->nr, static_cast<size_t>(sizeof(ps->param[0])), pcompar);

        /* remove doubles, keep the first one always. */
        j = 1;
        for (i = 1; (i < ps->nr); i++)
        {
            if ((ps->param[i].a[0] != ps->param[j-1].a[0]) ||
                (ps->param[i].a[1] != ps->param[j-1].a[1]) )
            {
                if (j != i)
                {
                    cp_param(&(ps->param[j]), &(ps->param[i]));
                }
                j++;
            }
        }
        fprintf(stderr, "Number of bonds was %d, now %d\n", ps->nr, j);
        ps->nr = j;
    }
    else
    {
        fprintf(stderr, "No bonds\n");
    }
}

void print_sums(t_atoms *atoms, bool bSystem)
{
    double      m, qtot;
    int         i;
    const char *where;

    if (bSystem)
    {
        where = " in system";
    }
    else
    {
        where = "";
    }

    m    = 0;
    qtot = 0;
    for (i = 0; (i < atoms->nr); i++)
    {
        m    += atoms->atom[i].m;
        qtot += atoms->atom[i].q;
    }

    fprintf(stderr, "Total mass%s %.3f a.m.u.\n", where, m);
    fprintf(stderr, "Total charge%s %.3f e\n", where, qtot);
}

static void check_restp_type(const char *name, int t1, int t2)
{
    if (t1 != t2)
    {
        gmx_fatal(FARGS, "Residues in one molecule have a different '%s' type: %d and %d", name, t1, t2);
    }
}

static void check_restp_types(t_restp *r0, t_restp *r1)
{
    int i;

    check_restp_type("all dihedrals", static_cast<int>(r0->bKeepAllGeneratedDihedrals), static_cast<int>(r1->bKeepAllGeneratedDihedrals));
    check_restp_type("nrexcl", r0->nrexcl, r1->nrexcl);
    check_restp_type("HH14", static_cast<int>(r0->bGenerateHH14Interactions), static_cast<int>(r1->bGenerateHH14Interactions));
    check_restp_type("remove dihedrals", static_cast<int>(r0->bRemoveDihedralIfWithImproper), static_cast<int>(r1->bRemoveDihedralIfWithImproper));

    for (i = 0; i < ebtsNR; i++)
    {
        check_restp_type(btsNames[i], r0->rb[i].type, r1->rb[i].type);
    }
}

static void add_atom_to_restp(t_restp *restp, int at_start, const t_hack *hack)
{
    char        buf[STRLEN];
    int         k;
    const char *Hnum = "123456";

    strcpy(buf, hack->nname);
    buf[strlen(buf)+1] = '\0';
    if (hack->nr > 1)
    {
        buf[strlen(buf)] = '-';
    }
    /* make space */
    restp->natom += hack->nr;
    srenew(restp->atom,     restp->natom);
    srenew(restp->atomname, restp->natom);
    srenew(restp->cgnr,     restp->natom);
    /* shift the rest */
    for (k = restp->natom-1; k > at_start+hack->nr; k--)
    {
        restp->atom[k] =
            restp->atom    [k - hack->nr];
        restp->atomname[k] =
            restp->atomname[k - hack->nr];
        restp->cgnr[k] =
            restp->cgnr    [k - hack->nr];
    }
    /* now add them */
    for (k = 0; k < hack->nr; k++)
    {
        /* set counter in atomname */
        if (hack->nr > 1)
        {
            buf[strlen(buf)-1] = Hnum[k];
        }
        snew( restp->atomname[at_start+1+k], 1);
        restp->atom     [at_start+1+k] = *hack->atom;
        *restp->atomname[at_start+1+k] = gmx_strdup(buf);
        if (hack->cgnr != NOTSET)
        {
            restp->cgnr   [at_start+1+k] = hack->cgnr;
        }
        else
        {
            restp->cgnr   [at_start+1+k] = restp->cgnr[at_start];
        }
    }
}

void get_hackblocks_rtp(t_hackblock **hb, t_restp **restp,
                        int nrtp, t_restp rtp[],
                        int nres, t_resinfo *resinfo,
                        int nterpairs,
                        t_hackblock **ntdb, t_hackblock **ctdb,
                        const int *rn, const int *rc,
                        bool bAllowMissing)
{
    int         i, j, k, l;
    char       *key;
    t_restp    *res;
    int         tern, terc;
    bool        bRM;

    snew(*hb, nres);
    snew(*restp, nres);
    /* first the termini */
    for (i = 0; i < nterpairs; i++)
    {
        if (rn[i] >= 0 && ntdb[i] != nullptr)
        {
            copy_t_hackblock(ntdb[i], &(*hb)[rn[i]]);
        }
        if (rc[i] >= 0 && ctdb[i] != nullptr)
        {
            merge_t_hackblock(ctdb[i], &(*hb)[rc[i]]);
        }
    }

    /* then the whole rtp */
    for (i = 0; i < nres; i++)
    {
        /* Here we allow a mismatch of one character when looking for the rtp entry.
         * For such a mismatch there should be only one mismatching name.
         * This is mainly useful for small molecules such as ions.
         * Note that this will usually not work for protein, DNA and RNA,
         * since there the residue names should be listed in residuetypes.dat
         * and an error will have been generated earlier in the process.
         */
        key = *resinfo[i].rtp;
        snew(resinfo[i].rtp, 1);
        *resinfo[i].rtp = search_rtp(key, nrtp, rtp);
        res             = get_restp(*resinfo[i].rtp, nrtp, rtp);
        copy_t_restp(res, &(*restp)[i]);

        /* Check that we do not have different bonded types in one molecule */
        check_restp_types(&(*restp)[0], &(*restp)[i]);

        tern = -1;
        for (j = 0; j < nterpairs && tern == -1; j++)
        {
            if (i == rn[j])
            {
                tern = j;
            }
        }
        terc = -1;
        for (j = 0; j < nterpairs && terc == -1; j++)
        {
            if (i == rc[j])
            {
                terc = j;
            }
        }
        bRM = merge_t_bondeds(res->rb, (*hb)[i].rb, tern >= 0, terc >= 0);

        if (bRM && ((tern >= 0 && ntdb[tern] == nullptr) ||
                    (terc >= 0 && ctdb[terc] == nullptr)))
        {
            const char *errString = "There is a dangling bond at at least one of the terminal ends and the force field does not provide terminal entries or files. Fix your terminal residues so that they match the residue database (.rtp) entries, or provide terminal database entries (.tdb).";
            if (bAllowMissing)
            {
                fprintf(stderr, "%s\n", errString);
            }
            else
            {
                gmx_fatal(FARGS, "%s", errString);
            }
        }
        else if (bRM && ((tern >= 0 && ntdb[tern]->nhack == 0) ||
                         (terc >= 0 && ctdb[terc]->nhack == 0)))
        {
            const char *errString = "There is a dangling bond at at least one of the terminal ends. Fix your coordinate file, add a new terminal database entry (.tdb), or select the proper existing terminal entry.";
            if (bAllowMissing)
            {
                fprintf(stderr, "%s\n", errString);
            }
            else
            {
                gmx_fatal(FARGS, "%s", errString);
            }
        }
    }

    /* now perform t_hack's on t_restp's,
       i.e. add's and deletes from termini database will be
       added to/removed from residue topology
       we'll do this on one big dirty loop, so it won't make easy reading! */
    for (i = 0; i < nres; i++)
    {
        for (j = 0; j < (*hb)[i].nhack; j++)
        {
            if ( (*hb)[i].hack[j].nr)
            {
                /* find atom in restp */
                for (l = 0; l < (*restp)[i].natom; l++)
                {
                    if ( ( (*hb)[i].hack[j].oname == nullptr &&
                           strcmp((*hb)[i].hack[j].a[0], *(*restp)[i].atomname[l]) == 0 ) ||
                         ( (*hb)[i].hack[j].oname != nullptr &&
                           strcmp((*hb)[i].hack[j].oname, *(*restp)[i].atomname[l]) == 0 ) )
                    {
                        break;
                    }
                }
                if (l == (*restp)[i].natom)
                {
                    /* If we are doing an atom rename only, we don't need
                     * to generate a fatal error if the old name is not found
                     * in the rtp.
                     */
                    /* Deleting can happen also only on the input atoms,
                     * not necessarily always on the rtp entry.
                     */
                    if (!((*hb)[i].hack[j].oname != nullptr &&
                          (*hb)[i].hack[j].nname != nullptr) &&
                        !((*hb)[i].hack[j].oname != nullptr &&
                          (*hb)[i].hack[j].nname == nullptr))
                    {
                        gmx_fatal(FARGS,
                                  "atom %s not found in buiding block %d%s "
                                  "while combining tdb and rtp",
                                  (*hb)[i].hack[j].oname != nullptr ?
                                  (*hb)[i].hack[j].oname : (*hb)[i].hack[j].a[0],
                                  i+1, *resinfo[i].rtp);
                    }
                }
                else
                {
                    if ( (*hb)[i].hack[j].oname == nullptr)
                    {
                        /* we're adding: */
                        add_atom_to_restp(&(*restp)[i], l, &(*hb)[i].hack[j]);
                    }
                    else
                    {
                        /* oname != NULL */
                        if ( (*hb)[i].hack[j].nname == nullptr)
                        {   /* we're deleting */
                            /* shift the rest */
                            (*restp)[i].natom--;
                            for (k = l; k < (*restp)[i].natom; k++)
                            {
                                (*restp)[i].atom    [k] = (*restp)[i].atom    [k+1];
                                (*restp)[i].atomname[k] = (*restp)[i].atomname[k+1];
                                (*restp)[i].cgnr    [k] = (*restp)[i].cgnr    [k+1];
                            }
                            /* give back space */
                            srenew((*restp)[i].atom,     (*restp)[i].natom);
                            srenew((*restp)[i].atomname, (*restp)[i].natom);
                            srenew((*restp)[i].cgnr,     (*restp)[i].natom);
                        }
                        else /* nname != NULL */
                        {    /* we're replacing */
                            snew( (*restp)[i].atomname[l], 1);
                            (*restp)[i].atom[l]      =       *(*hb)[i].hack[j].atom;
                            *(*restp)[i].atomname[l] = gmx_strdup((*hb)[i].hack[j].nname);
                            if ( (*hb)[i].hack[j].cgnr != NOTSET)
                            {
                                (*restp)[i].cgnr   [l] = (*hb)[i].hack[j].cgnr;
                            }
                        }
                    }
                }
            }
        }
    }
}

static bool atomname_cmp_nr(const char *anm, t_hack *hack, int *nr)
{

    if (hack->nr == 1)
    {
        *nr = 0;

        return (gmx_strcasecmp(anm, hack->nname) == 0);
    }
    else
    {
        if (isdigit(anm[strlen(anm)-1]))
        {
            *nr = anm[strlen(anm)-1] - '0';
        }
        else
        {
            *nr = 0;
        }
        if (*nr <= 0 || *nr > hack->nr)
        {
            return FALSE;
        }
        else
        {
            return (strlen(anm) == strlen(hack->nname) + 1 &&
                    gmx_strncasecmp(anm, hack->nname, strlen(hack->nname)) == 0);
        }
    }
}

static bool match_atomnames_with_rtp_atom(t_atoms *pdba, rvec *x, int atind,
                                          t_restp *rptr, t_hackblock *hbr,
                                          bool bVerbose)
{
    int      resnr;
    int      j, k;
    char    *oldnm, *newnm;
    int      anmnr;
    char    *start_at, buf[STRLEN];
    int      start_nr;
    bool     bReplaceReplace, bFoundInAdd;
    bool     bDeleted;

    oldnm = *pdba->atomname[atind];
    resnr = pdba->resinfo[pdba->atom[atind].resind].nr;

    bDeleted = FALSE;
    for (j = 0; j < hbr->nhack; j++)
    {
        if (hbr->hack[j].oname != nullptr && hbr->hack[j].nname != nullptr &&
            gmx_strcasecmp(oldnm, hbr->hack[j].oname) == 0)
        {
            /* This is a replace entry. */
            /* Check if we are not replacing a replaced atom. */
            bReplaceReplace = FALSE;
            for (k = 0; k < hbr->nhack; k++)
            {
                if (k != j &&
                    hbr->hack[k].oname != nullptr && hbr->hack[k].nname != nullptr &&
                    gmx_strcasecmp(hbr->hack[k].nname, hbr->hack[j].oname) == 0)
                {
                    /* The replace in hack[j] replaces an atom that
                     * was already replaced in hack[k], we do not want
                     * second or higher level replaces at this stage.
                     */
                    bReplaceReplace = TRUE;
                }
            }
            if (bReplaceReplace)
            {
                /* Skip this replace. */
                continue;
            }

            /* This atom still has the old name, rename it */
            newnm = hbr->hack[j].nname;
            for (k = 0; k < rptr->natom; k++)
            {
                if (gmx_strcasecmp(newnm, *rptr->atomname[k]) == 0)
                {
                    break;
                }
            }
            if (k == rptr->natom)
            {
                /* The new name is not present in the rtp.
                 * We need to apply the replace also to the rtp entry.
                 */

                /* We need to find the add hack that can add this atom
                 * to find out after which atom it should be added.
                 */
                bFoundInAdd = FALSE;
                for (k = 0; k < hbr->nhack; k++)
                {
                    if (hbr->hack[k].oname == nullptr &&
                        hbr->hack[k].nname != nullptr &&
                        atomname_cmp_nr(newnm, &hbr->hack[k], &anmnr))
                    {
                        if (anmnr <= 1)
                        {
                            start_at = hbr->hack[k].a[0];
                        }
                        else
                        {
                            sprintf(buf, "%s%d", hbr->hack[k].nname, anmnr-1);
                            start_at = buf;
                        }
                        for (start_nr = 0; start_nr < rptr->natom; start_nr++)
                        {
                            if (gmx_strcasecmp(start_at, (*rptr->atomname[start_nr])) == 0)
                            {
                                break;
                            }
                        }
                        if (start_nr == rptr->natom)
                        {
                            gmx_fatal(FARGS, "Could not find atom '%s' in residue building block '%s' to add atom '%s' to",
                                      start_at, rptr->resname, newnm);
                        }
                        /* We can add the atom after atom start_nr */
                        add_atom_to_restp(rptr, start_nr, &hbr->hack[j]);

                        bFoundInAdd = TRUE;
                    }
                }

                if (!bFoundInAdd)
                {
                    gmx_fatal(FARGS, "Could not find an 'add' entry for atom named '%s' corresponding to the 'replace' entry from atom name '%s' to '%s' for tdb or hdb database of residue type '%s'",
                              newnm,
                              hbr->hack[j].oname, hbr->hack[j].nname,
                              rptr->resname);
                }
            }

            if (bVerbose)
            {
                printf("Renaming atom '%s' in residue '%s' %d to '%s'\n",
                       oldnm, rptr->resname, resnr, newnm);
            }
            /* Rename the atom in pdba */
            snew(pdba->atomname[atind], 1);
            *pdba->atomname[atind] = gmx_strdup(newnm);
        }
        else if (hbr->hack[j].oname != nullptr && hbr->hack[j].nname == nullptr &&
                 gmx_strcasecmp(oldnm, hbr->hack[j].oname) == 0)
        {
            /* This is a delete entry, check if this atom is present
             * in the rtp entry of this residue.
             */
            for (k = 0; k < rptr->natom; k++)
            {
                if (gmx_strcasecmp(oldnm, *rptr->atomname[k]) == 0)
                {
                    break;
                }
            }
            if (k == rptr->natom)
            {
                /* This atom is not present in the rtp entry,
                 * delete is now from the input pdba.
                 */
                if (bVerbose)
                {
                    printf("Deleting atom '%s' in residue '%s' %d\n",
                           oldnm, rptr->resname, resnr);
                }
                /* We should free the atom name,
                 * but it might be used multiple times in the symtab.
                 * sfree(pdba->atomname[atind]);
                 */
                for (k = atind+1; k < pdba->nr; k++)
                {
                    pdba->atom[k-1]     = pdba->atom[k];
                    pdba->atomname[k-1] = pdba->atomname[k];
                    copy_rvec(x[k], x[k-1]);
                }
                pdba->nr--;
                bDeleted = TRUE;
            }
        }
    }

    return bDeleted;
}

void match_atomnames_with_rtp(t_restp restp[], t_hackblock hb[],
                              t_atoms *pdba, rvec *x,
                              bool bVerbose)
{
    int          i, j;
    char        *oldnm;
    t_restp     *rptr;

    for (i = 0; i < pdba->nr; i++)
    {
        oldnm = *pdba->atomname[i];
        rptr  = &restp[pdba->atom[i].resind];
        for (j = 0; (j < rptr->natom); j++)
        {
            if (gmx_strcasecmp(oldnm, *(rptr->atomname[j])) == 0)
            {
                break;
            }
        }
        if (j == rptr->natom)
        {
            /* Not found yet, check if we have to rename this atom */
            if (match_atomnames_with_rtp_atom(pdba, x, i,
                                              rptr, &(hb[pdba->atom[i].resind]),
                                              bVerbose))
            {
                /* We deleted this atom, decrease the atom counter by 1. */
                i--;
            }
        }
    }
}

#define NUM_CMAP_ATOMS 5
static void gen_cmap(t_params *psb, t_restp *restp, t_atoms *atoms)
{
    int         residx, i, j, k;
    const char *ptr;
    const char *pname;
    t_resinfo  *resinfo = atoms->resinfo;
    int         nres    = atoms->nres;
    bool        bAddCMAP;
    int         cmap_atomid[NUM_CMAP_ATOMS];
    int         cmap_chainnum = -1, this_residue_index;

    if (debug)
    {
        ptr = "cmap";
    }
    else
    {
        ptr = "check";
    }

    fprintf(stderr, "Making cmap torsions...\n");
    i = 0;
    /* Most cmap entries use the N atom from the next residue, so the last
     * residue should not have its CMAP entry in that case, but for things like
     * dipeptides we sometimes define a complete CMAP entry inside a residue,
     * and in this case we need to process everything through the last residue.
     */
    for (residx = 0; residx < nres; residx++)
    {
        /* Add CMAP terms from the list of CMAP interactions */
        for (j = 0; j < restp[residx].rb[ebtsCMAP].nb; j++)
        {
            bAddCMAP = TRUE;
            /* Loop over atoms in a candidate CMAP interaction and
             * check that they exist, are from the same chain and are
             * from residues labelled as protein. */
            for (k = 0; k < NUM_CMAP_ATOMS && bAddCMAP; k++)
            {
                /* Assign the pointer to the name of the next reference atom.
                 * This can use -/+ labels to refer to previous/next residue.
                 */
                pname = restp[residx].rb[ebtsCMAP].b[j].a[k];
                /* Skip this CMAP entry if it refers to residues before the
                 * first or after the last residue.
                 */
                if (((strchr(pname, '-') != nullptr) && (residx == 0)) ||
                    ((strchr(pname, '+') != nullptr) && (residx == nres-1)))
                {
                    bAddCMAP = FALSE;
                    break;
                }

                cmap_atomid[k] = search_atom(pname,
                                             i, atoms, ptr, TRUE);
                bAddCMAP = bAddCMAP && (cmap_atomid[k] != -1);
                if (!bAddCMAP)
                {
                    /* This break is necessary, because cmap_atomid[k]
                     * == -1 cannot be safely used as an index
                     * into the atom array. */
                    break;
                }
                this_residue_index = atoms->atom[cmap_atomid[k]].resind;
                if (0 == k)
                {
                    cmap_chainnum = resinfo[this_residue_index].chainnum;
                }
                else
                {
                    /* Does the residue for this atom have the same
                     * chain number as the residues for previous
                     * atoms? */
                    bAddCMAP = bAddCMAP &&
                        cmap_chainnum == resinfo[this_residue_index].chainnum;
                }
                /* Here we used to check that the residuetype was protein and
                 * disable bAddCMAP if that was not the case. However, some
                 * special residues (say, alanine dipeptides) might not adhere
                 * to standard naming, and if we start calling them normal
                 * protein residues the user will be bugged to select termini.
                 *
                 * Instead, I believe that the right course of action is to
                 * keep the CMAP interaction if it is present in the RTP file
                 * and we correctly identified all atoms (which is the case
                 * if we got here).
                 */
            }

            if (bAddCMAP)
            {
                add_cmap_param(psb, cmap_atomid[0], cmap_atomid[1], cmap_atomid[2], cmap_atomid[3], cmap_atomid[4], restp[residx].rb[ebtsCMAP].b[j].s);
            }
        }

        if (residx < nres-1)
        {
            while (atoms->atom[i].resind < residx+1)
            {
                i++;
            }
        }
    }
    /* Start the next residue */
}

static void
scrub_charge_groups(int *cgnr, int natoms)
{
    int i;

    for (i = 0; i < natoms; i++)
    {
        cgnr[i] = i+1;
    }
}


void pdb2top(FILE *top_file, const char *posre_fn, const char *molname,
             t_atoms *atoms, rvec **x, gpp_atomtype_t atype, t_symtab *tab,
             int nrtp, t_restp rtp[],
             t_restp *restp, t_hackblock *hb,
             bool bAllowMissing,
             bool bVsites, bool bVsiteAromatics,
             const char *ffdir,
             real mHmult,
             int nssbonds, t_ssbond *ssbonds,
             real long_bond_dist, real short_bond_dist,
             bool bDeuterate, bool bChargeGroups, bool bCmap,
             bool bRenumRes, bool bRTPresname)
{
    /*
       t_hackblock *hb;
       t_restp  *restp;
     */
    t_params          plist[F_NRE];
    t_excls          *excls;
    t_nextnb          nnb;
    int              *cgnr;
    int              *vsite_type;
    int               i, nmissat;
    int               bts[ebtsNR];
    gmx_residuetype_t*rt;

    init_plist(plist);
    gmx_residuetype_init(&rt);

    /* Make bonds */
    at2bonds(&(plist[F_BONDS]), hb,
             atoms, *x,
             long_bond_dist, short_bond_dist);

    /* specbonds: disulphide bonds & heme-his */
    do_ssbonds(&(plist[F_BONDS]),
               atoms, nssbonds, ssbonds,
               bAllowMissing);

    nmissat = name2type(atoms, &cgnr, restp, rt);
    if (nmissat)
    {
        if (bAllowMissing)
        {
            fprintf(stderr, "There were %d missing atoms in molecule %s\n",
                    nmissat, molname);
        }
        else
        {
            gmx_fatal(FARGS, "There were %d missing atoms in molecule %s, if you want to use this incomplete topology anyhow, use the option -missing",
                      nmissat, molname);
        }
    }

    /* Cleanup bonds (sort and rm doubles) */
    clean_bonds(&(plist[F_BONDS]));

    snew(vsite_type, atoms->nr);
    for (i = 0; i < atoms->nr; i++)
    {
        vsite_type[i] = NOTSET;
    }
    if (bVsites)
    {
        if (bVsiteAromatics)
        {
            fprintf(stdout, "The conversion of aromatic rings into virtual sites is deprecated "
                    "and may be removed in a future version of GROMACS");
        }
        /* determine which atoms will be vsites and add dummy masses
           also renumber atom numbers in plist[0..F_NRE]! */
        do_vsites(nrtp, rtp, atype, atoms, tab, x, plist,
                  &vsite_type, &cgnr, mHmult, bVsiteAromatics, ffdir);
    }

    /* Make Angles and Dihedrals */
    fprintf(stderr, "Generating angles, dihedrals and pairs...\n");
    snew(excls, atoms->nr);
    init_nnb(&nnb, atoms->nr, 4);
    gen_nnb(&nnb, plist);
    print_nnb(&nnb, "NNB");
    gen_pad(&nnb, atoms, restp, plist, excls, hb, bAllowMissing);
    done_nnb(&nnb);

    /* Make CMAP */
    if (bCmap)
    {
        gen_cmap(&(plist[F_CMAP]), restp, atoms);
        if (plist[F_CMAP].nr > 0)
        {
            fprintf(stderr, "There are %4d cmap torsion pairs\n",
                    plist[F_CMAP].nr);
        }
    }

    /* set mass of all remaining hydrogen atoms */
    if (mHmult != 1.0)
    {
        do_h_mass(&(plist[F_BONDS]), vsite_type, atoms, mHmult, bDeuterate);
    }
    sfree(vsite_type);

    /* Cleanup bonds (sort and rm doubles) */
    /* clean_bonds(&(plist[F_BONDS]));*/

    fprintf(stderr,
            "There are %4d dihedrals, %4d impropers, %4d angles\n"
            "          %4d pairs,     %4d bonds and  %4d virtual sites\n",
            plist[F_PDIHS].nr, plist[F_IDIHS].nr, plist[F_ANGLES].nr,
            plist[F_LJ14].nr, plist[F_BONDS].nr,
            plist[F_VSITE2].nr +
            plist[F_VSITE3].nr +
            plist[F_VSITE3FD].nr +
            plist[F_VSITE3FAD].nr +
            plist[F_VSITE3OUT].nr +
            plist[F_VSITE4FD].nr +
            plist[F_VSITE4FDN].nr );

    print_sums(atoms, FALSE);

    if (!bChargeGroups)
    {
        scrub_charge_groups(cgnr, atoms->nr);
    }

    if (bRenumRes)
    {
        for (i = 0; i < atoms->nres; i++)
        {
            atoms->resinfo[i].nr = i + 1;
            atoms->resinfo[i].ic = ' ';
        }
    }

    if (top_file)
    {
        fprintf(stderr, "Writing topology\n");
        /* We can copy the bonded types from the first restp,
         * since the types have to be identical for all residues in one molecule.
         */
        for (i = 0; i < ebtsNR; i++)
        {
            bts[i] = restp[0].rb[i].type;
        }
        write_top(top_file, posre_fn, molname,
                  atoms, bRTPresname,
                  bts, plist, excls, atype, cgnr, restp[0].nrexcl);
    }

    /* cleaning up */
    free_t_hackblock(atoms->nres, &hb);
    free_t_restp(atoms->nres, &restp);
    gmx_residuetype_destroy(rt);

    /* we should clean up hb and restp here, but that is a *L*O*T* of work! */
    sfree(cgnr);
    for (i = 0; i < F_NRE; i++)
    {
        sfree(plist[i].param);
    }
    for (i = 0; i < atoms->nr; i++)
    {
        sfree(excls[i].e);
    }
    sfree(excls);
}

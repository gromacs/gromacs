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
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/h_db.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/pgutil.h"
#include "gromacs/gmxpreprocess/specbond.h"
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

#include "hackblock.h"
#include "resall.h"

/* this must correspond to enum in pdb2top.h */
const char* hh[ehisNR] = { "HISD", "HISE", "HISH", "HIS1" };

static int missing_atoms(const PreprocessResidue* rp, int resind, t_atoms* at, int i0, int i)
{
    int nmiss = 0;
    for (int j = 0; j < rp->natom(); j++)
    {
        const char* name   = *(rp->atomname[j]);
        bool        bFound = false;
        for (int k = i0; k < i; k++)
        {
            bFound = (bFound || (gmx_strcasecmp(*(at->atomname[k]), name) == 0));
        }
        if (!bFound)
        {
            nmiss++;
            fprintf(stderr,
                    "\nWARNING: "
                    "atom %s is missing in residue %s %d in the pdb file\n",
                    name, *(at->resinfo[resind].name), at->resinfo[resind].nr);
            if (name[0] == 'H' || name[0] == 'h')
            {
                fprintf(stderr,
                        "         You might need to add atom %s to the hydrogen database of "
                        "building block %s\n"
                        "         in the file %s.hdb (see the manual)\n",
                        name, *(at->resinfo[resind].rtp), rp->filebase.c_str());
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

    return (fabs(x - ix) < tol);
}

static void choose_ff_impl(const char* ffsel, char* forcefield, int ff_maxlen, char* ffdir, int ffdir_maxlen)
{
    std::vector<gmx::DataFileInfo> ffdirs = fflib_enumerate_forcefields();
    const int                      nff    = ssize(ffdirs);

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
        ffs.push_back(gmx::stripSuffixIfPresent(ffdirs[i].name, fflib_forcefield_dir_ext()));
    }

    int sel;
    if (ffsel != nullptr)
    {
        sel        = -1;
        int cwdsel = -1;
        int nfound = 0;
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
                        "you would prefer a different one.\n",
                        ffsel, nfound);
            }
            else
            {
                std::string message = gmx::formatString(
                        "Force field '%s' occurs in %d places, but not in "
                        "the current directory.\n"
                        "Run without the -ff switch and select the force "
                        "field interactively.",
                        ffsel, nfound);
                GMX_THROW(gmx::InconsistentInputError(message));
            }
        }
        else if (nfound == 0)
        {
            std::string message = gmx::formatString(
                    "Could not find force field '%s' in current directory, "
                    "install tree or GMXLIB path.",
                    ffsel);
            GMX_THROW(gmx::InconsistentInputError(message));
        }
    }
    else if (nff > 1)
    {
        std::vector<std::string> desc;
        desc.reserve(ffdirs.size());
        for (int i = 0; i < nff; ++i)
        {
            std::string docFileName(gmx::Path::join(ffdirs[i].dir, ffdirs[i].name, fflib_forcefield_doc()));
            // TODO: Just try to open the file with a method that does not
            // throw/bail out with a fatal error instead of multiple checks.
            if (gmx::File::exists(docFileName, gmx::File::returnFalseOnError))
            {
                // TODO: Use a C++ API without such an intermediate/fixed-length buffer.
                char buf[STRLEN];
                /* We don't use fflib_open, because we don't want printf's */
                FILE* fp = gmx_ffopen(docFileName, "r");
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
                if (ffdirs[i].dir == ffdirs[j].dir
                    && ((desc[i][0] == '[' && desc[j][0] != '[')
                        || ((desc[i][0] == '[' || desc[j][0] != '[')
                            && gmx_strcasecmp(desc[i].c_str(), desc[j].c_str()) > 0)))
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
            if (i == 0 || ffdirs[i - 1].dir != ffdirs[i].dir)
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
            printf("%2d: %s\n", i + 1, desc[i].c_str());
        }

        sel = -1;
        // TODO: Add a C++ API for this.
        char  buf[STRLEN];
        char* pret;
        do
        {
            pret = fgets(buf, STRLEN, stdin);

            if (pret != nullptr)
            {
                sel = strtol(buf, nullptr, 10);
                sel--;
            }
        } while (pret == nullptr || (sel < 0) || (sel >= nff));

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
        std::string message = gmx::formatString("Length of force field name (%d) >= maxlen (%d)",
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
        std::string message = gmx::formatString("Length of force field dir (%d) >= maxlen (%d)",
                                                static_cast<int>(ffpath.length()), ffdir_maxlen);
        GMX_THROW(gmx::InvalidInputError(message));
    }
    strcpy(ffdir, ffpath.c_str());
}

void choose_ff(const char* ffsel, char* forcefield, int ff_maxlen, char* ffdir, int ffdir_maxlen)
{
    try
    {
        choose_ff_impl(ffsel, forcefield, ff_maxlen, ffdir, ffdir_maxlen);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
}

void choose_watermodel(const char* wmsel, const char* ffdir, char** watermodel)
{
    const char* fn_watermodels = "watermodels.dat";
    FILE*       fp;
    char        buf[STRLEN];
    int         nwm, sel, i;
    char**      model;
    char*       pret;

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
        fprintf(stderr, "No file '%s' found, will not include a water model\n", fn_watermodels);
        *watermodel = nullptr;

        return;
    }

    fp = fflib_open(filename);
    printf("\nSelect the Water Model:\n");
    nwm   = 0;
    model = nullptr;
    while (get_a_line(fp, buf, STRLEN))
    {
        srenew(model, nwm + 1);
        snew(model[nwm], STRLEN);
        sscanf(buf, "%s%n", model[nwm], &i);
        if (i > 0)
        {
            ltrim(buf + i);
            fprintf(stderr, "%2d: %s\n", nwm + 1, buf + i);
            nwm++;
        }
        else
        {
            sfree(model[nwm]);
        }
    }
    gmx_ffclose(fp);
    fprintf(stderr, "%2d: %s\n", nwm + 1, "None");

    sel = -1;
    do
    {
        pret = fgets(buf, STRLEN, stdin);

        if (pret != nullptr)
        {
            sel = strtol(buf, nullptr, 10);
            sel--;
        }
    } while (pret == nullptr || sel < 0 || sel > nwm);

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

static int name2type(t_atoms* at, int** cgnr, gmx::ArrayRef<const PreprocessResidue> usedPpResidues, ResidueType* rt)
{
    int    i, j, prevresind, i0, prevcg, cg, curcg;
    char*  name;
    bool   bNterm;
    double qt;
    int    nmissat;

    nmissat = 0;

    int resind = -1;
    bNterm     = false;
    i0         = 0;
    snew(*cgnr, at->nr);
    qt    = 0;
    curcg = 0;
    cg    = -1;

    for (i = 0; (i < at->nr); i++)
    {
        prevresind = resind;
        if (at->atom[i].resind != resind)
        {
            resind     = at->atom[i].resind;
            bool bProt = rt->namedResidueHasType(*(at->resinfo[resind].name), "Protein");
            bNterm     = bProt && (resind == 0);
            if (resind > 0)
            {
                nmissat += missing_atoms(&usedPpResidues[prevresind], prevresind, at, i0, i);
            }
            i0 = i;
        }
        if (at->atom[i].m == 0)
        {
            qt               = 0;
            prevcg           = cg;
            name             = *(at->atomname[i]);
            j                = search_jtype(usedPpResidues[resind], name, bNterm);
            at->atom[i].type = usedPpResidues[resind].atom[j].type;
            at->atom[i].q    = usedPpResidues[resind].atom[j].q;
            at->atom[i].m    = usedPpResidues[resind].atom[j].m;
            cg               = usedPpResidues[resind].cgnr[j];
            /* A charge group number -1 signals a separate charge group
             * for this atom.
             */
            if ((cg == -1) || (cg != prevcg) || (resind != prevresind))
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
    nmissat += missing_atoms(&usedPpResidues[resind], resind, at, i0, i);

    return nmissat;
}

static void print_top_heavy_H(FILE* out, real mHmult)
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
        fprintf(stderr,
                "WARNING: unsupported proton mass multiplier (%g) "
                "in pdb2top\n",
                mHmult);
    }
}

void print_top_comment(FILE* out, const char* filename, const char* ffdir, bool bITP)
{
    char  ffdir_parent[STRLEN];
    char* p;

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
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

    if (strchr(ffdir, '/') == nullptr)
    {
        fprintf(out, ";\tForce field was read from the standard GROMACS share directory.\n;\n\n");
    }
    else if (ffdir[0] == '.')
    {
        fprintf(out,
                ";\tForce field was read from current directory or a relative path - path "
                "added.\n;\n\n");
    }
    else
    {
        strncpy(ffdir_parent, ffdir, STRLEN - 1);
        ffdir_parent[STRLEN - 1] = '\0'; /*make sure it is 0-terminated even for long string*/
        p                        = strrchr(ffdir_parent, '/');

        *p = '\0';

        fprintf(out,
                ";\tForce field data was read from:\n"
                ";\t%s\n"
                ";\n"
                ";\tNote:\n"
                ";\tThis might be a non-standard force field location. When you use this topology, "
                "the\n"
                ";\tforce field must either be present in the current directory, or the location\n"
                ";\tspecified in the GMXLIB path variable or with the 'include' mdp file "
                "option.\n;\n\n",
                ffdir_parent);
    }
}

void print_top_header(FILE* out, const char* filename, bool bITP, const char* ffdir, real mHmult)
{
    const char* p;

    print_top_comment(out, filename, ffdir, bITP);

    print_top_heavy_H(out, mHmult);
    fprintf(out, "; Include forcefield parameters\n");

    p = strrchr(ffdir, '/');
    p = (ffdir[0] == '.' || p == nullptr) ? ffdir : p + 1;

    fprintf(out, "#include \"%s/%s\"\n\n", p, fflib_forcefield_itp());
}

static void print_top_posre(FILE* out, const char* pr)
{
    fprintf(out, "; Include Position restraint file\n");
    fprintf(out, "#ifdef POSRES\n");
    fprintf(out, "#include \"%s\"\n", pr);
    fprintf(out, "#endif\n\n");
}

static void print_top_water(FILE* out, const char* ffdir, const char* water)
{
    const char* p;
    char        buf[STRLEN];

    fprintf(out, "; Include water topology\n");

    p = strrchr(ffdir, '/');
    p = (ffdir[0] == '.' || p == nullptr) ? ffdir : p + 1;
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

static void print_top_system(FILE* out, const char* title)
{
    fprintf(out, "[ %s ]\n", dir2str(Directive::d_system));
    fprintf(out, "; Name\n");
    fprintf(out, "%s\n\n", title[0] ? title : "Protein");
}

void print_top_mols(FILE*                            out,
                    const char*                      title,
                    const char*                      ffdir,
                    const char*                      water,
                    gmx::ArrayRef<const std::string> incls,
                    gmx::ArrayRef<const t_mols>      mols)
{

    if (!incls.empty())
    {
        fprintf(out, "; Include chain topologies\n");
        for (const auto& incl : incls)
        {
            fprintf(out, "#include \"%s\"\n", gmx::Path::getFilename(incl).c_str());
        }
        fprintf(out, "\n");
    }

    if (water)
    {
        print_top_water(out, ffdir, water);
    }
    print_top_system(out, title);

    if (!mols.empty())
    {
        fprintf(out, "[ %s ]\n", dir2str(Directive::d_molecules));
        fprintf(out, "; %-15s %5s\n", "Compound", "#mols");
        for (const auto& mol : mols)
        {
            fprintf(out, "%-15s %5d\n", mol.name.c_str(), mol.nr);
        }
    }
}

void write_top(FILE*                                   out,
               const char*                             pr,
               const char*                             molname,
               t_atoms*                                at,
               bool                                    bRTPresname,
               int                                     bts[],
               gmx::ArrayRef<const InteractionsOfType> plist,
               t_excls                                 excls[],
               PreprocessingAtomTypes*                 atype,
               int*                                    cgnr,
               int                                     nrexcl)
/* NOTE: nrexcl is not the size of *excl! */
{
    if (at && atype && cgnr)
    {
        fprintf(out, "[ %s ]\n", dir2str(Directive::d_moleculetype));
        fprintf(out, "; %-15s %5s\n", "Name", "nrexcl");
        fprintf(out, "%-15s %5d\n\n", molname ? molname : "Protein", nrexcl);

        print_atoms(out, atype, at, cgnr, bRTPresname);
        print_bondeds(out, at->nr, Directive::d_bonds, F_BONDS, bts[ebtsBONDS], plist);
        print_bondeds(out, at->nr, Directive::d_constraints, F_CONSTR, 0, plist);
        print_bondeds(out, at->nr, Directive::d_constraints, F_CONSTRNC, 0, plist);
        print_bondeds(out, at->nr, Directive::d_pairs, F_LJ14, 0, plist);
        print_excl(out, at->nr, excls);
        print_bondeds(out, at->nr, Directive::d_angles, F_ANGLES, bts[ebtsANGLES], plist);
        print_bondeds(out, at->nr, Directive::d_dihedrals, F_PDIHS, bts[ebtsPDIHS], plist);
        print_bondeds(out, at->nr, Directive::d_dihedrals, F_IDIHS, bts[ebtsIDIHS], plist);
        print_bondeds(out, at->nr, Directive::d_cmap, F_CMAP, bts[ebtsCMAP], plist);
        print_bondeds(out, at->nr, Directive::d_polarization, F_POLARIZATION, 0, plist);
        print_bondeds(out, at->nr, Directive::d_thole_polarization, F_THOLE_POL, 0, plist);
        print_bondeds(out, at->nr, Directive::d_vsites2, F_VSITE2, 0, plist);
        print_bondeds(out, at->nr, Directive::d_vsites3, F_VSITE3, 0, plist);
        print_bondeds(out, at->nr, Directive::d_vsites3, F_VSITE3FD, 0, plist);
        print_bondeds(out, at->nr, Directive::d_vsites3, F_VSITE3FAD, 0, plist);
        print_bondeds(out, at->nr, Directive::d_vsites3, F_VSITE3OUT, 0, plist);
        print_bondeds(out, at->nr, Directive::d_vsites4, F_VSITE4FD, 0, plist);
        print_bondeds(out, at->nr, Directive::d_vsites4, F_VSITE4FDN, 0, plist);

        if (pr)
        {
            print_top_posre(out, pr);
        }
    }
}


static void do_ssbonds(InteractionsOfType*                ps,
                       t_atoms*                           atoms,
                       gmx::ArrayRef<const DisulfideBond> ssbonds,
                       bool                               bAllowMissing)
{
    for (const auto& bond : ssbonds)
    {
        int ri = bond.firstResidue;
        int rj = bond.secondResidue;
        int ai = search_res_atom(bond.firstAtom.c_str(), ri, atoms, "special bond", bAllowMissing);
        int aj = search_res_atom(bond.secondAtom.c_str(), rj, atoms, "special bond", bAllowMissing);
        if ((ai == -1) || (aj == -1))
        {
            gmx_fatal(FARGS, "Trying to make impossible special bond (%s-%s)!",
                      bond.firstAtom.c_str(), bond.secondAtom.c_str());
        }
        add_param(ps, ai, aj, {}, nullptr);
    }
}

static void at2bonds(InteractionsOfType*                  psb,
                     gmx::ArrayRef<MoleculePatchDatabase> globalPatches,
                     t_atoms*                             atoms,
                     gmx::ArrayRef<const gmx::RVec>       x,
                     real                                 long_bond_dist,
                     real                                 short_bond_dist)
{
    real        long_bond_dist2, short_bond_dist2;
    const char* ptr;

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
    int i = 0;
    for (int resind = 0; (resind < atoms->nres) && (i < atoms->nr); resind++)
    {
        /* add bonds from list of bonded interactions */
        for (const auto& patch : globalPatches[resind].rb[ebtsBONDS].b)
        {
            /* Unfortunately we can not issue errors or warnings
             * for missing atoms in bonds, as the hydrogens and terminal atoms
             * have not been added yet.
             */
            int ai = search_atom(patch.ai().c_str(), i, atoms, ptr, TRUE);
            int aj = search_atom(patch.aj().c_str(), i, atoms, ptr, TRUE);
            if (ai != -1 && aj != -1)
            {
                real dist2 = distance2(x[ai], x[aj]);
                if (dist2 > long_bond_dist2)

                {
                    fprintf(stderr, "Warning: Long Bond (%d-%d = %g nm)\n", ai + 1, aj + 1,
                            std::sqrt(dist2));
                }
                else if (dist2 < short_bond_dist2)
                {
                    fprintf(stderr, "Warning: Short Bond (%d-%d = %g nm)\n", ai + 1, aj + 1,
                            std::sqrt(dist2));
                }
                add_param(psb, ai, aj, {}, patch.s.c_str());
            }
        }
        /* add bonds from list of hacks (each added atom gets a bond) */
        while ((i < atoms->nr) && (atoms->atom[i].resind == resind))
        {
            for (const auto& patch : globalPatches[resind].hack)
            {
                if ((patch.tp > 0 || patch.type() == MoleculePatchType::Add)
                    && patch.a[0] == *(atoms->atomname[i]))
                {
                    switch (patch.tp)
                    {
                        case 9:                                        /* COOH terminus */
                            add_param(psb, i, i + 1, {}, nullptr);     /* C-O  */
                            add_param(psb, i, i + 2, {}, nullptr);     /* C-OA */
                            add_param(psb, i + 2, i + 3, {}, nullptr); /* OA-H */
                            break;
                        default:
                            for (int k = 0; (k < patch.nr); k++)
                            {
                                add_param(psb, i, i + k + 1, {}, nullptr);
                            }
                    }
                }
            }
            i++;
        }
        /* we're now at the start of the next residue */
    }
}

static bool pcompar(const InteractionOfType& a, const InteractionOfType& b)
{
    int d;

    if ((d = a.ai() - b.ai()) != 0)
    {
        return d < 0;
    }
    else if ((d = a.aj() - b.aj()) != 0)
    {
        return d < 0;
    }
    else
    {
        return a.interactionTypeName().length() > b.interactionTypeName().length();
    }
}

static void clean_bonds(InteractionsOfType* ps)
{
    if (ps->size() > 0)
    {
        /* Sort bonds */
        for (auto& bond : ps->interactionTypes)
        {
            bond.sortAtomIds();
        }
        std::sort(ps->interactionTypes.begin(), ps->interactionTypes.end(), pcompar);

        /* remove doubles, keep the first one always. */
        int oldNumber = ps->size();
        for (auto parm = ps->interactionTypes.begin() + 1; parm != ps->interactionTypes.end();)
        {
            auto prev = parm - 1;
            if (parm->ai() == prev->ai() && parm->aj() == prev->aj())
            {
                parm = ps->interactionTypes.erase(parm);
            }
            else
            {
                ++parm;
            }
        }
        fprintf(stderr, "Number of bonds was %d, now %zu\n", oldNumber, ps->size());
    }
    else
    {
        fprintf(stderr, "No bonds\n");
    }
}

void print_sums(const t_atoms* atoms, bool bSystem)
{
    double      m, qtot;
    int         i;
    const char* where;

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
        m += atoms->atom[i].m;
        qtot += atoms->atom[i].q;
    }

    fprintf(stderr, "Total mass%s %.3f a.m.u.\n", where, m);
    fprintf(stderr, "Total charge%s %.3f e\n", where, qtot);
}

static void check_restp_type(const char* name, int t1, int t2)
{
    if (t1 != t2)
    {
        gmx_fatal(FARGS, "Residues in one molecule have a different '%s' type: %d and %d", name, t1, t2);
    }
}

static void check_restp_types(const PreprocessResidue& r0, const PreprocessResidue& r1)
{
    check_restp_type("all dihedrals", static_cast<int>(r0.bKeepAllGeneratedDihedrals),
                     static_cast<int>(r1.bKeepAllGeneratedDihedrals));
    check_restp_type("nrexcl", r0.nrexcl, r1.nrexcl);
    check_restp_type("HH14", static_cast<int>(r0.bGenerateHH14Interactions),
                     static_cast<int>(r1.bGenerateHH14Interactions));
    check_restp_type("remove dihedrals", static_cast<int>(r0.bRemoveDihedralIfWithImproper),
                     static_cast<int>(r1.bRemoveDihedralIfWithImproper));

    for (int i = 0; i < ebtsNR; i++)
    {
        check_restp_type(btsNames[i], r0.rb[i].type, r1.rb[i].type);
    }
}

static void add_atom_to_restp(PreprocessResidue*   usedPpResidues,
                              t_symtab*            symtab,
                              int                  at_start,
                              const MoleculePatch* patch)
{
    /* now add them */
    for (int k = 0; k < patch->nr; k++)
    {
        /* set counter in atomname */
        std::string buf = patch->nname;
        if (patch->nr > 1)
        {
            buf.append(gmx::formatString("%d", k + 1));
        }
        usedPpResidues->atomname.insert(usedPpResidues->atomname.begin() + at_start + 1 + k,
                                        put_symtab(symtab, buf.c_str()));
        usedPpResidues->atom.insert(usedPpResidues->atom.begin() + at_start + 1 + k, patch->atom.back());
        if (patch->cgnr != NOTSET)
        {
            usedPpResidues->cgnr.insert(usedPpResidues->cgnr.begin() + at_start + 1 + k, patch->cgnr);
        }
        else
        {
            usedPpResidues->cgnr.insert(usedPpResidues->cgnr.begin() + at_start + 1 + k,
                                        usedPpResidues->cgnr[at_start]);
        }
    }
}

void get_hackblocks_rtp(std::vector<MoleculePatchDatabase>*    globalPatches,
                        std::vector<PreprocessResidue>*        usedPpResidues,
                        gmx::ArrayRef<const PreprocessResidue> rtpFFDB,
                        int                                    nres,
                        t_resinfo*                             resinfo,
                        int                                    nterpairs,
                        t_symtab*                              symtab,
                        gmx::ArrayRef<MoleculePatchDatabase*>  ntdb,
                        gmx::ArrayRef<MoleculePatchDatabase*>  ctdb,
                        gmx::ArrayRef<const int>               rn,
                        gmx::ArrayRef<const int>               rc,
                        bool                                   bAllowMissing)
{
    char* key;
    bool  bRM;

    globalPatches->resize(nres);
    usedPpResidues->clear();
    /* first the termini */
    for (int i = 0; i < nterpairs; i++)
    {
        if (rn[i] >= 0 && ntdb[i] != nullptr)
        {
            copyModificationBlocks(*ntdb[i], &globalPatches->at(rn[i]));
        }
        if (rc[i] >= 0 && ctdb[i] != nullptr)
        {
            mergeAtomAndBondModifications(*ctdb[i], &globalPatches->at(rc[i]));
        }
    }

    /* then the whole rtp */
    for (int i = 0; i < nres; i++)
    {
        /* Here we allow a mismatch of one character when looking for the rtp entry.
         * For such a mismatch there should be only one mismatching name.
         * This is mainly useful for small molecules such as ions.
         * Note that this will usually not work for protein, DNA and RNA,
         * since there the residue names should be listed in residuetypes.dat
         * and an error will have been generated earlier in the process.
         */
        key = *resinfo[i].rtp;

        resinfo[i].rtp = put_symtab(symtab, searchResidueDatabase(key, rtpFFDB).c_str());
        auto res       = getDatabaseEntry(*resinfo[i].rtp, rtpFFDB);
        usedPpResidues->push_back(PreprocessResidue());
        PreprocessResidue* newentry = &usedPpResidues->back();
        copyPreprocessResidues(*res, newentry, symtab);

        /* Check that we do not have different bonded types in one molecule */
        check_restp_types(usedPpResidues->front(), *newentry);

        int tern = -1;
        for (int j = 0; j < nterpairs && tern == -1; j++)
        {
            if (i == rn[j])
            {
                tern = j;
            }
        }
        int terc = -1;
        for (int j = 0; j < nterpairs && terc == -1; j++)
        {
            if (i == rc[j])
            {
                terc = j;
            }
        }
        bRM = mergeBondedInteractionList(res->rb, globalPatches->at(i).rb, tern >= 0, terc >= 0);

        if (bRM && ((tern >= 0 && ntdb[tern] == nullptr) || (terc >= 0 && ctdb[terc] == nullptr)))
        {
            const char* errString =
                    "There is a dangling bond at at least one of the terminal ends and the force "
                    "field does not provide terminal entries or files. Fix your terminal "
                    "residues so that they match the residue database (.rtp) entries, or provide "
                    "terminal database entries (.tdb).";
            if (bAllowMissing)
            {
                fprintf(stderr, "%s\n", errString);
            }
            else
            {
                gmx_fatal(FARGS, "%s", errString);
            }
        }
        else if (bRM && ((tern >= 0 && ntdb[tern]->nhack() == 0) || (terc >= 0 && ctdb[terc]->nhack() == 0)))
        {
            const char* errString =
                    "There is a dangling bond at at least one of the terminal ends. Fix your "
                    "coordinate file, add a new terminal database entry (.tdb), or select the "
                    "proper existing terminal entry.";
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

    /* Apply patchs to t_restp entries
       i.e. add's and deletes from termini database will be
       added to/removed from residue topology
       we'll do this on one big dirty loop, so it won't make easy reading! */
    for (auto modifiedResidue = globalPatches->begin(); modifiedResidue != globalPatches->end();
         modifiedResidue++)
    {
        const int          pos    = std::distance(globalPatches->begin(), modifiedResidue);
        PreprocessResidue* posres = &usedPpResidues->at(pos);
        for (auto patch = modifiedResidue->hack.begin(); patch != modifiedResidue->hack.end(); patch++)
        {
            if (patch->nr != 0)
            {
                /* find atom in restp */
                auto found = std::find_if(posres->atomname.begin(), posres->atomname.end(),
                                          [&patch](char** name) {
                                              return (patch->oname.empty() && patch->a[0] == *name)
                                                     || (patch->oname == *name);
                                          });

                if (found == posres->atomname.end())
                {
                    /* If we are doing an atom rename only, we don't need
                     * to generate a fatal error if the old name is not found
                     * in the rtp.
                     */
                    /* Deleting can happen also only on the input atoms,
                     * not necessarily always on the rtp entry.
                     */
                    if (patch->type() == MoleculePatchType::Add)
                    {
                        gmx_fatal(FARGS,
                                  "atom %s not found in buiding block %d%s "
                                  "while combining tdb and rtp",
                                  patch->oname.empty() ? patch->a[0].c_str() : patch->oname.c_str(),
                                  pos + 1, *resinfo[pos].rtp);
                    }
                }
                else
                {
                    int l = std::distance(posres->atomname.begin(), found);
                    switch (patch->type())
                    {
                        case MoleculePatchType::Add:
                        {
                            /* we're adding: */
                            add_atom_to_restp(posres, symtab, l, &(*patch));
                            break;
                        }
                        case MoleculePatchType::Delete:
                        { /* we're deleting */
                            posres->atom.erase(posres->atom.begin() + l);
                            posres->atomname.erase(posres->atomname.begin() + l);
                            posres->cgnr.erase(posres->cgnr.begin() + l);
                            break;
                        }
                        case MoleculePatchType::Replace:
                        {
                            /* we're replacing */
                            posres->atom[l]     = patch->atom.back();
                            posres->atomname[l] = put_symtab(symtab, patch->nname.c_str());
                            if (patch->cgnr != NOTSET)
                            {
                                posres->cgnr[l] = patch->cgnr;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
}

static bool atomname_cmp_nr(const char* anm, const MoleculePatch* patch, int* nr)
{

    if (patch->nr == 1)
    {
        *nr = 0;

        return (gmx::equalCaseInsensitive(anm, patch->nname));
    }
    else
    {
        if (isdigit(anm[strlen(anm) - 1]))
        {
            *nr = anm[strlen(anm) - 1] - '0';
        }
        else
        {
            *nr = 0;
        }
        if (*nr <= 0 || *nr > patch->nr)
        {
            return FALSE;
        }
        else
        {
            return (strlen(anm) == patch->nname.length() + 1
                    && gmx_strncasecmp(anm, patch->nname.c_str(), patch->nname.length()) == 0);
        }
    }
}

static bool match_atomnames_with_rtp_atom(t_atoms*                     pdba,
                                          gmx::ArrayRef<gmx::RVec>     x,
                                          t_symtab*                    symtab,
                                          int                          atind,
                                          PreprocessResidue*           localPpResidue,
                                          const MoleculePatchDatabase& singlePatch,
                                          bool                         bVerbose)
{
    int   resnr;
    char* oldnm;
    int   anmnr;
    bool  bDeleted;

    oldnm = *pdba->atomname[atind];
    resnr = pdba->resinfo[pdba->atom[atind].resind].nr;

    bDeleted = FALSE;
    for (auto patch = singlePatch.hack.begin(); patch != singlePatch.hack.end(); patch++)
    {
        if (patch->type() == MoleculePatchType::Replace && gmx::equalCaseInsensitive(oldnm, patch->oname))
        {
            /* This is a replace entry. */
            /* Check if we are not replacing a replaced atom. */
            bool bReplaceReplace = false;
            for (auto selfPatch = singlePatch.hack.begin(); selfPatch != singlePatch.hack.end(); selfPatch++)
            {
                if (patch != selfPatch && selfPatch->type() == MoleculePatchType::Replace
                    && gmx::equalCaseInsensitive(selfPatch->nname, patch->oname))
                {
                    /* The replace in patch replaces an atom that
                     * was already replaced in selfPatch, we do not want
                     * second or higher level replaces at this stage.
                     */
                    bReplaceReplace = true;
                }
            }
            if (bReplaceReplace)
            {
                /* Skip this replace. */
                continue;
            }

            /* This atom still has the old name, rename it */
            std::string newnm = patch->nname;
            auto        found = std::find_if(
                    localPpResidue->atomname.begin(), localPpResidue->atomname.end(),
                    [&newnm](char** name) { return gmx::equalCaseInsensitive(newnm, *name); });
            if (found == localPpResidue->atomname.end())
            {
                /* The new name is not present in the rtp.
                 * We need to apply the replace also to the rtp entry.
                 */

                /* We need to find the add hack that can add this atom
                 * to find out after which atom it should be added.
                 */
                bool bFoundInAdd = false;
                for (auto rtpModification = singlePatch.hack.begin();
                     rtpModification != singlePatch.hack.end(); rtpModification++)
                {
                    int         k = std::distance(localPpResidue->atomname.begin(), found);
                    std::string start_at;
                    if (rtpModification->type() == MoleculePatchType::Add
                        && atomname_cmp_nr(newnm.c_str(), &(*rtpModification), &anmnr))
                    {
                        if (anmnr <= 1)
                        {
                            start_at = singlePatch.hack[k].a[0];
                        }
                        else
                        {
                            start_at = gmx::formatString("%s%d", singlePatch.hack[k].nname.c_str(),
                                                         anmnr - 1);
                        }
                        auto found2 =
                                std::find_if(localPpResidue->atomname.begin(),
                                             localPpResidue->atomname.end(), [&start_at](char** name) {
                                                 return gmx::equalCaseInsensitive(start_at, *name);
                                             });
                        if (found2 == localPpResidue->atomname.end())
                        {
                            gmx_fatal(FARGS,
                                      "Could not find atom '%s' in residue building block '%s' to "
                                      "add atom '%s' to",
                                      start_at.c_str(), localPpResidue->resname.c_str(), newnm.c_str());
                        }
                        /* We can add the atom after atom start_nr */
                        add_atom_to_restp(localPpResidue, symtab,
                                          std::distance(localPpResidue->atomname.begin(), found2),
                                          &(*patch));

                        bFoundInAdd = true;
                    }
                }

                if (!bFoundInAdd)
                {
                    gmx_fatal(FARGS,
                              "Could not find an 'add' entry for atom named '%s' corresponding to "
                              "the 'replace' entry from atom name '%s' to '%s' for tdb or hdb "
                              "database of residue type '%s'",
                              newnm.c_str(), patch->oname.c_str(), patch->nname.c_str(),
                              localPpResidue->resname.c_str());
                }
            }

            if (bVerbose)
            {
                printf("Renaming atom '%s' in residue '%s' %d to '%s'\n", oldnm,
                       localPpResidue->resname.c_str(), resnr, newnm.c_str());
            }
            /* Rename the atom in pdba */
            pdba->atomname[atind] = put_symtab(symtab, newnm.c_str());
        }
        else if (patch->type() == MoleculePatchType::Delete
                 && gmx::equalCaseInsensitive(oldnm, patch->oname))
        {
            /* This is a delete entry, check if this atom is present
             * in the rtp entry of this residue.
             */
            auto found3 = std::find_if(
                    localPpResidue->atomname.begin(), localPpResidue->atomname.end(),
                    [&oldnm](char** name) { return gmx::equalCaseInsensitive(oldnm, *name); });
            if (found3 == localPpResidue->atomname.end())
            {
                /* This atom is not present in the rtp entry,
                 * delete is now from the input pdba.
                 */
                if (bVerbose)
                {
                    printf("Deleting atom '%s' in residue '%s' %d\n", oldnm,
                           localPpResidue->resname.c_str(), resnr);
                }
                /* We should free the atom name,
                 * but it might be used multiple times in the symtab.
                 * sfree(pdba->atomname[atind]);
                 */
                for (int k = atind + 1; k < pdba->nr; k++)
                {
                    pdba->atom[k - 1]     = pdba->atom[k];
                    pdba->atomname[k - 1] = pdba->atomname[k];
                    copy_rvec(x[k], x[k - 1]);
                }
                pdba->nr--;
                bDeleted = true;
            }
        }
    }

    return bDeleted;
}

void match_atomnames_with_rtp(gmx::ArrayRef<PreprocessResidue>     usedPpResidues,
                              gmx::ArrayRef<MoleculePatchDatabase> globalPatches,
                              t_atoms*                             pdba,
                              t_symtab*                            symtab,
                              gmx::ArrayRef<gmx::RVec>             x,
                              bool                                 bVerbose)
{
    for (int i = 0; i < pdba->nr; i++)
    {
        const char*        oldnm          = *pdba->atomname[i];
        PreprocessResidue* localPpResidue = &usedPpResidues[pdba->atom[i].resind];
        auto               found          = std::find_if(
                localPpResidue->atomname.begin(), localPpResidue->atomname.end(),
                [&oldnm](char** name) { return gmx::equalCaseInsensitive(oldnm, *name); });
        if (found == localPpResidue->atomname.end())
        {
            /* Not found yet, check if we have to rename this atom */
            if (match_atomnames_with_rtp_atom(pdba, x, symtab, i, localPpResidue,
                                              globalPatches[pdba->atom[i].resind], bVerbose))
            {
                /* We deleted this atom, decrease the atom counter by 1. */
                i--;
            }
        }
    }
}

#define NUM_CMAP_ATOMS 5
static void gen_cmap(InteractionsOfType* psb, gmx::ArrayRef<const PreprocessResidue> usedPpResidues, t_atoms* atoms)
{
    int         residx;
    const char* ptr;
    t_resinfo*  resinfo = atoms->resinfo;
    int         nres    = atoms->nres;
    int         cmap_atomid[NUM_CMAP_ATOMS];
    int         cmap_chainnum = -1;

    if (debug)
    {
        ptr = "cmap";
    }
    else
    {
        ptr = "check";
    }

    fprintf(stderr, "Making cmap torsions...\n");
    int i = 0;
    /* Most cmap entries use the N atom from the next residue, so the last
     * residue should not have its CMAP entry in that case, but for things like
     * dipeptides we sometimes define a complete CMAP entry inside a residue,
     * and in this case we need to process everything through the last residue.
     */
    for (residx = 0; residx < nres; residx++)
    {
        /* Add CMAP terms from the list of CMAP interactions */
        for (const auto& b : usedPpResidues[residx].rb[ebtsCMAP].b)
        {
            bool bAddCMAP = true;
            /* Loop over atoms in a candidate CMAP interaction and
             * check that they exist, are from the same chain and are
             * from residues labelled as protein. */
            for (int k = 0; k < NUM_CMAP_ATOMS && bAddCMAP; k++)
            {
                /* Assign the pointer to the name of the next reference atom.
                 * This can use -/+ labels to refer to previous/next residue.
                 */
                const char* pname = b.a[k].c_str();
                /* Skip this CMAP entry if it refers to residues before the
                 * first or after the last residue.
                 */
                if (((strchr(pname, '-') != nullptr) && (residx == 0))
                    || ((strchr(pname, '+') != nullptr) && (residx == nres - 1)))
                {
                    bAddCMAP = false;
                    break;
                }

                cmap_atomid[k] = search_atom(pname, i, atoms, ptr, TRUE);
                bAddCMAP       = bAddCMAP && (cmap_atomid[k] != -1);
                if (!bAddCMAP)
                {
                    /* This break is necessary, because cmap_atomid[k]
                     * == -1 cannot be safely used as an index
                     * into the atom array. */
                    break;
                }
                int this_residue_index = atoms->atom[cmap_atomid[k]].resind;
                if (0 == k)
                {
                    cmap_chainnum = resinfo[this_residue_index].chainnum;
                }
                else
                {
                    /* Does the residue for this atom have the same
                     * chain number as the residues for previous
                     * atoms? */
                    bAddCMAP = bAddCMAP && cmap_chainnum == resinfo[this_residue_index].chainnum;
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
                add_cmap_param(psb, cmap_atomid[0], cmap_atomid[1], cmap_atomid[2], cmap_atomid[3],
                               cmap_atomid[4], b.s.c_str());
            }
        }

        if (residx < nres - 1)
        {
            while (atoms->atom[i].resind < residx + 1)
            {
                i++;
            }
        }
    }
    /* Start the next residue */
}

static void scrub_charge_groups(int* cgnr, int natoms)
{
    int i;

    for (i = 0; i < natoms; i++)
    {
        cgnr[i] = i + 1;
    }
}


void pdb2top(FILE*                                  top_file,
             const char*                            posre_fn,
             const char*                            molname,
             t_atoms*                               atoms,
             std::vector<gmx::RVec>*                x,
             PreprocessingAtomTypes*                atype,
             t_symtab*                              tab,
             gmx::ArrayRef<const PreprocessResidue> rtpFFDB,
             gmx::ArrayRef<PreprocessResidue>       usedPpResidues,
             gmx::ArrayRef<MoleculePatchDatabase>   globalPatches,
             bool                                   bAllowMissing,
             bool                                   bVsites,
             bool                                   bVsiteAromatics,
             const char*                            ffdir,
             real                                   mHmult,
             gmx::ArrayRef<const DisulfideBond>     ssbonds,
             real                                   long_bond_dist,
             real                                   short_bond_dist,
             bool                                   bDeuterate,
             bool                                   bChargeGroups,
             bool                                   bCmap,
             bool                                   bRenumRes,
             bool                                   bRTPresname)
{
    std::array<InteractionsOfType, F_NRE> plist;
    t_excls*                              excls;
    int*                                  cgnr;
    int*                                  vsite_type;
    int                                   i, nmissat;
    int                                   bts[ebtsNR];

    ResidueType rt;

    /* Make bonds */
    at2bonds(&(plist[F_BONDS]), globalPatches, atoms, *x, long_bond_dist, short_bond_dist);

    /* specbonds: disulphide bonds & heme-his */
    do_ssbonds(&(plist[F_BONDS]), atoms, ssbonds, bAllowMissing);

    nmissat = name2type(atoms, &cgnr, usedPpResidues, &rt);
    if (nmissat)
    {
        if (bAllowMissing)
        {
            fprintf(stderr, "There were %d missing atoms in molecule %s\n", nmissat, molname);
        }
        else
        {
            gmx_fatal(FARGS,
                      "There were %d missing atoms in molecule %s, if you want to use this "
                      "incomplete topology anyhow, use the option -missing",
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
            fprintf(stdout,
                    "The conversion of aromatic rings into virtual sites is deprecated "
                    "and may be removed in a future version of GROMACS");
        }
        /* determine which atoms will be vsites and add dummy masses
           also renumber atom numbers in plist[0..F_NRE]! */
        do_vsites(rtpFFDB, atype, atoms, tab, x, plist, &vsite_type, &cgnr, mHmult, bVsiteAromatics, ffdir);
    }

    /* Make Angles and Dihedrals */
    fprintf(stderr, "Generating angles, dihedrals and pairs...\n");
    snew(excls, atoms->nr);
    gen_pad(atoms, usedPpResidues, plist, excls, globalPatches, bAllowMissing);

    /* Make CMAP */
    if (bCmap)
    {
        gen_cmap(&(plist[F_CMAP]), usedPpResidues, atoms);
        if (plist[F_CMAP].size() > 0)
        {
            fprintf(stderr, "There are %4zu cmap torsion pairs\n", plist[F_CMAP].size());
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
            "There are %4zu dihedrals, %4zu impropers, %4zu angles\n"
            "          %4zu pairs,     %4zu bonds and  %4zu virtual sites\n",
            plist[F_PDIHS].size(), plist[F_IDIHS].size(), plist[F_ANGLES].size(),
            plist[F_LJ14].size(), plist[F_BONDS].size(),
            plist[F_VSITE2].size() + plist[F_VSITE3].size() + plist[F_VSITE3FD].size()
                    + plist[F_VSITE3FAD].size() + plist[F_VSITE3OUT].size()
                    + plist[F_VSITE4FD].size() + plist[F_VSITE4FDN].size());

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
            bts[i] = usedPpResidues[0].rb[i].type;
        }
        write_top(top_file, posre_fn, molname, atoms, bRTPresname, bts, plist, excls, atype, cgnr,
                  usedPpResidues[0].nrexcl);
    }


    /* we should clean up hb and restp here, but that is a *L*O*T* of work! */
    sfree(cgnr);
    for (i = 0; i < atoms->nr; i++)
    {
        sfree(excls[i].e);
    }
    sfree(excls);
}

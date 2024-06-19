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

#include "ter_db.h"

#include <cctype>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <array>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxpreprocess/fflibutil.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/h_db.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringtoenumvalueconverter.h"
#include "gromacs/utility/stringutil.h"

#include "hackblock.h"
#include "resall.h"

/* use bonded types definitions in hackblock.h */
enum class ReplaceType : int
{
    Repl,
    Add,
    Del,
    Count
};

static const char* enumValueToString(ReplaceType enumValue)
{
    constexpr gmx::EnumerationArray<ReplaceType, const char*> replaceTypeNames = { "replace",
                                                                                   "add",
                                                                                   "delete" };
    return replaceTypeNames[enumValue];
}

template<typename EnumType>
static std::optional<EnumType> findTypeFromKeyword(char* keyw)
{
    gmx::StringToEnumValueConverter<EnumType, enumValueToString, gmx::StringCompareType::CaseInsensitive, gmx::StripStrings::Yes> converter;
    return converter.valueFrom(keyw);
}

static void read_atom(char* line, bool bAdd, std::string* nname, t_atom* a, PreprocessingAtomTypes* atype, int* cgnr)
{
    int    nr, i;
    char   buf[5][30];
    double m, q;

    /* This code is messy, because of support for different formats:
     * for replace: [new name] <atom type> <m> <q> [cgnr (old format)]
     * for add:                <atom type> <m> <q> [cgnr]
     */
    nr = sscanf(line, "%s %s %s %s %s", buf[0], buf[1], buf[2], buf[3], buf[4]);

    /* Here there an ambiguity due to the old replace format with cgnr,
     * which was read for years, but ignored in the rest of the code.
     * We have to assume that the atom type does not start with a digit
     * to make a line with 4 entries uniquely interpretable.
     */
    if (!bAdd && nr == 4 && isdigit(buf[1][0]))
    {
        nr = 3;
    }

    if (nr < 3 || nr > 4)
    {
        gmx_fatal(FARGS,
                  "Reading Termini Database: expected %d or %d items of atom data in stead of %d "
                  "on line\n%s",
                  3,
                  4,
                  nr,
                  line);
    }
    i = 0;
    if (!bAdd)
    {
        if (nr == 4)
        {
            *nname = buf[i++];
        }
        else
        {
            *nname = "";
        }
    }
    auto atomType = atype->atomTypeFromName(buf[i++]);
    if (atomType == std::nullopt)
    {
        GMX_THROW(gmx::InconsistentInputError(
                gmx::formatString("Atom type %s specified in terminal database has not been "
                                  "defined in the force field",
                                  buf[i - 1])));
    }
    else
    {
        a->type = atomType.value();
    }
    sscanf(buf[i++], "%lf", &m);
    a->m = m;
    sscanf(buf[i++], "%lf", &q);
    a->q = q;
    if (bAdd && nr == 4)
    {
        sscanf(buf[i++], "%d", cgnr);
    }
    else
    {
        *cgnr = NOTSET;
    }
}

static void print_atom(FILE* out, const t_atom& a, PreprocessingAtomTypes* atype)
{
    fprintf(out, "\t%s\t%g\t%g\n", atype->atomNameFromAtomType(a.type)->c_str(), a.m, a.q);
}

static void print_ter_db(const char*                                ff,
                         char                                       C,
                         gmx::ArrayRef<const MoleculePatchDatabase> tb,
                         PreprocessingAtomTypes*                    atype)
{
    std::string buf = gmx::formatString("%s-%c.tdb", ff, C);
    FILE*       out = gmx_fio_fopen(buf.c_str(), "w");

    for (const auto& modification : tb)
    {
        fprintf(out, "[ %s ]\n", modification.name.c_str());

        if (std::any_of(modification.hack.begin(), modification.hack.end(), [](const auto& mod) {
                return mod.type() == MoleculePatchType::Replace;
            }))
        {
            fprintf(out, "[ %s ]\n", enumValueToString(ReplaceType::Repl));
            for (const auto& hack : modification.hack)
            {
                if (hack.type() == MoleculePatchType::Replace)
                {
                    fprintf(out, "%s\t", hack.oname.c_str());
                    print_atom(out, hack.atom.back(), atype);
                }
            }
        }
        if (std::any_of(modification.hack.begin(), modification.hack.end(), [](const auto& mod) {
                return mod.type() == MoleculePatchType::Add;
            }))
        {
            fprintf(out, "[ %s ]\n", enumValueToString(ReplaceType::Add));
            for (const auto& hack : modification.hack)
            {
                if (hack.type() == MoleculePatchType::Add)
                {
                    print_ab(out, hack, hack.nname.c_str());
                    print_atom(out, hack.atom.back(), atype);
                }
            }
        }
        if (std::any_of(modification.hack.begin(), modification.hack.end(), [](const auto& mod) {
                return mod.type() == MoleculePatchType::Delete;
            }))
        {
            fprintf(out, "[ %s ]\n", enumValueToString(ReplaceType::Del));
            for (const auto& hack : modification.hack)
            {
                if (hack.type() == MoleculePatchType::Delete)
                {
                    fprintf(out, "%s\n", hack.oname.c_str());
                }
            }
        }
        for (auto bt : gmx::EnumerationWrapper<BondedTypes>{})
        {
            if (!modification.rb[bt].b.empty())
            {
                fprintf(out, "[ %s ]\n", enumValueToString(bt));
                for (const auto& b : modification.rb[bt].b)
                {
                    for (int k = 0; k < enumValueToNumIAtoms(bt); k++)
                    {
                        fprintf(out, "%s%s", k ? "\t" : "", b.a[k].c_str());
                    }
                    if (!b.s.empty())
                    {
                        fprintf(out, "\t%s", b.s.c_str());
                    }
                    fprintf(out, "\n");
                }
            }
        }
        fprintf(out, "\n");
    }
    gmx_fio_fclose(out);
}

static void read_ter_db_file(const std::filesystem::path&        fn,
                             std::vector<MoleculePatchDatabase>* tbptr,
                             PreprocessingAtomTypes*             atype)
{
    char header[STRLEN], buf[STRLEN], line[STRLEN];

    auto filebase = fflib_filename_base(fn);
    filebase.replace_extension();

    FILE* in = fflib_open(fn);

    std::optional<BondedTypes> btkw;
    std::optional<ReplaceType> rtkw;
    get_a_line(in, line, STRLEN);
    MoleculePatchDatabase* block = nullptr;
    while (!feof(in))
    {
        if (get_header(line, header))
        {
            /* this is a new block, or a new keyword */
            btkw = findTypeFromKeyword<BondedTypes>(header);
            rtkw = findTypeFromKeyword<ReplaceType>(header);

            if (!btkw.has_value() && !rtkw.has_value())
            {
                tbptr->emplace_back();
                block = &tbptr->back();
                clearModificationBlock(block);
                block->name     = header;
                block->filebase = filebase.string();
            }
        }
        else
        {
            if (block == nullptr)
            {
                gmx_fatal(FARGS,
                          "reading termini database: "
                          "directive expected before line:\n%s\n"
                          "This might be a file in an old format.",
                          line);
            }
            /* this is not a header, so it must be data */
            if (!btkw.has_value())
            {
                /* this is a hack: add/rename/delete atoms */
                /* make space for hacks */
                block->hack.emplace_back();
                MoleculePatch* hack = &block->hack.back();

                /* get data */
                int n = 0;
                if (*rtkw == ReplaceType::Repl || *rtkw == ReplaceType::Del)
                {
                    if (sscanf(line, "%s%n", buf, &n) != 1)
                    {
                        gmx_fatal(FARGS,
                                  "Reading Termini Database '%s': "
                                  "expected atom name on line\n%s",
                                  fn.string().c_str(),
                                  line);
                    }
                    hack->oname = buf;
                    /* we only replace or delete one atom at a time */
                    hack->nr = 1;
                }
                else if (*rtkw == ReplaceType::Add)
                {
                    read_ab(line, fn, hack);
                    get_a_line(in, line, STRLEN);
                }
                else
                {
                    gmx_fatal(FARGS,
                              "unimplemented keyword number %d (%s:%d)",
                              static_cast<int>(*rtkw),
                              __FILE__,
                              __LINE__);
                }
                if (*rtkw == ReplaceType::Repl || *rtkw == ReplaceType::Add)
                {
                    hack->atom.emplace_back();
                    read_atom(line + n, *rtkw == ReplaceType::Add, &hack->nname, &hack->atom.back(), atype, &hack->cgnr);
                    if (hack->nname.empty())
                    {
                        if (!hack->oname.empty())
                        {
                            hack->nname = hack->oname;
                        }
                        else
                        {
                            gmx_fatal(FARGS,
                                      "Reading Termini Database '%s': don't know which name the "
                                      "new atom should have on line\n%s",
                                      fn.string().c_str(),
                                      line);
                        }
                    }
                }
            }
            else if (*btkw >= BondedTypes::Bonds && *btkw < BondedTypes::Count)
            {
                /* this is bonded data: bonds, angles, dihedrals or impropers */
                int n = 0;
                block->rb[*btkw].b.emplace_back();
                BondedInteraction* newBond = &block->rb[*btkw].b.back();
                for (int j = 0; j < enumValueToNumIAtoms(*btkw); j++)
                {
                    int ni;
                    if (sscanf(line + n, "%s%n", buf, &ni) == 1)
                    {
                        newBond->a[j] = buf;
                    }
                    else
                    {
                        gmx_fatal(FARGS,
                                  "Reading Termini Database '%s': expected %d atom names (found "
                                  "%d) on line\n%s",
                                  fn.string().c_str(),
                                  enumValueToNumIAtoms(*btkw),
                                  j - 1,
                                  line);
                    }
                    n += ni;
                }
                strcpy(buf, "");
                sscanf(line + n, "%s", buf);
                newBond->s = buf;
            }
            else
            {
                gmx_fatal(FARGS,
                          "Reading Termini Database: Expecting a header at line\n"
                          "%s",
                          line);
            }
        }
        get_a_line(in, line, STRLEN);
    }

    gmx_ffclose(in);
}

int read_ter_db(const std::filesystem::path&        ffdir,
                char                                ter,
                std::vector<MoleculePatchDatabase>* tbptr,
                PreprocessingAtomTypes*             atype)
{
    std::string ext = gmx::formatString(".%c.tdb", ter);

    /* Search for termini database files.
     * Do not generate an error when none are found.
     */
    auto tdbf = fflib_search_file_end(ffdir, ext.c_str(), FALSE);
    tbptr->clear();
    for (const auto& filename : tdbf)
    {
        read_ter_db_file(filename, tbptr, atype);
    }

    if (debug)
    {
        print_ter_db("new", ter, *tbptr, atype);
    }

    return tbptr->size();
}

std::vector<MoleculePatchDatabase*> filter_ter(gmx::ArrayRef<MoleculePatchDatabase> tb, const char* resname)
{
    // TODO Four years later, no force fields have ever used this, so decide status of this feature
    /* Since some force fields (e.g. OPLS) needs different
     * atomtypes for different residues there could be a lot
     * of entries in the databases for specific residues
     * (e.g. GLY-NH3+, SER-NH3+, ALA-NH3+).
     *
     * To reduce the database size, we assume that a terminus specifier liker
     *
     * [ GLY|SER|ALA-NH3+ ]
     *
     * would cover all of the three residue types above.
     * Wildcards (*,?) are not OK. Don't worry about e.g. GLU vs. GLUH since
     * pdb2gmx only uses the first 3 letters when calling this routine.
     *
     * To automate this, this routines scans a list of termini
     * for the residue name "resname" and returns an allocated list of
     * pointers to the termini that could be applied to the
     * residue in question. The variable pointed to by nret will
     * contain the number of valid pointers in the list.
     * Remember to free the list when you are done with it...
     */

    auto                                none_idx = tb.end();
    std::vector<MoleculePatchDatabase*> list;

    for (auto it = tb.begin(); it != tb.end(); it++)
    {
        const char* s     = it->name.c_str();
        bool        found = false;
        do
        {
            if (gmx::equalCaseInsensitive(resname, s, 3))
            {
                found = true;
                list.push_back(&*it);
            }
            else
            {
                /* advance to next |-separated field */
                s = strchr(s, '|');
                if (s != nullptr)
                {
                    s++;
                }
            }
        } while (!found && s != nullptr);
    }

    /* All residue-specific termini have been added. We might have to fall
     * back on generic termini, which are characterized by not having
     * '-' in the name prior to the last position (which indicates charge).
     * The [ None ] alternative is special since we don't want that
     * to be the default, so we put it last in the list we return.
     */
    for (auto it = tb.begin(); it != tb.end(); it++)
    {
        const char* s = it->name.c_str();
        if (gmx::equalCaseInsensitive("None", it->name))
        {
            none_idx = it;
        }
        else
        {
            /* Time to see if there's a generic terminus that matches.
               Is there a hyphen? */
            const char* c = strchr(s, '-');

            /* A conjunction hyphen normally indicates a residue-specific
               terminus, which is named like "GLY-COOH". A generic terminus
               won't have a hyphen. */
            bool bFoundAnyHyphen = (c != nullptr);
            /* '-' as the last character indicates charge, so if that's
               the only one found e.g. "COO-", then it was not a conjunction
               hyphen, so this is a generic terminus */
            bool bOnlyFoundChargeHyphen = (bFoundAnyHyphen && *(c + 1) == '\0');
            /* Thus, "GLY-COO-" is not recognized as a generic terminus. */
            bool bFoundGenericTerminus = !bFoundAnyHyphen || bOnlyFoundChargeHyphen;
            if (bFoundGenericTerminus)
            {
                /* Check that we haven't already added a residue-specific version
                 * of this terminus.
                 */
                auto found = std::find_if(list.begin(), list.end(), [&s](const MoleculePatchDatabase* b) {
                    return strstr(b->name.c_str(), s) != nullptr;
                });
                if (found == list.end())
                {
                    list.push_back(&*it);
                }
            }
        }
    }
    if (none_idx != tb.end())
    {
        list.push_back(&*none_idx);
    }

    return list;
}


MoleculePatchDatabase* choose_ter(gmx::ArrayRef<MoleculePatchDatabase*> tb, const char* title)
{
    int sel, ret;

    printf("%s\n", title);
    int i = 0;
    for (const auto& modification : tb)
    {
        bool bIsZwitterion = (0 == gmx_wcmatch("*ZWITTERION*", modification->name.c_str()));
        printf("%2d: %s%s\n",
               i,
               modification->name.c_str(),
               bIsZwitterion ? " (only use with zwitterions containing exactly one residue)" : "");
        i++;
    }
    do
    {
        ret = fscanf(stdin, "%d", &sel);
    } while ((ret != 1) || (sel < 0) || (sel >= tb.ssize()));

    return tb[sel];
}

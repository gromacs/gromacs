/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020,2021, by the GROMACS development team, led by
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

#include "index.h"

#include <cassert>
#include <cctype>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <numeric>

#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/invblock.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

static gmx_bool gmx_ask_yesno(gmx_bool bASK)
{
    char c = 0;

    if (bASK)
    {
        do
        {
            c = toupper(fgetc(stdin));
        } while ((c != 'Y') && (c != 'N'));

        return (c == 'Y');
    }
    else
    {
        return FALSE;
    }
}

void write_index(const char* outf, t_blocka* b, char** gnames, gmx_bool bDuplicate, int natoms)
{
    FILE* out = gmx_ffopen(outf, "w");
    /* fprintf(out,"%5d  %5d\n",b->nr,b->nra); */
    for (int i = 0; (i < b->nr); i++)
    {
        fprintf(out, "[ %s ]", gnames[i]);
        for (int k = 0, j = b->index[i]; j < b->index[i + 1]; j++, k++)
        {
            const char sep = (k % 15 == 0 ? '\n' : ' ');
            fprintf(out, "%c%4d", sep, b->a[j] + 1);
        }
        fprintf(out, "\n");
    }

    /* Duplicate copy, useful for computational electrophysiology double-layer setups */
    if (bDuplicate)
    {
        fprintf(stderr, "Duplicating the whole system with an atom offset of %d atoms.\n", natoms);
        for (int i = 0; (i < b->nr); i++)
        {
            fprintf(out, "[ %s_copy ]", gnames[i]);
            for (int k = 0, j = b->index[i]; j < b->index[i + 1]; j++, k++)
            {
                const char sep = (k % 15 == 0 ? '\n' : ' ');
                fprintf(out, "%c%4d", sep, b->a[j] + 1 + natoms);
            }
            fprintf(out, "\n");
        }
    }

    gmx_ffclose(out);
}

void add_grp(t_blocka* b, char*** gnames, gmx::ArrayRef<const int> a, const std::string& name)
{
    srenew(b->index, b->nr + 2);
    srenew(*gnames, b->nr + 1);
    (*gnames)[b->nr] = gmx_strdup(name.c_str());

    srenew(b->a, b->nra + a.size());
    for (gmx::index i = 0; (i < a.ssize()); i++)
    {
        b->a[b->nra++] = a[i];
    }
    b->nr++;
    b->index[b->nr] = b->nra;
}

/*! \brief
 * Compare index in groups.
 *
 * Checks the index in \p a to the one in \p b at \p index.
 * If \p index is < 0, it is taken as relative to the end of \p b.
 *
 * \param[in] b Block with groups.
 * \param[in] a New group to compare to.
 * \param[in] index The index to check.
 * \returns True if groups are the same.
 */
static bool grp_cmp(t_blocka* b, gmx::ArrayRef<const int> a, int index)
{
    if (index < 0)
    {
        index = b->nr - 1 + index;
    }
    if (index >= b->nr)
    {
        gmx_fatal(FARGS, "no such index group %d in t_blocka (nr=%d)", index, b->nr);
    }
    /* compare sizes */
    if (a.ssize() != b->index[index + 1] - b->index[index])
    {
        return FALSE;
    }
    for (gmx::index i = 0; i < a.ssize(); i++)
    {
        if (a[i] != b->a[b->index[index] + i])
        {
            return false;
        }
    }
    return true;
}
//! Print out how many residues of a certain type are present.
static void p_status(gmx::ArrayRef<const std::string> restype, gmx::ArrayRef<const std::string> typenames)
{
    std::vector<int> counter(typenames.size(), 0);
    std::transform(typenames.begin(), typenames.end(), counter.begin(), [&restype](const std::string& typenm) {
        return std::count_if(restype.begin(), restype.end(), [&typenm](const std::string& res) {
            return gmx_strcasecmp(res.c_str(), typenm.c_str()) == 0;
        });
    });

    for (int i = 0; (i < typenames.ssize()); i++)
    {
        if (counter[i] > 0)
        {
            printf("There are: %5d %10s residues\n", counter[i], typenames[i].c_str());
        }
    }
}

/*! \brief
 * Return a new group of atoms with \p restype that match \p typestring,
 * or not matching if \p bMatch.
 *
 * \param[in] atoms Atoms data to use.
 * \param[in] restype Residuetypes to match.
 * \param[in] typestring Which type the residue should have.
 * \param[in] bMatch whether to return matching atoms or those that don't.
 * \returns Vector of atoms that match.
 */
static std::vector<int> mk_aid(const t_atoms*                   atoms,
                               gmx::ArrayRef<const std::string> restype,
                               const std::string&               typestring,
                               bool                             bMatch)
{
    std::vector<int> a;
    for (int i = 0; (i < atoms->nr); i++)
    {
        bool res = gmx_strcasecmp(restype[atoms->atom[i].resind].c_str(), typestring.c_str()) == 0;
        if (!bMatch)
        {
            res = !res;
        }
        if (res)
        {
            a.push_back(i);
        }
    }

    return a;
}

typedef struct
{
    char*    rname;
    gmx_bool bNeg;
    char*    gname;
} restp_t;

static void analyse_other(gmx::ArrayRef<std::string> restype,
                          const t_atoms*             atoms,
                          t_blocka*                  gb,
                          char***                    gn,
                          gmx_bool                   bASK,
                          gmx_bool                   bVerb)
{
    std::vector<restp_t> restp;
    int                  i = 0;

    for (; (i < atoms->nres); i++)
    {
        if (gmx_strcasecmp(restype[i].c_str(), "Protein")
            && gmx_strcasecmp(restype[i].c_str(), "DNA") && gmx_strcasecmp(restype[i].c_str(), "RNA")
            && gmx_strcasecmp(restype[i].c_str(), "Water"))
        {
            break;
        }
    }
    if (i < atoms->nres)
    {
        /* we have others */
        if (bVerb)
        {
            printf("Analysing residues not classified as Protein/DNA/RNA/Water and splitting into "
                   "groups...\n");
        }
        for (int k = 0; (k < atoms->nr); k++)
        {
            int         resind = atoms->atom[k].resind;
            const char* rname  = *atoms->resinfo[resind].name;
            if (gmx_strcasecmp(restype[resind].c_str(), "Protein")
                && gmx_strcasecmp(restype[resind].c_str(), "DNA")
                && gmx_strcasecmp(restype[resind].c_str(), "RNA")
                && gmx_strcasecmp(restype[resind].c_str(), "Water"))
            {
                auto found = std::find_if(restp.begin(), restp.end(), [rname](const auto& entry) {
                    return strcmp(entry.rname, rname) == 0;
                });
                if (found == restp.end())
                {
                    restp.emplace_back();
                    auto& last = restp.back();
                    last.rname = gmx_strdup(rname);
                    last.bNeg  = false;
                    last.gname = gmx_strdup(rname);
                }
            }
        }
        for (int i = 0; (i < gmx::ssize(restp)); i++)
        {
            std::vector<int> aid;
            for (int j = 0; (j < atoms->nr); j++)
            {
                const char* rname = *atoms->resinfo[atoms->atom[j].resind].name;
                if ((strcmp(restp[i].rname, rname) == 0 && !restp[i].bNeg)
                    || (strcmp(restp[i].rname, rname) != 0 && restp[i].bNeg))
                {
                    aid.push_back(j);
                }
            }
            add_grp(gb, gn, aid, restp[i].gname);
            if (bASK)
            {
                printf("split %s into atoms (y/n) ? ", restp[i].gname);
                fflush(stdout);
                if (gmx_ask_yesno(bASK))
                {
                    std::vector<const char*> attp;
                    for (size_t k = 0; (k < aid.size()); k++)
                    {
                        const char* aname = *atoms->atomname[aid[k]];
                        auto found = std::find_if(attp.begin(), attp.end(), [aname](const char* entry) {
                            return strcmp(aname, entry) == 0;
                        });
                        if (found == attp.end())
                        {
                            attp.emplace_back(aname);
                        }
                    }
                    if (attp.size() > 1)
                    {
                        const int natp = attp.size();
                        for (int l = 0; (l < natp); l++)
                        {
                            std::vector<int> aaid;
                            for (size_t k = 0; (k < aid.size()); k++)
                            {
                                const char* aname = *atoms->atomname[aid[k]];
                                if (strcmp(aname, attp[l]) == 0)
                                {
                                    aaid.push_back(aid[k]);
                                }
                            }
                            add_grp(gb, gn, aaid, attp[l]);
                        }
                    }
                }
            }
            sfree(restp[i].rname);
            sfree(restp[i].gname);
        }
    }
}

/*! \brief
 * Data necessary to construct a single (protein) index group in
 * analyse_prot().
 */
typedef struct gmx_help_make_index_group // NOLINT(clang-analyzer-optin.performance.Padding)
{
    /** The set of atom names that will be used to form this index group */
    const char** defining_atomnames;
    /** Size of the defining_atomnames array */
    int num_defining_atomnames;
    /** Name of this index group */
    const char* group_name;
    /** Whether the above atom names name the atoms in the group, or
        those not in the group */
    gmx_bool bTakeComplement;
    /** The index in wholename gives the first item in the arrays of
       atomnames that should be tested with 'gmx_strncasecmp' in stead of
       gmx_strcasecmp, or -1 if all items should be tested with strcasecmp
       This is comparable to using a '*' wildcard at the end of specific
       atom names, but that is more involved to implement...
     */
    int wholename;
    /** Only create this index group if it differs from the one specified in compareto,
       where -1 means to always create this group. */
    int compareto;
} t_gmx_help_make_index_group;

static void analyse_prot(gmx::ArrayRef<const std::string> restype,
                         const t_atoms*                   atoms,
                         t_blocka*                        gb,
                         char***                          gn,
                         gmx_bool                         bASK,
                         gmx_bool                         bVerb)
{
    /* lists of atomnames to be used in constructing index groups: */
    static const char* pnoh[]   = { "H", "HN" };
    static const char* pnodum[] = { "MN1",  "MN2",  "MCB1", "MCB2", "MCG1", "MCG2",
                                    "MCD1", "MCD2", "MCE1", "MCE2", "MNZ1", "MNZ2" };
    static const char* calpha[] = { "CA" };
    static const char* bb[]     = { "N", "CA", "C" };
    static const char* mc[]     = { "N", "CA", "C", "O", "O1", "O2", "OC1", "OC2", "OT", "OXT" };
    static const char* mcb[] = { "N", "CA", "CB", "C", "O", "O1", "O2", "OC1", "OC2", "OT", "OXT" };
    static const char* mch[] = { "N",  "CA",  "C",  "O",  "O1", "O2", "OC1", "OC2",
                                 "OT", "OXT", "H1", "H2", "H3", "H",  "HN" };

    static const t_gmx_help_make_index_group constructing_data[] = {
        { nullptr, 0, "Protein", TRUE, -1, -1 },
        { pnoh, asize(pnoh), "Protein-H", TRUE, 0, -1 },
        { calpha, asize(calpha), "C-alpha", FALSE, -1, -1 },
        { bb, asize(bb), "Backbone", FALSE, -1, -1 },
        { mc, asize(mc), "MainChain", FALSE, -1, -1 },
        { mcb, asize(mcb), "MainChain+Cb", FALSE, -1, -1 },
        { mch, asize(mch), "MainChain+H", FALSE, -1, -1 },
        { mch, asize(mch), "SideChain", TRUE, -1, -1 },
        { mch, asize(mch), "SideChain-H", TRUE, 11, -1 },
        { pnodum, asize(pnodum), "Prot-Masses", TRUE, -1, 0 },
    };
    const int num_index_groups = asize(constructing_data);

    char ndx_name[STRLEN];

    if (bVerb)
    {
        printf("Analysing Protein...\n");
    }
    std::vector<int> aid;

    /* calculate the number of protein residues */
    int npres = 0;
    for (int i = 0; (i < atoms->nres); i++)
    {
        if (0 == gmx_strcasecmp(restype[i].c_str(), "Protein"))
        {
            npres++;
        }
    }
    /* find matching or complement atoms */
    for (int i = 0; (i < num_index_groups); i++)
    {
        for (int n = 0; (n < atoms->nr); n++)
        {
            if (0 == gmx_strcasecmp(restype[atoms->atom[n].resind].c_str(), "Protein"))
            {
                bool match = false;
                for (int j = 0; (j < constructing_data[i].num_defining_atomnames); j++)
                {
                    /* skip digits at beginning of atomname, e.g. 1H */
                    char* atnm = *atoms->atomname[n];
                    while (isdigit(atnm[0]))
                    {
                        atnm++;
                    }
                    if ((constructing_data[i].wholename == -1) || (j < constructing_data[i].wholename))
                    {
                        if (0 == gmx_strcasecmp(constructing_data[i].defining_atomnames[j], atnm))
                        {
                            match = true;
                        }
                    }
                    else
                    {
                        if (0
                            == gmx_strncasecmp(constructing_data[i].defining_atomnames[j],
                                               atnm,
                                               strlen(constructing_data[i].defining_atomnames[j])))
                        {
                            match = true;
                        }
                    }
                }
                if (constructing_data[i].bTakeComplement != match)
                {
                    aid.push_back(n);
                }
            }
        }
        /* if we want to add this group always or it differs from previous
           group, add it: */
        if (-1 == constructing_data[i].compareto || !grp_cmp(gb, aid, constructing_data[i].compareto - i))
        {
            add_grp(gb, gn, aid, constructing_data[i].group_name);
        }
        aid.clear();
    }

    if (bASK)
    {
        for (int i = 0; (i < num_index_groups); i++)
        {
            printf("Split %12s into %5d residues (y/n) ? ", constructing_data[i].group_name, npres);
            if (gmx_ask_yesno(bASK))
            {
                aid.clear();
                for (int n = 0; ((atoms->atom[n].resind < npres) && (n < atoms->nr));)
                {
                    int resind = atoms->atom[n].resind;
                    for (; ((atoms->atom[n].resind == resind) && (n < atoms->nr)); n++)
                    {
                        bool match = false;
                        for (int j = 0; (j < constructing_data[i].num_defining_atomnames); j++)
                        {
                            if (0
                                == gmx_strcasecmp(constructing_data[i].defining_atomnames[j],
                                                  *atoms->atomname[n]))
                            {
                                match = true;
                            }
                        }
                        if (constructing_data[i].bTakeComplement != match)
                        {
                            aid.push_back(n);
                        }
                    }
                    /* copy the residuename to the tail of the groupname */
                    if (!aid.empty())
                    {
                        t_resinfo* ri = &atoms->resinfo[resind];
                        sprintf(ndx_name,
                                "%s_%s%d%c",
                                constructing_data[i].group_name,
                                *ri->name,
                                ri->nr,
                                ri->ic == ' ' ? '\0' : ri->ic);
                        add_grp(gb, gn, aid, ndx_name);
                        aid.clear();
                    }
                }
            }
        }
        printf("Make group with sidechain and C=O swapped (y/n) ? ");
        if (gmx_ask_yesno(bASK))
        {
            /* Make swap sidechain C=O index */
            aid.clear();
            for (int n = 0; ((atoms->atom[n].resind < npres) && (n < atoms->nr));)
            {
                int resind = atoms->atom[n].resind;
                int hold   = -1;
                for (; ((atoms->atom[n].resind == resind) && (n < atoms->nr)); n++)
                {
                    if (strcmp("CA", *atoms->atomname[n]) == 0)
                    {
                        aid.push_back(n);
                        hold = aid.size();
                        aid.resize(aid.size() + 3);
                    }
                    else if (strcmp("C", *atoms->atomname[n]) == 0)
                    {
                        if (hold == -1)
                        {
                            gmx_incons("Atom naming problem");
                        }
                        aid[hold] = n;
                    }
                    else if (strcmp("O", *atoms->atomname[n]) == 0)
                    {
                        if (hold == -1)
                        {
                            gmx_incons("Atom naming problem");
                        }
                        aid[hold + 1] = n;
                    }
                    else if (strcmp("O1", *atoms->atomname[n]) == 0)
                    {
                        if (hold == -1)
                        {
                            gmx_incons("Atom naming problem");
                        }
                        aid[hold + 1] = n;
                    }
                    else
                    {
                        aid.push_back(n);
                    }
                }
            }
            /* copy the residuename to the tail of the groupname */
            if (!aid.empty())
            {
                add_grp(gb, gn, aid, "SwapSC-CO");
            }
        }
    }
}


void analyse(const t_atoms* atoms, t_blocka* gb, char*** gn, gmx_bool bASK, gmx_bool bVerb)
{
    if (bVerb)
    {
        printf("Analysing residue names:\n");
    }
    /* Create system group, every single atom */
    std::vector<int> aid(atoms->nr);
    std::iota(aid.begin(), aid.end(), 0);
    add_grp(gb, gn, aid, "System");

    /* For every residue, get a pointer to the residue type name */
    ResidueTypeMap residueTypeMap = residueTypeMapFromLibraryFile("residuetypes.dat");

    std::vector<std::string> restype;
    std::vector<std::string> previousTypename;
    if (atoms->nres > 0)
    {
        const char* resnm = *atoms->resinfo[0].name;
        restype.emplace_back(typeOfNamedDatabaseResidue(residueTypeMap, resnm));
        previousTypename.push_back(restype.back());

        for (int i = 1; i < atoms->nres; i++)
        {
            const char* resnm = *atoms->resinfo[i].name;
            restype.emplace_back(typeOfNamedDatabaseResidue(residueTypeMap, resnm));

            /* Note that this does not lead to a N*N loop, but N*K, where
             * K is the number of residue _types_, which is small and independent of N.
             */
            bool found = false;
            for (size_t k = 0; k < previousTypename.size() && !found; k++)
            {
                found = strcmp(restype.back().c_str(), previousTypename[k].c_str()) == 0;
            }
            if (!found)
            {
                previousTypename.push_back(restype.back());
            }
        }
    }

    if (bVerb)
    {
        p_status(restype, previousTypename);
    }

    for (gmx::index k = 0; k < gmx::ssize(previousTypename); k++)
    {
        std::vector<int> aid = mk_aid(atoms, restype, previousTypename[k], TRUE);

        /* Check for special types to do fancy stuff with */

        if (!gmx_strcasecmp(previousTypename[k].c_str(), "Protein") && !aid.empty())
        {
            /* PROTEIN */
            analyse_prot(restype, atoms, gb, gn, bASK, bVerb);

            /* Create a Non-Protein group */
            std::vector<int> aid = mk_aid(atoms, restype, "Protein", FALSE);
            if ((!aid.empty()) && (gmx::ssize(aid) < atoms->nr))
            {
                add_grp(gb, gn, aid, "non-Protein");
            }
        }
        else if (!gmx_strcasecmp(previousTypename[k].c_str(), "Water") && !aid.empty())
        {
            add_grp(gb, gn, aid, previousTypename[k]);
            /* Add this group as 'SOL' too, for backward compatibility with older gromacs versions */
            add_grp(gb, gn, aid, "SOL");


            /* Solvent, create a negated group too */
            std::vector<int> aid = mk_aid(atoms, restype, "Water", FALSE);
            if ((!aid.empty()) && (gmx::ssize(aid) < atoms->nr))
            {
                add_grp(gb, gn, aid, "non-Water");
            }
        }
        else if (!aid.empty())
        {
            /* Other groups */
            add_grp(gb, gn, aid, previousTypename[k]);
            analyse_other(restype, atoms, gb, gn, bASK, bVerb);
        }
    }


    /* Create a merged water_and_ions group */
    int iwater = -1;
    int iion   = -1;
    int nwater = 0;
    int nion   = 0;

    for (int i = 0; i < gb->nr; i++)
    {
        if (!gmx_strcasecmp((*gn)[i], "Water"))
        {
            iwater = i;
            nwater = gb->index[i + 1] - gb->index[i];
        }
        else if (!gmx_strcasecmp((*gn)[i], "Ion"))
        {
            iion = i;
            nion = gb->index[i + 1] - gb->index[i];
        }
    }

    if (nwater > 0 && nion > 0)
    {
        srenew(gb->index, gb->nr + 2);
        srenew(*gn, gb->nr + 1);
        (*gn)[gb->nr] = gmx_strdup("Water_and_ions");
        srenew(gb->a, gb->nra + nwater + nion);
        if (nwater > 0)
        {
            for (int i = gb->index[iwater]; i < gb->index[iwater + 1]; i++)
            {
                gb->a[gb->nra++] = gb->a[i];
            }
        }
        if (nion > 0)
        {
            for (int i = gb->index[iion]; i < gb->index[iion + 1]; i++)
            {
                gb->a[gb->nra++] = gb->a[i];
            }
        }
        gb->nr++;
        gb->index[gb->nr] = gb->nra;
    }
}


void check_index(const char* gname, int n, int index[], const char* traj, int natoms)
{
    for (int i = 0; i < n; i++)
    {
        if (index[i] >= natoms)
        {
            gmx_fatal(FARGS,
                      "%s atom number (index[%d]=%d) is larger than the number of atoms in %s (%d)",
                      gname ? gname : "Index",
                      i + 1,
                      index[i] + 1,
                      traj ? traj : "the trajectory",
                      natoms);
        }
        else if (index[i] < 0)
        {
            gmx_fatal(FARGS,
                      "%s atom number (index[%d]=%d) is less than zero",
                      gname ? gname : "Index",
                      i + 1,
                      index[i] + 1);
        }
    }
}

t_blocka* init_index(const char* gfile, char*** grpname)
{
    t_blocka* b = nullptr;
    char      line[STRLEN], str[STRLEN];

    FILE* in = gmx_ffopen(gfile, "r");
    snew(b, 1);
    b->nr          = 0;
    b->index       = nullptr;
    b->nra         = 0;
    b->a           = nullptr;
    *grpname       = nullptr;
    int maxentries = 0;
    while (get_a_line(in, line, STRLEN))
    {
        if (get_header(line, str))
        {
            b->nr++;
            srenew(b->index, b->nr + 1);
            srenew(*grpname, b->nr);
            if (b->nr == 1)
            {
                b->index[0] = 0;
            }
            b->index[b->nr]       = b->index[b->nr - 1];
            (*grpname)[b->nr - 1] = gmx_strdup(str);
        }
        else
        {
            if (b->nr == 0)
            {
                gmx_fatal(FARGS, "The first header of your indexfile is invalid");
            }
            char* pt = line;
            while (sscanf(pt, "%s", str) == 1)
            {
                int i = b->index[b->nr];
                if (i >= maxentries)
                {
                    maxentries += 1024;
                    srenew(b->a, maxentries);
                }
                assert(b->a != nullptr); // for clang analyzer
                b->a[i] = strtol(str, nullptr, 10) - 1;
                b->index[b->nr]++;
                (b->nra)++;
                pt = strstr(pt, str) + strlen(str);
            }
        }
    }
    gmx_ffclose(in);

    for (int i = 0; (i < b->nr); i++)
    {
        assert(b->a != nullptr); // for clang analyzer
        for (int j = b->index[i]; (j < b->index[i + 1]); j++)
        {
            if (b->a[j] < 0)
            {
                fprintf(stderr, "\nWARNING: negative index %d in group %s\n\n", b->a[j], (*grpname)[i]);
            }
        }
    }

    return b;
}

static void minstring(char* str)
{
    for (int i = 0; (i < static_cast<int>(strlen(str))); i++)
    {
        if (str[i] == '-')
        {
            str[i] = '_';
        }
    }
}

int find_group(const char* s, int ngrps, char** grpname)
{
    char      string[STRLEN];
    bool      bMultiple = false;
    const int n         = strlen(s);
    int       aa        = -1;
    /* first look for whole name match */
    {
        for (int i = 0; i < ngrps; i++)
        {
            if (gmx_strcasecmp_min(s, grpname[i]) == 0)
            {
                if (aa != -1)
                {
                    bMultiple = true;
                }
                aa = i;
            }
        }
    }
    /* second look for first string match */
    if (aa == -1)
    {
        for (int i = 0; i < ngrps; i++)
        {
            if (gmx_strncasecmp_min(s, grpname[i], n) == 0)
            {
                if (aa != -1)
                {
                    bMultiple = TRUE;
                }
                aa = i;
            }
        }
    }
    /* last look for arbitrary substring match */
    if (aa == -1)
    {
        char key[STRLEN];
        strncpy(key, s, sizeof(key) - 1);
        key[STRLEN - 1] = '\0';
        upstring(key);
        minstring(key);
        for (int i = 0; i < ngrps; i++)
        {
            strncpy(string, grpname[i], STRLEN - 1);
            upstring(string);
            minstring(string);
            if (strstr(string, key) != nullptr)
            {
                if (aa != -1)
                {
                    bMultiple = TRUE;
                }
                aa = i;
            }
        }
    }
    if (bMultiple)
    {
        printf("Error: Multiple groups '%s' selected\n", s);
        aa = -1;
    }
    return aa;
}

static int qgroup(int* a, int ngrps, char** grpname)
{
    char  s[STRLEN];
    int   aa       = 0;
    bool  bInRange = false;
    char* end      = nullptr;

    do
    {
        fprintf(stderr, "Select a group: ");
        do
        {
            if (scanf("%s", s) != 1)
            {
                gmx_fatal(FARGS, "Cannot read from input");
            }
            trim(s); /* remove spaces */
        } while (strlen(s) == 0);
        aa = strtol(s, &end, 10);
        if (aa == 0 && end[0] != '\0') /* string entered */
        {
            aa = find_group(s, ngrps, grpname);
        }
        bInRange = (aa >= 0 && aa < ngrps);
        if (!bInRange)
        {
            printf("Error: No such group '%s'\n", s);
        }
    } while (!bInRange);
    printf("Selected %d: '%s'\n", aa, grpname[aa]);
    *a = aa;
    return aa;
}

static void
rd_groups(t_blocka* grps, char** grpname, char* gnames[], int ngrps, int isize[], int* index[], int grpnr[])
{
    if (grps->nr == 0)
    {
        gmx_fatal(FARGS, "Error: no groups in indexfile");
    }
    for (int i = 0; (i < grps->nr); i++)
    {
        fprintf(stderr, "Group %5d (%15s) has %5d elements\n", i, grpname[i], grps->index[i + 1] - grps->index[i]);
    }
    for (int i = 0; (i < ngrps); i++)
    {
        int gnr1 = 0;
        if (grps->nr > 1)
        {
            do
            {
                gnr1 = qgroup(&grpnr[i], grps->nr, grpname);
                if ((gnr1 < 0) || (gnr1 >= grps->nr))
                {
                    fprintf(stderr, "Select between %d and %d.\n", 0, grps->nr - 1);
                }
            } while ((gnr1 < 0) || (gnr1 >= grps->nr));
        }
        else
        {
            fprintf(stderr, "There is one group in the index\n");
            gnr1 = 0;
        }
        gnames[i] = gmx_strdup(grpname[gnr1]);
        isize[i]  = grps->index[gnr1 + 1] - grps->index[gnr1];
        snew(index[i], isize[i]);
        for (int j = 0; (j < isize[i]); j++)
        {
            index[i][j] = grps->a[grps->index[gnr1] + j];
        }
    }
}

void rd_index(const char* statfile, int ngrps, int isize[], int* index[], char* grpnames[])
{
    char** gnames = nullptr;
    int*   grpnr  = nullptr;

    snew(grpnr, ngrps);
    if (!statfile)
    {
        gmx_fatal(FARGS, "No index file specified");
    }
    t_blocka* grps = init_index(statfile, &gnames);
    rd_groups(grps, gnames, grpnames, ngrps, isize, index, grpnr);
    for (int i = 0; i < grps->nr; i++)
    {
        sfree(gnames[i]);
    }
    sfree(gnames);
    sfree(grpnr);
    done_blocka(grps);
    sfree(grps);
}

void get_index(const t_atoms* atoms, const char* fnm, int ngrps, int isize[], int* index[], char* grpnames[])
{
    char***   gnames = nullptr;
    t_blocka* grps   = nullptr;
    int*      grpnr  = nullptr;

    snew(grpnr, ngrps);
    snew(gnames, 1);
    if (fnm != nullptr)
    {
        grps = init_index(fnm, gnames);
    }
    else if (atoms)
    {
        snew(grps, 1);
        snew(grps->index, 1);
        analyse(atoms, grps, gnames, FALSE, FALSE);
    }
    else
    {
        gmx_incons("You need to supply a valid atoms structure or a valid index file name");
    }

    rd_groups(grps, *gnames, grpnames, ngrps, isize, index, grpnr);
    for (int i = 0; i < grps->nr; ++i)
    {
        sfree((*gnames)[i]);
    }
    sfree(*gnames);
    sfree(gnames);
    sfree(grpnr);
    done_blocka(grps);
    sfree(grps);
}

t_cluster_ndx* cluster_index(FILE* fplog, const char* ndx)
{
    t_cluster_ndx* c = nullptr;

    snew(c, 1);
    c->clust    = init_index(ndx, &c->grpname);
    c->maxframe = -1;
    for (int i = 0; (i < c->clust->nra); i++)
    {
        c->maxframe = std::max(c->maxframe, c->clust->a[i]);
    }
    fprintf(fplog ? fplog : stdout,
            "There are %d clusters containing %d structures, highest framenr is %d\n",
            c->clust->nr,
            c->clust->nra,
            c->maxframe);
    if (debug)
    {
        pr_blocka(debug, 0, "clust", c->clust, TRUE);
        for (int i = 0; (i < c->clust->nra); i++)
        {
            if ((c->clust->a[i] < 0) || (c->clust->a[i] > c->maxframe))
            {
                gmx_fatal(FARGS,
                          "Range check error for c->clust->a[%d] = %d\n"
                          "should be within 0 and %d",
                          i,
                          c->clust->a[i],
                          c->maxframe + 1);
            }
        }
    }
    c->inv_clust = make_invblocka(c->clust, c->maxframe);

    return c;
}

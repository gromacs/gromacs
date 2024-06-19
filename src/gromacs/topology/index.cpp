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

#include "gromacs/topology/index.h"

#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/invblock.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"
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

void write_index(const char* outf, gmx::ArrayRef<const IndexGroup> indexGroups, const bool duplicate, const int numAtoms)
{
    FILE* out = gmx_ffopen(outf, "w");

    for (const auto& indexGroup : indexGroups)
    {
        fprintf(out, "[ %s ]", indexGroup.name.c_str());
        int k = 0;
        for (const int particleIndex : indexGroup.particleIndices)
        {
            const char sep = (k % 15 == 0 ? '\n' : ' ');
            fprintf(out, "%c%4d", sep, particleIndex + 1);
            k++;
        }
        fprintf(out, "\n");
    }

    /* Duplicate copy, useful for computational electrophysiology double-layer setups */
    if (duplicate)
    {
        fprintf(stderr, "Duplicating the whole system with an atom offset of %d atoms.\n", numAtoms);
        for (const auto& indexGroup : indexGroups)
        {
            fprintf(out, "[ %s_copy ]", indexGroup.name.c_str());
            int k = 0;
            for (const int particleIndex : indexGroup.particleIndices)
            {
                const char sep = (k % 15 == 0 ? '\n' : ' ');
                fprintf(out, "%c%4d", sep, particleIndex + 1 + numAtoms);
                k++;
            }
            fprintf(out, "\n");
        }
    }

    gmx_ffclose(out);
}

/*! \brief
 * Compare index in groups.
 *
 * Checks the index in \p a to the one in \p b at \p index.
 * If \p index is < 0, it is taken as relative to the end of \p b.
 *
 * \param[in] indexGroups  List of index groups
 * \param[in] a New group to compare to.
 * \param[in] index The index to check, when negative this is the offset from the last index
 * \returns True if groups are the same.
 */
static bool grp_cmp(gmx::ArrayRef<const IndexGroup> indexGroups, gmx::ArrayRef<const int> a, int index)
{
    if (index >= indexGroups.ssize())
    {
        gmx_fatal(FARGS, "no such index group %d in index groups (nr=%td)", index, indexGroups.ssize());
    }

    const int indexToCompare = (index >= 0 ? index : indexGroups.ssize() - 1 + index);
    GMX_RELEASE_ASSERT(indexToCompare >= 0 && indexToCompare < indexGroups.ssize(),
                       "Index should be in range");
    gmx::ArrayRef<const int> indexGroup = indexGroups[indexToCompare].particleIndices;

    /* compare sizes */
    if (a.ssize() != indexGroup.ssize())
    {
        return false;
    }
    for (gmx::Index i = 0; i < a.ssize(); i++)
    {
        if (a[i] != indexGroup[i])
        {
            return false;
        }
    }
    return true;
}
//! Print out how many residues of a certain molecule category are present.
static void printMoleculeCategoryResidueCounts(gmx::ArrayRef<const std::pair<std::string, int>> moleculeCategoryCount)
{
    for (const auto& mcc : moleculeCategoryCount)
    {
        if (mcc.second > 0)
        {
            printf("There are: %5d %10s residues\n", mcc.second, mcc.first.c_str());
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
                          std::vector<IndexGroup>*   indexGroups,
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
            indexGroups->push_back({ restp[i].gname, aid });
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
                            indexGroups->push_back({ attp[l], aaid });
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
                         std::vector<IndexGroup>*         indexGroups,
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
        if (-1 == constructing_data[i].compareto
            || !grp_cmp(*indexGroups, aid, constructing_data[i].compareto - i))
        {
            indexGroups->push_back({ constructing_data[i].group_name, aid });
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
                for (int n = 0; ((n < atoms->nr) && (atoms->atom[n].resind < npres));)
                {
                    int resind = atoms->atom[n].resind;
                    for (; ((n < atoms->nr) && (atoms->atom[n].resind == resind)); n++)
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
                        indexGroups->push_back({ ndx_name, aid });
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
            for (int n = 0; ((n < atoms->nr) && (atoms->atom[n].resind < npres));)
            {
                int resind = atoms->atom[n].resind;
                int hold   = -1;
                for (; ((n < atoms->nr) && (atoms->atom[n].resind == resind)); n++)
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
                indexGroups->push_back({ "SwapSC-CO", aid });
            }
        }
    }
}


std::vector<IndexGroup> analyse(const t_atoms* atoms, gmx_bool bASK, gmx_bool bVerb)
{
    if (bVerb)
    {
        printf("Analysing residue names:\n");
    }

    std::vector<IndexGroup> indexGroups;

    /* Create system group, every single atom */
    std::vector<int> aid(atoms->nr);
    std::iota(aid.begin(), aid.end(), 0);
    indexGroups.push_back({ "System", aid });

    /* For every residue, get a pointer to the residue type name */
    ResidueTypeMap residueTypeMap = residueTypeMapFromLibraryFile("residuetypes.dat");

    std::vector<std::string>                 moleculeCategoryOfResidues;
    std::vector<std::pair<std::string, int>> moleculeCategoryCount;
    if (atoms->nres > 0)
    {
        for (int i = 0; i < atoms->nres; i++)
        {
            const char* resnm = *atoms->resinfo[i].name;
            moleculeCategoryOfResidues.emplace_back(typeOfNamedDatabaseResidue(residueTypeMap, resnm));

            /* Note that this does not lead to a N*N loop, but N*K, where
             * K is the number of residue _types_, which is small and independent of N.
             */
            bool found = false;
            for (auto& mcc : moleculeCategoryCount)
            {
                found = (moleculeCategoryOfResidues.back() == mcc.first);
                if (found)
                {
                    mcc.second++;
                    break;
                }
            }
            if (!found)
            {
                moleculeCategoryCount.emplace_back(moleculeCategoryOfResidues.back(), 1);
            }
        }
    }

    if (bVerb)
    {
        printMoleculeCategoryResidueCounts(moleculeCategoryCount);
    }

    bool haveAnalysedOther = false;
    for (const auto& mcc : moleculeCategoryCount)
    {
        const std::string& moleculeCategory = mcc.first;

        std::vector<int> aid = mk_aid(atoms, moleculeCategoryOfResidues, moleculeCategory, TRUE);

        /* Check for special types to do fancy stuff with */

        if (!gmx_strcasecmp(moleculeCategory.c_str(), "Protein") && !aid.empty())
        {
            /* PROTEIN */
            analyse_prot(moleculeCategoryOfResidues, atoms, &indexGroups, bASK, bVerb);

            /* Create a Non-Protein group */
            std::vector<int> aid = mk_aid(atoms, moleculeCategoryOfResidues, "Protein", FALSE);
            if ((!aid.empty()) && (gmx::ssize(aid) < atoms->nr))
            {
                indexGroups.push_back({ "non-Protein", aid });
            }
        }
        else if (!gmx_strcasecmp(moleculeCategory.c_str(), "Water") && !aid.empty())
        {
            indexGroups.push_back({ moleculeCategory, aid });
            /* Add this group as 'SOL' too, for backward compatibility with older gromacs versions */
            indexGroups.push_back({ "SOL", aid });


            /* Solvent, create a negated group too */
            std::vector<int> aid = mk_aid(atoms, moleculeCategoryOfResidues, "Water", FALSE);
            if ((!aid.empty()) && (gmx::ssize(aid) < atoms->nr))
            {
                indexGroups.push_back({ "non-Water", aid });
            }
        }
        else if (!gmx_strcasecmp(moleculeCategory.c_str(), "Ion") && !aid.empty())
        {
            indexGroups.push_back({ moleculeCategory, aid });
        }
        else if (!aid.empty() && !haveAnalysedOther)
        {
            /* Other groups */
            indexGroups.push_back({ moleculeCategory, aid });
            analyse_other(moleculeCategoryOfResidues, atoms, &indexGroups, bASK, bVerb);
            haveAnalysedOther = true;
        }
    }


    /* Create a merged water_and_ions group */
    int iwater = -1;
    int iion   = -1;
    int nwater = 0;
    int nion   = 0;

    for (gmx::Index i = 0; i < gmx::ssize(indexGroups); i++)
    {
        if (!gmx_strcasecmp(indexGroups[i].name.c_str(), "Water"))
        {
            iwater = i;
            nwater = indexGroups[i].particleIndices.size();
        }
        else if (!gmx_strcasecmp(indexGroups[i].name.c_str(), "Ion"))
        {
            iion = i;
            nion = indexGroups[i].particleIndices.size();
        }
    }

    if (nwater > 0 && nion > 0)
    {
        indexGroups.push_back({ "Water_and_ions", {} });
        std::vector<int>& a = indexGroups.back().particleIndices;
        a.insert(a.end(),
                 indexGroups[iwater].particleIndices.begin(),
                 indexGroups[iwater].particleIndices.end());
        a.insert(a.end(),
                 indexGroups[iion].particleIndices.begin(),
                 indexGroups[iion].particleIndices.end());
    }

    return indexGroups;
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

std::vector<IndexGroup> init_index(const std::filesystem::path& gfile)
{
    return init_index(gfile.string().c_str());
}

std::vector<IndexGroup> init_index(const char* gfile)
{
    std::vector<IndexGroup> indexGroups;

    char line[STRLEN], str[STRLEN];

    FILE* in = gmx_ffopen(gfile, "r");

    std::vector<int>* atomListPtr = nullptr;
    while (get_a_line(in, line, STRLEN))
    {
        if (get_header(line, str))
        {
            indexGroups.push_back({ str, {} });
            atomListPtr = &indexGroups.back().particleIndices;
        }
        else
        {
            if (indexGroups.empty())
            {
                gmx_fatal(FARGS, "The first header of your indexfile is invalid");
            }
            GMX_RELEASE_ASSERT(atomListPtr != nullptr,
                               "Here we should have a valid atom list pointer");
            char* pt = line;
            while (sscanf(pt, "%s", str) == 1)
            {
                atomListPtr->push_back(strtol(str, nullptr, 10) - 1);
                pt = strstr(pt, str) + strlen(str);
            }
        }
    }
    gmx_ffclose(in);

    for (const auto& indexGroup : indexGroups)
    {
        for (const int particleIndex : indexGroup.particleIndices)
        {
            if (particleIndex < 0)
            {
                fprintf(stderr,
                        "\nWARNING: negative index %d in group %s\n\n",
                        particleIndex,
                        indexGroup.name.c_str());
            }
        }
    }

    return indexGroups;
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

static const char* indexGroupName(const IndexGroup& indexGroup)
{
    return indexGroup.name.c_str();
}

static const char* indexGroupName(const char* indexGroup)
{
    return indexGroup;
}


template<typename T>
static int findGroupTemplated(const char* s, gmx::ArrayRef<T> indexGroups)
{
    char      string[STRLEN];
    bool      bMultiple = false;
    const int n         = strlen(s);
    int       aa        = -1;
    /* first look for whole name match */
    {
        for (gmx::Index i = 0; i < indexGroups.ssize(); i++)
        {
            if (gmx_strcasecmp_min(s, indexGroupName(indexGroups[i])) == 0)
            {
                if (aa != -1)
                {
                    bMultiple = true;
                }
                aa = i;
            }
        }
    }
    /* particleIndices look for first string match */
    if (aa == -1)
    {
        for (gmx::Index i = 0; i < indexGroups.ssize(); i++)
        {
            if (gmx_strncasecmp_min(s, indexGroupName(indexGroups[i]), n) == 0)
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
        for (gmx::Index i = 0; i < indexGroups.ssize(); i++)
        {
            strncpy(string, indexGroupName(indexGroups[i]), STRLEN - 1);
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

int find_group(const char* s, gmx::ArrayRef<const IndexGroup> indexGroups)
{
    return findGroupTemplated(s, indexGroups);
}

int find_group(const char* s, int ngrps, char** grpname)
{
    return findGroupTemplated(s, gmx::constArrayRefFromArray<const char*>(grpname, ngrps));
}

static int qgroup(gmx::ArrayRef<const IndexGroup> indexGroups)
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
            aa = find_group(s, indexGroups);
        }
        bInRange = (aa >= 0 && aa < indexGroups.ssize());
        if (!bInRange)
        {
            printf("Error: No such group '%s'\n", s);
        }
    } while (!bInRange);
    printf("Selected %d: '%s'\n", aa, indexGroups[aa].name.c_str());

    return aa;
}

static void rd_groups(gmx::ArrayRef<const IndexGroup> indexGroups,
                      char*                           gnames[],
                      int                             ngrps,
                      int                             isize[],
                      int*                            index[])
{
    if (indexGroups.empty())
    {
        gmx_fatal(FARGS, "Error: no groups in indexfile");
    }
    for (gmx::Index i = 0; i < indexGroups.ssize(); i++)
    {
        fprintf(stderr,
                "Group %5zd (%15s) has %5zd elements\n",
                i,
                indexGroups[i].name.c_str(),
                gmx::ssize(indexGroups[i].particleIndices));
    }
    for (int i = 0; (i < ngrps); i++)
    {
        int gnr1 = 0;
        if (indexGroups.size() > 1)
        {
            do
            {
                gnr1 = qgroup(indexGroups);
                if (gnr1 < 0 || gnr1 >= indexGroups.ssize())
                {
                    fprintf(stderr, "Select between %d and %zd.\n", 0, indexGroups.ssize() - 1);
                }
            } while (gnr1 < 0 || gnr1 >= indexGroups.ssize());
        }
        else
        {
            fprintf(stderr, "There is one group in the index\n");
            gnr1 = 0;
        }
        gnames[i] = gmx_strdup(indexGroups[gnr1].name.c_str());
        isize[i]  = gmx::ssize(indexGroups[gnr1].particleIndices);
        snew(index[i], isize[i]);
        for (int j = 0; (j < isize[i]); j++)
        {
            index[i][j] = indexGroups[gnr1].particleIndices[j];
        }
    }
}

void rd_index(const std::filesystem::path& statfile, int ngrps, int isize[], int* index[], char* grpnames[])
{
    rd_index(statfile.string().c_str(), ngrps, isize, index, grpnames);
}

void rd_index(const char* statfile, int ngrps, int isize[], int* index[], char* grpnames[])
{
    if (!statfile)
    {
        gmx_fatal(FARGS, "No index file specified");
    }
    const auto indexGroups = init_index(statfile);
    rd_groups(indexGroups, grpnames, ngrps, isize, index);
}

void get_index(const t_atoms*                              atoms,
               const std::optional<std::filesystem::path>& fnm,
               int                                         ngrps,
               int                                         isize[],
               int*                                        index[],
               char*                                       grpnames[])
{
    if (fnm)
    {
        get_index(atoms, fnm.value().string().c_str(), ngrps, isize, index, grpnames);
    }
    else
    {
        get_index(atoms, nullptr, ngrps, isize, index, grpnames);
    }
}

// Deprecated, remove and consolidate with the function above when there are no other callers
void get_index(const t_atoms* atoms, const char* fnm, int ngrps, int isize[], int* index[], char* grpnames[])
{
    std::vector<IndexGroup> indexGroups;
    if (fnm != nullptr)
    {
        indexGroups = init_index(fnm);
    }
    else if (atoms)
    {
        indexGroups = analyse(atoms, FALSE, FALSE);
    }
    else
    {
        gmx_incons("You need to supply a valid atoms structure or a valid index file name");
    }

    rd_groups(indexGroups, grpnames, ngrps, isize, index);
}

t_cluster_ndx cluster_index(FILE* fplog, const char* ndx)
{
    t_cluster_ndx c;

    c.clusters                = init_index(ndx);
    c.maxframe                = -1;
    gmx::Index totalNumFrames = 0;
    for (const auto& cluster : c.clusters)
    {
        for (const int frame : cluster.particleIndices)
        {
            c.maxframe = std::max(c.maxframe, frame);
        }
        totalNumFrames += gmx::ssize(cluster.particleIndices);
    }
    fprintf(fplog ? fplog : stdout,
            "There are %td clusters containing %td structures, highest framenr is %d\n",
            gmx::ssize(c.clusters),
            totalNumFrames,
            c.maxframe);

    if (debug)
    {
        pr_blocka(debug, 0, "clust", c.clusters, true);

        for (const auto& cluster : c.clusters)
        {
            int i = 0;
            for (const int frame : cluster.particleIndices)
            {
                if (frame < 0 || frame > c.maxframe)
                {
                    gmx_fatal(FARGS,
                              "Range check error for c.clust->a[%d] = %d\n"
                              "should be within 0 and %d",
                              i,
                              frame,
                              c.maxframe + 1);
                }
                i++;
            }
        }
    }

    gmx::ListOfLists<int> clusterList;
    for (const IndexGroup& cluster : c.clusters)
    {
        clusterList.pushBack(cluster.particleIndices);
    }
    c.inv_clust = make_invblock(clusterList, c.maxframe);

    return c;
}

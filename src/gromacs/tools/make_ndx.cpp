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

#include "make_ndx.h"

#include <cctype>
#include <cstring>

#include <algorithm>
#include <array>
#include <string>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

/* It's not nice to have size limits, but we should not spend more time
 * on this ancient tool, but instead use the new selection library.
 */
#define MAXNAMES 1024
#define NAME_LEN 1024
static const int NOTSET = -92637;

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static gmx_bool bCase = FALSE;

static int or_groups(int nr1, const int* at1, int nr2, const int* at2, int* nr, int* at)
{
    int      i1, i2, max = 0;
    gmx_bool bNotIncr;

    *nr = 0;

    bNotIncr = FALSE;
    for (i1 = 0; i1 < nr1; i1++)
    {
        if ((i1 > 0) && (at1[i1] <= max))
        {
            bNotIncr = TRUE;
        }
        max = at1[i1];
    }
    for (i1 = 0; i1 < nr2; i1++)
    {
        if ((i1 > 0) && (at2[i1] <= max))
        {
            bNotIncr = TRUE;
        }
        max = at2[i1];
    }

    if (bNotIncr)
    {
        printf("One of your groups is not ascending\n");
    }
    else
    {
        i1  = 0;
        i2  = 0;
        *nr = 0;
        while ((i1 < nr1) || (i2 < nr2))
        {
            if ((i2 == nr2) || ((i1 < nr1) && (at1[i1] < at2[i2])))
            {
                at[*nr] = at1[i1];
                (*nr)++;
                i1++;
            }
            else
            {
                if ((i2 < nr2) && ((i1 == nr1) || (at1[i1] > at2[i2])))
                {
                    at[*nr] = at2[i2];
                    (*nr)++;
                }
                i2++;
            }
        }

        printf("Merged two groups with OR: %d %d -> %d\n", nr1, nr2, *nr);
    }

    return *nr;
}

static int and_groups(int nr1, const int* at1, int nr2, const int* at2, int* nr, int* at)
{
    int i1, i2;

    *nr = 0;
    for (i1 = 0; i1 < nr1; i1++)
    {
        for (i2 = 0; i2 < nr2; i2++)
        {
            if (at1[i1] == at2[i2])
            {
                at[*nr] = at1[i1];
                (*nr)++;
            }
        }
    }

    printf("Merged two groups with AND: %d %d -> %d\n", nr1, nr2, *nr);

    return *nr;
}

static gmx_bool is_name_char(char c)
{
    /* This string should contain all characters that can not be
     * the first letter of a name due to the make_ndx syntax.
     */
    const char* spec = " !&|";

    return (c != '\0' && std::strchr(spec, c) == nullptr);
}

static int parse_names(char** string, int* n_names, gmx::ArrayRef<char*> names)
{
    int i;

    *n_names = 0;
    while (is_name_char((*string)[0]) || (*string)[0] == ' ')
    {
        if (is_name_char((*string)[0]))
        {
            if (*n_names >= MAXNAMES)
            {
                gmx_fatal(FARGS, "To many names: %d\n", *n_names + 1);
            }
            i = 0;
            while (is_name_char((*string)[i]))
            {
                names[*n_names][i] = (*string)[i];
                i++;
                if (i > NAME_LEN)
                {
                    printf("Name is too long, the maximum is %d characters\n", NAME_LEN);
                    return 0;
                }
            }
            names[*n_names][i] = '\0';
            if (!bCase)
            {
                upstring(names[*n_names]);
            }
            *string += i;
            (*n_names)++;
        }
        else
        {
            (*string)++;
        }
    }

    return *n_names;
}

static gmx_bool parse_int_char(char** string, int* nr, unsigned char* c)
{
    char*    orig;
    gmx_bool bRet;

    orig = *string;

    while ((*string)[0] == ' ')
    {
        (*string)++;
    }

    bRet = FALSE;

    *c = ' ';

    if (std::isdigit((*string)[0]))
    {
        *nr = (*string)[0] - '0';
        (*string)++;
        while (std::isdigit((*string)[0]))
        {
            *nr = (*nr) * 10 + (*string)[0] - '0';
            (*string)++;
        }
        if (std::isalpha((*string)[0]))
        {
            *c = (*string)[0];
            (*string)++;
        }
        /* Check if there is at most one non-digit character */
        if (!std::isalnum((*string)[0]))
        {
            bRet = TRUE;
        }
        else
        {
            *string = orig;
        }
    }
    else
    {
        *nr = NOTSET;
    }

    return bRet;
}

static gmx_bool parse_int(char** string, int* nr)
{
    char*         orig;
    unsigned char c;
    gmx_bool      bRet;

    orig = *string;
    bRet = parse_int_char(string, nr, &c);
    if (bRet && c != ' ')
    {
        *string = orig;
        bRet    = FALSE;
    }

    return bRet;
}

static gmx_bool isquote(char c)
{
    return (c == '\"');
}

static gmx_bool parse_string(char** string, int* nr, gmx::ArrayRef<const IndexGroup> indexGroups)
{
    char *s, *sp;
    char  c;

    while ((*string)[0] == ' ')
    {
        (*string)++;
    }

    (*nr) = NOTSET;
    if (isquote((*string)[0]))
    {
        c = (*string)[0];
        (*string)++;
        s  = gmx_strdup((*string));
        sp = std::strchr(s, c);
        if (sp != nullptr)
        {
            (*string) += sp - s + 1;
            sp[0] = '\0';
            (*nr) = find_group(s, indexGroups);
        }
    }

    return (*nr) != NOTSET;
}

static int select_atomnumbers(char** string, const t_atoms* atoms, int n1, int* nr, int* index, char* gname)
{
    char buf[STRLEN];
    int  i, up;

    *nr = 0;
    while ((*string)[0] == ' ')
    {
        (*string)++;
    }
    if ((*string)[0] == '-')
    {
        (*string)++;
        parse_int(string, &up);
        if ((n1 < 1) || (n1 > atoms->nr) || (up < 1) || (up > atoms->nr))
        {
            printf("Invalid atom range\n");
        }
        else
        {
            for (i = n1 - 1; i <= up - 1; i++)
            {
                index[*nr] = i;
                (*nr)++;
            }
            printf("Found %d atom%s in range %d-%d\n", *nr, (*nr == 1) ? "" : "s", n1, up);
            if (n1 == up)
            {
                sprintf(buf, "a_%d", n1);
            }
            else
            {
                sprintf(buf, "a_%d-%d", n1, up);
            }
            std::strcpy(gname, buf);
        }
    }
    else
    {
        i = n1;
        sprintf(gname, "a");
        do
        {
            if ((i - 1 >= 0) && (i - 1 < atoms->nr))
            {
                index[*nr] = i - 1;
                (*nr)++;
                sprintf(buf, "_%d", i);
                std::strcat(gname, buf);
            }
            else
            {
                printf("Invalid atom number %d\n", i);
                *nr = 0;
            }
        } while ((*nr != 0) && (parse_int(string, &i)));
    }

    return *nr;
}

static int
select_residuenumbers(char** string, const t_atoms* atoms, int n1, unsigned char c, int* nr, int* index, char* gname)
{
    char       buf[STRLEN];
    int        i, j, up;
    t_resinfo* ri;

    *nr = 0;
    while ((*string)[0] == ' ')
    {
        (*string)++;
    }
    if ((*string)[0] == '-')
    {
        /* Residue number range selection */
        if (c != ' ')
        {
            printf("Error: residue insertion codes can not be used with residue range selection\n");
            return 0;
        }
        (*string)++;
        parse_int(string, &up);

        for (i = 0; i < atoms->nr; i++)
        {
            ri = &atoms->resinfo[atoms->atom[i].resind];
            for (j = n1; (j <= up); j++)
            {
                if (ri->nr == j && (c == ' ' || ri->ic == c))
                {
                    index[*nr] = i;
                    (*nr)++;
                }
            }
        }
        printf("Found %d atom%s with res.nr. in range %d-%d\n", *nr, (*nr == 1) ? "" : "s", n1, up);
        if (n1 == up)
        {
            sprintf(buf, "r_%d", n1);
        }
        else
        {
            sprintf(buf, "r_%d-%d", n1, up);
        }
        std::strcpy(gname, buf);
    }
    else
    {
        /* Individual residue number/insertion code selection */
        j = n1;
        sprintf(gname, "r");
        do
        {
            for (i = 0; i < atoms->nr; i++)
            {
                ri = &atoms->resinfo[atoms->atom[i].resind];
                if (ri->nr == j && ri->ic == c)
                {
                    index[*nr] = i;
                    (*nr)++;
                }
            }
            sprintf(buf, "_%d", j);
            std::strcat(gname, buf);
        } while (parse_int_char(string, &j, &c));
    }

    return *nr;
}

static int
select_residueindices(char** string, const t_atoms* atoms, int n1, unsigned char c, int* nr, int* index, char* gname)
{
    /*this should be similar to select_residuenumbers except select by index (sequential numbering in file)*/
    /*resind+1 for 1-indexing*/
    char       buf[STRLEN];
    int        i, j, up;
    t_resinfo* ri;

    *nr = 0;
    while ((*string)[0] == ' ')
    {
        (*string)++;
    }
    if ((*string)[0] == '-')
    {
        /* Residue number range selection */
        if (c != ' ')
        {
            printf("Error: residue insertion codes can not be used with residue range selection\n");
            return 0;
        }
        (*string)++;
        parse_int(string, &up);

        for (i = 0; i < atoms->nr; i++)
        {
            ri = &atoms->resinfo[atoms->atom[i].resind];
            for (j = n1; (j <= up); j++)
            {
                if (atoms->atom[i].resind + 1 == j && (c == ' ' || ri->ic == c))
                {
                    index[*nr] = i;
                    (*nr)++;
                }
            }
        }
        printf("Found %d atom%s with resind.+1 in range %d-%d\n", *nr, (*nr == 1) ? "" : "s", n1, up);
        if (n1 == up)
        {
            sprintf(buf, "r_%d", n1);
        }
        else
        {
            sprintf(buf, "r_%d-%d", n1, up);
        }
        std::strcpy(gname, buf);
    }
    else
    {
        /* Individual residue number/insertion code selection */
        j = n1;
        sprintf(gname, "r");
        do
        {
            for (i = 0; i < atoms->nr; i++)
            {
                ri = &atoms->resinfo[atoms->atom[i].resind];
                if (atoms->atom[i].resind + 1 == j && ri->ic == c)
                {
                    index[*nr] = i;
                    (*nr)++;
                }
            }
            sprintf(buf, "_%d", j);
            std::strcat(gname, buf);
        } while (parse_int_char(string, &j, &c));
    }

    return *nr;
}


static gmx_bool atoms_from_residuenumbers(const t_atoms* atoms, const IndexGroup& indexGroup, int* nr, int* index)
{
    const int nres = atoms->nres;
    for (const int resnr : indexGroup.particleIndices)
    {
        if (resnr >= nres)
        {
            printf("Index %s contains number>nres (%d>%d)\n", indexGroup.name.c_str(), resnr + 1, nres);
            return FALSE;
        }
    }
    for (int i = 0; i < atoms->nr; i++)
    {
        const int selectedResnr = atoms->resinfo[atoms->atom[i].resind].nr;
        for (const int resnr : indexGroup.particleIndices)
        {
            if (resnr + 1 == selectedResnr)
            {
                index[*nr] = i;
                (*nr)++;
                break;
            }
        }
    }
    printf("Found %d atom%s in %td residues from group %s\n",
           *nr,
           (*nr == 1) ? "" : "s",
           gmx::ssize(indexGroup.particleIndices),
           indexGroup.name.c_str());
    return *nr != 0;
}

static gmx_bool comp_name(const char* name, const char* search)
{
    gmx_bool matches = TRUE;

    // Loop while name and search are not end-of-string and matches is true
    for (; *name && *search && matches; name++, search++)
    {
        if (*search == '?')
        {
            // still matching, continue to next character
            continue;
        }
        else if (*search == '*')
        {
            if (*(search + 1))
            {
                printf("WARNING: Currently '*' is only supported at the end of an expression\n");
            }
            // if * is the last char in search string, we have a match,
            // otherwise we just failed. Return in either case, we're done.
            return (*(search + 1) == '\0');
        }

        // Compare one character
        if (bCase)
        {
            matches = (*name == *search);
        }
        else
        {
            matches = (std::toupper(*name) == std::toupper(*search));
        }
    }

    matches = matches && (*name == '\0' && (*search == '\0' || *search == '*'));

    return matches;
}

static int select_chainnames(const t_atoms* atoms, int n_names, gmx::ArrayRef<char*> names, int* nr, int* index)
{
    char name[2];
    int  j;
    int  i;

    name[1] = 0;
    *nr     = 0;
    for (i = 0; i < atoms->nr; i++)
    {
        name[0] = atoms->resinfo[atoms->atom[i].resind].chainid;
        j       = 0;
        while (j < n_names && !comp_name(name, names[j]))
        {
            j++;
        }
        if (j < n_names)
        {
            index[*nr] = i;
            (*nr)++;
        }
    }
    printf("Found %d atom%s with chain identifier%s", *nr, (*nr == 1) ? "" : "s", (n_names == 1) ? "" : "s");
    for (j = 0; (j < n_names); j++)
    {
        printf(" %s", names[j]);
    }
    printf("\n");

    return *nr;
}

static int select_atomnames(const t_atoms* atoms, int n_names, gmx::ArrayRef<char*> names, int* nr, int* index, gmx_bool bType)
{
    char* name;
    int   j;
    int   i;

    *nr = 0;
    for (i = 0; i < atoms->nr; i++)
    {
        if (bType)
        {
            name = *(atoms->atomtype[i]);
        }
        else
        {
            name = *(atoms->atomname[i]);
        }
        j = 0;
        while (j < n_names && !comp_name(name, names[j]))
        {
            j++;
        }
        if (j < n_names)
        {
            index[*nr] = i;
            (*nr)++;
        }
    }
    printf("Found %d atoms with %s%s", *nr, bType ? "type" : "name", (n_names == 1) ? "" : "s");
    for (j = 0; (j < n_names); j++)
    {
        printf(" %s", names[j]);
    }
    printf("\n");

    return *nr;
}

static int select_residuenames(const t_atoms* atoms, int n_names, gmx::ArrayRef<char*> names, int* nr, int* index)
{
    char* name;
    int   j;
    int   i;

    *nr = 0;
    for (i = 0; i < atoms->nr; i++)
    {
        name = *(atoms->resinfo[atoms->atom[i].resind].name);
        j    = 0;
        while (j < n_names && !comp_name(name, names[j]))
        {
            j++;
        }
        if (j < n_names)
        {
            index[*nr] = i;
            (*nr)++;
        }
    }
    printf("Found %d atoms with residue name%s", *nr, (n_names == 1) ? "" : "s");
    for (j = 0; (j < n_names); j++)
    {
        printf(" %s", names[j]);
    }
    printf("\n");

    return *nr;
}

static void make_gname(int n, gmx::ArrayRef<char*> names, char* gname)
{
    int i;

    std::strcpy(gname, names[0]);
    for (i = 1; i < n; i++)
    {
        std::strcat(gname, "_");
        std::strcat(gname, names[i]);
    }
}

static void copy_group(gmx::ArrayRef<const int> group, int* nr, int* index)
{
    *nr = gmx::ssize(group);
    for (gmx::Index i = 0; i < *nr; i++)
    {
        index[i] = group[i];
    }
}

/* Remove groups in range [first, last] */
static void remove_group(const int first, const int last, std::vector<IndexGroup>* indexGroups)
{
    for (int j = 0; j <= last - first; j++)
    {
        if (first < 0 || first >= gmx::ssize(*indexGroups))
        {
            printf("Group %d does not exist\n", first + j);
        }
        else
        {
            printf("Removed group %d '%s'\n", first + j, (*indexGroups)[first].name.c_str());
            indexGroups->erase(indexGroups->begin() + first);
        }
    }
}

static void split_group(const t_atoms* atoms, int sel_nr, std::vector<IndexGroup>* indexGroups, gmx_bool bAtom)
{
    char buf[STRLEN];

    const std::string nameToSplit =
            (*indexGroups)[sel_nr].name; // Need a copy since indexGroups can reallocate
    printf("Splitting group %d '%s' into %s\n", sel_nr, nameToSplit.c_str(), bAtom ? "atoms" : "residues");

    gmx::ArrayRef<const int> groupToSplit = (*indexGroups)[sel_nr].particleIndices;
    int                      prevAtom     = -1;
    for (const int a : groupToSplit)
    {
        const int   resind = atoms->atom[a].resind;
        const char* name   = *(atoms->resinfo[resind].name);
        if (bAtom || (prevAtom == -1) || (atoms->atom[prevAtom].resind != resind))
        {
            if (bAtom)
            {
                sprintf(buf, "%s_%s_%d", nameToSplit.c_str(), *atoms->atomname[a], a + 1);
            }
            else
            {
                sprintf(buf, "%s_%s_%d", nameToSplit.c_str(), name, atoms->resinfo[resind].nr);
            }
            indexGroups->push_back({ buf, {} });
        }
        GMX_ASSERT(!indexGroups->empty(), "Should have already created an index group!");
        indexGroups->back().particleIndices.push_back(a);
        prevAtom = a;
    }
}

static int split_chain(const t_atoms* atoms, const rvec* x, int sel_nr, std::vector<IndexGroup>* indexGroups)
{
    char buf[STRLEN];
    int  j, nchain;
    int  i, natoms, *start = nullptr, *end = nullptr, ca_start, ca_end;
    rvec vec;

    natoms   = atoms->nr;
    nchain   = 0;
    ca_start = 0;

    while (ca_start < natoms)
    {
        while ((ca_start < natoms) && std::strcmp(*atoms->atomname[ca_start], "CA") != 0)
        {
            ca_start++;
        }
        if (ca_start < natoms)
        {
            srenew(start, nchain + 1);
            srenew(end, nchain + 1);
            start[nchain] = ca_start;
            while ((start[nchain] > 0)
                   && (atoms->atom[start[nchain] - 1].resind == atoms->atom[ca_start].resind))
            {
                start[nchain]--;
            }

            i = ca_start;
            do
            {
                ca_end = i;
                do
                {
                    i++;
                } while ((i < natoms) && std::strcmp(*atoms->atomname[i], "CA") != 0);
                if (i < natoms)
                {
                    rvec_sub(x[ca_end], x[i], vec);
                }
                else
                {
                    break;
                }
            } while (norm(vec) < 0.45);

            end[nchain] = ca_end;
            while ((end[nchain] + 1 < natoms)
                   && (atoms->atom[end[nchain] + 1].resind == atoms->atom[ca_end].resind))
            {
                end[nchain]++;
            }
            ca_start = end[nchain] + 1;
            nchain++;
        }
    }
    if (nchain == 1)
    {
        printf("Found 1 chain, will not split\n");
    }
    else
    {
        printf("Found %d chains\n", nchain);
    }
    for (j = 0; j < nchain; j++)
    {
        printf("%d:%6d atoms (%d to %d)\n", j + 1, end[j] - start[j] + 1, start[j] + 1, end[j] + 1);
    }

    if (nchain > 1)
    {
        for (j = 0; j < nchain; j++)
        {
            std::vector<int> particleIndices;
            for (const int a : (*indexGroups)[sel_nr].particleIndices)
            {
                if ((a >= start[j]) && (a <= end[j]))
                {
                    particleIndices.push_back(a);
                }
            }
            if (!particleIndices.empty())
            {
                sprintf(buf, "%s_chain%d", (*indexGroups)[sel_nr].name.c_str(), j + 1);
                indexGroups->push_back({ buf, particleIndices });
            }
        }
    }
    sfree(start);
    sfree(end);

    return nchain;
}

static gmx_bool check_have_atoms(const t_atoms* atoms, char* string)
{
    if (atoms == nullptr)
    {
        printf("Can not process '%s' without atom info, use option -f\n", string);
        return FALSE;
    }
    else
    {
        return TRUE;
    }
}

static gmx_bool parse_entry(char**                   string,
                            int                      natoms,
                            const t_atoms*           atoms,
                            std::vector<IndexGroup>* indexGroups,
                            int*                     nr,
                            int*                     index,
                            char*                    gname,
                            gmx::ArrayRef<char*>     entryNames)
{
    char*         ostring;
    int           j, n_names, sel_nr1;
    int           i, nr1, *index1;
    unsigned char c;
    gmx_bool      bRet, bCompl;

    bRet    = FALSE;
    sel_nr1 = NOTSET;

    while (*string[0] == ' ')
    {
        (*string)++;
    }

    if ((*string)[0] == '!')
    {
        bCompl = TRUE;
        (*string)++;
        while (*string[0] == ' ')
        {
            (*string)++;
        }
    }
    else
    {
        bCompl = FALSE;
    }

    ostring = *string;

    if (parse_int(string, &sel_nr1) || parse_string(string, &sel_nr1, *indexGroups))
    {
        if (sel_nr1 >= 0 && sel_nr1 < gmx::ssize(*indexGroups))
        {
            copy_group((*indexGroups)[sel_nr1].particleIndices, nr, index);
            std::strcpy(gname, (*indexGroups)[sel_nr1].name.c_str());
            printf("Copied index group %d '%s'\n", sel_nr1, (*indexGroups)[sel_nr1].name.c_str());
            bRet = TRUE;
        }
        else
        {
            printf("Group %d does not exist\n", sel_nr1);
        }
    }
    else if ((*string)[0] == 'a')
    {
        (*string)++;
        if (check_have_atoms(atoms, ostring))
        {
            if (parse_int(string, &sel_nr1))
            {
                bRet = (select_atomnumbers(string, atoms, sel_nr1, nr, index, gname) != 0);
            }
            else if (parse_names(string, &n_names, entryNames))
            {
                bRet = (select_atomnames(atoms, n_names, entryNames, nr, index, FALSE) != 0);
                make_gname(n_names, entryNames, gname);
            }
        }
    }
    else if ((*string)[0] == 't')
    {
        (*string)++;
        if (check_have_atoms(atoms, ostring) && parse_names(string, &n_names, entryNames))
        {
            if (!(atoms->haveType))
            {
                printf("Need a run input file to select atom types\n");
            }
            else
            {
                bRet = (select_atomnames(atoms, n_names, entryNames, nr, index, TRUE) != 0);
                make_gname(n_names, entryNames, gname);
            }
        }
    }
    else if (std::strncmp(*string, "res", 3) == 0)
    {
        (*string) += 3;
        if (check_have_atoms(atoms, ostring) && parse_int(string, &sel_nr1) && (sel_nr1 >= 0)
            && (sel_nr1 < gmx::ssize(*indexGroups)))
        {
            bRet = atoms_from_residuenumbers(atoms, (*indexGroups)[sel_nr1], nr, index);
            sprintf(gname, "atom_%s", (*indexGroups)[sel_nr1].name.c_str());
        }
    }
    else if (std::strncmp(*string, "ri", 2) == 0)
    {
        (*string) += 2;
        if (check_have_atoms(atoms, ostring) && parse_int_char(string, &sel_nr1, &c))
        {
            bRet = (select_residueindices(string, atoms, sel_nr1, c, nr, index, gname) != 0);
        }
    }
    else if ((*string)[0] == 'r')
    {
        (*string)++;
        if (check_have_atoms(atoms, ostring))
        {
            if (parse_int_char(string, &sel_nr1, &c))
            {
                bRet = (select_residuenumbers(string, atoms, sel_nr1, c, nr, index, gname) != 0);
            }
            else if (parse_names(string, &n_names, entryNames))
            {
                bRet = (select_residuenames(atoms, n_names, entryNames, nr, index) != 0);
                make_gname(n_names, entryNames, gname);
            }
        }
    }
    else if (std::strncmp(*string, "chain", 5) == 0)
    {
        (*string) += 5;
        if (check_have_atoms(atoms, ostring) && parse_names(string, &n_names, entryNames))
        {
            bRet = (select_chainnames(atoms, n_names, entryNames, nr, index) != 0);
            sprintf(gname, "ch%s", entryNames[0]);
            for (i = 1; i < n_names; i++)
            {
                std::strcat(gname, entryNames[i]);
            }
        }
    }
    if (bRet && bCompl)
    {
        snew(index1, natoms - *nr);
        nr1 = 0;
        for (i = 0; i < natoms; i++)
        {
            j = 0;
            while ((j < *nr) && (index[j] != i))
            {
                j++;
            }
            if (j == *nr)
            {
                if (nr1 >= natoms - *nr)
                {
                    printf("There are double atoms in your index group\n");
                    break;
                }
                index1[nr1] = i;
                nr1++;
            }
        }
        *nr = nr1;
        for (i = 0; i < nr1; i++)
        {
            index[i] = index1[i];
        }
        sfree(index1);

        for (i = std::strlen(gname) + 1; i > 0; i--)
        {
            gname[i] = gname[i - 1];
        }
        gname[0] = '!';
        printf("Complemented group: %d atoms\n", *nr);
    }

    return bRet;
}

static void list_residues(const t_atoms* atoms)
{
    int      i, j, start, end, prev_resind, resind;
    gmx_bool bDiff;

    /* Print all the residues, assuming continuous resnr count */
    start       = atoms->atom[0].resind;
    prev_resind = start;
    for (i = 0; i < atoms->nr; i++)
    {
        resind = atoms->atom[i].resind;
        if ((resind != prev_resind) || (i == atoms->nr - 1))
        {
            if ((bDiff = (std::strcmp(*atoms->resinfo[resind].name, *atoms->resinfo[start].name) != 0))
                || (i == atoms->nr - 1))
            {
                if (bDiff)
                {
                    end = prev_resind;
                }
                else
                {
                    end = resind;
                }
                if (end < start + 3)
                {
                    for (j = start; j <= end; j++)
                    {
                        printf("%4d %-5s", j + 1, *(atoms->resinfo[j].name));
                    }
                }
                else
                {
                    printf(" %4d - %4d %-5s  ", start + 1, end + 1, *(atoms->resinfo[start].name));
                }
                start = resind;
            }
        }
        prev_resind = resind;
    }
    printf("\n");
}

static void edit_index(int                      natoms,
                       const t_atoms*           atoms,
                       const rvec*              x,
                       std::vector<IndexGroup>* indexGroups,
                       gmx_bool                 bVerbose)
{
    char*    ostring;
    char     inp_string[STRLEN], *string;
    char     gname[STRLEN * 3], gname1[STRLEN], gname2[STRLEN];
    int      i, i0, i1, sel_nr, sel_nr2, newgroup;
    int      nr, nr1, nr2, *index, *index1, *index2;
    gmx_bool bAnd, bOr, bPrintOnce;

    string = nullptr;

    snew(index, natoms);
    snew(index1, natoms);
    snew(index2, natoms);

    newgroup   = NOTSET;
    bPrintOnce = TRUE;
    std::array<char*, MAXNAMES> entryNames;
    for (auto& name : entryNames)
    {
        snew(name, NAME_LEN + 1);
    }

    do
    {
        gname1[0] = '\0';
        if (bVerbose || bPrintOnce || newgroup != NOTSET)
        {
            printf("\n");
            if (bVerbose || bPrintOnce || newgroup == NOTSET)
            {
                i0 = 0;
                i1 = gmx::ssize(*indexGroups);
            }
            else
            {
                i0 = newgroup;
                i1 = newgroup + 1;
            }
            for (i = i0; i < i1; i++)
            {
                printf("%3d %-20s: %5td atoms\n",
                       i,
                       (*indexGroups)[i].name.c_str(),
                       gmx::ssize((*indexGroups)[i].particleIndices));
            }
            newgroup = NOTSET;
        }
        if (bVerbose || bPrintOnce)
        {
            printf("\n");
            printf(" nr : group      '!': not  'name' nr name   'splitch' nr    Enter: list "
                   "groups\n");
            printf(" 'a': atom       '&': and  'del' nr         'splitres' nr   'l': list "
                   "residues\n");
            printf(" 't': atom type  '|': or   'keep' nr        'splitat' nr    'h': help\n");
            printf(" 'r': residue              'res' nr         'chain' char\n");
            printf(" \"name\": group             'case': case %s         'q': save and quit\n",
                   bCase ? "insensitive" : "sensitive  ");
            printf(" 'ri': residue index\n");
            bPrintOnce = FALSE;
        }
        printf("\n");
        printf("> ");
        if (nullptr == fgets(inp_string, STRLEN, stdin))
        {
            gmx_fatal(FARGS, "Error reading user input");
        }
        inp_string[std::strlen(inp_string) - 1] = 0;
        printf("\n");
        string = inp_string;
        while (string[0] == ' ')
        {
            string++;
        }

        ostring = string;
        nr      = 0;
        if (string[0] == 'h')
        {
            printf(" nr                : selects an index group by number or quoted string.\n");
            printf("                     The string is first matched against the whole group "
                   "name,\n");
            printf("                     then against the beginning and finally against an\n");
            printf("                     arbitrary substring. A multiple match is an error.\n");

            printf(" 'a' nr1 [nr2 ...] : selects atoms, atom numbering starts at 1.\n");
            printf(" 'a' nr1 - nr2     : selects atoms in the range from nr1 to nr2.\n");
            printf(" 'a' name1[*] [name2[*] ...] : selects atoms by name(s), '?' matches any "
                   "char,\n");
            printf("                               wildcard '*' allowed at the end of a name.\n");
            printf(" 't' type1[*] [type2[*] ...] : as 'a', but for type, run input file "
                   "required.\n");
            printf(" 'r' nr1[ic1] [nr2[ic2] ...] : selects residues by number and insertion "
                   "code.\n");
            printf(" 'r' nr1 - nr2               : selects residues in the range from nr1 to "
                   "nr2.\n");
            printf(" 'r' name1[*] [name2[*] ...] : as 'a', but for residue names.\n");
            printf(" 'ri' nr1 - nr2              : selects residue indices, 1-indexed, (as opposed "
                   "to numbers) in the range from nr1 to nr2.\n");
            printf(" 'chain' ch1 [ch2 ...]       : selects atoms by chain identifier(s),\n");
            printf("                               not available with a .gro file as input.\n");
            printf(" !                 : takes the complement of a group with respect to all\n");
            printf("                     the atoms in the input file.\n");
            printf(" & |               : AND and OR, can be placed between any of the options\n");
            printf("                     above, the input is processed from left to right.\n");
            printf(" 'name' nr name    : rename group nr to name.\n");
            printf(" 'del' nr1 [- nr2] : deletes one group or groups in the range from nr1 to "
                   "nr2.\n");
            printf(" 'keep' nr         : deletes all groups except nr.\n");
            printf(" 'case'            : make all name compares case (in)sensitive.\n");
            printf(" 'splitch' nr      : split group into chains using CA distances.\n");
            printf(" 'splitres' nr     : split group into residues.\n");
            printf(" 'splitat' nr      : split group into atoms.\n");
            printf(" 'res' nr          : interpret numbers in group as residue numbers\n");
            printf(" Enter             : list the currently defined groups and commands\n");
            printf(" 'l'               : list the residues.\n");
            printf(" 'h'               : show this help.\n");
            printf(" 'q'               : save and quit.\n");
            printf("\n");
            printf(" Examples:\n");
            printf(" > 2 | 4 & r 3-5\n");
            printf(" selects all atoms from group 2 and 4 that have residue numbers 3, 4 or 5\n");
            printf(" > a C* & !a C CA\n");
            printf(" selects all atoms starting with 'C' but not the atoms 'C' and 'CA'\n");
            printf(" > \"protein\" & ! \"backb\"\n");
            printf(" selects all atoms that are in group 'protein' and not in group 'backbone'\n");
            if (bVerbose)
            {
                printf("\npress Enter ");
                getchar();
            }
        }
        else if (std::strncmp(string, "del", 3) == 0)
        {
            string += 3;
            if (parse_int(&string, &sel_nr))
            {
                while (string[0] == ' ')
                {
                    string++;
                }
                if (string[0] == '-')
                {
                    string++;
                    parse_int(&string, &sel_nr2);
                }
                else
                {
                    sel_nr2 = sel_nr;
                }
                while (string[0] == ' ')
                {
                    string++;
                }
                if (string[0] == '\0')
                {
                    remove_group(sel_nr, sel_nr2, indexGroups);
                }
                else
                {
                    printf("\nSyntax error: \"%s\"\n", string);
                }
            }
        }
        else if (std::strncmp(string, "keep", 4) == 0)
        {
            string += 4;
            if (parse_int(&string, &sel_nr))
            {
                remove_group(sel_nr + 1, gmx::ssize(*indexGroups) - 1, indexGroups);
                remove_group(0, sel_nr - 1, indexGroups);
            }
        }
        else if (std::strncmp(string, "name", 4) == 0)
        {
            string += 4;
            if (parse_int(&string, &sel_nr))
            {
                if ((sel_nr >= 0) && (sel_nr < gmx::ssize(*indexGroups)))
                {
                    sscanf(string, "%s", gname);
                    (*indexGroups)[sel_nr].name = gname;
                }
            }
        }
        else if (std::strncmp(string, "case", 4) == 0)
        {
            bCase = !bCase;
            printf("Switched to case %s\n", bCase ? "sensitive" : "insensitive");
        }
        else if (string[0] == 'v')
        {
            bVerbose = !bVerbose;
            printf("Turned verbose %s\n", bVerbose ? "on" : "off");
        }
        else if (string[0] == 'l')
        {
            if (check_have_atoms(atoms, ostring))
            {
                list_residues(atoms);
            }
        }
        else if (std::strncmp(string, "splitch", 7) == 0)
        {
            string += 7;
            if (check_have_atoms(atoms, ostring) && parse_int(&string, &sel_nr) && (sel_nr >= 0)
                && (sel_nr < gmx::ssize(*indexGroups)))
            {
                split_chain(atoms, x, sel_nr, indexGroups);
            }
        }
        else if (std::strncmp(string, "splitres", 8) == 0)
        {
            string += 8;
            if (check_have_atoms(atoms, ostring) && parse_int(&string, &sel_nr) && (sel_nr >= 0)
                && (sel_nr < gmx::ssize(*indexGroups)))
            {
                split_group(atoms, sel_nr, indexGroups, FALSE);
            }
        }
        else if (std::strncmp(string, "splitat", 7) == 0)
        {
            string += 7;
            if (check_have_atoms(atoms, ostring) && parse_int(&string, &sel_nr) && (sel_nr >= 0)
                && (sel_nr < gmx::ssize(*indexGroups)))
            {
                split_group(atoms, sel_nr, indexGroups, TRUE);
            }
        }
        else if (string[0] == '\0')
        {
            bPrintOnce = TRUE;
        }
        else if (string[0] != 'q')
        {
            nr2 = -1;
            if (parse_entry(&string, natoms, atoms, indexGroups, &nr, index, gname, entryNames))
            {
                do
                {
                    while (string[0] == ' ')
                    {
                        string++;
                    }

                    bAnd = FALSE;
                    bOr  = FALSE;
                    if (string[0] == '&')
                    {
                        bAnd = TRUE;
                    }
                    else if (string[0] == '|')
                    {
                        bOr = TRUE;
                    }

                    if (bAnd || bOr)
                    {
                        string++;
                        nr1 = nr;
                        for (i = 0; i < nr; i++)
                        {
                            index1[i] = index[i];
                        }
                        std::strcpy(gname1, gname);
                        if (parse_entry(&string, natoms, atoms, indexGroups, &nr2, index2, gname2, entryNames))
                        {
                            if (bOr)
                            {
                                or_groups(nr1, index1, nr2, index2, &nr, index);
                                sprintf(gname, "%s_%s", gname1, gname2);
                            }
                            else
                            {
                                and_groups(nr1, index1, nr2, index2, &nr, index);
                                sprintf(gname, "%s_&_%s", gname1, gname2);
                            }
                        }
                    }
                } while (bAnd || bOr);
            }
            while (string[0] == ' ')
            {
                string++;
            }
            if (string[0])
            {
                printf("\nSyntax error: \"%s\"\n", string);
            }
            else if (nr > 0)
            {
                indexGroups->push_back({ gname, { index, index + nr } });
            }
            else
            {
                printf("Group is empty\n");
            }
        }
    } while (string[0] != 'q');

    for (auto& name : entryNames)
    {
        sfree(name);
    }
    sfree(index);
    sfree(index1);
    sfree(index2);
}

//! Find the number of atoms in the system implied by the largest atom index
static int impliedNumberOfAtom(gmx::ArrayRef<const IndexGroup> indexGroups)
{
    int maxAtomIndex = -1;

    for (const auto& indexGroup : indexGroups)
    {
        for (const int a : indexGroup.particleIndices)
        {
            maxAtomIndex = std::max(maxAtomIndex, a);
        }
    }

    return maxAtomIndex + 1;
}

int gmx_make_ndx(int argc, char* argv[])
{
    const char* desc[] = { "Index groups are necessary for almost every GROMACS program.",
                           "All these programs can generate default index groups. You ONLY",
                           "have to use [THISMODULE] when you need SPECIAL index groups.",
                           "There is a default index group for the whole system, 9 default",
                           "index groups for proteins, and a default index group",
                           "is generated for every other residue name.",
                           "",
                           "When no index file is supplied, also [THISMODULE] will generate the",
                           "default groups.",
                           "With the index editor you can select on atom, residue and chain names",
                           "and numbers.",
                           "When a run input file is supplied you can also select on atom type.",
                           "You can use boolean operations, you can split groups",
                           "into chains, residues or atoms. You can delete and rename groups.",
                           "Type 'h' in the editor for more details.",
                           "",
                           "The atom numbering in the editor and the index file starts at 1.",
                           "",
                           "The [TT]-twin[tt] switch duplicates all index groups with an offset of",
                           "[TT]-natoms[tt], which is useful for Computational Electrophysiology",
                           "double-layer membrane setups.",
                           "",
                           "See also [gmx-select] [TT]-on[tt], which provides an alternative way",
                           "for constructing index groups.  It covers nearly all of [THISMODULE]",
                           "functionality, and in many cases much more." };

    static int      natoms     = 0;
    static gmx_bool bVerbose   = FALSE;
    static gmx_bool bDuplicate = FALSE;
    t_pargs         pa[]       = { { "-natoms",
                       FALSE,
                       etINT,
                       { &natoms },
                       "set number of atoms (default: read from coordinate or index file)" },
                     { "-twin",
                       FALSE,
                       etBOOL,
                       { &bDuplicate },
                       "Duplicate all index groups with an offset of -natoms" },
                     { "-verbose", FALSE, etBOOL, { &bVerbose }, "HIDDENVerbose output" } };
#define NPA asize(pa)

    gmx_output_env_t* oenv;
    const char*       stxfile;
    const char*       ndxoutfile;
    gmx_bool          bNatoms;
    t_atoms           atoms;
    rvec *            x, *v;
    PbcType           pbcType;
    matrix            box;
    t_filenm          fnm[] = { { efSTX, "-f", nullptr, ffOPTRD },
                       { efNDX, "-n", nullptr, ffOPTRDMULT },
                       { efNDX, "-o", nullptr, ffWRITE } };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, NPA, pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    stxfile                                     = ftp2fn_null(efSTX, NFILE, fnm);
    gmx::ArrayRef<const std::string> ndxInFiles = opt2fnsIfOptionSet("-n", NFILE, fnm);
    ndxoutfile                                  = opt2fn("-o", NFILE, fnm);
    bNatoms                                     = opt2parg_bSet("-natoms", NPA, pa);

    if (!stxfile && ndxInFiles.empty())
    {
        gmx_fatal(FARGS, "No input files (structure or index)");
    }

    gmx_mtop_t mtop;
    if (stxfile)
    {
        bool haveFullTopology = false;
        fprintf(stderr, "\nReading structure file\n");
        readConfAndTopology(stxfile, &haveFullTopology, &mtop, &pbcType, &x, &v, box);
        atoms = gmx_mtop_global_atoms(mtop);
        if (atoms.pdbinfo == nullptr)
        {
            snew(atoms.pdbinfo, atoms.nr);
        }
        natoms  = atoms.nr;
        bNatoms = TRUE;
    }
    else
    {
        x = nullptr;
    }

    /* read input file(s) */
    std::vector<IndexGroup> indexGroups;
    printf("Going to read %td old index file(s)\n", ndxInFiles.ssize());
    if (!ndxInFiles.empty())
    {
        for (const std::string& ndxInFile : ndxInFiles)
        {
            const auto indexGroups2 = init_index(ndxInFile.c_str());
            indexGroups.insert(indexGroups.begin(), indexGroups2.begin(), indexGroups2.end());
        }
    }
    else
    {
        indexGroups = analyse(&atoms, FALSE, TRUE);
    }

    if (!bNatoms)
    {
        natoms = impliedNumberOfAtom(indexGroups);
        printf("Deducing %d atoms in the system from indices in the index file\n", natoms);
    }
    edit_index(natoms, &atoms, x, &indexGroups, bVerbose);

    write_index(ndxoutfile, indexGroups, bDuplicate, natoms);

    if (stxfile)
    {
        sfree(v);
        sfree(x);
        done_atom(&atoms);
    }
    output_env_done(oenv);

    return 0;
}

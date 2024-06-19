/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2005- The GROMACS Authors
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

#include "espio.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static int get_espresso_word(FILE* fp, char word[])
{
    int ret, nc, i;

    ret = 0;
    nc  = 0;

    do
    {
        i = fgetc(fp);
        if (i != EOF)
        {
            if (i == ' ' || i == '\n' || i == '\t')
            {
                if (nc > 0)
                {
                    ret = 1;
                }
            }
            else if (i == '{')
            {
                if (nc == 0)
                {
                    word[nc++] = '{';
                }
                ret = 2;
            }
            else if (i == '}')
            {
                if (nc == 0)
                {
                    word[nc++] = '}';
                }
                ret = 3;
            }
            else
            {
                word[nc++] = static_cast<char>(i);
            }
        }
    } while (i != EOF && ret == 0);

    word[nc] = '\0';

    return ret;
}

static int check_open_parenthesis(FILE* fp, int r, const std::filesystem::path& infile, const char* keyword)
{
    int  level_inc;
    char word[STRLEN];

    level_inc = 0;
    if (r == 2)
    {
        level_inc++;
    }
    else
    {
        r = get_espresso_word(fp, word);
        if (r == 2)
        {
            level_inc++;
        }
        else
        {
            gmx_fatal(FARGS, "Expected '{' after '%s' in file '%s'", keyword, infile.string().c_str());
        }
    }

    return level_inc;
}

static int check_close_parenthesis(FILE* fp, int r, const std::filesystem::path& infile, const char* keyword)
{
    int  level_inc;
    char word[STRLEN];

    level_inc = 0;
    if (r == 3)
    {
        level_inc--;
    }
    else
    {
        r = get_espresso_word(fp, word);
        if (r == 3)
        {
            level_inc--;
        }
        else
        {
            gmx_fatal(
                    FARGS, "Expected '}' after section '%s' in file '%s'", keyword, infile.string().c_str());
        }
    }

    return level_inc;
}

enum
{
    espID,
    espPOS,
    espTYPE,
    espQ,
    espV,
    espF,
    espMOLECULE,
    espNR
};
static const char* const esp_prop[espNR] = { "id", "pos", "type", "q", "v", "f", "molecule" };

void gmx_espresso_read_conf(const std::filesystem::path& infile,
                            t_symtab*                    symtab,
                            char**                       name,
                            t_atoms*                     atoms,
                            rvec                         x[],
                            rvec*                        v,
                            matrix                       box)
{
    FILE*    fp;
    char     word[STRLEN], buf[STRLEN];
    int      level, r, nprop, p, i, m, molnr;
    int      prop[32];
    double   d;
    gmx_bool bFoundParticles, bFoundProp, bFoundVariable, bMol;

    if (name != nullptr)
    {
        // No title reading implemented for espresso files
        *name = gmx_strdup("");
    }

    clear_mat(box);

    atoms->haveMass    = FALSE;
    atoms->haveCharge  = FALSE;
    atoms->haveType    = FALSE;
    atoms->haveBState  = FALSE;
    atoms->havePdbInfo = FALSE;

    fp = gmx_fio_fopen(infile, "r");

    bFoundParticles = FALSE;
    bFoundVariable  = FALSE;
    bMol            = FALSE;
    level           = 0;
    while ((r = get_espresso_word(fp, word)))
    {
        if (level == 1 && std::strcmp(word, "particles") == 0 && !bFoundParticles)
        {
            bFoundParticles = TRUE;
            level += check_open_parenthesis(fp, r, infile, "particles");
            nprop = 0;
            while (level == 2 && (r = get_espresso_word(fp, word)))
            {
                bFoundProp = FALSE;
                for (p = 0; p < espNR; p++)
                {
                    if (strcmp(word, esp_prop[p]) == 0)
                    {
                        bFoundProp    = TRUE;
                        prop[nprop++] = p;
                        if (p == espQ)
                        {
                            atoms->haveCharge = TRUE;
                        }

                        if (debug)
                        {
                            fprintf(debug, "  prop[%d] = %s\n", nprop - 1, esp_prop[prop[nprop - 1]]);
                        }
                    }
                }
                if (!bFoundProp && word[0] != '}')
                {
                    gmx_fatal(FARGS, "Can not read Espresso files with particle property '%s'", word);
                }
                if (bFoundProp && p == espMOLECULE)
                {
                    bMol = TRUE;
                }
                if (r == 3)
                {
                    level--;
                }
            }

            i = 0;
            while (level > 0 && (r = get_espresso_word(fp, word)))
            {
                if (r == 2)
                {
                    level++;
                }
                else if (r == 3)
                {
                    level--;
                }
                if (level == 2)
                {
                    for (p = 0; p < nprop; p++)
                    {
                        switch (prop[p])
                        {
                            case espID:
                                r = get_espresso_word(fp, word);
                                /* Not used */
                                break;
                            case espPOS:
                                for (m = 0; m < 3; m++)
                                {
                                    r = get_espresso_word(fp, word);
                                    sscanf(word, "%lf", &d);
                                    x[i][m] = d;
                                }
                                break;
                            case espTYPE:
                                r                   = get_espresso_word(fp, word);
                                atoms->atom[i].type = std::strtol(word, nullptr, 10);
                                break;
                            case espQ:
                                r = get_espresso_word(fp, word);
                                sscanf(word, "%lf", &d);
                                atoms->atom[i].q = d;
                                break;
                            case espV:
                                for (m = 0; m < 3; m++)
                                {
                                    r = get_espresso_word(fp, word);
                                    sscanf(word, "%lf", &d);
                                    v[i][m] = d;
                                }
                                break;
                            case espF:
                                for (m = 0; m < 3; m++)
                                {
                                    r = get_espresso_word(fp, word);
                                    /* not used */
                                }
                                break;
                            case espMOLECULE:
                                r     = get_espresso_word(fp, word);
                                molnr = std::strtol(word, nullptr, 10);
                                if (i == 0 || atoms->resinfo[atoms->atom[i - 1].resind].nr != molnr)
                                {
                                    atoms->atom[i].resind = (i == 0 ? 0 : atoms->atom[i - 1].resind + 1);
                                    atoms->resinfo[atoms->atom[i].resind].nr      = molnr;
                                    atoms->resinfo[atoms->atom[i].resind].ic      = ' ';
                                    atoms->resinfo[atoms->atom[i].resind].chainid = ' ';
                                    atoms->resinfo[atoms->atom[i].resind].chainnum =
                                            molnr; /* Not sure if this is right? */
                                }
                                else
                                {
                                    atoms->atom[i].resind = atoms->atom[i - 1].resind;
                                }
                                break;
                        }
                    }
                    /* Generate an atom name from the particle type */
                    sprintf(buf, "T%hu", atoms->atom[i].type);
                    atoms->atomname[i] = put_symtab(symtab, buf);
                    if (bMol)
                    {
                        if (i == 0 || atoms->atom[i].resind != atoms->atom[i - 1].resind)
                        {
                            atoms->resinfo[atoms->atom[i].resind].name = put_symtab(symtab, "MOL");
                        }
                    }
                    else
                    {
                        /* Residue number is the atom number */
                        atoms->atom[i].resind = i;
                        /* Generate an residue name from the particle type */
                        if (atoms->atom[i].type < 26)
                        {
                            sprintf(buf, "T%c", 'A' + atoms->atom[i].type);
                        }
                        else
                        {
                            sprintf(buf, "T%c%c", 'A' + atoms->atom[i].type / 26, 'A' + atoms->atom[i].type % 26);
                        }
                        t_atoms_set_resinfo(atoms, i, symtab, buf, i, ' ', 0, ' ');
                    }

                    if (r == 3)
                    {
                        level--;
                    }
                    i++;
                }
            }
            atoms->nres = atoms->nr;

            if (i != atoms->nr)
            {
                gmx_fatal(FARGS,
                          "Internal inconsistency in Espresso routines, read %d atoms, expected %d "
                          "atoms",
                          i,
                          atoms->nr);
            }
        }
        else if (level == 1 && std::strcmp(word, "variable") == 0 && !bFoundVariable)
        {
            bFoundVariable = TRUE;
            level += check_open_parenthesis(fp, r, infile, "variable");
            while (level == 2 && (r = get_espresso_word(fp, word)))
            {
                if (level == 2 && std::strcmp(word, "box_l") == 0)
                {
                    for (m = 0; m < 3; m++)
                    {
                        r = get_espresso_word(fp, word);
                        sscanf(word, "%lf", &d);
                        box[m][m] = d;
                    }
                    level += check_close_parenthesis(fp, r, infile, "box_l");
                }
            }
        }
        else if (r == 2)
        {
            level++;
        }
        else if (r == 3)
        {
            level--;
        }
    }

    if (!bFoundParticles)
    {
        fprintf(stderr, "Did not find a particles section in Espresso file '%s'\n", infile.string().c_str());
    }

    gmx_fio_fclose(fp);
}

int get_espresso_coordnum(const std::filesystem::path& infile)
{
    FILE*    fp;
    char     word[STRLEN];
    int      natoms, level, r;
    gmx_bool bFoundParticles;

    natoms = 0;

    fp = gmx_fio_fopen(infile, "r");

    bFoundParticles = FALSE;
    level           = 0;
    while ((r = get_espresso_word(fp, word)) && !bFoundParticles)
    {
        if (level == 1 && strcmp(word, "particles") == 0 && !bFoundParticles)
        {
            bFoundParticles = TRUE;
            level += check_open_parenthesis(fp, r, infile, "particles");
            while (level > 0 && (r = get_espresso_word(fp, word)))
            {
                if (r == 2)
                {
                    level++;
                    if (level == 2)
                    {
                        natoms++;
                    }
                }
                else if (r == 3)
                {
                    level--;
                }
            }
        }
        else if (r == 2)
        {
            level++;
        }
        else if (r == 3)
        {
            level--;
        }
    }
    if (!bFoundParticles)
    {
        fprintf(stderr, "Did not find a particles section in Espresso file '%s'\n", infile.string().c_str());
    }

    gmx_fio_fclose(fp);

    return natoms;
}

void write_espresso_conf_indexed(FILE*          out,
                                 const char*    title,
                                 const t_atoms* atoms,
                                 int            nx,
                                 const int*     index,
                                 const rvec*    x,
                                 const rvec*    v,
                                 const matrix   box)
{
    int i, j;

    fprintf(out, "# %s\n", title);
    if (TRICLINIC(box))
    {
        gmx_warning("The Espresso format does not support triclinic unit-cells");
    }
    fprintf(out, "{variable {box_l %f %f %f}}\n", box[0][0], box[1][1], box[2][2]);

    fprintf(out, "{particles {id pos type q%s}\n", v ? " v" : "");
    for (i = 0; i < nx; i++)
    {
        if (index)
        {
            j = index[i];
        }
        else
        {
            j = i;
        }
        fprintf(out,
                "\t{%d %f %f %f %hu %g",
                j,
                x[j][XX],
                x[j][YY],
                x[j][ZZ],
                atoms->atom[j].type,
                atoms->atom[j].q);
        if (v)
        {
            fprintf(out, " %f %f %f", v[j][XX], v[j][YY], v[j][ZZ]);
        }
        fprintf(out, "}\n");
    }
    fprintf(out, "}\n");
}

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "confio.h"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/g96io.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static int get_espresso_word(FILE *fp, char word[])
{
    int  ret, nc, i;

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
                word[nc++] = (char)i;
            }
        }
    }
    while (i != EOF && ret == 0);

    word[nc] = '\0';

    /*  printf("'%s'\n",word); */

    return ret;
}

static int check_open_parenthesis(FILE *fp, int r,
                                  const char *infile, const char *keyword)
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
            gmx_fatal(FARGS, "Expected '{' after '%s' in file '%s'",
                      keyword, infile);
        }
    }

    return level_inc;
}

static int check_close_parenthesis(FILE *fp, int r,
                                   const char *infile, const char *keyword)
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
            gmx_fatal(FARGS, "Expected '}' after section '%s' in file '%s'",
                      keyword, infile);
        }
    }

    return level_inc;
}

enum {
    espID, espPOS, espTYPE, espQ, espV, espF, espMOLECULE, espNR
};
const char *esp_prop[espNR] = {
    "id", "pos", "type", "q", "v", "f",
    "molecule"
};

static void read_espresso_conf(const char *infile, char *title,
                               t_atoms *atoms, rvec x[], rvec *v, matrix box)
{
    t_symtab *symtab = NULL;
    FILE     *fp;
    char      word[STRLEN], buf[STRLEN];
    int       level, r, nprop, p, i, m, molnr;
    int       prop[32];
    double    d;
    gmx_bool  bFoundParticles, bFoundProp, bFoundVariable, bMol;

    if (!symtab)
    {
        snew(symtab, 1);
        open_symtab(symtab);
    }
    // TODO: The code does not understand titles it writes...
    title[0] = '\0';

    clear_mat(box);

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
            level          += check_open_parenthesis(fp, r, infile, "particles");
            nprop           = 0;
            while (level == 2 && (r = get_espresso_word(fp, word)))
            {
                bFoundProp = FALSE;
                for (p = 0; p < espNR; p++)
                {
                    if (strcmp(word, esp_prop[p]) == 0)
                    {
                        bFoundProp    = TRUE;
                        prop[nprop++] = p;
                        /* printf("  prop[%d] = %s\n",nprop-1,esp_prop[prop[nprop-1]]); */
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
                                atoms->atom[i].type = std::strtol(word, NULL, 10);
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
                                molnr = std::strtol(word, NULL, 10);
                                if (i == 0 ||
                                    atoms->resinfo[atoms->atom[i-1].resind].nr != molnr)
                                {
                                    atoms->atom[i].resind =
                                        (i == 0 ? 0 : atoms->atom[i-1].resind+1);
                                    atoms->resinfo[atoms->atom[i].resind].nr       = molnr;
                                    atoms->resinfo[atoms->atom[i].resind].ic       = ' ';
                                    atoms->resinfo[atoms->atom[i].resind].chainid  = ' ';
                                    atoms->resinfo[atoms->atom[i].resind].chainnum = molnr; /* Not sure if this is right? */
                                }
                                else
                                {
                                    atoms->atom[i].resind = atoms->atom[i-1].resind;
                                }
                                break;
                        }
                    }
                    /* Generate an atom name from the particle type */
                    sprintf(buf, "T%d", atoms->atom[i].type);
                    atoms->atomname[i] = put_symtab(symtab, buf);
                    if (bMol)
                    {
                        if (i == 0 || atoms->atom[i].resind != atoms->atom[i-1].resind)
                        {
                            atoms->resinfo[atoms->atom[i].resind].name =
                                put_symtab(symtab, "MOL");
                        }
                    }
                    else
                    {
                        /* Residue number is the atom number */
                        atoms->atom[i].resind = i;
                        /* Generate an residue name from the particle type */
                        if (atoms->atom[i].type < 26)
                        {
                            sprintf(buf, "T%c", 'A'+atoms->atom[i].type);
                        }
                        else
                        {
                            sprintf(buf, "T%c%c",
                                    'A'+atoms->atom[i].type/26, 'A'+atoms->atom[i].type%26);
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
                gmx_fatal(FARGS, "Internal inconsistency in Espresso routines, read %d atoms, expected %d atoms", i, atoms->nr);
            }
        }
        else if (level == 1 && std::strcmp(word, "variable") == 0 && !bFoundVariable)
        {
            bFoundVariable = TRUE;
            level         += check_open_parenthesis(fp, r, infile, "variable");
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
        fprintf(stderr, "Did not find a particles section in Espresso file '%s'\n",
                infile);
    }

    gmx_fio_fclose(fp);
}

static int get_espresso_coordnum(const char *infile)
{
    FILE    *fp;
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
            level          += check_open_parenthesis(fp, r, infile, "particles");
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
        fprintf(stderr, "Did not find a particles section in Espresso file '%s'\n",
                infile);
    }

    gmx_fio_fclose(fp);

    return natoms;
}

static void write_espresso_conf_indexed(FILE *out, const char *title,
                                        t_atoms *atoms, int nx, atom_id *index,
                                        rvec *x, rvec *v, matrix box)
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
        fprintf(out, "\t{%d %f %f %f %d %g",
                j, x[j][XX], x[j][YY], x[j][ZZ],
                atoms->atom[j].type, atoms->atom[j].q);
        if (v)
        {
            fprintf(out, " %f %f %f", v[j][XX], v[j][YY], v[j][ZZ]);
        }
        fprintf(out, "}\n");
    }
    fprintf(out, "}\n");
}

static void get_coordnum_fp (FILE *in, char *title, int *natoms)
{
    char line[STRLEN+1];

    fgets2 (title, STRLEN, in);
    fgets2 (line, STRLEN, in);
    if (sscanf (line, "%d", natoms) != 1)
    {
        gmx_fatal(FARGS, "gro file does not have the number of atoms on the second line");
    }
}

static void get_coordnum (const char *infile, int *natoms)
{
    FILE *in;
    char  title[STRLEN];

    in = gmx_fio_fopen(infile, "r");
    get_coordnum_fp(in, title, natoms);
    gmx_fio_fclose (in);
}

static gmx_bool get_w_conf(FILE *in, const char *infile, char *title,
                           t_symtab *symtab, t_atoms *atoms, int *ndec,
                           rvec x[], rvec *v, matrix box)
{
    char       name[6];
    char       resname[6], oldresname[6];
    char       line[STRLEN+1], *ptr;
    char       buf[256];
    double     x1, y1, z1, x2, y2, z2;
    rvec       xmin, xmax;
    int        natoms, i, m, resnr, newres, oldres, ddist, c;
    gmx_bool   bFirst, bVel;
    char      *p1, *p2, *p3;

    newres  = -1;
    oldres  = -12345; /* Unlikely number for the first residue! */
    ddist   = 0;

    /* Read the title and number of atoms */
    get_coordnum_fp(in, title, &natoms);

    if (natoms > atoms->nr)
    {
        gmx_fatal(FARGS, "gro file contains more atoms (%d) than expected (%d)",
                  natoms, atoms->nr);
    }
    else if (natoms <  atoms->nr)
    {
        fprintf(stderr, "Warning: gro file contains less atoms (%d) than expected"
                " (%d)\n", natoms, atoms->nr);
    }

    bFirst = TRUE;

    bVel = FALSE;

    resname[0]     = '\0';
    oldresname[0]  = '\0';

    /* just pray the arrays are big enough */
    for (i = 0; (i < natoms); i++)
    {
        if ((fgets2 (line, STRLEN, in)) == NULL)
        {
            gmx_fatal(FARGS, "Unexpected end of file in file %s at line %d",
                      infile, i+2);
        }
        if (strlen(line) < 39)
        {
            gmx_fatal(FARGS, "Invalid line in %s for atom %d:\n%s", infile, i+1, line);
        }

        /* determine read precision from distance between periods
           (decimal points) */
        if (bFirst)
        {
            bFirst = FALSE;
            p1     = strchr(line, '.');
            if (p1 == NULL)
            {
                gmx_fatal(FARGS, "A coordinate in file %s does not contain a '.'", infile);
            }
            p2 = strchr(&p1[1], '.');
            if (p2 == NULL)
            {
                gmx_fatal(FARGS, "A coordinate in file %s does not contain a '.'", infile);
            }
            ddist = p2 - p1;
            *ndec = ddist - 5;

            p3 = strchr(&p2[1], '.');
            if (p3 == NULL)
            {
                gmx_fatal(FARGS, "A coordinate in file %s does not contain a '.'", infile);
            }

            if (p3 - p2 != ddist)
            {
                gmx_fatal(FARGS, "The spacing of the decimal points in file %s is not consistent for x, y and z", infile);
            }
        }

        /* residue number*/
        memcpy(name, line, 5);
        name[5] = '\0';
        sscanf(name, "%d", &resnr);
        sscanf(line+5, "%5s", resname);

        if (resnr != oldres || strncmp(resname, oldresname, sizeof(resname)))
        {
            oldres = resnr;
            newres++;
            if (newres >= natoms)
            {
                gmx_fatal(FARGS, "More residues than atoms in %s (natoms = %d)",
                          infile, natoms);
            }
            atoms->atom[i].resind = newres;
            t_atoms_set_resinfo(atoms, i, symtab, resname, resnr, ' ', 0, ' ');
        }
        else
        {
            atoms->atom[i].resind = newres;
        }

        /* atomname */
        std::memcpy(name, line+10, 5);
        atoms->atomname[i] = put_symtab(symtab, name);

        /* Copy resname to oldresname after we are done with the sanity check above */
        std::strncpy(oldresname, resname, sizeof(oldresname));

        /* eventueel controle atomnumber met i+1 */

        /* coordinates (start after residue data) */
        ptr = line + 20;
        /* Read fixed format */
        for (m = 0; m < DIM; m++)
        {
            for (c = 0; (c < ddist && ptr[0]); c++)
            {
                buf[c] = ptr[0];
                ptr++;
            }
            buf[c] = '\0';
            if (sscanf (buf, "%lf %lf", &x1, &x2) != 1)
            {
                gmx_fatal(FARGS, "Something is wrong in the coordinate formatting of file %s. Note that gro is fixed format (see the manual)", infile);
            }
            else
            {
                x[i][m] = x1;
            }
        }

        /* velocities (start after residues and coordinates) */
        if (v)
        {
            /* Read fixed format */
            for (m = 0; m < DIM; m++)
            {
                for (c = 0; (c < ddist && ptr[0]); c++)
                {
                    buf[c] = ptr[0];
                    ptr++;
                }
                buf[c] = '\0';
                if (sscanf (buf, "%lf", &x1) != 1)
                {
                    v[i][m] = 0;
                }
                else
                {
                    v[i][m] = x1;
                    bVel    = TRUE;
                }
            }
        }
    }
    atoms->nres = newres + 1;

    /* box */
    fgets2 (line, STRLEN, in);
    if (sscanf (line, "%lf%lf%lf", &x1, &y1, &z1) != 3)
    {
        gmx_warning("Bad box in file %s", infile);

        /* Generate a cubic box */
        for (m = 0; (m < DIM); m++)
        {
            xmin[m] = xmax[m] = x[0][m];
        }
        for (i = 1; (i < atoms->nr); i++)
        {
            for (m = 0; (m < DIM); m++)
            {
                xmin[m] = std::min(xmin[m], x[i][m]);
                xmax[m] = std::max(xmax[m], x[i][m]);
            }
        }
        for (i = 0; i < DIM; i++)
        {
            for (m = 0; m < DIM; m++)
            {
                box[i][m] = 0.0;
            }
        }
        for (m = 0; (m < DIM); m++)
        {
            box[m][m] = (xmax[m]-xmin[m]);
        }
        fprintf(stderr, "Generated a cubic box %8.3f x %8.3f x %8.3f\n",
                box[XX][XX], box[YY][YY], box[ZZ][ZZ]);
    }
    else
    {
        /* We found the first three values, the diagonal elements */
        box[XX][XX] = x1;
        box[YY][YY] = y1;
        box[ZZ][ZZ] = z1;
        if (sscanf (line, "%*f%*f%*f%lf%lf%lf%lf%lf%lf",
                    &x1, &y1, &z1, &x2, &y2, &z2) != 6)
        {
            x1 = y1 = z1 = x2 = y2 = z2 = 0.0;
        }
        box[XX][YY] = x1;
        box[XX][ZZ] = y1;
        box[YY][XX] = z1;
        box[YY][ZZ] = x2;
        box[ZZ][XX] = y2;
        box[ZZ][YY] = z2;
    }

    return bVel;
}

static void read_whole_conf(const char *infile, char *title,
                            t_atoms *atoms, rvec x[], rvec *v, matrix box)
{
    FILE    *in;
    int      ndec;
    t_symtab symtab;

    /* open file */
    in = gmx_fio_fopen(infile, "r");

    open_symtab(&symtab);
    get_w_conf(in, infile, title, &symtab, atoms, &ndec, x, v, box);
    /* We can't free the symbols, as they are still used in atoms, so
     * the only choice is to leak them. */
    free_symtab(&symtab);

    gmx_fio_fclose(in);
}

static gmx_bool gmx_one_before_eof(FILE *fp)
{
    char     data[4];
    gmx_bool beof;

    if ((beof = fread(data, 1, 1, fp)) == 1)
    {
        gmx_fseek(fp, -1, SEEK_CUR);
    }
    return !beof;
}

gmx_bool gro_next_x_or_v(FILE *status, t_trxframe *fr)
{
    t_atoms  atoms;
    t_symtab symtab;
    char     title[STRLEN], *p;
    double   tt;
    int      ndec = 0, i;

    if (gmx_one_before_eof(status))
    {
        return FALSE;
    }

    open_symtab(&symtab);
    atoms.nr = fr->natoms;
    snew(atoms.atom, fr->natoms);
    atoms.nres = fr->natoms;
    snew(atoms.resinfo, fr->natoms);
    snew(atoms.atomname, fr->natoms);

    fr->bV    = get_w_conf(status, title, title, &symtab, &atoms, &ndec, fr->x, fr->v, fr->box);
    fr->bPrec = TRUE;
    fr->prec  = 1;
    /* prec = 10^ndec: */
    for (i = 0; i < ndec; i++)
    {
        fr->prec *= 10;
    }
    fr->title  = title;
    fr->bTitle = TRUE;
    fr->bX     = TRUE;
    fr->bBox   = TRUE;

    sfree(atoms.atom);
    sfree(atoms.resinfo);
    sfree(atoms.atomname);
    done_symtab(&symtab);

    if ((p = strstr(title, "t=")) != NULL)
    {
        p += 2;
        if (sscanf(p, "%lf", &tt) == 1)
        {
            fr->time  = tt;
            fr->bTime = TRUE;
        }
        else
        {
            fr->time  = 0;
            fr->bTime = FALSE;
        }
    }

    if (atoms.nr != fr->natoms)
    {
        gmx_fatal(FARGS, "Number of atoms in gro frame (%d) doesn't match the number in the previous frame (%d)", atoms.nr, fr->natoms);
    }

    return TRUE;
}

int gro_first_x_or_v(FILE *status, t_trxframe *fr)
{
    char title[STRLEN];

    frewind(status);
    fprintf(stderr, "Reading frames from gro file");
    get_coordnum_fp(status, title, &fr->natoms);
    frewind(status);
    fprintf(stderr, " '%s', %d atoms.\n", title, fr->natoms);
    fr->bTitle = TRUE;
    fr->title  = title;
    if (fr->natoms == 0)
    {
        gmx_file("No coordinates in gro file");
    }

    snew(fr->x, fr->natoms);
    snew(fr->v, fr->natoms);
    gro_next_x_or_v(status, fr);

    return fr->natoms;
}

static void make_hconf_format(int pr, gmx_bool bVel, char format[])
{
    int l, vpr;

    /* build format string for printing,
       something like "%8.3f" for x and "%8.4f" for v */
    if (pr < 0)
    {
        pr = 0;
    }
    if (pr > 30)
    {
        pr = 30;
    }
    l   = pr+5;
    vpr = pr+1;
    if (bVel)
    {
        sprintf(format, "%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df\n",
                l, pr, l, pr, l, pr, l, vpr, l, vpr, l, vpr);
    }
    else
    {
        sprintf(format, "%%%d.%df%%%d.%df%%%d.%df\n", l, pr, l, pr, l, pr);
    }

}

static void write_hconf_box(FILE *out, int pr, matrix box)
{
    char format[100];
    int  l;

    if (pr < 5)
    {
        pr = 5;
    }
    l = pr+5;

    if (box[XX][YY] || box[XX][ZZ] || box[YY][XX] || box[YY][ZZ] ||
        box[ZZ][XX] || box[ZZ][YY])
    {
        sprintf(format, "%%%d.%df%%%d.%df%%%d.%df"
                "%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df\n",
                l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr);
        fprintf(out, format,
                box[XX][XX], box[YY][YY], box[ZZ][ZZ],
                box[XX][YY], box[XX][ZZ], box[YY][XX],
                box[YY][ZZ], box[ZZ][XX], box[ZZ][YY]);
    }
    else
    {
        sprintf(format, "%%%d.%df%%%d.%df%%%d.%df\n", l, pr, l, pr, l, pr);
        fprintf(out, format,
                box[XX][XX], box[YY][YY], box[ZZ][ZZ]);
    }
}

void write_hconf_indexed_p(FILE *out, const char *title, t_atoms *atoms,
                           int nx, const atom_id index[], int pr,
                           rvec *x, rvec *v, matrix box)
{
    char resnm[6], nm[6], format[100];
    int  ai, i, resind, resnr;

    bromacs(format, 99);
    fprintf (out, "%s\n", (title && title[0]) ? title : format);
    fprintf (out, "%5d\n", nx);

    make_hconf_format(pr, v != NULL, format);

    for (i = 0; (i < nx); i++)
    {
        ai = index[i];

        resind = atoms->atom[ai].resind;
        std::strncpy(resnm, " ??? ", sizeof(resnm)-1);
        if (resind < atoms->nres)
        {
            std::strncpy(resnm, *atoms->resinfo[resind].name, sizeof(resnm)-1);
            resnr = atoms->resinfo[resind].nr;
        }
        else
        {
            std::strncpy(resnm, " ??? ", sizeof(resnm)-1);
            resnr = resind + 1;
        }

        if (atoms->atom)
        {
            std::strncpy(nm, *atoms->atomname[ai], sizeof(nm)-1);
        }
        else
        {
            std::strncpy(nm, " ??? ", sizeof(nm)-1);
        }

        fprintf(out, "%5d%-5.5s%5.5s%5d", resnr%100000, resnm, nm, (ai+1)%100000);
        /* next fprintf uses built format string */
        if (v)
        {
            fprintf(out, format,
                    x[ai][XX], x[ai][YY], x[ai][ZZ], v[ai][XX], v[ai][YY], v[ai][ZZ]);
        }
        else
        {
            fprintf(out, format,
                    x[ai][XX], x[ai][YY], x[ai][ZZ]);
        }
    }

    write_hconf_box(out, pr, box);

    fflush(out);
}

static void write_hconf_mtop(FILE *out, const char *title, gmx_mtop_t *mtop,
                             int pr,
                             rvec *x, rvec *v, matrix box)
{
    char                    format[100];
    int                     i, resnr;
    gmx_mtop_atomloop_all_t aloop;
    t_atom                 *atom;
    char                   *atomname, *resname;

    bromacs(format, 99);
    fprintf (out, "%s\n", (title && title[0]) ? title : format);
    fprintf (out, "%5d\n", mtop->natoms);

    make_hconf_format(pr, v != NULL, format);

    aloop = gmx_mtop_atomloop_all_init(mtop);
    while (gmx_mtop_atomloop_all_next(aloop, &i, &atom))
    {
        gmx_mtop_atomloop_all_names(aloop, &atomname, &resnr, &resname);

        fprintf(out, "%5d%-5.5s%5.5s%5d",
                resnr%100000, resname, atomname, (i+1)%100000);
        /* next fprintf uses built format string */
        if (v)
        {
            fprintf(out, format,
                    x[i][XX], x[i][YY], x[i][ZZ], v[i][XX], v[i][YY], v[i][ZZ]);
        }
        else
        {
            fprintf(out, format,
                    x[i][XX], x[i][YY], x[i][ZZ]);
        }
    }

    write_hconf_box(out, pr, box);

    fflush(out);
}

void write_hconf_p(FILE *out, const char *title, t_atoms *atoms, int pr,
                   rvec *x, rvec *v, matrix box)
{
    atom_id *aa;
    int      i;

    snew(aa, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        aa[i] = i;
    }
    write_hconf_indexed_p(out, title, atoms, atoms->nr, aa, pr, x, v, box);
    sfree(aa);
}

void write_conf_p(const char *outfile, const char *title,
                  t_atoms *atoms, int pr,
                  rvec *x, rvec *v, matrix box)
{
    FILE *out;

    out = gmx_fio_fopen(outfile, "w");
    write_hconf_p(out, title, atoms, pr, x, v, box);

    gmx_fio_fclose (out);
}

static void write_conf(const char *outfile, const char *title, t_atoms *atoms,
                       rvec *x, rvec *v, matrix box)
{
    write_conf_p(outfile, title, atoms, 3, x, v, box);
}

void write_sto_conf_indexed(const char *outfile, const char *title,
                            t_atoms *atoms,
                            rvec x[], rvec *v, int ePBC, matrix box,
                            atom_id nindex, atom_id index[])
{
    FILE       *out;
    int         ftp;
    t_trxframe  fr;

    ftp = fn2ftp(outfile);
    switch (ftp)
    {
        case efGRO:
            out = gmx_fio_fopen(outfile, "w");
            write_hconf_indexed_p(out, title, atoms, nindex, index, 3, x, v, box);
            gmx_fio_fclose(out);
            break;
        case efG96:
            clear_trxframe(&fr, TRUE);
            fr.bTitle = TRUE;
            fr.title  = title;
            fr.natoms = atoms->nr;
            fr.bAtoms = TRUE;
            fr.atoms  = atoms;
            fr.bX     = TRUE;
            fr.x      = x;
            if (v)
            {
                fr.bV = TRUE;
                fr.v  = v;
            }
            fr.bBox = TRUE;
            copy_mat(box, fr.box);
            out = gmx_fio_fopen(outfile, "w");
            write_g96_conf(out, &fr, nindex, index);
            gmx_fio_fclose(out);
            break;
        case efPDB:
        case efBRK:
        case efENT:
        case efPQR:
            out = gmx_fio_fopen(outfile, "w");
            write_pdbfile_indexed(out, title, atoms, x, ePBC, box, ' ', -1, nindex, index, NULL, TRUE);
            gmx_fio_fclose(out);
            break;
        case efESP:
            out = gmx_fio_fopen(outfile, "w");
            write_espresso_conf_indexed(out, title, atoms, nindex, index, x, v, box);
            gmx_fio_fclose(out);
            break;
        case efTPR:
            gmx_fatal(FARGS, "Sorry, can not write a topology to %s", outfile);
            break;
        default:
            gmx_incons("Not supported in write_sto_conf_indexed");
    }
}

void write_sto_conf(const char *outfile, const char *title, t_atoms *atoms,
                    rvec x[], rvec *v, int ePBC, matrix box)
{
    FILE       *out;
    int         ftp;
    t_trxframe  fr;

    ftp = fn2ftp(outfile);
    switch (ftp)
    {
        case efGRO:
            write_conf(outfile, title, atoms, x, v, box);
            break;
        case efG96:
            clear_trxframe(&fr, TRUE);
            fr.bTitle = TRUE;
            fr.title  = title;
            fr.natoms = atoms->nr;
            fr.bAtoms = TRUE;
            fr.atoms  = atoms;
            fr.bX     = TRUE;
            fr.x      = x;
            if (v)
            {
                fr.bV = TRUE;
                fr.v  = v;
            }
            fr.bBox = TRUE;
            copy_mat(box, fr.box);
            out = gmx_fio_fopen(outfile, "w");
            write_g96_conf(out, &fr, -1, NULL);
            gmx_fio_fclose(out);
            break;
        case efPDB:
        case efBRK:
        case efENT:
            out = gmx_fio_fopen(outfile, "w");
            write_pdbfile(out, title, atoms, x, ePBC, box, ' ', -1, NULL, TRUE);
            gmx_fio_fclose(out);
            break;
        case efESP:
            out = gmx_fio_fopen(outfile, "w");
            write_espresso_conf_indexed(out, title, atoms, atoms->nr, NULL, x, v, box);
            gmx_fio_fclose(out);
            break;
        case efTPR:
            gmx_fatal(FARGS, "Sorry, can not write a topology to %s", outfile);
            break;
        default:
            gmx_incons("Not supported in write_sto_conf");
    }
}

void write_sto_conf_mtop(const char *outfile, const char *title,
                         gmx_mtop_t *mtop,
                         rvec x[], rvec *v, int ePBC, matrix box)
{
    int     ftp;
    FILE   *out;
    t_atoms atoms;

    ftp = fn2ftp(outfile);
    switch (ftp)
    {
        case efGRO:
            out = gmx_fio_fopen(outfile, "w");
            write_hconf_mtop(out, title, mtop, 3, x, v, box);
            gmx_fio_fclose(out);
            break;
        default:
            /* This is a brute force approach which requires a lot of memory.
             * We should implement mtop versions of all writing routines.
             */
            atoms = gmx_mtop_global_atoms(mtop);

            write_sto_conf(outfile, title, &atoms, x, v, ePBC, box);

            done_atom(&atoms);
            break;
    }
}

void get_stx_coordnum(const char *infile, int *natoms)
{
    FILE      *in;
    int        ftp, tpxver, tpxgen;
    t_trxframe fr;
    char       g96_line[STRLEN+1];

    ftp = fn2ftp(infile);
    range_check(ftp, 0, efNR);
    switch (ftp)
    {
        case efGRO:
            get_coordnum(infile, natoms);
            break;
        case efG96:
            in        = gmx_fio_fopen(infile, "r");
            fr.title  = NULL;
            fr.natoms = -1;
            fr.atoms  = NULL;
            fr.x      = NULL;
            fr.v      = NULL;
            fr.f      = NULL;
            *natoms   = read_g96_conf(in, infile, &fr, g96_line);
            gmx_fio_fclose(in);
            break;
        case efPDB:
        case efBRK:
        case efENT:
            in = gmx_fio_fopen(infile, "r");
            get_pdb_coordnum(in, natoms);
            gmx_fio_fclose(in);
            break;
        case efESP:
            *natoms = get_espresso_coordnum(infile);
            break;
        case efTPR:
        {
            t_tpxheader tpx;

            read_tpxheader(infile, &tpx, TRUE, &tpxver, &tpxgen);
            *natoms = tpx.natoms;
            break;
        }
        default:
            gmx_fatal(FARGS, "File type %s not supported in get_stx_coordnum",
                      ftp2ext(ftp));
    }
}

static void tpx_make_chain_identifiers(t_atoms *atoms, t_block *mols)
{
    int  m, a, a0, a1, r;
    char c, chainid;
    int  chainnum;

    /* We always assign a new chain number, but save the chain id characters
     * for larger molecules.
     */
#define CHAIN_MIN_ATOMS 15

    chainnum = 0;
    chainid  = 'A';
    for (m = 0; m < mols->nr; m++)
    {
        a0 = mols->index[m];
        a1 = mols->index[m+1];
        if ((a1-a0 >= CHAIN_MIN_ATOMS) && (chainid <= 'Z'))
        {
            c = chainid;
            chainid++;
        }
        else
        {
            c = ' ';
        }
        for (a = a0; a < a1; a++)
        {
            atoms->resinfo[atoms->atom[a].resind].chainnum = chainnum;
            atoms->resinfo[atoms->atom[a].resind].chainid  = c;
        }
        chainnum++;
    }

    /* Blank out the chain id if there was only one chain */
    if (chainid == 'B')
    {
        for (r = 0; r < atoms->nres; r++)
        {
            atoms->resinfo[r].chainid = ' ';
        }
    }
}

void read_stx_conf(const char *infile, char *title, t_atoms *atoms,
                   rvec x[], rvec *v, int *ePBC, matrix box)
{
    FILE       *in;
    gmx_mtop_t *mtop;
    t_topology  top;
    t_trxframe  fr;
    int         i, ftp, natoms;
    char        g96_line[STRLEN+1];

    if (atoms->nr == 0)
    {
        fprintf(stderr, "Warning: Number of atoms in %s is 0\n", infile);
    }
    else if (atoms->atom == NULL)
    {
        gmx_mem("Uninitialized array atom");
    }

    if (ePBC)
    {
        *ePBC = -1;
    }

    ftp = fn2ftp(infile);
    switch (ftp)
    {
        case efGRO:
            read_whole_conf(infile, title, atoms, x, v, box);
            break;
        case efG96:
            fr.title  = NULL;
            fr.natoms = atoms->nr;
            fr.atoms  = atoms;
            fr.x      = x;
            fr.v      = v;
            fr.f      = NULL;
            in        = gmx_fio_fopen(infile, "r");
            read_g96_conf(in, infile, &fr, g96_line);
            gmx_fio_fclose(in);
            copy_mat(fr.box, box);
            std::strncpy(title, fr.title, STRLEN);
            break;
        case efPDB:
        case efBRK:
        case efENT:
            read_pdb_conf(infile, title, atoms, x, ePBC, box, TRUE, NULL);
            break;
        case efESP:
            read_espresso_conf(infile, title, atoms, x, v, box);
            break;
        case efTPR:
            snew(mtop, 1);
            i = read_tpx(infile, NULL, box, &natoms, x, v, NULL, mtop);
            if (ePBC)
            {
                *ePBC = i;
            }

            strcpy(title, *(mtop->name));

            /* Free possibly allocated memory */
            done_atom(atoms);

            *atoms = gmx_mtop_global_atoms(mtop);
            top    = gmx_mtop_t_to_t_topology(mtop);
            tpx_make_chain_identifiers(atoms, &top.mols);

            sfree(mtop);
            /* The strings in the symtab are still in use in the returned t_atoms
             * structure, so we should not free them. But there is no place to put the
             * symbols; the only choice is to leak the memory...
             * So we clear the symbol table before freeing the topology structure. */
            free_symtab(&top.symtab);
            done_top(&top);

            break;
        default:
            gmx_incons("Not supported in read_stx_conf");
    }
}

static void done_gmx_groups_t(gmx_groups_t *g)
{
    int i;

    for (i = 0; (i < egcNR); i++)
    {
        if (NULL != g->grps[i].nm_ind)
        {
            sfree(g->grps[i].nm_ind);
            g->grps[i].nm_ind = NULL;
        }
        if (NULL != g->grpnr[i])
        {
            sfree(g->grpnr[i]);
            g->grpnr[i] = NULL;
        }
    }
    /* The contents of this array is in symtab, don't free it here */
    sfree(g->grpname);
}

gmx_bool read_tps_conf(const char *infile, char *title, t_topology *top, int *ePBC,
                       rvec **x, rvec **v, matrix box, gmx_bool bMass)
{
    t_tpxheader      header;
    int              natoms, i, version, generation;
    gmx_bool         bTop, bXNULL = FALSE;
    gmx_mtop_t      *mtop;
    gmx_atomprop_t   aps;

    bTop  = fn2bTPX(infile);
    *ePBC = -1;
    if (bTop)
    {
        read_tpxheader(infile, &header, TRUE, &version, &generation);
        if (x)
        {
            snew(*x, header.natoms);
        }
        if (v)
        {
            snew(*v, header.natoms);
        }
        snew(mtop, 1);
        *ePBC = read_tpx(infile, NULL, box, &natoms,
                         (x == NULL) ? NULL : *x, (v == NULL) ? NULL : *v, NULL, mtop);
        *top = gmx_mtop_t_to_t_topology(mtop);
        /* In this case we need to throw away the group data too */
        done_gmx_groups_t(&mtop->groups);
        sfree(mtop);
        std::strcpy(title, *top->name);
        tpx_make_chain_identifiers(&top->atoms, &top->mols);
    }
    else
    {
        get_stx_coordnum(infile, &natoms);
        init_t_atoms(&top->atoms, natoms, (fn2ftp(infile) == efPDB));
        if (x == NULL)
        {
            snew(x, 1);
            bXNULL = TRUE;
        }
        snew(*x, natoms);
        if (v)
        {
            snew(*v, natoms);
        }
        read_stx_conf(infile, title, &top->atoms, *x, (v == NULL) ? NULL : *v, ePBC, box);
        if (bXNULL)
        {
            sfree(*x);
            sfree(x);
        }
        if (bMass)
        {
            aps = gmx_atomprop_init();
            for (i = 0; (i < natoms); i++)
            {
                if (!gmx_atomprop_query(aps, epropMass,
                                        *top->atoms.resinfo[top->atoms.atom[i].resind].name,
                                        *top->atoms.atomname[i],
                                        &(top->atoms.atom[i].m)))
                {
                    if (debug)
                    {
                        fprintf(debug, "Can not find mass for atom %s %d %s, setting to 1\n",
                                *top->atoms.resinfo[top->atoms.atom[i].resind].name,
                                top->atoms.resinfo[top->atoms.atom[i].resind].nr,
                                *top->atoms.atomname[i]);
                    }
                }
            }
            gmx_atomprop_destroy(aps);
        }
        top->idef.ntypes = -1;
    }

    return bTop;
}

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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

#include "groio.h"

#include <cstdio>
#include <cstring>

#include <algorithm>
#include <string>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static void get_coordnum_fp(FILE *in, char *title, int *natoms)
{
    char line[STRLEN+1];

    fgets2(title, STRLEN, in);
    fgets2(line, STRLEN, in);
    if (sscanf(line, "%d", natoms) != 1)
    {
        gmx_fatal(FARGS, "gro file does not have the number of atoms on the second line");
    }
}

void get_coordnum(const char *infile, int *natoms)
{
    FILE *in;
    char  title[STRLEN];

    in = gmx_fio_fopen(infile, "r");
    get_coordnum_fp(in, title, natoms);
    gmx_fio_fclose(in);
}

/* Note that the .gro reading routine still support variable precision
 * for backward compatibility with old .gro files.
 * We have removed writing of variable precision to avoid compatibility
 * issues with other software packages.
 */
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
    gmx_bool   bFirst, bVel, oldResFirst;
    char      *p1, *p2, *p3;

    oldres      = -1;
    newres      = -1;
    oldResFirst = FALSE;
    ddist       = 0;

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
        if ((fgets2(line, STRLEN, in)) == NULL)
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

        if (!oldResFirst || oldres != resnr || strncmp(resname, oldresname, sizeof(resname)))
        {
            oldres      = resnr;
            oldResFirst = TRUE;
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
            if (sscanf(buf, "%lf %lf", &x1, &x2) != 1)
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
                if (sscanf(buf, "%lf", &x1) != 1)
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
    fgets2(line, STRLEN, in);
    if (sscanf(line, "%lf%lf%lf", &x1, &y1, &z1) != 3)
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

void gmx_gro_read_conf(const char *infile,
                       t_topology *top, rvec x[], rvec *v, matrix box)
{
    FILE *in = gmx_fio_fopen(infile, "r");
    int   ndec;
    char  title[STRLEN];
    get_w_conf(in, infile, title, &top->symtab, &top->atoms, &ndec, x, v, box);
    top->name = put_symtab(&top->symtab, title);
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

static const char *get_hconf_format(bool haveVelocities)
{
    if (haveVelocities)
    {
        return "%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n";
    }
    else
    {
        return "%8.3f%8.3f%8.3f\n";
    }

}

static void write_hconf_box(FILE *out, const matrix box)
{
    if (box[XX][YY] || box[XX][ZZ] || box[YY][XX] || box[YY][ZZ] ||
        box[ZZ][XX] || box[ZZ][YY])
    {
        fprintf(out, "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",
                box[XX][XX], box[YY][YY], box[ZZ][ZZ],
                box[XX][YY], box[XX][ZZ], box[YY][XX],
                box[YY][ZZ], box[ZZ][XX], box[ZZ][YY]);
    }
    else
    {
        fprintf(out, "%10.5f%10.5f%10.5f\n",
                box[XX][XX], box[YY][YY], box[ZZ][ZZ]);
    }
}

void write_hconf_indexed_p(FILE *out, const char *title, const t_atoms *atoms,
                           int nx, const int index[],
                           const rvec *x, const rvec *v, const matrix box)
{
    char resnm[6], nm[6];
    int  ai, i, resind, resnr;

    fprintf(out, "%s\n", (title && title[0]) ? title : gmx::bromacs().c_str());
    fprintf(out, "%5d\n", nx);

    const char *format = get_hconf_format(v != NULL);

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

    write_hconf_box(out, box);

    fflush(out);
}

void write_hconf_mtop(FILE *out, const char *title, gmx_mtop_t *mtop,
                      const rvec *x, const rvec *v, const matrix box)
{
    int                     i, resnr;
    gmx_mtop_atomloop_all_t aloop;
    t_atom                 *atom;
    char                   *atomname, *resname;

    fprintf(out, "%s\n", (title && title[0]) ? title : gmx::bromacs().c_str());
    fprintf(out, "%5d\n", mtop->natoms);

    const char *format = get_hconf_format(v != NULL);

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

    write_hconf_box(out, box);

    fflush(out);
}

void write_hconf_p(FILE *out, const char *title, const t_atoms *atoms,
                   const rvec *x, const rvec *v, const matrix box)
{
    int     *aa;
    int      i;

    snew(aa, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        aa[i] = i;
    }
    write_hconf_indexed_p(out, title, atoms, atoms->nr, aa, x, v, box);
    sfree(aa);
}

void write_conf_p(const char *outfile, const char *title,
                  const t_atoms *atoms,
                  const rvec *x, const rvec *v, const matrix box)
{
    FILE *out;

    out = gmx_fio_fopen(outfile, "w");
    write_hconf_p(out, title, atoms, x, v, box);
    gmx_fio_fclose(out);
}

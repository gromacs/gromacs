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

#include "g96io.h"

#include <cinttypes>
#include <cstdio>
#include <cstring>

#include <string>

#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#define CHAR_SHIFT 24

static int read_g96_pos(char line[], t_symtab* symtab, FILE* fp, const std::filesystem::path& infile, t_trxframe* fr)
{
    t_atoms* atoms;
    gmx_bool bEnd;
    int      nwanted, natoms, atnr, resnr = 0, oldres, newres, shift;
    char     anm[STRLEN], resnm[STRLEN];
    char     c1, c2;
    double   db1, db2, db3;

    nwanted = fr->natoms;

    if (fr->atoms != nullptr)
    {
        GMX_RELEASE_ASSERT(symtab != nullptr,
                           "Reading a conformation from a g96 format with atom data requires a "
                           "valid symbol table");
    }
    atoms = fr->atoms;
    if (atoms != nullptr)
    {
        atoms->haveMass    = FALSE;
        atoms->haveCharge  = FALSE;
        atoms->haveType    = FALSE;
        atoms->haveBState  = FALSE;
        atoms->havePdbInfo = FALSE;
    }

    natoms = 0;

    if (fr->bX)
    {
        if (fr->bAtoms)
        {
            shift = CHAR_SHIFT;
        }
        else
        {
            shift = 0;
        }
        newres = -1;
        oldres = -666; /* Unlikely number for the first residue! */
        bEnd   = FALSE;
        while (!bEnd && fgets2(line, STRLEN, fp))
        {
            bEnd = (std::strncmp(line, "END", 3) == 0);
            if (!bEnd && (line[0] != '#'))
            {
                if (sscanf(line + shift, "%15lf%15lf%15lf", &db1, &db2, &db3) != 3)
                {
                    gmx_fatal(FARGS,
                              "Did not find 3 coordinates for atom %d in %s\n",
                              natoms + 1,
                              infile.string().c_str());
                }
                if ((nwanted != -1) && (natoms >= nwanted))
                {
                    gmx_fatal(FARGS,
                              "Found more coordinates (%d) in %s than expected %d\n",
                              natoms,
                              infile.string().c_str(),
                              nwanted);
                }
                if (atoms)
                {
                    if (fr->bAtoms
                        && (sscanf(line, "%5d%c%5s%c%5s%7d", &resnr, &c1, resnm, &c2, anm, &atnr) != 6))
                    {
                        if (oldres >= 0)
                        {
                            resnr = oldres;
                        }
                        else
                        {
                            resnr = 1;
                            strncpy(resnm, "???", sizeof(resnm) - 1);
                        }
                        strncpy(anm, "???", sizeof(anm) - 1);
                    }
                    atoms->atomname[natoms] = put_symtab(symtab, anm);
                    if (resnr != oldres)
                    {
                        oldres = resnr;
                        newres++;
                        if (newres >= atoms->nr)
                        {
                            gmx_fatal(FARGS,
                                      "More residues than atoms in %s (natoms = %d)",
                                      infile.string().c_str(),
                                      atoms->nr);
                        }
                        atoms->atom[natoms].resind = newres;
                        if (newres + 1 > atoms->nres)
                        {
                            atoms->nres = newres + 1;
                        }
                        t_atoms_set_resinfo(atoms, natoms, symtab, resnm, resnr, ' ', 0, ' ');
                    }
                    else
                    {
                        atoms->atom[natoms].resind = newres;
                    }
                }
                if (fr->x)
                {
                    fr->x[natoms][0] = db1;
                    fr->x[natoms][1] = db2;
                    fr->x[natoms][2] = db3;
                }
                natoms++;
            }
        }
        if ((nwanted != -1) && natoms != nwanted)
        {
            fprintf(stderr,
                    "Warning: found less coordinates (%d) in %s than expected %d\n",
                    natoms,
                    infile.string().c_str(),
                    nwanted);
        }
    }

    fr->natoms = natoms;

    return natoms;
}

static int read_g96_vel(char line[], FILE* fp, const std::filesystem::path& infile, t_trxframe* fr)
{
    gmx_bool bEnd;
    int      nwanted, natoms = -1, shift;
    double   db1, db2, db3;

    nwanted = fr->natoms;

    if (fr->v && fr->bV)
    {
        if (strcmp(line, "VELOCITYRED") == 0)
        {
            shift = 0;
        }
        else
        {
            shift = CHAR_SHIFT;
        }
        natoms = 0;
        bEnd   = FALSE;
        while (!bEnd && fgets2(line, STRLEN, fp))
        {
            bEnd = (strncmp(line, "END", 3) == 0);
            if (!bEnd && (line[0] != '#'))
            {
                if (sscanf(line + shift, "%15lf%15lf%15lf", &db1, &db2, &db3) != 3)
                {
                    gmx_fatal(FARGS,
                              "Did not find 3 velocities for atom %d in %s",
                              natoms + 1,
                              infile.string().c_str());
                }
                if ((nwanted != -1) && (natoms >= nwanted))
                {
                    gmx_fatal(FARGS,
                              "Found more velocities (%d) in %s than expected %d\n",
                              natoms,
                              infile.string().c_str(),
                              nwanted);
                }
                if (fr->v)
                {
                    fr->v[natoms][0] = db1;
                    fr->v[natoms][1] = db2;
                    fr->v[natoms][2] = db3;
                }
                natoms++;
            }
        }
        if ((nwanted != -1) && (natoms != nwanted))
        {
            fprintf(stderr,
                    "Warning: found less velocities (%d) in %s than expected %d\n",
                    natoms,
                    infile.string().c_str(),
                    nwanted);
        }
    }

    return natoms;
}

int read_g96_conf(FILE* fp, const std::filesystem::path& infile, char** name, t_trxframe* fr, t_symtab* symtab, char* line)
{
    gmx_bool bAtStart, bTime, bAtoms, bPos, bVel, bBox, bEnd, bFinished;
    int      natoms, nbp;
    double   db1, db2, db3, db4, db5, db6, db7, db8, db9;

    bAtStart = (ftell(fp) == 0);

    clear_trxframe(fr, FALSE);

    natoms = 0;

    if (bAtStart)
    {
        bool foundTitle = false;
        while (!foundTitle && fgets2(line, STRLEN, fp))
        {
            foundTitle = (std::strcmp(line, "TITLE") == 0);
        }
        fgets2(line, STRLEN, fp);
        if (name != nullptr)
        {
            *name = gmx_strdup(line);
        }
        bEnd = FALSE;
        while (!bEnd && fgets2(line, STRLEN, fp))
        {
            bEnd = (std::strcmp(line, "END") == 0);
        }
        fgets2(line, STRLEN, fp);
    }

    /* Do not get a line if we are not at the start of the file, *
     * because without a parameter file we don't know what is in *
     * the trajectory and we have already read the line in the   *
     * previous call (VERY DIRTY).                               */
    bFinished = FALSE;
    do
    {
        bTime  = (std::strcmp(line, "TIMESTEP") == 0);
        bAtoms = (std::strcmp(line, "POSITION") == 0);
        bPos   = (bAtoms || (strcmp(line, "POSITIONRED") == 0));
        bVel   = (std::strncmp(line, "VELOCITY", 8) == 0);
        bBox   = (std::strcmp(line, "BOX") == 0);
        if (bTime)
        {
            if (!fr->bTime && !fr->bX)
            {
                fr->bStep = bTime;
                fr->bTime = bTime;
                do
                {
                    bFinished = (fgets2(line, STRLEN, fp) == nullptr);
                } while (!bFinished && (line[0] == '#'));
                sscanf(line, "%15" SCNd64 "%15lf", &(fr->step), &db1);
                fr->time = db1;
            }
            else
            {
                bFinished = TRUE;
            }
        }
        if (bPos)
        {
            if (!fr->bX)
            {
                fr->bAtoms = bAtoms;
                fr->bX     = bPos;
                natoms     = read_g96_pos(line, symtab, fp, infile, fr);
            }
            else
            {
                bFinished = TRUE;
            }
        }
        if (fr->v && bVel)
        {
            fr->bV = bVel;
            natoms = read_g96_vel(line, fp, infile, fr);
        }
        if (bBox)
        {
            fr->bBox = bBox;
            clear_mat(fr->box);
            bEnd = FALSE;
            while (!bEnd && fgets2(line, STRLEN, fp))
            {
                bEnd = (strncmp(line, "END", 3) == 0);
                if (!bEnd && (line[0] != '#'))
                {
                    nbp = sscanf(line,
                                 "%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf",
                                 &db1,
                                 &db2,
                                 &db3,
                                 &db4,
                                 &db5,
                                 &db6,
                                 &db7,
                                 &db8,
                                 &db9);
                    if (nbp < 3)
                    {
                        gmx_fatal(FARGS, "Found a BOX line, but no box in %s", infile.string().c_str());
                    }
                    fr->box[XX][XX] = db1;
                    fr->box[YY][YY] = db2;
                    fr->box[ZZ][ZZ] = db3;
                    if (nbp == 9)
                    {
                        fr->box[XX][YY] = db4;
                        fr->box[XX][ZZ] = db5;
                        fr->box[YY][XX] = db6;
                        fr->box[YY][ZZ] = db7;
                        fr->box[ZZ][XX] = db8;
                        fr->box[ZZ][YY] = db9;
                    }
                }
            }
            bFinished = TRUE;
        }
    } while (!bFinished && (fgets2(line, STRLEN, fp) != nullptr));

    fr->natoms = natoms;

    return natoms;
}

void write_g96_conf(FILE* out, const char* title, const t_trxframe* fr, int nindex, const int* index)
{
    t_atoms* atoms;
    int      nout, i, a;

    atoms = fr->atoms;

    if (index)
    {
        nout = nindex;
    }
    else
    {
        nout = fr->natoms;
    }

    fprintf(out, "TITLE\n%s\nEND\n", title);
    if (fr->bStep || fr->bTime)
    {
        /* Officially the time format is %15.9, which is not enough for 10 ns */
        fprintf(out, "TIMESTEP\n%15" PRId64 "%15.6f\nEND\n", fr->step, fr->time);
    }
    if (fr->bX)
    {
        if (fr->bAtoms)
        {
            fprintf(out, "POSITION\n");
            for (i = 0; i < nout; i++)
            {
                if (index)
                {
                    a = index[i];
                }
                else
                {
                    a = i;
                }
                std::string residueName = *atoms->resinfo[atoms->atom[a].resind].name;
                std::string atomName    = *atoms->atomname[a];
                // g96 is fixed format, so we need to limit strings to max 5
                // characters.
                residueName.resize(5);
                atomName.resize(5);
                fprintf(out,
                        "%5d %-5s %-5s%7d%15.9f%15.9f%15.9f\n",
                        (atoms->resinfo[atoms->atom[a].resind].nr) % 100000,
                        residueName.c_str(),
                        atomName.c_str(),
                        (i + 1) % 10000000,
                        fr->x[a][XX],
                        fr->x[a][YY],
                        fr->x[a][ZZ]);
            }
        }
        else
        {
            fprintf(out, "POSITIONRED\n");
            for (i = 0; i < nout; i++)
            {
                if (index)
                {
                    a = index[i];
                }
                else
                {
                    a = i;
                }
                fprintf(out, "%15.9f%15.9f%15.9f\n", fr->x[a][XX], fr->x[a][YY], fr->x[a][ZZ]);
            }
        }
        fprintf(out, "END\n");
    }
    if (fr->bV)
    {
        if (fr->bAtoms)
        {
            fprintf(out, "VELOCITY\n");
            for (i = 0; i < nout; i++)
            {
                if (index)
                {
                    a = index[i];
                }
                else
                {
                    a = i;
                }
                fprintf(out,
                        "%5d %-5s %-5s%7d%15.9f%15.9f%15.9f\n",
                        (atoms->resinfo[atoms->atom[a].resind].nr) % 100000,
                        *atoms->resinfo[atoms->atom[a].resind].name,
                        *atoms->atomname[a],
                        (i + 1) % 10000000,
                        fr->v[a][XX],
                        fr->v[a][YY],
                        fr->v[a][ZZ]);
            }
        }
        else
        {
            fprintf(out, "VELOCITYRED\n");
            for (i = 0; i < nout; i++)
            {
                if (index)
                {
                    a = index[i];
                }
                else
                {
                    a = i;
                }
                fprintf(out, "%15.9f%15.9f%15.9f\n", fr->v[a][XX], fr->v[a][YY], fr->v[a][ZZ]);
            }
        }
        fprintf(out, "END\n");
    }
    if (fr->bBox)
    {
        fprintf(out, "BOX\n");
        fprintf(out, "%15.9f%15.9f%15.9f", fr->box[XX][XX], fr->box[YY][YY], fr->box[ZZ][ZZ]);
        if ((fr->box[XX][YY] != 0.0F) || (fr->box[XX][ZZ] != 0.0F) || (fr->box[YY][XX] != 0.0F)
            || (fr->box[YY][ZZ] != 0.0F) || (fr->box[ZZ][XX] != 0.0F) || (fr->box[ZZ][YY] != 0.0F))
        {
            fprintf(out,
                    "%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f",
                    fr->box[XX][YY],
                    fr->box[XX][ZZ],
                    fr->box[YY][XX],
                    fr->box[YY][ZZ],
                    fr->box[ZZ][XX],
                    fr->box[ZZ][YY]);
        }
        fprintf(out, "\n");
        fprintf(out, "END\n");
    }
}

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

#ifndef _fgrid_h
#define _fgrid_h

#include <stdio.h>

#include "xdlgitem.h"

typedef struct
{
    edlgitem edlg;
    bool     bDef;
    int      nname;
    char**   name;
    char *   set, *get, *def, *help;
} t_fitem;

typedef struct
{
    char*     name;
    int       x, y, w, h;
    int       nfitem;
    t_fitem** fitem;
} t_fgroup;

typedef struct
{
    int      x, y, w, h;
    t_fitem* fitem;
} t_fsimple;

typedef struct
{
    int         w, h;
    int         nfgroup;
    t_fgroup**  fgroup;
    int         nfsimple;
    t_fsimple** fsimple;
} t_fgrid;

typedef enum
{
    eGRIDEXP,
    eACCOEXP,
    eACCCEXP,
    eGRPEXP,
    eITEMEXP,
    eSAMEPOINT,
    eTOOWIDE,
    eTOOHIGH,
    eQUOTE,
    eNOVALS
} eDLGERR;

void ReadDlgErr(const char* infile, eDLGERR err, const char* s);

t_fgrid* FGridFromFile(const char* infile);

void DoneFGrid(t_fgrid* fgrid);

void DumpFGrid(t_fgrid* fgrid);

void ReadQuoteString(const char* infile, FILE* in, char* buf);

#endif /* _fgrid_h */

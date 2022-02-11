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

#ifndef _xdlghi_h
#define _xdlghi_h

#include <stdarg.h>

#include "Xstuff.h"
#include "x11.h"
#include "xdlg.h"

typedef struct
{
    int         nitem;
    int         w, h;
    t_dlgitem** list;
} t_dlgitemlist;

extern t_dlgitem** CreateRadioButtonGroup(t_x11* x11,
                                          char*  szTitle,
                                          t_id   GroupID,
                                          int    nrb,
                                          t_id   rb[],
                                          int    nSelect,
                                          char*  szRB[],
                                          int    x0,
                                          int    y0);
/* This routine creates a radio button group at the
 * specified position. The return values is a pointer to an
 * array of dlgitems, the array has length (nrb+1) with the +1
 * because of the groupbox.
 * nSelect is the ordinal of the selected button.
 */

extern t_dlgitem**
CreateDlgitemGroup(t_x11* x11, const char* szTitle, t_id GroupID, int x0, int y0, int nitem, ...);
/* This routine creates a dlgitem group at the
 * specified position. The return values is a pointer to an
 * array of dlgitems, the array has length (nitem+1) with the +1
 * because of the groupbox.
 */

extern t_dlg* ReadDlg(t_x11*       x11,
                      Window       Parent,
                      const char*  title,
                      const char*  infile,
                      int          x0,
                      int          y0,
                      bool         bAutoPosition,
                      bool         bUseMon,
                      DlgCallback* cb,
                      void*        data);
/* Read a dialog box from a template file */

#endif /* _xdlghi_h */

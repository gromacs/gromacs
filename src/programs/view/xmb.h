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

#ifndef _xmb_h
#define _xmb_h

#include "manager.h"
#include "x11.h"
#include "xdlg.h"
#include "xmb.h"

#define MB_OK 1
#define MB_CANCEL (1 << 1)
#define MB_OKCANCEL (MB_OK | MB_CANCEL)
#define MB_YES (1 << 2)
#define MB_NO (1 << 3)
#define MB_YESNO (MB_YES | MB_NO)
#define MB_ICONSTOP (1 << 16)
#define MB_ICONINFORMATION (1 << 17)
#define MB_ICONEXCLAMATION (1 << 18)
#define MB_ICONGMX (1 << 19)
#define MB_SYSTEMMODAL (1 << 20)
#define MB_APPLMODAL (1 << 21)
#define MB_DONTSHOW (1 << 22)

t_dlg* MessageBox(t_x11*             x11,
                  Window             Parent,
                  const char*        title,
                  int                nlines,
                  const char* const* lines,
                  unsigned long      Flags,
                  DlgCallback*       cb,
                  void*              data);

#endif /* _xmb_h */

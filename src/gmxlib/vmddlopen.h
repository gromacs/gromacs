/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * This file is part of Gromacs        Copyright (c) 1991-2008
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

/***************************************************************************
 * cr
 * cr            (C) Copyright 1995-2009 The Board of Trustees of the
 * cr                        University of Illinois
 * cr                         All Rights Reserved
 * cr
   Developed by:           Theoretical and Computational Biophysics Group
                        University of Illinois at Urbana-Champaign
                        http://www.ks.uiuc.edu/

   Permission is hereby granted, free of charge, to any person obtaining a copy of
   this software and associated documentation files (the Software), to deal with
   the Software without restriction, including without limitation the rights to
   use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
   of the Software, and to permit persons to whom the Software is furnished to
   do so, subject to the following conditions:

   Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimers.

   Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimers in the documentation
   and/or other materials provided with the distribution.

   Neither the names of Theoretical and Computational Biophysics Group,
   University of Illinois at Urbana-Champaign, nor the names of its contributors
   may be used to endorse or promote products derived from this Software without
   specific prior written permission.

   THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
   THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
   OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
   OTHER DEALINGS WITH THE SOFTWARE.
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vmddlopen.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.9 $      $Date: 2009/07/07 02:40:05 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Routines for loading dynamic link libraries and shared object files
 *   on various platforms, abstracting from machine dependent APIs.
 *
 * LICENSE:
 *   UIUC Open Source License
 *   http://www.ks.uiuc.edu/Research/vmd/plugins/pluginlicense.html
 *
 ***************************************************************************/

/*
 * vmddlopen: thin multi-platform wrapper around dlopen/LoadLibrary
 */

#ifndef VMD_DLOPEN__

#ifdef __cplusplus
extern "C" {
#endif

/* Try to open the specified library.  All symbols must be resolved or the
 * load will fail (RTLD_NOW).
 */
void *vmddlopen(const char *fname);

/* Try to load the specified symbol using the given handle.  Returns NULL if
 * the symbol cannot be loaded.
 */
void *vmddlsym(void *h, const char *sym);

/* Unload the library.  Return 0 on success, nonzero on error.
 */
int vmddlclose(void *h);

/* Return last error from any of the above functions.  Not thread-safe on
 * Windows due to static buffer in our code.
 */
const char *vmddlerror(void);

#ifdef __cplusplus
}
#endif

#endif

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
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

#ifndef _copyrite_h
#define _copyrite_h
#include "visibility.h"

#include <stdio.h>
#include "types/simple.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Has to be a function, so we can get version number from autoconf */
GMX_LIBGMX_EXPORT
const char *GromacsVersion(void);

void
gmx_print_version_info(FILE *fp);


static const char *
CopyrightText[] = {
    "Contributions from Mark Abraham, Emile Apol, Rossen Apostolov, ",
    "Herman J.C. Berendsen, Aldert van Buuren, Pär Bjelkmar,  ",
    "Rudi van Drunen, Anton Feenstra, Gerrit Groenhof, Christoph Junghans, ",
    "Peter Kasson, Carsten Kutzner, Per Larsson, Pieter Meulenhoff, ",
    "Teemu Murtola, Szilard Pall, Sander Pronk, Roland Schulz, ",
    "Michael Shirts, Alfons Sijbers, Peter Tieleman,\n",
    "Berk Hess, David van der Spoel, and Erik Lindahl.\n",
    "Copyright (c) 1991-2000, University of Groningen, The Netherlands.",
    "Copyright (c) 2001-2012,2013, The GROMACS development team at",
    "Uppsala University & The Royal Institute of Technology, Sweden.",
    "check out http://www.gromacs.org for more information.\n"
};

static const char *
LicenseText[] = {
    "This program is free software; you can redistribute it and/or",
    "modify it under the terms of the GNU Lesser General Public License",
    "as published by the Free Software Foundation; either version 2.1",
    "of the License, or (at your option) any later version."
};


void
pr_difftime(FILE *out, double dt);

GMX_LIBGMX_EXPORT
void
CopyRight(FILE *out, const char *szProgram);


/* For both bromacs() and cool_quote() you have to provide a pointer to
 * a string of reasonable length (say 256) and the string length. This
 * is necessary to make the routines threadsafe and avoid allocating
 * a new string each time. The retstring pointer will be the return value.
 */
GMX_LIBGMX_EXPORT
void
bromacs(char *retstring, int retsize);

/* For cool_quote, the number of the quote used will be returned in cqnum
 * if it is non-NULL.
 */
GMX_LIBGMX_EXPORT
void
cool_quote(char *retstring, int retsize, int *cqnum);

gmx_bool
be_cool(void);
/* Return TRUE when the user is COOL, FALSE otherwise */

GMX_LIBGMX_EXPORT
void
thanx(FILE *fp);

enum {
    eCITEGMX, eCITEBATH, eCITESHAKE, eCITESETTLE, eCITESOR,
    eCITEDISRE, eCITERF, eCITELINCS, eCITENR
};

GMX_LIBGMX_EXPORT
void
please_cite(FILE *fp, const char *key);
/* Print a message asking to cite something... */

#ifdef __cplusplus
}
#endif

#endif  /* _copyright_h */

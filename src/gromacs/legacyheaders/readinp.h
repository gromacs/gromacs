/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#ifndef _readinp_h
#define _readinp_h

#include <string.h>

#include "gromacs/legacyheaders/warninp.h"
#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int      count;     /* sort order for output  */
    gmx_bool bObsolete; /* whether it is an obsolete param value */
    gmx_bool bSet;      /* whether it it has been read out */
    char    *name;      /* name of the parameter */
    char    *value;     /* parameter value string */
    int      inp_count; /* number of einps read. Only valid for the first item
                                                 in the inpfile list. */
} t_inpfile;
/* entry in input files (like .mdp files).
   Initally read in with read_inpfile, then filled in with missing values
   through get_eint, get_ereal, etc. */

t_inpfile *read_inpfile(const char *fn, int *ninp,
                        warninp_t wi);
/* Create & populate a t_inpfile struct from values in file fn.
   fn = the file name
   ninp = the number of read parameters
   cppopts = the cpp-style options for #include paths and #defines */

void write_inpfile(const char *fn, int ninp, t_inpfile inp[],
                   gmx_bool bHaltOnUnknown,
                   warninp_t wi);

void replace_inp_entry(int ninp, t_inpfile *inp,
                       const char *old_entry, const char *new_entry);

int search_einp(int ninp, const t_inpfile *inp, const char *name);
/* Return the index of an .mdp field with the given name within the
 * inp array, if it exists. Return -1 if it does not exist. */

int get_eint(int *ninp, t_inpfile **inp, const char *name, int def,
             warninp_t wi);

gmx_int64_t get_eint64(int *ninp, t_inpfile **inp,
                       const char *name, gmx_int64_t def,
                       warninp_t);

double get_ereal(int *ninp, t_inpfile **inp, const char *name, double def,
                 warninp_t wi);

const char *get_estr(int *ninp, t_inpfile **inp, const char *name, const char *def);

int get_eeenum(int *ninp, t_inpfile **inp, const char *name, const char **defs,
               warninp_t wi);
/* defs must be NULL terminated */

int get_eenum(int *ninp, t_inpfile **inp, const char *name, const char **defs);
/* defs must be NULL terminated */

/* Here are some dirty macros to extract data from the inp structures.
 * Most macros assume the variables ninp, inp and wi are present.
 * Elements are removed from the list after reading.
 */
#define REM_TYPE(name)       replace_inp_entry(ninp, inp, name, NULL)
#define REPL_TYPE(old, new)   replace_inp_entry(ninp, inp, old, new)
#define STYPE(name, var, def)  if ((tmp = get_estr(&ninp, &inp, name, def)) != NULL) strcpy(var, tmp)
#define STYPENC(name, def) get_estr(&ninp, &inp, name, def)
#define ITYPE(name, var, def)  var    = get_eint(&ninp, &inp, name, def, wi)
#define STEPTYPE(name, var, def)  var = get_eint64(&ninp, &inp, name, def, wi)
#define RTYPE(name, var, def)  var    = get_ereal(&ninp, &inp, name, def, wi)
#define ETYPE(name, var, defs) var    = get_eenum(&ninp, &inp, name, defs)
#define EETYPE(name, var, defs) var   = get_eeenum(&ninp, &inp, name, defs, wi)
#define CCTYPE(s) STYPENC("\n; " s, NULL)
#define CTYPE(s)  STYPENC("; " s, NULL)
/* This last one prints a comment line where you can add some explanation */

#ifdef __cplusplus
}
#endif

#endif

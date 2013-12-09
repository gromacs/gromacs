/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013, by the GROMACS development team, led by
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
#ifndef GMX_COMMANDLINE_PARGS_H
#define GMX_COMMANDLINE_PARGS_H

#include "../legacyheaders/types/simple.h"
#include "../legacyheaders/oenv.h"
#include "../fileio/filenm.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* This structure is used for parsing arguments off the comand line */
enum
{
    etINT, etGMX_INT64, etREAL, etTIME, etSTR,    etBOOL, etRVEC,   etENUM, etNR
};

typedef struct
{
    const char *option;
    gmx_bool    bSet;
    int         type;
    union
    {
        void            *v; /* This is a nasty workaround, to be able to use initialized */
        int             *i; /* arrays */
        gmx_int64_t     *is;
        real            *r;
        const char     **c; /* Must be pointer to string (when type == etSTR)         */
        /* or null terminated list of enums (when type == etENUM) */
        gmx_bool        *b;
        rvec            *rv;
    }           u;
    const char *desc;
} t_pargs;

gmx_bool is_hidden(t_pargs *pa);
/* Return TRUE when the option is a secret one */

int nenum(const char *const enumc[]);
/* returns ordinal number of selected enum from args
 * depends on enumc[0] pointing to one of the other elements
 * array must be terminated by a NULL pointer
 */

int opt2parg_int(const char *option, int nparg, t_pargs pa[]);

gmx_bool opt2parg_gmx_bool(const char *option, int nparg, t_pargs pa[]);

real opt2parg_real(const char *option, int nparg, t_pargs pa[]);

const char *opt2parg_str(const char *option, int nparg, t_pargs pa[]);

const char *opt2parg_enum(const char *option, int nparg, t_pargs pa[]);

gmx_bool opt2parg_bSet(const char *option, int nparg, t_pargs pa[]);


#define PCA_CAN_VIEW       (1<<5)
/* add option -w to view output files (must be implemented in program) */
#define PCA_CAN_BEGIN      (1<<6)
#define PCA_CAN_END        (1<<7)
#define PCA_CAN_DT         (1<<14)
#define PCA_CAN_TIME       (PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_DT)
/* adds options -b and -e for begin and end time for reading trajectories */
#define PCA_TIME_UNIT      (1<<15)
/* set time unit for output */
#define PCA_KEEP_ARGS      (1<<8)
/* keep parsed args in argv (doesn't make sense without NOEXIT_ON_ARGS) */
#define PCA_CAN_SET_DEFFNM (1<<10)
/* does something for non-master mdrun nodes */
#define PCA_NOEXIT_ON_ARGS (1<<11)
/* no fatal_error when invalid options are encountered */
#define PCA_QUIET          (1<<12)
/* does something for non-master mdrun nodes */
#define PCA_BE_NICE        (1<<13)
/* Default to low priority, unless configured with --disable-nice */
#define PCA_NOT_READ_NODE  (1<<16)
/* Is this node not reading: for parallel all nodes but the master */

gmx_bool parse_common_args(int *argc, char *argv[], unsigned long Flags,
                           int nfile, t_filenm fnm[], int npargs, t_pargs *pa,
                           int ndesc, const char **desc,
                           int nbugs, const char **bugs,
                           output_env_t *oenv);
/* Get arguments from the arg-list. The arguments extracted
 * are removed from the list. If manual is NULL a default message is displayed
 * when errors are encountered. The Flags argument, when non-0 enables
 * some input checks. Using this routine also means that the arguments
 * -b and -e will be used for begin and end time, whether this is
 * appropriate or not!
 */

#ifdef __cplusplus
}
#endif

#endif

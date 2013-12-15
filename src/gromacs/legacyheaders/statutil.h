/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _statutil_h
#define _statutil_h

#include "../fileio/filenm.h"
#include "readinp.h"
#include "oenv.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0 /* avoid screwing up indentation */
}
#endif


/* The code below is to facilitate controlled begin and end of
   trajectory reading. Corresponding routines in
   src/gmxlib/tcontrol.c
 */
enum {
    TBEGIN, TEND, TDELTA, TNR
};

gmx_bool bTimeSet(int tcontrol);

real rTimeValue(int tcontrol);

void setTimeValue(int tcontrol, real value);

/* End trajectory time control */

/* LEGACY FUNCTIONS

   The program names, command lines, etc. are now also set in the output_env
   structure. That is now the preferred location, but the functions here
   are still available as legacy functions. Because they all act on inherently
   global informaion, their existence in a multi-threaded environment is not
   a real problem. */

/* set the program name to the provided string, but note
 * that it must be a real file - we determine the library
 * directory from its location!
 */
const char *Program(void);
/* Id. without leading directory */
const char *ShortProgram(void);

/*****************************************************
 *         Some command line parsing routines
 *****************************************************/

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

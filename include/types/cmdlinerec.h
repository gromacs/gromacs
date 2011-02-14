/*  -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.5
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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifndef _cmdlinerec_h_
#define _cmdlinerec_h_

#include "simple.h"

#ifdef __cplusplus
extern "C" {
#endif

#if 0
}
#endif

/** \file include/types/cmdlinerec.h
 *
 *  \brief Structure for mdrun command-line arguments
 *
 *  This structure packages the command-line arguments to mdrun to
 *  enable them to be passed into the lower-level routines in a way
 *  that enhances future maintainability better than the old style of
 *  passing many different arguments to many different functions. In
 *  particular, maintaining several mdrun-work-alike utility programs
 *  was becoming a nightmare.
 *
 *  \param fnm                   Structure for managing input filenames
 *  \param nfile                 Size of the fnm structure
 *  \param oenv                  Controls the output environment
 *  \param nstglobalcomm         Number of steps between global communication when running in parallel
 *  \param dd_ncells             Number of DD cells in each direction
 *  \param dd_node_order         Controls the mapping of MPMD processes to hardware when using DD-PME
 *  \param dd_comm_distance_min  Maximum distance over which bonded interactions can occur when using DD
 *  \param rconstr_max           Maximum distance over which P-LINCS constraints can act when using DD
 *  \param dd_dlb_opt            Controls the use of dynamic load balancing with DD
 *  \param dd_dlb_scale          Minimum allowed scaling of the DD cell size with dynamic load balancing
 *  \param dd_cell_size          DD cell sizes for static load balancing
 *  \param nstepout              Number of steps between writing remaining runtime, etc.
 *  \param resetstep             Step at which to reset performance counters
 *  \param repl_ex_nst           Number of steps between replica-exchange attempts
 *  \param repl_ex_seed          Replica-exchange random-number generator seed
 *  \param pforce                Minimum force that will be dumped to stderr when debugging
 *  \param cpt_period            Minutes between writing of checkpoint file
 *  \param max_hours             Terminate after 0.99 of this time (hours)
 *  \param deviceOptions         Device option string
 */

typedef struct
{
    t_filenm *fnm;
    int nfile;
    output_env_t oenv;
    int nstglobalcomm;
    ivec dd_ncells;
    int dd_node_order;
    real dd_comm_distance_min;
    real rconstr_max;
    const char *dd_dlb_opt;
    real dd_dlb_scale;
    const char *dd_cell_size[DIM];
    int nstepout;
    int resetstep;
    int repl_ex_nst;
    int repl_ex_seed;
    real pforce;
    real cpt_period;
    real max_hours;
#ifdef GMX_OPENMM
    const char *deviceOptions;
#endif
} gmx_cmdlinerec;

/* We can't actually have an abstract type, but for consistency here's
 * a similar typedef */
typedef gmx_cmdlinerec *gmx_cmdlinerec_t;

#ifdef __cplusplus
}
#endif

#endif

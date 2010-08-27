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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "statutil.h"
#include "gmx_fatal.h"

#ifdef GMX_THREADS
#include "thread_mpi.h"
#endif

/* The source code in this file should be thread-safe. 
         Please keep it that way. */

/* Globals for trajectory input */
typedef struct {
  real t;
  gmx_bool bSet;
} t_timecontrol;

static t_timecontrol timecontrol[TNR] = {
  { 0, FALSE },
  { 0, FALSE },
  { 0, FALSE }
};

#ifdef GMX_THREADS
static tMPI_Thread_mutex_t tc_mutex=TMPI_THREAD_MUTEX_INITIALIZER;
#endif

typedef struct {
  real tfactor;
  const char *tstr,*xvgstr;
} t_timeconvert;

static const t_timeconvert timeconvert[] = {
    { 0,                   NULL,  NULL       },
    { 1e3,  		   "fs",  "fs"       },
    { 1,    		   "ps",  "ps"       },
    { 1e-3, 		   "ns",  "ns"       },
    { 1e-6, 		   "us",  "\\mus"    }, 
    { 1e-9, 		   "ms",  "ms"       },
    { 1e-12, 		   "s",   "s"        },
    { (1.0/60.0)*1e-12,    "m",   "m"        },
    { (1.0/3600.0)*1e-12,  "h",   "h"        },
    { 0,		   NULL,  NULL       }
};

gmx_bool bTimeSet(int tcontrol)
{
    gmx_bool ret;

#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&tc_mutex);
#endif
    range_check(tcontrol,0,TNR);
    ret=timecontrol[tcontrol].bSet;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&tc_mutex);
#endif

    return ret;
}
  
real rTimeValue(int tcontrol)
{
    real ret;

#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&tc_mutex);
#endif
    range_check(tcontrol,0,TNR);
    ret=timecontrol[tcontrol].t;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&tc_mutex);
#endif
    return ret;
}
  
void setTimeValue(int tcontrol,real value)
{
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&tc_mutex);
#endif
    range_check(tcontrol,0,TNR);
    timecontrol[tcontrol].t = value;
    timecontrol[tcontrol].bSet = TRUE;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&tc_mutex);
#endif
}



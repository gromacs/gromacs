/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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


#include <ctype.h>
#include "sysstuff.h"
#include "macros.h"
#include "string2.h"
#include "smalloc.h"
#include "pbc.h"
#include "statutil.h"
#include "names.h"
#include "vec.h"
#include "futil.h"
#include "wman.h"
#include "tpxio.h"
#include "gmx_fatal.h"
#include "network.h"
#include "vec.h"
#include "mtop_util.h"
#include "gmxfio.h"
#include "oenv.h"

#ifdef GMX_THREADS
#include "thread_mpi.h"
#endif



/* The source code in this file should be thread-safe. 
      Please keep it that way. */

/******************************************************************
 *
 *             T R A J E C T O R Y   S T U F F
 *
 ******************************************************************/

/* read only time names */
/* These must correspond to the time units type time_unit_t in statutil.h */
static const real timefactors[] =   { 0,  1e3,  1, 1e-3, 1e-6, 1e-9, 1e-12, 0 };
static const real timeinvfactors[] ={ 0, 1e-3,  1,  1e3,  1e6,  1e9,  1e12, 0 };
static const char *time_units_str[] = { NULL, "fs", "ps", "ns", "us", 
                                        "\\mus", "ms", "s" };
static const char *time_units_xvgr[] = { NULL, "fs", "ps", "ns",  
                                        "ms", "s", NULL };



/***** OUTPUT_ENV MEMBER FUNCTIONS ******/

void output_env_init(output_env_t oenv,  int argc, char *argv[],
                     time_unit_t tmu, gmx_bool view, xvg_format_t xvg_format,
                     int verbosity, int debug_level)
{
    int i;
    int cmdlength=0;
    char *argvzero=NULL;

    oenv->time_unit  = tmu;
    oenv->view=view;
    oenv->xvg_format = xvg_format;
    oenv->verbosity=verbosity;
    oenv->debug_level=debug_level;
    oenv->program_name=NULL;

    if (argv)
        argvzero=argv[0];

    /* set program name */
    /* When you run a dynamically linked program before installing
     * it, libtool uses wrapper scripts and prefixes the name with "lt-".
     * Until libtool is fixed to set argv[0] right, rip away the prefix:
     */
    if (argvzero)
    {
        if(strlen(argvzero)>3 && !strncmp(argvzero,"lt-",3))
            oenv->program_name=strdup(argvzero+3);
        else
            oenv->program_name=strdup(argvzero);
    }
    if (oenv->program_name == NULL)
        oenv->program_name = strdup("GROMACS");
   
    /* copy command line */ 
    if (argv) 
    {
        cmdlength = strlen(argvzero);
        for (i=1; i<argc; i++) 
        {
            cmdlength += strlen(argv[i]);
        }
    }
        
    /* Fill the cmdline string */
    snew(oenv->cmd_line,cmdlength+argc+1);
    for (i=0; i<argc; i++) 
    {
        strcat(oenv->cmd_line,argv[i]);
        strcat(oenv->cmd_line," ");
    }
}


void output_env_init_default(output_env_t oenv)
{
    output_env_init(oenv, 0, NULL, time_ps, FALSE, exvgNONE, 0, 0);
}


void output_env_done(output_env_t oenv)
{
    sfree(oenv->program_name);
    sfree(oenv->cmd_line);
    sfree(oenv);
}



int output_env_get_verbosity(const output_env_t oenv)
{
    return oenv->verbosity;
}

int output_env_get_debug_level(const output_env_t oenv)
{
    return oenv->debug_level;
}


const char *output_env_get_time_unit(const output_env_t oenv)
{
    return time_units_str[oenv->time_unit];
}

const char *output_env_get_time_label(const output_env_t oenv)
{
    char *label;
    snew(label, 20);
    
    sprintf(label,"Time (%s)",time_units_str[oenv->time_unit] ? 
            time_units_str[oenv->time_unit]: "ps");
    
    return label;
}

const char *output_env_get_xvgr_tlabel(const output_env_t oenv)
{
    char *label;
    snew(label, 20);
    
    sprintf(label,"Time (%s)", time_units_xvgr[oenv->time_unit] ?
            time_units_xvgr[oenv->time_unit] : "ps");
    
    return label;
}


real output_env_get_time_factor(const output_env_t oenv)
{
    return timefactors[oenv->time_unit];
}

real output_env_get_time_invfactor(const output_env_t oenv)
{
    return timeinvfactors[oenv->time_unit];
}

real output_env_conv_time(const output_env_t oenv, real time)
{
    return time*timefactors[oenv->time_unit];
}


void output_env_conv_times(const output_env_t oenv, int n, real *time)
{
    int i;
    double fact=timefactors[oenv->time_unit];
    
    if (fact!=1.)
        for(i=0; i<n; i++)
            time[i] *= fact;
}

gmx_bool output_env_get_view(const output_env_t oenv)
{
    return oenv->view;
}

xvg_format_t output_env_get_xvg_format(const output_env_t oenv)
{
    return oenv->xvg_format;
}

const char *output_env_get_program_name(const output_env_t oenv)
{
    return oenv->program_name;
}

const char *output_env_get_short_program_name(const output_env_t oenv)
{
    const char *pr,*ret;
    pr=ret=oenv->program_name; 
    if ((pr=strrchr(ret,DIR_SEPARATOR)) != NULL)
        ret=pr+1;
    /* Strip away the libtool prefix if it's still there. */
    if(strlen(ret) > 3 && !strncmp(ret, "lt-", 3))
        ret = ret + 3;
    return ret;
}



const char *output_env_get_cmd_line(const output_env_t oenv)
{
    return oenv->cmd_line;
}



/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: molprop_sqlite3.c,v 1.23 2009/06/01 06:13:18 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#ifdef HAVE_SQLITE3
#include <sqlite3.h>
#endif
#include "gmx_fatal.h"
#include "futil.h"
#include "macros.h"
#include "grompp.h"
#include "smalloc.h"
#include "string2.h"
#include "molprop_sqlite3.h"

void gmx_molprop_read_sqlite3(int np,gmx_molprop_t mp[],const char *sqlite_file)
{
    int rc;
#ifdef HAVE_SQLITE3
    sqlite3 *db;

    if (NULL == sqlite_file)
        return;
            
    rc = sqlite3_initialize();
    if (SQLITE_OK != rc)
        gmx_fatal(FARGS,"Initializing sqlite. Sqlite3 code %d.",rc);
        
    db = NULL;
    rc = sqlite3_open_v2(sqlite_file,&db,SQLITE_OPEN_READONLY,NULL);
    if (SQLITE_OK != rc)
        gmx_fatal(FARGS,"Opening sqlite database %s in read-only mode. Sqlite3 code %d.",
                  sqlite_file,rc);
    /* Now database is open and everything is Hunky Dory */
    
    /* Seems like we're done, close down and say goodbye */
    rc = sqlite3_close(db);
    if (SQLITE_OK != rc)
        gmx_fatal(FARGS,"Closing sqlite database %s. Sqlite3 code %d.",
                  sqlite_file,rc);
    
    rc = sqlite3_shutdown();
    if (SQLITE_OK != rc)
        gmx_fatal(FARGS,"Shutting down sqlite. Sqlite3 code %d.",rc);
#else
    fprintf(stderr,"No support for sqlite3 database in this executable.\n");
    fprintf(stderr,"Please rebuild gromacs with cmake flag -DGMX_SQLITE3=ON set.\n");
#endif
}

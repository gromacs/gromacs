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
#ifdef HAVE_LIBSQLITE3
#include <sqlite3.h>
#endif
#include "gmx_fatal.h"
#include "futil.h"
#include "macros.h"
#include "grompp.h"
#include "smalloc.h"
#include "string2.h"
#include "molprop_sqlite3.h"

static void my_memcallback(void *ptr)
{
    ;
}

#ifdef HAVE_LIBSQLITE3
static void check_sqlite3(sqlite3 *db,char *extra,int rc)
{
    const char *msg;
    
    if (NULL != db) 
    {
        if (SQLITE_OK != sqlite3_errcode(db)) 
        {
            msg = sqlite3_errmsg(db);
            sqlite3_close(db);
            sqlite3_shutdown();
            gmx_fatal(FARGS,"%s: %s",extra,msg);  
        }
    }
    else if (SQLITE_OK != rc)
    {
        gmx_fatal(FARGS,"%s",extra);  
    }
}
#endif

void gmx_molprop_read_sqlite3(int np,gmx_molprop_t mp[],const char *sqlite_file)
{
#ifdef HAVE_LIBSQLITE3
    sqlite3 *db = NULL;
    sqlite3_stmt *stmt=NULL;
    char sql_str[1024];
    char *iupac,*iupac2,*prop,*ref;
    double value,error;
    int i,cidx,expref,rc;
    
    if (NULL == sqlite_file)
        return;
            
    check_sqlite3(NULL,"Initializing sqlite",
                  sqlite3_initialize());
        
    check_sqlite3(NULL,"Opening sqlite database n read-only mode",
                  sqlite3_open_v2(sqlite_file,&db,SQLITE_OPEN_READONLY,NULL));
    
    /* Now database is open and everything is Hunky Dory */
    fprintf(stderr,"Opened SQLite3 database %s\n",sqlite_file);
    
    /* Now present a query statement */
    sprintf(sql_str,"SELECT mol.iupac,pt.prop,gp.value,gp.error,ref.ref FROM molecules as mol,gasproperty as gp,proptypes as pt,reference as ref WHERE (mol.molid = gp.molid) AND (gp.propid = pt.propid) AND (gp.refid = ref.refid) AND (mol.iupac = \"?\")");
    check_sqlite3(db,"Preparing statement",
                  sqlite3_prepare_v2(db,sql_str,strlen(sql_str),&stmt,NULL));
    for(i=0; (i<np); i++) 
    {
        iupac = gmx_molprop_get_iupac(mp[i]);
        if (NULL != iupac)
        {
            do 
            {
                check_sqlite3(db,"Binding text",
                              sqlite3_bind_text(stmt,1,iupac,-1,NULL));
                rc = sqlite3_step(stmt);
                check_sqlite3(db,"Stepping",rc);
                cidx   = 0;
                iupac2 = sqlite3_column_text(stmt,cidx++);
                prop   = sqlite3_column_text(stmt,cidx++);
                value  = sqlite3_column_double(stmt,cidx++);
                error  = sqlite3_column_double(stmt,cidx++);
                ref    = sqlite3_column_text(stmt,cidx++);
                if (strcasecmp(iupac,iupac2) != 0)
                    gmx_fatal(FARGS,"Selected '%s' from database but got '%s'. WTF?!",
                              iupac,iupac2);
                gmx_molprop_add_experiment(mp[i],ref,"minimum",&expref);
                if (strcasecmp(prop,"Polarizability") == 0)
                    gmx_molprop_add_polar(mp[i],expref,"Experiment","Unit",
                                          0,0,0,value,error);
                
            } while (SQLITE_ROW == rc);
            check_sqlite3(db,"Resetting sqlite3 statement",
                          sqlite3_reset(stmt));
        }
    }
    check_sqlite3(db,"Finalizing sqlite3 statement",
                  sqlite3_finalize(stmt));
    
    /* Seems like we're done, close down and say goodbye */
    check_sqlite3(NULL,"Closing sqlite database",
                  sqlite3_close(db));
    
    check_sqlite3(NULL,"Shutting down sqlite. Sqlite3 code %d.",
                  sqlite3_shutdown());
                  
#else
    fprintf(stderr,"No support for sqlite3 database in this executable.\n");
    fprintf(stderr,"Please rebuild gromacs with cmake flag -DGMX_SQLITE3=ON set.\n");
#endif
}

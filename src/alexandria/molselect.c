/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: molselect.c,v 1.51 2009/06/01 06:13:18 spoel Exp $
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
#include <strdb.h>
#include <string2.h>
#include <smalloc.h>
#include <strdb.h>
#include "molselect.h"

const char *ims_names[imsNR] = { "Train", "Test", "Ignore", "Unknown" };

typedef struct {
    const char *iupac;
    int  status,index;
} t_ims;

typedef struct gmx_molselect {
    int nmol;
    t_ims *ims;
} gmx_molselect;

static int ims_comp(const void *a,const void *b)
{
    t_ims *ia = (t_ims *)a;
    t_ims *ib = (t_ims *)b;
  
    return strcasecmp(ia->iupac,ib->iupac);
}

gmx_molselect_t gmx_molselect_init(const char *fn)
{
    gmx_molselect *gms;
    char **strings;
    char **tmp;
    int i,j;
  
    snew(gms,1);
    gms->nmol = get_file(fn,&strings);
    snew(gms->ims,gms->nmol);
    for(i=0; (i<gms->nmol); i++) {
        tmp = split('|',strings[i]);
        if (NULL != tmp[0])
        {
            gms->ims[i].iupac = strdup(tmp[0]);
            gms->ims[i].index   = i+1;
            if (NULL != tmp[1]) {
                for(j=0; (j<imsNR); j++) 
                    if (strcasecmp(ims_names[j],tmp[1]) == 0)
                        break;
                if (j < imsNR)
                    gms->ims[i].status = j;
                else
                {
                    gms->ims[i].status = imsUnknown;
                    fprintf(stderr,"Unknown status '%s' for molecule %s on line %d in file %s\n",
                            tmp[1],tmp[0],i,fn);
                }
            }
            else {
                gms->ims[i].status = imsUnknown;
                fprintf(stderr,"No status field for molecule %s on line %d in file %s\n",
                        tmp[0],i,fn);
            }
        }
        sfree(strings[i]);
        j = 0;
        while (NULL != tmp[j])
        {
            sfree(tmp[j]);
            j++;
        }
    }
    sfree(strings);
  
    /* Sort the molecules for faster searching down below */
    qsort(gms->ims,gms->nmol,sizeof(gms->ims[0]),ims_comp);
  
    return (gmx_molselect_t) gms;
}

void gmx_molselect_done(gmx_molselect_t gms)
{
    gmx_molselect *g = (gmx_molselect *)gms;
    int i;
  
    for(i=0; (i<g->nmol); i++)
        sfree(g->ims[i].iupac);
    sfree(g->ims);
    g->ims = NULL;
    sfree(g);
}

int gmx_molselect_status(gmx_molselect_t gms,const char *iupac)
{
    gmx_molselect *g = (gmx_molselect *)gms;
    t_ims key;
    t_ims *ims;

    if (NULL == iupac)
        return imsUnknown;
    key.iupac = iupac;
    ims = (t_ims *)bsearch(&key,g->ims,g->nmol,sizeof(g->ims[0]),ims_comp);
    if (NULL != ims)
        return ims->status;
    else
        return imsUnknown;
}

int gmx_molselect_index(gmx_molselect_t gms,const char *iupac)
{
    gmx_molselect *g = (gmx_molselect *)gms;
    t_ims key;
    t_ims *ims;

    if (NULL == iupac)
        return g->nmol;
    key.iupac = iupac;
    key.status = 1;
    key.index = 0;
    ims = (t_ims *)bsearch(&key,g->ims,g->nmol,sizeof(g->ims[0]),ims_comp);
    if (NULL != ims)
        return ims->index;
    else
        return g->nmol;
}


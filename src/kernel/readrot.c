/*
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vec.h"
#include "smalloc.h"
#include "readir.h"
#include "names.h"

static char s_vec[STRLEN];


static void string2dvec(char buf[], dvec nums)
{
    if (sscanf(buf,"%lf%lf%lf",&nums[0],&nums[1],&nums[2]) != 3)
        gmx_fatal(FARGS,"Expected three numbers at input line %s",buf);
}


extern char **read_rotparams(int *ninp_p,t_inpfile **inp_p,t_rot *rot) 
{
    int  ninp,g,m;
    t_inpfile *inp;
    const char *tmp;
    char **grpbuf;
    char buf[STRLEN];
    dvec vec;
    t_rotgrp *rotg;

    ninp   = *ninp_p;
    inp    = *inp_p;
    
    /* read rotation parameters */
    ITYPE("rot_nstrout",     rot->nstrout, 100);
    ITYPE("rot_nsttout",     rot->nsttout, 1000);
    CTYPE("Number of rotation groups");
    ITYPE("rot_ngroups",     rot->ngrp,1);
    
    if (rot->ngrp < 1) {
        gmx_fatal(FARGS,"rot_ngroups should be >= 1");
    }
    
    snew(rot->grp,rot->ngrp);
    
    /* Read the rotation groups */
    snew(grpbuf,rot->ngrp);
    for(g=0; g<rot->ngrp; g++) {
        rotg = &rot->grp[g];
        snew(grpbuf[g],STRLEN);
        CTYPE("Group name");
        sprintf(buf,"rot_group%d",g);
        STYPE(buf,              grpbuf[g], "");
        
        CTYPE("Rotation type can be Fixed, FixedPlane, Flexible1 or Flexible2");
        sprintf(buf,"rot_type%d",g);
        ETYPE(buf,              rotg->eType, erotg_names);
        
        CTYPE("Rotation vector");
        sprintf(buf,"rot_vec%d",g);
        STYPE(buf,              s_vec, "1.0 0.0 0.0");
        string2dvec(s_vec,vec);
        /* Normalize the rotation vector */
        if (dnorm(vec) != 0)
            dsvmul(1.0/dnorm(vec),vec,vec);
        else
            gmx_fatal(FARGS,"rot_vec%d must not be 0!", g);
        fprintf(stderr, "Enforced rotation: Group %d (%s) normalized rot. vector: %f %f %f\n", 
                g, erotg_names[rotg->eType], vec[0], vec[1], vec[2]);
        for(m=0; m<DIM; m++)
            rotg->vec[m] = vec[m];
        
        CTYPE("Emission point for the fixed axis");
        sprintf(buf,"rot_offset%d",g);
        STYPE(buf,              s_vec, "0.0 0.0 0.0");
        string2dvec(s_vec,vec);
        for(m=0; m<DIM; m++)
            rotg->offset[m] = vec[m];

        CTYPE("Rotation rate (nm/ps) and force constant [kJ/(mol*nm^2)]");
        sprintf(buf,"rot_rate%d",g);
        RTYPE(buf,              rotg->rate, 0.0);
        sprintf(buf,"rot_k%d",g);
        RTYPE(buf,              rotg->k, 0.0);
        CTYPE("Slab distance for flexible rotation (nm)");
        sprintf(buf,"rot_slab_distance%d",g);
        RTYPE(buf,              rotg->slab_dist, 1.5);
        CTYPE("Minimum value of Gaussian for the force to be evaluated");
        sprintf(buf,"rot_min_gaussian%d",g);
        RTYPE(buf,              rotg->min_gaussian, 1e-3);
    }
    
    *ninp_p   = ninp;
    *inp_p    = inp;
    
    return grpbuf;
}


extern void make_rotation_groups(t_rot *rot,char **rotgnames,t_blocka *grps,char **gnames)
{
    int      g,ig=-1,i;
    t_rotgrp *rotgrp;
    
    
    for (g=0; g<rot->ngrp; g++)
    {
        rotgrp = &rot->grp[g];
        if (g == 0 && strcmp(rotgnames[g],"") == 0)
        {
            rotgrp->nat = 0;
        }
        else
        {
            ig = search_string(rotgnames[g],grps->nr,gnames);
            rotgrp->nat = grps->index[ig+1] - grps->index[ig];
        }
        
        if (rotgrp->nat > 0)
        {
            fprintf(stderr,"Rotation group %d '%s' has %d atoms\n",g,rotgnames[g],rotgrp->nat);
            snew(rotgrp->ind,rotgrp->nat);
            for(i=0; i<rotgrp->nat; i++)
                rotgrp->ind[i] = grps->a[grps->index[ig]+i];            
        }
        else
            gmx_fatal(FARGS,"Rotation group %d '%s' is empty",g,rotgnames[g]);
    }
}

/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * $Id: nonbonded.c,v 1.36 2008/12/03 16:07:05 hess Exp $
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
 
#include "adress.h"
#include <math.h>
#include "types/simple.h"
#include "typedefs.h"
#include "vec.h"

real 
adress_weight(real            x,
              real            y,
              real            z,
              int             adresstype,
              real            adressr,
              real            adressw,
              real            refx,
              real            refy,
              real            refz)
{
    real l2 = adressr+adressw;
    real sqr_dl,dl;
    real tmp;
    
    switch(adresstype)
    {
    case 1:              
        /* constant value for weighting function = adressw */
        return adressw;
    case 2:              
        /* plane through center of box, varies in x direction */
        sqr_dl = (x-refx)*(x-refx);
        break;
    case 3:
        /* point at center of box, assuming cubic geometry */
        sqr_dl = (x-refx)*(x-refx)+(y-refy)*(y-refy)+(z-refz)*(z-refz);
        break;
    case 4:
        /* get reference from reference solute, still need to figure out how to read from index */
        sqr_dl = (x-refx)*(x-refx)+(y-refy)*(y-refy)+(z-refz)*(z-refz);
        break;
    default:
        /* default to explicit simulation */
        return 1;
    }
    
    dl=sqrt(sqr_dl);
    
    /* molecule is coarse grained */
    if (dl > l2)
    {
        return 0;
    }
    /* molecule is explicit */
    else if (dl < adressr)
    {
        return 1;
    }
    /* hybrid region */
    else
    {
        tmp=cos((dl-adressr)*M_PI/2/l2);
        return tmp*tmp;
    }
}

void
update_adress_weights(t_forcerec *         fr,
                      t_mdatoms *          mdatoms,
                      rvec                 x[],
                      matrix               box)
{
    int            i,nr;
    int            adresstype;
    real           adressr;
    real           adressw;
    real           ix,iy,iz;
    real           refx,refy,refz;
    real *         wf;
    unsigned short * ptype;
    nr                 = mdatoms->homenr;
    adresstype         = fr->userint1;
    adressr            = fr->userreal1;
    adressw            = fr->userreal2;
    wf                 = mdatoms->wf;
    ptype              = mdatoms->ptype;

    if(adresstype == 4)
    {
        /* get refx,refy,refz from reference solute */
    }
    else
    {
        /* reference is the center of the box */
        refx           = box[XX][XX]/2.0;
        refy           = box[YY][YY]/2.0;
        refz           = box[ZZ][ZZ]/2.0;
    }

    for(i=0;i<nr;i++)
    {
        if(ptype[i] == eptVSite)
        {
            ix             = x[i][0];
            iy             = x[i][1];
            iz             = x[i][2];
            wf[i]          = adress_weight(ix,iy,iz,adresstype,adressr,adressw,refx,refy,refz);
        }
    }
}

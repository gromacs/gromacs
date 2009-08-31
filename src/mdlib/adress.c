/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.5
 * Written by Christoph Junghans, Brad Lambeth, and possibly others.
 * Copyright (c) 2009 Christoph Junghans, Brad Lambeth.
 * All rights reserved.

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
adress_weight(rvec            x,
              int             adresstype,
              real            adressr,
              real            adressw,
              rvec            ref,
              rvec            box2,
              matrix          box)
{
    int  i;
    real l2 = adressr+adressw;
    real dx,dx2;
    real sqr_dl,dl;
    real tmp;

    sqr_dl = 0.0;

    switch(adresstype)
    {
    case 1:              
        /* constant value for weighting function = adressw */
        return adressw;
    case 2:              
        /* plane through center of box, varies in x direction */
        dx             = x[0]-ref[0];
        sqr_dl         = dx*dx;
        break;
    case 3:
        /* point at center of box, assuming cubic geometry */
        for(i=0;i<3;i++){
            dx         = x[i]-ref[i];
            sqr_dl    += dx*dx;
        }
        break;
    case 4:
        /* get reference from shortest distance to reference solute */
        for(i=0;i<3;i++){
            dx         = x[i]-ref[i];
            dx2        = dx*dx;
            if(dx2 > box2[i]){
                while(dx2 > box2[i]){
                    if(dx<0){
                        dx += box[i][i];
                    }
                    else{
                        dx -= box[i][i];
                    }
                    dx2    = dx*dx;
                }
            }
            sqr_dl    += dx2;
        }
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
#ifndef ADRESS_SWITCHFCT_NEW
        tmp=cos((dl-adressr)*M_PI/2/l2);
#else	
        /* shift the weight past the long flat part.  this makes 
         * the correction look like the old cos^2 function, and 
         * is approximately correct (f(0,1)=1E-8,f'(0,1)=1E-10).
         * To use unstretched correction, change prefactor below and 
         * adjust prefactor in src/gmxlib/nonbonded/nb_generic_cg.c */
        tmp=1.0-(dl-adressr)/l2;
#endif
        return tmp;
    }
}

void
update_adress_weights(t_forcerec *         fr,
                      t_mdatoms *          mdatoms,
                      rvec                 x[],
                      matrix               box)
{
    int            i,j,k,nr;
    int            adresstype;
    real           adressr;
    real           adressw;
    rvec           ix;
    rvec           ref;
    rvec           box2;
    real *         wf;
    nr                 = mdatoms->homenr;
    adresstype         = fr->userint1;
    adressr            = fr->userreal1;
    adressw            = fr->userreal2;
    wf                 = mdatoms->wf;

    if(adresstype == 4)
    {
        /* get refx,refy,refz from reference solute 
         * for now, assume its the last molecule */  
        j              = nr-1;
        for(i=0;i<3;i++) {
            /* need square of half the box length for shortest distance to solute */
            box2[i]    = box[i][i]/2.0;
            box2[i]   *= box2[i];
            ref[i]     = x[j][i];
        }
    }
    else
    {
        /* reference is the center of the box */
        for(i=0;i<3;i++){
            ref[i]     = box[i][i]/2.0;
        }
        /* avoid compiler warning */
        for(i=0;i<3;i++){
            box2[i]    = 0.0;
        }
    }

    k=0;
    for(i=0;i<nr;i++)
    {
        /* only calculate wf for virtual particles */
        if(mdatoms->ptype[i] == eptVSite) 
        {
            for(j=0;j<3;j++){
                ix[j]      = x[i][j];
            }
            wf[i]          = adress_weight(ix,adresstype,adressr,adressw,ref,box2,box);
            /* Assign wf value to explicit atoms of the molecule 
             * this requires that every molecule end in a virtual 
             * site which determines the weight of the molecule. 
             * The molecule may contain other virtual sites, in 
             * that case we'll need to be smarter about assigning
             * the weights. This will work for now. */
            for(j=k;j<i;j++){
                wf[j]      = wf[i];
            }
            k=i+1;
        }
    }
}

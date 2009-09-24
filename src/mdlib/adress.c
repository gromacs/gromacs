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
    case eAdressOff:
        /* default to explicit simulation */
        return 1;
    case eAdressConst:              
        /* constant value for weighting function = adressw */
        return adressw;
    case eAdressXSplit:              
        /* plane through center of box, varies in x direction */
        dx             = x[0]-ref[0];
        sqr_dl         = dx*dx;
        break;
    case eAdressSphere:
        /* point at center of box, assuming cubic geometry */
        for(i=0;i<3;i++){
            dx         = x[i]-ref[i];
            sqr_dl    += dx*dx;
        }
        break;
    case eAdressRefMol:
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
        tmp=cos((dl-adressr)*M_PI/2/adressw);
        return tmp*tmp;
#else	
        /* shift the weight past the long flat part.  this makes 
         * the correction look like the old cos^2 function, and 
         * is approximately correct (f(0,1)=1E-8,f'(0,1)=1E-10).
         * To use unstretched correction, change prefactor below and 
         * adjust prefactor in src/gmxlib/nonbonded/nb_generic_cg.c */
        tmp=1.0-(dl-adressr)/adressw;
        return tmp;
#endif
    }
}

void
get_adress_ref(int             adresstype,
               matrix          box,
               rvec            box2,
               rvec            ref)
{
    int i;

    if(adresstype == eAdressRefMol)
    {
        /* get refx,refy,refz from reference solute 
         * for now, assume its the last molecule */  
//        j              = nr-1;
        for(i=0;i<3;i++) {
            /* need square of half the box length for shortest distance to solute */
            box2[i]    = box[i][i]/2.0;
            box2[i]   *= box2[i];
//            ref[i]     = x[j][i];
            ref[i]     = box[i][i]/2.0;
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
}

void
update_adress_weights_com(FILE *               fplog,
                          int                  cg0,
                          int                  cg1,
                          t_block *            cgs,
                          rvec                 x[],
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          matrix               box)
{
    int            icg,k,k0,k1,d;
    real           nrcg,inv_ncg,mtot,inv_mtot;
    atom_id *      cgindex;
    int            adresstype;
    real           adressr,adressw;
    rvec           ix,ref,box2;
    real *         massT;
    real *         wf;

    adresstype         = fr->adress_type;
    adressr            = fr->adress_ex_width;
    adressw            = fr->adress_hy_width;
    massT              = mdatoms->massT;
    wf                 = mdatoms->wf;

    get_adress_ref(adresstype,box,box2,ref);

    /* Since this is center of mass AdResS, the vsite is not guaranteed
     * to be on the same node as the constructing atoms.  Therefore we 
     * loop over the charge groups, calculate their center of mass,
     * then use this to calculate wf for each atom.  This wastes vsite
     * construction, but it's the only way to assure that the explicit
     * atoms have the same wf as their vsite. */

#ifdef DEBUG
    fprintf(fplog,"Calculating center of mass for charge groups %d to %d\n",
            cg0,cg1);
#endif
    cgindex = cgs->index;
    
    /* Compute the center of mass for all charge groups */
    for(icg=cg0; (icg<cg1); icg++) 
    {
        k0      = cgindex[icg];
        k1      = cgindex[icg+1];
        nrcg    = k1-k0;
        if (nrcg == 1)
        {
            wf[k0] = adress_weight(x[k0],adresstype,adressr,adressw,ref,box2,box);
        }
        else
        {
            mtot = 0.0;
            for(k=k0; (k<k1); k++)
            {
                mtot += massT[k];
            }
            if (mtot > 0.0)
            {
                inv_mtot = 1.0/mtot;
                
                clear_rvec(ix);
                for(k=k0; (k<k1); k++)
                {
                    for(d=0; (d<DIM); d++)
                    {
                        ix[d] += x[k][d]*massT[k];
                    }
                }
                for(d=0; (d<DIM); d++)
                {
                    ix[d] *= inv_mtot;
                }
            }
            /* Calculate the center of gravity if the charge group mtot=0 (only vsites) */
            else
            {
                inv_ncg = 1.0/nrcg;

                clear_rvec(ix);
                for(k=k0; (k<k1); k++)
                {
                    for(d=0; (d<DIM); d++)
                    {
                        ix[d] += x[k][d];
                    }
                }
                for(d=0; (d<DIM); d++)
                {
                    ix[d] *= inv_ncg;
                }
            }

            /* Set wf of all atoms in charge group equal to wf of com */
            wf[k0] = adress_weight(ix,adresstype,adressr,adressw,ref,box2,box);
            for(k=(k0+1); (k<k1); k++)
            {
                wf[k] = wf[k0];
            }
        }
    }
}
        
void
update_adress_weights_cog(t_iparams            ip[],
                          t_ilist              ilist[],
                          rvec                 x[],
                          t_forcerec *         fr,
                          t_mdatoms *          mdatoms,
                          matrix               box)
{
    int            i,j,k,nr,nra,inc;
    int            ftype,adresstype;
    t_iatom        avsite,ai,aj,ak,al;
    t_iatom *      ia;
    real           adressr,adressw;
    rvec           ref,box2;
    real *         wf;

    adresstype         = fr->adress_type;
    adressr            = fr->adress_ex_width;
    adressw            = fr->adress_hy_width;
    wf                 = mdatoms->wf;

    get_adress_ref(adresstype,box,box2,ref);

    /* Since this is center of gravity AdResS, we know the vsite
     * is in the same charge group as the constructing atoms.
     * Loop over vsite types, calculate the weight of the vsite,
     * then assign that weight to the constructing atoms.  We
     * shouldn't need to worry about pbc since this was taken
     * care of during vsite construction, which necessarily comes
     * before this. */

    for(ftype=0; (ftype<F_NRE); ftype++) 
    {
        if (interaction_function[ftype].flags & IF_VSITE) 
        {
            nra    = interaction_function[ftype].nratoms;
            nr     = ilist[ftype].nr;
            ia     = ilist[ftype].iatoms;
            
            for(i=0; (i<nr); ) 
            {
                /* The vsite and first constructing atom */
                avsite     = ia[1];
                ai         = ia[2];
                wf[avsite] = adress_weight(x[avsite],adresstype,adressr,adressw,ref,box2,box);
                wf[ai]     = wf[avsite];

                /* Assign the vsite wf to rest of constructing atoms depending on type */
                inc = nra+1;
                switch (ftype) {
                case F_VSITE2:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    break;
                case F_VSITE3:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    ak     = ia[4];
                    wf[ak] = wf[avsite];
                    break;
                case F_VSITE3FD:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    ak     = ia[4];
                    wf[ak] = wf[avsite];
                    break;
                case F_VSITE3FAD:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    ak     = ia[4];
                    wf[ak] = wf[avsite];
                    break;
                case F_VSITE3OUT:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    ak     = ia[4];
                    wf[ak] = wf[avsite];
                    break;
                case F_VSITE4FD:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    ak     = ia[4];
                    wf[ak] = wf[avsite];
                    al     = ia[5];
                    wf[al] = wf[avsite];
                    break;
                case F_VSITE4FDN:
                    aj     = ia[3];
                    wf[aj] = wf[avsite];
                    ak     = ia[4];
                    wf[ak] = wf[avsite];
                    al     = ia[5];
                    wf[al] = wf[avsite];
                    break;
                case F_VSITEN:
                    inc    = 3*ip[ia[0]].vsiten.n;
                    for(j=3; j<inc; j+=3) 
                    {
                        ai = ia[j+2];
                        wf[ai] = wf[avsite];
                    }
                    break;
                default:
                    gmx_fatal(FARGS,"No such vsite type %d in %s, line %d",
                              ftype,__FILE__,__LINE__);
                }

                /* Increment loop variables */
                i  += inc;
                ia += inc;
            }
        }
    }
}


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
#include "maths.h"
#include "pbc.h"
#include "types/simple.h"
#include "typedefs.h"
#include "vec.h"

real 
adress_weight(rvec            x,
              int             adresstype,
              real            adressr,
              real            adressw,
              bool            bnew_wf,
              rvec *          ref,
              t_pbc *         pbc)
{
    int  i;
    real l2 = adressr+adressw;
    real sqr_dl,dl;
    real tmp;
    rvec dx;

    sqr_dl = 0.0;

    if (pbc) 
    {
        pbc_dx(pbc,(*ref),x,dx);
    } 
    else 
    {
        rvec_sub((*ref),x,dx);
    }

    switch(adresstype)
    {
    case eAdressOff:
        /* default to explicit simulation */
        return 1;
    case eAdressConst:              
        /* constant value for weighting function = adressw */
        return adressw;
    case eAdressXSplit:              
        /* plane through center of ref, varies in x direction */
        sqr_dl         = dx[0]*dx[0];
        break;
    case eAdressSphere:
        /* point at center of ref, assuming cubic geometry */
        for(i=0;i<3;i++){
            sqr_dl    += dx[i]*dx[i];
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
    else if (bnew_wf)
    {
        tmp=1.0-(dl-adressr)/adressw;
        return tmp;
    }
    else
    {
        tmp=cos((dl-adressr)*M_PI/2/adressw);
        return tmp*tmp;
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
                          t_pbc *              pbc)
{
    int            icg,k,k0,k1,d;
    real           nrcg,inv_ncg,mtot,inv_mtot;
    atom_id *      cgindex;
    rvec           ix;
    int            adresstype;
    real           adressr,adressw;
    bool           bnew_wf;
    rvec *         ref;
    real *         massT;
    real *         wf;

    adresstype         = fr->adress_type;
    adressr            = fr->adress_ex_width;
    adressw            = fr->adress_hy_width;
    bnew_wf            = fr->badress_new_wf;
    massT              = mdatoms->massT;
    wf                 = mdatoms->wf;
    ref                = &(fr->adress_refmol);

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
            wf[k0] = adress_weight(x[k0],adresstype,adressr,adressw,bnew_wf,ref,pbc);
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
            wf[k0] = adress_weight(ix,adresstype,adressr,adressw,bnew_wf,ref,pbc);
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
                          t_pbc *              pbc)
{
    int            i,j,k,nr,nra,inc;
    int            ftype,adresstype;
    t_iatom        avsite,ai,aj,ak,al;
    t_iatom *      ia;
    real           adressr,adressw;
    bool           bnew_wf;
    rvec *         ref;
    real *         wf;

    adresstype         = fr->adress_type;
    adressr            = fr->adress_ex_width;
    adressw            = fr->adress_hy_width;
    bnew_wf            = fr->badress_new_wf;
    wf                 = mdatoms->wf;
    ref                = &(fr->adress_refmol);

    /* Since this is center of geometry AdResS, we know the vsite
     * is in the same charge group node as the constructing atoms.
     * Loop over vsite types, calculate the weight of the vsite,
     * then assign that weight to the constructing atoms. */

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
                wf[avsite] = adress_weight(x[avsite],adresstype,adressr,adressw,bnew_wf,ref,pbc);
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

void
adress_thermo_force(int                  cg0,
                    int                  cg1,
                    t_block *            cgs,
                    rvec                 x[],
                    rvec                 f[],
                    t_forcerec *         fr,
                    t_mdatoms *          mdatoms,
                    t_pbc *              pbc)
{
    int              icg,k0,k1,n0,nnn,nrcg;
    int              adresstype;
    real             adressw;
    bool             bnew_wf;
    atom_id *        cgindex;
    unsigned short * ptype;
    rvec *           ref;
    real *           wf;
    real             tabscale;
    real *           ATFtab;
    rvec             dr;
    real             w,wsq,wmin1,wmin1sq,wp,wt,rinv;
    real             eps,eps2,F,Geps,Heps2,Fp,dmu_dwp,dwp_dr,fscal;

    adresstype       = fr->adress_type;
    adressw          = fr->adress_hy_width;
    bnew_wf          = fr->badress_new_wf;
    cgindex          = cgs->index;
    ptype            = mdatoms->ptype;
    ref              = &(fr->adress_refmol);
    wf               = mdatoms->wf;
    tabscale         = fr->atf_tab.scale;
    ATFtab           = fr->atf_tab.tab;

    for(icg=cg0; (icg<cg1); icg++)
    {
        k0           = cgindex[icg];
        k1           = cgindex[icg+1];
        nrcg         = k1-k0;
        /* avoid confusing TIP4P vsite with CG vsite */
        if (nrcg == 1)
        {
            if (ptype[k0] == eptVSite)
            {
                w    = wf[k0];
                /* is it hybrid? */
                if (w > 0 && w < 1)
                {
                    fscal            = 0;
                    if (pbc)
                    {
                        pbc_dx(pbc,(*ref),x[k0],dr);
                    }
                    else
                    {
                        rvec_sub((*ref),x[k0],dr);
                    }

                    rinv             = gmx_invsqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
                    wsq              = w*w;
                    wmin1            = w - 1.0;
                    wmin1sq          = wmin1*wmin1;
                    dwp_dr           = 1.0/adressw;

                    /* get the dwp_dr part of the force */
                    if(bnew_wf)
                    {
                        wp           = 1.0-exp(-M_LN2*wsq/wmin1sq);
                        dwp_dr      *= M_LN2*w*exp(M_LN2*(1.0-2.0*w)/wmin1sq)/(wmin1*wmin1sq);
                    }
                    else
                    {
                        wp           = wsq;
                        dwp_dr      *= -2.0*M_PI*sqrt(wsq*w*(1.0-w));
                    }

                    /* now get the dmu/dwp part from the table */
                    wt               = wp*tabscale;
                    n0               = wt;
                    eps              = wt-n0;
                    eps2             = eps*eps;
                    nnn              = 4*n0;
                    F                = ATFtab[nnn+1];
                    Geps             = eps*ATFtab[nnn+2];
                    Heps2            = eps2*ATFtab[nnn+3];
                    Fp               = F+Geps+Heps2;
                    dmu_dwp          = -(Fp+Geps+2.0*Heps2)*tabscale;

                    fscal            = dmu_dwp*dwp_dr*rinv;

                    /* now add thermo force to f_novirsum */
                    f[k0][0]        += fscal*dr[0];
                    if (adresstype != eAdressXSplit)
                    {
                        f[k0][1]    += fscal*dr[1];
                        f[k0][2]    += fscal*dr[2];
                    }
                }
            }
        }
    }
}

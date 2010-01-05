/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
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

#include <math.h>
#include "typedefs.h"
#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "physics.h"
#include "vec.h"
#include "futil.h"
#include "xvgr.h"
#include "gmx_fatal.h"
#include "bondf.h"
#include "copyrite.h"
#include "disre.h"
#include "main.h"
#include "mtop_util.h"
#include "dihre.h"

void init_dihres(FILE *fplog,gmx_mtop_t *mtop,t_inputrec *ir,t_fcdata *fcd)
{
  int count;

  fcd->dihre_fc = ir->dihre_fc;

  count = gmx_mtop_ftype_count(mtop,F_DIHRES);

  if (fplog && count) {
    fprintf(fplog,"There are %d dihedral restraints\n",count);
  }
}

real ta_dihres(int nfa,const t_iatom forceatoms[],const t_iparams ip[],
               const rvec x[],rvec f[],rvec fshift[],
               const t_pbc *pbc,const t_graph *g,
               real lambda,real *dvdl,
               const t_mdatoms *md,t_fcdata *fcd,
               int *ddgatindex)
{
    real vtot = 0;
    int  ai,aj,ak,al,i,k,type,typep,label,power,t1,t2,t3;
    real phi0A,phi0B,dphiA,dphiB,kfacA,kfacB,phi0,dphi,kfac;
    real phi,ddphi,ddp,ddp2,dp,sign,d2r,fc,L1;
    rvec r_ij,r_kj,r_kl,m,n;
    
    L1 = 1.0-lambda;
    
    fc  = fcd->dihre_fc;
    d2r = DEG2RAD;
    k   = 0;
    
    for(i=0; (i<nfa); ) {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        ak   = forceatoms[i++];
        al   = forceatoms[i++];
        
        phi0A  = ip[type].dihres.phiA*d2r;
        dphiA  = ip[type].dihres.dphiA*d2r;
        kfacA  = ip[type].dihres.kfacA*fc; 
        
        phi0B  = ip[type].dihres.phiB*d2r;
        dphiB  = ip[type].dihres.dphiB*d2r;
        kfacB  = ip[type].dihres.kfacB*fc; 
        
        phi0  = L1*phi0A + lambda*phi0B;
        dphi  = L1*dphiA + lambda*dphiB;
        kfac = L1*kfacA + lambda*kfacB;
        
        power = ip[type].dihres.power;
        label = ip[type].dihres.label;
        
        phi = dih_angle(x[ai],x[aj],x[ak],x[al],pbc,r_ij,r_kj,r_kl,m,n,
                        &sign,&t1,&t2,&t3);	  
        /* 84 flops */
        
        if (debug)
            fprintf(debug,"dihres[%d]: %d %d %d %d : phi=%f, dphi=%f, kfac=%f, power=%d, label=%d\n",
                    k++,ai,aj,ak,al,phi0,dphi,kfac,power,label);
        
        /* phi can jump if phi0 is close to Pi/-Pi, which will cause huge
         * force changes if we just apply a normal harmonic.
         * Instead, we first calculate phi-phi0 and take it modulo (-Pi,Pi).
         * This means we will never have the periodicity problem, unless
         * the dihedral is Pi away from phiO, which is very unlikely due to
         * the potential.
         */
        dp = phi-phi0;
        make_dp_periodic(&dp);

        if (dp > dphi) 
        {
            ddp = dp-dphi;
        }
        else if (dp < -dphi)
        { 
            ddp = dp+dphi;
        }
        else 
        {
            ddp = 0;
        }
        
        if (ddp != 0.0) 
        {
            ddp2 = ddp*ddp;
            vtot += 0.5*kfac*ddp2;
            ddphi = kfac*ddp;
            
            *dvdl += 0.5*(kfacB - kfacA)*ddp2; 	
            /* lambda dependence from changing restraint distances */
            if (ddp > 0)  
            {
                *dvdl -= kfac*ddp*((dphiB - dphiA)+(phi0B - phi0A));  
                
            } 
            else if (ddp < 0 )
            {
                *dvdl += kfac*ddp*((dphiB - dphiA)-(phi0B - phi0A)); 
            }
            
            do_dih_fup(ai,aj,ak,al,ddphi,r_ij,r_kj,r_kl,m,n,
                       f,fshift,pbc,g,x,t1,t2,t3);		/* 112		*/
        }
    }
    return vtot;
}

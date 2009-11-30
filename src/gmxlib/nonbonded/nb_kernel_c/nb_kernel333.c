/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.  
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "vec.h"
#ifdef GMX_THREAD_SHM_FDECOMP
#include "thread_mpi.h"
#endif

#include "nb_kernel333.h"

/*
 * Gromacs nonbonded kernel nb_kernel333
 * Coulomb interaction:     Tabulated
 * VdW interaction:         Tabulated
 * water optimization:      TIP4P - other atoms
 * Calculate forces:        yes
 */
void nb_kernel333(
                    int *           p_nri,
                    int *           iinr,
                    int *           jindex,
                    int *           jjnr,
                    int *           shift,
                    real *          shiftvec,
                    real *          fshift,
                    int *           gid,
                    real *          pos,
                    real *          faction,
                    real *          charge,
                    real *          p_facel,
                    real *          p_krf,
                    real *          p_crf,
                    real *          Vc,
                    int *           type,
                    int *           p_ntype,
                    real *          vdwparam,
                    real *          Vvdw,
                    real *          p_tabscale,
                    real *          VFtab,
                    real *          invsqrta,
                    real *          dvda,
                    real *          p_gbtabscale,
                    real *          GBtab,
                    int *           p_nthreads,
                    int *           count,
                    void *          mtx,
                    int *           outeriter,
                    int *           inneriter,
                    real *          work)
{
    int           nri,ntype,nthreads;
    real          facel,krf,crf,tabscale,gbtabscale;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    int           nn0,nn1,nouter,ninner;
    real          shX,shY,shZ;
    real          fscal,tx,ty,tz;
    real          jq;
    real          qq,vcoul,vctot;
    int           nti;
    int           tj;
    real          Vvdw6,Vvdwtot;
    real          Vvdw12;
    real          r,rt,eps,eps2;
    int           n0,nnn;
    real          Y,F,Geps,Heps2,Fp,VV;
    real          FF;
    real          fijC;
    real          fijD,fijR;
    real          ix1,iy1,iz1,fix1,fiy1,fiz1;
    real          ix2,iy2,iz2,fix2,fiy2,fiz2;
    real          ix3,iy3,iz3,fix3,fiy3,fiz3;
    real          ix4,iy4,iz4,fix4,fiy4,fiz4;
    real          jx1,jy1,jz1,fjx1,fjy1,fjz1;
    real          dx11,dy11,dz11,rsq11,rinv11;
    real          dx21,dy21,dz21,rsq21,rinv21;
    real          dx31,dy31,dz31,rsq31,rinv31;
    real          dx41,dy41,dz41,rsq41,rinv41;
    real          qH,qM;
    real          c6,c12;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    

    /* Initialize water data */
    ii               = iinr[0];        
    qH               = facel*charge[ii+1];
    qM               = facel*charge[ii+3];
    nti              = 2*ntype*type[ii];


    /* Reset outer and inner iteration counters */
    nouter           = 0;              
    ninner           = 0;              

    /* Loop over thread workunits */
    
    do
    {
#ifdef GMX_THREAD_SHM_FDECOMP
        tMPI_Thread_mutex_lock((tMPI_Thread_mutex_t *)mtx);
        nn0              = *count;         
		
        /* Take successively smaller chunks (at least 10 lists) */
        nn1              = nn0+(nri-nn0)/(2*nthreads)+10;
        *count           = nn1;            
        tMPI_Thread_mutex_unlock((tMPI_Thread_mutex_t *)mtx);
        if(nn1>nri) nn1=nri;
#else
	    nn0 = 0;
		nn1 = nri;
#endif
        /* Start outer loop over neighborlists */
        
        for(n=nn0; (n<nn1); n++)
        {

            /* Load shift vector for this list */
            is3              = 3*shift[n];     
            shX              = shiftvec[is3];  
            shY              = shiftvec[is3+1];
            shZ              = shiftvec[is3+2];

            /* Load limits for loop over neighbors */
            nj0              = jindex[n];      
            nj1              = jindex[n+1];    

            /* Get outer coordinate index */
            ii               = iinr[n];        
            ii3              = 3*ii;           

            /* Load i atom data, add shift vector */
            ix1              = shX + pos[ii3+0];
            iy1              = shY + pos[ii3+1];
            iz1              = shZ + pos[ii3+2];
            ix2              = shX + pos[ii3+3];
            iy2              = shY + pos[ii3+4];
            iz2              = shZ + pos[ii3+5];
            ix3              = shX + pos[ii3+6];
            iy3              = shY + pos[ii3+7];
            iz3              = shZ + pos[ii3+8];
            ix4              = shX + pos[ii3+9];
            iy4              = shY + pos[ii3+10];
            iz4              = shZ + pos[ii3+11];

            /* Zero the potential energy for this list */
            vctot            = 0;              
            Vvdwtot          = 0;              

            /* Clear i atom forces */
            fix1             = 0;              
            fiy1             = 0;              
            fiz1             = 0;              
            fix2             = 0;              
            fiy2             = 0;              
            fiz2             = 0;              
            fix3             = 0;              
            fiy3             = 0;              
            fiz3             = 0;              
            fix4             = 0;              
            fiy4             = 0;              
            fiz4             = 0;              
            
            for(k=nj0; (k<nj1); k++)
            {

                /* Get j neighbor index, and coordinate index */
                jnr              = jjnr[k];        
                j3               = 3*jnr;          

                /* load j atom coordinates */
                jx1              = pos[j3+0];      
                jy1              = pos[j3+1];      
                jz1              = pos[j3+2];      

                /* Calculate distance */
                dx11             = ix1 - jx1;      
                dy11             = iy1 - jy1;      
                dz11             = iz1 - jz1;      
                rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
                dx21             = ix2 - jx1;      
                dy21             = iy2 - jy1;      
                dz21             = iz2 - jz1;      
                rsq21            = dx21*dx21+dy21*dy21+dz21*dz21;
                dx31             = ix3 - jx1;      
                dy31             = iy3 - jy1;      
                dz31             = iz3 - jz1;      
                rsq31            = dx31*dx31+dy31*dy31+dz31*dz31;
                dx41             = ix4 - jx1;      
                dy41             = iy4 - jy1;      
                dz41             = iz4 - jz1;      
                rsq41            = dx41*dx41+dy41*dy41+dz41*dz41;

                /* Calculate 1/r and 1/r2 */
                rinv11           = gmx_invsqrt(rsq11);
                rinv21           = gmx_invsqrt(rsq21);
                rinv31           = gmx_invsqrt(rsq31);
                rinv41           = gmx_invsqrt(rsq41);

                /* Load parameters for j atom */
                tj               = nti+2*type[jnr];
                c6               = vdwparam[tj];   
                c12              = vdwparam[tj+1]; 

                /* Calculate table index */
                r                = rsq11*rinv11;   

                /* Calculate table index */
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 12*n0;          

                /* Tabulated VdW interaction - dispersion */
                nnn              = nnn+4;          
                Y                = VFtab[nnn];     
                F                = VFtab[nnn+1];   
                Geps             = eps*VFtab[nnn+2];
                Heps2            = eps2*VFtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                FF               = Fp+Geps+2.0*Heps2;
                Vvdw6            = c6*VV;          
                fijD             = c6*FF;          

                /* Tabulated VdW interaction - repulsion */
                nnn              = nnn+4;          
                Y                = VFtab[nnn];     
                F                = VFtab[nnn+1];   
                Geps             = eps*VFtab[nnn+2];
                Heps2            = eps2*VFtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                FF               = Fp+Geps+2.0*Heps2;
                Vvdw12           = c12*VV;         
                fijR             = c12*FF;         
                Vvdwtot          = Vvdwtot+ Vvdw6 + Vvdw12;
                fscal            = -((fijD+fijR)*tabscale)*rinv11;

                /* Calculate temporary vectorial force */
                tx               = fscal*dx11;     
                ty               = fscal*dy11;     
                tz               = fscal*dz11;     

                /* Increment i atom force */
                fix1             = fix1 + tx;      
                fiy1             = fiy1 + ty;      
                fiz1             = fiz1 + tz;      

                /* Decrement j atom force */
                fjx1             = faction[j3+0] - tx;
                fjy1             = faction[j3+1] - ty;
                fjz1             = faction[j3+2] - tz;

                /* Load parameters for j atom */
                jq               = charge[jnr+0];  
                qq               = qH*jq;          

                /* Calculate table index */
                r                = rsq21*rinv21;   

                /* Calculate table index */
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 12*n0;          

                /* Tabulated coulomb interaction */
                Y                = VFtab[nnn];     
                F                = VFtab[nnn+1];   
                Geps             = eps*VFtab[nnn+2];
                Heps2            = eps2*VFtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                FF               = Fp+Geps+2.0*Heps2;
                vcoul            = qq*VV;          
                fijC             = qq*FF;          
                vctot            = vctot + vcoul;  
                fscal            = -((fijC)*tabscale)*rinv21;

                /* Calculate temporary vectorial force */
                tx               = fscal*dx21;     
                ty               = fscal*dy21;     
                tz               = fscal*dz21;     

                /* Increment i atom force */
                fix2             = fix2 + tx;      
                fiy2             = fiy2 + ty;      
                fiz2             = fiz2 + tz;      

                /* Decrement j atom force */
                fjx1             = fjx1 - tx;      
                fjy1             = fjy1 - ty;      
                fjz1             = fjz1 - tz;      

                /* Load parameters for j atom */

                /* Calculate table index */
                r                = rsq31*rinv31;   

                /* Calculate table index */
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 12*n0;          

                /* Tabulated coulomb interaction */
                Y                = VFtab[nnn];     
                F                = VFtab[nnn+1];   
                Geps             = eps*VFtab[nnn+2];
                Heps2            = eps2*VFtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                FF               = Fp+Geps+2.0*Heps2;
                vcoul            = qq*VV;          
                fijC             = qq*FF;          
                vctot            = vctot + vcoul;  
                fscal            = -((fijC)*tabscale)*rinv31;

                /* Calculate temporary vectorial force */
                tx               = fscal*dx31;     
                ty               = fscal*dy31;     
                tz               = fscal*dz31;     

                /* Increment i atom force */
                fix3             = fix3 + tx;      
                fiy3             = fiy3 + ty;      
                fiz3             = fiz3 + tz;      

                /* Decrement j atom force */
                fjx1             = fjx1 - tx;      
                fjy1             = fjy1 - ty;      
                fjz1             = fjz1 - tz;      

                /* Load parameters for j atom */
                qq               = qM*jq;          

                /* Calculate table index */
                r                = rsq41*rinv41;   

                /* Calculate table index */
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 12*n0;          

                /* Tabulated coulomb interaction */
                Y                = VFtab[nnn];     
                F                = VFtab[nnn+1];   
                Geps             = eps*VFtab[nnn+2];
                Heps2            = eps2*VFtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                FF               = Fp+Geps+2.0*Heps2;
                vcoul            = qq*VV;          
                fijC             = qq*FF;          
                vctot            = vctot + vcoul;  
                fscal            = -((fijC)*tabscale)*rinv41;

                /* Calculate temporary vectorial force */
                tx               = fscal*dx41;     
                ty               = fscal*dy41;     
                tz               = fscal*dz41;     

                /* Increment i atom force */
                fix4             = fix4 + tx;      
                fiy4             = fiy4 + ty;      
                fiz4             = fiz4 + tz;      

                /* Decrement j atom force */
                faction[j3+0]    = fjx1 - tx;      
                faction[j3+1]    = fjy1 - ty;      
                faction[j3+2]    = fjz1 - tz;      

                /* Inner loop uses 179 flops/iteration */
            }
            

            /* Add i forces to mem and shifted force list */
            faction[ii3+0]   = faction[ii3+0] + fix1;
            faction[ii3+1]   = faction[ii3+1] + fiy1;
            faction[ii3+2]   = faction[ii3+2] + fiz1;
            faction[ii3+3]   = faction[ii3+3] + fix2;
            faction[ii3+4]   = faction[ii3+4] + fiy2;
            faction[ii3+5]   = faction[ii3+5] + fiz2;
            faction[ii3+6]   = faction[ii3+6] + fix3;
            faction[ii3+7]   = faction[ii3+7] + fiy3;
            faction[ii3+8]   = faction[ii3+8] + fiz3;
            faction[ii3+9]   = faction[ii3+9] + fix4;
            faction[ii3+10]  = faction[ii3+10] + fiy4;
            faction[ii3+11]  = faction[ii3+11] + fiz4;
            fshift[is3]      = fshift[is3]+fix1+fix2+fix3+fix4;
            fshift[is3+1]    = fshift[is3+1]+fiy1+fiy2+fiy3+fiy4;
            fshift[is3+2]    = fshift[is3+2]+fiz1+fiz2+fiz3+fiz4;

            /* Add potential energies to the group for this list */
            ggid             = gid[n];         
            Vc[ggid]         = Vc[ggid] + vctot;
            Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;

            /* Increment number of inner iterations */
            ninner           = ninner + nj1 - nj0;

            /* Outer loop uses 38 flops/iteration */
        }
        

        /* Increment number of outer iterations */
        nouter           = nouter + nn1 - nn0;
    }
    while (nn1<nri);
    

    /* Write outer/inner iteration count to pointers */
    *outeriter       = nouter;         
    *inneriter       = ninner;         
}





/*
 * Gromacs nonbonded kernel nb_kernel333nf
 * Coulomb interaction:     Tabulated
 * VdW interaction:         Tabulated
 * water optimization:      TIP4P - other atoms
 * Calculate forces:        no
 */
void nb_kernel333nf(
                    int *           p_nri,
                    int *           iinr,
                    int *           jindex,
                    int *           jjnr,
                    int *           shift,
                    real *          shiftvec,
                    real *          fshift,
                    int *           gid,
                    real *          pos,
                    real *          faction,
                    real *          charge,
                    real *          p_facel,
                    real *          p_krf,
                    real *          p_crf,
                    real *          Vc,
                    int *           type,
                    int *           p_ntype,
                    real *          vdwparam,
                    real *          Vvdw,
                    real *          p_tabscale,
                    real *          VFtab,
                    real *          invsqrta,
                    real *          dvda,
                    real *          p_gbtabscale,
                    real *          GBtab,
                    int *           p_nthreads,
                    int *           count,
                    void *          mtx,
                    int *           outeriter,
                    int *           inneriter,
                    real *          work)
{
    int           nri,ntype,nthreads;
    real          facel,krf,crf,tabscale,gbtabscale;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    int           nn0,nn1,nouter,ninner;
    real          shX,shY,shZ;
    real          jq;
    real          qq,vcoul,vctot;
    int           nti;
    int           tj;
    real          Vvdw6,Vvdwtot;
    real          Vvdw12;
    real          r,rt,eps,eps2;
    int           n0,nnn;
    real          Y,F,Geps,Heps2,Fp,VV;
    real          ix1,iy1,iz1;
    real          ix2,iy2,iz2;
    real          ix3,iy3,iz3;
    real          ix4,iy4,iz4;
    real          jx1,jy1,jz1;
    real          dx11,dy11,dz11,rsq11,rinv11;
    real          dx21,dy21,dz21,rsq21,rinv21;
    real          dx31,dy31,dz31,rsq31,rinv31;
    real          dx41,dy41,dz41,rsq41,rinv41;
    real          qH,qM;
    real          c6,c12;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    

    /* Initialize water data */
    ii               = iinr[0];        
    qH               = facel*charge[ii+1];
    qM               = facel*charge[ii+3];
    nti              = 2*ntype*type[ii];


    /* Reset outer and inner iteration counters */
    nouter           = 0;              
    ninner           = 0;              

    /* Loop over thread workunits */
    
    do
    {
#ifdef GMX_THREAD_SHM_FDECOMP
        tMPI_Thread_mutex_lock((tMPI_Thread_mutex_t *)mtx);
        nn0              = *count;         
		
        /* Take successively smaller chunks (at least 10 lists) */
        nn1              = nn0+(nri-nn0)/(2*nthreads)+10;
        *count           = nn1;            
        tMPI_Thread_mutex_unlock((tMPI_Thread_mutex_t *)mtx);
        if(nn1>nri) nn1=nri;
#else
	    nn0 = 0;
		nn1 = nri;
#endif
        /* Start outer loop over neighborlists */
        
        for(n=nn0; (n<nn1); n++)
        {

            /* Load shift vector for this list */
            is3              = 3*shift[n];     
            shX              = shiftvec[is3];  
            shY              = shiftvec[is3+1];
            shZ              = shiftvec[is3+2];

            /* Load limits for loop over neighbors */
            nj0              = jindex[n];      
            nj1              = jindex[n+1];    

            /* Get outer coordinate index */
            ii               = iinr[n];        
            ii3              = 3*ii;           

            /* Load i atom data, add shift vector */
            ix1              = shX + pos[ii3+0];
            iy1              = shY + pos[ii3+1];
            iz1              = shZ + pos[ii3+2];
            ix2              = shX + pos[ii3+3];
            iy2              = shY + pos[ii3+4];
            iz2              = shZ + pos[ii3+5];
            ix3              = shX + pos[ii3+6];
            iy3              = shY + pos[ii3+7];
            iz3              = shZ + pos[ii3+8];
            ix4              = shX + pos[ii3+9];
            iy4              = shY + pos[ii3+10];
            iz4              = shZ + pos[ii3+11];

            /* Zero the potential energy for this list */
            vctot            = 0;              
            Vvdwtot          = 0;              

            /* Clear i atom forces */
            
            for(k=nj0; (k<nj1); k++)
            {

                /* Get j neighbor index, and coordinate index */
                jnr              = jjnr[k];        
                j3               = 3*jnr;          

                /* load j atom coordinates */
                jx1              = pos[j3+0];      
                jy1              = pos[j3+1];      
                jz1              = pos[j3+2];      

                /* Calculate distance */
                dx11             = ix1 - jx1;      
                dy11             = iy1 - jy1;      
                dz11             = iz1 - jz1;      
                rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
                dx21             = ix2 - jx1;      
                dy21             = iy2 - jy1;      
                dz21             = iz2 - jz1;      
                rsq21            = dx21*dx21+dy21*dy21+dz21*dz21;
                dx31             = ix3 - jx1;      
                dy31             = iy3 - jy1;      
                dz31             = iz3 - jz1;      
                rsq31            = dx31*dx31+dy31*dy31+dz31*dz31;
                dx41             = ix4 - jx1;      
                dy41             = iy4 - jy1;      
                dz41             = iz4 - jz1;      
                rsq41            = dx41*dx41+dy41*dy41+dz41*dz41;

                /* Calculate 1/r and 1/r2 */
                rinv11           = gmx_invsqrt(rsq11);
                rinv21           = gmx_invsqrt(rsq21);
                rinv31           = gmx_invsqrt(rsq31);
                rinv41           = gmx_invsqrt(rsq41);

                /* Load parameters for j atom */
                tj               = nti+2*type[jnr];
                c6               = vdwparam[tj];   
                c12              = vdwparam[tj+1]; 

                /* Calculate table index */
                r                = rsq11*rinv11;   

                /* Calculate table index */
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 12*n0;          

                /* Tabulated VdW interaction - dispersion */
                nnn              = nnn+4;          
                Y                = VFtab[nnn];     
                F                = VFtab[nnn+1];   
                Geps             = eps*VFtab[nnn+2];
                Heps2            = eps2*VFtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                Vvdw6            = c6*VV;          

                /* Tabulated VdW interaction - repulsion */
                nnn              = nnn+4;          
                Y                = VFtab[nnn];     
                F                = VFtab[nnn+1];   
                Geps             = eps*VFtab[nnn+2];
                Heps2            = eps2*VFtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                Vvdw12           = c12*VV;         
                Vvdwtot          = Vvdwtot+ Vvdw6 + Vvdw12;

                /* Load parameters for j atom */
                jq               = charge[jnr+0];  
                qq               = qH*jq;          

                /* Calculate table index */
                r                = rsq21*rinv21;   

                /* Calculate table index */
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 12*n0;          

                /* Tabulated coulomb interaction */
                Y                = VFtab[nnn];     
                F                = VFtab[nnn+1];   
                Geps             = eps*VFtab[nnn+2];
                Heps2            = eps2*VFtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                vcoul            = qq*VV;          
                vctot            = vctot + vcoul;  

                /* Load parameters for j atom */

                /* Calculate table index */
                r                = rsq31*rinv31;   

                /* Calculate table index */
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 12*n0;          

                /* Tabulated coulomb interaction */
                Y                = VFtab[nnn];     
                F                = VFtab[nnn+1];   
                Geps             = eps*VFtab[nnn+2];
                Heps2            = eps2*VFtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                vcoul            = qq*VV;          
                vctot            = vctot + vcoul;  

                /* Load parameters for j atom */
                qq               = qM*jq;          

                /* Calculate table index */
                r                = rsq41*rinv41;   

                /* Calculate table index */
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 12*n0;          

                /* Tabulated coulomb interaction */
                Y                = VFtab[nnn];     
                F                = VFtab[nnn+1];   
                Geps             = eps*VFtab[nnn+2];
                Heps2            = eps2*VFtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                vcoul            = qq*VV;          
                vctot            = vctot + vcoul;  

                /* Inner loop uses 110 flops/iteration */
            }
            

            /* Add i forces to mem and shifted force list */

            /* Add potential energies to the group for this list */
            ggid             = gid[n];         
            Vc[ggid]         = Vc[ggid] + vctot;
            Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;

            /* Increment number of inner iterations */
            ninner           = ninner + nj1 - nj0;

            /* Outer loop uses 14 flops/iteration */
        }
        

        /* Increment number of outer iterations */
        nouter           = nouter + nn1 - nn0;
    }
    while (nn1<nri);
    

    /* Write outer/inner iteration count to pointers */
    *outeriter       = nouter;         
    *inneriter       = ninner;         
}



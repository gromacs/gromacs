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

#include "nb_kernel102.h"

/*
 * Gromacs nonbonded kernel nb_kernel102
 * Coulomb interaction:     Normal Coulomb
 * VdW interaction:         Not calculated
 * water optimization:      pairs of SPC/TIP3P interactions
 * Calculate forces:        yes
 */
void nb_kernel102(
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
    real          rinvsq;
    real          qq,vcoul,vctot;
    real          ix1,iy1,iz1,fix1,fiy1,fiz1;
    real          ix2,iy2,iz2,fix2,fiy2,fiz2;
    real          ix3,iy3,iz3,fix3,fiy3,fiz3;
    real          jx1,jy1,jz1,fjx1,fjy1,fjz1;
    real          jx2,jy2,jz2,fjx2,fjy2,fjz2;
    real          jx3,jy3,jz3,fjx3,fjy3,fjz3;
    real          dx11,dy11,dz11,rsq11,rinv11;
    real          dx12,dy12,dz12,rsq12,rinv12;
    real          dx13,dy13,dz13,rsq13,rinv13;
    real          dx21,dy21,dz21,rsq21,rinv21;
    real          dx22,dy22,dz22,rsq22,rinv22;
    real          dx23,dy23,dz23,rsq23,rinv23;
    real          dx31,dy31,dz31,rsq31,rinv31;
    real          dx32,dy32,dz32,rsq32,rinv32;
    real          dx33,dy33,dz33,rsq33,rinv33;
    real          qO,qH,qqOO,qqOH,qqHH;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    

    /* Initialize water data */
    ii               = iinr[0];        
    qO               = charge[ii];     
    qH               = charge[ii+1];   
    qqOO             = facel*qO*qO;    
    qqOH             = facel*qO*qH;    
    qqHH             = facel*qH*qH;    


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

            /* Zero the potential energy for this list */
            vctot            = 0;              

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
            
            for(k=nj0; (k<nj1); k++)
            {

                /* Get j neighbor index, and coordinate index */
                jnr              = jjnr[k];        
                j3               = 3*jnr;          

                /* load j atom coordinates */
                jx1              = pos[j3+0];      
                jy1              = pos[j3+1];      
                jz1              = pos[j3+2];      
                jx2              = pos[j3+3];      
                jy2              = pos[j3+4];      
                jz2              = pos[j3+5];      
                jx3              = pos[j3+6];      
                jy3              = pos[j3+7];      
                jz3              = pos[j3+8];      

                /* Calculate distance */
                dx11             = ix1 - jx1;      
                dy11             = iy1 - jy1;      
                dz11             = iz1 - jz1;      
                rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
                dx12             = ix1 - jx2;      
                dy12             = iy1 - jy2;      
                dz12             = iz1 - jz2;      
                rsq12            = dx12*dx12+dy12*dy12+dz12*dz12;
                dx13             = ix1 - jx3;      
                dy13             = iy1 - jy3;      
                dz13             = iz1 - jz3;      
                rsq13            = dx13*dx13+dy13*dy13+dz13*dz13;
                dx21             = ix2 - jx1;      
                dy21             = iy2 - jy1;      
                dz21             = iz2 - jz1;      
                rsq21            = dx21*dx21+dy21*dy21+dz21*dz21;
                dx22             = ix2 - jx2;      
                dy22             = iy2 - jy2;      
                dz22             = iz2 - jz2;      
                rsq22            = dx22*dx22+dy22*dy22+dz22*dz22;
                dx23             = ix2 - jx3;      
                dy23             = iy2 - jy3;      
                dz23             = iz2 - jz3;      
                rsq23            = dx23*dx23+dy23*dy23+dz23*dz23;
                dx31             = ix3 - jx1;      
                dy31             = iy3 - jy1;      
                dz31             = iz3 - jz1;      
                rsq31            = dx31*dx31+dy31*dy31+dz31*dz31;
                dx32             = ix3 - jx2;      
                dy32             = iy3 - jy2;      
                dz32             = iz3 - jz2;      
                rsq32            = dx32*dx32+dy32*dy32+dz32*dz32;
                dx33             = ix3 - jx3;      
                dy33             = iy3 - jy3;      
                dz33             = iz3 - jz3;      
                rsq33            = dx33*dx33+dy33*dy33+dz33*dz33;

                /* Calculate 1/r and 1/r2 */
                rinv11           = gmx_invsqrt(rsq11);
                rinv12           = gmx_invsqrt(rsq12);
                rinv13           = gmx_invsqrt(rsq13);
                rinv21           = gmx_invsqrt(rsq21);
                rinv22           = gmx_invsqrt(rsq22);
                rinv23           = gmx_invsqrt(rsq23);
                rinv31           = gmx_invsqrt(rsq31);
                rinv32           = gmx_invsqrt(rsq32);
                rinv33           = gmx_invsqrt(rsq33);

                /* Load parameters for j atom */
                qq               = qqOO;           
                rinvsq           = rinv11*rinv11;  

                /* Coulomb interaction */
                vcoul            = qq*rinv11;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

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
                qq               = qqOH;           
                rinvsq           = rinv12*rinv12;  

                /* Coulomb interaction */
                vcoul            = qq*rinv12;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

                /* Calculate temporary vectorial force */
                tx               = fscal*dx12;     
                ty               = fscal*dy12;     
                tz               = fscal*dz12;     

                /* Increment i atom force */
                fix1             = fix1 + tx;      
                fiy1             = fiy1 + ty;      
                fiz1             = fiz1 + tz;      

                /* Decrement j atom force */
                fjx2             = faction[j3+3] - tx;
                fjy2             = faction[j3+4] - ty;
                fjz2             = faction[j3+5] - tz;

                /* Load parameters for j atom */
                qq               = qqOH;           
                rinvsq           = rinv13*rinv13;  

                /* Coulomb interaction */
                vcoul            = qq*rinv13;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

                /* Calculate temporary vectorial force */
                tx               = fscal*dx13;     
                ty               = fscal*dy13;     
                tz               = fscal*dz13;     

                /* Increment i atom force */
                fix1             = fix1 + tx;      
                fiy1             = fiy1 + ty;      
                fiz1             = fiz1 + tz;      

                /* Decrement j atom force */
                fjx3             = faction[j3+6] - tx;
                fjy3             = faction[j3+7] - ty;
                fjz3             = faction[j3+8] - tz;

                /* Load parameters for j atom */
                qq               = qqOH;           
                rinvsq           = rinv21*rinv21;  

                /* Coulomb interaction */
                vcoul            = qq*rinv21;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

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
                qq               = qqHH;           
                rinvsq           = rinv22*rinv22;  

                /* Coulomb interaction */
                vcoul            = qq*rinv22;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

                /* Calculate temporary vectorial force */
                tx               = fscal*dx22;     
                ty               = fscal*dy22;     
                tz               = fscal*dz22;     

                /* Increment i atom force */
                fix2             = fix2 + tx;      
                fiy2             = fiy2 + ty;      
                fiz2             = fiz2 + tz;      

                /* Decrement j atom force */
                fjx2             = fjx2 - tx;      
                fjy2             = fjy2 - ty;      
                fjz2             = fjz2 - tz;      

                /* Load parameters for j atom */
                qq               = qqHH;           
                rinvsq           = rinv23*rinv23;  

                /* Coulomb interaction */
                vcoul            = qq*rinv23;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

                /* Calculate temporary vectorial force */
                tx               = fscal*dx23;     
                ty               = fscal*dy23;     
                tz               = fscal*dz23;     

                /* Increment i atom force */
                fix2             = fix2 + tx;      
                fiy2             = fiy2 + ty;      
                fiz2             = fiz2 + tz;      

                /* Decrement j atom force */
                fjx3             = fjx3 - tx;      
                fjy3             = fjy3 - ty;      
                fjz3             = fjz3 - tz;      

                /* Load parameters for j atom */
                qq               = qqOH;           
                rinvsq           = rinv31*rinv31;  

                /* Coulomb interaction */
                vcoul            = qq*rinv31;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

                /* Calculate temporary vectorial force */
                tx               = fscal*dx31;     
                ty               = fscal*dy31;     
                tz               = fscal*dz31;     

                /* Increment i atom force */
                fix3             = fix3 + tx;      
                fiy3             = fiy3 + ty;      
                fiz3             = fiz3 + tz;      

                /* Decrement j atom force */
                faction[j3+0]    = fjx1 - tx;      
                faction[j3+1]    = fjy1 - ty;      
                faction[j3+2]    = fjz1 - tz;      

                /* Load parameters for j atom */
                qq               = qqHH;           
                rinvsq           = rinv32*rinv32;  

                /* Coulomb interaction */
                vcoul            = qq*rinv32;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

                /* Calculate temporary vectorial force */
                tx               = fscal*dx32;     
                ty               = fscal*dy32;     
                tz               = fscal*dz32;     

                /* Increment i atom force */
                fix3             = fix3 + tx;      
                fiy3             = fiy3 + ty;      
                fiz3             = fiz3 + tz;      

                /* Decrement j atom force */
                faction[j3+3]    = fjx2 - tx;      
                faction[j3+4]    = fjy2 - ty;      
                faction[j3+5]    = fjz2 - tz;      

                /* Load parameters for j atom */
                qq               = qqHH;           
                rinvsq           = rinv33*rinv33;  

                /* Coulomb interaction */
                vcoul            = qq*rinv33;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

                /* Calculate temporary vectorial force */
                tx               = fscal*dx33;     
                ty               = fscal*dy33;     
                tz               = fscal*dz33;     

                /* Increment i atom force */
                fix3             = fix3 + tx;      
                fiy3             = fiy3 + ty;      
                fiz3             = fiz3 + tz;      

                /* Decrement j atom force */
                faction[j3+6]    = fjx3 - tx;      
                faction[j3+7]    = fjy3 - ty;      
                faction[j3+8]    = fjz3 - tz;      

                /* Inner loop uses 234 flops/iteration */
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
            fshift[is3]      = fshift[is3]+fix1+fix2+fix3;
            fshift[is3+1]    = fshift[is3+1]+fiy1+fiy2+fiy3;
            fshift[is3+2]    = fshift[is3+2]+fiz1+fiz2+fiz3;

            /* Add potential energies to the group for this list */
            ggid             = gid[n];         
            Vc[ggid]         = Vc[ggid] + vctot;

            /* Increment number of inner iterations */
            ninner           = ninner + nj1 - nj0;

            /* Outer loop uses 28 flops/iteration */
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
 * Gromacs nonbonded kernel nb_kernel102nf
 * Coulomb interaction:     Normal Coulomb
 * VdW interaction:         Not calculated
 * water optimization:      pairs of SPC/TIP3P interactions
 * Calculate forces:        no
 */
void nb_kernel102nf(
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
    real          qq,vcoul,vctot;
    real          ix1,iy1,iz1;
    real          ix2,iy2,iz2;
    real          ix3,iy3,iz3;
    real          jx1,jy1,jz1;
    real          jx2,jy2,jz2;
    real          jx3,jy3,jz3;
    real          dx11,dy11,dz11,rsq11,rinv11;
    real          dx12,dy12,dz12,rsq12,rinv12;
    real          dx13,dy13,dz13,rsq13,rinv13;
    real          dx21,dy21,dz21,rsq21,rinv21;
    real          dx22,dy22,dz22,rsq22,rinv22;
    real          dx23,dy23,dz23,rsq23,rinv23;
    real          dx31,dy31,dz31,rsq31,rinv31;
    real          dx32,dy32,dz32,rsq32,rinv32;
    real          dx33,dy33,dz33,rsq33,rinv33;
    real          qO,qH,qqOO,qqOH,qqHH;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    

    /* Initialize water data */
    ii               = iinr[0];        
    qO               = charge[ii];     
    qH               = charge[ii+1];   
    qqOO             = facel*qO*qO;    
    qqOH             = facel*qO*qH;    
    qqHH             = facel*qH*qH;    


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

            /* Zero the potential energy for this list */
            vctot            = 0;              

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
                jx2              = pos[j3+3];      
                jy2              = pos[j3+4];      
                jz2              = pos[j3+5];      
                jx3              = pos[j3+6];      
                jy3              = pos[j3+7];      
                jz3              = pos[j3+8];      

                /* Calculate distance */
                dx11             = ix1 - jx1;      
                dy11             = iy1 - jy1;      
                dz11             = iz1 - jz1;      
                rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
                dx12             = ix1 - jx2;      
                dy12             = iy1 - jy2;      
                dz12             = iz1 - jz2;      
                rsq12            = dx12*dx12+dy12*dy12+dz12*dz12;
                dx13             = ix1 - jx3;      
                dy13             = iy1 - jy3;      
                dz13             = iz1 - jz3;      
                rsq13            = dx13*dx13+dy13*dy13+dz13*dz13;
                dx21             = ix2 - jx1;      
                dy21             = iy2 - jy1;      
                dz21             = iz2 - jz1;      
                rsq21            = dx21*dx21+dy21*dy21+dz21*dz21;
                dx22             = ix2 - jx2;      
                dy22             = iy2 - jy2;      
                dz22             = iz2 - jz2;      
                rsq22            = dx22*dx22+dy22*dy22+dz22*dz22;
                dx23             = ix2 - jx3;      
                dy23             = iy2 - jy3;      
                dz23             = iz2 - jz3;      
                rsq23            = dx23*dx23+dy23*dy23+dz23*dz23;
                dx31             = ix3 - jx1;      
                dy31             = iy3 - jy1;      
                dz31             = iz3 - jz1;      
                rsq31            = dx31*dx31+dy31*dy31+dz31*dz31;
                dx32             = ix3 - jx2;      
                dy32             = iy3 - jy2;      
                dz32             = iz3 - jz2;      
                rsq32            = dx32*dx32+dy32*dy32+dz32*dz32;
                dx33             = ix3 - jx3;      
                dy33             = iy3 - jy3;      
                dz33             = iz3 - jz3;      
                rsq33            = dx33*dx33+dy33*dy33+dz33*dz33;

                /* Calculate 1/r and 1/r2 */
                rinv11           = gmx_invsqrt(rsq11);
                rinv12           = gmx_invsqrt(rsq12);
                rinv13           = gmx_invsqrt(rsq13);
                rinv21           = gmx_invsqrt(rsq21);
                rinv22           = gmx_invsqrt(rsq22);
                rinv23           = gmx_invsqrt(rsq23);
                rinv31           = gmx_invsqrt(rsq31);
                rinv32           = gmx_invsqrt(rsq32);
                rinv33           = gmx_invsqrt(rsq33);

                /* Load parameters for j atom */
                qq               = qqOO;           

                /* Coulomb interaction */
                vcoul            = qq*rinv11;      
                vctot            = vctot+vcoul;    

                /* Load parameters for j atom */
                qq               = qqOH;           

                /* Coulomb interaction */
                vcoul            = qq*rinv12;      
                vctot            = vctot+vcoul;    

                /* Load parameters for j atom */
                qq               = qqOH;           

                /* Coulomb interaction */
                vcoul            = qq*rinv13;      
                vctot            = vctot+vcoul;    

                /* Load parameters for j atom */
                qq               = qqOH;           

                /* Coulomb interaction */
                vcoul            = qq*rinv21;      
                vctot            = vctot+vcoul;    

                /* Load parameters for j atom */
                qq               = qqHH;           

                /* Coulomb interaction */
                vcoul            = qq*rinv22;      
                vctot            = vctot+vcoul;    

                /* Load parameters for j atom */
                qq               = qqHH;           

                /* Coulomb interaction */
                vcoul            = qq*rinv23;      
                vctot            = vctot+vcoul;    

                /* Load parameters for j atom */
                qq               = qqOH;           

                /* Coulomb interaction */
                vcoul            = qq*rinv31;      
                vctot            = vctot+vcoul;    

                /* Load parameters for j atom */
                qq               = qqHH;           

                /* Coulomb interaction */
                vcoul            = qq*rinv32;      
                vctot            = vctot+vcoul;    

                /* Load parameters for j atom */
                qq               = qqHH;           

                /* Coulomb interaction */
                vcoul            = qq*rinv33;      
                vctot            = vctot+vcoul;    

                /* Inner loop uses 135 flops/iteration */
            }
            

            /* Add i forces to mem and shifted force list */

            /* Add potential energies to the group for this list */
            ggid             = gid[n];         
            Vc[ggid]         = Vc[ggid] + vctot;

            /* Increment number of inner iterations */
            ninner           = ninner + nj1 - nj0;

            /* Outer loop uses 10 flops/iteration */
        }
        

        /* Increment number of outer iterations */
        nouter           = nouter + nn1 - nn0;
    }
    while (nn1<nri);
    

    /* Write outer/inner iteration count to pointers */
    *outeriter       = nouter;         
    *inneriter       = ninner;         
}



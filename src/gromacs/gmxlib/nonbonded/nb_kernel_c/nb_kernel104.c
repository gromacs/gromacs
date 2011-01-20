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

#include "nb_kernel104.h"

/*
 * Gromacs nonbonded kernel nb_kernel104
 * Coulomb interaction:     Normal Coulomb
 * VdW interaction:         Not calculated
 * water optimization:      pairs of TIP4P interactions
 * Calculate forces:        yes
 */
void nb_kernel104(
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
    real          ix2,iy2,iz2,fix2,fiy2,fiz2;
    real          ix3,iy3,iz3,fix3,fiy3,fiz3;
    real          ix4,iy4,iz4,fix4,fiy4,fiz4;
    real          jx2,jy2,jz2,fjx2,fjy2,fjz2;
    real          jx3,jy3,jz3,fjx3,fjy3,fjz3;
    real          jx4,jy4,jz4,fjx4,fjy4,fjz4;
    real          dx22,dy22,dz22,rsq22,rinv22;
    real          dx23,dy23,dz23,rsq23,rinv23;
    real          dx24,dy24,dz24,rsq24,rinv24;
    real          dx32,dy32,dz32,rsq32,rinv32;
    real          dx33,dy33,dz33,rsq33,rinv33;
    real          dx34,dy34,dz34,rsq34,rinv34;
    real          dx42,dy42,dz42,rsq42,rinv42;
    real          dx43,dy43,dz43,rsq43,rinv43;
    real          dx44,dy44,dz44,rsq44,rinv44;
    real          qH,qM,qqMM,qqMH,qqHH;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    

    /* Initialize water data */
    ii               = iinr[0];        
    qH               = charge[ii+1];   
    qM               = charge[ii+3];   
    qqMM             = facel*qM*qM;    
    qqMH             = facel*qM*qH;    
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

            /* Clear i atom forces */
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
                jx2              = pos[j3+3];      
                jy2              = pos[j3+4];      
                jz2              = pos[j3+5];      
                jx3              = pos[j3+6];      
                jy3              = pos[j3+7];      
                jz3              = pos[j3+8];      
                jx4              = pos[j3+9];      
                jy4              = pos[j3+10];     
                jz4              = pos[j3+11];     

                /* Calculate distance */
                dx22             = ix2 - jx2;      
                dy22             = iy2 - jy2;      
                dz22             = iz2 - jz2;      
                rsq22            = dx22*dx22+dy22*dy22+dz22*dz22;
                dx23             = ix2 - jx3;      
                dy23             = iy2 - jy3;      
                dz23             = iz2 - jz3;      
                rsq23            = dx23*dx23+dy23*dy23+dz23*dz23;
                dx24             = ix2 - jx4;      
                dy24             = iy2 - jy4;      
                dz24             = iz2 - jz4;      
                rsq24            = dx24*dx24+dy24*dy24+dz24*dz24;
                dx32             = ix3 - jx2;      
                dy32             = iy3 - jy2;      
                dz32             = iz3 - jz2;      
                rsq32            = dx32*dx32+dy32*dy32+dz32*dz32;
                dx33             = ix3 - jx3;      
                dy33             = iy3 - jy3;      
                dz33             = iz3 - jz3;      
                rsq33            = dx33*dx33+dy33*dy33+dz33*dz33;
                dx34             = ix3 - jx4;      
                dy34             = iy3 - jy4;      
                dz34             = iz3 - jz4;      
                rsq34            = dx34*dx34+dy34*dy34+dz34*dz34;
                dx42             = ix4 - jx2;      
                dy42             = iy4 - jy2;      
                dz42             = iz4 - jz2;      
                rsq42            = dx42*dx42+dy42*dy42+dz42*dz42;
                dx43             = ix4 - jx3;      
                dy43             = iy4 - jy3;      
                dz43             = iz4 - jz3;      
                rsq43            = dx43*dx43+dy43*dy43+dz43*dz43;
                dx44             = ix4 - jx4;      
                dy44             = iy4 - jy4;      
                dz44             = iz4 - jz4;      
                rsq44            = dx44*dx44+dy44*dy44+dz44*dz44;

                /* Calculate 1/r and 1/r2 */
                rinv22           = gmx_invsqrt(rsq22);
                rinv23           = gmx_invsqrt(rsq23);
                rinv24           = gmx_invsqrt(rsq24);
                rinv32           = gmx_invsqrt(rsq32);
                rinv33           = gmx_invsqrt(rsq33);
                rinv34           = gmx_invsqrt(rsq34);
                rinv42           = gmx_invsqrt(rsq42);
                rinv43           = gmx_invsqrt(rsq43);
                rinv44           = gmx_invsqrt(rsq44);

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
                fjx2             = faction[j3+3] - tx;
                fjy2             = faction[j3+4] - ty;
                fjz2             = faction[j3+5] - tz;

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
                fjx3             = faction[j3+6] - tx;
                fjy3             = faction[j3+7] - ty;
                fjz3             = faction[j3+8] - tz;

                /* Load parameters for j atom */
                qq               = qqMH;           
                rinvsq           = rinv24*rinv24;  

                /* Coulomb interaction */
                vcoul            = qq*rinv24;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

                /* Calculate temporary vectorial force */
                tx               = fscal*dx24;     
                ty               = fscal*dy24;     
                tz               = fscal*dz24;     

                /* Increment i atom force */
                fix2             = fix2 + tx;      
                fiy2             = fiy2 + ty;      
                fiz2             = fiz2 + tz;      

                /* Decrement j atom force */
                fjx4             = faction[j3+9] - tx;
                fjy4             = faction[j3+10] - ty;
                fjz4             = faction[j3+11] - tz;

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
                fjx2             = fjx2 - tx;      
                fjy2             = fjy2 - ty;      
                fjz2             = fjz2 - tz;      

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
                fjx3             = fjx3 - tx;      
                fjy3             = fjy3 - ty;      
                fjz3             = fjz3 - tz;      

                /* Load parameters for j atom */
                qq               = qqMH;           
                rinvsq           = rinv34*rinv34;  

                /* Coulomb interaction */
                vcoul            = qq*rinv34;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

                /* Calculate temporary vectorial force */
                tx               = fscal*dx34;     
                ty               = fscal*dy34;     
                tz               = fscal*dz34;     

                /* Increment i atom force */
                fix3             = fix3 + tx;      
                fiy3             = fiy3 + ty;      
                fiz3             = fiz3 + tz;      

                /* Decrement j atom force */
                fjx4             = fjx4 - tx;      
                fjy4             = fjy4 - ty;      
                fjz4             = fjz4 - tz;      

                /* Load parameters for j atom */
                qq               = qqMH;           
                rinvsq           = rinv42*rinv42;  

                /* Coulomb interaction */
                vcoul            = qq*rinv42;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

                /* Calculate temporary vectorial force */
                tx               = fscal*dx42;     
                ty               = fscal*dy42;     
                tz               = fscal*dz42;     

                /* Increment i atom force */
                fix4             = fix4 + tx;      
                fiy4             = fiy4 + ty;      
                fiz4             = fiz4 + tz;      

                /* Decrement j atom force */
                faction[j3+3]    = fjx2 - tx;      
                faction[j3+4]    = fjy2 - ty;      
                faction[j3+5]    = fjz2 - tz;      

                /* Load parameters for j atom */
                qq               = qqMH;           
                rinvsq           = rinv43*rinv43;  

                /* Coulomb interaction */
                vcoul            = qq*rinv43;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

                /* Calculate temporary vectorial force */
                tx               = fscal*dx43;     
                ty               = fscal*dy43;     
                tz               = fscal*dz43;     

                /* Increment i atom force */
                fix4             = fix4 + tx;      
                fiy4             = fiy4 + ty;      
                fiz4             = fiz4 + tz;      

                /* Decrement j atom force */
                faction[j3+6]    = fjx3 - tx;      
                faction[j3+7]    = fjy3 - ty;      
                faction[j3+8]    = fjz3 - tz;      

                /* Load parameters for j atom */
                qq               = qqMM;           
                rinvsq           = rinv44*rinv44;  

                /* Coulomb interaction */
                vcoul            = qq*rinv44;      
                vctot            = vctot+vcoul;    
                fscal            = (vcoul)*rinvsq; 

                /* Calculate temporary vectorial force */
                tx               = fscal*dx44;     
                ty               = fscal*dy44;     
                tz               = fscal*dz44;     

                /* Increment i atom force */
                fix4             = fix4 + tx;      
                fiy4             = fiy4 + ty;      
                fiz4             = fiz4 + tz;      

                /* Decrement j atom force */
                faction[j3+9]    = fjx4 - tx;      
                faction[j3+10]   = fjy4 - ty;      
                faction[j3+11]   = fjz4 - tz;      

                /* Inner loop uses 234 flops/iteration */
            }
            

            /* Add i forces to mem and shifted force list */
            faction[ii3+3]   = faction[ii3+3] + fix2;
            faction[ii3+4]   = faction[ii3+4] + fiy2;
            faction[ii3+5]   = faction[ii3+5] + fiz2;
            faction[ii3+6]   = faction[ii3+6] + fix3;
            faction[ii3+7]   = faction[ii3+7] + fiy3;
            faction[ii3+8]   = faction[ii3+8] + fiz3;
            faction[ii3+9]   = faction[ii3+9] + fix4;
            faction[ii3+10]  = faction[ii3+10] + fiy4;
            faction[ii3+11]  = faction[ii3+11] + fiz4;
            fshift[is3]      = fshift[is3]+fix2+fix3+fix4;
            fshift[is3+1]    = fshift[is3+1]+fiy2+fiy3+fiy4;
            fshift[is3+2]    = fshift[is3+2]+fiz2+fiz3+fiz4;

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
 * Gromacs nonbonded kernel nb_kernel104nf
 * Coulomb interaction:     Normal Coulomb
 * VdW interaction:         Not calculated
 * water optimization:      pairs of TIP4P interactions
 * Calculate forces:        no
 */
void nb_kernel104nf(
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
    real          ix2,iy2,iz2;
    real          ix3,iy3,iz3;
    real          ix4,iy4,iz4;
    real          jx2,jy2,jz2;
    real          jx3,jy3,jz3;
    real          jx4,jy4,jz4;
    real          dx22,dy22,dz22,rsq22,rinv22;
    real          dx23,dy23,dz23,rsq23,rinv23;
    real          dx24,dy24,dz24,rsq24,rinv24;
    real          dx32,dy32,dz32,rsq32,rinv32;
    real          dx33,dy33,dz33,rsq33,rinv33;
    real          dx34,dy34,dz34,rsq34,rinv34;
    real          dx42,dy42,dz42,rsq42,rinv42;
    real          dx43,dy43,dz43,rsq43,rinv43;
    real          dx44,dy44,dz44,rsq44,rinv44;
    real          qH,qM,qqMM,qqMH,qqHH;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    

    /* Initialize water data */
    ii               = iinr[0];        
    qH               = charge[ii+1];   
    qM               = charge[ii+3];   
    qqMM             = facel*qM*qM;    
    qqMH             = facel*qM*qH;    
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

            /* Clear i atom forces */
            
            for(k=nj0; (k<nj1); k++)
            {

                /* Get j neighbor index, and coordinate index */
                jnr              = jjnr[k];        
                j3               = 3*jnr;          

                /* load j atom coordinates */
                jx2              = pos[j3+3];      
                jy2              = pos[j3+4];      
                jz2              = pos[j3+5];      
                jx3              = pos[j3+6];      
                jy3              = pos[j3+7];      
                jz3              = pos[j3+8];      
                jx4              = pos[j3+9];      
                jy4              = pos[j3+10];     
                jz4              = pos[j3+11];     

                /* Calculate distance */
                dx22             = ix2 - jx2;      
                dy22             = iy2 - jy2;      
                dz22             = iz2 - jz2;      
                rsq22            = dx22*dx22+dy22*dy22+dz22*dz22;
                dx23             = ix2 - jx3;      
                dy23             = iy2 - jy3;      
                dz23             = iz2 - jz3;      
                rsq23            = dx23*dx23+dy23*dy23+dz23*dz23;
                dx24             = ix2 - jx4;      
                dy24             = iy2 - jy4;      
                dz24             = iz2 - jz4;      
                rsq24            = dx24*dx24+dy24*dy24+dz24*dz24;
                dx32             = ix3 - jx2;      
                dy32             = iy3 - jy2;      
                dz32             = iz3 - jz2;      
                rsq32            = dx32*dx32+dy32*dy32+dz32*dz32;
                dx33             = ix3 - jx3;      
                dy33             = iy3 - jy3;      
                dz33             = iz3 - jz3;      
                rsq33            = dx33*dx33+dy33*dy33+dz33*dz33;
                dx34             = ix3 - jx4;      
                dy34             = iy3 - jy4;      
                dz34             = iz3 - jz4;      
                rsq34            = dx34*dx34+dy34*dy34+dz34*dz34;
                dx42             = ix4 - jx2;      
                dy42             = iy4 - jy2;      
                dz42             = iz4 - jz2;      
                rsq42            = dx42*dx42+dy42*dy42+dz42*dz42;
                dx43             = ix4 - jx3;      
                dy43             = iy4 - jy3;      
                dz43             = iz4 - jz3;      
                rsq43            = dx43*dx43+dy43*dy43+dz43*dz43;
                dx44             = ix4 - jx4;      
                dy44             = iy4 - jy4;      
                dz44             = iz4 - jz4;      
                rsq44            = dx44*dx44+dy44*dy44+dz44*dz44;

                /* Calculate 1/r and 1/r2 */
                rinv22           = gmx_invsqrt(rsq22);
                rinv23           = gmx_invsqrt(rsq23);
                rinv24           = gmx_invsqrt(rsq24);
                rinv32           = gmx_invsqrt(rsq32);
                rinv33           = gmx_invsqrt(rsq33);
                rinv34           = gmx_invsqrt(rsq34);
                rinv42           = gmx_invsqrt(rsq42);
                rinv43           = gmx_invsqrt(rsq43);
                rinv44           = gmx_invsqrt(rsq44);

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
                qq               = qqMH;           

                /* Coulomb interaction */
                vcoul            = qq*rinv24;      
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

                /* Load parameters for j atom */
                qq               = qqMH;           

                /* Coulomb interaction */
                vcoul            = qq*rinv34;      
                vctot            = vctot+vcoul;    

                /* Load parameters for j atom */
                qq               = qqMH;           

                /* Coulomb interaction */
                vcoul            = qq*rinv42;      
                vctot            = vctot+vcoul;    

                /* Load parameters for j atom */
                qq               = qqMH;           

                /* Coulomb interaction */
                vcoul            = qq*rinv43;      
                vctot            = vctot+vcoul;    

                /* Load parameters for j atom */
                qq               = qqMM;           

                /* Coulomb interaction */
                vcoul            = qq*rinv44;      
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



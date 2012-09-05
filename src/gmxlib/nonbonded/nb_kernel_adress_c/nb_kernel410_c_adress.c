/*
 * Copyright (c) Erik Lindahl, David van der Spoel 2003
 * 
 * This file is generated automatically at compile time
 * by the program mknb in the Gromacs distribution.
 *
 * Options used when generation this file:
 * Language:         c
 * Precision:        single
 * Threads:          yes
 * Software invsqrt: no
 * PowerPC invsqrt:  no
 * Prefetch forces:  no
 * Adress kernel:  yes
 * Comments:         no
 */
#ifdef HAVE_CONFIG_H
#include<config.h>
#endif
#ifdef GMX_THREAD_SHM_FDECOMP
#include<thread_mpi.h>
#endif
#define ALMOST_ZERO 1e-30
#define ALMOST_ONE 1-(1e-30)
#include<math.h>

#include "nb_kernel410_adress.h"



/*
 * Gromacs nonbonded kernel nb_kernel410_adress_cg
 * Coulomb interaction:     Generalized-Born
 * VdW interaction:         Lennard-Jones
 * water optimization:      No
 * Calculate forces:        yes
 */
void nb_kernel410_adress_cg(
                    int *           p_nri,
                    int *           iinr,
                    int *           jindex,
                    int *           jjnr,
                    int *           shift,
                    real *         shiftvec,
                    real *         fshift,
                    int *           gid,
                    real *         pos,
                    real *         faction,
                    real *         charge,
                    real *         p_facel,
                    real *         p_krf,
                    real *         p_crf,
                    real *         Vc,
                    int *           type,
                    int *           p_ntype,
                    real *         vdwparam,
                    real *         Vvdw,
                    real *         p_tabscale,
                    real *         VFtab,
                    real *         invsqrta,
                    real *         dvda,
                    real *         p_gbtabscale,
                    real *         GBtab,
                    int *           p_nthreads,
                    int *           count,
                    void *          mtx,
                    int *           outeriter,
                    int *           inneriter,
                    real           force_cap,
                    real *         wf)
{
    int           nri,ntype,nthreads;
    real         facel,krf,crf,tabscale,gbtabscale;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    int           nn0,nn1,nouter,ninner;
    real         shX,shY,shZ;
    real         fscal,tx,ty,tz;
    real         rinvsq;
    real         iq;
    real         qq,vcoul,vctot;
    int           nti;
    int           tj;
    real         rinvsix;
    real         Vvdw6,Vvdwtot;
    real         Vvdw12;
    real         r,rt,eps,eps2;
    int           n0,nnn;
    real         Y,F,Geps,Heps2,Fp,VV;
    real         FF;
    real         fijC;
    real         isai,isaj,isaprod,gbscale,vgb;
    real         dvdasum,dvdatmp,dvdaj,fgb;
    real         ix1,iy1,iz1,fix1,fiy1,fiz1;
    real         jx1,jy1,jz1;
    real         dx11,dy11,dz11,rsq11,rinv11;
    real         c6,c12;
    real         weight_cg1, weight_cg2, weight_product;
    real         hybscal;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    gbtabscale       = *p_gbtabscale;  
    nouter           = 0;              
    ninner           = 0;              
    
    do
    {
        #ifdef GMX_THREAD_SHM_FDECOMP
        tMPI_Thread_mutex_lock((tMPI_Thread_mutex_t *)mtx);
        nn0              = *count;         
        nn1              = nn0+(nri-nn0)/(2*nthreads)+10;
        *count           = nn1;            
        tMPI_Thread_mutex_unlock((tMPI_Thread_mutex_t *)mtx);
        if(nn1>nri) nn1=nri;
        #else
        nn0 = 0;
        nn1 = nri;
        #endif
        
        for(n=nn0; (n<nn1); n++)
        {
            is3              = 3*shift[n];     
            shX              = shiftvec[is3];  
            shY              = shiftvec[is3+1];
            shZ              = shiftvec[is3+2];
            nj0              = jindex[n];      
            nj1              = jindex[n+1];    
            ii               = iinr[n];        
            ii3              = 3*ii;           
            ix1              = shX + pos[ii3+0];
            iy1              = shY + pos[ii3+1];
            iz1              = shZ + pos[ii3+2];
            iq               = facel*charge[ii];
            isai             = invsqrta[ii];   
            nti              = 2*ntype*type[ii];
            weight_cg1       = wf[ii];         
            vctot            = 0;              
            Vvdwtot          = 0;              
            dvdasum          = 0;              
            fix1             = 0;              
            fiy1             = 0;              
            fiz1             = 0;              
            
            for(k=nj0; (k<nj1); k++)
            {
                jnr              = jjnr[k];        
                weight_cg2       = wf[jnr];        
                weight_product   = weight_cg1*weight_cg2;
                if (weight_product < ALMOST_ZERO) {
                       hybscal = 1.0;
                }
                else if (weight_product >= ALMOST_ONE)
                {
                  /* force is zero, skip this molecule */
                       continue;
                }
                else
                {
                   hybscal = 1.0 - weight_product;
                }
                j3               = 3*jnr;          
                jx1              = pos[j3+0];      
                jy1              = pos[j3+1];      
                jz1              = pos[j3+2];      
                dx11             = ix1 - jx1;      
                dy11             = iy1 - jy1;      
                dz11             = iz1 - jz1;      
                rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
                rinv11           = 1.0/sqrt(rsq11);
                isaj             = invsqrta[jnr];  
                isaprod          = isai*isaj;      
                qq               = iq*charge[jnr]; 
                vcoul            = qq*rinv11;      
                fscal            = vcoul*rinv11;   
                qq               = isaprod*(-qq);  
                gbscale          = isaprod*gbtabscale;
                tj               = nti+2*type[jnr];
                c6               = vdwparam[tj];   
                c12              = vdwparam[tj+1]; 
                rinvsq           = rinv11*rinv11;  
                dvdaj            = dvda[jnr];      
                r                = rsq11*rinv11;   
                rt               = r*gbscale;      
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 4*n0;           
                Y                = GBtab[nnn];     
                F                = GBtab[nnn+1];   
                Geps             = eps*GBtab[nnn+2];
                Heps2            = eps2*GBtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                FF               = Fp+Geps+2.0*Heps2;
                vgb              = qq*VV;          
                fijC             = qq*FF*gbscale;  
                dvdatmp          = -0.5*(vgb+fijC*r);
                dvdasum          = dvdasum + dvdatmp;
                dvda[jnr]        = dvdaj+dvdatmp*isaj*isaj;
                vctot            = vctot + vcoul;  
                rinvsix          = rinvsq*rinvsq*rinvsq;
                Vvdw6            = c6*rinvsix;     
                Vvdw12           = c12*rinvsix*rinvsix;
                Vvdwtot          = Vvdwtot+Vvdw12-Vvdw6;
                fscal            = (12.0*Vvdw12-6.0*Vvdw6)*rinvsq-(fijC-fscal)*rinv11;
                fscal *= hybscal;
                tx               = fscal*dx11;     
                ty               = fscal*dy11;     
                tz               = fscal*dz11;     
                fix1             = fix1 + tx;      
                fiy1             = fiy1 + ty;      
                fiz1             = fiz1 + tz;      
                faction[j3+0]    = faction[j3+0] - tx;
                faction[j3+1]    = faction[j3+1] - ty;
                faction[j3+2]    = faction[j3+2] - tz;
            }
            
            faction[ii3+0]   = faction[ii3+0] + fix1;
            faction[ii3+1]   = faction[ii3+1] + fiy1;
            faction[ii3+2]   = faction[ii3+2] + fiz1;
            fshift[is3]      = fshift[is3]+fix1;
            fshift[is3+1]    = fshift[is3+1]+fiy1;
            fshift[is3+2]    = fshift[is3+2]+fiz1;
            ggid             = gid[n];         
            Vc[ggid]         = Vc[ggid] + vctot;
            Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;
            dvda[ii]         = dvda[ii] + dvdasum*isai*isai;
            ninner           = ninner + nj1 - nj0;
        }
        
        nouter           = nouter + nn1 - nn0;
    }
    while (nn1<nri);
    
    *outeriter       = nouter;         
    *inneriter       = ninner;         
}





/*
 * Gromacs nonbonded kernel nb_kernel410_adress_ex
 * Coulomb interaction:     Generalized-Born
 * VdW interaction:         Lennard-Jones
 * water optimization:      No
 * Calculate forces:        yes
 */
void nb_kernel410_adress_ex(
                    int *           p_nri,
                    int *           iinr,
                    int *           jindex,
                    int *           jjnr,
                    int *           shift,
                    real *         shiftvec,
                    real *         fshift,
                    int *           gid,
                    real *         pos,
                    real *         faction,
                    real *         charge,
                    real *         p_facel,
                    real *         p_krf,
                    real *         p_crf,
                    real *         Vc,
                    int *           type,
                    int *           p_ntype,
                    real *         vdwparam,
                    real *         Vvdw,
                    real *         p_tabscale,
                    real *         VFtab,
                    real *         invsqrta,
                    real *         dvda,
                    real *         p_gbtabscale,
                    real *         GBtab,
                    int *           p_nthreads,
                    int *           count,
                    void *          mtx,
                    int *           outeriter,
                    int *           inneriter,
                    real           force_cap,
                    real *         wf)
{
    int           nri,ntype,nthreads;
    real         facel,krf,crf,tabscale,gbtabscale;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    int           nn0,nn1,nouter,ninner;
    real         shX,shY,shZ;
    real         fscal,tx,ty,tz;
    real         rinvsq;
    real         iq;
    real         qq,vcoul,vctot;
    int           nti;
    int           tj;
    real         rinvsix;
    real         Vvdw6,Vvdwtot;
    real         Vvdw12;
    real         r,rt,eps,eps2;
    int           n0,nnn;
    real         Y,F,Geps,Heps2,Fp,VV;
    real         FF;
    real         fijC;
    real         isai,isaj,isaprod,gbscale,vgb;
    real         dvdasum,dvdatmp,dvdaj,fgb;
    real         ix1,iy1,iz1,fix1,fiy1,fiz1;
    real         jx1,jy1,jz1;
    real         dx11,dy11,dz11,rsq11,rinv11;
    real         c6,c12;
    real         weight_cg1, weight_cg2, weight_product;
    real         hybscal;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    gbtabscale       = *p_gbtabscale;  
    nouter           = 0;              
    ninner           = 0;              
    
    do
    {
        #ifdef GMX_THREAD_SHM_FDECOMP
        tMPI_Thread_mutex_lock((tMPI_Thread_mutex_t *)mtx);
        nn0              = *count;         
        nn1              = nn0+(nri-nn0)/(2*nthreads)+10;
        *count           = nn1;            
        tMPI_Thread_mutex_unlock((tMPI_Thread_mutex_t *)mtx);
        if(nn1>nri) nn1=nri;
        #else
        nn0 = 0;
        nn1 = nri;
        #endif
        
        for(n=nn0; (n<nn1); n++)
        {
            is3              = 3*shift[n];     
            shX              = shiftvec[is3];  
            shY              = shiftvec[is3+1];
            shZ              = shiftvec[is3+2];
            nj0              = jindex[n];      
            nj1              = jindex[n+1];    
            ii               = iinr[n];        
            ii3              = 3*ii;           
            ix1              = shX + pos[ii3+0];
            iy1              = shY + pos[ii3+1];
            iz1              = shZ + pos[ii3+2];
            iq               = facel*charge[ii];
            isai             = invsqrta[ii];   
            nti              = 2*ntype*type[ii];
            weight_cg1       = wf[ii];         
            vctot            = 0;              
            Vvdwtot          = 0;              
            dvdasum          = 0;              
            fix1             = 0;              
            fiy1             = 0;              
            fiz1             = 0;              
            
            for(k=nj0; (k<nj1); k++)
            {
                jnr              = jjnr[k];        
                weight_cg2       = wf[jnr];        
                weight_product   = weight_cg1*weight_cg2;
                if (weight_product < ALMOST_ZERO) {
                /* force is zero, skip this molecule */
                 continue;
                }
                else if (weight_product >= ALMOST_ONE)
                {
                       hybscal = 1.0;
                }
                else
                {
                   hybscal = weight_product;
                }
                j3               = 3*jnr;          
                jx1              = pos[j3+0];      
                jy1              = pos[j3+1];      
                jz1              = pos[j3+2];      
                dx11             = ix1 - jx1;      
                dy11             = iy1 - jy1;      
                dz11             = iz1 - jz1;      
                rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
                rinv11           = 1.0/sqrt(rsq11);
                isaj             = invsqrta[jnr];  
                isaprod          = isai*isaj;      
                qq               = iq*charge[jnr]; 
                vcoul            = qq*rinv11;      
                fscal            = vcoul*rinv11;   
                qq               = isaprod*(-qq);  
                gbscale          = isaprod*gbtabscale;
                tj               = nti+2*type[jnr];
                c6               = vdwparam[tj];   
                c12              = vdwparam[tj+1]; 
                rinvsq           = rinv11*rinv11;  
                dvdaj            = dvda[jnr];      
                r                = rsq11*rinv11;   
                rt               = r*gbscale;      
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 4*n0;           
                Y                = GBtab[nnn];     
                F                = GBtab[nnn+1];   
                Geps             = eps*GBtab[nnn+2];
                Heps2            = eps2*GBtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                FF               = Fp+Geps+2.0*Heps2;
                vgb              = qq*VV;          
                fijC             = qq*FF*gbscale;  
                dvdatmp          = -0.5*(vgb+fijC*r);
                dvdasum          = dvdasum + dvdatmp;
                dvda[jnr]        = dvdaj+dvdatmp*isaj*isaj;
                vctot            = vctot + vcoul;  
                rinvsix          = rinvsq*rinvsq*rinvsq;
                Vvdw6            = c6*rinvsix;     
                Vvdw12           = c12*rinvsix*rinvsix;
                Vvdwtot          = Vvdwtot+Vvdw12-Vvdw6;
                fscal            = (12.0*Vvdw12-6.0*Vvdw6)*rinvsq-(fijC-fscal)*rinv11;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx11;     
                ty               = fscal*dy11;     
                tz               = fscal*dz11;     
                fix1             = fix1 + tx;      
                fiy1             = fiy1 + ty;      
                fiz1             = fiz1 + tz;      
                faction[j3+0]    = faction[j3+0] - tx;
                faction[j3+1]    = faction[j3+1] - ty;
                faction[j3+2]    = faction[j3+2] - tz;
            }
            
            faction[ii3+0]   = faction[ii3+0] + fix1;
            faction[ii3+1]   = faction[ii3+1] + fiy1;
            faction[ii3+2]   = faction[ii3+2] + fiz1;
            fshift[is3]      = fshift[is3]+fix1;
            fshift[is3+1]    = fshift[is3+1]+fiy1;
            fshift[is3+2]    = fshift[is3+2]+fiz1;
            ggid             = gid[n];         
            Vc[ggid]         = Vc[ggid] + vctot;
            Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;
            dvda[ii]         = dvda[ii] + dvdasum*isai*isai;
            ninner           = ninner + nj1 - nj0;
        }
        
        nouter           = nouter + nn1 - nn0;
    }
    while (nn1<nri);
    
    *outeriter       = nouter;         
    *inneriter       = ninner;         
}



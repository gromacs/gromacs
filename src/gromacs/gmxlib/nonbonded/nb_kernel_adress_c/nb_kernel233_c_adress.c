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



/*
 * Gromacs nonbonded kernel nb_kernel233_adress_cg
 * Coulomb interaction:     Reaction field
 * VdW interaction:         Tabulated
 * water optimization:      TIP4P - other atoms
 * Calculate forces:        yes
 */
void nb_kernel233_adress_cg(
                    int *           p_nri,
                    int *           iinr,
                    int *           jindex,
                    int *           jjnr,
                    int *           shift,
                    float *         shiftvec,
                    float *         fshift,
                    int *           gid,
                    float *         pos,
                    float *         faction,
                    float *         charge,
                    float *         p_facel,
                    float *         p_krf,
                    float *         p_crf,
                    float *         Vc,
                    int *           type,
                    int *           p_ntype,
                    float *         vdwparam,
                    float *         Vvdw,
                    float *         p_tabscale,
                    float *         VFtab,
                    float *         invsqrta,
                    float *         dvda,
                    float *         p_gbtabscale,
                    float *         GBtab,
                    int *           p_nthreads,
                    int *           count,
                    void *          mtx,
                    int *           outeriter,
                    int *           inneriter,
                    float           force_cap,
                    float *         wf)
{
    int           nri,ntype,nthreads;
    float         facel,krf,crf,tabscale,gbtabscale;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    int           nn0,nn1,nouter,ninner;
    float         shX,shY,shZ;
    float         fscal,tx,ty,tz;
    float         rinvsq;
    float         jq;
    float         qq,vcoul,vctot;
    int           nti;
    int           tj;
    float         Vvdw6,Vvdwtot;
    float         Vvdw12;
    float         r,rt,eps,eps2;
    int           n0,nnn;
    float         Y,F,Geps,Heps2,Fp,VV;
    float         FF;
    float         fijD,fijR;
    float         krsq;
    float         ix1,iy1,iz1,fix1,fiy1,fiz1;
    float         ix2,iy2,iz2,fix2,fiy2,fiz2;
    float         ix3,iy3,iz3,fix3,fiy3,fiz3;
    float         ix4,iy4,iz4,fix4,fiy4,fiz4;
    float         jx1,jy1,jz1,fjx1,fjy1,fjz1;
    float         dx11,dy11,dz11,rsq11,rinv11;
    float         dx21,dy21,dz21,rsq21,rinv21;
    float         dx31,dy31,dz31,rsq31,rinv31;
    float         dx41,dy41,dz41,rsq41,rinv41;
    float         qH,qM;
    float         c6,c12;
    float         weight_cg1, weight_cg2, weight_product;
    float         hybscal;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    ii               = iinr[0];        
    qH               = facel*charge[ii+1];
    qM               = facel*charge[ii+3];
    nti              = 2*ntype*type[ii];

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
            ix2              = shX + pos[ii3+3];
            iy2              = shY + pos[ii3+4];
            iz2              = shZ + pos[ii3+5];
            ix3              = shX + pos[ii3+6];
            iy3              = shY + pos[ii3+7];
            iz3              = shZ + pos[ii3+8];
            ix4              = shX + pos[ii3+9];
            iy4              = shY + pos[ii3+10];
            iz4              = shZ + pos[ii3+11];
            weight_cg1       = wf[ii];         
            vctot            = 0;              
            Vvdwtot          = 0;              
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
                rinv11           = 1.0/sqrt(rsq11);
                rinv21           = 1.0/sqrt(rsq21);
                rinv31           = 1.0/sqrt(rsq31);
                rinv41           = 1.0/sqrt(rsq41);
                tj               = nti+2*type[jnr];
                c6               = vdwparam[tj];   
                c12              = vdwparam[tj+1]; 
                r                = rsq11*rinv11;   
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 8*n0;           
                Y                = VFtab[nnn];     
                F                = VFtab[nnn+1];   
                Geps             = eps*VFtab[nnn+2];
                Heps2            = eps2*VFtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                FF               = Fp+Geps+2.0*Heps2;
                Vvdw6            = c6*VV;          
                fijD             = c6*FF;          
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
                fscal *= hybscal;
                tx               = fscal*dx11;     
                ty               = fscal*dy11;     
                tz               = fscal*dz11;     
                fix1             = fix1 + tx;      
                fiy1             = fiy1 + ty;      
                fiz1             = fiz1 + tz;      
                fjx1             = faction[j3+0] - tx;
                fjy1             = faction[j3+1] - ty;
                fjz1             = faction[j3+2] - tz;
                jq               = charge[jnr+0];  
                qq               = qH*jq;          
                rinvsq           = rinv21*rinv21;  
                krsq             = krf*rsq21;      
                vcoul            = qq*(rinv21+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv21-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx21;     
                ty               = fscal*dy21;     
                tz               = fscal*dz21;     
                fix2             = fix2 + tx;      
                fiy2             = fiy2 + ty;      
                fiz2             = fiz2 + tz;      
                fjx1             = fjx1 - tx;      
                fjy1             = fjy1 - ty;      
                fjz1             = fjz1 - tz;      
                rinvsq           = rinv31*rinv31;  
                krsq             = krf*rsq31;      
                vcoul            = qq*(rinv31+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv31-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx31;     
                ty               = fscal*dy31;     
                tz               = fscal*dz31;     
                fix3             = fix3 + tx;      
                fiy3             = fiy3 + ty;      
                fiz3             = fiz3 + tz;      
                fjx1             = fjx1 - tx;      
                fjy1             = fjy1 - ty;      
                fjz1             = fjz1 - tz;      
                qq               = qM*jq;          
                rinvsq           = rinv41*rinv41;  
                krsq             = krf*rsq41;      
                vcoul            = qq*(rinv41+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv41-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx41;     
                ty               = fscal*dy41;     
                tz               = fscal*dz41;     
                fix4             = fix4 + tx;      
                fiy4             = fiy4 + ty;      
                fiz4             = fiz4 + tz;      
                faction[j3+0]    = fjx1 - tx;      
                faction[j3+1]    = fjy1 - ty;      
                faction[j3+2]    = fjz1 - tz;      
            }
            
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
            ggid             = gid[n];         
            Vc[ggid]         = Vc[ggid] + vctot;
            Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;
            ninner           = ninner + nj1 - nj0;
        }
        
        nouter           = nouter + nn1 - nn0;
    }
    while (nn1<nri);
    
    *outeriter       = nouter;         
    *inneriter       = ninner;         
}





/*
 * Gromacs nonbonded kernel nb_kernel233_adress_ex
 * Coulomb interaction:     Reaction field
 * VdW interaction:         Tabulated
 * water optimization:      TIP4P - other atoms
 * Calculate forces:        yes
 */
void nb_kernel233_adress_ex(
                    int *           p_nri,
                    int *           iinr,
                    int *           jindex,
                    int *           jjnr,
                    int *           shift,
                    float *         shiftvec,
                    float *         fshift,
                    int *           gid,
                    float *         pos,
                    float *         faction,
                    float *         charge,
                    float *         p_facel,
                    float *         p_krf,
                    float *         p_crf,
                    float *         Vc,
                    int *           type,
                    int *           p_ntype,
                    float *         vdwparam,
                    float *         Vvdw,
                    float *         p_tabscale,
                    float *         VFtab,
                    float *         invsqrta,
                    float *         dvda,
                    float *         p_gbtabscale,
                    float *         GBtab,
                    int *           p_nthreads,
                    int *           count,
                    void *          mtx,
                    int *           outeriter,
                    int *           inneriter,
                    float           force_cap,
                    float *         wf)
{
    int           nri,ntype,nthreads;
    float         facel,krf,crf,tabscale,gbtabscale;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    int           nn0,nn1,nouter,ninner;
    float         shX,shY,shZ;
    float         fscal,tx,ty,tz;
    float         rinvsq;
    float         jq;
    float         qq,vcoul,vctot;
    int           nti;
    int           tj;
    float         Vvdw6,Vvdwtot;
    float         Vvdw12;
    float         r,rt,eps,eps2;
    int           n0,nnn;
    float         Y,F,Geps,Heps2,Fp,VV;
    float         FF;
    float         fijD,fijR;
    float         krsq;
    float         ix1,iy1,iz1,fix1,fiy1,fiz1;
    float         ix2,iy2,iz2,fix2,fiy2,fiz2;
    float         ix3,iy3,iz3,fix3,fiy3,fiz3;
    float         ix4,iy4,iz4,fix4,fiy4,fiz4;
    float         jx1,jy1,jz1,fjx1,fjy1,fjz1;
    float         dx11,dy11,dz11,rsq11,rinv11;
    float         dx21,dy21,dz21,rsq21,rinv21;
    float         dx31,dy31,dz31,rsq31,rinv31;
    float         dx41,dy41,dz41,rsq41,rinv41;
    float         qH,qM;
    float         c6,c12;
    float         weight_cg1, weight_cg2, weight_product;
    float         hybscal;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    ii               = iinr[0];        
    qH               = facel*charge[ii+1];
    qM               = facel*charge[ii+3];
    nti              = 2*ntype*type[ii];

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
            ix2              = shX + pos[ii3+3];
            iy2              = shY + pos[ii3+4];
            iz2              = shZ + pos[ii3+5];
            ix3              = shX + pos[ii3+6];
            iy3              = shY + pos[ii3+7];
            iz3              = shZ + pos[ii3+8];
            ix4              = shX + pos[ii3+9];
            iy4              = shY + pos[ii3+10];
            iz4              = shZ + pos[ii3+11];
            weight_cg1       = wf[ii];         
            vctot            = 0;              
            Vvdwtot          = 0;              
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
                rinv11           = 1.0/sqrt(rsq11);
                rinv21           = 1.0/sqrt(rsq21);
                rinv31           = 1.0/sqrt(rsq31);
                rinv41           = 1.0/sqrt(rsq41);
                tj               = nti+2*type[jnr];
                c6               = vdwparam[tj];   
                c12              = vdwparam[tj+1]; 
                r                = rsq11*rinv11;   
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 8*n0;           
                Y                = VFtab[nnn];     
                F                = VFtab[nnn+1];   
                Geps             = eps*VFtab[nnn+2];
                Heps2            = eps2*VFtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                FF               = Fp+Geps+2.0*Heps2;
                Vvdw6            = c6*VV;          
                fijD             = c6*FF;          
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
                fjx1             = faction[j3+0] - tx;
                fjy1             = faction[j3+1] - ty;
                fjz1             = faction[j3+2] - tz;
                jq               = charge[jnr+0];  
                qq               = qH*jq;          
                rinvsq           = rinv21*rinv21;  
                krsq             = krf*rsq21;      
                vcoul            = qq*(rinv21+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv21-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx21;     
                ty               = fscal*dy21;     
                tz               = fscal*dz21;     
                fix2             = fix2 + tx;      
                fiy2             = fiy2 + ty;      
                fiz2             = fiz2 + tz;      
                fjx1             = fjx1 - tx;      
                fjy1             = fjy1 - ty;      
                fjz1             = fjz1 - tz;      
                rinvsq           = rinv31*rinv31;  
                krsq             = krf*rsq31;      
                vcoul            = qq*(rinv31+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv31-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx31;     
                ty               = fscal*dy31;     
                tz               = fscal*dz31;     
                fix3             = fix3 + tx;      
                fiy3             = fiy3 + ty;      
                fiz3             = fiz3 + tz;      
                fjx1             = fjx1 - tx;      
                fjy1             = fjy1 - ty;      
                fjz1             = fjz1 - tz;      
                qq               = qM*jq;          
                rinvsq           = rinv41*rinv41;  
                krsq             = krf*rsq41;      
                vcoul            = qq*(rinv41+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv41-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx41;     
                ty               = fscal*dy41;     
                tz               = fscal*dz41;     
                fix4             = fix4 + tx;      
                fiy4             = fiy4 + ty;      
                fiz4             = fiz4 + tz;      
                faction[j3+0]    = fjx1 - tx;      
                faction[j3+1]    = fjy1 - ty;      
                faction[j3+2]    = fjz1 - tz;      
            }
            
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
            ggid             = gid[n];         
            Vc[ggid]         = Vc[ggid] + vctot;
            Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;
            ninner           = ninner + nj1 - nj0;
        }
        
        nouter           = nouter + nn1 - nn0;
    }
    while (nn1<nri);
    
    *outeriter       = nouter;         
    *inneriter       = ninner;         
}



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
 * Gromacs nonbonded kernel nb_kernel212_adress_cg
 * Coulomb interaction:     Reaction field
 * VdW interaction:         Lennard-Jones
 * water optimization:      pairs of SPC/TIP3P interactions
 * Calculate forces:        yes
 */
void nb_kernel212_adress_cg(
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
    float         qq,vcoul,vctot;
    int           tj;
    float         rinvsix;
    float         Vvdw6,Vvdwtot;
    float         Vvdw12;
    float         krsq;
    float         ix1,iy1,iz1,fix1,fiy1,fiz1;
    float         ix2,iy2,iz2,fix2,fiy2,fiz2;
    float         ix3,iy3,iz3,fix3,fiy3,fiz3;
    float         jx1,jy1,jz1,fjx1,fjy1,fjz1;
    float         jx2,jy2,jz2,fjx2,fjy2,fjz2;
    float         jx3,jy3,jz3,fjx3,fjy3,fjz3;
    float         dx11,dy11,dz11,rsq11,rinv11;
    float         dx12,dy12,dz12,rsq12,rinv12;
    float         dx13,dy13,dz13,rsq13,rinv13;
    float         dx21,dy21,dz21,rsq21,rinv21;
    float         dx22,dy22,dz22,rsq22,rinv22;
    float         dx23,dy23,dz23,rsq23,rinv23;
    float         dx31,dy31,dz31,rsq31,rinv31;
    float         dx32,dy32,dz32,rsq32,rinv32;
    float         dx33,dy33,dz33,rsq33,rinv33;
    float         qO,qH,qqOO,qqOH,qqHH;
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
    qO               = charge[ii];     
    qH               = charge[ii+1];   
    qqOO             = facel*qO*qO;    
    qqOH             = facel*qO*qH;    
    qqHH             = facel*qH*qH;    
    tj               = 2*(ntype+1)*type[ii];
    c6               = vdwparam[tj];   
    c12              = vdwparam[tj+1]; 

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
                jx2              = pos[j3+3];      
                jy2              = pos[j3+4];      
                jz2              = pos[j3+5];      
                jx3              = pos[j3+6];      
                jy3              = pos[j3+7];      
                jz3              = pos[j3+8];      
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
                rinv11           = 1.0/sqrt(rsq11);
                rinv12           = 1.0/sqrt(rsq12);
                rinv13           = 1.0/sqrt(rsq13);
                rinv21           = 1.0/sqrt(rsq21);
                rinv22           = 1.0/sqrt(rsq22);
                rinv23           = 1.0/sqrt(rsq23);
                rinv31           = 1.0/sqrt(rsq31);
                rinv32           = 1.0/sqrt(rsq32);
                rinv33           = 1.0/sqrt(rsq33);
                qq               = qqOO;           
                rinvsq           = rinv11*rinv11;  
                krsq             = krf*rsq11;      
                vcoul            = qq*(rinv11+krsq-crf);
                vctot            = vctot+vcoul;    
                rinvsix          = rinvsq*rinvsq*rinvsq;
                Vvdw6            = c6*rinvsix;     
                Vvdw12           = c12*rinvsix*rinvsix;
                Vvdwtot          = Vvdwtot+Vvdw12-Vvdw6;
                fscal            = (qq*(rinv11-2.0*krsq)+12.0*Vvdw12-6.0*Vvdw6)*rinvsq;
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
                qq               = qqOH;           
                rinvsq           = rinv12*rinv12;  
                krsq             = krf*rsq12;      
                vcoul            = qq*(rinv12+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv12-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx12;     
                ty               = fscal*dy12;     
                tz               = fscal*dz12;     
                fix1             = fix1 + tx;      
                fiy1             = fiy1 + ty;      
                fiz1             = fiz1 + tz;      
                fjx2             = faction[j3+3] - tx;
                fjy2             = faction[j3+4] - ty;
                fjz2             = faction[j3+5] - tz;
                qq               = qqOH;           
                rinvsq           = rinv13*rinv13;  
                krsq             = krf*rsq13;      
                vcoul            = qq*(rinv13+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv13-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx13;     
                ty               = fscal*dy13;     
                tz               = fscal*dz13;     
                fix1             = fix1 + tx;      
                fiy1             = fiy1 + ty;      
                fiz1             = fiz1 + tz;      
                fjx3             = faction[j3+6] - tx;
                fjy3             = faction[j3+7] - ty;
                fjz3             = faction[j3+8] - tz;
                qq               = qqOH;           
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
                qq               = qqHH;           
                rinvsq           = rinv22*rinv22;  
                krsq             = krf*rsq22;      
                vcoul            = qq*(rinv22+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv22-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx22;     
                ty               = fscal*dy22;     
                tz               = fscal*dz22;     
                fix2             = fix2 + tx;      
                fiy2             = fiy2 + ty;      
                fiz2             = fiz2 + tz;      
                fjx2             = fjx2 - tx;      
                fjy2             = fjy2 - ty;      
                fjz2             = fjz2 - tz;      
                qq               = qqHH;           
                rinvsq           = rinv23*rinv23;  
                krsq             = krf*rsq23;      
                vcoul            = qq*(rinv23+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv23-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx23;     
                ty               = fscal*dy23;     
                tz               = fscal*dz23;     
                fix2             = fix2 + tx;      
                fiy2             = fiy2 + ty;      
                fiz2             = fiz2 + tz;      
                fjx3             = fjx3 - tx;      
                fjy3             = fjy3 - ty;      
                fjz3             = fjz3 - tz;      
                qq               = qqOH;           
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
                faction[j3+0]    = fjx1 - tx;      
                faction[j3+1]    = fjy1 - ty;      
                faction[j3+2]    = fjz1 - tz;      
                qq               = qqHH;           
                rinvsq           = rinv32*rinv32;  
                krsq             = krf*rsq32;      
                vcoul            = qq*(rinv32+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv32-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx32;     
                ty               = fscal*dy32;     
                tz               = fscal*dz32;     
                fix3             = fix3 + tx;      
                fiy3             = fiy3 + ty;      
                fiz3             = fiz3 + tz;      
                faction[j3+3]    = fjx2 - tx;      
                faction[j3+4]    = fjy2 - ty;      
                faction[j3+5]    = fjz2 - tz;      
                qq               = qqHH;           
                rinvsq           = rinv33*rinv33;  
                krsq             = krf*rsq33;      
                vcoul            = qq*(rinv33+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv33-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx33;     
                ty               = fscal*dy33;     
                tz               = fscal*dz33;     
                fix3             = fix3 + tx;      
                fiy3             = fiy3 + ty;      
                fiz3             = fiz3 + tz;      
                faction[j3+6]    = fjx3 - tx;      
                faction[j3+7]    = fjy3 - ty;      
                faction[j3+8]    = fjz3 - tz;      
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
            fshift[is3]      = fshift[is3]+fix1+fix2+fix3;
            fshift[is3+1]    = fshift[is3+1]+fiy1+fiy2+fiy3;
            fshift[is3+2]    = fshift[is3+2]+fiz1+fiz2+fiz3;
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
 * Gromacs nonbonded kernel nb_kernel212_adress_ex
 * Coulomb interaction:     Reaction field
 * VdW interaction:         Lennard-Jones
 * water optimization:      pairs of SPC/TIP3P interactions
 * Calculate forces:        yes
 */
void nb_kernel212_adress_ex(
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
    float         qq,vcoul,vctot;
    int           tj;
    float         rinvsix;
    float         Vvdw6,Vvdwtot;
    float         Vvdw12;
    float         krsq;
    float         ix1,iy1,iz1,fix1,fiy1,fiz1;
    float         ix2,iy2,iz2,fix2,fiy2,fiz2;
    float         ix3,iy3,iz3,fix3,fiy3,fiz3;
    float         jx1,jy1,jz1,fjx1,fjy1,fjz1;
    float         jx2,jy2,jz2,fjx2,fjy2,fjz2;
    float         jx3,jy3,jz3,fjx3,fjy3,fjz3;
    float         dx11,dy11,dz11,rsq11,rinv11;
    float         dx12,dy12,dz12,rsq12,rinv12;
    float         dx13,dy13,dz13,rsq13,rinv13;
    float         dx21,dy21,dz21,rsq21,rinv21;
    float         dx22,dy22,dz22,rsq22,rinv22;
    float         dx23,dy23,dz23,rsq23,rinv23;
    float         dx31,dy31,dz31,rsq31,rinv31;
    float         dx32,dy32,dz32,rsq32,rinv32;
    float         dx33,dy33,dz33,rsq33,rinv33;
    float         qO,qH,qqOO,qqOH,qqHH;
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
    qO               = charge[ii];     
    qH               = charge[ii+1];   
    qqOO             = facel*qO*qO;    
    qqOH             = facel*qO*qH;    
    qqHH             = facel*qH*qH;    
    tj               = 2*(ntype+1)*type[ii];
    c6               = vdwparam[tj];   
    c12              = vdwparam[tj+1]; 

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
                jx2              = pos[j3+3];      
                jy2              = pos[j3+4];      
                jz2              = pos[j3+5];      
                jx3              = pos[j3+6];      
                jy3              = pos[j3+7];      
                jz3              = pos[j3+8];      
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
                rinv11           = 1.0/sqrt(rsq11);
                rinv12           = 1.0/sqrt(rsq12);
                rinv13           = 1.0/sqrt(rsq13);
                rinv21           = 1.0/sqrt(rsq21);
                rinv22           = 1.0/sqrt(rsq22);
                rinv23           = 1.0/sqrt(rsq23);
                rinv31           = 1.0/sqrt(rsq31);
                rinv32           = 1.0/sqrt(rsq32);
                rinv33           = 1.0/sqrt(rsq33);
                qq               = qqOO;           
                rinvsq           = rinv11*rinv11;  
                krsq             = krf*rsq11;      
                vcoul            = qq*(rinv11+krsq-crf);
                vctot            = vctot+vcoul;    
                rinvsix          = rinvsq*rinvsq*rinvsq;
                Vvdw6            = c6*rinvsix;     
                Vvdw12           = c12*rinvsix*rinvsix;
                Vvdwtot          = Vvdwtot+Vvdw12-Vvdw6;
                fscal            = (qq*(rinv11-2.0*krsq)+12.0*Vvdw12-6.0*Vvdw6)*rinvsq;
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
                qq               = qqOH;           
                rinvsq           = rinv12*rinv12;  
                krsq             = krf*rsq12;      
                vcoul            = qq*(rinv12+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv12-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx12;     
                ty               = fscal*dy12;     
                tz               = fscal*dz12;     
                fix1             = fix1 + tx;      
                fiy1             = fiy1 + ty;      
                fiz1             = fiz1 + tz;      
                fjx2             = faction[j3+3] - tx;
                fjy2             = faction[j3+4] - ty;
                fjz2             = faction[j3+5] - tz;
                qq               = qqOH;           
                rinvsq           = rinv13*rinv13;  
                krsq             = krf*rsq13;      
                vcoul            = qq*(rinv13+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv13-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx13;     
                ty               = fscal*dy13;     
                tz               = fscal*dz13;     
                fix1             = fix1 + tx;      
                fiy1             = fiy1 + ty;      
                fiz1             = fiz1 + tz;      
                fjx3             = faction[j3+6] - tx;
                fjy3             = faction[j3+7] - ty;
                fjz3             = faction[j3+8] - tz;
                qq               = qqOH;           
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
                qq               = qqHH;           
                rinvsq           = rinv22*rinv22;  
                krsq             = krf*rsq22;      
                vcoul            = qq*(rinv22+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv22-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx22;     
                ty               = fscal*dy22;     
                tz               = fscal*dz22;     
                fix2             = fix2 + tx;      
                fiy2             = fiy2 + ty;      
                fiz2             = fiz2 + tz;      
                fjx2             = fjx2 - tx;      
                fjy2             = fjy2 - ty;      
                fjz2             = fjz2 - tz;      
                qq               = qqHH;           
                rinvsq           = rinv23*rinv23;  
                krsq             = krf*rsq23;      
                vcoul            = qq*(rinv23+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv23-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx23;     
                ty               = fscal*dy23;     
                tz               = fscal*dz23;     
                fix2             = fix2 + tx;      
                fiy2             = fiy2 + ty;      
                fiz2             = fiz2 + tz;      
                fjx3             = fjx3 - tx;      
                fjy3             = fjy3 - ty;      
                fjz3             = fjz3 - tz;      
                qq               = qqOH;           
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
                faction[j3+0]    = fjx1 - tx;      
                faction[j3+1]    = fjy1 - ty;      
                faction[j3+2]    = fjz1 - tz;      
                qq               = qqHH;           
                rinvsq           = rinv32*rinv32;  
                krsq             = krf*rsq32;      
                vcoul            = qq*(rinv32+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv32-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx32;     
                ty               = fscal*dy32;     
                tz               = fscal*dz32;     
                fix3             = fix3 + tx;      
                fiy3             = fiy3 + ty;      
                fiz3             = fiz3 + tz;      
                faction[j3+3]    = fjx2 - tx;      
                faction[j3+4]    = fjy2 - ty;      
                faction[j3+5]    = fjz2 - tz;      
                qq               = qqHH;           
                rinvsq           = rinv33*rinv33;  
                krsq             = krf*rsq33;      
                vcoul            = qq*(rinv33+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv33-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx33;     
                ty               = fscal*dy33;     
                tz               = fscal*dz33;     
                fix3             = fix3 + tx;      
                fiy3             = fiy3 + ty;      
                fiz3             = fiz3 + tz;      
                faction[j3+6]    = fjx3 - tx;      
                faction[j3+7]    = fjy3 - ty;      
                faction[j3+8]    = fjz3 - tz;      
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
            fshift[is3]      = fshift[is3]+fix1+fix2+fix3;
            fshift[is3+1]    = fshift[is3+1]+fiy1+fiy2+fiy3;
            fshift[is3+2]    = fshift[is3+2]+fiz1+fiz2+fiz3;
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



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
 * Gromacs nonbonded kernel nb_kernel234_adress_cg
 * Coulomb interaction:     Reaction field
 * VdW interaction:         Tabulated
 * water optimization:      pairs of TIP4P interactions
 * Calculate forces:        yes
 */
void nb_kernel234_adress_cg(
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
    float         jx1,jy1,jz1;
    float         jx2,jy2,jz2,fjx2,fjy2,fjz2;
    float         jx3,jy3,jz3,fjx3,fjy3,fjz3;
    float         jx4,jy4,jz4,fjx4,fjy4,fjz4;
    float         dx11,dy11,dz11,rsq11,rinv11;
    float         dx22,dy22,dz22,rsq22,rinv22;
    float         dx23,dy23,dz23,rsq23,rinv23;
    float         dx24,dy24,dz24,rsq24,rinv24;
    float         dx32,dy32,dz32,rsq32,rinv32;
    float         dx33,dy33,dz33,rsq33,rinv33;
    float         dx34,dy34,dz34,rsq34,rinv34;
    float         dx42,dy42,dz42,rsq42,rinv42;
    float         dx43,dy43,dz43,rsq43,rinv43;
    float         dx44,dy44,dz44,rsq44,rinv44;
    float         qH,qM,qqMM,qqMH,qqHH;
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
    qH               = charge[ii+1];   
    qM               = charge[ii+3];   
    qqMM             = facel*qM*qM;    
    qqMH             = facel*qM*qH;    
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
                jx2              = pos[j3+3];      
                jy2              = pos[j3+4];      
                jz2              = pos[j3+5];      
                jx3              = pos[j3+6];      
                jy3              = pos[j3+7];      
                jz3              = pos[j3+8];      
                jx4              = pos[j3+9];      
                jy4              = pos[j3+10];     
                jz4              = pos[j3+11];     
                dx11             = ix1 - jx1;      
                dy11             = iy1 - jy1;      
                dz11             = iz1 - jz1;      
                rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
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
                rinv11           = 1.0/sqrt(rsq11);
                rinv22           = 1.0/sqrt(rsq22);
                rinv23           = 1.0/sqrt(rsq23);
                rinv24           = 1.0/sqrt(rsq24);
                rinv32           = 1.0/sqrt(rsq32);
                rinv33           = 1.0/sqrt(rsq33);
                rinv34           = 1.0/sqrt(rsq34);
                rinv42           = 1.0/sqrt(rsq42);
                rinv43           = 1.0/sqrt(rsq43);
                rinv44           = 1.0/sqrt(rsq44);
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
                faction[j3+0]    = faction[j3+0] - tx;
                faction[j3+1]    = faction[j3+1] - ty;
                faction[j3+2]    = faction[j3+2] - tz;
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
                fjx2             = faction[j3+3] - tx;
                fjy2             = faction[j3+4] - ty;
                fjz2             = faction[j3+5] - tz;
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
                fjx3             = faction[j3+6] - tx;
                fjy3             = faction[j3+7] - ty;
                fjz3             = faction[j3+8] - tz;
                qq               = qqMH;           
                rinvsq           = rinv24*rinv24;  
                krsq             = krf*rsq24;      
                vcoul            = qq*(rinv24+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv24-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx24;     
                ty               = fscal*dy24;     
                tz               = fscal*dz24;     
                fix2             = fix2 + tx;      
                fiy2             = fiy2 + ty;      
                fiz2             = fiz2 + tz;      
                fjx4             = faction[j3+9] - tx;
                fjy4             = faction[j3+10] - ty;
                fjz4             = faction[j3+11] - tz;
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
                fjx2             = fjx2 - tx;      
                fjy2             = fjy2 - ty;      
                fjz2             = fjz2 - tz;      
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
                fjx3             = fjx3 - tx;      
                fjy3             = fjy3 - ty;      
                fjz3             = fjz3 - tz;      
                qq               = qqMH;           
                rinvsq           = rinv34*rinv34;  
                krsq             = krf*rsq34;      
                vcoul            = qq*(rinv34+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv34-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx34;     
                ty               = fscal*dy34;     
                tz               = fscal*dz34;     
                fix3             = fix3 + tx;      
                fiy3             = fiy3 + ty;      
                fiz3             = fiz3 + tz;      
                fjx4             = fjx4 - tx;      
                fjy4             = fjy4 - ty;      
                fjz4             = fjz4 - tz;      
                qq               = qqMH;           
                rinvsq           = rinv42*rinv42;  
                krsq             = krf*rsq42;      
                vcoul            = qq*(rinv42+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv42-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx42;     
                ty               = fscal*dy42;     
                tz               = fscal*dz42;     
                fix4             = fix4 + tx;      
                fiy4             = fiy4 + ty;      
                fiz4             = fiz4 + tz;      
                faction[j3+3]    = fjx2 - tx;      
                faction[j3+4]    = fjy2 - ty;      
                faction[j3+5]    = fjz2 - tz;      
                qq               = qqMH;           
                rinvsq           = rinv43*rinv43;  
                krsq             = krf*rsq43;      
                vcoul            = qq*(rinv43+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv43-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx43;     
                ty               = fscal*dy43;     
                tz               = fscal*dz43;     
                fix4             = fix4 + tx;      
                fiy4             = fiy4 + ty;      
                fiz4             = fiz4 + tz;      
                faction[j3+6]    = fjx3 - tx;      
                faction[j3+7]    = fjy3 - ty;      
                faction[j3+8]    = fjz3 - tz;      
                qq               = qqMM;           
                rinvsq           = rinv44*rinv44;  
                krsq             = krf*rsq44;      
                vcoul            = qq*(rinv44+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv44-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                tx               = fscal*dx44;     
                ty               = fscal*dy44;     
                tz               = fscal*dz44;     
                fix4             = fix4 + tx;      
                fiy4             = fiy4 + ty;      
                fiz4             = fiz4 + tz;      
                faction[j3+9]    = fjx4 - tx;      
                faction[j3+10]   = fjy4 - ty;      
                faction[j3+11]   = fjz4 - tz;      
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
 * Gromacs nonbonded kernel nb_kernel234_adress_ex
 * Coulomb interaction:     Reaction field
 * VdW interaction:         Tabulated
 * water optimization:      pairs of TIP4P interactions
 * Calculate forces:        yes
 */
void nb_kernel234_adress_ex(
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
    float         jx1,jy1,jz1;
    float         jx2,jy2,jz2,fjx2,fjy2,fjz2;
    float         jx3,jy3,jz3,fjx3,fjy3,fjz3;
    float         jx4,jy4,jz4,fjx4,fjy4,fjz4;
    float         dx11,dy11,dz11,rsq11,rinv11;
    float         dx22,dy22,dz22,rsq22,rinv22;
    float         dx23,dy23,dz23,rsq23,rinv23;
    float         dx24,dy24,dz24,rsq24,rinv24;
    float         dx32,dy32,dz32,rsq32,rinv32;
    float         dx33,dy33,dz33,rsq33,rinv33;
    float         dx34,dy34,dz34,rsq34,rinv34;
    float         dx42,dy42,dz42,rsq42,rinv42;
    float         dx43,dy43,dz43,rsq43,rinv43;
    float         dx44,dy44,dz44,rsq44,rinv44;
    float         qH,qM,qqMM,qqMH,qqHH;
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
    qH               = charge[ii+1];   
    qM               = charge[ii+3];   
    qqMM             = facel*qM*qM;    
    qqMH             = facel*qM*qH;    
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
                jx2              = pos[j3+3];      
                jy2              = pos[j3+4];      
                jz2              = pos[j3+5];      
                jx3              = pos[j3+6];      
                jy3              = pos[j3+7];      
                jz3              = pos[j3+8];      
                jx4              = pos[j3+9];      
                jy4              = pos[j3+10];     
                jz4              = pos[j3+11];     
                dx11             = ix1 - jx1;      
                dy11             = iy1 - jy1;      
                dz11             = iz1 - jz1;      
                rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
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
                rinv11           = 1.0/sqrt(rsq11);
                rinv22           = 1.0/sqrt(rsq22);
                rinv23           = 1.0/sqrt(rsq23);
                rinv24           = 1.0/sqrt(rsq24);
                rinv32           = 1.0/sqrt(rsq32);
                rinv33           = 1.0/sqrt(rsq33);
                rinv34           = 1.0/sqrt(rsq34);
                rinv42           = 1.0/sqrt(rsq42);
                rinv43           = 1.0/sqrt(rsq43);
                rinv44           = 1.0/sqrt(rsq44);
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
                faction[j3+0]    = faction[j3+0] - tx;
                faction[j3+1]    = faction[j3+1] - ty;
                faction[j3+2]    = faction[j3+2] - tz;
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
                fjx2             = faction[j3+3] - tx;
                fjy2             = faction[j3+4] - ty;
                fjz2             = faction[j3+5] - tz;
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
                fjx3             = faction[j3+6] - tx;
                fjy3             = faction[j3+7] - ty;
                fjz3             = faction[j3+8] - tz;
                qq               = qqMH;           
                rinvsq           = rinv24*rinv24;  
                krsq             = krf*rsq24;      
                vcoul            = qq*(rinv24+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv24-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx24;     
                ty               = fscal*dy24;     
                tz               = fscal*dz24;     
                fix2             = fix2 + tx;      
                fiy2             = fiy2 + ty;      
                fiz2             = fiz2 + tz;      
                fjx4             = faction[j3+9] - tx;
                fjy4             = faction[j3+10] - ty;
                fjz4             = faction[j3+11] - tz;
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
                fjx2             = fjx2 - tx;      
                fjy2             = fjy2 - ty;      
                fjz2             = fjz2 - tz;      
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
                fjx3             = fjx3 - tx;      
                fjy3             = fjy3 - ty;      
                fjz3             = fjz3 - tz;      
                qq               = qqMH;           
                rinvsq           = rinv34*rinv34;  
                krsq             = krf*rsq34;      
                vcoul            = qq*(rinv34+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv34-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx34;     
                ty               = fscal*dy34;     
                tz               = fscal*dz34;     
                fix3             = fix3 + tx;      
                fiy3             = fiy3 + ty;      
                fiz3             = fiz3 + tz;      
                fjx4             = fjx4 - tx;      
                fjy4             = fjy4 - ty;      
                fjz4             = fjz4 - tz;      
                qq               = qqMH;           
                rinvsq           = rinv42*rinv42;  
                krsq             = krf*rsq42;      
                vcoul            = qq*(rinv42+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv42-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx42;     
                ty               = fscal*dy42;     
                tz               = fscal*dz42;     
                fix4             = fix4 + tx;      
                fiy4             = fiy4 + ty;      
                fiz4             = fiz4 + tz;      
                faction[j3+3]    = fjx2 - tx;      
                faction[j3+4]    = fjy2 - ty;      
                faction[j3+5]    = fjz2 - tz;      
                qq               = qqMH;           
                rinvsq           = rinv43*rinv43;  
                krsq             = krf*rsq43;      
                vcoul            = qq*(rinv43+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv43-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx43;     
                ty               = fscal*dy43;     
                tz               = fscal*dz43;     
                fix4             = fix4 + tx;      
                fiy4             = fiy4 + ty;      
                fiz4             = fiz4 + tz;      
                faction[j3+6]    = fjx3 - tx;      
                faction[j3+7]    = fjy3 - ty;      
                faction[j3+8]    = fjz3 - tz;      
                qq               = qqMM;           
                rinvsq           = rinv44*rinv44;  
                krsq             = krf*rsq44;      
                vcoul            = qq*(rinv44+krsq-crf);
                vctot            = vctot+vcoul;    
                fscal            = (qq*(rinv44-2.0*krsq))*rinvsq;
                fscal *= hybscal;
                if(force_cap>0 && (fabs(fscal)> force_cap)){
                fscal=force_cap*fscal/fabs(fscal);
                }
                tx               = fscal*dx44;     
                ty               = fscal*dy44;     
                tz               = fscal*dz44;     
                fix4             = fix4 + tx;      
                fiy4             = fiy4 + ty;      
                fiz4             = fiz4 + tz;      
                faction[j3+9]    = fjx4 - tx;      
                faction[j3+10]   = fjy4 - ty;      
                faction[j3+11]   = fjz4 - tz;      
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



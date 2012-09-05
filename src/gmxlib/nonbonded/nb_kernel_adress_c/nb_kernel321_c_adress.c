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

#include "nb_kernel321_adress.h"



/*
 * Gromacs nonbonded kernel nb_kernel321_adress_cg
 * Coulomb interaction:     Tabulated
 * VdW interaction:         Buckingham
 * water optimization:      SPC/TIP3P - other atoms
 * Calculate forces:        yes
 */
void nb_kernel321_adress_cg(
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
    real         jq;
    real         qq,vcoul,vctot;
    int           nti;
    int           tj;
    real         rinvsix;
    real         Vvdw6,Vvdwtot;
    real         r,rt,eps,eps2;
    int           n0,nnn;
    real         Y,F,Geps,Heps2,Fp,VV;
    real         FF;
    real         fijC;
    real         Vvdwexp,br;
    real         ix1,iy1,iz1,fix1,fiy1,fiz1;
    real         ix2,iy2,iz2,fix2,fiy2,fiz2;
    real         ix3,iy3,iz3,fix3,fiy3,fiz3;
    real         jx1,jy1,jz1,fjx1,fjy1,fjz1;
    real         dx11,dy11,dz11,rsq11,rinv11;
    real         dx21,dy21,dz21,rsq21,rinv21;
    real         dx31,dy31,dz31,rsq31,rinv31;
    real         qO,qH;
    real         c6,cexp1,cexp2;
    real         weight_cg1, weight_cg2, weight_product;
    real         hybscal;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    ii               = iinr[0];        
    qO               = facel*charge[ii];
    qH               = facel*charge[ii+1];
    nti              = 3*ntype*type[ii];

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
                rinv11           = 1.0/sqrt(rsq11);
                rinv21           = 1.0/sqrt(rsq21);
                rinv31           = 1.0/sqrt(rsq31);
                jq               = charge[jnr+0];  
                qq               = qO*jq;          
                tj               = nti+3*type[jnr];
                c6               = vdwparam[tj];   
                cexp1            = vdwparam[tj+1]; 
                cexp2            = vdwparam[tj+2]; 
                rinvsq           = rinv11*rinv11;  
                r                = rsq11*rinv11;   
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 4*n0;           
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
                rinvsix          = rinvsq*rinvsq*rinvsq;
                Vvdw6            = c6*rinvsix;     
                br               = cexp2*rsq11*rinv11;
                Vvdwexp          = cexp1*exp(-br); 
                Vvdwtot          = Vvdwtot+Vvdwexp-Vvdw6;
                fscal            = (br*Vvdwexp-6.0*Vvdw6)*rinvsq-((fijC)*tabscale)*rinv11;
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
                qq               = qH*jq;          
                r                = rsq21*rinv21;   
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 4*n0;           
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
                r                = rsq31*rinv31;   
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 4*n0;           
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
 * Gromacs nonbonded kernel nb_kernel321_adress_ex
 * Coulomb interaction:     Tabulated
 * VdW interaction:         Buckingham
 * water optimization:      SPC/TIP3P - other atoms
 * Calculate forces:        yes
 */
void nb_kernel321_adress_ex(
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
    real         jq;
    real         qq,vcoul,vctot;
    int           nti;
    int           tj;
    real         rinvsix;
    real         Vvdw6,Vvdwtot;
    real         r,rt,eps,eps2;
    int           n0,nnn;
    real         Y,F,Geps,Heps2,Fp,VV;
    real         FF;
    real         fijC;
    real         Vvdwexp,br;
    real         ix1,iy1,iz1,fix1,fiy1,fiz1;
    real         ix2,iy2,iz2,fix2,fiy2,fiz2;
    real         ix3,iy3,iz3,fix3,fiy3,fiz3;
    real         jx1,jy1,jz1,fjx1,fjy1,fjz1;
    real         dx11,dy11,dz11,rsq11,rinv11;
    real         dx21,dy21,dz21,rsq21,rinv21;
    real         dx31,dy31,dz31,rsq31,rinv31;
    real         qO,qH;
    real         c6,cexp1,cexp2;
    real         weight_cg1, weight_cg2, weight_product;
    real         hybscal;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    ii               = iinr[0];        
    qO               = facel*charge[ii];
    qH               = facel*charge[ii+1];
    nti              = 3*ntype*type[ii];

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
                rinv11           = 1.0/sqrt(rsq11);
                rinv21           = 1.0/sqrt(rsq21);
                rinv31           = 1.0/sqrt(rsq31);
                jq               = charge[jnr+0];  
                qq               = qO*jq;          
                tj               = nti+3*type[jnr];
                c6               = vdwparam[tj];   
                cexp1            = vdwparam[tj+1]; 
                cexp2            = vdwparam[tj+2]; 
                rinvsq           = rinv11*rinv11;  
                r                = rsq11*rinv11;   
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 4*n0;           
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
                rinvsix          = rinvsq*rinvsq*rinvsq;
                Vvdw6            = c6*rinvsix;     
                br               = cexp2*rsq11*rinv11;
                Vvdwexp          = cexp1*exp(-br); 
                Vvdwtot          = Vvdwtot+Vvdwexp-Vvdw6;
                fscal            = (br*Vvdwexp-6.0*Vvdw6)*rinvsq-((fijC)*tabscale)*rinv11;
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
                qq               = qH*jq;          
                r                = rsq21*rinv21;   
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 4*n0;           
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
                r                = rsq31*rinv31;   
                rt               = r*tabscale;     
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 4*n0;           
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



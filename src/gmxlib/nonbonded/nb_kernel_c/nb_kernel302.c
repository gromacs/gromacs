/*
 * Copyright (c) Erik Lindahl, David van der Spoel 2003
 * 
 * This file is generated automatically at compile time
 * by the program mknb in the Gromacs distribution.
 *
 * Options used when generation this file:
 * Language:         c
 * Precision:        single
 * Threads:          no
 * Software invsqrt: no
 * PowerPC invsqrt:  no
 * Prefetch forces:  no
 * Comments:         no
 */
#ifdef HAVE_CONFIG_H
#include<config.h>
#endif
#include<math.h>
#include "types/simple.h"
#include "localpressure.h"



/*
 * Gromacs nonbonded kernel nb_kernel302
 * Coulomb interaction:     Tabulated
 * VdW interaction:         Not calculated
 * water optimization:      pairs of SPC/TIP3P interactions
 * Calculate forces:        yes
 */
void nb_kernel302(
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
      gmx_localp_grid_t * localp_grid = (gmx_localp_grid_t *)work;
    int           nri,ntype,nthreads;
    real          facel,krf,crf,tabscale,gbtabscale;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    real          shX,shY,shZ;
    real          fscal,tx,ty,tz;
    real          qq,vcoul,vctot;
    real          r,rt,eps,eps2;
    int           n0,nnn;
    real          Y,F,Geps,Heps2,Fp,VV;
    real          FF;
    real          fijC;
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
    ii               = iinr[0];        
    qO               = charge[ii];     
    qH               = charge[ii+1];   
    qqOO             = facel*qO*qO;    
    qqOH             = facel*qO*qH;    
    qqHH             = facel*qH*qH;    

    nj1              = 0;              
    
    for(n=0; (n<nri); n++)
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
        vctot            = 0;              
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
            fscal            = -((fijC)*tabscale)*rinv11;
            tx               = fscal*dx11;     
            ty               = fscal*dy11;     
            tz               = fscal*dz11;     
            gmx_spread_local_virial_on_grid(localp_grid,

            ix1,iy1,iz1,jx1,jy1,jz1,tx,ty,tz);

            fix1             = fix1 + tx;      
            fiy1             = fiy1 + ty;      
            fiz1             = fiz1 + tz;      
            fjx1             = faction[j3+0] - tx;
            fjy1             = faction[j3+1] - ty;
            fjz1             = faction[j3+2] - tz;
            qq               = qqOH;           
            r                = rsq12*rinv12;   
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
            fscal            = -((fijC)*tabscale)*rinv12;
            tx               = fscal*dx12;     
            ty               = fscal*dy12;     
            tz               = fscal*dz12;     
            gmx_spread_local_virial_on_grid(localp_grid,

            ix1,iy1,iz1,jx2,jy2,jz2,tx,ty,tz);

            fix1             = fix1 + tx;      
            fiy1             = fiy1 + ty;      
            fiz1             = fiz1 + tz;      
            fjx2             = faction[j3+3] - tx;
            fjy2             = faction[j3+4] - ty;
            fjz2             = faction[j3+5] - tz;
            qq               = qqOH;           
            r                = rsq13*rinv13;   
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
            fscal            = -((fijC)*tabscale)*rinv13;
            tx               = fscal*dx13;     
            ty               = fscal*dy13;     
            tz               = fscal*dz13;     
            gmx_spread_local_virial_on_grid(localp_grid,

            ix1,iy1,iz1,jx3,jy3,jz3,tx,ty,tz);

            fix1             = fix1 + tx;      
            fiy1             = fiy1 + ty;      
            fiz1             = fiz1 + tz;      
            fjx3             = faction[j3+6] - tx;
            fjy3             = faction[j3+7] - ty;
            fjz3             = faction[j3+8] - tz;
            qq               = qqOH;           
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
            tx               = fscal*dx21;     
            ty               = fscal*dy21;     
            tz               = fscal*dz21;     
            gmx_spread_local_virial_on_grid(localp_grid,

            ix2,iy2,iz2,jx1,jy1,jz1,tx,ty,tz);

            fix2             = fix2 + tx;      
            fiy2             = fiy2 + ty;      
            fiz2             = fiz2 + tz;      
            fjx1             = fjx1 - tx;      
            fjy1             = fjy1 - ty;      
            fjz1             = fjz1 - tz;      
            qq               = qqHH;           
            r                = rsq22*rinv22;   
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
            fscal            = -((fijC)*tabscale)*rinv22;
            tx               = fscal*dx22;     
            ty               = fscal*dy22;     
            tz               = fscal*dz22;     
            gmx_spread_local_virial_on_grid(localp_grid,

            ix2,iy2,iz2,jx2,jy2,jz2,tx,ty,tz);

            fix2             = fix2 + tx;      
            fiy2             = fiy2 + ty;      
            fiz2             = fiz2 + tz;      
            fjx2             = fjx2 - tx;      
            fjy2             = fjy2 - ty;      
            fjz2             = fjz2 - tz;      
            qq               = qqHH;           
            r                = rsq23*rinv23;   
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
            fscal            = -((fijC)*tabscale)*rinv23;
            tx               = fscal*dx23;     
            ty               = fscal*dy23;     
            tz               = fscal*dz23;     
            gmx_spread_local_virial_on_grid(localp_grid,

            ix2,iy2,iz2,jx3,jy3,jz3,tx,ty,tz);

            fix2             = fix2 + tx;      
            fiy2             = fiy2 + ty;      
            fiz2             = fiz2 + tz;      
            fjx3             = fjx3 - tx;      
            fjy3             = fjy3 - ty;      
            fjz3             = fjz3 - tz;      
            qq               = qqOH;           
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
            tx               = fscal*dx31;     
            ty               = fscal*dy31;     
            tz               = fscal*dz31;     
            gmx_spread_local_virial_on_grid(localp_grid,

            ix3,iy3,iz3,jx1,jy1,jz1,tx,ty,tz);

            fix3             = fix3 + tx;      
            fiy3             = fiy3 + ty;      
            fiz3             = fiz3 + tz;      
            faction[j3+0]    = fjx1 - tx;      
            faction[j3+1]    = fjy1 - ty;      
            faction[j3+2]    = fjz1 - tz;      
            qq               = qqHH;           
            r                = rsq32*rinv32;   
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
            fscal            = -((fijC)*tabscale)*rinv32;
            tx               = fscal*dx32;     
            ty               = fscal*dy32;     
            tz               = fscal*dz32;     
            gmx_spread_local_virial_on_grid(localp_grid,

            ix3,iy3,iz3,jx2,jy2,jz2,tx,ty,tz);

            fix3             = fix3 + tx;      
            fiy3             = fiy3 + ty;      
            fiz3             = fiz3 + tz;      
            faction[j3+3]    = fjx2 - tx;      
            faction[j3+4]    = fjy2 - ty;      
            faction[j3+5]    = fjz2 - tz;      
            qq               = qqHH;           
            r                = rsq33*rinv33;   
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
            fscal            = -((fijC)*tabscale)*rinv33;
            tx               = fscal*dx33;     
            ty               = fscal*dy33;     
            tz               = fscal*dz33;     
            gmx_spread_local_virial_on_grid(localp_grid,

            ix3,iy3,iz3,jx3,jy3,jz3,tx,ty,tz);

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
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}





/*
 * Gromacs nonbonded kernel nb_kernel302nf
 * Coulomb interaction:     Tabulated
 * VdW interaction:         Not calculated
 * water optimization:      pairs of SPC/TIP3P interactions
 * Calculate forces:        no
 */
void nb_kernel302nf(
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
      gmx_localp_grid_t * localp_grid = (gmx_localp_grid_t *)work;
    int           nri,ntype,nthreads;
    real          facel,krf,crf,tabscale,gbtabscale;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    real          shX,shY,shZ;
    real          qq,vcoul,vctot;
    real          r,rt,eps,eps2;
    int           n0,nnn;
    real          Y,F,Geps,Heps2,Fp,VV;
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
    ii               = iinr[0];        
    qO               = charge[ii];     
    qH               = charge[ii+1];   
    qqOO             = facel*qO*qO;    
    qqOH             = facel*qO*qH;    
    qqHH             = facel*qH*qH;    

    nj1              = 0;              
    
    for(n=0; (n<nri); n++)
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
        vctot            = 0;              
        
        for(k=nj0; (k<nj1); k++)
        {
            jnr              = jjnr[k];        
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
            vcoul            = qq*VV;          
            vctot            = vctot + vcoul;  
            qq               = qqOH;           
            r                = rsq12*rinv12;   
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
            vcoul            = qq*VV;          
            vctot            = vctot + vcoul;  
            qq               = qqOH;           
            r                = rsq13*rinv13;   
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
            vcoul            = qq*VV;          
            vctot            = vctot + vcoul;  
            qq               = qqOH;           
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
            vcoul            = qq*VV;          
            vctot            = vctot + vcoul;  
            qq               = qqHH;           
            r                = rsq22*rinv22;   
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
            vcoul            = qq*VV;          
            vctot            = vctot + vcoul;  
            qq               = qqHH;           
            r                = rsq23*rinv23;   
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
            vcoul            = qq*VV;          
            vctot            = vctot + vcoul;  
            qq               = qqOH;           
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
            vcoul            = qq*VV;          
            vctot            = vctot + vcoul;  
            qq               = qqHH;           
            r                = rsq32*rinv32;   
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
            vcoul            = qq*VV;          
            vctot            = vctot + vcoul;  
            qq               = qqHH;           
            r                = rsq33*rinv33;   
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
            vcoul            = qq*VV;          
            vctot            = vctot + vcoul;  
        }
        
        ggid             = gid[n];         
        Vc[ggid]         = Vc[ggid] + vctot;
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}



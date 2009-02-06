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
 * Software invsqrt: yes
 * Prefetch forces:  no
 * Comments:         no
 */
#ifdef HAVE_CONFIG_H
#include<config.h>
#endif
#include<math.h>
#include<vec.h>



/*
 * Gromacs nonbonded kernel nb_kernel430
 * Coulomb interaction:     Generalized-Born
 * VdW interaction:         Tabulated
 * water optimization:      No
 * Calculate forces:        yes
 */
void nb_kernel430_x86_64_sse(
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
                    float *         work)
{
    int           nri,ntype,nthreads;
    float         facel,krf,crf,tabscale,gbtabscale,vgb,fgb;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    float         shX,shY,shZ;
    float         fscal,tx,ty,tz;
    float         iq;
    float         qq,vcoul,vctot;
    int           nti;
    int           tj;
    float         Vvdw6,Vvdwtot;
    float         Vvdw12;
    float         r,rt,eps,eps2;
    int           n0,nnn;
    float         Y,F,Geps,Heps2,Fp,VV;
    float         FF;
    float         fijC;
    float         fijD,fijR;
    float         isai,isaj,isaprod,gbscale;
    float         dvdasum,dvdatmp,dvdaj;
    float         ix1,iy1,iz1,fix1,fiy1,fiz1;
    float         jx1,jy1,jz1;
    float         dx11,dy11,dz11,rsq11,rinv11;
    float         c6,c12;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    gbtabscale       = *p_gbtabscale;  
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
        iq               = facel*charge[ii];
        isai             = invsqrta[ii];   
        nti              = 2*ntype*type[ii];
        vctot            = 0;              
        Vvdwtot          = 0;              
        dvdasum          = 0;              
        fix1             = 0;              
        fiy1             = 0;              
        fiz1             = 0;              
        
        for(k=nj0; (k<nj1); k++)
        {
            jnr              = jjnr[k];        
            j3               = 3*jnr;          
            jx1              = pos[j3+0];      
            jy1              = pos[j3+1];      
            jz1              = pos[j3+2];      
            dx11             = ix1 - jx1;      
            dy11             = iy1 - jy1;      
            dz11             = iz1 - jz1;      
            rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
			rinv11           = invsqrt(rsq11);
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
            dvdatmp          = -(vgb+fijC*r)*0.5;
            dvdasum          = dvdasum + dvdatmp;
            dvda[jnr]        = dvdaj+dvdatmp*isaj*isaj;
            vctot            = vctot + vcoul;  
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
            fscal            = -((fijD+fijR)*tabscale+fijC-fscal)*rinv11;
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
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}





/*
 * Gromacs nonbonded kernel nb_kernel430nf
 * Coulomb interaction:     Generalized-Born
 * VdW interaction:         Tabulated
 * water optimization:      No
 * Calculate forces:        no
 */
void nb_kernel430nf_x86_64_sse(
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
                    float *         work)
{
    int           nri,ntype,nthreads;
    float         facel,krf,crf,tabscale,gbtabscale,vgb,fgb;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    float         shX,shY,shZ;
    float         iq;
    float         qq,vcoul,vctot;
    int           nti;
    int           tj;
    float         Vvdw6,Vvdwtot;
    float         Vvdw12;
    float         r,rt,eps,eps2;
    int           n0,nnn;
    float         Y,F,Geps,Heps2,Fp,VV;
    float         isai,isaj,isaprod,gbscale;
    float         ix1,iy1,iz1;
    float         jx1,jy1,jz1;
    float         dx11,dy11,dz11,rsq11,rinv11;
    float         c6,c12;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    gbtabscale       = *p_gbtabscale;  
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
        iq               = facel*charge[ii];
        isai             = invsqrta[ii];   
        nti              = 2*ntype*type[ii];
        vctot            = 0;              
        Vvdwtot          = 0;              
        
        for(k=nj0; (k<nj1); k++)
        {
            jnr              = jjnr[k];        
            j3               = 3*jnr;          
            jx1              = pos[j3+0];      
            jy1              = pos[j3+1];      
            jz1              = pos[j3+2];      
            dx11             = ix1 - jx1;      
            dy11             = iy1 - jy1;      
            dz11             = iz1 - jz1;      
            rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
            rinv11           = invsqrt(rsq11);
            isaj             = invsqrta[jnr];  
            isaprod          = isai*isaj;      
            qq               = iq*charge[jnr]; 
            vcoul            = qq*rinv11;      
            qq               = isaprod*(-qq);  
            gbscale          = isaprod*gbtabscale;
            tj               = nti+2*type[jnr];
            c6               = vdwparam[tj];   
            c12              = vdwparam[tj+1]; 
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
            vgb              = qq*VV;          
            vctot            = vctot + vcoul;  
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
            Vvdw6            = c6*VV;          
            nnn              = nnn+4;          
            Y                = VFtab[nnn];     
            F                = VFtab[nnn+1];   
            Geps             = eps*VFtab[nnn+2];
            Heps2            = eps2*VFtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            Vvdw12           = c12*VV;         
            Vvdwtot          = Vvdwtot+ Vvdw6 + Vvdw12;
        }
        
        ggid             = gid[n];         
        Vc[ggid]         = Vc[ggid] + vctot;
        Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}



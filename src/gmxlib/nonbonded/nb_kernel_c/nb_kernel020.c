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
 * Gromacs nonbonded kernel nb_kernel020
 * Coulomb interaction:     Not calculated
 * VdW interaction:         Buckingham
 * water optimization:      No
 * Calculate forces:        yes
 */
void nb_kernel020(
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
    real          rinvsq;
    int           nti;
    int           tj;
    real          rinvsix;
    real          Vvdw6,Vvdwtot;
    real          Vvdwexp,br;
    real          ix1,iy1,iz1,fix1,fiy1,fiz1;
    real          jx1,jy1,jz1;
    real          dx11,dy11,dz11,rsq11,rinv11;
    real          c6,cexp1,cexp2;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
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
        nti              = 3*ntype*type[ii];
        Vvdwtot          = 0;              
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
            rinv11           = 1.0/sqrt(rsq11);
            tj               = nti+3*type[jnr];
            c6               = vdwparam[tj];   
            cexp1            = vdwparam[tj+1]; 
            cexp2            = vdwparam[tj+2]; 
            rinvsq           = rinv11*rinv11;  
            rinvsix          = rinvsq*rinvsq*rinvsq;
            Vvdw6            = c6*rinvsix;     
            br               = cexp2*rsq11*rinv11;
            Vvdwexp          = cexp1*exp(-br); 
            Vvdwtot          = Vvdwtot+Vvdwexp-Vvdw6;
            fscal            = (br*Vvdwexp-6.0*Vvdw6)*rinvsq;
            tx               = fscal*dx11;     
            ty               = fscal*dy11;     
            tz               = fscal*dz11;     
            gmx_spread_local_virial_on_grid(localp_grid,

            ix1,iy1,iz1,jx1,jy1,jz1,tx,ty,tz);

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
        Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}





/*
 * Gromacs nonbonded kernel nb_kernel020nf
 * Coulomb interaction:     Not calculated
 * VdW interaction:         Buckingham
 * water optimization:      No
 * Calculate forces:        no
 */
void nb_kernel020nf(
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
    real          rinvsq;
    int           nti;
    int           tj;
    real          rinvsix;
    real          Vvdw6,Vvdwtot;
    real          Vvdwexp,br;
    real          ix1,iy1,iz1;
    real          jx1,jy1,jz1;
    real          dx11,dy11,dz11,rsq11,rinv11;
    real          c6,cexp1,cexp2;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
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
        nti              = 3*ntype*type[ii];
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
            rinv11           = 1.0/sqrt(rsq11);
            tj               = nti+3*type[jnr];
            c6               = vdwparam[tj];   
            cexp1            = vdwparam[tj+1]; 
            cexp2            = vdwparam[tj+2]; 
            rinvsq           = rinv11*rinv11;  
            rinvsix          = rinvsq*rinvsq*rinvsq;
            Vvdw6            = c6*rinvsix;     
            br               = cexp2*rsq11*rinv11;
            Vvdwexp          = cexp1*exp(-br); 
            Vvdwtot          = Vvdwtot+Vvdwexp-Vvdw6;
        }
        
        ggid             = gid[n];         
        Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}



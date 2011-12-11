C
C                This source code is part of
C
C                 G   R   O   M   A   C   S
C
C Copyright (c) 1991-2000, University of Groningen, The Netherlands.
C Copyright (c) 2001-2009, The GROMACS Development Team
C
C Gromacs is a library for molecular simulation and trajectory analysis,
C written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
C a full list of developers and information, check out http://www.gromacs.org
C
C This program is free software; you can redistribute it and/or modify it under
C the terms of the GNU Lesser General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version.
C As a special exception, you may use this file as part of a free software
C library without restriction.  Specifically, if other files instantiate
C templates or use macros or inline functions from this file, or you compile
C this file and link it with other files to produce an executable, this
C file does not by itself cause the resulting executable to be covered by
C the GNU Lesser General Public License.
C
C In plain-speak: do not worry about classes/macros/templates either - only
C changes to the library have to be LGPL, not an application linking with it.
C
C To help fund GROMACS development, we humbly ask that you cite
C the papers people have written on it - you can find them on the website!
C

#ifdef HAVE_CONFIG_H
#  include<config.h>
#endif

#ifdef GMX_DOUBLE
#  define gmxreal real*8
#else
#  define gmxreal real*4
#endif



C
C Gromacs nonbonded kernel pwr6kernel334
C Coulomb interaction:     Tabulated
C VdW interaction:         Tabulated
C water optimization:      pairs of TIP4P interactions
C Calculate forces:        yes
C
      subroutine pwr6kernel334(
     &                          nri,
     &                          iinr,
     &                          jindex,
     &                          jjnr,
     &                          shift,
     &                          shiftvec,
     &                          fshift,
     &                          gid,
     &                          pos,
     &                          faction,
     &                          charge,
     &                          facel,
     &                          krf,
     &                          crf,
     &                          Vc,
     &                          type,
     &                          ntype,
     &                          vdwparam,
     &                          Vvdw,
     &                          tabscale,
     &                          VFtab,
     &                          invsqrta,
     &                          dvda,
     &                          gbtabscale,
     &                          GBtab,
     &                          nthreads,
     &                          count,
     &                          mtx,
     &                          outeriter,
     &                          inneriter,
     &                          work)
      implicit      none
      integer*4     nri,iinr(*),jindex(*),jjnr(*),shift(*)
      gmxreal       shiftvec(*),fshift(*),pos(*),faction(*)
      integer*4     gid(*),type(*),ntype
      gmxreal       charge(*),facel,krf,crf,Vc(*),vdwparam(*)
      gmxreal       Vvdw(*),tabscale,VFtab(*)
      gmxreal       invsqrta(*),dvda(*),gbtabscale,GBtab(*)
      integer*4     nthreads,count,mtx,outeriter,inneriter
      gmxreal       work(*)

      integer*4     n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid
      integer*4     nn0,nn1,nouter,ninner
      gmxreal       shX,shY,shZ
      gmxreal       fscal,tx,ty,tz
      gmxreal       qq,vcoul,vctot
      integer*4     tj
      gmxreal       Vvdw6,Vvdwtot
      gmxreal       Vvdw12
      gmxreal       r,rt,eps,eps2
      integer*4     n0,nnn
      gmxreal       Y,F,Geps,Heps2,Fp,VV
      gmxreal       FF
      gmxreal       fijC
      gmxreal       fijD,fijR
      gmxreal       ix1,iy1,iz1,fix1,fiy1,fiz1
      gmxreal       ix2,iy2,iz2,fix2,fiy2,fiz2
      gmxreal       ix3,iy3,iz3,fix3,fiy3,fiz3
      gmxreal       ix4,iy4,iz4,fix4,fiy4,fiz4
      gmxreal       jx1,jy1,jz1
      gmxreal       jx2,jy2,jz2,fjx2,fjy2,fjz2
      gmxreal       jx3,jy3,jz3,fjx3,fjy3,fjz3
      gmxreal       jx4,jy4,jz4,fjx4,fjy4,fjz4
      gmxreal       dx11,dy11,dz11,rsq11,rinv11
      gmxreal       dx22,dy22,dz22,rsq22,rinv22
      gmxreal       dx23,dy23,dz23,rsq23,rinv23
      gmxreal       dx24,dy24,dz24,rsq24,rinv24
      gmxreal       dx32,dy32,dz32,rsq32,rinv32
      gmxreal       dx33,dy33,dz33,rsq33,rinv33
      gmxreal       dx34,dy34,dz34,rsq34,rinv34
      gmxreal       dx42,dy42,dz42,rsq42,rinv42
      gmxreal       dx43,dy43,dz43,rsq43,rinv43
      gmxreal       dx44,dy44,dz44,rsq44,rinv44
      gmxreal       qH,qM,qqMM,qqMH,qqHH
      gmxreal       c6,c12


C    Initialize water data
      ii               = iinr(1)+1       
      qH               = charge(ii+1)    
      qM               = charge(ii+3)    
      qqMM             = facel*qM*qM     
      qqMH             = facel*qM*qH     
      qqHH             = facel*qH*qH     
      tj               = 2*(ntype+1)*type(ii)+1
      c6               = vdwparam(tj)    
      c12              = vdwparam(tj+1)  


C    Reset outer and inner iteration counters
      nouter           = 0               
      ninner           = 0               

C    Loop over thread workunits
   10 call pwr6kernelsync(mtx,count,nri,nthreads,nn0,nn1)
        if(nn1.gt.nri) nn1=nri

C      Start outer loop over neighborlists
        
        do n=nn0+1,nn1

C        Load shift vector for this list
          is3              = 3*shift(n)+1    
          shX              = shiftvec(is3)   
          shY              = shiftvec(is3+1) 
          shZ              = shiftvec(is3+2) 

C        Load limits for loop over neighbors
          nj0              = jindex(n)+1     
          nj1              = jindex(n+1)     

C        Get outer coordinate index
          ii               = iinr(n)+1       
          ii3              = 3*ii-2          

C        Load i atom data, add shift vector
          ix1              = shX + pos(ii3+0)
          iy1              = shY + pos(ii3+1)
          iz1              = shZ + pos(ii3+2)
          ix2              = shX + pos(ii3+3)
          iy2              = shY + pos(ii3+4)
          iz2              = shZ + pos(ii3+5)
          ix3              = shX + pos(ii3+6)
          iy3              = shY + pos(ii3+7)
          iz3              = shZ + pos(ii3+8)
          ix4              = shX + pos(ii3+9)
          iy4              = shY + pos(ii3+10)
          iz4              = shZ + pos(ii3+11)

C        Zero the potential energy for this list
          vctot            = 0               
          Vvdwtot          = 0               

C        Clear i atom forces
          fix1             = 0               
          fiy1             = 0               
          fiz1             = 0               
          fix2             = 0               
          fiy2             = 0               
          fiz2             = 0               
          fix3             = 0               
          fiy3             = 0               
          fiz3             = 0               
          fix4             = 0               
          fiy4             = 0               
          fiz4             = 0               
          
          do k=nj0,nj1

C          Get j neighbor index, and coordinate index
            jnr              = jjnr(k)+1       
            j3               = 3*jnr-2         

C          load j atom coordinates
            jx1              = pos(j3+0)       
            jy1              = pos(j3+1)       
            jz1              = pos(j3+2)       
            jx2              = pos(j3+3)       
            jy2              = pos(j3+4)       
            jz2              = pos(j3+5)       
            jx3              = pos(j3+6)       
            jy3              = pos(j3+7)       
            jz3              = pos(j3+8)       
            jx4              = pos(j3+9)       
            jy4              = pos(j3+10)      
            jz4              = pos(j3+11)      

C          Calculate distance
            dx11             = ix1 - jx1       
            dy11             = iy1 - jy1       
            dz11             = iz1 - jz1       
            rsq11            = dx11*dx11+dy11*dy11+dz11*dz11
            dx22             = ix2 - jx2       
            dy22             = iy2 - jy2       
            dz22             = iz2 - jz2       
            rsq22            = dx22*dx22+dy22*dy22+dz22*dz22
            dx23             = ix2 - jx3       
            dy23             = iy2 - jy3       
            dz23             = iz2 - jz3       
            rsq23            = dx23*dx23+dy23*dy23+dz23*dz23
            dx24             = ix2 - jx4       
            dy24             = iy2 - jy4       
            dz24             = iz2 - jz4       
            rsq24            = dx24*dx24+dy24*dy24+dz24*dz24
            dx32             = ix3 - jx2       
            dy32             = iy3 - jy2       
            dz32             = iz3 - jz2       
            rsq32            = dx32*dx32+dy32*dy32+dz32*dz32
            dx33             = ix3 - jx3       
            dy33             = iy3 - jy3       
            dz33             = iz3 - jz3       
            rsq33            = dx33*dx33+dy33*dy33+dz33*dz33
            dx34             = ix3 - jx4       
            dy34             = iy3 - jy4       
            dz34             = iz3 - jz4       
            rsq34            = dx34*dx34+dy34*dy34+dz34*dz34
            dx42             = ix4 - jx2       
            dy42             = iy4 - jy2       
            dz42             = iz4 - jz2       
            rsq42            = dx42*dx42+dy42*dy42+dz42*dz42
            dx43             = ix4 - jx3       
            dy43             = iy4 - jy3       
            dz43             = iz4 - jz3       
            rsq43            = dx43*dx43+dy43*dy43+dz43*dz43
            dx44             = ix4 - jx4       
            dy44             = iy4 - jy4       
            dz44             = iz4 - jz4       
            rsq44            = dx44*dx44+dy44*dy44+dz44*dz44

C          Calculate 1/r and 1/r2

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv11           = frsqrtes(rsq11) 
#else
            rinv11           = frsqrte(dble(rsq11)) 
#endif
            rinv11           = (0.5*rinv11*(3.0-((rsq11*rinv11)
     &  *rinv11)))
#ifdef GMX_DOUBLE
            rinv11           = (0.5*rinv11*(3.0-((rsq11*rinv11)
     &  *rinv11)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv22           = frsqrtes(rsq22) 
#else
            rinv22           = frsqrte(dble(rsq22)) 
#endif
            rinv22           = (0.5*rinv22*(3.0-((rsq22*rinv22)
     &  *rinv22)))
#ifdef GMX_DOUBLE
            rinv22           = (0.5*rinv22*(3.0-((rsq22*rinv22)
     &  *rinv22)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv23           = frsqrtes(rsq23) 
#else
            rinv23           = frsqrte(dble(rsq23)) 
#endif
            rinv23           = (0.5*rinv23*(3.0-((rsq23*rinv23)
     &  *rinv23)))
#ifdef GMX_DOUBLE
            rinv23           = (0.5*rinv23*(3.0-((rsq23*rinv23)
     &  *rinv23)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv24           = frsqrtes(rsq24) 
#else
            rinv24           = frsqrte(dble(rsq24)) 
#endif
            rinv24           = (0.5*rinv24*(3.0-((rsq24*rinv24)
     &  *rinv24)))
#ifdef GMX_DOUBLE
            rinv24           = (0.5*rinv24*(3.0-((rsq24*rinv24)
     &  *rinv24)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv32           = frsqrtes(rsq32) 
#else
            rinv32           = frsqrte(dble(rsq32)) 
#endif
            rinv32           = (0.5*rinv32*(3.0-((rsq32*rinv32)
     &  *rinv32)))
#ifdef GMX_DOUBLE
            rinv32           = (0.5*rinv32*(3.0-((rsq32*rinv32)
     &  *rinv32)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv33           = frsqrtes(rsq33) 
#else
            rinv33           = frsqrte(dble(rsq33)) 
#endif
            rinv33           = (0.5*rinv33*(3.0-((rsq33*rinv33)
     &  *rinv33)))
#ifdef GMX_DOUBLE
            rinv33           = (0.5*rinv33*(3.0-((rsq33*rinv33)
     &  *rinv33)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv34           = frsqrtes(rsq34) 
#else
            rinv34           = frsqrte(dble(rsq34)) 
#endif
            rinv34           = (0.5*rinv34*(3.0-((rsq34*rinv34)
     &  *rinv34)))
#ifdef GMX_DOUBLE
            rinv34           = (0.5*rinv34*(3.0-((rsq34*rinv34)
     &  *rinv34)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv42           = frsqrtes(rsq42) 
#else
            rinv42           = frsqrte(dble(rsq42)) 
#endif
            rinv42           = (0.5*rinv42*(3.0-((rsq42*rinv42)
     &  *rinv42)))
#ifdef GMX_DOUBLE
            rinv42           = (0.5*rinv42*(3.0-((rsq42*rinv42)
     &  *rinv42)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv43           = frsqrtes(rsq43) 
#else
            rinv43           = frsqrte(dble(rsq43)) 
#endif
            rinv43           = (0.5*rinv43*(3.0-((rsq43*rinv43)
     &  *rinv43)))
#ifdef GMX_DOUBLE
            rinv43           = (0.5*rinv43*(3.0-((rsq43*rinv43)
     &  *rinv43)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv44           = frsqrtes(rsq44) 
#else
            rinv44           = frsqrte(dble(rsq44)) 
#endif
            rinv44           = (0.5*rinv44*(3.0-((rsq44*rinv44)
     &  *rinv44)))
#ifdef GMX_DOUBLE
            rinv44           = (0.5*rinv44*(3.0-((rsq44*rinv44)
     &  *rinv44)))
#endif

C          Load parameters for j atom

C          Calculate table index
            r                = rsq11*rinv11    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated VdW interaction - dispersion
            nnn              = nnn+4           
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            FF               = Fp+Geps+2.0*Heps2
            Vvdw6            = c6*VV           
            fijD             = c6*FF           

C          Tabulated VdW interaction - repulsion
            nnn              = nnn+4           
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            FF               = Fp+Geps+2.0*Heps2
            Vvdw12           = c12*VV          
            fijR             = c12*FF          
            Vvdwtot          = Vvdwtot+ Vvdw6 + Vvdw12
            fscal            = -((fijD+fijR)*tabscale)*rinv11

C          Calculate temporary vectorial force
            tx               = fscal*dx11      
            ty               = fscal*dy11      
            tz               = fscal*dz11      

C          Increment i atom force
            fix1             = fix1 + tx       
            fiy1             = fiy1 + ty       
            fiz1             = fiz1 + tz       

C          Decrement j atom force
            faction(j3+0)    = faction(j3+0) - tx
            faction(j3+1)    = faction(j3+1) - ty
            faction(j3+2)    = faction(j3+2) - tz

C          Load parameters for j atom
            qq               = qqHH            

C          Calculate table index
            r                = rsq22*rinv22    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            FF               = Fp+Geps+2.0*Heps2
            vcoul            = qq*VV           
            fijC             = qq*FF           
            vctot            = vctot + vcoul   
            fscal            = -((fijC)*tabscale)*rinv22

C          Calculate temporary vectorial force
            tx               = fscal*dx22      
            ty               = fscal*dy22      
            tz               = fscal*dz22      

C          Increment i atom force
            fix2             = fix2 + tx       
            fiy2             = fiy2 + ty       
            fiz2             = fiz2 + tz       

C          Decrement j atom force
            fjx2             = faction(j3+3) - tx
            fjy2             = faction(j3+4) - ty
            fjz2             = faction(j3+5) - tz

C          Load parameters for j atom
            qq               = qqHH            

C          Calculate table index
            r                = rsq23*rinv23    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            FF               = Fp+Geps+2.0*Heps2
            vcoul            = qq*VV           
            fijC             = qq*FF           
            vctot            = vctot + vcoul   
            fscal            = -((fijC)*tabscale)*rinv23

C          Calculate temporary vectorial force
            tx               = fscal*dx23      
            ty               = fscal*dy23      
            tz               = fscal*dz23      

C          Increment i atom force
            fix2             = fix2 + tx       
            fiy2             = fiy2 + ty       
            fiz2             = fiz2 + tz       

C          Decrement j atom force
            fjx3             = faction(j3+6) - tx
            fjy3             = faction(j3+7) - ty
            fjz3             = faction(j3+8) - tz

C          Load parameters for j atom
            qq               = qqMH            

C          Calculate table index
            r                = rsq24*rinv24    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            FF               = Fp+Geps+2.0*Heps2
            vcoul            = qq*VV           
            fijC             = qq*FF           
            vctot            = vctot + vcoul   
            fscal            = -((fijC)*tabscale)*rinv24

C          Calculate temporary vectorial force
            tx               = fscal*dx24      
            ty               = fscal*dy24      
            tz               = fscal*dz24      

C          Increment i atom force
            fix2             = fix2 + tx       
            fiy2             = fiy2 + ty       
            fiz2             = fiz2 + tz       

C          Decrement j atom force
            fjx4             = faction(j3+9) - tx
            fjy4             = faction(j3+10) - ty
            fjz4             = faction(j3+11) - tz

C          Load parameters for j atom
            qq               = qqHH            

C          Calculate table index
            r                = rsq32*rinv32    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            FF               = Fp+Geps+2.0*Heps2
            vcoul            = qq*VV           
            fijC             = qq*FF           
            vctot            = vctot + vcoul   
            fscal            = -((fijC)*tabscale)*rinv32

C          Calculate temporary vectorial force
            tx               = fscal*dx32      
            ty               = fscal*dy32      
            tz               = fscal*dz32      

C          Increment i atom force
            fix3             = fix3 + tx       
            fiy3             = fiy3 + ty       
            fiz3             = fiz3 + tz       

C          Decrement j atom force
            fjx2             = fjx2 - tx       
            fjy2             = fjy2 - ty       
            fjz2             = fjz2 - tz       

C          Load parameters for j atom
            qq               = qqHH            

C          Calculate table index
            r                = rsq33*rinv33    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            FF               = Fp+Geps+2.0*Heps2
            vcoul            = qq*VV           
            fijC             = qq*FF           
            vctot            = vctot + vcoul   
            fscal            = -((fijC)*tabscale)*rinv33

C          Calculate temporary vectorial force
            tx               = fscal*dx33      
            ty               = fscal*dy33      
            tz               = fscal*dz33      

C          Increment i atom force
            fix3             = fix3 + tx       
            fiy3             = fiy3 + ty       
            fiz3             = fiz3 + tz       

C          Decrement j atom force
            fjx3             = fjx3 - tx       
            fjy3             = fjy3 - ty       
            fjz3             = fjz3 - tz       

C          Load parameters for j atom
            qq               = qqMH            

C          Calculate table index
            r                = rsq34*rinv34    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            FF               = Fp+Geps+2.0*Heps2
            vcoul            = qq*VV           
            fijC             = qq*FF           
            vctot            = vctot + vcoul   
            fscal            = -((fijC)*tabscale)*rinv34

C          Calculate temporary vectorial force
            tx               = fscal*dx34      
            ty               = fscal*dy34      
            tz               = fscal*dz34      

C          Increment i atom force
            fix3             = fix3 + tx       
            fiy3             = fiy3 + ty       
            fiz3             = fiz3 + tz       

C          Decrement j atom force
            fjx4             = fjx4 - tx       
            fjy4             = fjy4 - ty       
            fjz4             = fjz4 - tz       

C          Load parameters for j atom
            qq               = qqMH            

C          Calculate table index
            r                = rsq42*rinv42    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            FF               = Fp+Geps+2.0*Heps2
            vcoul            = qq*VV           
            fijC             = qq*FF           
            vctot            = vctot + vcoul   
            fscal            = -((fijC)*tabscale)*rinv42

C          Calculate temporary vectorial force
            tx               = fscal*dx42      
            ty               = fscal*dy42      
            tz               = fscal*dz42      

C          Increment i atom force
            fix4             = fix4 + tx       
            fiy4             = fiy4 + ty       
            fiz4             = fiz4 + tz       

C          Decrement j atom force
            faction(j3+3)    = fjx2 - tx       
            faction(j3+4)    = fjy2 - ty       
            faction(j3+5)    = fjz2 - tz       

C          Load parameters for j atom
            qq               = qqMH            

C          Calculate table index
            r                = rsq43*rinv43    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            FF               = Fp+Geps+2.0*Heps2
            vcoul            = qq*VV           
            fijC             = qq*FF           
            vctot            = vctot + vcoul   
            fscal            = -((fijC)*tabscale)*rinv43

C          Calculate temporary vectorial force
            tx               = fscal*dx43      
            ty               = fscal*dy43      
            tz               = fscal*dz43      

C          Increment i atom force
            fix4             = fix4 + tx       
            fiy4             = fiy4 + ty       
            fiz4             = fiz4 + tz       

C          Decrement j atom force
            faction(j3+6)    = fjx3 - tx       
            faction(j3+7)    = fjy3 - ty       
            faction(j3+8)    = fjz3 - tz       

C          Load parameters for j atom
            qq               = qqMM            

C          Calculate table index
            r                = rsq44*rinv44    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            FF               = Fp+Geps+2.0*Heps2
            vcoul            = qq*VV           
            fijC             = qq*FF           
            vctot            = vctot + vcoul   
            fscal            = -((fijC)*tabscale)*rinv44

C          Calculate temporary vectorial force
            tx               = fscal*dx44      
            ty               = fscal*dy44      
            tz               = fscal*dz44      

C          Increment i atom force
            fix4             = fix4 + tx       
            fiy4             = fiy4 + ty       
            fiz4             = fiz4 + tz       

C          Decrement j atom force
            faction(j3+9)    = fjx4 - tx       
            faction(j3+10)   = fjy4 - ty       
            faction(j3+11)   = fjz4 - tz       

C          Inner loop uses 433 flops/iteration
          end do
          

C        Add i forces to mem and shifted force list
          faction(ii3+0)   = faction(ii3+0) + fix1
          faction(ii3+1)   = faction(ii3+1) + fiy1
          faction(ii3+2)   = faction(ii3+2) + fiz1
          faction(ii3+3)   = faction(ii3+3) + fix2
          faction(ii3+4)   = faction(ii3+4) + fiy2
          faction(ii3+5)   = faction(ii3+5) + fiz2
          faction(ii3+6)   = faction(ii3+6) + fix3
          faction(ii3+7)   = faction(ii3+7) + fiy3
          faction(ii3+8)   = faction(ii3+8) + fiz3
          faction(ii3+9)   = faction(ii3+9) + fix4
          faction(ii3+10)  = faction(ii3+10) + fiy4
          faction(ii3+11)  = faction(ii3+11) + fiz4
          fshift(is3)      = fshift(is3)+fix1+fix2+fix3+fix4
          fshift(is3+1)    = fshift(is3+1)+fiy1+fiy2+fiy3+fiy4
          fshift(is3+2)    = fshift(is3+2)+fiz1+fiz2+fiz3+fiz4

C        Add potential energies to the group for this list
          ggid             = gid(n)+1        
          Vc(ggid)         = Vc(ggid) + vctot
          Vvdw(ggid)       = Vvdw(ggid) + Vvdwtot

C        Increment number of inner iterations
          ninner           = ninner + nj1 - nj0

C        Outer loop uses 38 flops/iteration
        end do
        

C      Increment number of outer iterations
        nouter           = nouter + nn1 - nn0
      if(nn1.lt.nri) goto 10

C    Write outer/inner iteration count to pointers
      outeriter        = nouter          
      inneriter        = ninner          
      return
      end






C
C Gromacs nonbonded kernel pwr6kernel334nf
C Coulomb interaction:     Tabulated
C VdW interaction:         Tabulated
C water optimization:      pairs of TIP4P interactions
C Calculate forces:        no
C
      subroutine pwr6kernel334nf(
     &                          nri,
     &                          iinr,
     &                          jindex,
     &                          jjnr,
     &                          shift,
     &                          shiftvec,
     &                          fshift,
     &                          gid,
     &                          pos,
     &                          faction,
     &                          charge,
     &                          facel,
     &                          krf,
     &                          crf,
     &                          Vc,
     &                          type,
     &                          ntype,
     &                          vdwparam,
     &                          Vvdw,
     &                          tabscale,
     &                          VFtab,
     &                          invsqrta,
     &                          dvda,
     &                          gbtabscale,
     &                          GBtab,
     &                          nthreads,
     &                          count,
     &                          mtx,
     &                          outeriter,
     &                          inneriter,
     &                          work)
      implicit      none
      integer*4     nri,iinr(*),jindex(*),jjnr(*),shift(*)
      gmxreal       shiftvec(*),fshift(*),pos(*),faction(*)
      integer*4     gid(*),type(*),ntype
      gmxreal       charge(*),facel,krf,crf,Vc(*),vdwparam(*)
      gmxreal       Vvdw(*),tabscale,VFtab(*)
      gmxreal       invsqrta(*),dvda(*),gbtabscale,GBtab(*)
      integer*4     nthreads,count,mtx,outeriter,inneriter
      gmxreal       work(*)

      integer*4     n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid
      integer*4     nn0,nn1,nouter,ninner
      gmxreal       shX,shY,shZ
      gmxreal       qq,vcoul,vctot
      integer*4     tj
      gmxreal       Vvdw6,Vvdwtot
      gmxreal       Vvdw12
      gmxreal       r,rt,eps,eps2
      integer*4     n0,nnn
      gmxreal       Y,F,Geps,Heps2,Fp,VV
      gmxreal       ix1,iy1,iz1
      gmxreal       ix2,iy2,iz2
      gmxreal       ix3,iy3,iz3
      gmxreal       ix4,iy4,iz4
      gmxreal       jx1,jy1,jz1
      gmxreal       jx2,jy2,jz2
      gmxreal       jx3,jy3,jz3
      gmxreal       jx4,jy4,jz4
      gmxreal       dx11,dy11,dz11,rsq11,rinv11
      gmxreal       dx22,dy22,dz22,rsq22,rinv22
      gmxreal       dx23,dy23,dz23,rsq23,rinv23
      gmxreal       dx24,dy24,dz24,rsq24,rinv24
      gmxreal       dx32,dy32,dz32,rsq32,rinv32
      gmxreal       dx33,dy33,dz33,rsq33,rinv33
      gmxreal       dx34,dy34,dz34,rsq34,rinv34
      gmxreal       dx42,dy42,dz42,rsq42,rinv42
      gmxreal       dx43,dy43,dz43,rsq43,rinv43
      gmxreal       dx44,dy44,dz44,rsq44,rinv44
      gmxreal       qH,qM,qqMM,qqMH,qqHH
      gmxreal       c6,c12


C    Initialize water data
      ii               = iinr(1)+1       
      qH               = charge(ii+1)    
      qM               = charge(ii+3)    
      qqMM             = facel*qM*qM     
      qqMH             = facel*qM*qH     
      qqHH             = facel*qH*qH     
      tj               = 2*(ntype+1)*type(ii)+1
      c6               = vdwparam(tj)    
      c12              = vdwparam(tj+1)  


C    Reset outer and inner iteration counters
      nouter           = 0               
      ninner           = 0               

C    Loop over thread workunits
   10 call pwr6kernelsync(mtx,count,nri,nthreads,nn0,nn1)
        if(nn1.gt.nri) nn1=nri

C      Start outer loop over neighborlists
        
        do n=nn0+1,nn1

C        Load shift vector for this list
          is3              = 3*shift(n)+1    
          shX              = shiftvec(is3)   
          shY              = shiftvec(is3+1) 
          shZ              = shiftvec(is3+2) 

C        Load limits for loop over neighbors
          nj0              = jindex(n)+1     
          nj1              = jindex(n+1)     

C        Get outer coordinate index
          ii               = iinr(n)+1       
          ii3              = 3*ii-2          

C        Load i atom data, add shift vector
          ix1              = shX + pos(ii3+0)
          iy1              = shY + pos(ii3+1)
          iz1              = shZ + pos(ii3+2)
          ix2              = shX + pos(ii3+3)
          iy2              = shY + pos(ii3+4)
          iz2              = shZ + pos(ii3+5)
          ix3              = shX + pos(ii3+6)
          iy3              = shY + pos(ii3+7)
          iz3              = shZ + pos(ii3+8)
          ix4              = shX + pos(ii3+9)
          iy4              = shY + pos(ii3+10)
          iz4              = shZ + pos(ii3+11)

C        Zero the potential energy for this list
          vctot            = 0               
          Vvdwtot          = 0               

C        Clear i atom forces
          
          do k=nj0,nj1

C          Get j neighbor index, and coordinate index
            jnr              = jjnr(k)+1       
            j3               = 3*jnr-2         

C          load j atom coordinates
            jx1              = pos(j3+0)       
            jy1              = pos(j3+1)       
            jz1              = pos(j3+2)       
            jx2              = pos(j3+3)       
            jy2              = pos(j3+4)       
            jz2              = pos(j3+5)       
            jx3              = pos(j3+6)       
            jy3              = pos(j3+7)       
            jz3              = pos(j3+8)       
            jx4              = pos(j3+9)       
            jy4              = pos(j3+10)      
            jz4              = pos(j3+11)      

C          Calculate distance
            dx11             = ix1 - jx1       
            dy11             = iy1 - jy1       
            dz11             = iz1 - jz1       
            rsq11            = dx11*dx11+dy11*dy11+dz11*dz11
            dx22             = ix2 - jx2       
            dy22             = iy2 - jy2       
            dz22             = iz2 - jz2       
            rsq22            = dx22*dx22+dy22*dy22+dz22*dz22
            dx23             = ix2 - jx3       
            dy23             = iy2 - jy3       
            dz23             = iz2 - jz3       
            rsq23            = dx23*dx23+dy23*dy23+dz23*dz23
            dx24             = ix2 - jx4       
            dy24             = iy2 - jy4       
            dz24             = iz2 - jz4       
            rsq24            = dx24*dx24+dy24*dy24+dz24*dz24
            dx32             = ix3 - jx2       
            dy32             = iy3 - jy2       
            dz32             = iz3 - jz2       
            rsq32            = dx32*dx32+dy32*dy32+dz32*dz32
            dx33             = ix3 - jx3       
            dy33             = iy3 - jy3       
            dz33             = iz3 - jz3       
            rsq33            = dx33*dx33+dy33*dy33+dz33*dz33
            dx34             = ix3 - jx4       
            dy34             = iy3 - jy4       
            dz34             = iz3 - jz4       
            rsq34            = dx34*dx34+dy34*dy34+dz34*dz34
            dx42             = ix4 - jx2       
            dy42             = iy4 - jy2       
            dz42             = iz4 - jz2       
            rsq42            = dx42*dx42+dy42*dy42+dz42*dz42
            dx43             = ix4 - jx3       
            dy43             = iy4 - jy3       
            dz43             = iz4 - jz3       
            rsq43            = dx43*dx43+dy43*dy43+dz43*dz43
            dx44             = ix4 - jx4       
            dy44             = iy4 - jy4       
            dz44             = iz4 - jz4       
            rsq44            = dx44*dx44+dy44*dy44+dz44*dz44

C          Calculate 1/r and 1/r2

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv11           = frsqrtes(rsq11) 
#else
            rinv11           = frsqrte(dble(rsq11)) 
#endif
            rinv11           = (0.5*rinv11*(3.0-((rsq11*rinv11)
     &  *rinv11)))
#ifdef GMX_DOUBLE
            rinv11           = (0.5*rinv11*(3.0-((rsq11*rinv11)
     &  *rinv11)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv22           = frsqrtes(rsq22) 
#else
            rinv22           = frsqrte(dble(rsq22)) 
#endif
            rinv22           = (0.5*rinv22*(3.0-((rsq22*rinv22)
     &  *rinv22)))
#ifdef GMX_DOUBLE
            rinv22           = (0.5*rinv22*(3.0-((rsq22*rinv22)
     &  *rinv22)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv23           = frsqrtes(rsq23) 
#else
            rinv23           = frsqrte(dble(rsq23)) 
#endif
            rinv23           = (0.5*rinv23*(3.0-((rsq23*rinv23)
     &  *rinv23)))
#ifdef GMX_DOUBLE
            rinv23           = (0.5*rinv23*(3.0-((rsq23*rinv23)
     &  *rinv23)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv24           = frsqrtes(rsq24) 
#else
            rinv24           = frsqrte(dble(rsq24)) 
#endif
            rinv24           = (0.5*rinv24*(3.0-((rsq24*rinv24)
     &  *rinv24)))
#ifdef GMX_DOUBLE
            rinv24           = (0.5*rinv24*(3.0-((rsq24*rinv24)
     &  *rinv24)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv32           = frsqrtes(rsq32) 
#else
            rinv32           = frsqrte(dble(rsq32)) 
#endif
            rinv32           = (0.5*rinv32*(3.0-((rsq32*rinv32)
     &  *rinv32)))
#ifdef GMX_DOUBLE
            rinv32           = (0.5*rinv32*(3.0-((rsq32*rinv32)
     &  *rinv32)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv33           = frsqrtes(rsq33) 
#else
            rinv33           = frsqrte(dble(rsq33)) 
#endif
            rinv33           = (0.5*rinv33*(3.0-((rsq33*rinv33)
     &  *rinv33)))
#ifdef GMX_DOUBLE
            rinv33           = (0.5*rinv33*(3.0-((rsq33*rinv33)
     &  *rinv33)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv34           = frsqrtes(rsq34) 
#else
            rinv34           = frsqrte(dble(rsq34)) 
#endif
            rinv34           = (0.5*rinv34*(3.0-((rsq34*rinv34)
     &  *rinv34)))
#ifdef GMX_DOUBLE
            rinv34           = (0.5*rinv34*(3.0-((rsq34*rinv34)
     &  *rinv34)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv42           = frsqrtes(rsq42) 
#else
            rinv42           = frsqrte(dble(rsq42)) 
#endif
            rinv42           = (0.5*rinv42*(3.0-((rsq42*rinv42)
     &  *rinv42)))
#ifdef GMX_DOUBLE
            rinv42           = (0.5*rinv42*(3.0-((rsq42*rinv42)
     &  *rinv42)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv43           = frsqrtes(rsq43) 
#else
            rinv43           = frsqrte(dble(rsq43)) 
#endif
            rinv43           = (0.5*rinv43*(3.0-((rsq43*rinv43)
     &  *rinv43)))
#ifdef GMX_DOUBLE
            rinv43           = (0.5*rinv43*(3.0-((rsq43*rinv43)
     &  *rinv43)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv44           = frsqrtes(rsq44) 
#else
            rinv44           = frsqrte(dble(rsq44)) 
#endif
            rinv44           = (0.5*rinv44*(3.0-((rsq44*rinv44)
     &  *rinv44)))
#ifdef GMX_DOUBLE
            rinv44           = (0.5*rinv44*(3.0-((rsq44*rinv44)
     &  *rinv44)))
#endif

C          Load parameters for j atom

C          Calculate table index
            r                = rsq11*rinv11    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated VdW interaction - dispersion
            nnn              = nnn+4           
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            Vvdw6            = c6*VV           

C          Tabulated VdW interaction - repulsion
            nnn              = nnn+4           
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            Vvdw12           = c12*VV          
            Vvdwtot          = Vvdwtot+ Vvdw6 + Vvdw12

C          Load parameters for j atom
            qq               = qqHH            

C          Calculate table index
            r                = rsq22*rinv22    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            vcoul            = qq*VV           
            vctot            = vctot + vcoul   

C          Load parameters for j atom
            qq               = qqHH            

C          Calculate table index
            r                = rsq23*rinv23    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            vcoul            = qq*VV           
            vctot            = vctot + vcoul   

C          Load parameters for j atom
            qq               = qqMH            

C          Calculate table index
            r                = rsq24*rinv24    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            vcoul            = qq*VV           
            vctot            = vctot + vcoul   

C          Load parameters for j atom
            qq               = qqHH            

C          Calculate table index
            r                = rsq32*rinv32    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            vcoul            = qq*VV           
            vctot            = vctot + vcoul   

C          Load parameters for j atom
            qq               = qqHH            

C          Calculate table index
            r                = rsq33*rinv33    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            vcoul            = qq*VV           
            vctot            = vctot + vcoul   

C          Load parameters for j atom
            qq               = qqMH            

C          Calculate table index
            r                = rsq34*rinv34    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            vcoul            = qq*VV           
            vctot            = vctot + vcoul   

C          Load parameters for j atom
            qq               = qqMH            

C          Calculate table index
            r                = rsq42*rinv42    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            vcoul            = qq*VV           
            vctot            = vctot + vcoul   

C          Load parameters for j atom
            qq               = qqMH            

C          Calculate table index
            r                = rsq43*rinv43    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            vcoul            = qq*VV           
            vctot            = vctot + vcoul   

C          Load parameters for j atom
            qq               = qqMM            

C          Calculate table index
            r                = rsq44*rinv44    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 12*n0+1         

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            vcoul            = qq*VV           
            vctot            = vctot + vcoul   

C          Inner loop uses 268 flops/iteration
          end do
          

C        Add i forces to mem and shifted force list

C        Add potential energies to the group for this list
          ggid             = gid(n)+1        
          Vc(ggid)         = Vc(ggid) + vctot
          Vvdw(ggid)       = Vvdw(ggid) + Vvdwtot

C        Increment number of inner iterations
          ninner           = ninner + nj1 - nj0

C        Outer loop uses 14 flops/iteration
        end do
        

C      Increment number of outer iterations
        nouter           = nouter + nn1 - nn0
      if(nn1.lt.nri) goto 10

C    Write outer/inner iteration count to pointers
      outeriter        = nouter          
      inneriter        = ninner          
      return
      end




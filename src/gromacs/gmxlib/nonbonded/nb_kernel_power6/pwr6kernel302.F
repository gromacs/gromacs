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
C Gromacs nonbonded kernel pwr6kernel302
C Coulomb interaction:     Tabulated
C VdW interaction:         Not calculated
C water optimization:      pairs of SPC/TIP3P interactions
C Calculate forces:        yes
C
      subroutine pwr6kernel302(
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
      gmxreal       r,rt,eps,eps2
      integer*4     n0,nnn
      gmxreal       Y,F,Geps,Heps2,Fp,VV
      gmxreal       FF
      gmxreal       fijC
      gmxreal       ix1,iy1,iz1,fix1,fiy1,fiz1
      gmxreal       ix2,iy2,iz2,fix2,fiy2,fiz2
      gmxreal       ix3,iy3,iz3,fix3,fiy3,fiz3
      gmxreal       jx1,jy1,jz1,fjx1,fjy1,fjz1
      gmxreal       jx2,jy2,jz2,fjx2,fjy2,fjz2
      gmxreal       jx3,jy3,jz3,fjx3,fjy3,fjz3
      gmxreal       dx11,dy11,dz11,rsq11,rinv11
      gmxreal       dx12,dy12,dz12,rsq12,rinv12
      gmxreal       dx13,dy13,dz13,rsq13,rinv13
      gmxreal       dx21,dy21,dz21,rsq21,rinv21
      gmxreal       dx22,dy22,dz22,rsq22,rinv22
      gmxreal       dx23,dy23,dz23,rsq23,rinv23
      gmxreal       dx31,dy31,dz31,rsq31,rinv31
      gmxreal       dx32,dy32,dz32,rsq32,rinv32
      gmxreal       dx33,dy33,dz33,rsq33,rinv33
      gmxreal       qO,qH,qqOO,qqOH,qqHH


C    Initialize water data
      ii               = iinr(1)+1       
      qO               = charge(ii)      
      qH               = charge(ii+1)    
      qqOO             = facel*qO*qO     
      qqOH             = facel*qO*qH     
      qqHH             = facel*qH*qH     


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

C        Zero the potential energy for this list
          vctot            = 0               

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

C          Calculate distance
            dx11             = ix1 - jx1       
            dy11             = iy1 - jy1       
            dz11             = iz1 - jz1       
            rsq11            = dx11*dx11+dy11*dy11+dz11*dz11
            dx12             = ix1 - jx2       
            dy12             = iy1 - jy2       
            dz12             = iz1 - jz2       
            rsq12            = dx12*dx12+dy12*dy12+dz12*dz12
            dx13             = ix1 - jx3       
            dy13             = iy1 - jy3       
            dz13             = iz1 - jz3       
            rsq13            = dx13*dx13+dy13*dy13+dz13*dz13
            dx21             = ix2 - jx1       
            dy21             = iy2 - jy1       
            dz21             = iz2 - jz1       
            rsq21            = dx21*dx21+dy21*dy21+dz21*dz21
            dx22             = ix2 - jx2       
            dy22             = iy2 - jy2       
            dz22             = iz2 - jz2       
            rsq22            = dx22*dx22+dy22*dy22+dz22*dz22
            dx23             = ix2 - jx3       
            dy23             = iy2 - jy3       
            dz23             = iz2 - jz3       
            rsq23            = dx23*dx23+dy23*dy23+dz23*dz23
            dx31             = ix3 - jx1       
            dy31             = iy3 - jy1       
            dz31             = iz3 - jz1       
            rsq31            = dx31*dx31+dy31*dy31+dz31*dz31
            dx32             = ix3 - jx2       
            dy32             = iy3 - jy2       
            dz32             = iz3 - jz2       
            rsq32            = dx32*dx32+dy32*dy32+dz32*dz32
            dx33             = ix3 - jx3       
            dy33             = iy3 - jy3       
            dz33             = iz3 - jz3       
            rsq33            = dx33*dx33+dy33*dy33+dz33*dz33

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
            rinv12           = frsqrtes(rsq12) 
#else
            rinv12           = frsqrte(dble(rsq12)) 
#endif
            rinv12           = (0.5*rinv12*(3.0-((rsq12*rinv12)
     &  *rinv12)))
#ifdef GMX_DOUBLE
            rinv12           = (0.5*rinv12*(3.0-((rsq12*rinv12)
     &  *rinv12)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv13           = frsqrtes(rsq13) 
#else
            rinv13           = frsqrte(dble(rsq13)) 
#endif
            rinv13           = (0.5*rinv13*(3.0-((rsq13*rinv13)
     &  *rinv13)))
#ifdef GMX_DOUBLE
            rinv13           = (0.5*rinv13*(3.0-((rsq13*rinv13)
     &  *rinv13)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv21           = frsqrtes(rsq21) 
#else
            rinv21           = frsqrte(dble(rsq21)) 
#endif
            rinv21           = (0.5*rinv21*(3.0-((rsq21*rinv21)
     &  *rinv21)))
#ifdef GMX_DOUBLE
            rinv21           = (0.5*rinv21*(3.0-((rsq21*rinv21)
     &  *rinv21)))
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
            rinv31           = frsqrtes(rsq31) 
#else
            rinv31           = frsqrte(dble(rsq31)) 
#endif
            rinv31           = (0.5*rinv31*(3.0-((rsq31*rinv31)
     &  *rinv31)))
#ifdef GMX_DOUBLE
            rinv31           = (0.5*rinv31*(3.0-((rsq31*rinv31)
     &  *rinv31)))
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

C          Load parameters for j atom
            qq               = qqOO            

C          Calculate table index
            r                = rsq11*rinv11    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            fscal            = -((fijC)*tabscale)*rinv11

C          Calculate temporary vectorial force
            tx               = fscal*dx11      
            ty               = fscal*dy11      
            tz               = fscal*dz11      

C          Increment i atom force
            fix1             = fix1 + tx       
            fiy1             = fiy1 + ty       
            fiz1             = fiz1 + tz       

C          Decrement j atom force
            fjx1             = faction(j3+0) - tx
            fjy1             = faction(j3+1) - ty
            fjz1             = faction(j3+2) - tz

C          Load parameters for j atom
            qq               = qqOH            

C          Calculate table index
            r                = rsq12*rinv12    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            fscal            = -((fijC)*tabscale)*rinv12

C          Calculate temporary vectorial force
            tx               = fscal*dx12      
            ty               = fscal*dy12      
            tz               = fscal*dz12      

C          Increment i atom force
            fix1             = fix1 + tx       
            fiy1             = fiy1 + ty       
            fiz1             = fiz1 + tz       

C          Decrement j atom force
            fjx2             = faction(j3+3) - tx
            fjy2             = faction(j3+4) - ty
            fjz2             = faction(j3+5) - tz

C          Load parameters for j atom
            qq               = qqOH            

C          Calculate table index
            r                = rsq13*rinv13    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            fscal            = -((fijC)*tabscale)*rinv13

C          Calculate temporary vectorial force
            tx               = fscal*dx13      
            ty               = fscal*dy13      
            tz               = fscal*dz13      

C          Increment i atom force
            fix1             = fix1 + tx       
            fiy1             = fiy1 + ty       
            fiz1             = fiz1 + tz       

C          Decrement j atom force
            fjx3             = faction(j3+6) - tx
            fjy3             = faction(j3+7) - ty
            fjz3             = faction(j3+8) - tz

C          Load parameters for j atom
            qq               = qqOH            

C          Calculate table index
            r                = rsq21*rinv21    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            fscal            = -((fijC)*tabscale)*rinv21

C          Calculate temporary vectorial force
            tx               = fscal*dx21      
            ty               = fscal*dy21      
            tz               = fscal*dz21      

C          Increment i atom force
            fix2             = fix2 + tx       
            fiy2             = fiy2 + ty       
            fiz2             = fiz2 + tz       

C          Decrement j atom force
            fjx1             = fjx1 - tx       
            fjy1             = fjy1 - ty       
            fjz1             = fjz1 - tz       

C          Load parameters for j atom
            qq               = qqHH            

C          Calculate table index
            r                = rsq22*rinv22    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            fjx2             = fjx2 - tx       
            fjy2             = fjy2 - ty       
            fjz2             = fjz2 - tz       

C          Load parameters for j atom
            qq               = qqHH            

C          Calculate table index
            r                = rsq23*rinv23    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            fjx3             = fjx3 - tx       
            fjy3             = fjy3 - ty       
            fjz3             = fjz3 - tz       

C          Load parameters for j atom
            qq               = qqOH            

C          Calculate table index
            r                = rsq31*rinv31    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            fscal            = -((fijC)*tabscale)*rinv31

C          Calculate temporary vectorial force
            tx               = fscal*dx31      
            ty               = fscal*dy31      
            tz               = fscal*dz31      

C          Increment i atom force
            fix3             = fix3 + tx       
            fiy3             = fiy3 + ty       
            fiz3             = fiz3 + tz       

C          Decrement j atom force
            faction(j3+0)    = fjx1 - tx       
            faction(j3+1)    = fjy1 - ty       
            faction(j3+2)    = fjz1 - tz       

C          Load parameters for j atom
            qq               = qqHH            

C          Calculate table index
            r                = rsq32*rinv32    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            faction(j3+3)    = fjx2 - tx       
            faction(j3+4)    = fjy2 - ty       
            faction(j3+5)    = fjz2 - tz       

C          Load parameters for j atom
            qq               = qqHH            

C          Calculate table index
            r                = rsq33*rinv33    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            faction(j3+6)    = fjx3 - tx       
            faction(j3+7)    = fjy3 - ty       
            faction(j3+8)    = fjz3 - tz       

C          Inner loop uses 378 flops/iteration
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
          fshift(is3)      = fshift(is3)+fix1+fix2+fix3
          fshift(is3+1)    = fshift(is3+1)+fiy1+fiy2+fiy3
          fshift(is3+2)    = fshift(is3+2)+fiz1+fiz2+fiz3

C        Add potential energies to the group for this list
          ggid             = gid(n)+1        
          Vc(ggid)         = Vc(ggid) + vctot

C        Increment number of inner iterations
          ninner           = ninner + nj1 - nj0

C        Outer loop uses 28 flops/iteration
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
C Gromacs nonbonded kernel pwr6kernel302nf
C Coulomb interaction:     Tabulated
C VdW interaction:         Not calculated
C water optimization:      pairs of SPC/TIP3P interactions
C Calculate forces:        no
C
      subroutine pwr6kernel302nf(
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
      gmxreal       r,rt,eps,eps2
      integer*4     n0,nnn
      gmxreal       Y,F,Geps,Heps2,Fp,VV
      gmxreal       ix1,iy1,iz1
      gmxreal       ix2,iy2,iz2
      gmxreal       ix3,iy3,iz3
      gmxreal       jx1,jy1,jz1
      gmxreal       jx2,jy2,jz2
      gmxreal       jx3,jy3,jz3
      gmxreal       dx11,dy11,dz11,rsq11,rinv11
      gmxreal       dx12,dy12,dz12,rsq12,rinv12
      gmxreal       dx13,dy13,dz13,rsq13,rinv13
      gmxreal       dx21,dy21,dz21,rsq21,rinv21
      gmxreal       dx22,dy22,dz22,rsq22,rinv22
      gmxreal       dx23,dy23,dz23,rsq23,rinv23
      gmxreal       dx31,dy31,dz31,rsq31,rinv31
      gmxreal       dx32,dy32,dz32,rsq32,rinv32
      gmxreal       dx33,dy33,dz33,rsq33,rinv33
      gmxreal       qO,qH,qqOO,qqOH,qqHH


C    Initialize water data
      ii               = iinr(1)+1       
      qO               = charge(ii)      
      qH               = charge(ii+1)    
      qqOO             = facel*qO*qO     
      qqOH             = facel*qO*qH     
      qqHH             = facel*qH*qH     


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

C        Zero the potential energy for this list
          vctot            = 0               

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

C          Calculate distance
            dx11             = ix1 - jx1       
            dy11             = iy1 - jy1       
            dz11             = iz1 - jz1       
            rsq11            = dx11*dx11+dy11*dy11+dz11*dz11
            dx12             = ix1 - jx2       
            dy12             = iy1 - jy2       
            dz12             = iz1 - jz2       
            rsq12            = dx12*dx12+dy12*dy12+dz12*dz12
            dx13             = ix1 - jx3       
            dy13             = iy1 - jy3       
            dz13             = iz1 - jz3       
            rsq13            = dx13*dx13+dy13*dy13+dz13*dz13
            dx21             = ix2 - jx1       
            dy21             = iy2 - jy1       
            dz21             = iz2 - jz1       
            rsq21            = dx21*dx21+dy21*dy21+dz21*dz21
            dx22             = ix2 - jx2       
            dy22             = iy2 - jy2       
            dz22             = iz2 - jz2       
            rsq22            = dx22*dx22+dy22*dy22+dz22*dz22
            dx23             = ix2 - jx3       
            dy23             = iy2 - jy3       
            dz23             = iz2 - jz3       
            rsq23            = dx23*dx23+dy23*dy23+dz23*dz23
            dx31             = ix3 - jx1       
            dy31             = iy3 - jy1       
            dz31             = iz3 - jz1       
            rsq31            = dx31*dx31+dy31*dy31+dz31*dz31
            dx32             = ix3 - jx2       
            dy32             = iy3 - jy2       
            dz32             = iz3 - jz2       
            rsq32            = dx32*dx32+dy32*dy32+dz32*dz32
            dx33             = ix3 - jx3       
            dy33             = iy3 - jy3       
            dz33             = iz3 - jz3       
            rsq33            = dx33*dx33+dy33*dy33+dz33*dz33

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
            rinv12           = frsqrtes(rsq12) 
#else
            rinv12           = frsqrte(dble(rsq12)) 
#endif
            rinv12           = (0.5*rinv12*(3.0-((rsq12*rinv12)
     &  *rinv12)))
#ifdef GMX_DOUBLE
            rinv12           = (0.5*rinv12*(3.0-((rsq12*rinv12)
     &  *rinv12)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv13           = frsqrtes(rsq13) 
#else
            rinv13           = frsqrte(dble(rsq13)) 
#endif
            rinv13           = (0.5*rinv13*(3.0-((rsq13*rinv13)
     &  *rinv13)))
#ifdef GMX_DOUBLE
            rinv13           = (0.5*rinv13*(3.0-((rsq13*rinv13)
     &  *rinv13)))
#endif

C          PowerPC intrinsics 1/sqrt lookup table
#ifndef GMX_BLUEGENE
            rinv21           = frsqrtes(rsq21) 
#else
            rinv21           = frsqrte(dble(rsq21)) 
#endif
            rinv21           = (0.5*rinv21*(3.0-((rsq21*rinv21)
     &  *rinv21)))
#ifdef GMX_DOUBLE
            rinv21           = (0.5*rinv21*(3.0-((rsq21*rinv21)
     &  *rinv21)))
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
            rinv31           = frsqrtes(rsq31) 
#else
            rinv31           = frsqrte(dble(rsq31)) 
#endif
            rinv31           = (0.5*rinv31*(3.0-((rsq31*rinv31)
     &  *rinv31)))
#ifdef GMX_DOUBLE
            rinv31           = (0.5*rinv31*(3.0-((rsq31*rinv31)
     &  *rinv31)))
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

C          Load parameters for j atom
            qq               = qqOO            

C          Calculate table index
            r                = rsq11*rinv11    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            qq               = qqOH            

C          Calculate table index
            r                = rsq12*rinv12    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            qq               = qqOH            

C          Calculate table index
            r                = rsq13*rinv13    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            qq               = qqOH            

C          Calculate table index
            r                = rsq21*rinv21    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            r                = rsq22*rinv22    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            nnn              = 4*n0+1          

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
            qq               = qqOH            

C          Calculate table index
            r                = rsq31*rinv31    

C          Calculate table index
            rt               = r*tabscale      
            n0               = rt              
            eps              = rt-n0           
            eps2             = eps*eps         
            nnn              = 4*n0+1          

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
            nnn              = 4*n0+1          

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
            nnn              = 4*n0+1          

C          Tabulated coulomb interaction
            Y                = VFtab(nnn)      
            F                = VFtab(nnn+1)    
            Geps             = eps*VFtab(nnn+2)
            Heps2            = eps2*VFtab(nnn+3)
            Fp               = F+Geps+Heps2    
            VV               = Y+eps*Fp        
            vcoul            = qq*VV           
            vctot            = vctot + vcoul   

C          Inner loop uses 234 flops/iteration
          end do
          

C        Add i forces to mem and shifted force list

C        Add potential energies to the group for this list
          ggid             = gid(n)+1        
          Vc(ggid)         = Vc(ggid) + vctot

C        Increment number of inner iterations
          ninner           = ninner + nj1 - nj0

C        Outer loop uses 10 flops/iteration
        end do
        

C      Increment number of outer iterations
        nouter           = nouter + nn1 - nn0
      if(nn1.lt.nri) goto 10

C    Write outer/inner iteration count to pointers
      outeriter        = nouter          
      inneriter        = ninner          
      return
      end




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

C
C Gromacs nonbonded kernel f77dkernel203
C Coulomb interaction:     Reaction field
C VdW interaction:         Not calculated
C water optimization:      TIP4P - other atoms
C Calculate forces:        yes
C
      subroutine f77dkernel203(
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
      real*8        shiftvec(*),fshift(*),pos(*),faction(*)
      integer*4     gid(*),type(*),ntype
      real*8        charge(*),facel,krf,crf,Vc(*),vdwparam(*)
      real*8        Vvdw(*),tabscale,VFtab(*)
      real*8        invsqrta(*),dvda(*),gbtabscale,GBtab(*)
      integer*4     nthreads,count,mtx,outeriter,inneriter
      real*8        work(*)

      integer*4     n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid
      integer*4     nn0,nn1,nouter,ninner
      real*8        shX,shY,shZ
      real*8        fscal,tx,ty,tz
      real*8        rinvsq
      real*8        jq
      real*8        qq,vcoul,vctot
      real*8        krsq
      real*8        ix2,iy2,iz2,fix2,fiy2,fiz2
      real*8        ix3,iy3,iz3,fix3,fiy3,fiz3
      real*8        ix4,iy4,iz4,fix4,fiy4,fiz4
      real*8        jx1,jy1,jz1,fjx1,fjy1,fjz1
      real*8        dx21,dy21,dz21,rsq21,rinv21
      real*8        dx31,dy31,dz31,rsq31,rinv31
      real*8        dx41,dy41,dz41,rsq41,rinv41
      real*8        qH,qM


C    Initialize water data
      ii               = iinr(1)+1       
      qH               = facel*charge(ii+1)
      qM               = facel*charge(ii+3)


C    Reset outer and inner iteration counters
      nouter           = 0               
      ninner           = 0               

C    Loop over thread workunits
   10 call f77kernelsync(mtx,count,nri,nthreads,nn0,nn1)
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

C        Clear i atom forces
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

C          Calculate distance
            dx21             = ix2 - jx1       
            dy21             = iy2 - jy1       
            dz21             = iz2 - jz1       
            rsq21            = dx21*dx21+dy21*dy21+dz21*dz21
            dx31             = ix3 - jx1       
            dy31             = iy3 - jy1       
            dz31             = iz3 - jz1       
            rsq31            = dx31*dx31+dy31*dy31+dz31*dz31
            dx41             = ix4 - jx1       
            dy41             = iy4 - jy1       
            dz41             = iz4 - jz1       
            rsq41            = dx41*dx41+dy41*dy41+dz41*dz41

C          Calculate 1/r and 1/r2
            rinv21           = 1.0/sqrt(rsq21) 
            rinv31           = 1.0/sqrt(rsq31) 
            rinv41           = 1.0/sqrt(rsq41) 

C          Load parameters for j atom
            jq               = charge(jnr+0)   
            qq               = qH*jq           
            rinvsq           = rinv21*rinv21   

C          Coulomb reaction-field interaction
            krsq             = krf*rsq21       
            vcoul            = qq*(rinv21+krsq-crf)
            vctot            = vctot+vcoul     
            fscal            = (qq*(rinv21-2.0*krsq))*rinvsq

C          Calculate temporary vectorial force
            tx               = fscal*dx21      
            ty               = fscal*dy21      
            tz               = fscal*dz21      

C          Increment i atom force
            fix2             = fix2 + tx       
            fiy2             = fiy2 + ty       
            fiz2             = fiz2 + tz       

C          Decrement j atom force
            fjx1             = faction(j3+0) - tx
            fjy1             = faction(j3+1) - ty
            fjz1             = faction(j3+2) - tz

C          Load parameters for j atom
            rinvsq           = rinv31*rinv31   

C          Coulomb reaction-field interaction
            krsq             = krf*rsq31       
            vcoul            = qq*(rinv31+krsq-crf)
            vctot            = vctot+vcoul     
            fscal            = (qq*(rinv31-2.0*krsq))*rinvsq

C          Calculate temporary vectorial force
            tx               = fscal*dx31      
            ty               = fscal*dy31      
            tz               = fscal*dz31      

C          Increment i atom force
            fix3             = fix3 + tx       
            fiy3             = fiy3 + ty       
            fiz3             = fiz3 + tz       

C          Decrement j atom force
            fjx1             = fjx1 - tx       
            fjy1             = fjy1 - ty       
            fjz1             = fjz1 - tz       

C          Load parameters for j atom
            qq               = qM*jq           
            rinvsq           = rinv41*rinv41   

C          Coulomb reaction-field interaction
            krsq             = krf*rsq41       
            vcoul            = qq*(rinv41+krsq-crf)
            vctot            = vctot+vcoul     
            fscal            = (qq*(rinv41-2.0*krsq))*rinvsq

C          Calculate temporary vectorial force
            tx               = fscal*dx41      
            ty               = fscal*dy41      
            tz               = fscal*dz41      

C          Increment i atom force
            fix4             = fix4 + tx       
            fiy4             = fiy4 + ty       
            fiz4             = fiz4 + tz       

C          Decrement j atom force
            faction(j3+0)    = fjx1 - tx       
            faction(j3+1)    = fjy1 - ty       
            faction(j3+2)    = fjz1 - tz       

C          Inner loop uses 113 flops/iteration
          end do
          

C        Add i forces to mem and shifted force list
          faction(ii3+3)   = faction(ii3+3) + fix2
          faction(ii3+4)   = faction(ii3+4) + fiy2
          faction(ii3+5)   = faction(ii3+5) + fiz2
          faction(ii3+6)   = faction(ii3+6) + fix3
          faction(ii3+7)   = faction(ii3+7) + fiy3
          faction(ii3+8)   = faction(ii3+8) + fiz3
          faction(ii3+9)   = faction(ii3+9) + fix4
          faction(ii3+10)  = faction(ii3+10) + fiy4
          faction(ii3+11)  = faction(ii3+11) + fiz4
          fshift(is3)      = fshift(is3)+fix2+fix3+fix4
          fshift(is3+1)    = fshift(is3+1)+fiy2+fiy3+fiy4
          fshift(is3+2)    = fshift(is3+2)+fiz2+fiz3+fiz4

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
C Gromacs nonbonded kernel f77dkernel203nf
C Coulomb interaction:     Reaction field
C VdW interaction:         Not calculated
C water optimization:      TIP4P - other atoms
C Calculate forces:        no
C
      subroutine f77dkernel203nf(
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
      real*8        shiftvec(*),fshift(*),pos(*),faction(*)
      integer*4     gid(*),type(*),ntype
      real*8        charge(*),facel,krf,crf,Vc(*),vdwparam(*)
      real*8        Vvdw(*),tabscale,VFtab(*)
      real*8        invsqrta(*),dvda(*),gbtabscale,GBtab(*)
      integer*4     nthreads,count,mtx,outeriter,inneriter
      real*8        work(*)

      integer*4     n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid
      integer*4     nn0,nn1,nouter,ninner
      real*8        shX,shY,shZ
      real*8        jq
      real*8        qq,vcoul,vctot
      real*8        krsq
      real*8        ix2,iy2,iz2
      real*8        ix3,iy3,iz3
      real*8        ix4,iy4,iz4
      real*8        jx1,jy1,jz1
      real*8        dx21,dy21,dz21,rsq21,rinv21
      real*8        dx31,dy31,dz31,rsq31,rinv31
      real*8        dx41,dy41,dz41,rsq41,rinv41
      real*8        qH,qM


C    Initialize water data
      ii               = iinr(1)+1       
      qH               = facel*charge(ii+1)
      qM               = facel*charge(ii+3)


C    Reset outer and inner iteration counters
      nouter           = 0               
      ninner           = 0               

C    Loop over thread workunits
   10 call f77kernelsync(mtx,count,nri,nthreads,nn0,nn1)
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

C        Clear i atom forces
          
          do k=nj0,nj1

C          Get j neighbor index, and coordinate index
            jnr              = jjnr(k)+1       
            j3               = 3*jnr-2         

C          load j atom coordinates
            jx1              = pos(j3+0)       
            jy1              = pos(j3+1)       
            jz1              = pos(j3+2)       

C          Calculate distance
            dx21             = ix2 - jx1       
            dy21             = iy2 - jy1       
            dz21             = iz2 - jz1       
            rsq21            = dx21*dx21+dy21*dy21+dz21*dz21
            dx31             = ix3 - jx1       
            dy31             = iy3 - jy1       
            dz31             = iz3 - jz1       
            rsq31            = dx31*dx31+dy31*dy31+dz31*dz31
            dx41             = ix4 - jx1       
            dy41             = iy4 - jy1       
            dz41             = iz4 - jz1       
            rsq41            = dx41*dx41+dy41*dy41+dz41*dz41

C          Calculate 1/r and 1/r2
            rinv21           = 1.0/sqrt(rsq21) 
            rinv31           = 1.0/sqrt(rsq31) 
            rinv41           = 1.0/sqrt(rsq41) 

C          Load parameters for j atom
            jq               = charge(jnr+0)   
            qq               = qH*jq           

C          Coulomb reaction-field interaction
            krsq             = krf*rsq21       
            vcoul            = qq*(rinv21+krsq-crf)
            vctot            = vctot+vcoul     

C          Load parameters for j atom

C          Coulomb reaction-field interaction
            krsq             = krf*rsq31       
            vcoul            = qq*(rinv31+krsq-crf)
            vctot            = vctot+vcoul     

C          Load parameters for j atom
            qq               = qM*jq           

C          Coulomb reaction-field interaction
            krsq             = krf*rsq41       
            vcoul            = qq*(rinv41+krsq-crf)
            vctot            = vctot+vcoul     

C          Inner loop uses 71 flops/iteration
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




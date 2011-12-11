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
C Gromacs nonbonded kernel f77dkernel020
C Coulomb interaction:     Not calculated
C VdW interaction:         Buckingham
C water optimization:      No
C Calculate forces:        yes
C
      subroutine f77dkernel020(
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
      integer*4     nti
      integer*4     tj
      real*8        rinvsix
      real*8        Vvdw6,Vvdwtot
      real*8        Vvdwexp,br
      real*8        ix1,iy1,iz1,fix1,fiy1,fiz1
      real*8        jx1,jy1,jz1
      real*8        dx11,dy11,dz11,rsq11,rinv11
      real*8        c6,cexp1,cexp2


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
          ix1              = shX + pos(ii3+0)
          iy1              = shY + pos(ii3+1)
          iz1              = shZ + pos(ii3+2)

C        Load parameters for i atom
          nti              = 3*ntype*type(ii)

C        Zero the potential energy for this list
          Vvdwtot          = 0               

C        Clear i atom forces
          fix1             = 0               
          fiy1             = 0               
          fiz1             = 0               
          
          do k=nj0,nj1

C          Get j neighbor index, and coordinate index
            jnr              = jjnr(k)+1       
            j3               = 3*jnr-2         

C          load j atom coordinates
            jx1              = pos(j3+0)       
            jy1              = pos(j3+1)       
            jz1              = pos(j3+2)       

C          Calculate distance
            dx11             = ix1 - jx1       
            dy11             = iy1 - jy1       
            dz11             = iz1 - jz1       
            rsq11            = dx11*dx11+dy11*dy11+dz11*dz11

C          Calculate 1/r and 1/r2
            rinv11           = 1.0/sqrt(rsq11) 

C          Load parameters for j atom
            tj               = nti+3*type(jnr)+1
            c6               = vdwparam(tj)    
            cexp1            = vdwparam(tj+1)  
            cexp2            = vdwparam(tj+2)  
            rinvsq           = rinv11*rinv11   

C          Buckingham interaction
            rinvsix          = rinvsq*rinvsq*rinvsq
            Vvdw6            = c6*rinvsix      
            br               = cexp2*rsq11*rinv11
            Vvdwexp          = cexp1*exp(-br)  
            Vvdwtot          = Vvdwtot+Vvdwexp-Vvdw6
            fscal            = (br*Vvdwexp-6.0*Vvdw6)*rinvsq

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

C          Inner loop uses 66 flops/iteration
          end do
          

C        Add i forces to mem and shifted force list
          faction(ii3+0)   = faction(ii3+0) + fix1
          faction(ii3+1)   = faction(ii3+1) + fiy1
          faction(ii3+2)   = faction(ii3+2) + fiz1
          fshift(is3)      = fshift(is3)+fix1
          fshift(is3+1)    = fshift(is3+1)+fiy1
          fshift(is3+2)    = fshift(is3+2)+fiz1

C        Add potential energies to the group for this list
          ggid             = gid(n)+1        
          Vvdw(ggid)       = Vvdw(ggid) + Vvdwtot

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






C
C Gromacs nonbonded kernel f77dkernel020nf
C Coulomb interaction:     Not calculated
C VdW interaction:         Buckingham
C water optimization:      No
C Calculate forces:        no
C
      subroutine f77dkernel020nf(
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
      real*8        rinvsq
      integer*4     nti
      integer*4     tj
      real*8        rinvsix
      real*8        Vvdw6,Vvdwtot
      real*8        Vvdwexp,br
      real*8        ix1,iy1,iz1
      real*8        jx1,jy1,jz1
      real*8        dx11,dy11,dz11,rsq11,rinv11
      real*8        c6,cexp1,cexp2


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
          ix1              = shX + pos(ii3+0)
          iy1              = shY + pos(ii3+1)
          iz1              = shZ + pos(ii3+2)

C        Load parameters for i atom
          nti              = 3*ntype*type(ii)

C        Zero the potential energy for this list
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

C          Calculate distance
            dx11             = ix1 - jx1       
            dy11             = iy1 - jy1       
            dz11             = iz1 - jz1       
            rsq11            = dx11*dx11+dy11*dy11+dz11*dz11

C          Calculate 1/r and 1/r2
            rinv11           = 1.0/sqrt(rsq11) 

C          Load parameters for j atom
            tj               = nti+3*type(jnr)+1
            c6               = vdwparam(tj)    
            cexp1            = vdwparam(tj+1)  
            cexp2            = vdwparam(tj+2)  
            rinvsq           = rinv11*rinv11   

C          Buckingham interaction
            rinvsix          = rinvsq*rinvsq*rinvsq
            Vvdw6            = c6*rinvsix      
            br               = cexp2*rsq11*rinv11
            Vvdwexp          = cexp1*exp(-br)  
            Vvdwtot          = Vvdwtot+Vvdwexp-Vvdw6

C          Inner loop uses 53 flops/iteration
          end do
          

C        Add i forces to mem and shifted force list

C        Add potential energies to the group for this list
          ggid             = gid(n)+1        
          Vvdw(ggid)       = Vvdw(ggid) + Vvdwtot

C        Increment number of inner iterations
          ninner           = ninner + nj1 - nj0

C        Outer loop uses 4 flops/iteration
        end do
        

C      Increment number of outer iterations
        nouter           = nouter + nn1 - nn0
      if(nn1.lt.nri) goto 10

C    Write outer/inner iteration count to pointers
      outeriter        = nouter          
      inneriter        = ninner          
      return
      end




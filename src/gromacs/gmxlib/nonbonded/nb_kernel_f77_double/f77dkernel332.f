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
C Gromacs nonbonded kernel f77dkernel332
C Coulomb interaction:     Tabulated
C VdW interaction:         Tabulated
C water optimization:      pairs of SPC/TIP3P interactions
C Calculate forces:        yes
C
      subroutine f77dkernel332(
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
      real*8        qq,vcoul,vctot
      integer*4     tj
      real*8        Vvdw6,Vvdwtot
      real*8        Vvdw12
      real*8        r,rt,eps,eps2
      integer*4     n0,nnn
      real*8        Y,F,Geps,Heps2,Fp,VV
      real*8        FF
      real*8        fijC
      real*8        fijD,fijR
      real*8        ix1,iy1,iz1,fix1,fiy1,fiz1
      real*8        ix2,iy2,iz2,fix2,fiy2,fiz2
      real*8        ix3,iy3,iz3,fix3,fiy3,fiz3
      real*8        jx1,jy1,jz1,fjx1,fjy1,fjz1
      real*8        jx2,jy2,jz2,fjx2,fjy2,fjz2
      real*8        jx3,jy3,jz3,fjx3,fjy3,fjz3
      real*8        dx11,dy11,dz11,rsq11,rinv11
      real*8        dx12,dy12,dz12,rsq12,rinv12
      real*8        dx13,dy13,dz13,rsq13,rinv13
      real*8        dx21,dy21,dz21,rsq21,rinv21
      real*8        dx22,dy22,dz22,rsq22,rinv22
      real*8        dx23,dy23,dz23,rsq23,rinv23
      real*8        dx31,dy31,dz31,rsq31,rinv31
      real*8        dx32,dy32,dz32,rsq32,rinv32
      real*8        dx33,dy33,dz33,rsq33,rinv33
      real*8        qO,qH,qqOO,qqOH,qqHH
      real*8        c6,c12


C    Initialize water data
      ii               = iinr(1)+1       
      qO               = charge(ii)      
      qH               = charge(ii+1)    
      qqOO             = facel*qO*qO     
      qqOH             = facel*qO*qH     
      qqHH             = facel*qH*qH     
      tj               = 2*(ntype+1)*type(ii)+1
      c6               = vdwparam(tj)    
      c12              = vdwparam(tj+1)  


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
          ix2              = shX + pos(ii3+3)
          iy2              = shY + pos(ii3+4)
          iz2              = shZ + pos(ii3+5)
          ix3              = shX + pos(ii3+6)
          iy3              = shY + pos(ii3+7)
          iz3              = shZ + pos(ii3+8)

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
            rinv11           = 1.0/sqrt(rsq11) 
            rinv12           = 1.0/sqrt(rsq12) 
            rinv13           = 1.0/sqrt(rsq13) 
            rinv21           = 1.0/sqrt(rsq21) 
            rinv22           = 1.0/sqrt(rsq22) 
            rinv23           = 1.0/sqrt(rsq23) 
            rinv31           = 1.0/sqrt(rsq31) 
            rinv32           = 1.0/sqrt(rsq32) 
            rinv33           = 1.0/sqrt(rsq33) 

C          Load parameters for j atom
            qq               = qqOO            

C          Calculate table index
            r                = rsq11*rinv11    

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
            fscal            = -((fijC+fijD+fijR)*tabscale)*rinv11

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
            faction(j3+6)    = fjx3 - tx       
            faction(j3+7)    = fjy3 - ty       
            faction(j3+8)    = fjz3 - tz       

C          Inner loop uses 440 flops/iteration
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
          Vvdw(ggid)       = Vvdw(ggid) + Vvdwtot

C        Increment number of inner iterations
          ninner           = ninner + nj1 - nj0

C        Outer loop uses 29 flops/iteration
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
C Gromacs nonbonded kernel f77dkernel332nf
C Coulomb interaction:     Tabulated
C VdW interaction:         Tabulated
C water optimization:      pairs of SPC/TIP3P interactions
C Calculate forces:        no
C
      subroutine f77dkernel332nf(
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
      real*8        qq,vcoul,vctot
      integer*4     tj
      real*8        Vvdw6,Vvdwtot
      real*8        Vvdw12
      real*8        r,rt,eps,eps2
      integer*4     n0,nnn
      real*8        Y,F,Geps,Heps2,Fp,VV
      real*8        ix1,iy1,iz1
      real*8        ix2,iy2,iz2
      real*8        ix3,iy3,iz3
      real*8        jx1,jy1,jz1
      real*8        jx2,jy2,jz2
      real*8        jx3,jy3,jz3
      real*8        dx11,dy11,dz11,rsq11,rinv11
      real*8        dx12,dy12,dz12,rsq12,rinv12
      real*8        dx13,dy13,dz13,rsq13,rinv13
      real*8        dx21,dy21,dz21,rsq21,rinv21
      real*8        dx22,dy22,dz22,rsq22,rinv22
      real*8        dx23,dy23,dz23,rsq23,rinv23
      real*8        dx31,dy31,dz31,rsq31,rinv31
      real*8        dx32,dy32,dz32,rsq32,rinv32
      real*8        dx33,dy33,dz33,rsq33,rinv33
      real*8        qO,qH,qqOO,qqOH,qqHH
      real*8        c6,c12


C    Initialize water data
      ii               = iinr(1)+1       
      qO               = charge(ii)      
      qH               = charge(ii+1)    
      qqOO             = facel*qO*qO     
      qqOH             = facel*qO*qH     
      qqHH             = facel*qH*qH     
      tj               = 2*(ntype+1)*type(ii)+1
      c6               = vdwparam(tj)    
      c12              = vdwparam(tj+1)  


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
          ix2              = shX + pos(ii3+3)
          iy2              = shY + pos(ii3+4)
          iz2              = shZ + pos(ii3+5)
          ix3              = shX + pos(ii3+6)
          iy3              = shY + pos(ii3+7)
          iz3              = shZ + pos(ii3+8)

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
            rinv11           = 1.0/sqrt(rsq11) 
            rinv12           = 1.0/sqrt(rsq12) 
            rinv13           = 1.0/sqrt(rsq13) 
            rinv21           = 1.0/sqrt(rsq21) 
            rinv22           = 1.0/sqrt(rsq22) 
            rinv23           = 1.0/sqrt(rsq23) 
            rinv31           = 1.0/sqrt(rsq31) 
            rinv32           = 1.0/sqrt(rsq32) 
            rinv33           = 1.0/sqrt(rsq33) 

C          Load parameters for j atom
            qq               = qqOO            

C          Calculate table index
            r                = rsq11*rinv11    

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
            qq               = qqOH            

C          Calculate table index
            r                = rsq12*rinv12    

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
            qq               = qqOH            

C          Calculate table index
            r                = rsq13*rinv13    

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
            qq               = qqOH            

C          Calculate table index
            r                = rsq21*rinv21    

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
            qq               = qqOH            

C          Calculate table index
            r                = rsq31*rinv31    

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

C          Inner loop uses 286 flops/iteration
          end do
          

C        Add i forces to mem and shifted force list

C        Add potential energies to the group for this list
          ggid             = gid(n)+1        
          Vc(ggid)         = Vc(ggid) + vctot
          Vvdw(ggid)       = Vvdw(ggid) + Vvdwtot

C        Increment number of inner iterations
          ninner           = ninner + nj1 - nj0

C        Outer loop uses 11 flops/iteration
        end do
        

C      Increment number of outer iterations
        nouter           = nouter + nn1 - nn0
      if(nn1.lt.nri) goto 10

C    Write outer/inner iteration count to pointers
      outeriter        = nouter          
      inneriter        = ninner          
      return
      end




C
C        @(#) in_loopf.f 1.17 18 Aug 1996
C 
C        This  source-code  is  part  of
C 
C        G    R    O    M    A    C    S
C 
C  GROningen MAchine for Chemical Simulations
C 
C  Copyright (c) 1990-1995,
C  BIOSON Research Institute, Dept. of Biophysical Chemistry,
C  University of Groningen, The Netherlands
C 
C  Please refer to:
C  GROMACS: A Message Passing Parallel Molecular Dynamics Implementation
C  H.J.C. Berendsen, D. van der Spoel and R. van Drunen
C  Comp. Phys. Comm. 1995 (in press)
C 
C  Also check out our WWW page:
C  http://rugmd0.chem.rug.nl/~gmx/gmx.cgi
C  or e-mail to:
C  gromacs@chem.rug.nl
C 
C  And Hey:
C  Gyas ROwers Mature At Cryogenic Speed
C

C
C     This code is meant to be called from C routines.
C     Therefore all indices start at 0, although the arrays
C     start at 1, if an array contains an index we must add 1 to it.
C     EG: jjnr points to particles starting at 0
C         type is indexed from 1 to ...
C

      subroutine FORCOUL(ix,iy,iz,qi,
     $     pos,nj,jjnr,charge,
     $     faction,fip,
     $     Vc)
      
      real      ix,iy,iz,qi,Vc
      real      pos(*),charge(*),faction(*),fip(3)
      integer*4 nj,jjnr(*)
      
      integer   k,jnr,j3
      real      fX,fY,fZ
      real      rijX,rijY,rijZ
      real      fijscal,rsq,vijcoul,vctot
      real      rinv1,rinv2
      real      fjx,fjy,fjz
      real      tx,ty,tz

      fX    = 0
      fY    = 0
      fZ    = 0
      vctot = 0
      
cray compiler directive ignore vector dependencies      
c$dir ivdep
      do k=1,nj
         jnr         = jjnr(k)+1
         j3          = 3*jnr-2
         fjx         = faction(j3)
         fjy         = faction(j3+1)
         fjz         = faction(j3+2)
         rijX        = ix - pos(j3)
         rijY        = iy - pos(j3+1)
         rijZ        = iz - pos(j3+2)
         rsq         = (rijX*rijX)+(rijY*rijY)+(rijZ*rijZ)

         rinv1       = 1.0/sqrt(rsq)
         vijcoul     = qi*charge(jnr)*rinv1
         rinv2       = rinv1*rinv1
         vctot       = vctot+vijcoul
         fijscal     = vijcoul*rinv2
         
         tx       = rijX*fijscal
         ty       = rijY*fijscal
         tz       = rijZ*fijscal
         fX       = fX+tx
         fY       = fY+ty
         fZ       = fZ+tz
         faction(j3)   = fjx-tx
         faction(j3+1) = fjy-ty
         faction(j3+2) = fjz-tz
      end do

      fip(1) = fX
      fip(2) = fY
      fip(3) = fZ
      Vc     = Vc+vctot
      
      return
      
      end
      
      subroutine FORLJC(ix,iy,iz,qi,
     $     pos,nj,type,jjnr,charge,nbfp,
     $     faction,fip,
     $     Vc,Vnb)
      
      implicit none
      
      real      ix,iy,iz,qi
      real      pos(*),charge(*),faction(*),fip(3)
      integer*4 nj,jjnr(*),type(*)
      real      Vc,Vnb,nbfp(*)
      
      integer   k,jnr,j3,tj
      real      twelve,six
      real      fX,fY,fZ
      real      rijX,rijY,rijZ
      real      fijscal,vijcoul
      real      vctot,vnbtot
      real      rinv1,rinv2,rinv6
      real      fjx,fjy,fjz
      real      tx,ty,tz,vnb6,vnb12

      parameter(twelve=12.0,six=6.0)
            
      fX     = 0
      fY     = 0
      fZ     = 0
      vctot  = 0
      vnbtot = 0
      
cray compiler directive ignore vector dependencies      
c$dir ivdep
      do k=1,nj
         jnr   = jjnr(k)+1
         j3    = 3*jnr-2
         rijX  = ix - pos(j3)
         rijY  = iy - pos(j3+1)
         rijZ  = iz - pos(j3+2)

         rinv1       = 1.0/sqrt((rijX*rijX)+(rijY*rijY)+(rijZ*rijZ))
         rinv2       = rinv1*rinv1
         rinv6       = rinv2*rinv2*rinv2
         
         tj          = 2*type(jnr)+1
         vnb6        = nbfp(tj)*rinv6
         vnb12       = nbfp(tj+1)*rinv6*rinv6
         vijcoul     = qi*charge(jnr)*rinv1
         
         vctot       = vctot+vijcoul
         vnbtot      = vnbtot+vnb12-vnb6
         fijscal     = (twelve*vnb12-six*vnb6+vijcoul)*rinv2
         
         fjx           = faction(j3)
         tx            = rijX*fijscal
         fX            = fX + tx
         faction(j3)   = fjx - tx
         fjy           = faction(j3+1)
         ty            = rijY*fijscal
         fY            = fY + ty
         faction(j3+1) = fjy - ty
         fjz           = faction(j3+2)
         tz            = rijZ*fijscal
         fZ            = fZ + tz
         faction(j3+2) = fjz - tz
         
      end do
 
      fip(1) = fX
      fip(2) = fY
      fip(3) = fZ
      Vc     = Vc  + vctot
      Vnb    = Vnb + vnbtot

      return
      
      end
      
      subroutine FORWATER(i0,xw,fudge,pos,nj,type,jjnr,
     &     charge,nbfp,faction,fw,Vc,Vnb)

      integer*4 i0,nj,type(*),jjnr(*)
      real      xw(9),fudge
      real      pos(*),charge(*),nbfp(*)
      real      faction(*),fw(*),Vc,Vnb
     
      real      qO,qH
      
      qO   = charge(i0+1)*fudge
      qH   = charge(i0+2)*fudge
      
      call FORLJC(xw(1),xw(2),xw(3),qO,
     $     pos,nj,type,jjnr,charge,nbfp,
     $     faction,fw(1),
     $     Vc,Vnb)

      call FORCOUL(xw(4),xw(5),xw(6),qH,
     $     pos,nj,jjnr,charge,
     $     faction,fw(4),
     $     Vc)
     
      call FORCOUL(xw(7),xw(8),xw(9),qH,
     $     pos,nj,jjnr,charge,
     $     faction,fw(7),
     $     Vc)
      
      return
      
      end
      
      subroutine FORWCOUL(i0,xw,fudge,pos,nj,jjnr,
     &     charge,faction,fw,Vc)
      
      integer*4 i0,nj,jjnr(*)
      real      xw(9),fudge
      real      pos(*),charge(*)
      real      faction(*),fw(*),Vc
      
      integer*4 m,i3
      real      qi
      
c
c     Call the coulomb routine three times with the same neighbor list
c      
      do m=1,3
         qi = charge(i0+m)*fudge
         i3 = 3*m-2
         call FORCOUL(xw(i3),xw(i3+1),xw(i3+2),qi,
     $        pos,nj,jjnr,charge,
     $        faction,fw(i3),
     $        Vc)
      end do
      
      return
      
      end
      
      subroutine FORFREE(ix,iy,iz,inr,
     &     pos,nj,jjnr,typeA,typeB,fudge,chargeA,chargeB,
     &     nbfpA,nbfpB,faction,fip,Vc,Vnb,lambda,dvdlambda,
     &     krf,crf,tfac,trunctab)
               
      real      ix, iy, iz, pos(*), fudge
      integer*4 inr, nj, jjnr(*), typeA(*), typeB(*)
      real      chargeA(*), chargeB(*)
      real      nbfpA(*), nbfpB(*), faction(*), fip(*)
      real      Vc, Vnb, lambda, dvdlambda,krf,crf,tfac,trunctab(*)

      real       twelve, six
      integer    k,jnr,j3,tjA,tjB,m
      real       fX,fY,fZ
      real       rijX,rijY,rijZ
      real       fijscal,rsq,vijcoul,vctot,vnbtot
      real       rinv1,rinv2,rinv6,rinv12,fjx,fjy,fjz
      real       tx,ty,tz,vnb6,vnb12,tt,ffj,L1,dvdl
      real       qi,qiA,qiB,qj,qjA,qjB,dqi,c6,c12,c6a,c6b,c12a,c12b
  
      twelve = 12.0
      six    =  6.0
  
      fX     = 0.0
      fY     = 0.0
      fZ     = 0.0
      vctot  = 0.0
      vnbtot = 0.0
      L1     = 1.0-lambda
      dvdl   = 0.0
  
      qiA    = fudge*chargeA(inr+1)
      qiB    = fudge*chargeB(inr+1)
      qi     = L1*qiA+lambda*qiB
      dqi    = qiB-qiA
  
      do k=1,nj 
         jnr            = jjnr(k)+1
         j3             = 3*jnr-2
         rijX           = ix - pos(j3)
         rijY           = iy - pos(j3+1)
         rijZ           = iz - pos(j3+2)
         
         rinv1          = 1.0/sqrt((rijX*rijX)+(rijY*rijY)+(rijZ*rijZ))
         rinv2          = rinv1*rinv1
         rinv6          = rinv2*rinv2*rinv2
    
         qjA            = chargeA(jnr)
         qjB            = chargeB(jnr)
         qj             = L1*qjA+lambda*qjB
         
         dvdl           = dvdl + ((qjB-qjA)*qi+(dqi*qj))*rinv1
         
         tjA            = 2*typeA(jnr)+1
         tjB            = 2*typeB(jnr)+1
         
         c6a            = nbfpA(tjA)
         c6b            = nbfpB(tjB)
         c12a           = nbfpA(tjA+1)
         c12b           = nbfpB(tjB+1)
         c6             = L1*c6a  + lambda*c6b
         c12            = L1*c12a + lambda*c12b
         vnb6           = c6  * rinv6
         rinv12         = rinv6*rinv6
         vnb12          = c12 * rinv12
         dvdl           = dvdl + (c12b-c12a)*rinv12 - (c6b-c6a)*rinv6
         
         vijcoul        = qi*qj*rinv1
         vctot          = vctot+vijcoul
         vnbtot         = vnbtot+vnb12-vnb6
         fijscal        = (twelve*vnb12-six*vnb6+vijcoul)*rinv2
         
         fjx            = faction(j3+XX)
         tx             = rijX*fijscal
         fX             = fX + tx
         faction(j3+XX) = fjx - tx
         fjy            = faction(j3+YY)
         ty             = rijY*fijscal
         fY             = fY + ty
         faction(j3+YY) = fjy - ty
         fjz            = faction(j3+ZZ)
         tz             = rijZ*fijscal
         fZ             = fZ + tz
         faction(j3+ZZ) = fjz - tz
      end do
      
      fip(1)     = fX
      fip(2)     = fY
      fip(3)     = fZ
      Vc         = Vc+vctot
      Vnb        = Vnb+vnbtot
      dvdlambda  = dvdlambda+dvdl

      return
      
      end

      subroutine EXTRACT(nn,VFtab,eps,eps2,VV,FF)
      
      integer*4 nn
      real      VFtab(*),eps,eps2,VV,FF
      
      real      Y,F,Geps,Heps2,Fp
      
      Y     = VFtab(nn)
      F     = VFtab(nn+1)
      Geps  = VFtab(nn+2)*eps
      Heps2 = VFtab(nn+3)*eps2
      Fp    = F+Geps+Heps2
      
      VV    = Y+eps*Fp
      FF    = Fp+Geps+2.0*Heps2
      
      return
      end
      
      subroutine FORTAB(ix,iy,iz,qi,
     $     pos,nj,type,jjnr,charge,nbfp,
     $     faction,fip,
     $     Vc,Vnb,ntab,tabscale,
     $     VFtab)
      
      implicit none
      
      real      ix,iy,iz,qi
      real      pos(*),charge(*),faction(*),fip(3)
      integer*4 nj,jjnr(*),type(*),ntab
      real      Vc,Vnb,nbfp(*),tabscale
      real      VFtab(*)
      
      integer     k,jnr,j3,tj
      real        fX,fY,fZ
      real        rijX,rijY,rijZ
      real        vijcoul,fijD,fijR,fijC,fijscal
      real        fjx,fjy,fjz
      real        tx,ty,tz,vnb6,vnb12
      real        vctot,vnbtot
      real        qq,c6,c12,rsq
      real        r1,r1t,h_1,eps,VV,FF
      real        Y,F,G,epsH,eGeH,FeGeH
      integer     n0,n1,nn
  
      fX     = 0
      fY     = 0
      fZ     = 0
      vctot  = 0
      vnbtot = 0
      h_1    = tabscale
      
      do k=1,nj
         jnr            = jjnr(k)+1
         j3             = 3*jnr-2
         rijX           = ix - pos(j3)
         rijY           = iy - pos(j3+1)
         rijZ           = iz - pos(j3+2)
         
         rsq            = (rijX*rijX)+(rijY*rijY)+(rijZ*rijZ)
         r1             = sqrt(rsq)
         r1t            = r1*tabscale
         n0             = r1t
         n1             = 12*n0+1
         eps            = (r1t-n0)

C     Coulomb
C     call EXTRACT(n1,VFtab,eps,eps2,VV,FF)
         nn    = n1
         Y     = VFtab(nn)
         F     = VFtab(nn+1)
         G     = VFtab(nn+2)
         epsH  = VFtab(nn+3)*eps
         eGeH  = eps*(G+epsH)
         FeGeH = F+eGeH

         VV    = Y+eps*FeGeH
         FF    = FeGeH+eGeH+epsH

         qq       = qi*charge(jnr)
         vijcoul  = qq*VV
         fijC     = qq*FF
         vctot    = vctot  + vijcoul
         
C     Dispersion 
C     call EXTRACT(n1+4,VFtab,eps,eps2,VV,FF)
         nn    = n1+4
         Y     = VFtab(nn)
         F     = VFtab(nn+1)
         G     = VFtab(nn+2)
         epsH  = VFtab(nn+3)*eps
         eGeH  = eps*(G+epsH)
         FeGeH = F+eGeH
         
         VV    = Y+eps*FeGeH
         FF    = FeGeH+eGeH+epsH

         tj    = 2*type(jnr)+1
         c6    = nbfp(tj)
         vnb6  = c6*VV
         fijD  = c6*FF
         
C     Repulsion 
C     call EXTRACT(n1+8,VFtab,eps,eps2,VV,FF)
         nn    = n1+8
         Y     = VFtab(nn)
         F     = VFtab(nn+1)
         G     = VFtab(nn+2)
         epsH  = VFtab(nn+3)*eps
         eGeH  = eps*(G+epsH)
         FeGeH = F+eGeH

         VV    = Y+eps*FeGeH
         FF    = FeGeH+eGeH+epsH

         c12   = nbfp(tj+1)
         vnb12 = c12*VV
         fijR  = c12*FF
         vnbtot= vnbtot + vnb12 + vnb6
         
C     Total force
         fijscal        = -(fijD + fijR + fijC)*h_1/r1
         
         fjx            = faction(j3)
         tx             = rijX*fijscal
         fX             = fX + tx
         faction(j3)    = fjx - tx
         fjy            = faction(j3+1)
         ty             = rijY*fijscal
         fY             = fY + ty
         faction(j3+1) = fjy - ty
         fjz            = faction(j3+2)
         tz             = rijZ*fijscal
         fZ             = fZ + tz
         faction(j3+2) = fjz - tz
      end do
      
      fip(1) = fX
      fip(2) = fY
      fip(3) = fZ
      Vc   = Vc+vctot
      Vnb  = Vnb+vnbtot
         
      return
      end

      subroutine FORCOULTAB(ix,iy,iz,qi,
     $     pos,nj,type,jjnr,charge,nbfp,
     $     faction,fip,
     $     Vc,Vnb,ntab,tabscale,
     $     VFtab)
      
      implicit none
      
      real      ix,iy,iz,qi
      real      pos(*),charge(*),faction(*),fip(3)
      integer*4 nj,jjnr(*),type(*),ntab
      real      Vc,Vnb,nbfp(*),tabscale
      real      VFtab(*)
      
      integer     k,jnr,j3,tj
      real        fX,fY,fZ
      real        rijX,rijY,rijZ
      real        vijcoul,fijD,fijR,fijC,fijscal
      real        fjx,fjy,fjz
      real        tx,ty,tz,vnb6,vnb12
      real        vctot,vnbtot
      real        qq,c6,c12,rsq
      real        r1,r1t,h_1,eps,VV,FF
      real        Y,F,G,epsH,eGeH,FeGeH
      integer     n0,n1,nn
  
      fX     = 0
      fY     = 0
      fZ     = 0
      vctot  = 0
      h_1    = tabscale
      
      do k=1,nj
         jnr            = jjnr(k)+1
         j3             = 3*jnr-2
         rijX           = ix - pos(j3)
         rijY           = iy - pos(j3+1)
         rijZ           = iz - pos(j3+2)
         
         rsq            = (rijX*rijX)+(rijY*rijY)+(rijZ*rijZ)
         r1             = sqrt(rsq)
         r1t            = r1*tabscale
         n0             = r1t
         n1             = 12*n0+1
         eps            = (r1t-n0)

C     Coulomb
C         call EXTRACT(n1,VFtab,eps,eps2,VV,FF)
         nn    = n1
         Y     = VFtab(nn)
         F     = VFtab(nn+1)
         G     = VFtab(nn+2)
         epsH  = VFtab(nn+3)*eps
         eGeH  = eps*(G+epsH)
         FeGeH = F+eGeH

         VV    = Y+eps*FeGeH
         FF    = FeGeH+eGeH+epsH

         qq             = qi*charge(jnr)
         vijcoul        = qq*VV
         fijC           = qq*FF
         vctot          = vctot  + vijcoul
         
C     Total force
         fijscal        = -fijC*h_1/r1
         
         fjx            = faction(j3)
         tx             = rijX*fijscal
         fX             = fX + tx
         faction(j3) = fjx - tx
         fjy            = faction(j3+1)
         ty             = rijY*fijscal
         fY             = fY + ty
         faction(j3+1) = fjy - ty
         fjz            = faction(j3+2)
         tz             = rijZ*fijscal
         fZ             = fZ + tz
         faction(j3+2) = fjz - tz
      end do
      
      fip(1) = fX
      fip(2) = fY
      fip(3) = fZ
      Vc   = Vc+vctot
         
      return
      end


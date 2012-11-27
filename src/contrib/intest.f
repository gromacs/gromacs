C
C This file is part of the GROMACS molecular simulation package.
C
C Copyright (c) 1991-2001
C BIOSON Research Institute, Dept. of Biophysical Chemistry
C University of Groningen, The Netherlands
C Copyright (c) 2012, by the GROMACS development team, led by
C David van der Spoel, Berk Hess, Erik Lindahl, and including many
C others, as listed in the AUTHORS file in the top-level source
C directory and at http://www.gromacs.org.
C
C GROMACS is free software; you can redistribute it and/or
C modify it under the terms of the GNU Lesser General Public License
C as published by the Free Software Foundation; either version 2.1
C of the License, or (at your option) any later version.
C
C GROMACS is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C Lesser General Public License for more details.
C
C You should have received a copy of the GNU Lesser General Public
C License along with GROMACS; if not, see
C http://www.gnu.org/licenses, or write to the Free Software Foundation,
C Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
C
C If you want to redistribute modifications to GROMACS, please
C consider that scientific software is very special. Version
C control is crucial - bugs must be traceable. We will be happy to
C consider code for inclusion in the official distribution, but
C derived work must not be called official GROMACS. Details are found
C in the README & COPYING files - if they are missing, get the
C official version at http://www.gromacs.org.
C
C To help us fund GROMACS development, we humbly ask that you cite
C the research papers on the package. Check out http://www.gromacs.org.

C
C     This code is meant to be called from C routines.
C     Therefore all indices start at 0, although the arrays
C     start at 1, if an array contains an index we must add 1 to it.
C     EG: jjnr points to particles starting at 0
C         type is indexed from 1 to ...
C

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
      
      program main
      
      implicit none
      
      integer maxatom,maxx,maxlist,maxtype
      parameter(maxatom=1000,maxx=3*maxatom,maxlist=100)
      parameter(maxtype=19)
      
      real*4 ix,iy,iz,qi,pos(maxx),faction(maxx),fip(3)
      real*4 charge(maxatom),nbfp(2*maxtype),Vc,Vnb
      integer type(maxatom),jjnr(maxlist)
      integer i,i3,j
      
c     setup benchmark
      do i=1,maxtype
         nbfp(2*i-1) = 1e-3
         nbfp(2*i)   = 1e-6
      end do
      
      type(1)=1
      do i=2,maxatom
         type(i)=1+mod(type(i-1)+91,maxtype)
      end do
      
      do i=1,maxatom
         i3=3*(i-1)
         pos(i3+1) = i
         pos(i3+2) = i
         pos(i3+3) = 1
         
         charge(i) = mod(i,2)-0.5
      end do

      jjnr(1) = 13
      do i=2,maxlist
         jjnr(i) = mod(jjnr(i-1)+13,maxatom)
      end do      
         
      ix = 0.0
      iy = 0.0
      iz = 0.0
      qi = 1.0
      Vc = 0.0
      Vnb= 0.0
      
c     run it


      do j=1,100
         do i=1,maxatom
            
            call FORLJC(ix,iy,iz,qi,
     $           pos,maxlist,type,jjnr,charge,nbfp,
     $           faction,fip,
     $           Vc,Vnb)
            
         end do
      end do
      
      print *,Vc, Vnb
      
      stop
      end

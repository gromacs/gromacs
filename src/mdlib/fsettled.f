
C
C                This source code is part of
C 
C                 G   R   O   M   A   C   S
C 
C          GROningen MAchine for Chemical Simulations
C 
C                        VERSION 3.0
C 
C Copyright (c) 1991-2001
C BIOSON Research Institute, Dept. of Biophysical Chemistry
C University of Groningen, The Netherlands
C 
C This program is free software; you can redistribute it and/or
C modify it under the terms of the GNU General Public License
C as published by the Free Software Foundation; either version 2
C of the License, or (at your option) any later version.
C 
C If you want to redistribute modifications, please consider that
C scientific software is very special. Version control is crucial -
C bugs must be traceable. We will be happy to consider code for
C inclusion in the official distribution, but derived work must not
C be called official GROMACS. Details are found in the README & COPYING
C files - if they are missing, get the official version at www.gromacs.org.
C 
C To help us fund GROMACS development, we humbly ask that you cite
C the papers on the package - you can find them in the top README file.
C 
C Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
C 
C And Hey:
C GROup of MAchos and Cynical Suckers

c     IMPORTANT IMPORTANT IMPORTANT IMPORTANT !
c     Note that this file comes in two flavours -
c     fshake.f for single precision and fshaked.f 
c     for double precision. The only difference is 
c     the size of the real variables.
c     This is an unfortunate, but necessary setup      
c     since not all f77 compilers (e.g. g77) have
c     switches to change the real size, and neither
c     do all f77 compilers support preprocessing.
c     Thus, if you edit one of the files, make sure  
c      to change to other similarly!
      
      subroutine fsettle(nshake,owptr,b4,after,
     &     dOH,dHH,mO,mH,error)
      implicit none
c*****************************************************************
c                                                               **
c    Subroutine : setlep - reset positions of TIP3P waters      **
c    Author : Shuichi Miyamoto                                  **
c    Date of last update : Oct. 1, 1992                         **
c                                                               **
c    Reference for the SETTLE algorithm                         **
c           S. Miyamoto et al., J. Comp. Chem., 13, 952 (1992). **
c                                                               **
c*****************************************************************
      integer nshake,owptr(*),error

      real*8  b4(*),after(*),mO,mH,dOH,dHH
      real*8  wo,wh,wohh,ra,rb,rc,rc2,tmp,tmp2
      real*8 gama, beta, alpa, xcom, ycom, zcom, al2be2
      real*8 axlng, aylng, azlng, trns11, trns21, trns31, trns12, 
     &     trns22, trns32, trns13, trns23, trns33, cosphi, 
     &     costhe, sinphi, sinthe, cospsi, xaksxd, yaksxd, 
     &     xakszd, yakszd, zakszd, zaksxd, xaksyd, 
     &     xb0, yb0, zb0, xc0, yc0, zc0, xa1
      real*8 ya1, za1, xb1, yb1
      real*8 zb1, xc1, yc1, zc1, yaksyd, zaksyd, sinpsi, 
     &     xa3, ya3, za3, xb3, yb3, zb3, xc3, yc3, zc3, 
     &     xb0d, yb0d, xc0d, yc0d, xa1d, ya1d, 
     &     za1d, xb1d, yb1d, zb1d, xc1d, yc1d, zc1d, ya2d, 
     &     xb2d, yb2d, yc2d, 
     &     xa3d, ya3d, za3d, xb3d, yb3d, zb3d, xc3d, yc3d, zc3d
      integer i,ow1,hw2,hw3

      error= -1
      wo   = mO
      wh   = mH
      wohh = mO+2.0*mH
      rc   = dHH/2.0
      ra   = 2.0*wh*dsqrt(dOH*dOH-rc*rc)/wohh
      rb   = dsqrt(dOH*dOH-rc*rc)-ra
      rc2  = dHH
      wo   = wo/wohh
      wh   = wh/wohh
c
cray compiler directive ignore vector dependencies      
c$dir ivdep
      do i=1,nshake
c                                                --- Step1  A1' ---
c
         ow1 = 3*owptr(i)+1
         hw2 = ow1+3
         hw3 = ow1+6
         xb0 = b4(hw2)   - b4(ow1)
         yb0 = b4(hw2+1) - b4(ow1+1)
         zb0 = b4(hw2+2) - b4(ow1+2)
         xc0 = b4(hw3)   - b4(ow1)
         yc0 = b4(hw3+1) - b4(ow1+1)
         zc0 = b4(hw3+2) - b4(ow1+2)
c
         xcom = ( after(ow1  )*wo + 
     &        (after(hw2  ) + after(hw3  )) * wh )
         ycom = ( after(ow1+1)*wo + 
     &        (after(hw2+1) + after(hw3+1)) * wh )
         zcom = ( after(ow1+2)*wo +
     &        (after(hw2+2) + after(hw3+2)) * wh )
c
         xa1 = after(ow1  ) - xcom
         ya1 = after(ow1+1) - ycom
         za1 = after(ow1+2) - zcom
         xb1 = after(hw2  ) - xcom
         yb1 = after(hw2+1) - ycom
         zb1 = after(hw2+2) - zcom
         xc1 = after(hw3  ) - xcom
         yc1 = after(hw3+1) - ycom
         zc1 = after(hw3+2) - zcom
c
         xaksZd = yb0*zc0 - zb0*yc0
         yaksZd = zb0*xc0 - xb0*zc0
         zaksZd = xb0*yc0 - yb0*xc0
         xaksXd = ya1*zaksZd - za1*yaksZd
         yaksXd = za1*xaksZd - xa1*zaksZd
         zaksXd = xa1*yaksZd - ya1*xaksZd
         xaksYd = yaksZd*zaksXd - zaksZd*yaksXd
         yaksYd = zaksZd*xaksXd - xaksZd*zaksXd
         zaksYd = xaksZd*yaksXd - yaksZd*xaksXd
c
         axlng = 1.0/dsqrt ( xaksXd * xaksXd + yaksXd * yaksXd
     &                                  + zaksXd * zaksXd )
         aylng = 1.0/dsqrt ( xaksYd * xaksYd + yaksYd * yaksYd
     &                                  + zaksYd * zaksYd )
         azlng = 1.0/dsqrt ( xaksZd * xaksZd + yaksZd * yaksZd
     &                                  + zaksZd * zaksZd )
         trns11 = xaksXd * axlng
         trns21 = yaksXd * axlng
         trns31 = zaksXd * axlng
         trns12 = xaksYd * aylng
         trns22 = yaksYd * aylng
         trns32 = zaksYd * aylng
         trns13 = xaksZd * azlng
         trns23 = yaksZd * azlng
         trns33 = zaksZd * azlng
c 
         xb0d = trns11*xb0 + trns21*yb0 + trns31*zb0
         yb0d = trns12*xb0 + trns22*yb0 + trns32*zb0
         xc0d = trns11*xc0 + trns21*yc0 + trns31*zc0
         yc0d = trns12*xc0 + trns22*yc0 + trns32*zc0
         xa1d = trns11*xa1 + trns21*ya1 + trns31*za1
         ya1d = trns12*xa1 + trns22*ya1 + trns32*za1
         za1d = trns13*xa1 + trns23*ya1 + trns33*za1
         xb1d = trns11*xb1 + trns21*yb1 + trns31*zb1
         yb1d = trns12*xb1 + trns22*yb1 + trns32*zb1
         zb1d = trns13*xb1 + trns23*yb1 + trns33*zb1
         xc1d = trns11*xc1 + trns21*yc1 + trns31*zc1
         yc1d = trns12*xc1 + trns22*yc1 + trns32*zc1
         zc1d = trns13*xc1 + trns23*yc1 + trns33*zc1
c
         sinphi = za1d / ra
         tmp    = 1.0 - sinphi*sinphi
         if ( tmp .le. 0.0) then
            error = i-1
            cosphi = 0
         else
            cosphi = dsqrt (tmp)
         endif
         sinpsi = ( zb1d - zc1d ) / (rc2 * cosphi)
         tmp2   = 1.0 - sinpsi*sinpsi
         if ( tmp2 .le. 0.0 ) then
            error = i-1
            cospsi = 0
         else
            cospsi = dsqrt (tmp2)
         endif
c 
         ya2d =   ra * cosphi
         xb2d = - rc * cospsi
c        xc2d =   rc * cospsi
         yb2d = - rb * cosphi - rc *sinpsi * sinphi
         yc2d = - rb * cosphi + rc *sinpsi * sinphi
c        xb2d2 = xb2d * xb2d
c        hh2 = 4.d0 * xb2d2 + (yb2d-yc2d) * (yb2d-yc2d)
c    &                      + (zb1d-zc1d) * (zb1d-zc1d)
c        deltx = 2.d0 * xb2d + dsqrt ( 4.d0 * xb2d2 - hh2 + hhhh )
c        xb2d = xb2d - deltx * 0.5d0
c        xc2d = xc2d + deltx * 0.5d0
c
c                                                --- Step3  al,be,ga ---
c
         alpa = ( xb2d * (xb0d-xc0d) + yb0d * yb2d + yc0d * yc2d )
         beta = ( xb2d * (yc0d-yb0d) + xb0d * yb2d + xc0d * yc2d )
         gama = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d
c
         al2be2 = alpa * alpa + beta * beta
         sinthe = ( alpa*gama - beta * dsqrt ( al2be2 - gama * gama ) )
     &            / al2be2
c
c                                                --- Step4  A3' ---
c
         costhe = dsqrt (1.d0 - sinthe * sinthe )
         xa3d = - ya2d * sinthe
         ya3d =   ya2d * costhe
         za3d = za1d
         xb3d =   xb2d * costhe - yb2d * sinthe
         yb3d =   xb2d * sinthe + yb2d * costhe
         zb3d = zb1d
         xc3d = - xb2d * costhe - yc2d * sinthe
         yc3d = - xb2d * sinthe + yc2d * costhe
         zc3d = zc1d
c                                                --- Step5  A3 ---
c
c
         xa3 = trns11*xa3d + trns12*ya3d + trns13*za3d
         ya3 = trns21*xa3d + trns22*ya3d + trns23*za3d
         za3 = trns31*xa3d + trns32*ya3d + trns33*za3d
         xb3 = trns11*xb3d + trns12*yb3d + trns13*zb3d
         yb3 = trns21*xb3d + trns22*yb3d + trns23*zb3d
         zb3 = trns31*xb3d + trns32*yb3d + trns33*zb3d
         xc3 = trns11*xc3d + trns12*yc3d + trns13*zc3d
         yc3 = trns21*xc3d + trns22*yc3d + trns23*zc3d
         zc3 = trns31*xc3d + trns32*yc3d + trns33*zc3d
c
         after(ow1  ) = xcom + xa3
         after(ow1+1) = ycom + ya3
         after(ow1+2) = zcom + za3
         after(hw2  ) = xcom + xb3
         after(hw2+1) = ycom + yb3
         after(hw2+2) = zcom + zb3
         after(hw3  ) = xcom + xc3
         after(hw3+1) = ycom + yc3
         after(hw3+2) = zcom + zc3
c
      end do
      return
      end


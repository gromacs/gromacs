      SUBROUTINE fournf(data,nn,ndim,isign)
      
      implicit none
      
      INTEGER isign,ndim,nn(ndim)
      REAL data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
     &     k2,n,nprev,nrem,ntot
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
 11   continue
      nprev=1
      do 18 idim=1,ndim
         n=nn(idim)
         nrem=ntot/(n*nprev)
         ip1=2*nprev
         ip2=ip1*n
         ip3=ip2*nrem
         i2rev=1
         do 14 i2=1,ip2,ip1
            if(i2.lt.i2rev)then
               do 13 i1=i2,i2+ip1-2,2
                  do 12 i3=i1,ip3,ip2
                     i3rev=i2rev+i3-i2
                     tempr=data(i3)
                     tempi=data(i3+1)
                     data(i3)=data(i3rev)
                     data(i3+1)=data(i3rev+1)
                     data(i3rev)=tempr
                     data(i3rev+1)=tempi
 12               continue
 13            continue
            endif
       ibit=ip2/2
 1     if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
          i2rev=i2rev-ibit
          ibit=ibit/2
          goto 1
       endif
       i2rev=i2rev+ibit
 14   continue
      ifp1=ip1
 2    if(ifp1.lt.ip2)then
         ifp2=2*ifp1
         theta=isign*6.28318530717959d0/(ifp2/ip1)
         wpr=-2.d0*sin(0.5d0*theta)**2
         wpi=sin(theta)
         wr=1.d0
         wi=0.d0
          do 17 i3=1,ifp1,ip1
             do 16 i1=i3,i3+ip1-2,2
                do 15 i2=i1,ip3,ifp2
                   k1=i2
                   k2=k1+ifp1
                   tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                   tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                   data(k2)=data(k1)-tempr
                   data(k2+1)=data(k1+1)-tempi
                   data(k1)=data(k1)+tempr
                   data(k1+1)=data(k1+1)+tempi
 15             continue
 16          continue
             wtemp=wr
             wr=wr*wpr-wi*wpi+wr
             wi=wi*wpr+wtemp*wpi+wi
 17       continue
          ifp1=ifp2
          goto 2
       endif
       nprev=n*nprev
 18   continue
      return
      END
C     (C) Copr. 1986-92 Numerical Recipes Software +).
      

      SUBROUTINE rlft3f(data,speq,nn1,nn2,nn3,isign)
      
      implicit none
      
      INTEGER isign,nn1,nn2,nn3
      COMPLEX data(nn1/2,nn2,nn3),speq(nn2,nn3)
      
CU    USES fourn

      INTEGER i1,i2,i3,j1,j2,j3,nn(3)
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      COMPLEX c1,c2,h1,h2,w
      
      c1    = cmplx(0.5,0.0)
      c2    = cmplx(0.0,-0.5*isign)
      theta = 6.28318530717959d0/dble(isign*nn1)
      wpr   = -2.0d0*sin(0.5d0*theta)**2
      wpi   = sin(theta)
      nn(1) = nn1/2
      nn(2) = nn2
      nn(3) = nn3
      
      if (isign.eq.1) then
         call fournf(data,nn,3,isign)
         do  i3=1,nn3
            do  i2=1,nn2
               speq(i2,i3)=data(1,i2,i3)
            end do
         end do
      endif
      do i3=1,nn3
         j3=1
         if (i3.ne.1) j3=nn3-i3+2
         wr=1.0d0
         wi=0.0d0
         do i1=1,nn1/4+1
            j1=nn1/2-i1+2
            do i2=1,nn2
               j2=1
               if (i2.ne.1) j2=nn2-i2+2
               if (i1.eq.1) then
                  h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
                  h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
                  data(1,i2,i3)=h1+h2
                  speq(j2,j3)=conjg(h1-h2)
               else
                  h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
                  h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
                  data(i1,i2,i3)=h1+w*h2
                  data(j1,j2,j3)=conjg(h1-w*h2)
               endif
            end do
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
            w=cmplx(sngl(wr),sngl(wi))
         end do
      end do
      if (isign.eq.-1) then
         call fournf(data,nn,3,isign)
      endif
      return
      END
C     (C) Copr. 1986-92 Numerical Recipes Software +).

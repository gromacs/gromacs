c     IMPORTANT!
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
 
      subroutine flincsp(x,f,fp,ncons,
     $     bla1,bla2,blnr,blbnb,blc,blcc,blm,
     $     nrec,invmass,r,rhs1,rhs2,sol)

      implicit none

      real*8  x(*),f(*),fp(*),blc(*),blcc(*),blm(*),invmass(*)
      real*8  r(*),rhs1(*),rhs2(*),sol(*)
      integer*4 ncons,nrec,bla1(*),bla2(*),blnr(*),blbnb(*)

      integer*4 b,i,j,k,n,b3,i3,j3,it,rec
      real*8  tmp0,tmp1,tmp2,im1,im2,mvb,rlen 
      real*8  u0,u1,u2,v0,v1,v2


      do b=1,ncons
         b3=3*(b-1)
         i=3*bla1(b)
         j=3*bla2(b)
         tmp0=x(i+1)-x(j+1)
         tmp1=x(i+2)-x(j+2)
         tmp2=x(i+3)-x(j+3)
         rlen=1.0/sqrt(tmp0*tmp0+tmp1*tmp1+tmp2*tmp2)
         r(b3+1)=rlen*tmp0
         r(b3+2)=rlen*tmp1
         r(b3+3)=rlen*tmp2
      enddo

      do b=1,ncons
         b3=3*(b-1)
         tmp0=r(b3+1)
         tmp1=r(b3+2)
         tmp2=r(b3+3)
         i=3*bla1(b)
         j=3*bla2(b)
         do n=blnr(b)+1,blnr(b+1)
            k=3*blbnb(n)
            blm(n)=blcc(n)*(tmp0*r(k+1)+tmp1*r(k+2)+tmp2*r(k+3))
         enddo
         mvb=blc(b)*(tmp0*(f(i+1)-f(j+1))+
     $               tmp1*(f(i+2)-f(j+2))+    
     $               tmp2*(f(i+3)-f(j+3)))
         rhs1(b)=mvb
         sol(b)=mvb
      enddo

      
      do rec=1,nrec,2
         do b=1,ncons
            mvb=0.
            do n=blnr(b)+1,blnr(b+1)
               j=blbnb(n)+1
               mvb=mvb+blm(n)*rhs1(j)
            enddo
            rhs2(b)=mvb
            sol(b)=sol(b)+mvb
         enddo 
         if (rec .lt. nrec) then
            do b=1,ncons
               mvb=0.
               do n=blnr(b)+1,blnr(b+1)
                  j=blbnb(n)+1
                  mvb=mvb+blm(n)*rhs2(j)
               enddo
               rhs1(b)=mvb
               sol(b)=sol(b)+mvb
            enddo 
            endif
         enddo
         
         
      do b=1,ncons
         b3=3*(b-1)
         i=bla1(b)
         j=bla2(b)
         i3=3*i
         j3=3*j
         mvb=blc(b)*sol(b)
         im1=invmass(i+1)
         im2=invmass(j+1)
         tmp0=r(b3+1)*mvb
         tmp1=r(b3+2)*mvb
         tmp2=r(b3+3)*mvb
         u0=fp(i3+1)-tmp0*im1
         u1=fp(i3+2)-tmp1*im1
         u2=fp(i3+3)-tmp2*im1
         v0=fp(j3+1)+tmp0*im2
         v1=fp(j3+2)+tmp1*im2
         v2=fp(j3+3)+tmp2*im2
         fp(i3+1)=u0
         fp(i3+2)=u1
         fp(i3+3)=u2
         fp(j3+1)=v0
         fp(j3+2)=v1
         fp(j3+3)=v2
      enddo

      end


      subroutine flincs(x,xp,ncons,
     $     bla1,bla2,blnr,blbnb,bllen,blc,blcc,blm,
     $     nit,nrec,invmass,r,rhs1,rhs2,sol,wangle,warn,
     $     lambda)

      implicit none

      real*4  x(*),xp(*),bllen(*),blc(*),blcc(*),blm(*),invmass(*)
      real*4  r(*),rhs1(*),rhs2(*),sol(*),wangle,lambda(*)
      real*4  tmp0,tmp1,tmp2,im1,im2,mvb,rlen,len,wfac,lam 
      real*4  u0,u1,u2,v0,v1,v2

      integer*4 ncons,nit,nrec,bla1(*),bla2(*),blnr(*),blbnb(*)
      integer*4 warn
      integer*4 b,i,j,k,n,b3,i3,j3,it,rec


      warn=0

      do b=1,ncons
         b3=3*(b-1)
         i=3*bla1(b)
         j=3*bla2(b)
         tmp0=x(i+1)-x(j+1)
         tmp1=x(i+2)-x(j+2)
         tmp2=x(i+3)-x(j+3)
         rlen=1.0/sqrt(tmp0*tmp0+tmp1*tmp1+tmp2*tmp2)
         r(b3+1)=rlen*tmp0
         r(b3+2)=rlen*tmp1
         r(b3+3)=rlen*tmp2
      enddo

      do b=1,ncons
         b3=3*(b-1)
         tmp0=r(b3+1)
         tmp1=r(b3+2)
         tmp2=r(b3+3)
         len=bllen(b)
         i=3*bla1(b)
         j=3*bla2(b)
         do n=blnr(b)+1,blnr(b+1)
            k=3*blbnb(n)
            blm(n)=blcc(n)*(tmp0*r(k+1)+tmp1*r(k+2)+tmp2*r(k+3))
         enddo
         mvb=blc(b)*(tmp0*(xp(i+1)-xp(j+1))+
     $               tmp1*(xp(i+2)-xp(j+2))+    
     $               tmp2*(xp(i+3)-xp(j+3))-len)
         rhs1(b)=mvb
         sol(b)=mvb
      enddo

      
      do rec=1,nrec,2
         do b=1,ncons
            mvb=0.
            do n=blnr(b)+1,blnr(b+1)
               j=blbnb(n)+1
               mvb=mvb+blm(n)*rhs1(j)
            enddo
            rhs2(b)=mvb
            sol(b)=sol(b)+mvb
         enddo 
         if (rec .lt. nrec) then
            do b=1,ncons
               mvb=0.
               do n=blnr(b)+1,blnr(b+1)
                  j=blbnb(n)+1
                  mvb=mvb+blm(n)*rhs2(j)
               enddo
               rhs1(b)=mvb
               sol(b)=sol(b)+mvb
            enddo 
            endif
         enddo
         
         
      do b=1,ncons
         b3=3*(b-1)
         i=bla1(b)
         j=bla2(b)
         i3=3*i
         j3=3*j
         mvb=blc(b)*sol(b)
         lambda(b)=-mvb
         im1=invmass(i+1)
         im2=invmass(j+1)
         tmp0=r(b3+1)*mvb
         tmp1=r(b3+2)*mvb
         tmp2=r(b3+3)*mvb
         u0=xp(i3+1)-tmp0*im1
         u1=xp(i3+2)-tmp1*im1
         u2=xp(i3+3)-tmp2*im1
         v0=xp(j3+1)+tmp0*im2
         v1=xp(j3+2)+tmp1*im2
         v2=xp(j3+3)+tmp2*im2
         xp(i3+1)=u0
         xp(i3+2)=u1
         xp(i3+3)=u2
         xp(j3+1)=v0
         xp(j3+2)=v1
         xp(j3+3)=v2
      enddo



c     ********  Correction for centripetal effects  ********

      wfac=cos(0.01745*wangle)
      wfac=wfac*wfac

      do it=1,nit

         do b=1,ncons
            len=bllen(b)
            i=3*bla1(b)
            j=3*bla2(b)
            tmp0=xp(i+1)-xp(j+1)
            tmp1=xp(i+2)-xp(j+2)
            tmp2=xp(i+3)-xp(j+3)
            u1=len*len
            u0=2.*u1-(tmp0*tmp0+tmp1*tmp1+tmp2*tmp2)
            if (u0 .lt. wfac*u1) warn=b  
            if (u0 .lt. 0.) u0=0.
            mvb=blc(b)*(len-sqrt(u0))
            rhs1(b)=mvb
            sol(b)=mvb
         enddo


         do rec=1,nrec,2
            do b=1,ncons
               mvb=0.
               do n=blnr(b)+1,blnr(b+1)
                  j=blbnb(n)+1
                  mvb=mvb+blm(n)*rhs1(j)
               enddo
               rhs2(b)=mvb
               sol(b)=sol(b)+mvb
            enddo 
            if (rec .lt. nrec) then
               do b=1,ncons
                  mvb=0.
                  do n=blnr(b)+1,blnr(b+1)
                     j=blbnb(n)+1
                     mvb=mvb+blm(n)*rhs2(j)
                  enddo
                  rhs1(b)=mvb
                  sol(b)=sol(b)+mvb
               enddo 
            endif
         enddo
         
         
         do b=1,ncons
            b3=3*(b-1)
            i=bla1(b)
            j=bla2(b)
            i3=3*i
            j3=3*j
            lam=lambda(b)
            mvb=blc(b)*sol(b)
            lambda(b)=lam-mvb
            im1=invmass(i+1)
            im2=invmass(j+1)
            tmp0=r(b3+1)*mvb
            tmp1=r(b3+2)*mvb
            tmp2=r(b3+3)*mvb
            u0=xp(i3+1)-tmp0*im1
            u1=xp(i3+2)-tmp1*im1
            u2=xp(i3+3)-tmp2*im1
            v0=xp(j3+1)+tmp0*im2
            v1=xp(j3+2)+tmp1*im2
            v2=xp(j3+3)+tmp2*im2
            xp(i3+1)=u0
            xp(i3+2)=u1
            xp(i3+3)=u2
            xp(j3+1)=v0
            xp(j3+2)=v1
            xp(j3+3)=v2
         enddo

      enddo

      return

      end

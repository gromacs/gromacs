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
 
      subroutine fshake(iatom,ncon,nit,maxnit,
     $     dist2,xp,rij,m2,invmass,tt,lagr,error)
      
      implicit none
      
      integer*4 iatom(*)
      integer*4 ncon,nit,maxnit
      integer*4 error
      real*4    dist2(*),xp(*),rij(*),m2(*),invmass(*),tt(*),lagr(*)
      
      integer*4 ll,i,j,i3,j3,l3,nconv,iconv
      integer*4 ix,iy,iz,jx,jy,jz
      real*4  toler,rpij2,rrpr,tx,ty,tz,diff,acor,im,jm
      real*4  xh,yh,zh,rijx,rijy,rijz
      real*4  tix,tiy,tiz
      real*4  tjx,tjy,tjz
      real*4  mytol
      
      parameter(mytol=1e-6)
      
      error=0
      do nit=1,maxnit
         nconv=0
         do ll=1,ncon
            l3    = 3*(ll-1)
            rijx  = rij(l3+1)
            rijy  = rij(l3+2)
            rijz  = rij(l3+3)
            i     = iatom(l3+2)
            j     = iatom(l3+3)
            i3    = 3*i
            j3    = 3*j
            ix    = i3+1
            iy    = i3+2
            iz    = i3+3
            jx    = j3+1
            jy    = j3+2
            jz    = j3+3
            
            tx      = xp(ix)-xp(jx)
            ty      = xp(iy)-xp(jy)
            tz      = xp(iz)-xp(jz)
            rpij2   = tx*tx+ty*ty+tz*tz
            
            toler   = dist2(ll)
            diff    = toler-rpij2
            iconv   = abs(diff)*tt(ll)

            if (iconv .ne. 0) then
               nconv   = nconv + iconv
               rrpr    = rijx*tx+rijy*ty+rijz*tz
            
               if (rrpr .lt. mytol*toler) then
                  error   = ll
               else
                  acor    = diff*m2(ll)/rrpr
                  lagr(ll) = lagr(ll) + acor
                  xh      = rijx*acor
                  yh      = rijy*acor
                  zh      = rijz*acor
                  im      = invmass(i)
                  jm      = invmass(j)
                  if ((im .ne. 0) .and. (jm .ne. 0)) then
                     xp(ix) = xp(ix) + xh*im
                     xp(iy) = xp(iy) + yh*im
                     xp(iz) = xp(iz) + zh*im
                     xp(jx) = xp(jx) - xh*jm
                     xp(jy) = xp(jy) - yh*jm
                     xp(jz) = xp(jz) - zh*jm
                  else
                     if ((im .eq. 0) .and. (jm .ne. 0)) then
                        xp(ix) = xp(ix) + xh*jm
                        xp(iy) = xp(iy) + yh*jm
                        xp(iz) = xp(iz) + zh*jm
                     else if ((jm .eq. 0) .and. (im .ne. 0)) then
                        xp(jx) = xp(jx) - xh*im
                        xp(jy) = xp(jy) - yh*im
                        xp(jz) = xp(jz) - zh*im
                     else
                        error = ll
                     endif
                  endif
               end if
            end if
         end do
         
         if (nconv .eq. 0) goto 10
         if (error .ne. 0) goto 10
      end do
      
 10   return
 
      end
      

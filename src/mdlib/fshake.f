      subroutine FORSHAKE(iatom,ncon,nit,maxnit,
     $     dist2,xp,rij,m2,invmass,tt,error)
      
      implicit none
      
      integer*4 iatom(*)
      integer*4 ncon,nit,maxnit
      integer*4 error
      real      dist2(*),xp(*),rij(*),m2(*),invmass(*),tt(*)
      
      integer*4 ll,i,j,i3,j3,l3,nconv,iconv
      integer*4 ix,iy,iz,jx,jy,jz
      real    toler,rpij2,rrpr,tx,ty,tz,diff,acor,im,jm
      real    xh,yh,zh,rijx,rijy,rijz
      real    tix,tiy,tiz
      real    tjx,tjy,tjz
      real    mytol
      
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
                  im      = invmass(i+1)
                  jm      = invmass(j+1)
                  xh      = rijx*acor
                  yh      = rijy*acor
                  zh      = rijz*acor
                  tix     = xp(ix) + xh*im
                  tiy     = xp(iy) + yh*im
                  tiz     = xp(iz) + zh*im
                  tjx     = xp(jx) - xh*jm
                  tjy     = xp(jy) - yh*jm
                  tjz     = xp(jz) - zh*jm
                  xp(ix) = tix
                  xp(iy) = tiy
                  xp(iz) = tiz
                  xp(jx) = tjx
                  xp(jy) = tjy
                  xp(jz) = tjz
               end if
            end if
         end do
         
         if (nconv .eq. 0) goto 10
         if (error .ne. 0) goto 10
      end do
      
 10   return
 
      end
      

      subroutine doconv(n,c,g)
      
      integer*4 n
      real      g(*),c(*)
      
      integer*4 i,i2
      real      gg

C     Do an optimized convolution
      do i=1,n
         i2 = 2*i
         gg = g(i)
         c(i2)   = gg*c(i2)
         c(i2+1) = gg*c(i2+1)
      end do
      
      return
      end
      

      program tstsqrtf
      
      integer  i
      real     x,y,z,dy,rdy,invsqrt
      
      call fillbuf

      print '(5(a12,2x))','X','Y','Z','DY','DY/Z'
      do i=1,1000
         x = 1.0*i
         y = invsqrt(x)
         z = 1.0/sqrt(x)
         dy = y-z
         rdy = dy/z
         print '(5(e12.5,2x))',x,y,z,dy,rdy
      end do
      
      stop
      
      end

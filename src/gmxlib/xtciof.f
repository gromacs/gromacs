      subroutine xtccheck(ret,ivar)
c
c     Passed variables
c
      integer ret,ivar
      
      if ( ret .eq. 0 ) then
         if ( ivar .eq. 1 ) print *,'XTC-Error reading/writing natoms'
         if ( ivar .eq. 2 ) print *,'XTC-Error reading/writing step'
         if ( ivar .eq. 3 ) print *,'XTC-Error reading/writing time'
         if ( ivar .eq. 4 ) print *,'XTC-Error reading/writing box'
         if ( ivar .eq. 5 ) print *,'XTC-Error reading/writing x'
         stop
      endif
      
      return
      end
        
      subroutine xtcheader(xd,magic,natoms,step,time,ret)
c
c     Passed variables
c
      integer xd,natoms,step,ret
      real time
      
      call xdrfint(xd,magic,ret)
      if ( ret .eq. 0 ) then
         return
      endif
      
      call xdrfint(xd,natoms,ret)
      call xtccheck(ret,1)
      
      call xdrfint(xd,step,ret)
      call xtccheck(ret,2)
      
      call xdrffloat(xd,time,ret)
      call xtccheck(ret,3)

      return
      end
      
      subroutine xtccoord(xd,natoms,box,x,prec,ret)
c
c     Passed variables
c
      integer xd,natoms,ret
      real box(9),x(*),prec
c     
c     local
c
      integer i
      
      do i=1,9
         call xdrffloat(xd,box(i),ret)
         call xtccheck(ret,4)
      end do
  
      call xdrf3dfcoord(xd,x,natoms,prec,ret)
      call xtccheck(ret,5)
      
      return
      end
              
      subroutine xtcio(xd,natoms,step,time,box,x,prec,mode,ret)
c
c     Passed variables
c
      integer xd,natoms,step,mode,ret
      real box(9)
      real x(*)
      real time,prec
c
c     local variables
c
      integer xtcmagic
      integer magic
      
      xtcmagic=1995
      
      if (mode .eq. 0) magic=xtcmagic
      
      call xtcheader(xd,magic,natoms,step,time,ret)

      if (ret .eq. 0) return
            
      if (mode .eq. 1) then
         if (magic .ne. xtcmagic) then
            print *,'Fatal error, magic number read as ',magic,
     &           ' should be ',xtcmagic
         endif
      endif

      call xtccoord(xd,natoms,box,x,prec,ret)
      
      return 
      end

      subroutine readxtc(xd,natoms,step,time,box,x,prec,ret)

      integer xd,natoms,step,mode,ret
      real    box(9),time,prec
      
      call xtcio(xd,natoms,step,time,box,x,prec,1,ret)
      
      return 
      end

      subroutine writextc(xd,natoms,step,time,box,x,prec,ret)

      integer xd,natoms,step,mode,ret
      real    box(9),time,prec
      
      call xtcio(xd,natoms,step,time,box,x,prec,0,ret)
      
      return 
      end

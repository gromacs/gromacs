      program testxtc
      
      parameter (maxatom=10000,maxx=3*maxatom)
      integer xd,xd2,natoms,step,ret,i
      real    time,box(9),x(maxx)
      
      call xdrfopen(xd,"test.xtc","r",ret)
      print *,'opened test.xtc, ret=',ret
      call xdrfopen(xd2,"testout.xtc","w",ret)
      print *,'opened testout.xtc, ret=',ret
      
 10   call readxtc(xd,natoms,step,time,box,x,prec,ret)
 
      if ( ret .eq. 1 ) then
         call writextc(xd2,natoms,step,time,box,x,prec,ret)
         goto 10
      endif
         
      stop
      end

      program conveig
      
      parameter (maxatom=10000,maxdim=3*maxatom)
      
      implicit none
      
      integer      dim,i,j,j3,nvec,natom
      real*8       ev(maxdim),evnorm
      character*80 fn
      
      print *,'Number of ED atoms ?'
      read *,natom
      
      if (natom .gt. maxatom) then
         print *,'Too big!'
         stop
      endif
      
      dim=3*natom
      
      print *,'Number of eigenvectors ?'
      read *,nvec
      
      print *,'Using ',nvec,' dimensions for vectors'

      call getarg(1,fn)
            
      open (unit=10,file=fn,form='formatted')
      open (unit=11,file='eigenvec.dat',form='formatted')
      do i=1,nvec
         read (10,'(6e12.5)') (ev(j),j=1,dim)
         
         do j=1,natom
            j3=3*(j-1)+1
            evnorm=sqrt(ev(j3)**2+ev(j3+1)**2+ev(j3+2)**2)
            write (11,'(i5,2x,i5,2x,e12.5,2x,e12.5,2x,e12.5,2x,e12.5)') 
     &           i,j,evnorm,ev(j3),ev(j3+1),ev(j3+2)
         end do
         
      end do
      
      stop
      end
      

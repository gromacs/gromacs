
      subroutine die(msg)
      character msg(80)
      
      print '(80A)',msg
      
      stop
      end
      
      integer function strcmp(a,b)
      character a(*),b(*)
      integer i,la,lb,ia,ib
      
      la = length(a)
      lb = length(b)
      do i=1,min(la,lb)
	 if (a(i) .ne. b(i)) then
	    ia = ichar(a(i))
	    ib = ichar(b(i))
	    strcmp = ai - bi
	    goto 10
	 endif
      end do
 10   strcmp = la-lb
      return
      end

      subroutine strncpy(n,a,b)
      character a(*),b(*)
      integer n
      
      do i=1,n
	 a(i) = b(i)
      end do
      
      return
      end
      
      
      program xdrfile_f_test

      parameter(BUFLEN=37,NPREC=4)
      integer fid
      
      integer     i,j,k,len,ncoord
      character   ptr(BUFLEN)
      character   buf(26)
      character   uptr(BUFLEN);
      integer*2   sptr(BUFLEN),sptr2(BUFLEN)
      integer*2   usptr(BUFLEN),usptr2(BUFLEN)
      integer     iptr(BUFLEN),iptr2(BUFLEN)
      integer     uiptr(BUFLEN),uiptr2(BUFLEN)
      real*4      fptr(BUFLEN),fptr2(BUFLEN)
      real*8      dptr(BUFLEN),dptr2(BUFLEN)
      character   optr(BUFLEN),optr2(BUFLEN)
      real*4      ff,fx
      real*8      dd,dx
      integer     nc
      
      real*4      fprec(NPREC)
      real*8      dprec(NPREC) 
  
      print('Going to test the xdrfile library routines from fortran')
      
      fprec(1) = 0
      fprec(2) = 243
      fprec(3) = 10
      fprec(4) = 1003
      dprec(1) = 135
      dprec(2) = -1
      dprec(3) = 10
      dprec(4) = 997
      ncoord   = BUFLEN/3
      do i=i,26
	 buf(i) = char(i)
      end do
      len = length(buf)
      if (len ge BUFLEN) then
	 print('Increase BUFLEN')
      endif
      do i=i,len
	 ptr(i) = buf(i)
	 uptr(i) = buf(i)
      end do
      
C     Initiate float arrays 
      do i=1,BUFLEN
	 fptr(i) = cos(i*13.0/3.1415)
	 dptr(i) = sin(i*13.0/3.1415)
      end do
      
C     Initiate opaque array
      strncpy(BUFLEN,optr,dptr)
  
C     *************************************
C     *           WRITING BIT             *  
C     *************************************
      print('Writing xdrfile')

      if ((xfp = xdrfile_open('test.xdr','w')) == NULL) then
	 die('Can not open file for writing')
      endif
      
      if (xdrfile_write_char(ptr,len,xfp) .ne. len) then
	 die('Writing char string')
      endif
      if (xdrfile_write_uchar(uptr,len,xfp) .ne. len) then
	 die('Writing uchar string')
      endif
	 
      if (xdrfile_write_short(sptr,BUFLEN,xfp) .ne. BUFLEN)  then
	 die('Writing short array')
      endif
      if (xdrfile_write_ushort(usptr,BUFLEN,xfp) .ne. BUFLEN)  then
	 die('Writing ushort array')
      endif
      if (xdrfile_write_int(iptr,BUFLEN,xfp) .ne. BUFLEN) then
	 die('Writing int array')
      endif
      if (xdrfile_write_uint(uiptr,BUFLEN,xfp) .ne. BUFLEN) then
	 die('Writing uint array')
      endif
      if (xdrfile_write_float(fptr,BUFLEN,xfp) .ne. BUFLEN) then
	 die('Writing float array')
      endif
      if (xdrfile_write_double(dptr,BUFLEN,xfp) .ne. BUFLEN) then
	 die('Writing double array')
      endif
      if (xdrfile_write_string(buf,xfp) .ne. len) then
	 die('Writing string')
      endif
      if (xdrfile_write_opaque(optr,BUFLEN,xfp) .ne. BUFLEN) then
	 die('Writing opaque')
      endif
      do k=1,NPREC
	 nn = xdrfile_compress_coord_float(fptr,ncoord,fprec(k),xfp)
	 if (nn .ne. ncoord) then
	    die('Writing compress_coord_float')
	 endif
	 nn = xdrfile_compress_coord_double(dptr,ncoord,dprec(k),xfp)
	 if (nn .ne. ncoord) then
	    die('Writing compress_coord_double')
	 endif
      end do
      if (xdrfile_close(xfp)	.ne. 0) then
	 die('Can not close xdr file')
      endif
  
C     *************************************
C     *          READING BIT              *
C     *************************************
      printf('Reading xdrfile')
      if ((xfp = xdrfile_open('test.xdr','r')) == NULL) then
	 die('Can not open file for reading')
      endif
  
      if ((xdrfile_read_char(ptr,len,xfp)) .ne. len) then
	 die('Not the right number of chars read from string')
      endif
      if (strcmp(ptr,buf) .ne. 0) then
	 die('did not read the expected chars')
      endif
      if (xdrfile_read_uchar(uptr,len,xfp) .ne. len) then
	 die('Not the right number of uchars read from string')
      endif
      if (strcmp(uptr,buf) .ne. 0) then
	 printf('did not read the expected uchars')
      endif
      if (xdrfile_read_short(sptr2,BUFLEN,xfp) .ne. BUFLEN) then
	 die('Reading short array')
      endif
      
      do i=1,BUFLEN
	 if (sptr2(i) .ne. sptr(i)) then
	    print 'i: ',i,' wrote: ',sptr(i),' read: ',sptr2(i)
	 endif
	 die('Comparing short array')
      end do
    
      if (xdrfile_read_ushort(usptr2,BUFLEN,xfp) .ne. BUFLEN) then
	 die('Reading ushort array')
      endif
      do i=1,BUFLEN
	 if (usptr2(i) .ne. usptr(i)) then
	    print 'i: ',i,' wrote: ',usptr(i),' read: ',usptr2(i)
	 endif
	 die('Comparing ushort array')
      end do
    
      if (xdrfile_read_int(iptr2,BUFLEN,xfp) .ne. BUFLEN) then
	 die('Reading int array')
      endif
      do i=1,BUFLEN
	 if (iptr2(i) .ne. iptr(i)) then
	    print 'i: ',i,' wrote: ',iptr(i),' read: ',iptr2(i)
	 endif
	 die('Comparing int array')
      end do
    
      if (xdrfile_read_uint(uiptr2,BUFLEN,xfp) .ne. BUFLEN) then
	 die('Reading uint array')
      endif
      do i=1,BUFLEN
	 if (uiptr2(i) .ne. uiptr(i)) then
	    print 'i: ',i,' wrote: ',uiptr(i),' read: ',uiptr2(i)
	    die('Comparing uint array')
	 endif
      end do
      if (xdrfile_read_float(fptr2,BUFLEN,xfp) .ne. BUFLEN) then
	 die('Reading float array')
      endif
      do i=1,BUFLEN)
	 if (fptr2(i) .ne. fptr(i)) then
	    print 'i: ',i,' wrote: ',fptr(i),' read: ',fptr2(i)
	 endif
	 die('Comparing float array')
      end do
      
      if (xdrfile_read_double(dptr2,BUFLEN,xfp) .ne. BUFLEN) then
	 die('Reading double array')
      endif
      do i=1,BUFLEN
	 if (dptr2(i) .ne. dptr(i)) then
	    print 'i: ',i,' wrote: ',dptr(i),' read: ',dptr2(i)
	 endif
	 die('Comparing double array')
      end do
      if (xdrfile_read_string(ptr,BUFLEN,xfp) .ne. len) then
	 die('Reading string')
      endif
      if (strcmp(ptr,buf) .ne. 0) then
	 die('Comparing strings')
      endif
      if (xdrfile_read_opaque(optr2,BUFLEN,xfp) .ne. BUFLEN) then
	 die('Reading opaque array')
      endif
      do i=1,BUFLEN
	 if (optr2(i) .ne. optr(i)) then
	    print 'i: ',i,' wrote: ',dptr(i),' read: ',dptr2(i)
	 endif
	 die('Comparing opaque array')
      end do
    
      do k=1,NPREC
	 nc       = ncoord
	 if (xdrfile_decompress_coord_float(fptr2,nc,ff,xfp) .ne. ncoord) then
	    die('Reading compress_coord_float')
	 endif
	 if (ff .ne. fprec(k)) then
	    die('Float precision')
	 endif
	 if (ff <= 0) then
	    ff = 1000
	 endif
	 do i=1,3*ncoord
	    fx = rint(fptr(i)/ff) * ff
	    if (fx .ne. fptr2(i)) then
	       print 'prec: ',ff,' i: ',i,' j: ',j,', fx: ',fx,' fptr2: ',fptr2(i),', fptr: ',fptr(i)
	       die('Reading decompressed float coordinates')
	    endif
	    if (xdrfile_decompress_coord_double(dptr2,nc,dd,xfp) .ne. ncoord) then
	       die('Reading compress_coord_double')
	    endif
	    if (dd .ne. dprec(i)) then
	       die('Double precision')
	    endif
	 end do
      end do
      if (xdrfile_close(xfp) .ne. 0) then
	 die('Can not close xdr file')
      endif
      
      print 'No errors'
      
      return 0
      end
      

	program xdr_test
	
	implicit none
	integer xd, xd2

	integer i, j, num_of_coord, ret, framecnt
	real coord(10000)
	real coord2(10000)
	real d0, d1, d2, d3
	real prec
	
	prec = 1000.0
c	/* open file containing 3d coordinates */
	open(7, file="test.gmx", status="old")

c	/* read first line which contains number of coordinates */
	read(7, *) num_of_coord, d0, d1, d2, d3

c	/* open xdr file to which compressed coordinates will be written */
	call xdrfopen(xd, "test.xdr", "w", ret)
	
c	/* just as test write the first line using normal xdr routine */
	call xdrfint(xd, num_of_coord, ret)
	call xdrffloat(xd, d0, ret)
	call xdrffloat(xd, d1, ret)
	call xdrffloat(xd, d2, ret)
	call xdrffloat(xd, d3, ret)

	framecnt = 0
 10	continue
	call read_frame(7, num_of_coord, coord, ret)
	if (ret .eq. 0) goto 20
	call xdrf3dfcoord(xd, coord, num_of_coord, prec, ret)
	framecnt = framecnt + 1
	goto 10
 20	continue

	call xdrfclose(xd, ret)
	close(7)

c	/* Now do the inverse ! */

c	/* open file to write decompressed data */
	open(8, file="test.out", status="unknown")

	call xdrfopen(xd2, "test.xdr", "r", ret)
	call xdrfint(xd2, num_of_coord, ret)
	call xdrffloat(xd2, d0, ret)
	call xdrffloat(xd2, d1, ret)
	call xdrffloat(xd2, d2, ret)
	call xdrffloat(xd2, d3, ret)
	write(8,'(i5,f8.3,f8.3,f8.3,f8.3)') num_of_coord, d0, d1,
     &		d2, d3

	do 30, i=1,framecnt
	    call xdrf3dfcoord(xd2, coord2, num_of_coord, prec, ret)
	    do 40, j=1, num_of_coord * 3, 3
	    	write (8,'(f8.3,f9.3,f9.3)') coord2(j), coord2(j+1),
     &			coord2(j+2)
 40	    continue
 30	continue
	call xdrfclose(xd2, ret)
	close(8)
	end

	
	subroutine read_frame(in, num_of_coord,  coord, ret)

	implicit none
	integer in
	integer i, num_of_coord
	real coord(*)
	integer ret

	do 100, i=1, num_of_coord*3, 3
	    read (in, *, err = 120, end=120) coord(i), coord(i+1),
     &		coord(i+2)
 100	continue
	ret = 1
	return
 120	ret = 0
	return
	end

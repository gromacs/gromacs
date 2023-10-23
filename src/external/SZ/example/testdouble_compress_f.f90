program p
        use sz
        use rw
        implicit none
        character(len=32) :: arg
        real(kind=8), dimension(:,:,:), allocatable :: grid
        integer(kind=C_SIZE_T) :: gridsize1,gridsize2,gridsize3
        real(kind=8) :: res=0
        integer :: i,j,k
        integer(kind=4) :: ierr;
        integer(kind=C_SIZE_T) outSize !the size of the compressed stream
        INTEGER(kind=1), DIMENSION(:), allocatable :: Bytes
        gridsize1 = 10 
        gridsize2 = 10 
        gridsize3 = 10 
        
        write (6,*) 'start....'
        allocate(grid(gridsize1,gridsize2,gridsize3))
        DO i=1,gridsize1
                DO j=1,gridsize2
                        DO k=1,gridsize3
                                grid(i,j,k)=i+j+k
                        END DO
                END DO
        END DO
     
        call getarg(1, arg)
        call SZ_Init(arg,ierr)
       
        call SZ_Compress(grid, Bytes, outSize)
        call writeData(Bytes, outSize, 'test_f.sz') 
         
        ! Free memory
        deallocate(grid)
        deallocate(Bytes)
        call SZ_Finalize()
        write (6,*) 'done.'
end program p

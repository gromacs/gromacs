program p
        use sz
        use rw
        implicit none
        character(len=32) :: arg
        integer(kind=1), dimension(:), allocatable :: bytes
        real(kind=8), dimension(:,:,:), allocatable :: grid
        integer(kind=C_SIZE_T) :: gridsize1,gridsize2,gridsize3
        real(kind=8) :: res=0
        integer :: i,j,k
        integer(kind=4) :: ierr; 
        integer(kind=C_SIZE_T) outSize !the size of the compressed stream
        gridsize1 = 10 
        gridsize2 = 10 
        gridsize3 = 10 
        write (6,*) 'start....'
     
        call getarg(1, arg)
        call SZ_Init(arg,ierr)
      
        call readData('test_f.sz', bytes, outSize)
 
        call SZ_Decompress(bytes, grid, gridsize1, gridsize2, gridsize3)
        open(unit=10,file='test_f.txt')
        DO i=1,gridsize3
                DO j=1,gridsize2
                        DO k=1,gridsize1
                                write (10,*) grid(k,j,i)
                        END DO
                END DO
        END DO
        
        deallocate(grid)
        deallocate(bytes) 
        write (6,*) 'done.'
        call SZ_Finalize()
end program p

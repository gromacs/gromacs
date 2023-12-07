      program main

c*********************************************************************72
c
cc MAIN is the main program for MD_OPENMP.
c
c  Discussion:
c
c    The program implements a simple molecular dynamics simulation.
c
c    The program uses Open MP directives to allow parallel computation.
c
c    The velocity Verlet time integration scheme is used.
c
c    The particles interact with a central pair potential.
c
c    Output of the program is saved in the TNG format, which is why this
c    code is included in the TNG API release.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    8 Jan 2013
c
c  Author:
c
c    Original FORTRAN90 version by Bill Magro.
c    FORTRAN77 version by John Burkardt.
c    TNG trajectory output by Magnus Lundborg.
c
c  Parameters:
c
c    None
c
      use omp_lib

      implicit none

      integer nd
      parameter ( nd = 3 )
      integer np
      parameter ( np = 250 )
      integer step_num
      parameter ( step_num = 1000 )

      double precision acc(nd,np)
      double precision box(nd)
      double precision box_shape(9)
      double precision dt
      parameter ( dt = 0.0001D+00 )
      double precision e0
      double precision force(nd,np)
      integer i
      integer id
      double precision kinetic
      double precision mass
      parameter ( mass = 1.0D+00 )
      double precision pos(nd,np)
      double precision potential
      integer proc_num
      integer seed
      integer step
      integer step_print
      integer step_print_index
      integer step_print_num
      integer step_save
      integer*8 sparse_save
      integer thread_num
      double precision vel(nd,np)
      double precision wtime

c
c  Cray pointers are not standard fortran 77, but must be used to allocate
c  memory properly.
c
      pointer (traj_p, traj)
      pointer (molecule_p, molecule)
      pointer (chain_p, chain)
      pointer (residue_p, residue)
      pointer (atom_p, atom)
      byte traj
      byte molecule
      byte chain
      byte residue
      byte atom

c
c  The TNG functions expect 64 bit integers
c
      integer*8 n_frames_per_frame_set
      integer*8 frames_saved_cnt
      integer*8 tng_n_particles

c
c  These constants are also defined in tng_io.h, but need to
c  set in fortran as well. This can be copied to any fortran
c  source code using the tng_io library.
c
      integer*8 TNG_UNCOMPRESSED
      parameter ( TNG_UNCOMPRESSED = 0)
      integer TNG_NON_TRAJECTORY_BLOCK
      parameter ( TNG_NON_TRAJECTORY_BLOCK = 0)
      integer TNG_TRAJECTORY_BLOCK
      parameter ( TNG_TRAJECTORY_BLOCK = 1)
      integer*8 TNG_GENERAL_INFO
      parameter ( TNG_GENERAL_INFO = 0 )
      integer*8 TNG_MOLECULES
      parameter ( TNG_MOLECULES = 1 )
      integer*8 TNG_TRAJECTORY_FRAME_SET
      parameter ( TNG_TRAJECTORY_FRAME_SET = 2 )
      integer*8 TNG_PARTICLE_MAPPING
      parameter ( TNG_PARTICLE_MAPPING = 3 )
      integer*8 TNG_TRAJ_BOX_SHAPE
      parameter ( TNG_TRAJ_BOX_SHAPE = 10000 )
      integer*8 TNG_TRAJ_POSITIONS
      parameter ( TNG_TRAJ_POSITIONS = 10001 )
      integer*8 TNG_TRAJ_VELOCITIES
      parameter ( TNG_TRAJ_VELOCITIES = 10002 )
      integer*8 TNG_TRAJ_FORCES
      parameter ( TNG_TRAJ_FORCES = 10003 )
      integer TNG_SKIP_HASH
      parameter ( TNG_SKIP_HASH = 0 )
      integer TNG_USE_HASH
      parameter ( TNG_USE_HASH = 1 )
      integer TNG_CHAR_DATA
      parameter ( TNG_CHAR_DATA = 0 )
      integer TNG_INT_DATA
      parameter ( TNG_INT_DATA = 1 )
      integer TNG_FLOAT_DATA
      parameter ( TNG_FLOAT_DATA = 2 )
      integer TNG_DOUBLE_DATA
      parameter ( TNG_DOUBLE_DATA = 3 )

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MD_OPENMP'
      write ( *, '(a)' ) '  FORTRAN77/OpenMP version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A molecular dynamics program.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  NP, the number of particles in the simulation is ', np
      write ( *, '(a,i8)' )
     &  '  STEP_NUM, the number of time steps, is ', step_num
      write ( *, '(a,g14.6)' )
     &  '  DT, the size of each time step, is ', dt
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  The number of processors = ', omp_get_num_procs ( )
      write ( *, '(a,i8)' )
     &  '  The number of threads    = ', omp_get_max_threads ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Initializing trajectory storage.'
      call tng_trajectory_init(traj_p)

c
c  N.B. The TNG output file should be modified according to needs
c
      call tng_output_file_set(traj, TNG_EXAMPLE_FILES_DIR"tng_md_out_f77.tng")

      write ( *, '(a)' ) '  Creating molecules in trajectory.'
      tng_n_particles = np
      call tng_molecule_add(traj, "water", molecule_p)
      call tng_molecule_chain_add(traj, molecule, "W", chain_p)
      call tng_chain_residue_add(traj, chain, "WAT", residue_p)
      call tng_residue_atom_add(traj, residue, "O", "O", atom_p)
      call tng_molecule_cnt_set(traj, molecule, tng_n_particles)
c
c  Set the dimensions of the box.
c
      do i = 1, 9
        box_shape(i) = 0.0
      end do
      do i = 1, nd
        box(i) = 10.0D+00
        box_shape(i*nd + i) = box(i)
      end do

c
c  Add the box shape data block
c
      call tng_data_block_add(traj, TNG_TRAJ_BOX_SHAPE, "BOX SHAPE",
     &  TNG_DOUBLE_DATA, TNG_NON_TRAJECTORY_BLOCK, int(1, 8),
     &  int(9, 8), int(1, 8), TNG_UNCOMPRESSED, box_shape)

c
c  Write the file headers
c
      call tng_file_headers_write(traj, TNG_USE_HASH)

c
c  Set initial positions, velocities, and accelerations.
c
      write ( *, '(a)' )
     &  '  Initializing positions, velocities, and accelerations.'
      seed = 123456789
      call initialize ( np, nd, box, seed, pos, vel, acc )
c
c  Compute the forces and energies.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Computing initial forces and energies.'

      call compute ( np, nd, pos, vel, mass, force, potential,
     &  kinetic )

      e0 = potential + kinetic

c
c  Saving frequency
c
      step_save = 5

      step_print = 0
      step_print_index = 0
      step_print_num = 10
      sparse_save = 10

      frames_saved_cnt = 0

c
c  This is the main time stepping loop.
c
      write ( *, '(a,i4,a)' ) '  Every', step_save,
     &  ' steps particle positions, velocities and forces are'
      write ( *, '(a)' ) '  saved to a TNG trajectory file.'
      write ( *, '(a)' )
      write ( *, '(a)' )
     &  '  At each step, we report the potential and kinetic energies.'
      write ( *, '(a)' )
     &  '  The sum of these energies should be a constant.'
      write ( *, '(a)' )
     &  '  As an accuracy check, we also print the relative error'
      write ( *, '(a)' ) '  in the total energy.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '      Step      Potential       Kinetic        (P+K-E0)/E0'
      write ( *, '(a)' )
     &  '                Energy P        Energy K       ' //
     &  'Relative Energy Error'
      write ( *, '(a)' ) ' '

      step = 0
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' )
     &  step, potential, kinetic, ( potential + kinetic - e0 ) / e0
      step_print_index = step_print_index + 1
      step_print = ( step_print_index * step_num ) / step_print_num

c
c  Create a frame set for writing data
c
      call tng_num_frames_per_frame_set_get(traj,
     &  n_frames_per_frame_set)
      call tng_frame_set_new(traj, int(0, 8), n_frames_per_frame_set)

c
c  Add empty data blocks.
c
      call tng_particle_data_block_add(traj, TNG_TRAJ_POSITIONS,
     &  "POSITIONS", TNG_DOUBLE_DATA, TNG_TRAJECTORY_BLOCK,
     &  n_frames_per_frame_set, int(3, 8), int(1, 8), int(0, 8),
     &  tng_n_particles, TNG_UNCOMPRESSED, %VAL(int(0, 8)))

      call tng_particle_data_block_add(traj, TNG_TRAJ_VELOCITIES,
     &  "VELOCITIES", TNG_DOUBLE_DATA, TNG_TRAJECTORY_BLOCK,
     &  n_frames_per_frame_set, int(3, 8), int(1, 8), int(0, 8),
     &  tng_n_particles, TNG_UNCOMPRESSED, %VAL(int(0, 8)))

      call tng_particle_data_block_add(traj, TNG_TRAJ_FORCES,
     &  "FORCES", TNG_DOUBLE_DATA, TNG_TRAJECTORY_BLOCK,
     &  n_frames_per_frame_set, int(3, 8), int(1, 8), int(0, 8),
     &  tng_n_particles, TNG_UNCOMPRESSED, %VAL(int(0, 8)))

c
c  The potential energy data block is saved sparsely.
c
      call tng_data_block_add(traj, int(10101, 8),
     &  "POTENTIAL ENERGY", TNG_DOUBLE_DATA, TNG_TRAJECTORY_BLOCK,
     &  n_frames_per_frame_set, int(1, 8), sparse_save,
     &  TNG_UNCOMPRESSED, %VAL(int(0, 8)))


c
c  Write the frame set to disk
c
      call tng_frame_set_write(traj, TNG_USE_HASH)

      wtime = omp_get_wtime ( )

      do step = 1, step_num

        call compute ( np, nd, pos, vel, mass, force, potential,
     &    kinetic )

        if ( step .eq. step_print ) then

          write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' )
     &      step, potential, kinetic, ( potential + kinetic - e0 ) / e0

          step_print_index = step_print_index + 1
          step_print = ( step_print_index * step_num ) / step_print_num

        end if

c
c  Output to TNG file at regular intervals, specified by step_save
c
        if ( step_save .EQ. 0 .OR. mod(step, step_save) .EQ. 0 ) then
          call tng_frame_particle_data_write(traj, frames_saved_cnt,
     &    TNG_TRAJ_POSITIONS, int(0, 8), tng_n_particles, pos,
     &    TNG_USE_HASH)
          call tng_frame_particle_data_write(traj, frames_saved_cnt,
     &    TNG_TRAJ_VELOCITIES, int(0, 8), tng_n_particles, vel,
     &    TNG_USE_HASH)
          call tng_frame_particle_data_write(traj, frames_saved_cnt,
     &    TNG_TRAJ_FORCES, int(0, 8), tng_n_particles, force,
     &    TNG_USE_HASH)
          frames_saved_cnt = frames_saved_cnt + 1

          if (mod(step, step_save * sparse_save) .EQ. 0) then
            call tng_frame_data_write(traj, frames_saved_cnt,
     &      int(10101, 8), potential, TNG_USE_HASH)
          end if

        end if

        call update ( np, nd, pos, vel, force, acc, mass, dt )

      end do

      wtime = omp_get_wtime ( ) - wtime

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Elapsed time for main computation:'
      write ( *, '(2x,g14.6,a)' ) wtime, ' seconds'
c
c  Terminate.
c
      call tng_trajectory_destroy(traj_p)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MD_OPENMP'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine compute ( np, nd, pos, vel, mass, f, pot, kin )

c*********************************************************************72
c
cc COMPUTE computes the forces and energies.
c
c  Discussion:
c
c    The computation of forces and energies is fully parallel.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 July 2009
c
c  Author:
c
c    Original FORTRAN90 version by Bill Magro.
c    FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer NP, the number of particles.
c
c    Input, integer ND, the number of spatial dimensions.
c
c    Input, double precision POS(ND,NP), the position of each particle.
c
c    Input, double precision VEL(ND,NP), the velocity of each particle.
c
c    Input, double precision MASS, the mass of each particle.
c
      implicit none

      integer np
      integer nd

      double precision d
      double precision d2
      double precision dv
      double precision f(nd,np)
      integer i
      integer j
      integer k
      double precision kin
      double precision mass
      double precision PI2
      parameter ( PI2 = 3.141592653589793D+00 / 2.0D+00 )
      double precision pos(nd,np)
      double precision pot
      double precision rij(nd)
      double precision v
      double precision vel(nd,np)

      pot = 0.0D+00
      kin = 0.0D+00

c$omp parallel
c$omp& shared ( f, nd, np, pos, vel )
c$omp& private ( d, d2, i, j, k, rij )

c$omp do reduction ( + : pot, kin )
      do i = 1, np
c
c  Compute the potential energy and forces.
c
        do k = 1, nd
          f(k,i) = 0.0D+00
        end do

        do j = 1, np

          if ( i .ne. j ) then

            call dist ( nd, pos(1,i), pos(1,j), rij, d )
c
c  Attribute half of the potential energy to particle J.
c
            d2 = min ( d, pi2 )

            pot = pot + 0.5D+00 * ( sin ( d2 ) )**2

            do k = 1, nd
              f(k,i) = f(k,i) - rij(k) * sin ( 2.0D+00 * d2 ) / d
            end do

          end if

        end do
c
c  Compute the kinetic energy.
c
        do k = 1, nd
          kin = kin + vel(k,i)**2
        end do

      end do
c$omp end do

c$omp end parallel

      kin = kin * 0.5D+00 * mass

      return
      end
      subroutine dist ( nd, r1, r2, dr, d )

c*********************************************************************72
c
cc DIST computes the displacement (and its norm) between two particles.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2007
c
c  Author:
c
c    Original FORTRAN90 version by Bill Magro.
c    FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer ND, the number of spatial dimensions.
c
c    Input, double precision R1(ND), R2(ND), the positions of the particles.
c
c    Output, double precision DR(ND), the displacement vector.
c
c    Output, double precision D, the Euclidean norm of the displacement.
c
      implicit none

      integer nd

      double precision d
      double precision dr(nd)
      integer i
      double precision r1(nd)
      double precision r2(nd)

      do i = 1, nd
        dr(i) = r1(i) - r2(i)
      end do

      d = 0.0D+00
      do i = 1, nd
        d = d + dr(i)**2
      end do
      d = sqrt ( d )

      return
      end
      subroutine initialize ( np, nd, box, seed, pos, vel, acc )

c*********************************************************************72
c
cc INITIALIZE initializes the positions, velocities, and accelerations.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2007
c
c  Author:
c
c    Original FORTRAN90 version by Bill Magro.
c    FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer NP, the number of particles.
c
c    Input, integer ND, the number of spatial dimensions.
c
c    Input, double precision BOX(ND), specifies the maximum position
c    of particles in each dimension.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision POS(ND,NP), the position of each particle.
c
c    Output, double precision VEL(ND,NP), the velocity of each particle.
c
c    Output, double precision ACC(ND,NP), the acceleration of each particle.
c
      implicit none

      integer np
      integer nd

      double precision acc(nd,np)
      double precision box(nd)
      integer i
      integer j
      double precision pos(nd,np)
      double precision r8_uniform_01
      integer seed
      double precision vel(nd,np)
c
c  Give the particles random positions within the box.
c
      do i = 1, nd
        do j = 1, np
          pos(i,j) = r8_uniform_01 ( seed )
        end do
      end do

c$omp parallel
c$omp& shared ( acc, box, nd, np, pos, vel )
c$omp& private ( i, j )

c$omp do
      do j = 1, np
        do i = 1, nd
          pos(i,j) = box(i) * pos(i,j)
          vel(i,j) = 0.0D+00
          acc(i,j) = 0.0D+00
        end do
      end do
c$omp end do

c$omp end parallel

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision r8_uniform_01
      integer k
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
      subroutine update ( np, nd, pos, vel, f, acc, mass, dt )

c*********************************************************************72
c
cc UPDATE performs the time integration, using a velocity Verlet algorithm.
c
c  Discussion:
c
c    The time integration is fully parallel.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2007
c
c  Author:
c
c    Original FORTRAN90 version by Bill Magro.
c    FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer NP, the number of particles.
c
c    Input, integer ND, the number of spatial dimensions.
c
c    Input/output, double precision POS(ND,NP), the position of each particle.
c
c    Input/output, double precision VEL(ND,NP), the velocity of each particle.
c
c    Input, double precision MASS, the mass of each particle.
c
c    Input/output, double precision ACC(ND,NP), the acceleration of each
c    particle.
c
      implicit none

      integer np
      integer nd

      double precision acc(nd,np)
      double precision dt
      double precision f(nd,np)
      integer i
      integer j
      double precision mass
      double precision pos(nd,np)
      double precision rmass
      double precision vel(nd,np)

      rmass = 1.0D+00 / mass

c$omp parallel
c$omp& shared ( acc, dt, f, nd, np, pos, rmass, vel )
c$omp& private ( i, j )

c$omp do
      do j = 1, np
        do i = 1, nd

          pos(i,j) = pos(i,j)
     &      + vel(i,j) * dt + 0.5D+00 * acc(i,j) * dt * dt

          vel(i,j) = vel(i,j)
     &      + 0.5D+00 * dt * ( f(i,j) * rmass + acc(i,j) )

          acc(i,j) = f(i,j) * rmass

        end do
      end do
c$omp end do

c$omp end parallel

      return
      end

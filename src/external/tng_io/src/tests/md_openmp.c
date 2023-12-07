/*
 * This code is part of the tng binary trajectory format.
 *
 * Copyright (c) 2020, by the GROMACS development team.
 * TNG was orginally written by Magnus Lundborg, Daniel Spångberg and
 * Rossen Apostolov. The API is implemented mainly by Magnus Lundborg,
 * Daniel Spångberg and Anders Gärdenäs.
 *
 * Please see the AUTHORS file for more information.
 *
 * The TNG library is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 *
 * To help us fund future development, we humbly ask that you cite
 * the research papers on the package.
 *
 * Check out http://www.gromacs.org for more information.
 */
#ifdef TNG_BUILD_OPENMP_EXAMPLES

#    include "tng/tng_io.h"
#    include <stdlib.h>
#    include <stdio.h>
#    include <time.h>
#    include <math.h>
#    include <omp.h>

int    main();
void   compute(int np, int nd, double pos[], const double vel[], double mass, double f[], double* pot, double* kin);
double dist(int nd, const double r1[], const double r2[], double dr[]);
void initialize(int np, int nd, const double box[], int* seed, double pos[], double vel[], double acc[]);
double r8_uniform_01(int* seed);
void   timestamp(void);
void   update(int np, int nd, double pos[], double vel[], const double f[], double acc[], double mass, double dt);

/******************************************************************************/

int main()

/******************************************************************************/
/*
    Purpose:

        MAIN is the main program for MD_OPENMP.

    Discussion:

        MD implements a simple molecular dynamics simulation.

        The program uses Open MP directives to allow parallel computation.

        The velocity Verlet time integration scheme is used.

        The particles interact with a central pair potential.

        Output of the program is saved in the TNG format, which is why this
        code is included in the TNG API release.

    Licensing:

        This code is distributed under the GNU LGPL license.

    Modified:

        8 Jan 2013

    Author:

        Original FORTRAN77 version by Bill Magro.
        C version by John Burkardt.
        TNG trajectory output by Magnus Lundborg.

    Parameters:

        None
*/
{
    double*          acc;
    double*          box;
    double*          box_shape;
    double           dt = 0.0002;
    double           e0;
    double*          force;
    int              i;
    double           kinetic;
    double           mass = 1.0;
    int              nd   = 3;
    int              np   = 50;
    double*          pos;
    double           potential;
    int              proc_num;
    int              seed = 123456789;
    int              step;
    int              step_num = 50000;
    int              step_print;
    int              step_print_index;
    int              step_print_num;
    int              step_save;
    int64_t          sparse_save;
    double*          vel;
    double           wtime;
    tng_trajectory_t traj;
    tng_molecule_t   molecule;
    tng_chain_t      chain;
    tng_residue_t    residue;
    tng_atom_t       atom;
    int64_t          n_frames_per_frame_set;
    int              frames_saved_cnt = 0;

    timestamp();

    proc_num = omp_get_num_procs();

    acc       = (double*)malloc(nd * np * sizeof(double));
    box       = (double*)malloc(nd * sizeof(double));
    box_shape = (double*)malloc(9 * sizeof(double));
    force     = (double*)malloc(nd * np * sizeof(double));
    pos       = (double*)malloc(nd * np * sizeof(double));
    vel       = (double*)malloc(nd * np * sizeof(double));

    printf("\n");
    printf("MD_OPENMP\n");
    printf("  C/OpenMP version\n");
    printf("\n");
    printf("  A molecular dynamics program.\n");

    printf("\n");
    printf("  NP, the number of particles in the simulation is %d\n", np);
    printf("  STEP_NUM, the number of time steps, is %d\n", step_num);
    printf("  DT, the size of each time step, is %f\n", dt);

    printf("\n");
    printf("  Number of processors available = %d\n", proc_num);
    printf("  Number of threads =              %d\n", omp_get_max_threads());


    printf("\n");
    printf("  Initializing trajectory storage.\n");
    if (tng_trajectory_init(&traj) != TNG_SUCCESS)
    {
        tng_trajectory_destroy(&traj);
        printf("  Cannot init trajectory.\n");
        exit(1);
    }
    tng_output_file_set(traj, TNG_EXAMPLE_FILES_DIR "tng_md_out.tng");


    /* Set molecules data */
    printf("  Creating molecules in trajectory.\n");
    tng_molecule_add(traj, "water", &molecule);
    tng_molecule_chain_add(traj, molecule, "W", &chain);
    tng_chain_residue_add(traj, chain, "WAT", &residue);
    if (tng_residue_atom_add(traj, residue, "O", "O", &atom) == TNG_CRITICAL)
    {
        tng_trajectory_destroy(&traj);
        printf("  Cannot create molecules.\n");
        exit(1);
    }
    tng_molecule_cnt_set(traj, molecule, np);


    /*
        Set the dimensions of the box.
    */
    for (i = 0; i < 9; i++)
    {
        box_shape[i] = 0.0;
    }
    for (i = 0; i < nd; i++)
    {
        box[i]                = 10.0;
        box_shape[i * nd + i] = box[i];
    }


    /* Add the box shape data block and write the file headers */
    if (tng_data_block_add(traj, TNG_TRAJ_BOX_SHAPE, "BOX SHAPE", TNG_DOUBLE_DATA,
                           TNG_NON_TRAJECTORY_BLOCK, 1, 9, 1, TNG_UNCOMPRESSED, box_shape)
                == TNG_CRITICAL
        || tng_file_headers_write(traj, TNG_USE_HASH) == TNG_CRITICAL)
    {
        free(box_shape);
        tng_trajectory_destroy(&traj);
        printf("  Cannot write trajectory headers and box shape.\n");
        exit(1);
    }
    free(box_shape);

    printf("\n");
    printf("  Initializing positions, velocities, and accelerations.\n");
    /*
        Set initial positions, velocities, and accelerations.
    */
    initialize(np, nd, box, &seed, pos, vel, acc);
    /*
        Compute the forces and energies.
    */
    printf("\n");
    printf("  Computing initial forces and energies.\n");

    compute(np, nd, pos, vel, mass, force, &potential, &kinetic);

    e0 = potential + kinetic;

    /* Saving frequency */
    step_save = 500;

    step_print       = 0;
    step_print_index = 0;
    step_print_num   = 10;
    sparse_save      = 100;

    /*
        This is the main time stepping loop:
            Compute forces and energies,
            Update positions, velocities, accelerations.
    */
    printf("  Every %d steps particle positions, velocities and forces are\n", step_save);
    printf("  saved to a TNG trajectory file.\n");
    printf("\n");
    printf("  At certain step intervals, we report the potential and kinetic energies.\n");
    printf("  The sum of these energies should be a constant.\n");
    printf("  As an accuracy check, we also print the relative error\n");
    printf("  in the total energy.\n");
    printf("\n");
    printf("      Step      Potential       Kinetic        (P+K-E0)/E0\n");
    printf("                Energy P        Energy K       Relative Energy Error\n");
    printf("\n");

    step = 0;
    printf("  %8d  %14f  %14f  %14e\n", step, potential, kinetic, (potential + kinetic - e0) / e0);
    step_print_index++;
    step_print = (step_print_index * step_num) / step_print_num;

    /* Create a frame set for writing data */
    tng_num_frames_per_frame_set_get(traj, &n_frames_per_frame_set);
    if (tng_frame_set_new(traj, 0, n_frames_per_frame_set) != TNG_SUCCESS)
    {
        printf("Error creating frame set %d. %s: %d\n", i, __FILE__, __LINE__);
        exit(1);
    }

    /* Add empty data blocks */
    if (tng_particle_data_block_add(traj, TNG_TRAJ_POSITIONS, "POSITIONS", TNG_DOUBLE_DATA, TNG_TRAJECTORY_BLOCK,
                                    n_frames_per_frame_set, 3, 1, 0, np, TNG_UNCOMPRESSED, 0)
        != TNG_SUCCESS)
    {
        printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }
    if (tng_particle_data_block_add(traj, TNG_TRAJ_VELOCITIES, "VELOCITIES", TNG_DOUBLE_DATA, TNG_TRAJECTORY_BLOCK,
                                    n_frames_per_frame_set, 3, 1, 0, np, TNG_UNCOMPRESSED, 0)
        != TNG_SUCCESS)
    {
        printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }
    if (tng_particle_data_block_add(traj, TNG_TRAJ_FORCES, "FORCES", TNG_DOUBLE_DATA, TNG_TRAJECTORY_BLOCK,
                                    n_frames_per_frame_set, 3, 1, 0, np, TNG_UNCOMPRESSED, 0)
        != TNG_SUCCESS)
    {
        printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }

    /* There is no standard ID for potential energy. Pick one. The
       potential energy will not be saved every frame - it is sparsely
       saved. */
    if (tng_data_block_add(traj, 10101, "POTENTIAL ENERGY", TNG_DOUBLE_DATA, TNG_TRAJECTORY_BLOCK,
                           n_frames_per_frame_set, 1, sparse_save, TNG_UNCOMPRESSED, 0)
        != TNG_SUCCESS)
    {
        printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }

    /* Write the frame set to disk */
    if (tng_frame_set_write(traj, TNG_USE_HASH) != TNG_SUCCESS)
    {
        printf("Error writing frame set. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }

    wtime = omp_get_wtime();

    if (tng_frame_particle_data_write(traj, frames_saved_cnt, TNG_TRAJ_POSITIONS, 0, np, pos, TNG_USE_HASH)
        != TNG_SUCCESS)
    {
        printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }
    if (tng_frame_particle_data_write(traj, frames_saved_cnt, TNG_TRAJ_VELOCITIES, 0, np, vel, TNG_USE_HASH)
        != TNG_SUCCESS)
    {
        printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }
    if (tng_frame_particle_data_write(traj, frames_saved_cnt, TNG_TRAJ_FORCES, 0, np, force, TNG_USE_HASH)
        != TNG_SUCCESS)
    {
        printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }
    if (step % (step_save * sparse_save) == 0)
    {
        if (tng_frame_data_write(traj, frames_saved_cnt, 10101, &potential, TNG_USE_HASH) != TNG_SUCCESS)
        {
            printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
            exit(1);
        }
    }
    frames_saved_cnt++;

    for (step = 1; step < step_num; step++)
    {
        compute(np, nd, pos, vel, mass, force, &potential, &kinetic);

        if (step == step_print)
        {
            printf("  %8d  %14f  %14f  %14e\n", step, potential, kinetic, (potential + kinetic - e0) / e0);
            step_print_index++;
            step_print = (step_print_index * step_num) / step_print_num;
        }
        if (step % step_save == 0)
        {
            if (tng_frame_particle_data_write(traj, frames_saved_cnt, TNG_TRAJ_POSITIONS, 0, np, pos, TNG_USE_HASH)
                != TNG_SUCCESS)
            {
                printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
                exit(1);
            }
            if (tng_frame_particle_data_write(traj, frames_saved_cnt, TNG_TRAJ_VELOCITIES, 0, np, vel, TNG_USE_HASH)
                != TNG_SUCCESS)
            {
                printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
                exit(1);
            }
            if (tng_frame_particle_data_write(traj, frames_saved_cnt, TNG_TRAJ_FORCES, 0, np, force, TNG_USE_HASH)
                != TNG_SUCCESS)
            {
                printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
                exit(1);
            }
            if (step % (step_save * sparse_save) == 0)
            {
                if (tng_frame_data_write(traj, frames_saved_cnt, 10101, &potential, TNG_USE_HASH) != TNG_SUCCESS)
                {
                    printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
                    exit(1);
                }
            }
            frames_saved_cnt++;
        }
        update(np, nd, pos, vel, force, acc, mass, dt);
    }
    wtime = omp_get_wtime() - wtime;

    printf("\n");
    printf("  Elapsed time for main computation:\n");
    printf("  %f seconds.\n", wtime);

    free(acc);
    free(box);
    free(force);
    free(pos);
    free(vel);
    /*
        Terminate.
    */
    tng_trajectory_destroy(&traj);

    printf("\n");
    printf("MD_OPENMP\n");
    printf("  Normal end of execution.\n");

    printf("\n");
    timestamp();

    return 0;
}
/******************************************************************************/

void compute(int np, int nd, double pos[], const double vel[], double mass, double f[], double* pot, double* kin)

/******************************************************************************/
/*
    Purpose:

        COMPUTE computes the forces and energies.

    Discussion:

        The computation of forces and energies is fully parallel.

        The potential function V(X) is a harmonic well which smoothly
        saturates to a maximum value at PI/2:

            v(x) = ( sin ( min ( x, PI2 ) ) )**2

        The derivative of the potential is:

            dv(x) = 2.0 * sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
                        = sin ( 2.0 * min ( x, PI2 ) )

    Licensing:

        This code is distributed under the GNU LGPL license.

    Modified:

        21 November 2007

    Author:

        Original FORTRAN77 version by Bill Magro.
        C version by John Burkardt.

    Parameters:

        Input, int NP, the number of particles.

        Input, int ND, the number of spatial dimensions.

        Input, double POS[ND*NP], the position of each particle.

        Input, double VEL[ND*NP], the velocity of each particle.

        Input, double MASS, the mass of each particle.

        Output, double F[ND*NP], the forces.

        Output, double *POT, the total potential energy.

        Output, double *KIN, the total kinetic energy.
*/
{
    double d;
    double d2;
    int    i;
    int    j;
    int    k = 0;
    double ke;
    double pe;
    double PI2 = 3.141592653589793 / 2.0;
    double rij[3];

    pe = 0.0;
    ke = 0.0;

#    pragma omp parallel shared(f, nd, np, pos, vel) private(i, j, k, rij, d, d2)


#    pragma omp for reduction(+ : pe, ke)
    for (k = 0; k < np; k++)
    {
        /*
            Compute the potential energy and forces.
        */
        for (i = 0; i < nd; i++)
        {
            f[i + k * nd] = 0.0;
        }

        for (j = 0; j < np; j++)
        {
            if (k != j)
            {
                d = dist(nd, pos + k * nd, pos + j * nd, rij);
                /*
                    Attribute half of the potential energy to particle J.
                */
                if (d < PI2)
                {
                    d2 = d;
                }
                else
                {
                    d2 = PI2;
                }

                pe = pe + 0.5 * pow(sin(d2), 2);

                for (i = 0; i < nd; i++)
                {
                    f[i + k * nd] = f[i + k * nd] - rij[i] * sin(2.0 * d2) / d;
                }
            }
        }
        /*
            Compute the kinetic energy.
        */
        for (i = 0; i < nd; i++)
        {
            ke = ke + vel[i + k * nd] * vel[i + k * nd];
        }
    }

    ke = ke * 0.5 * mass;

    *pot = pe;
    *kin = ke;
}
/******************************************************************************/

double dist(int nd, const double r1[], const double r2[], double dr[])

/******************************************************************************/
/*
    Purpose:

        DIST computes the displacement (and its norm) between two particles.

    Licensing:

        This code is distributed under the GNU LGPL license.

    Modified:

        21 November 2007

    Author:

        Original FORTRAN77 version by Bill Magro.
        C version by John Burkardt.

    Parameters:

        Input, int ND, the number of spatial dimensions.

        Input, double R1[ND], R2[ND], the positions of the particles.

        Output, double DR[ND], the displacement vector.

        Output, double D, the Euclidean norm of the displacement.
*/
{
    double d;
    int    i;

    d = 0.0;
    for (i = 0; i < nd; i++)
    {
        dr[i] = r1[i] - r2[i];
        d     = d + dr[i] * dr[i];
    }
    d = sqrt(d);

    return d;
}
/******************************************************************************/

void initialize(int np, int nd, const double box[], int* seed, double pos[], double vel[], double acc[])

/******************************************************************************/
/*
    Purpose:

        INITIALIZE initializes the positions, velocities, and accelerations.

    Licensing:

        This code is distributed under the GNU LGPL license.

    Modified:

        21 November 2007

    Author:

        Original FORTRAN77 version by Bill Magro.
        C version by John Burkardt.

    Parameters:

        Input, int NP, the number of particles.

        Input, int ND, the number of spatial dimensions.

        Input, double BOX[ND], specifies the maximum position
        of particles in each dimension.

        Input, int *SEED, a seed for the random number generator.

        Output, double POS[ND*NP], the position of each particle.

        Output, double VEL[ND*NP], the velocity of each particle.

        Output, double ACC[ND*NP], the acceleration of each particle.
*/
{
    int i;
    int j;
    /*
        Give the particles random positions within the box.
    */
    for (i = 0; i < nd; i++)
    {
        for (j = 0; j < np; j++)
        {
            pos[i + j * nd] = box[i] * r8_uniform_01(seed);
        }
    }

    for (j = 0; j < np; j++)
    {
        for (i = 0; i < nd; i++)
        {
            vel[i + j * nd] = 0.0;
        }
    }
    for (j = 0; j < np; j++)
    {
        for (i = 0; i < nd; i++)
        {
            acc[i + j * nd] = 0.0;
        }
    }
}
/******************************************************************************/

double r8_uniform_01(int* seed)

/******************************************************************************/
/*
    Purpose:

        R8_UNIFORM_01 is a unit pseudorandom R8.

    Discussion:

        This routine implements the recursion

            seed = 16807 * seed mod ( 2**31 - 1 )
            unif = seed / ( 2**31 - 1 )

        The integer arithmetic never requires more than 32 bits,
        including a sign bit.

    Licensing:

        This code is distributed under the GNU LGPL license.

    Modified:

        11 August 2004

    Author:

        John Burkardt

    Reference:

        Paul Bratley, Bennett Fox, Linus Schrage,
        A Guide to Simulation,
        Springer Verlag, pages 201-202, 1983.

        Bennett Fox,
        Algorithm 647:
        Implementation and Relative Efficiency of Quasirandom
        Sequence Generators,
        ACM Transactions on Mathematical Software,
        Volume 12, Number 4, pages 362-376, 1986.

    Parameters:

        Input/output, int *SEED, a seed for the random number generator.

        Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
        0 and 1.
*/
{
    int    k;
    double r;

    k = *seed / 127773;

    *seed = 16807 * (*seed - k * 127773) - k * 2836;

    if (*seed < 0)
    {
        *seed = *seed + 2147483647;
    }

    r = (double)(*seed) * 4.656612875E-10;

    return r;
}
/******************************************************************************/

void timestamp(void)

/******************************************************************************/
/*
    Purpose:

        TIMESTAMP prints the current YMDHMS date as a time stamp.

    Example:

        31 May 2001 09:45:54 AM

    Licensing:

        This code is distributed under the GNU LGPL license.

    Modified:

        24 September 2003

    Author:

        John Burkardt

    Parameters:

        None
*/
{
#    define TIME_SIZE 40

    static char      time_buffer[TIME_SIZE];
    const struct tm* tm;
    time_t           now;

    now = time(NULL);
    tm  = localtime(&now);

    strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

    printf("%s\n", time_buffer);
#    undef TIME_SIZE
}
/******************************************************************************/

void update(int np, int nd, double pos[], double vel[], const double f[], double acc[], double mass, double dt)

/******************************************************************************/
/*
    Purpose:

        UPDATE updates positions, velocities and accelerations.

    Discussion:

        The time integration is fully parallel.

        A velocity Verlet algorithm is used for the updating.

        x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
        v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
        a(t+dt) = f(t) / m

    Licensing:

        This code is distributed under the GNU LGPL license.

    Modified:

        17 April 2009

    Author:

        Original FORTRAN77 version by Bill Magro.
        C version by John Burkardt.

    Parameters:

        Input, int NP, the number of particles.

        Input, int ND, the number of spatial dimensions.

        Input/output, double POS[ND*NP], the position of each particle.

        Input/output, double VEL[ND*NP], the velocity of each particle.

        Input, double F[ND*NP], the force on each particle.

        Input/output, double ACC[ND*NP], the acceleration of each particle.

        Input, double MASS, the mass of each particle.

        Input, double DT, the time step.
*/
{
    int    i;
    int    j = 0;
    double rmass;

    rmass = 1.0 / mass;

#    pragma omp parallel shared(acc, dt, f, nd, np, pos, rmass, vel) private(i, j)

#    pragma omp for
    for (j = 0; j < np; j++)
    {
        for (i = 0; i < nd; i++)
        {
            pos[i + j * nd] = pos[i + j * nd] + vel[i + j * nd] * dt + 0.5 * acc[i + j * nd] * dt * dt;
            vel[i + j * nd] = vel[i + j * nd] + 0.5 * dt * (f[i + j * nd] * rmass + acc[i + j * nd]);
            acc[i + j * nd] = f[i + j * nd] * rmass;
        }
    }
}

#endif

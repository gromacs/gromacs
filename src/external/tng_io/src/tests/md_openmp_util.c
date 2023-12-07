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

int   main();
void  compute(int np, int nd, float pos[], const float vel[], float mass, float f[], float* pot, float* kin);
float dist(int nd, const float r1[], const float r2[], float dr[]);
void initialize(int np, int nd, const float box[], int* seed, float pos[], float vel[], float acc[]);
float r8_uniform_01(int* seed);
void  timestamp(void);
void  update(int np, int nd, float pos[], float vel[], const float f[], float acc[], float mass, float dt);

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
        code is included in the TNG API release. The high-level API of the
        TNG API is used where appropriate.

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
    float*           acc;
    float*           box;
    float*           box_shape;
    float            dt = 0.0002;
    float            e0;
    float*           force;
    int              i;
    float            kinetic;
    float            mass = 1.0;
    int              nd   = 3;
    int              np   = 50;
    float*           pos;
    float            potential;
    int              proc_num;
    int              seed = 123456789;
    int              step;
    int              step_num = 50000;
    int              step_print;
    int              step_print_index;
    int              step_print_num;
    int              step_save;
    float*           vel;
    float            wtime;
    tng_trajectory_t traj;
    tng_molecule_t   molecule;
    tng_chain_t      chain;
    tng_residue_t    residue;
    tng_atom_t       atom;

    timestamp();

    proc_num = omp_get_num_procs();

    acc       = (float*)malloc(nd * np * sizeof(float));
    box       = (float*)malloc(nd * sizeof(float));
    box_shape = (float*)malloc(9 * sizeof(float));
    force     = (float*)malloc(nd * np * sizeof(float));
    pos       = (float*)malloc(nd * np * sizeof(float));
    vel       = (float*)malloc(nd * np * sizeof(float));

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
    /* Initialize the TNG trajectory */
    tng_util_trajectory_open(TNG_EXAMPLE_FILES_DIR "tng_md_out.tng", 'w', &traj);


    /* Set molecules data */
    /* N.B. This is still not done using utility functions. The low-level API
     * is used. */
    printf("  Creating molecules in trajectory.\n");
    tng_molecule_add(traj, "water", &molecule);
    tng_molecule_chain_add(traj, molecule, "W", &chain);
    tng_chain_residue_add(traj, chain, "WAT", &residue);
    if (tng_residue_atom_add(traj, residue, "O", "O", &atom) == TNG_CRITICAL)
    {
        tng_util_trajectory_close(&traj);
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
        box[i] = 10.0;
        /* box_shape stores 9 values according to the TNG specs */
        box_shape[i * nd + i] = box[i];
    }


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
    step_save = 400;

    step_print       = 0;
    step_print_index = 0;
    step_print_num   = 10;

    /*
        This is the main time stepping loop:
            Compute forces and energies,
            Update positions, velocities, accelerations.
    */
    printf("  Every %d steps box shape, particle positions, velocities and forces are\n", step_save);
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

    /* Set the output frequency of box shape, positions, velocities and forces */
    if (tng_util_box_shape_write_frequency_set(traj, step_save) != TNG_SUCCESS)
    {
        printf("Error setting writing frequency data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }
    if (tng_util_pos_write_frequency_set(traj, step_save) != TNG_SUCCESS)
    {
        printf("Error setting writing frequency data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }
    if (tng_util_vel_write_frequency_set(traj, step_save) != TNG_SUCCESS)
    {
        printf("Error setting writing frequency data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }
    if (tng_util_force_write_frequency_set(traj, step_save) != TNG_SUCCESS)
    {
        printf("Error setting writing frequency data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }

    /* Write the first frame of box shape, positions, velocities and forces */
    if (tng_util_box_shape_write(traj, 0, box_shape) != TNG_SUCCESS)
    {
        printf("Error writing box shape. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }
    if (tng_util_pos_write(traj, 0, pos) != TNG_SUCCESS)
    {
        printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }
    if (tng_util_vel_write(traj, 0, vel) != TNG_SUCCESS)
    {
        printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }
    if (tng_util_force_write(traj, 0, force) != TNG_SUCCESS)
    {
        printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
        exit(1);
    }

    wtime = omp_get_wtime();

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
            /* Write box shape, positions, velocities and forces */
            if (tng_util_box_shape_write(traj, step, box_shape) != TNG_SUCCESS)
            {
                printf("Error writing box shape. %s: %d\n", __FILE__, __LINE__);
                exit(1);
            }
            if (tng_util_pos_write(traj, step, pos) != TNG_SUCCESS)
            {
                printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
                break;
            }
            if (tng_util_vel_write(traj, step, vel) != TNG_SUCCESS)
            {
                printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
                break;
            }
            if (tng_util_force_write(traj, step, force) != TNG_SUCCESS)
            {
                printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
                break;
            }
        }
        update(np, nd, pos, vel, force, acc, mass, dt);
    }
    wtime = omp_get_wtime() - wtime;

    printf("\n");
    printf("  Elapsed time for main computation:\n");
    printf("  %f seconds.\n", wtime);

    free(acc);
    free(box);
    free(box_shape);
    free(force);
    free(pos);
    free(vel);

    /* Close the TNG output. */
    tng_util_trajectory_close(&traj);

    printf("\n");
    printf("MD_OPENMP\n");
    printf("  Normal end of execution.\n");

    printf("\n");
    timestamp();

    return 0;
}
/******************************************************************************/

void compute(int np, int nd, float pos[], const float vel[], float mass, float f[], float* pot, float* kin)

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

        Input, float POS[ND*NP], the position of each particle.

        Input, float VEL[ND*NP], the velocity of each particle.

        Input, float MASS, the mass of each particle.

        Output, float F[ND*NP], the forces.

        Output, float *POT, the total potential energy.

        Output, float *KIN, the total kinetic energy.
*/
{
    float d;
    float d2;
    int   i;
    int   j;
    int   k = 0;
    float ke;
    float pe;
    float PI2 = 3.141592653589793 / 2.0;
    float rij[3];

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

                pe = pe + 0.5 * pow(sinf(d2), 2);

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

float dist(int nd, const float r1[], const float r2[], float dr[])

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

        Input, float R1[ND], R2[ND], the positions of the particles.

        Output, float DR[ND], the displacement vector.

        Output, float D, the Euclidean norm of the displacement.
*/
{
    float d;
    int   i;

    d = 0.0;
    for (i = 0; i < nd; i++)
    {
        dr[i] = r1[i] - r2[i];
        d     = d + dr[i] * dr[i];
    }
    d = sqrtf(d);

    return d;
}
/******************************************************************************/

void initialize(int np, int nd, const float box[], int* seed, float pos[], float vel[], float acc[])

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

        Input, float BOX[ND], specifies the maximum position
        of particles in each dimension.

        Input, int *SEED, a seed for the random number generator.

        Output, float POS[ND*NP], the position of each particle.

        Output, float VEL[ND*NP], the velocity of each particle.

        Output, float ACC[ND*NP], the acceleration of each particle.
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

float r8_uniform_01(int* seed)

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

        Output, float R8_UNIFORM_01, a new pseudorandom variate, strictly between
        0 and 1.
*/
{
    int   k;
    float r;

    k = *seed / 127773;

    *seed = 16807 * (*seed - k * 127773) - k * 2836;

    if (*seed < 0)
    {
        *seed = *seed + 2147483647;
    }

    r = (float)(*seed) * 4.656612875E-10;

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

void update(int np, int nd, float pos[], float vel[], const float f[], float acc[], float mass, float dt)

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

        Input/output, float POS[ND*NP], the position of each particle.

        Input/output, float VEL[ND*NP], the velocity of each particle.

        Input, float F[ND*NP], the force on each particle.

        Input/output, float ACC[ND*NP], the acceleration of each particle.

        Input, float MASS, the mass of each particle.

        Input, float DT, the time step.
*/
{
    int   i;
    int   j = 0;
    float rmass;

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

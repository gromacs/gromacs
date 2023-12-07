#ifdef TNG_BUILD_OPENMP_EXAMPLES

/* This code is part of the tng binary trajectory format.
 *
 * Written by Magnus Lundborg
 * Copyright (c) 2012-2013, The GROMACS development team.
 * Check out http://www.gromacs.org for more information.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */

#    include "tng/tng_io.h"

#    include <stdlib.h>
#    include <stdio.h>


/* N.B. this code is for testing parallel reading of trajectory frame sets. The
 * performance is not improved very much and is to a large extent limited by
 * disk i/o. It can however be used as inspiration for writing parallel code
 * using the TNG library. The code is NOT fully tested and may behave strangely. */

int main(int argc, char** argv)
{
    tng_trajectory_t           traj, local_traj = 0;
    union data_values***       local_positions = 0; // A 3-dimensional array to be populated
    union data_values**        particle_pos    = 0;
    int64_t                    n_particles, n_values_per_frame, n_frame_sets, n_frames;
    int64_t                    n_frames_per_frame_set, tot_n_frames = 0;
    char                       data_type;
    int                        i, j, fail;
    int64_t                    particle = 0, local_first_frame, local_last_frame;
    char                       atom_name[64], res_name[64];
    tng_trajectory_frame_set_t frame_set = 0;

    if (argc <= 1)
    {
        printf("No file specified\n");
        printf("Usage:\n");
        printf("tng_parallel_read <tng_file> [particle number = %" PRId64 "]\n", particle);
        exit(1);
    }

    // A reference must be passed to allocate memory
    if (tng_trajectory_init(&traj) != TNG_SUCCESS)
    {
        tng_trajectory_destroy(&traj);
        exit(1);
    }
    tng_input_file_set(traj, argv[1]);

    tng_current_frame_set_get(traj, &frame_set);

    // Read the file headers
    tng_file_headers_read(traj, TNG_USE_HASH);

    if (argc >= 3)
    {
        particle = strtoll(argv[2], 0, 10);
    }

    tng_num_frame_sets_get(traj, &n_frame_sets);
    tng_num_frames_per_frame_set_get(traj, &n_frames_per_frame_set);

    particle_pos = malloc(sizeof(union data_values*) * n_frame_sets * n_frames_per_frame_set);
    for (i = n_frame_sets * n_frames_per_frame_set; i--;)
    {
        /* Assume 3 values per frame even if it's not determined yet */
        particle_pos[i] = malloc(sizeof(union data_values) * 3);
    }

    printf("%" PRId64 " frame sets\n", n_frame_sets);

    if (tng_atom_name_of_particle_nr_get(traj, particle, atom_name, sizeof(atom_name)) == TNG_SUCCESS
        && tng_residue_name_of_particle_nr_get(traj, particle, res_name, sizeof(res_name)) == TNG_SUCCESS)
    {
        printf("Particle: %s (%s)\n", atom_name, res_name);
    }
    else
    {
        printf("Particle name not found\n");
    }

    fail = 0;

#    pragma omp parallel private(n_frames, n_particles, n_values_per_frame, local_first_frame, \
                                 local_last_frame, j,                                          \
                                 fail) firstprivate(local_traj, local_positions, frame_set)    \
            shared(data_type, traj, n_frame_sets, particle_pos, particle, i, tot_n_frames) default(none)
    {
        /* Each tng_trajectory_t keeps its own file pointers and i/o positions.
         * Therefore there must be a copy for each thread. */
        tng_trajectory_init_from_src(traj, &local_traj);
#    pragma omp for
        for (i = 0; i < n_frame_sets; i++)
        {
            if (tng_frame_set_nr_find(local_traj, i) != TNG_SUCCESS)
            {
                printf("FAILED finding frame set %d!\n", i);
                tot_n_frames = 0;
                fail         = 1;
            }
            if (tng_particle_data_get(local_traj, TNG_TRAJ_POSITIONS, &local_positions, &n_frames,
                                      &n_particles, &n_values_per_frame, &data_type)
                != TNG_SUCCESS)
            {
                printf("FAILED getting particle data\n");
                tot_n_frames = 0;
                fail         = 1;
            }
            if (!fail)
            {
                tng_current_frame_set_get(local_traj, &frame_set);
                tng_frame_set_frame_range_get(local_traj, frame_set, &local_first_frame, &local_last_frame);
                //         printf("Frame %"PRId64"-%"PRId64":\n", local_first_frame, local_last_frame);
                //         printf("%"PRId64" %"PRId64" %"PRId64"\n", n_frames, n_particles, n_values_per_frame);
                tot_n_frames += n_frames;
                for (j = 0; j < n_frames; j++)
                {
                    particle_pos[local_first_frame + j][0] = local_positions[j][particle][0];
                    particle_pos[local_first_frame + j][1] = local_positions[j][particle][1];
                    particle_pos[local_first_frame + j][2] = local_positions[j][particle][2];
                }
            }
        }

        // Free memory
        if (local_positions)
        {
            tng_particle_data_values_free(local_traj, local_positions, n_frames, n_particles,
                                          n_values_per_frame, data_type);
        }
        tng_trajectory_destroy(&local_traj);
    }
    switch (data_type)
    {
        case TNG_INT_DATA:
            for (j = 0; j < tot_n_frames; j++)
            {
                printf("\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\n", particle_pos[j][0].i,
                       particle_pos[j][1].i, particle_pos[j][2].i);
            }
            break;
        case TNG_FLOAT_DATA:
            for (j = 0; j < tot_n_frames; j++)
            {
                printf("\t%f\t%f\t%f\n", particle_pos[j][0].f, particle_pos[j][1].f,
                       particle_pos[j][2].f);
            }
            break;
        case TNG_DOUBLE_DATA:
            for (j = 0; j < tot_n_frames; j++)
            {
                printf("\t%f\t%f\t%f\n", particle_pos[j][0].d, particle_pos[j][1].d,
                       particle_pos[j][2].d);
            }
            break;
        default: break;
    }

    /* Free more memory */
    for (i = n_frame_sets * n_frames_per_frame_set; i--;)
    {
        free(particle_pos[i]);
    }
    free(particle_pos);

    tng_trajectory_destroy(&traj);

    return (0);
}

#endif

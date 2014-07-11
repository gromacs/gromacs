/* This code is part of the tng binary trajectory format.
 *
 * Written by Magnus Lundborg
 * Copyright (c) 2012-2014, The GROMACS development team.
 * Check out http://www.gromacs.org for more information.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */

#include "tng/tng_io.h"

#ifdef USE_STD_INTTYPES_H
#include <inttypes.h>
#endif

#include <stdlib.h>
#include <string.h>
#include "tng/version.h"

static tng_function_status tng_test_setup_molecules(tng_trajectory_t traj)
{
    tng_molecule_t molecule;
    tng_chain_t chain;
    tng_residue_t residue;
    tng_atom_t atom;
    int64_t cnt;

    tng_molecule_add(traj, "water", &molecule);
    tng_molecule_chain_add(traj, molecule, "W", &chain);
    tng_chain_residue_add(traj, chain, "WAT", &residue);
    if(tng_residue_atom_add(traj, residue, "O", "O", &atom) == TNG_CRITICAL)
    {
        return(TNG_CRITICAL);
    }
    if(tng_residue_atom_add(traj, residue, "HO1", "H", &atom) == TNG_CRITICAL)
    {
        return(TNG_CRITICAL);
    }
    if(tng_residue_atom_add(traj, residue, "HO2", "H", &atom) == TNG_CRITICAL)
    {
        return(TNG_CRITICAL);
    }
    tng_molecule_cnt_set(traj, molecule, 200);
    tng_molecule_cnt_get(traj, molecule, &cnt);
/*     printf("Created %"PRId64" %s molecules.\n", cnt, molecule->name); */

/*     traj->molecule_cnt_list[traj->n_molecules-1] = 5;
//     tng_molecule_name_set(traj, &traj->molecules[1], "ligand");
//     tng_molecule_name_set(traj, &traj->molecules[2], "water");
//     tng_molecule_name_set(traj, &traj->molecules[3], "dummy");
//     traj->molecules[0].id = 0;
//     traj->molecules[1].id = 1;
//     traj->molecules[2].id = 2;
//     traj->molecules[3].id = 3;

//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom1", "type1") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom2", "type1") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom3", "type1") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom4", "type2") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom5", "type2") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom6", "type2") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom7", "type3") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "C1", "C") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "O1", "O") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "H11", "H") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "H12", "H") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "H13", "H") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
*/
    return(TNG_SUCCESS);
}

static tng_function_status tng_test_read_and_write_file
                (tng_trajectory_t traj, const char hash_mode)
{
    tng_function_status stat;

    stat = tng_file_headers_read(traj, hash_mode);
    if(stat == TNG_CRITICAL)
    {
        return(stat);
    }
    stat = tng_file_headers_write(traj, hash_mode);
    if(stat == TNG_CRITICAL)
    {
        return(stat);
    }

    while(stat == TNG_SUCCESS)
    {
        stat = tng_frame_set_read_next(traj, hash_mode);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
        stat = tng_frame_set_write(traj, hash_mode);
    }

    return(stat);
}

static tng_function_status tng_test_write_and_read_traj(tng_trajectory_t *traj,
                                                        const char hash_mode)
{
    int i, j, k, nr, cnt;
    float *data, *molpos, *charges;
    int64_t mapping[300], n_particles, n_frames_per_frame_set, tot_n_mols;
    int64_t codec_id;
//     int64_t frame_nr;
    double box_shape[9];
    char atom_type[16], annotation[128];
    tng_function_status stat = TNG_SUCCESS;

    tng_medium_stride_length_set(*traj, 10);
    tng_long_stride_length_set(*traj, 100);

    tng_first_user_name_set(*traj, "User1");
    tng_first_program_name_set(*traj, "tng_testing");

    /* Create molecules */
    if(tng_test_setup_molecules(*traj) == TNG_CRITICAL)
    {
        return(TNG_CRITICAL);
    }

    /* Set the box shape */
    box_shape[1] = box_shape[2] = box_shape[3] = box_shape[5] = box_shape[6] =
    box_shape[7] = 0;
    box_shape[0] = 150.0;
    box_shape[4] = 145.5;
    box_shape[8] = 155.5;
    if(tng_data_block_add(*traj, TNG_TRAJ_BOX_SHAPE, "BOX SHAPE", TNG_DOUBLE_DATA,
                       TNG_NON_TRAJECTORY_BLOCK, 1, 9, 1, TNG_UNCOMPRESSED,
                       box_shape) == TNG_CRITICAL)
    {
        tng_trajectory_destroy(traj);
        printf("Cannot write trajectory box shape.\n");
        exit(1);
    }

    /* Set partial charges (treat the water as TIP3P. */
    tng_num_particles_get(*traj, &n_particles);
    charges = malloc(sizeof(float) * n_particles);
    for(i = 0; i < n_particles; i++)
    {
        stat = tng_atom_type_of_particle_nr_get(*traj, i, atom_type,
                                                sizeof(atom_type));
        if(stat == TNG_CRITICAL)
        {
            break;
        }
        if(atom_type[0] == 'O')
        {
            charges[i] = -0.834;
        }
        else if(atom_type[0] == 'H')
        {
            charges[i] = 0.417;
        }
    }
    if(stat == TNG_CRITICAL)
    {
        free(charges);
        printf("Failed setting partial charges.\n");
        return(TNG_CRITICAL);
    }

    stat = tng_particle_data_block_add(*traj, TNG_TRAJ_PARTIAL_CHARGES, "PARTIAL CHARGES",
                                       TNG_FLOAT_DATA, TNG_NON_TRAJECTORY_BLOCK,
                                       1, 1, 1, 0, n_particles,
                                       TNG_UNCOMPRESSED, charges);
    free(charges);
    if(stat != TNG_SUCCESS)
    {
        printf("Failed adding partial charges\n");
        return(TNG_CRITICAL);
    }


    /* Generate a custom annotation data block */
    strcpy(annotation, "This trajectory was generated from tng_io_testing. "
                       "It is not a real MD trajectory.");
    if(tng_data_block_add(*traj, TNG_TRAJ_GENERAL_COMMENTS, "COMMENTS", TNG_CHAR_DATA,
                          TNG_NON_TRAJECTORY_BLOCK, 1, 1, 1, TNG_UNCOMPRESSED,
                          annotation) != TNG_SUCCESS)
    {
        printf("Failed adding details annotation data block.\n");
        return(TNG_CRITICAL);
    }

    /* Write file headers (includes non trajectory data blocks */
    if(tng_file_headers_write(*traj, hash_mode) == TNG_CRITICAL)
    {
        printf("Cannot write file headers.\n");
    }


    tng_num_frames_per_frame_set_get(*traj, &n_frames_per_frame_set);

    data = malloc(sizeof(float) * n_particles *
                  n_frames_per_frame_set * 3);
    if(!data)
    {
        printf("Cannot allocate memory. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    tng_num_molecules_get(*traj, &tot_n_mols);
    molpos = malloc(sizeof(float) * tot_n_mols * 3);

    /* Set initial coordinates */
    for(i = 0; i < tot_n_mols; i++)
    {
        nr = i * 3;
        /* Somewhat random coordinates (between 0 and 100),
         * but not specifying a random seed */
        molpos[nr] = 100.0 * rand() / (RAND_MAX + 1.0);
        molpos[nr+1] = 100.0 * rand() / (RAND_MAX + 1.0);
        molpos[nr+2] = 100.0 * rand() / (RAND_MAX + 1.0);
    }

    /* Generate 200 frame sets - each with 100 frames (by default) */
    for(i = 0; i < 200; i++)
    {
        cnt = 0;
        if(i < 100)
        {
            codec_id = TNG_GZIP_COMPRESSION;
        }
        else
        {
            codec_id = TNG_TNG_COMPRESSION;
        }
        for(j = 0; j < n_frames_per_frame_set; j++)
        {
            for(k = 0; k < tot_n_mols; k++)
            {
                nr = k * 3;
                /* Move -1 to 1 */
                molpos[nr] += 2 * (rand() / (RAND_MAX + 1.0)) - 1;
                molpos[nr+1] += 2 * (rand() / (RAND_MAX + 1.0)) - 1;
                molpos[nr+2] += 2 * (rand() / (RAND_MAX + 1.0)) - 1;

                data[cnt++] = molpos[nr];
                data[cnt++] = molpos[nr + 1];
                data[cnt++] = molpos[nr + 2];
                data[cnt++] = molpos[nr] + 1;
                data[cnt++] = molpos[nr + 1] + 1;
                data[cnt++] = molpos[nr + 2] + 1;
                data[cnt++] = molpos[nr] - 1;
                data[cnt++] = molpos[nr + 1] - 1;
                data[cnt++] = molpos[nr + 2] - 1;
            }
        }
        if(tng_frame_set_new(*traj, i * n_frames_per_frame_set,
            n_frames_per_frame_set) != TNG_SUCCESS)
        {
            printf("Error creating frame set %d. %s: %d\n",
                   i, __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TNG_CRITICAL);
        }

        tng_frame_set_particle_mapping_free(*traj);

        /* Setup particle mapping. Use 4 different mapping blocks with arbitrary
         * mappings. */
        for(k=0; k<150; k++)
        {
            mapping[k]=k;
        }
        if(tng_particle_mapping_add(*traj, 0, 150, mapping) != TNG_SUCCESS)
        {
            printf("Error creating particle mapping. %s: %d\n",
                   __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TNG_CRITICAL);
        }
        for(k=0; k<150; k++)
        {
            mapping[k]=599-k;
        }
        if(tng_particle_mapping_add(*traj, 150, 150, mapping) != TNG_SUCCESS)
        {
            printf("Error creating particle mapping. %s: %d\n",
                   __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TNG_CRITICAL);
        }
        for(k=0; k<150; k++)
        {
            mapping[k]=k+150;
        }
        if(tng_particle_mapping_add(*traj, 300, 150, mapping) != TNG_SUCCESS)
        {
            printf("Error creating particle mapping. %s: %d\n",
                   __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TNG_CRITICAL);
        }
        for(k=0; k<150; k++)
        {
            mapping[k]=449-k;
        }
        if(tng_particle_mapping_add(*traj, 450, 150, mapping) != TNG_SUCCESS)
        {
            printf("Error creating particle mapping. %s: %d\n",
                   __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TNG_CRITICAL);
        }

        /* Add the positions in a data block */
        if(tng_particle_data_block_add(*traj, TNG_TRAJ_POSITIONS,
                                       "POSITIONS",
                                       TNG_FLOAT_DATA,
                                       TNG_TRAJECTORY_BLOCK,
                                       n_frames_per_frame_set, 3,
                                       1, 0, n_particles,
/*                                        TNG_UNCOMPRESSED, */
                                       codec_id,
                                       data) != TNG_SUCCESS)
        {
            printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TNG_CRITICAL);
        }
        /* Write the frame set */
        if(tng_frame_set_write(*traj, hash_mode) != TNG_SUCCESS)
        {
            printf("Error writing frame set. %s: %d\n", __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TNG_CRITICAL);
        }
    }

    /* Write two more frame sets one frame at a time */

    /* Make a new frame set - if always using the same mapping blocks
     * it is not necessary to explicitly add a new frame set - it will
     * be added automatically when adding data for a frame */
/*    if(tng_frame_set_new(*traj, i * n_frames_per_frame_set,
        n_frames_per_frame_set) != TNG_SUCCESS)
    {
        printf("Error creating frame set %d. %s: %d\n",
                i, __FILE__, __LINE__);
        free(molpos);
        free(data);
        return(TNG_CRITICAL);
    }

    frame_nr = i * n_frames_per_frame_set;

    for(k=0; k<300; k++)
    {
        mapping[k]=k;
    }
    *//* Just use two particle mapping blocks in this frame set *//*
    if(tng_particle_mapping_add(*traj, 0, 300, mapping) != TNG_SUCCESS)
    {
        printf("Error creating particle mapping. %s: %d\n",
                __FILE__, __LINE__);
        free(molpos);
        free(data);
        return(TNG_CRITICAL);
    }
    for(k=0; k<300; k++)
    {
        mapping[k]=599-k;
    }
    if(tng_particle_mapping_add(*traj, 300, 300, mapping) != TNG_SUCCESS)
    {
        printf("Error creating particle mapping. %s: %d\n",
                __FILE__, __LINE__);
        free(molpos);
        free(data);
        return(TNG_CRITICAL);
    }

    *//* Add the data block to the current frame set *//*
    if(tng_particle_data_block_add(*traj, TNG_TRAJ_POSITIONS,
                                    "POSITIONS",
                                    TNG_FLOAT_DATA,
                                    TNG_TRAJECTORY_BLOCK,
                                    n_frames_per_frame_set, 3,
                                    1, 0, n_particles,
                                    TNG_UNCOMPRESSED,
                                    0) != TNG_SUCCESS)
    {
        printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
        free(molpos);
        free(data);
        return(TNG_CRITICAL);
    }

    *//* Write the frame set to disk *//*
    if(tng_frame_set_write(*traj, hash_mode) != TNG_SUCCESS)
    {
        printf("Error writing frame set. %s: %d\n", __FILE__, __LINE__);
        free(molpos);
        free(data);
        return(TNG_CRITICAL);
    }

    *//* Write particle data to disk - one frame at a time *//*
    for(i = 0; i < n_frames_per_frame_set * 2; i++)
    {
        for(j = 0; j < 2; j++)
        {
            cnt = 0;
            for(k = 0; k < tot_n_mols/2; k++)
            {
                nr = k * 3;
                *//* Move -1 to 1 *//*
                molpos[nr] += 2 * (rand() / (RAND_MAX + 1.0)) - 1;
                molpos[nr+1] += 2 * (rand() / (RAND_MAX + 1.0)) - 1;
                molpos[nr+2] += 2 * (rand() / (RAND_MAX + 1.0)) - 1;

                data[cnt++] = molpos[nr];
                data[cnt++] = molpos[nr + 1];
                data[cnt++] = molpos[nr + 2];
                data[cnt++] = molpos[nr] + 1;
                data[cnt++] = molpos[nr + 1] + 1;
                data[cnt++] = molpos[nr + 2] + 1;
                data[cnt++] = molpos[nr] - 1;
                data[cnt++] = molpos[nr + 1] - 1;
                data[cnt++] = molpos[nr + 2] - 1;
            }
            if(tng_frame_particle_data_write(*traj, frame_nr + i,
                                          TNG_TRAJ_POSITIONS, j * 300, 300,
                                          data, hash_mode) != TNG_SUCCESS)
            {
                printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
                free(molpos);
                free(data);
                return(TNG_CRITICAL);
            }
        }
    }
*/
    free(molpos);
    free(data);

    tng_trajectory_destroy(traj);
    tng_trajectory_init(traj);
    tng_input_file_set(*traj, TNG_EXAMPLE_FILES_DIR "tng_test.tng");

    stat = tng_file_headers_read(*traj, hash_mode);

    while(stat == TNG_SUCCESS)
    {
        stat = tng_frame_set_read_next(*traj, hash_mode);
    }

    return(stat);
}

/* This test relies on knowing that the box shape is stored as double */
tng_function_status tng_test_get_box_data(tng_trajectory_t traj)
{
    int64_t n_frames, n_values_per_frame;
    union data_values **values = 0;
    char type;

    if(tng_data_get(traj, TNG_TRAJ_BOX_SHAPE, &values, &n_frames,
                    &n_values_per_frame, &type) != TNG_SUCCESS)
    {
        printf("Failed getting box shape. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

/*
//     int64_t i, j;
//     printf("Box shape:");
//     for(i=0; i<n_frames; i++)
//     {
//         for(j=0; j<n_values_per_frame; j++)
//         {
//             printf("\t%f", (values[i][j]).d);
//         }
//         printf("\n");
//     }
*/
    tng_data_values_free(traj, values, n_frames, n_values_per_frame, type);

    return(TNG_SUCCESS);
}

/* This test relies on knowing that the positions are stored as float
 * and that the data is not sparse (i.e. as many frames in the data
 * as in the frame set */
tng_function_status tng_test_get_positions_data(tng_trajectory_t traj,
                                                const char hash_mode)
{
    int64_t n_frames, n_particles, n_values_per_frame;
    union data_values ***values = 0;
    char type;

    if(tng_particle_data_get(traj, TNG_TRAJ_POSITIONS, &values, &n_frames,
                             &n_particles, &n_values_per_frame, &type) !=
       TNG_SUCCESS)
    {
        printf("Failed getting particle positions. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

/*
//     int64_t i, j, k;
//     struct tng_trajectory_frame_set *frame_set =
//     &traj->current_trajectory_frame_set;
//     for(i = 0; i<n_frames; i++)
//     {
//         printf("Frame %"PRId64"\n", frame_set->first_frame + i);
//         for(j = 0; j<n_particles; j++)
//         {
//             printf("Particle %"PRId64":", j);
//             for(k=0; k<n_values_per_frame; k++)
//             {
//                 printf("\t%f", (values[i][j][k]).f);
//             }
//             printf("\n");
//         }
//     }
*/
    tng_particle_data_values_free(traj, values, n_frames, n_particles,
                                  n_values_per_frame, type);

    values = 0;

    tng_particle_data_interval_get(traj, TNG_TRAJ_POSITIONS, 11000, 11499,
                                   hash_mode, &values, &n_particles,
                                   &n_values_per_frame, &type);

    /* Here the particle positions can be printed */

    tng_particle_data_values_free(traj, values, 500, n_particles,
                                  n_values_per_frame, type);

    return(TNG_SUCCESS);
}


tng_function_status tng_test_append(tng_trajectory_t traj, const char hash_mode)
{
    tng_function_status stat;

    stat = tng_util_trajectory_open(TNG_EXAMPLE_FILES_DIR "tng_test.tng", 'a', &traj);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    tng_last_user_name_set(traj, "User2");
    tng_last_program_name_set(traj, "tng_testing");
    tng_file_headers_write(traj, hash_mode);

    stat = tng_util_trajectory_close(&traj);

    return(stat);
}

int main()
{
    tng_trajectory_t traj;
    tng_function_status stat;
    char time_str[TNG_MAX_DATE_STR_LEN];
    char version_str[TNG_MAX_STR_LEN];
    char hash_mode = TNG_USE_HASH;

    tng_version(traj, version_str, TNG_MAX_STR_LEN);
    if(strncmp(TNG_VERSION, version_str, TNG_MAX_STR_LEN) == 0)
    {
        printf("Test version control: \t\t\t\tSucceeded.\n");
    }
    else
    {
        printf("Test version control: \t\t\t\tFailed.\n");
    }

    if(tng_trajectory_init(&traj) != TNG_SUCCESS)
    {
        tng_trajectory_destroy(&traj);
        printf("Test Init trajectory:\t\t\t\tFailed. %s: %d.\n",
               __FILE__, __LINE__);
        exit(1);
    }
    printf("Test Init trajectory:\t\t\t\tSucceeded.\n");

    tng_time_get_str(traj, time_str);

    printf("Creation time: %s\n", time_str);

    tng_input_file_set(traj, TNG_EXAMPLE_FILES_DIR "tng_example.tng");
    tng_output_file_set(traj, TNG_EXAMPLE_FILES_DIR "tng_example_out.tng");


    if(tng_test_read_and_write_file(traj, hash_mode) == TNG_CRITICAL)
    {
        printf("Test Read and write file:\t\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Read and write file:\t\t\tSucceeded.\n");
    }

    if(tng_test_get_box_data(traj) != TNG_SUCCESS)
    {
        printf("Test Get data:\t\t\t\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Get data:\t\t\t\t\tSucceeded.\n");
    }

    if(tng_trajectory_destroy(&traj) == TNG_CRITICAL ||
       tng_trajectory_init(&traj) == TNG_CRITICAL)
    {
        printf("Test Destroy and init trajectory:\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Destroy and init trajectory:\t\tSucceeded.\n");
    }


    tng_output_file_set(traj, TNG_EXAMPLE_FILES_DIR "tng_test.tng");

    if(tng_test_write_and_read_traj(&traj, hash_mode) == TNG_CRITICAL)
    {
        printf("Test Write and read file:\t\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Write and read file:\t\t\tSucceeded.\n");
    }

    if(tng_test_get_positions_data(traj, hash_mode) != TNG_SUCCESS)
    {
        printf("Test Get particle data:\t\t\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Get particle data:\t\t\t\tSucceeded.\n");
    }

    if(tng_trajectory_destroy(&traj) == TNG_CRITICAL)
    {
        printf("Test Destroy trajectory:\t\t\tFailed. %s: %d.\n",
               __FILE__, __LINE__);
        exit(1);
    }
    else
    {
        printf("Test Destroy trajectory:\t\t\tSucceeded.\n");
    }


    stat = tng_util_trajectory_open(TNG_EXAMPLE_FILES_DIR "tng_test.tng", 'r', &traj);

    if(stat != TNG_SUCCESS)
    {
        printf("Test Utility function open:\t\t\tFailed. %s: %d.\n",
               __FILE__, __LINE__);
        exit(1);
    }
    else
    {
        printf("Test Utility function open:\t\t\tSucceeded.\n");
    }

    stat = tng_util_trajectory_close(&traj);
    if(stat != TNG_SUCCESS)
    {
        printf("Test Utility function close:\t\t\tFailed. %s: %d.\n",
               __FILE__, __LINE__);
        exit(1);
    }
    else
    {
        printf("Test Utility function close:\t\t\tSucceeded.\n");
    }

    if(tng_test_append(traj, hash_mode) != TNG_SUCCESS)
    {
        printf("Test Append:\t\t\t\t\tFailed. %s: %d.\n",
               __FILE__, __LINE__);
        exit(1);
    }
    else
    {
        printf("Test Append:\t\t\t\t\tSucceeded.\n");
    }

    printf("Tests finished\n");

    exit(0);
}

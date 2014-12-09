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
#include <math.h>
#include "tng/version.h"

#define BOX_SHAPE_X 150.0
#define BOX_SHAPE_Y 145.5
#define BOX_SHAPE_Z 155.5
#define N_FRAME_SETS 100
#define MEDIUM_STRIDE_LEN 5
#define LONG_STRIDE_LEN 25
#define TIME_PER_FRAME 2e-15
#define COMPRESSION_PRECISION 1000
#define USER_NAME "USER 1"
#define PROGRAM_NAME "tng_testing"
#define COMPUTER_NAME "Unknown computer"
#define FORCEFIELD_NAME "No forcefield"

static tng_function_status tng_test_setup_molecules(tng_trajectory_t traj)
{
    tng_molecule_t molecule;
    tng_chain_t chain;
    tng_residue_t residue;
    tng_atom_t atom;
    tng_bond_t bond;
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
    tng_molecule_bond_add(traj, molecule, 0, 1, &bond);
    tng_molecule_bond_add(traj, molecule, 0, 2, &bond);
    tng_molecule_cnt_set(traj, molecule, 200);
    tng_molecule_cnt_get(traj, molecule, &cnt);

    if(cnt != 200)
    {
        return(TNG_CRITICAL);
    }
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

static tng_function_status tng_test_molecules(tng_trajectory_t traj)
{
    tng_molecule_t molecule, molecule_new;
    tng_chain_t chain;
    tng_residue_t residue;
    tng_atom_t atom;
    int64_t cnt, *bonds_to, *bonds_from;
    char var_atoms, str[TNG_MAX_STR_LEN];
    tng_function_status stat;

    stat = tng_num_molecule_types_get(traj, &cnt);
    if(stat != TNG_SUCCESS || cnt != 1)
    {
        printf("Molecule reading error. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    stat = tng_num_molecules_get(traj, &cnt);
    if(stat != TNG_SUCCESS || cnt != 200)
    {
        printf("Molecule reading error. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    stat = tng_num_particles_variable_get(traj, &var_atoms);
    if(stat != TNG_SUCCESS || var_atoms)
    {
        printf("Molecule reading error. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    stat = tng_molecule_of_index_get(traj, 0, &molecule);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot find molecule. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_molecule_find(traj, "water", -1, &molecule);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot find molecule. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat =tng_molecule_name_get(traj, molecule, str, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get molecule name. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_molecule_num_chains_get(traj, molecule, &cnt);
    if(stat != TNG_SUCCESS || cnt != 1)
    {
        printf("Cannot get number of chains in molecule. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_molecule_chain_of_index_get(traj, molecule, 0, &chain);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get chain in molecule. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_molecule_chain_find(traj, molecule, "W", -1, &chain);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get chain in molecule. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_molecule_num_residues_get(traj, molecule, &cnt);
    if(stat != TNG_SUCCESS || cnt != 1)
    {
        printf("Cannot get number of residues in molecule. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_molecule_residue_of_index_get(traj, molecule, 0, &residue);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get residue in molecule. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_molecule_num_atoms_get(traj, molecule, &cnt);
    if(stat != TNG_SUCCESS || cnt != 3)
    {
        printf("Cannot get number of atoms in molecule. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_molecule_atom_of_index_get(traj, molecule, 0, &atom);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get atom in molecule. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_molecule_atom_find(traj, molecule, "O", -1, &atom);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get atom in molecule. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }

    stat =tng_chain_name_get(traj, chain, str, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get name of chain. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_chain_num_residues_get(traj, chain, &cnt);
    if(stat != TNG_SUCCESS || cnt != 1)
    {
        printf("Cannot get number of residues in chain. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_chain_residue_of_index_get(traj, chain, 0, &residue);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get residue in chain. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_chain_residue_find(traj, chain, "WAT", -1, &residue);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get residue in chain. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }

    stat = tng_residue_name_get(traj, residue, str, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get name of residue. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_residue_num_atoms_get(traj, residue, &cnt);
    if(stat != TNG_SUCCESS || cnt != 3)
    {
        printf("Cannot get number of atoms in residue. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_residue_atom_of_index_get(traj, residue, 0, &atom);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get residue of atom. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }

    stat = tng_atom_name_get(traj, atom, str, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get name of atom. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_atom_type_get(traj, atom, str, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get atom type of atom. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }

    stat = tng_molecule_id_of_particle_nr_get(traj, 0, &cnt);
    if(stat != TNG_SUCCESS || cnt != 1)
    {
        printf("Cannot get molecule id of atom. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_residue_id_of_particle_nr_get(traj, 0, &cnt);
    if(stat != TNG_SUCCESS || cnt != 0)
    {
        printf("Cannot get residue id of atom. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_global_residue_id_of_particle_nr_get(traj, 599, &cnt);
    if(stat != TNG_SUCCESS || cnt != 199)
    {
        printf("Cannot get global residue id of atom. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_molecule_name_of_particle_nr_get(traj, 0, str, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get molecule name of atom. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_chain_name_of_particle_nr_get(traj, 0, str, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get chain name of atom. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_residue_name_of_particle_nr_get(traj, 0, str, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get residue name of atom. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_atom_name_of_particle_nr_get(traj, 0, str, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get atom name of atom. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }

    stat = tng_molecule_alloc(traj, &molecule_new);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot setup new molecule. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_molecule_name_set(traj, molecule_new, "TEST");
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot set name of new molecule. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_molecule_existing_add(traj, &molecule_new);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot add new molecule to molecule system. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }

    stat = tng_molsystem_bonds_get(traj, &cnt, &bonds_from, &bonds_to);
    if(stat != TNG_SUCCESS || cnt != 400)
    {
        printf("Cannot get bonds in molecule system. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }

    free(bonds_from);
    free(bonds_to);

    return(TNG_SUCCESS);
}

static tng_function_status tng_test_read_and_write_file
                (tng_trajectory_t traj, const char hash_mode)
{
    char file_name[TNG_MAX_STR_LEN];
    tng_function_status stat;

    stat = tng_input_file_get(traj, file_name, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Could not get name of input file. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_output_file_get(traj, file_name, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Could not get name of output file. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }

    stat = tng_file_headers_read(traj, hash_mode);
    if(stat != TNG_SUCCESS)
    {
        printf("Could not read headers. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_file_headers_write(traj, hash_mode);
    if(stat != TNG_SUCCESS)
    {
        printf("Could not write headers. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }

    while(stat == TNG_SUCCESS)
    {
        stat = tng_frame_set_read_next(traj, hash_mode);
        if(stat == TNG_CRITICAL)
        {
            printf("Could not read frame set. %s: %d\n",
                   __FILE__, __LINE__);
            return(stat);
        }
        if(stat == TNG_FAILURE)
        {
            return(TNG_SUCCESS);
        }
        stat = tng_frame_set_write(traj, hash_mode);
    }

    return(stat);
}

static tng_function_status tng_test_write_and_read_traj(tng_trajectory_t *traj,
                                                        const char hash_mode)
{
    int i, j, k, nr, cnt, dependency;
    float *data, *molpos, *charges;
    int64_t mapping[300], n_particles, n_frames_per_frame_set, tot_n_mols;
    int64_t codec_id;
    int64_t dist_exp = -9, temp_int, temp_int2;
//     int64_t frame_nr;
    double box_shape[9], temp_double;
    char atom_type[16], annotation[128];
    char temp_str[TNG_MAX_STR_LEN];
    tng_trajectory_frame_set_t frame_set;
    tng_function_status stat = TNG_SUCCESS;

    tng_medium_stride_length_set(*traj, MEDIUM_STRIDE_LEN);
    tng_long_stride_length_set(*traj, LONG_STRIDE_LEN);

    tng_first_user_name_set(*traj, USER_NAME);
    tng_first_program_name_set(*traj, PROGRAM_NAME);
    tng_first_computer_name_set(*traj, COMPUTER_NAME);
    tng_forcefield_name_set(*traj, FORCEFIELD_NAME);

    tng_compression_precision_set(*traj, COMPRESSION_PRECISION);

    tng_distance_unit_exponential_set(*traj, dist_exp);

    tng_time_per_frame_set(*traj, TIME_PER_FRAME);

    /* Create molecules */
    if(tng_test_setup_molecules(*traj) == TNG_CRITICAL)
    {
        return(TNG_CRITICAL);
    }

    /* Set the box shape */
    box_shape[1] = box_shape[2] = box_shape[3] = box_shape[5] = box_shape[6] =
    box_shape[7] = 0;
    box_shape[0] = BOX_SHAPE_X;
    box_shape[4] = BOX_SHAPE_Y;
    box_shape[8] = BOX_SHAPE_Z;
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
        printf("Failed setting partial charges. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    stat = tng_particle_data_block_add(*traj, TNG_TRAJ_PARTIAL_CHARGES, "PARTIAL CHARGES",
                                       TNG_FLOAT_DATA, TNG_NON_TRAJECTORY_BLOCK,
                                       1, 1, 1, 0, n_particles,
                                       TNG_UNCOMPRESSED, charges);
    free(charges);
    if(stat != TNG_SUCCESS)
    {
        printf("Failed adding partial charges. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }


    /* Generate a custom annotation data block */
    strcpy(annotation, "This trajectory was generated from tng_io_testing. "
                       "It is not a real MD trajectory.");
    if(tng_data_block_add(*traj, TNG_TRAJ_GENERAL_COMMENTS, "COMMENTS", TNG_CHAR_DATA,
                          TNG_NON_TRAJECTORY_BLOCK, 1, 1, 1, TNG_UNCOMPRESSED,
                          annotation) != TNG_SUCCESS)
    {
        printf("Failed adding details annotation data block. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* Write file headers (includes non trajectory data blocks */
    if(tng_file_headers_write(*traj, hash_mode) == TNG_CRITICAL)
    {
        printf("Cannot write file headers. %s: %d\n",
               __FILE__, __LINE__);
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

    /* Generate frame sets - each with 100 frames (by default) */
    for(i = 0; i < N_FRAME_SETS; i++)
    {
        cnt = 0;
        if(i < N_FRAME_SETS/2)
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
        if(tng_frame_set_with_time_new(*traj, i * n_frames_per_frame_set,
                                       n_frames_per_frame_set, 2e-15 * (i*n_frames_per_frame_set)) != TNG_SUCCESS)
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

    free(molpos);
    free(data);

    tng_trajectory_destroy(traj);
    tng_trajectory_init(traj);
    tng_input_file_set(*traj, TNG_EXAMPLE_FILES_DIR "tng_test.tng");

    stat = tng_file_headers_read(*traj, hash_mode);

    tng_first_user_name_get(*traj, temp_str, TNG_MAX_STR_LEN);
    if(strcmp(USER_NAME, temp_str) != 0)
    {
        printf("User name does not match when reading written file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    tng_first_program_name_get(*traj, temp_str, TNG_MAX_STR_LEN);
    if(strcmp(PROGRAM_NAME, temp_str) != 0)
    {
        printf("Program name does not match when reading written file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    tng_first_computer_name_get(*traj, temp_str, TNG_MAX_STR_LEN);
    if(strcmp(COMPUTER_NAME, temp_str) != 0)
    {
        printf("Computer name does not match when reading written file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    tng_forcefield_name_get(*traj, temp_str, TNG_MAX_STR_LEN);
    if(strcmp(FORCEFIELD_NAME, temp_str) != 0)
    {
        printf("Forcefield name does not match when reading written file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    tng_medium_stride_length_get(*traj, &temp_int);
    if(temp_int != MEDIUM_STRIDE_LEN)
    {
        printf("Stride length does not match when reading written file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    tng_long_stride_length_get(*traj, &temp_int);
    if(temp_int != LONG_STRIDE_LEN)
    {
        printf("Stride length does not match when reading written file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    tng_compression_precision_get(*traj, &temp_double);
    if(temp_double != COMPRESSION_PRECISION)
    {
        printf("Compression precision does not match when reading written file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    tng_distance_unit_exponential_get(*traj, &temp_int);
    if(temp_int != dist_exp)
    {
        printf("Distance unit exponential does not match when reading written file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    stat = tng_test_molecules(*traj);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    i = 0;
    while(stat == TNG_SUCCESS)
    {
        stat = tng_frame_set_read_next(*traj, hash_mode);
        tng_current_frame_set_get(*traj, &frame_set);
        tng_frame_set_prev_frame_set_file_pos_get(*traj, frame_set, &temp_int);
        tng_frame_set_next_frame_set_file_pos_get(*traj, frame_set, &temp_int2);
        if(i > 0)
        {
            if(temp_int == -1)
            {
                printf("File position of previous frame set not correct. %s: %d\n",
                       __FILE__, __LINE__);
                return(TNG_FAILURE);
            }
        }
        else if(temp_int != -1)
        {
            printf("File position of previous frame set not correct. %s: %d\n",
                    __FILE__, __LINE__);
            return(TNG_FAILURE);
        }
        if(i < N_FRAME_SETS -1)
        {
            if(temp_int2 == -1)
            {
                printf("File position of next frame set not correct. %s: %d\n",
                       __FILE__, __LINE__);
                return(TNG_FAILURE);
            }
        }
        else if(temp_int2 != -1)
        {
            printf("File position of previous next set not correct. %s: %d\n",
                    __FILE__, __LINE__);
            return(TNG_FAILURE);
        }
        i++;
    }
    if(stat == TNG_CRITICAL)
    {
        return(stat);
    }

    tng_time_per_frame_get(*traj, &temp_double);
    if(fabs(TIME_PER_FRAME - temp_double) > 0.000001)
    {
        printf("Time per frame does not match when reading written file. %s: %d\n",
               __FILE__, __LINE__);
        printf("Value: %e, expected value: %e\n", temp_double, TIME_PER_FRAME);
        return(TNG_FAILURE);
    }

    stat = tng_frame_set_nr_find(*traj, (int64_t)(0.30*N_FRAME_SETS));
    if(stat != TNG_SUCCESS)
    {
        printf("Could not find frame set %"PRId64". %s: %d\n", (int64_t)0.30*N_FRAME_SETS,
               __FILE__, __LINE__);
        return(stat);
    }

    stat = tng_frame_set_nr_find(*traj, (int64_t)(0.75*N_FRAME_SETS));
    if(stat != TNG_SUCCESS)
    {
        printf("Could not find frame set %"PRId64". %s: %d\n", (int64_t)0.75*N_FRAME_SETS,
               __FILE__, __LINE__);
        return(stat);
    }

    tng_current_frame_set_get(*traj, &frame_set);
    tng_frame_set_frame_range_get(*traj, frame_set, &temp_int, &temp_int2);
    if(temp_int !=  75 * n_frames_per_frame_set)
    {
        printf("Unexpected first frame in frame set. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    stat = tng_frame_set_read_current_only_data_from_block_id(*traj, hash_mode, TNG_TRAJ_POSITIONS);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot read positions in current frame set. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_frame_set_read_next_only_data_from_block_id(*traj, hash_mode, TNG_TRAJ_POSITIONS);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot read positions in next frame set. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_data_block_name_get(*traj, TNG_TRAJ_POSITIONS, temp_str, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS || strcmp("POSITIONS", temp_str) != 0)
    {
        printf("Cannot get name of data block or unexpected name. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_data_block_name_get(*traj, TNG_TRAJ_FORCES, temp_str, TNG_MAX_STR_LEN);
    if(stat != TNG_FAILURE)
    {
        printf("Trying to retrieve name of non-existent data block did not return failure. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_data_block_dependency_get(*traj, TNG_TRAJ_POSITIONS, &dependency);
    if(stat != TNG_SUCCESS || dependency != TNG_FRAME_DEPENDENT + TNG_PARTICLE_DEPENDENT)
    {
        printf("Cannot get dependency of data block or unexpected dependency. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_data_block_num_values_per_frame_get(*traj, TNG_TRAJ_POSITIONS, &temp_int);
    if(stat != TNG_SUCCESS || temp_int != 3)
    {
        printf("Cannot get number of values per frame of data block or unexpected value. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    stat = tng_data_get_stride_length(*traj, TNG_TRAJ_POSITIONS, 100, &temp_int);
    if(stat != TNG_SUCCESS || temp_int != 1)
    {
        printf("Cannot get stride length of data block or unexpected value. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    return(TNG_SUCCESS);
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

    /* The X dimension in the example file is 50 */
    if(fabs(values[0][0].d - 50) > 0.000001)
    {
        printf("Unexpected value in box shape. %s: %d\n", __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    tng_data_values_free(traj, values, n_frames, n_values_per_frame, type);

    return(TNG_SUCCESS);
}

/* This test relies on knowing that the positions are stored as float
 * and that the data is not sparse (i.e. as many frames in the data
 * as in the frame set */
tng_function_status tng_test_get_positions_data(tng_trajectory_t traj,
                                                const char hash_mode)
{
    int64_t i, j, k, n_frames, n_particles, n_values_per_frame;
    union data_values ***values = 0;
    char type;

    if(tng_particle_data_get(traj, TNG_TRAJ_POSITIONS, &values, &n_frames,
                             &n_particles, &n_values_per_frame, &type) !=
       TNG_SUCCESS)
    {
        printf("Failed getting particle positions. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(n_values_per_frame != 3)
    {
        printf("Number of values per frame does not match expected value. %s: %d\n",
               __FILE__, __LINE__);
        tng_particle_data_values_free(traj, values, n_frames, n_particles,
                                      n_values_per_frame, type);
        return(TNG_FAILURE);
    }

    for(i = 0; i < n_frames; i++)
    {
//         printf("%"PRId64"\n", i);
        for(j = 0; j < n_particles; j++)
        {
            for(k = 0; k < n_values_per_frame; k++)
            {
//                 printf("%f ", values[i][j][k].f);
                if(values[i][j][k].f < -500 || values[i][j][k].f > 500)
                {
                    printf("Coordinates not in range. %s: %d\n",
                           __FILE__, __LINE__);
                    tng_particle_data_values_free(traj, values, n_frames, n_particles,
                                                  n_values_per_frame, type);
                    return(TNG_FAILURE);
                }
            }
//             printf("\n");
        }
    }

    if(tng_particle_data_interval_get(traj, TNG_TRAJ_POSITIONS, 111000, 111499,
                                      hash_mode, &values, &n_particles,
                                      &n_values_per_frame, &type) == TNG_SUCCESS)
    {
        printf("Getting particle positions succeeded when it should have failed. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(tng_particle_data_interval_get(traj, TNG_TRAJ_POSITIONS, 1000, 1050,
                                      hash_mode, &values, &n_particles,
                                      &n_values_per_frame, &type) != TNG_SUCCESS)
    {
        printf("Failed getting particle positions. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    for(i = 0; i < 50; i++)
    {
//         printf("%"PRId64"\n", i);
        for(j = 0; j < n_particles; j++)
        {
            for(k = 0; k < n_values_per_frame; k++)
            {
//                 printf("%f ", values[i][j][k].f);
                if(values[i][j][k].f < -500 || values[i][j][k].f > 500)
                {
                    printf("Coordinates not in range. %s: %d\n",
                           __FILE__, __LINE__);
                    tng_particle_data_values_free(traj, values, n_frames, n_particles,
                                                  n_values_per_frame, type);
                    return(TNG_FAILURE);
                }
            }
//             printf("\n");
        }
    }

    tng_particle_data_values_free(traj, values, n_frames, n_particles,
                                  n_values_per_frame, type);

    return(TNG_SUCCESS);
}

tng_function_status tng_test_utility_functions(tng_trajectory_t traj, const char hash_mode)
{
    tng_function_status stat;
    int64_t n_particles, i, j, k, codec_id, n_frames, n_frames_per_frame_set;
    int64_t n_frames_to_read=30, stride_len, next_frame, n_blocks, *block_ids = 0;
    double time, multiplier;
    float *positions = 0;

    stat = tng_util_trajectory_open(TNG_EXAMPLE_FILES_DIR "tng_test.tng", 'r', &traj);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    stat = tng_util_time_of_frame_get(traj, 50, &time);
    if(stat != TNG_SUCCESS || fabs(time - 100e-13) > 0.000001)
    {
        printf("Unexpected time at frame 50. %s: %d\n", __FILE__, __LINE__);
        printf("Value: %e, expected value: %e\n", time, 100e-13);
        return(stat);
    }
    stat = tng_util_time_of_frame_get(traj, 100, &time);
    if(stat != TNG_SUCCESS || fabs(time - 200e-13) > 0.000001)
    {
        printf("Unexpected time at frame 100. %s: %d\n", __FILE__, __LINE__);
        printf("Value: %e, expected value: %e\n", time, 100e-13);
        return(stat);
    }

    tng_num_frames_per_frame_set_get(traj, &n_frames_per_frame_set);

    stat = tng_util_num_frames_with_data_of_block_id_get(traj, TNG_TRAJ_POSITIONS, &n_frames);
    if(stat != TNG_SUCCESS || n_frames != n_frames_per_frame_set * N_FRAME_SETS)
    {
        printf("Unexpected number of frames with positions data. %s: %d\n",
               __FILE__, __LINE__);
        printf("Value: %"PRId64", expected value: %"PRId64"\n", n_frames,
               n_frames_per_frame_set * N_FRAME_SETS);
        return(stat);
    }

    tng_num_frames_per_frame_set_get(traj, &n_frames_per_frame_set);

    stat = tng_util_num_frames_with_data_of_block_id_get(traj, TNG_TRAJ_POSITIONS, &n_frames);
    if(stat != TNG_SUCCESS || n_frames != n_frames_per_frame_set * N_FRAME_SETS)
    {
        return(stat);
    }

    tng_num_particles_get(traj, &n_particles);

    stat = tng_util_pos_read_range(traj, 1, n_frames_to_read, &positions, &stride_len);
    if(stat != TNG_SUCCESS)
    {
        if(positions)
        {
            free(positions);
        }
        return(stat);
    }

    for(i = 0; i < n_frames_to_read / stride_len; i++)
    {
        for(j = 0; j < n_particles; j++)
        {
            for(k = 0; k < 3; k++)
            {
                if(positions[i*n_particles + j*3 + k] < -500 || positions[i*n_particles + j*3 + k] > 500)
                {
                    printf("Coordinates not in range. %s: %d\n",
                           __FILE__, __LINE__);
                    free(positions);
                    return(TNG_FAILURE);
                }
            }
        }
    }

    free(positions);

    stat=tng_util_trajectory_next_frame_present_data_blocks_find(traj, n_frames_to_read,
                                                                 0, 0, &next_frame,
                                                                 &n_blocks, &block_ids);
    if(block_ids)
    {
        free(block_ids);
    }
    if(stat != TNG_SUCCESS || n_blocks != 1 || next_frame != n_frames_to_read + stride_len)
    {
        printf("Unexpected data blocks in next frame. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    stat = tng_util_frame_current_compression_get(traj, TNG_TRAJ_POSITIONS, &codec_id, &multiplier);
    if(stat != TNG_SUCCESS || codec_id != TNG_GZIP_COMPRESSION)
    {
        printf("Could not get compression. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    stat = tng_util_trajectory_close(&traj);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    return(TNG_SUCCESS);
}


tng_function_status tng_test_append(tng_trajectory_t traj, const char hash_mode)
{
    char str[TNG_MAX_STR_LEN];
    int64_t n_frames, n_particles, i;
    double time, *velocities;
    tng_function_status stat;

    stat = tng_util_trajectory_open(TNG_EXAMPLE_FILES_DIR "tng_test.tng", 'a', &traj);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot open trajectory. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }

    stat = tng_last_user_name_set(traj, USER_NAME);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot set last user name. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_last_user_name_get(traj, str, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get last user name. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_last_program_name_set(traj, PROGRAM_NAME);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot set last program name. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_last_program_name_get(traj, str, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get last program name. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_last_computer_name_set(traj, "Still " COMPUTER_NAME);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot set last computer name. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }
    stat = tng_last_computer_name_get(traj, str, TNG_MAX_STR_LEN);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot get last computer name. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }

    stat = tng_file_headers_write(traj, hash_mode);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot write file headers. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }

    tng_num_frames_get(traj, &n_frames);
    tng_frame_set_of_frame_find(traj, n_frames - 1);
    tng_util_time_of_frame_get(traj, n_frames - 1, &time);
    time += TIME_PER_FRAME;
    tng_num_particles_get(traj, &n_particles);

    velocities = malloc(sizeof(double) * n_particles * 3);
    if(!velocities)
    {
        printf("Cannot allocate memory. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    for(i = 0; i < n_particles * 3; i++)
    {
        velocities[i] = i;
    }

    stat = tng_util_vel_with_time_double_write(traj, n_frames, time, velocities);

    free(velocities);

    stat = tng_util_trajectory_close(&traj);

    return(stat);
}

tng_function_status tng_test_copy_container(tng_trajectory_t traj, const char hash_mode)
{
    tng_trajectory_t dest;
    tng_function_status stat;

    stat = tng_util_trajectory_open(TNG_EXAMPLE_FILES_DIR "tng_test.tng", 'r', &traj);
    if(stat != TNG_SUCCESS)
    {
        printf("Cannot open trajectory. %s: %d\n",
               __FILE__, __LINE__);
        return(stat);
    }

    stat = tng_trajectory_init_from_src(traj, &dest);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    stat = tng_molecule_system_copy(traj, dest);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    stat = tng_util_trajectory_close(&traj);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }
    stat = tng_util_trajectory_close(&dest);

    return(stat);
}

int main()
{
    tng_trajectory_t traj;
    char time_str[TNG_MAX_DATE_STR_LEN];
    char version_str[TNG_MAX_STR_LEN];
    char hash_mode = TNG_USE_HASH;

    tng_version(traj, version_str, TNG_MAX_STR_LEN);
    printf("Test version control:\t\t\t\t");
    if(strncmp(TNG_VERSION, version_str, TNG_MAX_STR_LEN) == 0)
    {
        printf("Succeeded.\n");
    }
    else
    {
        printf("Failed.\n");
    }

    printf("Test Init trajectory:\t\t\t\t");
    if(tng_trajectory_init(&traj) != TNG_SUCCESS)
    {
        tng_trajectory_destroy(&traj);
        printf("Failed. %s: %d.\n", __FILE__, __LINE__);
        exit(1);
    }
    printf("Succeeded.\n");

    tng_time_get_str(traj, time_str);

    printf("Creation time: %s\n", time_str);

    tng_input_file_set(traj, TNG_EXAMPLE_FILES_DIR "tng_example.tng");
    tng_output_file_set(traj, TNG_EXAMPLE_FILES_DIR "tng_example_out.tng");

    printf("Test Read and write file:\t\t\t");
    if(tng_test_read_and_write_file(traj, hash_mode) != TNG_SUCCESS)
    {
        printf("Failed. %s: %d\n", __FILE__, __LINE__);
    }
    else
    {
        printf("Succeeded.\n");
    }

    printf("Test Get data:\t\t\t\t\t");
    if(tng_test_get_box_data(traj) != TNG_SUCCESS)
    {
        printf("Failed. %s: %d\n", __FILE__, __LINE__);
    }
    else
    {
        printf("Succeeded.\n");
    }

    printf("Test Destroy and init trajectory:\t\t");
    if(tng_trajectory_destroy(&traj) != TNG_SUCCESS ||
       tng_trajectory_init(&traj) != TNG_SUCCESS)
    {
        printf("Failed. %s: %d\n", __FILE__, __LINE__);
    }
    else
    {
        printf("Succeeded.\n");
    }


    tng_output_file_set(traj, TNG_EXAMPLE_FILES_DIR "tng_test.tng");

    printf("Test Write and read file:\t\t\t");
    if(tng_test_write_and_read_traj(&traj, hash_mode) != TNG_SUCCESS)
    {
        printf("Failed. %s: %d\n", __FILE__, __LINE__);
    }
    else
    {
        printf("Succeeded.\n");
    }

    printf("Test Get particle data:\t\t\t\t");
    if(tng_test_get_positions_data(traj, hash_mode) != TNG_SUCCESS)
    {
        printf("Failed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Succeeded.\n");
    }

    printf("Test Destroy trajectory:\t\t\t");
    if(tng_trajectory_destroy(&traj) != TNG_SUCCESS)
    {
        printf("Failed. %s: %d.\n", __FILE__, __LINE__);
        exit(1);
    }
    else
    {
        printf("Succeeded.\n");
    }

    printf("Test Utility functions:\t\t\t\t");
    if(tng_test_utility_functions(traj, hash_mode) != TNG_SUCCESS)
    {
        printf("Failed. %s: %d.\n", __FILE__, __LINE__);
        exit(1);
    }
    else
    {
        printf("Succeeded.\n");
    }

    printf("Test Append:\t\t\t\t\t");
    if(tng_test_append(traj, hash_mode) != TNG_SUCCESS)
    {
        printf("Failed. %s: %d.\n", __FILE__, __LINE__);
    }
    else
    {
        printf("Succeeded.\n");
    }

    printf("Test Copy trajectory container:\t\t\t");
    if(tng_test_copy_container(traj, hash_mode) != TNG_SUCCESS)
    {
        printf("Failed. %s: %d.\n", __FILE__, __LINE__);
    }
    else
    {
        printf("Succeeded.\n");
    }

    printf("Tests finished\n");

    exit(0);
}

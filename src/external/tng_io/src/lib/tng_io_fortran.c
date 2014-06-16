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

#include "../../include/tng_io.h"

/* The following is for calling the library from fortran */

tng_function_status DECLSPECDLLEXPORT tng_trajectory_init_(tng_trajectory_t *tng_data_p)
{
    return(tng_trajectory_init(tng_data_p));
}

tng_function_status DECLSPECDLLEXPORT tng_trajectory_destroy_(tng_trajectory_t *tng_data_p)
{
    return(tng_trajectory_destroy(tng_data_p));
}

tng_function_status DECLSPECDLLEXPORT tng_trajectory_init_from_src_(tng_trajectory_t src,
                                                   tng_trajectory_t *dest_p)
{
    return(tng_trajectory_init_from_src(src, dest_p));
}

tng_function_status DECLSPECDLLEXPORT tng_input_file_get_(const tng_trajectory_t tng_data,
                                        char *file_name, const int max_len)
{
    return(tng_input_file_get(tng_data, file_name, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_input_file_set_(tng_trajectory_t tng_data,
                                        const char *file_name, int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, file_name, name_len);
    name[name_len] = 0;
    stat = tng_input_file_set(tng_data, name);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_output_file_get_(const tng_trajectory_t tng_data,
                                         char *file_name, const int max_len)
{
    return(tng_output_file_get(tng_data, file_name, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_output_file_set_(tng_trajectory_t tng_data,
                                         const char *file_name, int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, file_name, name_len);
    name[name_len] = 0;
    stat = tng_output_file_set(tng_data, name);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_first_program_name_get_(const tng_trajectory_t tng_data,
                                                char *name, const int max_len)
{
    return(tng_first_program_name_get(tng_data, name, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_output_file_endianness_get_
                (tng_trajectory_t tng_data, tng_file_endianness *endianness)
{
    return(tng_output_file_endianness_get(tng_data, endianness));
}

tng_function_status DECLSPECDLLEXPORT tng_output_file_endianness_set_
                (tng_trajectory_t tng_data, const tng_file_endianness *endianness)
{
    return(tng_output_file_endianness_set(tng_data, *endianness));
}

tng_function_status DECLSPECDLLEXPORT tng_first_program_name_set_(tng_trajectory_t tng_data,
                                                const char *new_name,
                                                int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, new_name, name_len);
    name[name_len] = 0;
    stat = tng_first_program_name_set(tng_data, name);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_last_program_name_get_(const tng_trajectory_t tng_data,
                                                char *name, const int max_len)
{
    return(tng_last_program_name_get(tng_data, name, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_last_program_name_set_(tng_trajectory_t tng_data,
                                               const char *new_name,
                                               int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, new_name, name_len);
    name[name_len] = 0;
    stat = tng_last_program_name_set(tng_data, name);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_first_user_name_get_(const tng_trajectory_t tng_data,
                                             char *name, const int max_len)
{
    return(tng_first_user_name_get(tng_data, name, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_first_user_name_set_(tng_trajectory_t tng_data,
                                             const char *new_name,
                                             int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, new_name, name_len);
    name[name_len] = 0;
    stat = tng_first_user_name_set(tng_data, name);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_last_user_name_get_(const tng_trajectory_t tng_data,
                                            char *name, const int max_len)
{
    return(tng_last_user_name_get(tng_data, name, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_last_user_name_set_(tng_trajectory_t tng_data,
                                            const char *new_name,
                                            int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, new_name, name_len);
    name[name_len] = 0;
    stat = tng_last_user_name_set(tng_data, name);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_first_computer_name_get_(const tng_trajectory_t tng_data,
                                                 char *name, const int max_len)
{
    return(tng_first_computer_name_get(tng_data, name, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_first_computer_name_set_(tng_trajectory_t tng_data,
                                                 const char *new_name,
                                                 int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, new_name, name_len);
    name[name_len] = 0;
    stat = tng_first_computer_name_set(tng_data, name);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_last_computer_name_get_(const tng_trajectory_t tng_data,
                                                 char *name, const int max_len)
{
    return(tng_last_computer_name_get(tng_data, name, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_last_computer_name_set_(tng_trajectory_t tng_data,
                                                const char *new_name,
                                                int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, new_name, name_len);
    name[name_len] = 0;
    stat = tng_last_computer_name_set(tng_data, name);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_first_signature_get_
                (const tng_trajectory_t tng_data,
                 char *signature, const int max_len)
{
    return(tng_first_signature_get(tng_data, signature, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_first_signature_set_(tng_trajectory_t tng_data,
                                             const char *signature,
                                             int sign_len)
{
    char *sign = malloc(sign_len + 1);
    tng_function_status stat;

    strncpy(sign, signature, sign_len);
    sign[sign_len] = 0;
    stat = tng_first_signature_set(tng_data, sign);
    free(sign);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_last_signature_get_
                (const tng_trajectory_t tng_data,
                 char *signature, const int max_len)
{
    return(tng_last_signature_get(tng_data, signature, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_last_signature_set_
                (tng_trajectory_t tng_data,
                 const char *signature,
                 int sign_len)
{
    char *sign = malloc(sign_len + 1);
    tng_function_status stat;

    strncpy(sign, signature, sign_len);
    sign[sign_len] = 0;
    stat = tng_last_signature_set(tng_data, sign);
    free(sign);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_forcefield_name_get_
                (const tng_trajectory_t tng_data,
                 char *name, const int max_len)
{
    return(tng_forcefield_name_get(tng_data, name, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_forcefield_name_set_
                (tng_trajectory_t tng_data,
                 const char *new_name,
                 int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, new_name, name_len);
    name[name_len] = 0;
    stat = tng_forcefield_name_set(tng_data, name);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_medium_stride_length_get_
                (const tng_trajectory_t tng_data,
                 int64_t *len)
{
    return(tng_medium_stride_length_get(tng_data, len));
}

tng_function_status DECLSPECDLLEXPORT tng_medium_stride_length_set_
                (tng_trajectory_t tng_data,
                 const int64_t *len)
{
    return(tng_medium_stride_length_set(tng_data, *len));
}

tng_function_status DECLSPECDLLEXPORT tng_long_stride_length_get_
                (const tng_trajectory_t tng_data,
                 int64_t *len)
{
    return(tng_long_stride_length_get(tng_data, len));
}

tng_function_status DECLSPECDLLEXPORT tng_long_stride_length_set_
                (tng_trajectory_t tng_data,
                 const int64_t *len)
{
    return(tng_long_stride_length_set(tng_data, *len));
}

tng_function_status DECLSPECDLLEXPORT tng_time_per_frame_get_
                (const tng_trajectory_t tng_data,
                 double *time)
{
    return(tng_time_per_frame_get(tng_data, time));
}

tng_function_status DECLSPECDLLEXPORT tng_time_per_frame_set_
                (tng_trajectory_t tng_data,
                 const double *time)
{
    return(tng_time_per_frame_set(tng_data, *time));
}

tng_function_status DECLSPECDLLEXPORT tng_input_file_len_get_
                (const tng_trajectory_t tng_data,
                 int64_t *len)
{
    return(tng_input_file_len_get(tng_data, len));
}

tng_function_status DECLSPECDLLEXPORT tng_num_frames_get_
                (const tng_trajectory_t tng_data,
                 int64_t *n)
{
    return(tng_num_frames_get(tng_data, n));
}

tng_function_status DECLSPECDLLEXPORT tng_num_particles_get_
                (const tng_trajectory_t tng_data,
                 int64_t *n)
{
    return(tng_num_particles_get(tng_data, n));
}

tng_function_status DECLSPECDLLEXPORT tng_num_molecules_get_
                (const tng_trajectory_t tng_data,
                 int64_t *n)
{
    return(tng_num_molecules_get(tng_data, n));
}

tng_function_status DECLSPECDLLEXPORT tng_distance_unit_exponential_get_
                (const tng_trajectory_t tng_data,
                 int64_t *exp)
{
    return(tng_distance_unit_exponential_get(tng_data, exp));
}

tng_function_status DECLSPECDLLEXPORT tng_distance_unit_exponential_set_
                (const tng_trajectory_t tng_data,
                 const int64_t *exp)
{
    return(tng_distance_unit_exponential_set(tng_data, *exp));
}

tng_function_status DECLSPECDLLEXPORT tng_num_frames_per_frame_set_get_
                (const tng_trajectory_t tng_data,
                 int64_t *n)
{
    return(tng_num_frames_per_frame_set_get(tng_data, n));
}

tng_function_status DECLSPECDLLEXPORT tng_num_frames_per_frame_set_set_
                (const tng_trajectory_t tng_data,
                 int64_t *n)
{
    return(tng_num_frames_per_frame_set_set(tng_data, *n));
}

tng_function_status DECLSPECDLLEXPORT tng_num_frame_sets_get_
                (const tng_trajectory_t tng_data,
                 int64_t *n)
{
    return(tng_num_frame_sets_get(tng_data, n));
}

tng_function_status DECLSPECDLLEXPORT tng_current_frame_set_get_
                (tng_trajectory_t tng_data,
                 tng_trajectory_frame_set_t *frame_set_p)
{
    return(tng_current_frame_set_get(tng_data, frame_set_p));
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_nr_find_(tng_trajectory_t tng_data,
                                                             const int64_t *nr)
{
    return(tng_frame_set_nr_find(tng_data, *nr));
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_of_frame_find_
                (tng_trajectory_t tng_data,
                 const int64_t *frame)
{
    return(tng_frame_set_of_frame_find(tng_data, *frame));
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_next_frame_set_file_pos_get_
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *pos)
{
    return(tng_frame_set_next_frame_set_file_pos_get(tng_data, frame_set, pos));
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_prev_frame_set_file_pos_get_
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *pos)
{
    return(tng_frame_set_prev_frame_set_file_pos_get(tng_data, frame_set, pos));
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_frame_range_get_
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *first_frame,
                 int64_t *last_frame)
{
    return(tng_frame_set_frame_range_get(tng_data, frame_set, first_frame,
                                         last_frame));
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_init_(const tng_trajectory_t tng_data,
                                       tng_molecule_t molecule)
{
    return(tng_molecule_init(tng_data, molecule));
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_destroy_
                (const tng_trajectory_t tng_data,
                 tng_molecule_t molecule)
{
    return(tng_molecule_destroy(tng_data, molecule));
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_add_(tng_trajectory_t tng_data,
                                                        const char *name,
                                                        tng_molecule_t *molecule,
                                                        int name_len)
{
    char *n = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(n, name, name_len);
    n[name_len] = 0;
    stat = tng_molecule_add(tng_data, n, molecule);
    free(n);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_name_set_(tng_trajectory_t tng_data,
                                                             tng_molecule_t molecule,
                                                             const char *new_name,
                                                             int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, new_name, name_len);
    name[name_len] = 0;
    stat = tng_molecule_name_set(tng_data, molecule, name);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_cnt_get_(tng_trajectory_t tng_data,
                                                            tng_molecule_t molecule,
                                                            int64_t *cnt)
{
    return(tng_molecule_cnt_get(tng_data, molecule, cnt));
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_cnt_set_(tng_trajectory_t tng_data,
                                                            tng_molecule_t molecule,
                                                            int64_t *cnt)
{
    return(tng_molecule_cnt_set(tng_data, molecule, *cnt));
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_find_(tng_trajectory_t tng_data,
                                                         const char *name,
                                                         int64_t nr,
                                                         tng_molecule_t *molecule,
                                                         int name_len)
{
    char *n = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(n, name, name_len);
    n[name_len] = 0;
    stat = tng_molecule_find(tng_data, n, nr, molecule);
    free(n);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_chain_find_(tng_trajectory_t tng_data,
                                                               tng_molecule_t molecule,
                                                               const char *name,
                                                               int64_t id,
                                                               tng_chain_t *chain,
                                                               int name_len)
{
    char *n = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(n, name, name_len);
    n[name_len] = 0;
    stat = tng_molecule_chain_find(tng_data, molecule, n, id, chain);
    free(n);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_chain_add_(tng_trajectory_t tng_data,
                                                              tng_molecule_t molecule,
                                                              const char *name,
                                                              tng_chain_t *chain,
                                                              int name_len)
{
    char *n = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(n, name, name_len);
    n[name_len] = 0;
    stat = tng_molecule_chain_add(tng_data, molecule, n, chain);
    free(n);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_chain_name_set_(tng_trajectory_t tng_data,
                                                          tng_chain_t chain,
                                                          const char *new_name,
                                                          int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, new_name, name_len);
    name[name_len] = 0;
    stat = tng_chain_name_set(tng_data, chain, name);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_chain_residue_add_(tng_trajectory_t tng_data,
                                                             tng_chain_t chain,
                                                             const char *name,
                                                             tng_residue_t *residue,
                                                             int name_len)
{
    char *n = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(n, name, name_len);
    n[name_len] = 0;
    stat = tng_chain_residue_add(tng_data, chain, n, residue);
    free(n);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_residue_name_set_(tng_trajectory_t tng_data,
                                                            tng_residue_t residue,
                                                            const char *new_name,
                                                            int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, new_name, name_len);
    name[name_len] = 0;
    stat = tng_residue_name_set(tng_data, residue, name);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_residue_atom_add_(tng_trajectory_t tng_data,
                                                            tng_residue_t residue,
                                                            const char *atom_name,
                                                            const char *atom_type,
                                                            tng_atom_t *atom,
                                                            int name_len,
                                                            int type_len)
{
    char *name = malloc(name_len + 1);
    char *type = malloc(type_len + 1);
    tng_function_status stat;

    strncpy(name, atom_name, name_len);
    strncpy(type, atom_type, type_len);
    name[name_len] = 0;
    type[type_len] = 0;
    stat = tng_residue_atom_add(tng_data, residue, name, type, atom);
    free(name);
    free(type);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_atom_name_set_(tng_trajectory_t tng_data,
                                                         tng_atom_t atom,
                                                         const char *new_name,
                                                         int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, new_name, name_len);
    name[name_len] = 0;
    stat = tng_atom_name_set(tng_data, atom, name);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_atom_type_set_(tng_trajectory_t tng_data,
                                                         tng_atom_t atom,
                                                         const char *new_type,
                                                         int type_len)
{
    char *type = malloc(type_len + 1);
    tng_function_status stat;

    strncpy(type, new_type, type_len);
    type[type_len] = 0;
    stat = tng_atom_type_set(tng_data, atom, type);
    free(type);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_name_of_particle_nr_get_
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    return(tng_molecule_name_of_particle_nr_get(tng_data, nr, name, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_chain_name_of_particle_nr_get_
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    return(tng_chain_name_of_particle_nr_get(tng_data, nr, name, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_residue_name_of_particle_nr_get_
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    return(tng_residue_name_of_particle_nr_get(tng_data, nr, name, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_residue_id_of_particle_nr_get_
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 int64_t *id)
{
    return(tng_residue_id_of_particle_nr_get(tng_data, nr, id));
}

tng_function_status DECLSPECDLLEXPORT tng_atom_name_of_particle_nr_get_
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    return(tng_atom_name_of_particle_nr_get(tng_data, nr, name, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_atom_type_of_particle_nr_get_
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *type,
                 int max_len)
{
    return(tng_atom_type_of_particle_nr_get(tng_data, nr, type, max_len));
}

tng_function_status DECLSPECDLLEXPORT tng_particle_mapping_add_
                (tng_trajectory_t tng_data,
                 const int64_t *first_particle_number,
                 const int64_t *n_particles,
                 const int64_t *mapping_table)
{
    return(tng_particle_mapping_add(tng_data, *first_particle_number,
                                    *n_particles, mapping_table));
}

tng_function_status DECLSPECDLLEXPORT tng_file_headers_read_(tng_trajectory_t tng_data,
                                                             const char *hash_mode)
{
    return(tng_file_headers_read(tng_data, *hash_mode));
}

tng_function_status DECLSPECDLLEXPORT tng_file_headers_write_
                (tng_trajectory_t tng_data,
                 const char *hash_mode)
{
    return(tng_file_headers_write(tng_data, *hash_mode));
}

tng_function_status DECLSPECDLLEXPORT tng_block_read_next_
                (tng_trajectory_t tng_data,
                 tng_gen_block_t block_data,
                 const char *hash_mode)
{
    return(tng_block_read_next(tng_data, block_data, *hash_mode));
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_read_next_
                (tng_trajectory_t tng_data,
                 const char *hash_mode)
{
    return(tng_frame_set_read_next(tng_data, *hash_mode));
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_write_(tng_trajectory_t tng_data,
                                                           const char *hash_mode)
{
    return(tng_frame_set_write(tng_data, *hash_mode));
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_new_(tng_trajectory_t tng_data,
                                                         const int64_t *first_frame,
                                                         const int64_t *n_frames)
{
    return(tng_frame_set_new(tng_data, *first_frame, *n_frames));
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_with_time_new_
                (tng_trajectory_t tng_data,
                 const int64_t *first_frame,
                 const int64_t *n_frames,
                 const double *first_frame_time)
{
    return(tng_frame_set_with_time_new(tng_data, *first_frame, *n_frames,
                             *first_frame_time));
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_first_frame_time_set_
                (tng_trajectory_t tng_data,
                 const double *first_frame_time)
{
    return(tng_frame_set_first_frame_time_set(tng_data, *first_frame_time));
}

tng_function_status DECLSPECDLLEXPORT tng_data_block_add_
                (tng_trajectory_t tng_data,
                 const int64_t *id,
                 const char *block_name,
                 const char *datatype,
                 const char *block_type_flag,
                 int64_t *n_frames,
                 const int64_t *n_values_per_frame,
                 const int64_t *stride_length,
                 const int64_t *codec_id,
                 void *new_data,
                 int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, block_name, name_len);
    name[name_len] = 0;
    stat = tng_data_block_add(tng_data, *id, name, *datatype, *block_type_flag,
                              *n_frames, *n_values_per_frame, *stride_length,
                              *codec_id, new_data);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_particle_data_block_add_
                (tng_trajectory_t tng_data,
                 const int64_t *id,
                 const char *block_name,
                 const char *datatype,
                 const char *block_type_flag,
                 int64_t *n_frames,
                 const int64_t *n_values_per_frame,
                 const int64_t *stride_length,
                 const int64_t *first_particle_number,
                 const int64_t *n_particles,
                 const int64_t *codec_id,
                 void *new_data,
                 int name_len)
{
    char *name = malloc(name_len + 1);
    tng_function_status stat;

    strncpy(name, block_name, name_len);
    name[name_len] = 0;
    stat = tng_particle_data_block_add(tng_data, *id, name, *datatype,
                                       *block_type_flag, *n_frames,
                                       *n_values_per_frame, *stride_length,
                                       *first_particle_number, *n_particles,
                                       *codec_id, new_data);
    free(name);
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_frame_data_write_
                (tng_trajectory_t tng_data,
                 const int64_t *frame_nr,
                 const int64_t *block_id,
                 const void *data,
                 const char *hash_mode)
{
    return(tng_frame_data_write(tng_data, *frame_nr, *block_id, data,
                                *hash_mode));
}

tng_function_status DECLSPECDLLEXPORT tng_frame_particle_data_write_
                (tng_trajectory_t tng_data,
                 const int64_t *frame_nr,
                 const int64_t *block_id,
                 const int64_t *val_first_particle,
                 const int64_t *val_n_particles,
                 const void *data,
                 const char *hash_mode)
{
    return(tng_frame_particle_data_write(tng_data, *frame_nr, *block_id,
                                         *val_first_particle, *val_n_particles,
                                         data, *hash_mode));
}

tng_function_status DECLSPECDLLEXPORT tng_data_values_free_
                (const tng_trajectory_t tng_data,
                 union data_values **values,
                 const int64_t *n_frames,
                 const int64_t *n_values_per_frame,
                 const char *type)
{
    return(tng_data_values_free(tng_data, values, *n_frames,
                                *n_values_per_frame, *type));
}

tng_function_status DECLSPECDLLEXPORT tng_particle_data_values_free_
                (const tng_trajectory_t tng_data,
                 union data_values ***values,
                 const int64_t *n_frames,
                 const int64_t *n_particles,
                 const int64_t *n_values_per_frame,
                 const char *type)
{
    return(tng_particle_data_values_free(tng_data, values, *n_frames, *n_particles,
                                         *n_values_per_frame, *type));
}

tng_function_status DECLSPECDLLEXPORT tng_data_get_
                (tng_trajectory_t tng_data,
                 const int64_t *block_id,
                 union data_values ***values,
                 int64_t *n_frames,
                 int64_t *n_values_per_frame,
                 char *type)
{
    return(tng_data_get(tng_data, *block_id, values, n_frames,
                        n_values_per_frame, type));
}

tng_function_status DECLSPECDLLEXPORT tng_data_interval_get_
                (tng_trajectory_t tng_data,
                 const int64_t *block_id,
                 const int64_t *start_frame_nr,
                 const int64_t *end_frame_nr,
                 const char *hash_mode,
                 union data_values ***values,
                 int64_t *n_values_per_frame,
                 char *type)
{
    return(tng_data_interval_get(tng_data, *block_id, *start_frame_nr,
                                 *end_frame_nr, *hash_mode, values,
                                 n_values_per_frame, type));
}

tng_function_status DECLSPECDLLEXPORT tng_particle_data_get_
                (tng_trajectory_t tng_data,
                 const int64_t *block_id,
                 union data_values ****values,
                 int64_t *n_frames,
                 int64_t *n_particles,
                 int64_t *n_values_per_frame,
                 char *type)
{
    return(tng_particle_data_get(tng_data, *block_id, values, n_frames,
                                 n_particles, n_values_per_frame, type));
}

tng_function_status DECLSPECDLLEXPORT tng_particle_data_interval_get_
                (tng_trajectory_t tng_data,
                 const int64_t *block_id,
                 const int64_t *start_frame_nr,
                 const int64_t *end_frame_nr,
                 const char *hash_mode,
                 union data_values ****values,
                 int64_t *n_particles,
                 int64_t *n_values_per_frame,
                 char *type)
{
    return(tng_particle_data_interval_get(tng_data, *block_id, *start_frame_nr,
                                          *end_frame_nr, *hash_mode, values,
                                          n_particles, n_values_per_frame,
                                          type));
}

tng_function_status DECLSPECDLLEXPORT tng_time_get_str_
                (const tng_trajectory_t tng_data,
                 char *time, int64_t str_len)
{
    return(tng_time_get_str(tng_data, time));
}

tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_open_
                (const char *filename, const char *mode,
                 tng_trajectory_t *tng_data_p)
{
    return(tng_util_trajectory_open(filename, *mode, tng_data_p));
}

tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_close_
                (tng_trajectory_t *tng_data_p)
{
    return(tng_util_trajectory_close(tng_data_p));
}

tng_function_status DECLSPECDLLEXPORT tng_util_pos_read_
                (tng_trajectory_t tng_data,
                 float **positions,
                 int64_t *stride_length)
{
    return(tng_util_pos_read(tng_data, positions, stride_length));
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_read_
                (tng_trajectory_t tng_data,
                 float **velocities,
                 int64_t *stride_length)
{
    return(tng_util_vel_read(tng_data, velocities, stride_length));
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_read_
                (tng_trajectory_t tng_data,
                 float **forces,
                 int64_t *stride_length)
{
    return(tng_util_force_read(tng_data, forces, stride_length));
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_read_
                (tng_trajectory_t tng_data,
                 float **box_shape,
                 int64_t *stride_length)
{
    return(tng_util_box_shape_read(tng_data, box_shape, stride_length));
}

tng_function_status DECLSPECDLLEXPORT tng_util_pos_read_range_
                (tng_trajectory_t tng_data,
                 const int64_t *first_frame,
                 const int64_t *last_frame,
                 float **positions,
                 int64_t *stride_length)
{
    return(tng_util_pos_read_range(tng_data, *first_frame, *last_frame,
                                         positions, stride_length));
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_read_range_
                (tng_trajectory_t tng_data,
                 const int64_t *first_frame,
                 const int64_t *last_frame,
                 float **velocities,
                 int64_t *stride_length)
{
    return(tng_util_vel_read_range(tng_data, *first_frame, *last_frame,
                                         velocities, stride_length));
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_read_range_
                (tng_trajectory_t tng_data,
                 const int64_t *first_frame,
                 const int64_t *last_frame,
                 float **forces,
                 int64_t *stride_length)
{
    return(tng_util_force_read_range(tng_data, *first_frame, *last_frame,
                                         forces, stride_length));
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_read_range_
                (tng_trajectory_t tng_data,
                 const int64_t *first_frame,
                 const int64_t *last_frame,
                 float **box_shape,
                 int64_t *stride_length)
{
    return(tng_util_box_shape_read_range(tng_data, *first_frame, *last_frame,
                                         box_shape, stride_length));
}

tng_function_status DECLSPECDLLEXPORT tng_util_generic_write_frequency_set_
                (tng_trajectory_t tng_data,
                 const int64_t *f,
                 const int64_t *n_values_per_frame,
                 const int64_t *block_id,
                 const char *block_name,
                 const char *particle_dependency,
                 const char *compression)
{
    return(tng_util_generic_write_frequency_set(tng_data, *f,
                                                *n_values_per_frame, *block_id,
                                                block_name,
                                                *particle_dependency,
                                                *compression));
}

tng_function_status DECLSPECDLLEXPORT tng_util_pos_write_frequency_set_
                (tng_trajectory_t tng_data,
                 const int64_t *f)
{
    return(tng_util_pos_write_frequency_set(tng_data, *f));
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_write_frequency_set_
                (tng_trajectory_t tng_data,
                 const int64_t *f)
{
    return(tng_util_vel_write_frequency_set(tng_data, *f));
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_write_frequency_set_
                (tng_trajectory_t tng_data,
                 const int64_t *f)
{
    return(tng_util_force_write_frequency_set(tng_data, *f));
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_write_frequency_set_
                (tng_trajectory_t tng_data,
                 const int64_t *f)
{
    return(tng_util_box_shape_write_frequency_set(tng_data, *f));
}

tng_function_status DECLSPECDLLEXPORT tng_util_generic_write_
                (tng_trajectory_t tng_data,
                 const int64_t *frame_nr,
                 const float *values,
                 const int64_t *n_values_per_frame,
                 const int64_t *block_id,
                 const char *block_name,
                 const char *particle_dependency,
                 const char *compression)
{
    return(tng_util_generic_write(tng_data, *frame_nr, values,
                                  *n_values_per_frame, *block_id, block_name,
                                  *particle_dependency, *compression));
}

tng_function_status DECLSPECDLLEXPORT tng_util_pos_write_
                (tng_trajectory_t tng_data,
                 const int64_t *frame_nr,
                 const float *positions)
{
    return(tng_util_vel_write(tng_data, *frame_nr, positions));
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_write_
                (tng_trajectory_t tng_data,
                 const int64_t *frame_nr,
                 const float *velocities)
{
    return(tng_util_vel_write(tng_data, *frame_nr, velocities));
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_write_
                (tng_trajectory_t tng_data,
                 const int64_t *frame_nr,
                 const float *forces)
{
    return(tng_util_force_write(tng_data, *frame_nr, forces));
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_write_
                (tng_trajectory_t tng_data,
                 const int64_t *frame_nr,
                 const float *box_shape)
{
    return(tng_util_box_shape_write(tng_data, *frame_nr, box_shape));
}

tng_function_status DECLSPECDLLEXPORT tng_util_pos_with_time_write_
                (tng_trajectory_t tng_data,
                 const int64_t *frame_nr,
                 const int64_t *time,
                 const float *positions)
{
    return(tng_util_pos_with_time_write(tng_data, *frame_nr, *time,
                                        positions));
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_with_time_write_
                (tng_trajectory_t tng_data,
                 const int64_t *frame_nr,
                 const int64_t *time,
                 const float *velocities)
{
    return(tng_util_vel_with_time_write(tng_data, *frame_nr, *time,
                                        velocities));
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_with_time_write_
                (tng_trajectory_t tng_data,
                 const int64_t *frame_nr,
                 const int64_t *time,
                 const float *forces)
{
    return(tng_util_force_with_time_write(tng_data, *frame_nr, *time,
                                          forces));
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_with_time_write_
                (tng_trajectory_t tng_data,
                 const int64_t *frame_nr,
                 const int64_t *time,
                 const float *box_shape)
{
    return(tng_util_box_shape_with_time_write(tng_data, *frame_nr, *time,
                                              box_shape));
}

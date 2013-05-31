/* This code is part of the tng binary trajectory format.
 *
 *                      VERSION 1.0
 *
 * Written by Magnus Lundborg
 * Copyright (c) 2012, The GROMACS development team.
 * Check out http://www.gromacs.org for more information.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 */

#ifdef USE_STD_INTTYPES_H
#include <inttypes.h>
#endif

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef USE_ZLIB
#include <zlib.h>
#endif

#include "tng_io.h"
#include "md5.h"


struct tng_bond {
    /** One of the atoms of the bond */
    int64_t from_atom_id;
    /** The other atom of the bond */
    int64_t to_atom_id;
};

struct tng_atom {
    /** The residue containing this atom */
    tng_residue_t residue;
    /** A unique (per molecule) ID number of the atom */
    int64_t id;
    /** The atom_type (depending on the forcefield) */
    char *atom_type;
    /** The name of the atom */
    char *name;
};

struct tng_residue {
    /** The chain containing this residue */
    tng_chain_t chain;
    /** A unique (per chain) ID number of the residue */
    int64_t id;
    /** The name of the residue */
    char *name;
    /** The number of atoms in the residue */
    int64_t n_atoms;
    /** A list of atoms in the residue */
    int64_t atoms_offset;
};

struct tng_chain {
    /** The molecule containing this chain */
    tng_molecule_t molecule;
    /** A unique (per molecule) ID number of the chain */
    int64_t id;
    /** The name of the chain */
    char *name;
    /** The number of residues in the chain */
    int64_t n_residues;
    /** A list of residues in the chain */
    tng_residue_t residues;
};

struct tng_molecule {
    /** A unique ID number of the molecule */
    int64_t id;
    /** Quaternary structure of the molecule.
     *  1 => monomeric
     *  2 => dimeric
     *  3 => trimeric
     *  etc */
    int64_t quaternary_str;
    /** The number of chains in the molecule */
    int64_t n_chains;
    /** The number of residues in the molecule */
    int64_t n_residues;
    /** The number of atoms in the molecule */
    int64_t n_atoms;
    /** The number of bonds in the molecule. If the bonds are not specified this
     * value can be 0. */
    int64_t n_bonds;
    /** The name of the molecule */
    char *name;
    /** A list of chains in the molecule */
    tng_chain_t chains;
    /** A list of residues in the molecule */
    tng_residue_t residues;
    /** A list of the atoms in the molecule */
    tng_atom_t atoms;
    /** A list of the bonds in the molecule */
    tng_bond_t bonds;
};

struct tng_gen_block {
    /** The size of the block header in bytes */
    int64_t header_contents_size;
    /** The size of the block contents in bytes */
    int64_t block_contents_size;
    /** The ID of the block to determine its type */
    int64_t id;
    /** The MD5 hash of the block to verify integrity */
    char hash[TNG_HASH_LEN];
    /** The name of the block */
    char *name;
    /** The library version used to write the block */
    int64_t block_version;
    /** The full block header contents */
    char *header_contents;
    /** The full block contents */
    char *block_contents;
};

struct tng_particle_mapping {
    /** The index number of the first particle in this mapping block */
    int64_t num_first_particle;
    /** The number of particles list in this mapping block */
    int64_t n_particles;
    /** the mapping of index numbers to the real particle numbers in the
     * trajectory. real_particle_numbers[0] is the real particle number
     * (as it is numbered in the molecular system) of the first particle
     * in the data blocks covered by this particle mapping block */
    int64_t *real_particle_numbers;
};

struct tng_trajectory_frame_set {
    /** The number of different particle mapping blocks present. */
    int64_t n_mapping_blocks;
    /** The atom mappings of this frame set */
    struct tng_particle_mapping *mappings;
    /** The first frame of this frame set */
    int64_t first_frame;
    /** The number of frames in this frame set */
    int64_t n_frames;
    /** The number of written frames in this frame set (used when writing one
     * frame at a time. */
    int64_t n_written_frames;

    /** A list of the number of each molecule type - only used when using
     * variable number of atoms */
    int64_t *molecule_cnt_list;
    /** The number of particles/atoms - only used when using variable number
     * of atoms */
    int64_t n_particles;
    /** The file position of the next frame set */
    int64_t next_frame_set_file_pos;
    /** The file position of the previous frame set */
    int64_t prev_frame_set_file_pos;
    /** The file position of the frame set one long stride step ahead */
    int64_t medium_stride_next_frame_set_file_pos;
    /** The file position of the frame set one long stride step behind */
    int64_t medium_stride_prev_frame_set_file_pos;
    /** The file position of the frame set one long stride step ahead */
    int64_t long_stride_next_frame_set_file_pos;
    /** The file position of the frame set one long stride step behind */
    int64_t long_stride_prev_frame_set_file_pos;

    /* The data blocks in a frame set are trajectory data blocks */
    /** The number of trajectory data blocks of particle dependent data */
    int n_particle_data_blocks;
    /** A list of data blocks containing particle dependent data */
    struct tng_particle_data *tr_particle_data;
    /** The number of trajectory data blocks independent of particles */
    int n_data_blocks;
    /** A list of data blocks containing particle indepdendent data */
    struct tng_non_particle_data *tr_data;
};

/* FIXME: Should there be a pointer to a tng_gen_block from each data block? */
struct tng_particle_data {
    /** The block ID of the data block containing this particle data.
     *  This is used to determine the kind of data that is stored */
    int64_t block_id;
    /** The name of the data block. This is used to determine the kind of
     *  data that is stored */
    char *block_name;
    /** The type of data stored. */
    tng_data_type datatype;
    /** The frame number of the first data value */
    int64_t first_frame_with_data;
    /** The number of frames in this frame set */
    int64_t n_frames;
    /** The number of values stored per frame */
    int64_t n_values_per_frame;
    /** The number of frames between each data point - e.g. when
     *  storing sparse data. */
    int64_t stride_length;
    /** ID of the CODEC used for compression 0 == no compression. */
    int64_t codec_id;
    /** The multiplier used for getting integer values for compression */
    double compression_multiplier;
    /** A 3-dimensional array of values, sized
     *  n_frames * n_particles * n_values_per_frame */
    union data_values ***values;
};

struct tng_non_particle_data {
    /** The ID of the data block */
    int64_t block_id;
    /** The name of the data block. This is used to determine the kind of
     *  data that is stored */
    char *block_name;
    /** The type of data stored. */
    tng_data_type datatype;
    /** The first frame number of the first data value */
    int64_t first_frame_with_data;
    /** The number of frames in this data block */
    int64_t n_frames;
    /** The number of values stored per frame */
    int64_t n_values_per_frame;
    /** The number of frames between each data value, e.g. if storing data
     *  that is not saved every frame. */
    int64_t stride_length;
    /** ID of the CODEC used for compression. 0 == no compression. */
    int64_t codec_id;
    /** Compressed data is stored as integers. This compression multiplier is
     *  the multiplication factor to convert from integer to float/double */
    double compression_multiplier;
    /** A 2-dimensional array of values, sized n_frames * n_values_per_frame */
    union data_values **values;
};



struct tng_trajectory {
    /** The path of the input trajectory file */
    char *input_file_path;
    /** A handle to the input file */
    FILE *input_file;
    /** The length of the input file */
    long int input_file_len;
    /** The path of the output trajectory file */
    char *output_file_path;
    /** A handle to the output file */
    FILE *output_file;
    /** Function to swap 32 bit values to and from the endianness of the
     * input file */
    tng_function_status (*input_endianness_swap_func_32)(const tng_trajectory_t, int32_t *);
    /** Function to swap 64 bit values to and from the endianness of the
     * input file */
    tng_function_status (*input_endianness_swap_func_64)(const tng_trajectory_t, int64_t *);
    /** Function to swap 32 bit values to and from the endianness of the
     * input file */
    tng_function_status (*output_endianness_swap_func_32)(const tng_trajectory_t, int32_t *);
    /** Function to swap 64 bit values to and from the endianness of the
     * input file */
    tng_function_status (*output_endianness_swap_func_64)(const tng_trajectory_t, int64_t *);
    /** The endianness of 32 bit values of the current computer */
    tng_endianness_32 endianness_32;
    /** The endianness of 64 bit values of the current computer */
    tng_endianness_64 endianness_64;

    /** The name of the program producing this trajectory */
    char *first_program_name;
    /** The forcefield used in the simulations */
    char *forcefield_name;
    /** The name of the user running the simulations */
    char *first_user_name;
    /** The name of the computer on which the simulations were performed */
    char *first_computer_name;
    /** The PGP signature of the user creating the file. */
    char *first_pgp_signature;
    /** The name of the program used when making last modifications to the
     *  file */
    char *last_program_name;
    /** The name of the user making the last modifications to the file */
    char *last_user_name;
    /** The name of the computer on which the last modifications were made */
    char *last_computer_name;
    /** The PGP signature of the user making the last modifications to the
     *  file. */
    char *last_pgp_signature;
    /** The time (n seconds since 1970) when the file was created */
    int64_t time;

    /** A flag indicating if the number of atoms can vary throughout the
     *  simulation, e.g. using a grand canonical ensemble */
    char var_num_atoms_flag;
    /** The number of frames in a frame set. It is allowed to have frame sets
     *  with fewer frames, but this will help searching for specific frames */
    int64_t frame_set_n_frames;
    /** The number of frame sets in a medium stride step */
    int64_t medium_stride_length;
    /** The number of frame sets in a long stride step */
    int64_t long_stride_length;

    /** The number of different kinds of molecules in the trajectory */
    int64_t n_molecules;
    /** A list of molecules in the trajectory */
    tng_molecule_t molecules;
    /** A list of the count of each molecule - if using variable number of
     *  particles this will be specified in each frame set */
    int64_t *molecule_cnt_list;
    /** The total number of particles/atoms. If using variable number of
     *  particles this will be specified in each frame set */
    int64_t n_particles;

     /** The pos in the src file of the first frame set */
    int64_t first_trajectory_frame_set_input_file_pos;
    /** The pos in the dest file of the first frame set */
    int64_t first_trajectory_frame_set_output_file_pos;
    /** The pos in the src file of the last frame set */
    int64_t last_trajectory_frame_set_input_file_pos;
    /** The pos in the dest file of the last frame set */
    int64_t last_trajectory_frame_set_output_file_pos;
    /** The currently active frame set */
    struct tng_trajectory_frame_set current_trajectory_frame_set;
    /** The pos in the src file of the current frame set */
    long int current_trajectory_frame_set_input_file_pos;
    /** The pos in the dest file of the current frame set */
    long int current_trajectory_frame_set_output_file_pos;
    /** The number of frame sets in the trajectory N.B. Not saved in file and
     *  cannot be trusted to be up-to-date */
    int64_t n_trajectory_frame_sets;

    /** The number of trajectory blocks in the file */
    int64_t n_trajectory_blocks;

    /* These data blocks are non-trajectory data blocks */
    /** The number of non-frame dependent particle dependent data blocks */
    int n_particle_data_blocks;
    /** A list of data blocks containing particle dependent data */
    struct tng_particle_data *non_tr_particle_data;

    /** The number of frame and particle independent data blocks */
    int n_data_blocks;
    /** A list of frame and particle indepdendent data blocks */
    struct tng_non_particle_data *non_tr_data;
};




/** This function swaps the byte order of a 32 bit numerical variable
 * to big endian.
 * It does not only work with integer, but e.g. floats need casting.
 * If the byte order is already big endian no change is needed.
 * @param tng_data is a trajectory data container.
 * @param v is a pointer to a 32 bit numerical value (float or integer).
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the current
 * byte order is not recognised.
 */
static tng_function_status tng_swap_byte_order_big_endian_32
                (const tng_trajectory_t tng_data, int32_t *v)
{
    switch(tng_data->endianness_32)
    {
    case TNG_LITTLE_ENDIAN_32: // Byte order is reversed.
        *v = ((*v & 0xFF000000) >> 24) | // Move first byte to end
             ((*v & 0x00FF0000) >> 8) |  // Move 2nd byte to pos 3
             ((*v & 0x0000FF00) << 8) |  // Move 3rd byte to pos 2
             ((*v & 0x000000FF) << 24);  // Move last byte to first

        return(TNG_SUCCESS);

    case TNG_BYTE_PAIR_SWAP_32: // byte pair swap
        *v = ((*v & 0xFFFF0000) >> 16) |
             ((*v & 0x0000FFFF) << 16);

        return(TNG_SUCCESS);

    case TNG_BIG_ENDIAN_32: // Already correct
        return(TNG_SUCCESS);

    default:
        return(TNG_FAILURE);
    }
}

/** This function swaps the byte order of a 64 bit numerical variable
 * to big endian.
 * It does not only work with integer, but e.g. floats need casting.
 * The byte order swapping routine can convert four different byte
 * orders to big endian.
 * If the byte order is already big endian no change is needed.
 * @param tng_data is a trajectory data container.
 * @param v is a pointer to a 64 bit numerical value (double or integer).
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the current
 * byte order is not recognised.
 */
static tng_function_status tng_swap_byte_order_big_endian_64
                (const tng_trajectory_t tng_data, int64_t *v)
{
    switch(tng_data->endianness_64)
    {
    case TNG_LITTLE_ENDIAN_64: // Byte order is reversed.
        *v = ((*v & 0xFF00000000000000) >> 56) | // Move first byte to end
             ((*v & 0x00FF000000000000) >> 40) | // Move 2nd byte to pos 7
             ((*v & 0x0000FF0000000000) >> 24) | // Move 3rd byte to pos 6
             ((*v & 0x000000FF00000000) >> 8 ) | // Move 4th byte to pos 5
             ((*v & 0x00000000FF000000) << 8 ) | // Move 5th byte to pos 4
             ((*v & 0x0000000000FF0000) << 24) | // Move 6th byte to pos 3
             ((*v & 0x000000000000FF00) << 40) | // Move 7th byte to pos 2
             ((*v & 0x00000000000000FF) << 56);  // Move last byte to first

        return(TNG_SUCCESS);

    case TNG_QUAD_SWAP_64: // Byte quad swap
        *v = ((*v & 0xFFFFFFFF00000000) >> 32) |
             ((*v & 0x00000000FFFFFFFF) << 32);

        return(TNG_SUCCESS);

    case TNG_BYTE_PAIR_SWAP_64: // Byte pair swap
        *v = ((*v & 0xFFFF0000FFFF0000) >> 16) |
             ((*v & 0x0000FFFF0000FFFF) << 16);

        return(TNG_SUCCESS);

    case TNG_BYTE_SWAP_64: // Byte swap
        *v = ((*v & 0xFF00FF00FF00FF00) >> 8) |
             ((*v & 0x00FF00FF00FF00FF) << 8);

        return(TNG_SUCCESS);

    case TNG_BIG_ENDIAN_64: // Already correct
        return(TNG_SUCCESS);

    default:
        return(TNG_FAILURE);
    }
}

/** This function swaps the byte order of a 32 bit numerical variable
 * to little endian.
 * It does not only work with integer, but e.g. floats need casting.
 * If the byte order is already little endian no change is needed.
 * @param tng_data is a trajectory data container.
 * @param v is a pointer to a 32 bit numerical value (float or integer).
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the current
 * byte order is not recognised.
 */
static tng_function_status tng_swap_byte_order_little_endian_32
                (const tng_trajectory_t tng_data, int32_t *v)
{
    switch(tng_data->endianness_32)
    {
    case TNG_LITTLE_ENDIAN_32: // Already correct
        return(TNG_SUCCESS);

    case TNG_BYTE_PAIR_SWAP_32: // byte pair swapped big endian to little endian
        *v = ((*v & 0xFF00FF00) >> 8) |
             ((*v & 0x00FF00FF) << 8);

        return(TNG_SUCCESS);

    case TNG_BIG_ENDIAN_32: // Byte order is reversed.
        *v = ((*v & 0xFF000000) >> 24) | // Move first byte to end
             ((*v & 0x00FF0000) >> 8) |  // Move 2nd byte to pos 3
             ((*v & 0x0000FF00) << 8) |  // Move 3rd byte to pos 2
             ((*v & 0x000000FF) << 24);  // Move last byte to first

        return(TNG_SUCCESS);

    default:
        return(TNG_FAILURE);
    }
}

/** This function swaps the byte order of a 64 bit numerical variable
 * to little endian.
 * It does not only work with integer, but e.g. floats need casting.
 * The byte order swapping routine can convert four different byte
 * orders to little endian.
 * If the byte order is already little endian no change is needed.
 * @param tng_data is a trajectory data container.
 * @param v is a pointer to a 64 bit numerical value (double or integer).
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the current
 * byte order is not recognised.
 */
static tng_function_status tng_swap_byte_order_little_endian_64
                (const tng_trajectory_t tng_data, int64_t *v)
{
    switch(tng_data->endianness_64)
    {
    case TNG_LITTLE_ENDIAN_64: // Already correct
        return(TNG_SUCCESS);

    case TNG_QUAD_SWAP_64: // Byte quad swapped big endian to little endian
        *v = ((*v & 0xFF000000FF000000) >> 24) |
             ((*v & 0x00FF000000FF0000) >> 8) |
             ((*v & 0x0000FF000000FF00) << 8) |
             ((*v & 0x000000FF000000FF) << 24);

        return(TNG_SUCCESS);

    case TNG_BYTE_PAIR_SWAP_64: // Byte pair swapped big endian to little endian
        *v = ((*v & 0xFF00FF0000000000) >> 40) |
             ((*v & 0x00FF00FF00000000) >> 24) |
             ((*v & 0x00000000FF00FF00) << 24) |
             ((*v & 0x0000000000FF00FF) << 40);

        return(TNG_SUCCESS);

    case TNG_BYTE_SWAP_64: // Byte swapped big endian to little endian
        *v = ((*v & 0xFFFF000000000000) >> 48) |
             ((*v & 0x0000FFFF00000000) >> 16) |
             ((*v & 0x00000000FFFF0000) << 16) |
             ((*v & 0x000000000000FFFF) << 48);

        return(TNG_SUCCESS);

    case TNG_BIG_ENDIAN_64: // Byte order is reversed.
        *v = ((*v & 0xFF00000000000000) >> 56) | // Move first byte to end
             ((*v & 0x00FF000000000000) >> 40) | // Move 2nd byte to pos 7
             ((*v & 0x0000FF0000000000) >> 24) | // Move 3rd byte to pos 6
             ((*v & 0x000000FF00000000) >> 8 ) | // Move 4th byte to pos 5
             ((*v & 0x00000000FF000000) << 8 ) | // Move 5th byte to pos 4
             ((*v & 0x0000000000FF0000) << 24) | // Move 6th byte to pos 3
             ((*v & 0x000000000000FF00) << 40) | // Move 7th byte to pos 2
             ((*v & 0x00000000000000FF) << 56);  // Move last byte to first

        return(TNG_SUCCESS);

    default:
        return(TNG_FAILURE);
    }
}
/** Generate the md5 hash of a block.
 * The hash is created based on the actual block contents.
 * @param block is a general block container.
 * @return TNG_SUCCESS (0) if successful.
 */
static tng_function_status tng_block_hash_generate(tng_gen_block_t block)
{
    md5_state_t md5_state;

    md5_init(&md5_state);
    md5_append(&md5_state, (md5_byte_t *)block->block_contents,
               block->block_contents_size);
    md5_finish(&md5_state, (md5_byte_t *)block->hash);

    return(TNG_SUCCESS);
}

/** Compare the current block hash (e.g. read from file) with the hash
 * calculated from the current contents.
 * If the current hash is not set skip the comparison.
 * @param block is a general block container.
 * @param results If the hashes match results is set to TNG_TRUE, otherwise it is
 * set to TNG_FALSE. If the hash was not set results is set to TNG_TRUE.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (1) if the hash was not
 * set.
 */
static tng_function_status hash_match_verify(tng_gen_block_t block,
                                             tng_bool *results)
{
    md5_state_t md5_state;
    char hash[TNG_HASH_LEN];

    if(strncmp(block->hash, "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", 16) == 0)
    {
        *results = TNG_TRUE;
        return(TNG_FAILURE);
    }
    md5_init(&md5_state);
    md5_append(&md5_state, (md5_byte_t *)block->block_contents,
               block->block_contents_size);
    md5_finish(&md5_state, (md5_byte_t *)hash);

    if(strncmp(block->hash, hash, 16) != 0)
    {
        *results = TNG_FALSE;
//         int i;
//         printf("Hash in file:  ");
//         for(i = 0; i<16; i++)
//         {
//             printf("%hhX ", block->hash[i]);
//         }
//         printf("\n");
//         printf("Expected hash: ");
//         for(i = 0; i<16; i++)
//         {
//             printf("%hhX ", hash[i]);
//         }
//         printf("\n");
//         if(block->block_contents_size < 100)
//         {
//             for(i = 0; i < block->block_contents_size; i++)
//             {
//                 printf("%hhX ", block->block_contents[i]);
//             }
//         }
    }
    else
    {
        *results = TNG_TRUE;
    }

    return(TNG_SUCCESS);
}

/** Open the input file if it is not already opened.
 * @param tng_data is a trajectory data container.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_input_file_init(tng_trajectory_t tng_data)
{
    if(!tng_data->input_file)
    {
        if(!tng_data->input_file_path)
        {
            printf("No file specified for reading. %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        tng_data->input_file = fopen(tng_data->input_file_path, "r");
        if(!tng_data->input_file)
        {
            printf("Cannot open file %s. %s: %d\n",
                   tng_data->input_file_path, __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }
    return(TNG_SUCCESS);
}

/** Open the output file if it is not already opened
 * @param tng_data is a trajectory data container.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_output_file_init(tng_trajectory_t tng_data)
{
    if(!tng_data->output_file)
    {
        if(!tng_data->output_file_path)
        {
            printf("No file specified for writing. %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }

        tng_data->output_file = fopen(tng_data->output_file_path, "w+");

        if(!tng_data->output_file)
        {
            printf("Cannot open file %s. %s: %d\n",
                   tng_data->output_file_path, __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }
    return(TNG_SUCCESS);
}

/** Setup a file block container.
 * @param block_p a pointer to memory to initialise as a file block container.
 * @details Memory is allocated during initialisation.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_block_init(struct tng_gen_block **block_p)
{
//     printf("In tng_block_init\n");

    tng_gen_block_t block;

    *block_p = malloc(sizeof(struct tng_gen_block));

    if(!*block_p)
    {
        printf("Cannot allocate memory (%lu bytes). %s: %d\n",
               sizeof(struct tng_gen_block), __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    block = *block_p;

    block->id = -1;
    /* Reset the hash */
    memcpy(block->hash, "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", TNG_HASH_LEN);
    block->name = 0;
    block->block_version = TNG_VERSION;
    block->header_contents = 0;
    block->header_contents_size = 0;
    block->block_contents = 0;
    block->block_contents_size = 0;

    return(TNG_SUCCESS);
}

/**
 * @brief Clean up a file block container.
 * @param block_p a pointer to the file block container to destroy.
 * @details All allocated memory in the data structure is freed, as well as
 * block_p itself.
 * @return TNG_SUCCESS (0) if successful.
 */
static tng_function_status tng_block_destroy(struct tng_gen_block **block_p)
{
    tng_gen_block_t block = *block_p;

    if(!*block_p)
    {
        return(TNG_SUCCESS);
    }

//     printf("Destroying block\n");
    if(block->name)
    {
        free(block->name);
        block->name = 0;
    }
    if(block->header_contents)
    {
        free(block->header_contents);
        block->header_contents = 0;
    }
    if(block->block_contents)
    {
        free(block->block_contents);
        block->block_contents = 0;
    }

    free(*block_p);
    *block_p = 0;

    return(TNG_SUCCESS);
}

/** Read the header of a data block, regardless of its type
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_block_header_read
                (tng_trajectory_t tng_data, tng_gen_block_t block)
{
    int len, offset = 0;

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    /* First read the header size to be able to read the whole header. */
    if(fread(&block->header_contents_size, sizeof(block->header_contents_size),
        1, tng_data->input_file) == 0)
    {
        printf("Cannot read header size. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* If this was the size of the general info block check the endianness */
    if(ftell(tng_data->input_file) < 9)
    {
        /* File is little endian */
        if ( *(const uint8_t*)&block->header_contents_size != 0x00 &&
             *(const uint8_t*)(&block->header_contents_size + 7) == 0x00)
        {
            /* If the architecture endianness is little endian no byte swap
             * will be needed. Otherwise use the functions to swap to little
             * endian */
            if(tng_data->endianness_32 == TNG_LITTLE_ENDIAN_32)
            {
                tng_data->input_endianness_swap_func_32 = 0;
            }
            else
            {
                tng_data->input_endianness_swap_func_32 =
                &tng_swap_byte_order_little_endian_32;
            }
            if(tng_data->endianness_64 == TNG_LITTLE_ENDIAN_64)
            {
                tng_data->input_endianness_swap_func_64 = 0;
            }
            else
            {
                tng_data->input_endianness_swap_func_64 =
                &tng_swap_byte_order_little_endian_64;
            }
        }
        /* File is big endian */
        else
        {
            /* If the architecture endianness is big endian no byte swap
             * will be needed. Otherwise use the functions to swap to big
             * endian */
            if(tng_data->endianness_32 == TNG_BIG_ENDIAN_32)
            {
                tng_data->input_endianness_swap_func_32 = 0;
            }
            else
            {
                tng_data->input_endianness_swap_func_32 =
                &tng_swap_byte_order_big_endian_32;
            }
            if(tng_data->endianness_64 == TNG_BIG_ENDIAN_64)
            {
                tng_data->input_endianness_swap_func_64 = 0;
            }
            else
            {
                tng_data->input_endianness_swap_func_64 =
                &tng_swap_byte_order_big_endian_64;
            }
        }
    }

    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &block->header_contents_size)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    /* Move the reading position to the beginning of the header. */
    fseek(tng_data->input_file, -sizeof(block->header_contents_size),
          SEEK_CUR);

    /* If there is already memory allocated for the contents free it (we do not
     * know if the size is correct). */
    if(block->header_contents)
    {
        free(block->header_contents);
    }

    block->header_contents = malloc(block->header_contents_size);
    if(!block->header_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->header_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* Read the whole header into header_contents. This way it can be saved
     * even if it cannot be interpreted
     * for one reason or another. */
    if(fread(block->header_contents, block->header_contents_size, 1,
        tng_data->input_file) == 0)
    {
        printf("Cannot read header. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* The header contents size has already been read. Skip ahead. */
    offset = sizeof(block->header_contents_size);


    /* Copy the respective parameters from the header contents block */
    memcpy(&block->block_contents_size, block->header_contents+offset,
           sizeof(block->block_contents_size));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &block->block_contents_size)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    offset += sizeof(block->block_contents_size);

    memcpy(&block->id, block->header_contents+offset, sizeof(block->id));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &block->id)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    offset += sizeof(block->id);

    memcpy(block->hash, block->header_contents+offset, TNG_HASH_LEN);
    offset += TNG_HASH_LEN;

    if(block->name && strcmp(block->name, block->header_contents+offset) != 0)
    {
        free(block->name);
        block->name = 0;
    }
    len = tng_min(strlen(block->header_contents+offset) + 1, TNG_MAX_STR_LEN);
    if(!block->name)
    {
        block->name = malloc(len);
        if(!block->name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                    __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        strncpy(block->name, block->header_contents+offset, len);
    }
    offset += len;

    memcpy(&block->block_version, block->header_contents+offset,
           sizeof(block->block_version));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &block->block_version)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(block->block_version);

    return(TNG_SUCCESS);
}

/** Write a whole block, both header and contents, regardless of it type
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
static tng_function_status tng_block_verbatim_write(tng_trajectory_t tng_data,
                                                    tng_gen_block_t block)
{
    if(!block->header_contents)
    {
        printf("No contents to write. %s: %d\n", __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    if(fwrite(block->header_contents, block->header_contents_size, 1,
                tng_data->output_file) != 1)
    {
        printf("Could not write all header data. %s: %d\n",
                __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(!block->block_contents)
    {
        printf("No block data to write. %s: %d\n",
                __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    if(fwrite(block->block_contents, block->block_contents_size, 1,
                tng_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n",
                __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }
    return(TNG_SUCCESS);
}

/** Write the header of a data block, regardless of its type
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH an md5 hash will be generated and written.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_block_header_write
                (tng_trajectory_t tng_data,
                 tng_gen_block_t block,
                 const tng_hash_mode hash_mode)
{
    int name_len, offset = 0;

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        printf("Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }


    if(!block->name)
    {
        block->name = malloc(1);
        if(!block->name)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        block->name[0] = 0;
    }

    name_len = tng_min(strlen(block->name) + 1, TNG_MAX_STR_LEN);

    if(hash_mode == TNG_USE_HASH)
    {
        tng_block_hash_generate(block);
    }

    /* Calculate the size of the header to write */
    block->header_contents_size = sizeof(block->header_contents_size) +
                                  sizeof(block->block_contents_size) +
                                  sizeof(block->id) +
                                  sizeof(block->block_version) +
                                  TNG_HASH_LEN +
                                  name_len;

    if(block->header_contents)
    {
        free(block->header_contents);
    }

    block->header_contents = malloc(block->header_contents_size);
    if(!block->header_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->header_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* First copy all data into the header_contents block and finally write
     * the whole block at once. */
    memcpy(block->header_contents, &block->header_contents_size,
           sizeof(block->header_contents_size));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(block->header_contents_size);

    memcpy(block->header_contents+offset, &block->block_contents_size,
           sizeof(block->block_contents_size));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(block->block_contents_size);

    memcpy(block->header_contents+offset, &block->id, sizeof(block->id));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(block->id);

    memcpy(block->header_contents+offset, block->hash, TNG_HASH_LEN);
    offset += TNG_HASH_LEN;

    strncpy(block->header_contents+offset, block->name, name_len);
    offset += name_len;

    memcpy(block->header_contents+offset, &block->block_version,
           sizeof(block->block_version));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(block->block_version);

    if(fwrite(block->header_contents, block->header_contents_size,
       1, tng_data->output_file) != 1)
    {
        printf("Could not write all header data. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }
    return(TNG_SUCCESS);
}

/** Read a general info block. This is the first block of a TNG file.
 *  Populate the fields in tng_data.
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_general_info_block_read
                (tng_trajectory_t tng_data, tng_gen_block_t block,
                 const tng_hash_mode hash_mode)
{
    int len, offset = 0;
    tng_bool same_hash;

    void *temp;

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    temp = realloc(block->block_contents, block->block_contents_size);
    if(!temp)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        free(block->block_contents);
        return(TNG_CRITICAL);
    }
    block->block_contents = temp;

    /* Read the whole block into block_contents to be able to write it to disk
     * even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
             tng_data->input_file) == 0)
    {
        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */

    if(hash_mode == TNG_USE_HASH)
    {
        hash_match_verify(block, &same_hash);
        if(same_hash != TNG_TRUE)
        {
            printf("General info block contents corrupt. Hashes do not match. "
                "%s: %d\n",
                __FILE__, __LINE__);
    //         return(TNG_FAILURE);
        }
    }

    len = tng_min(strlen(block->block_contents) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->first_program_name, len);
    if(!temp)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->first_program_name);
        return(TNG_CRITICAL);
    }
    tng_data->first_program_name = temp;
    strncpy(tng_data->first_program_name, block->block_contents, len);
    offset += len;

    len = tng_min(strlen(block->block_contents + offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->last_program_name, len);
    if(!temp)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->last_program_name);
        return(TNG_CRITICAL);
    }
    tng_data->last_program_name = temp;
    strncpy(tng_data->last_program_name, block->block_contents, len);
    offset += len;

    len = tng_min(strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->first_user_name, len);
    if(!temp)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->first_user_name);
        return(TNG_CRITICAL);
    }
    tng_data->first_user_name = temp;
    strncpy(tng_data->first_user_name, block->block_contents+offset, len);
    offset += len;

    len = tng_min(strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->last_user_name, len);
    if(!temp)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->last_user_name);
        return(TNG_CRITICAL);
    }
    tng_data->last_user_name = temp;
    strncpy(tng_data->last_user_name, block->block_contents+offset, len);
    offset += len;

    len = tng_min(strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->first_computer_name, len);
    if(!temp)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->first_computer_name);
        return(TNG_CRITICAL);
    }
    tng_data->first_computer_name = temp;
    strncpy(tng_data->first_computer_name, block->block_contents+offset, len);
    offset += len;

    len = tng_min(strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->last_computer_name, len);
    if(!temp)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->last_computer_name);
        return(TNG_CRITICAL);
    }
    tng_data->last_computer_name = temp;
    strncpy(tng_data->last_computer_name, block->block_contents+offset, len);
    offset += len;

    len = tng_min(strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->first_pgp_signature, len);
    if(!temp)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->first_pgp_signature);
        return(TNG_CRITICAL);
    }
    tng_data->first_pgp_signature = temp;
    strncpy(tng_data->first_pgp_signature, block->block_contents+offset, len);
    offset += len;

    len = tng_min(strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->last_pgp_signature, len);
    if(!temp)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->last_pgp_signature);
        return(TNG_CRITICAL);
    }
    tng_data->last_pgp_signature = temp;
    strncpy(tng_data->last_pgp_signature, block->block_contents+offset, len);
    offset += len;

    len = tng_min(strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->forcefield_name, len);
    if(!temp)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->forcefield_name);
        return(TNG_CRITICAL);
    }
    tng_data->forcefield_name = temp;
    strncpy(tng_data->forcefield_name, block->block_contents+offset, len);
    offset += len;

    memcpy(&tng_data->time, block->block_contents+offset,
           sizeof(tng_data->time));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &tng_data->time)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->time);

    memcpy(&tng_data->var_num_atoms_flag, block->block_contents+offset,
           sizeof(tng_data->var_num_atoms_flag));
    offset += sizeof(tng_data->var_num_atoms_flag);

    memcpy(&tng_data->frame_set_n_frames, block->block_contents+offset,
           sizeof(tng_data->frame_set_n_frames));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                 &tng_data->frame_set_n_frames)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->frame_set_n_frames);

    memcpy(&tng_data->first_trajectory_frame_set_input_file_pos,
           block->block_contents+offset,
           sizeof(tng_data->first_trajectory_frame_set_input_file_pos));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                          &tng_data->first_trajectory_frame_set_input_file_pos)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->first_trajectory_frame_set_input_file_pos);

    tng_data->current_trajectory_frame_set.next_frame_set_file_pos =
    tng_data->first_trajectory_frame_set_input_file_pos;


    memcpy(&tng_data->last_trajectory_frame_set_input_file_pos,
           block->block_contents+offset,
           sizeof(tng_data->last_trajectory_frame_set_input_file_pos));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                          &tng_data->last_trajectory_frame_set_input_file_pos)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->last_trajectory_frame_set_input_file_pos);

    memcpy(&tng_data->medium_stride_length, block->block_contents+offset,
           sizeof(tng_data->medium_stride_length));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                               &tng_data->medium_stride_length)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->medium_stride_length);

    memcpy(&tng_data->long_stride_length, block->block_contents+offset,
           sizeof(tng_data->long_stride_length));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                 &tng_data->long_stride_length)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    return(TNG_SUCCESS);
}

/** Write a general info block. This is the first block of a TNG file.
 * @param tng_data is a trajectory data container.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH an md5 hash will be generated and written.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_general_info_block_write
                (tng_trajectory_t tng_data,
                 const tng_hash_mode hash_mode)
{
    int first_program_name_len, first_user_name_len;
    int first_computer_name_len, first_pgp_signature_len;
    int last_program_name_len, last_user_name_len;
    int last_computer_name_len, last_pgp_signature_len;
    int forcefield_name_len, name_len;
    int offset = 0;
    tng_gen_block_t block;

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    /* If the strings are unallocated allocate memory for just string
     * termination */
    if(!tng_data->first_program_name)
    {
        tng_data->first_program_name = malloc(1);
        if(!tng_data->first_program_name)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        tng_data->first_program_name[0] = 0;
    }
    if(!tng_data->last_program_name)
    {
        tng_data->last_program_name = malloc(1);
        if(!tng_data->last_program_name)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        tng_data->last_program_name[0] = 0;
    }
    if(!tng_data->first_user_name)
    {
        tng_data->first_user_name = malloc(1);
        if(!tng_data->first_user_name)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        tng_data->first_user_name[0] = 0;
    }
    if(!tng_data->last_user_name)
    {
        tng_data->last_user_name = malloc(1);
        if(!tng_data->last_user_name)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        tng_data->last_user_name[0] = 0;
    }
    if(!tng_data->first_computer_name)
    {
        tng_data->first_computer_name = malloc(1);
        if(!tng_data->first_computer_name)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        tng_data->first_computer_name[0] = 0;
    }
    if(!tng_data->last_computer_name)
    {
        tng_data->last_computer_name = malloc(1);
        if(!tng_data->last_computer_name)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        tng_data->last_computer_name[0] = 0;
    }
    if(!tng_data->first_pgp_signature)
    {
        tng_data->first_pgp_signature = malloc(1);
        if(!tng_data->first_pgp_signature)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        tng_data->first_pgp_signature[0] = 0;
    }
    if(!tng_data->last_pgp_signature)
    {
        tng_data->last_pgp_signature = malloc(1);
        if(!tng_data->last_pgp_signature)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        tng_data->last_pgp_signature[0] = 0;
    }
    if(!tng_data->forcefield_name)
    {
        tng_data->forcefield_name = malloc(1);
        if(!tng_data->forcefield_name)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        tng_data->forcefield_name[0] = 0;
    }

    tng_block_init(&block);

    name_len = strlen("GENERAL INFO");

    block->name = malloc(name_len + 1);
    if(!block->name)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n",
                name_len+1, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    strcpy(block->name, "GENERAL INFO");
    block->id = TNG_GENERAL_INFO;

    first_program_name_len = tng_min(strlen(tng_data->first_program_name) + 1,
                           TNG_MAX_STR_LEN);
    last_program_name_len = tng_min(strlen(tng_data->last_program_name) + 1,
                           TNG_MAX_STR_LEN);
    first_user_name_len = tng_min(strlen(tng_data->first_user_name) + 1,
                        TNG_MAX_STR_LEN);
    last_user_name_len = tng_min(strlen(tng_data->last_user_name) + 1,
                        TNG_MAX_STR_LEN);
    first_computer_name_len = tng_min(strlen(tng_data->first_computer_name) + 1,
                            TNG_MAX_STR_LEN);
    last_computer_name_len = tng_min(strlen(tng_data->last_computer_name) + 1,
                            TNG_MAX_STR_LEN);
    first_pgp_signature_len = tng_min(strlen(tng_data->first_pgp_signature) + 1,
                            TNG_MAX_STR_LEN);
    last_pgp_signature_len = tng_min(strlen(tng_data->last_pgp_signature) + 1,
                            TNG_MAX_STR_LEN);
    forcefield_name_len = tng_min(strlen(tng_data->forcefield_name) + 1,
                              TNG_MAX_STR_LEN);

    block->block_contents_size = sizeof(tng_data->time) +
                sizeof(tng_data->var_num_atoms_flag) +
                sizeof(tng_data->frame_set_n_frames) +
                sizeof(tng_data->first_trajectory_frame_set_input_file_pos) +
                sizeof(tng_data->last_trajectory_frame_set_input_file_pos) +
                sizeof(tng_data->medium_stride_length) +
                sizeof(tng_data->long_stride_length) +
                first_program_name_len +
                last_program_name_len +
                first_user_name_len +
                last_user_name_len +
                first_computer_name_len +
                last_computer_name_len +
                first_pgp_signature_len +
                last_pgp_signature_len +
                forcefield_name_len;

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    strncpy(block->block_contents, tng_data->first_program_name, first_program_name_len);
    offset += first_program_name_len;

    strncpy(block->block_contents+offset, tng_data->last_program_name, last_program_name_len);
    offset += last_program_name_len;

    strncpy(block->block_contents+offset, tng_data->first_user_name, first_user_name_len);
    offset += first_user_name_len;

    strncpy(block->block_contents+offset, tng_data->last_user_name, last_user_name_len);
    offset += last_user_name_len;

    strncpy(block->block_contents+offset, tng_data->first_computer_name,
            first_computer_name_len);
    offset += first_computer_name_len;

    strncpy(block->block_contents+offset, tng_data->last_computer_name,
            last_computer_name_len);
    offset += last_computer_name_len;

    strncpy(block->block_contents+offset, tng_data->first_pgp_signature,
            first_pgp_signature_len);
    offset += first_pgp_signature_len;

    strncpy(block->block_contents+offset, tng_data->last_pgp_signature,
            last_pgp_signature_len);
    offset += last_pgp_signature_len;

    strncpy(block->block_contents+offset, tng_data->forcefield_name,
            forcefield_name_len);
    offset += forcefield_name_len;

    memcpy(block->block_contents+offset, &tng_data->time,
           sizeof(tng_data->time));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->time);

    memcpy(block->block_contents+offset, &tng_data->var_num_atoms_flag,
           sizeof(tng_data->var_num_atoms_flag));
    offset += sizeof(tng_data->var_num_atoms_flag);

    memcpy(block->block_contents+offset, &tng_data->frame_set_n_frames,
           sizeof(tng_data->frame_set_n_frames));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->frame_set_n_frames);

    memcpy(block->block_contents+offset,
           &tng_data->first_trajectory_frame_set_input_file_pos,
           sizeof(tng_data->first_trajectory_frame_set_input_file_pos));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->first_trajectory_frame_set_input_file_pos);

    memcpy(block->block_contents+offset,
           &tng_data->last_trajectory_frame_set_input_file_pos,
           sizeof(tng_data->last_trajectory_frame_set_input_file_pos));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->last_trajectory_frame_set_input_file_pos);

    memcpy(block->block_contents+offset, &tng_data->medium_stride_length,
           sizeof(tng_data->medium_stride_length));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->medium_stride_length);

    memcpy(block->block_contents+offset, &tng_data->long_stride_length,
           sizeof(tng_data->long_stride_length));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(tng_block_header_write(tng_data, block, hash_mode) != TNG_SUCCESS)
    {
        printf("Cannot write header of file %s. %s: %d\n",
               tng_data->output_file_path, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
        tng_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

/** Read the chain data of a molecules block.
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @param chain is the chain data container.
 * @param offset is the offset of the block input and is updated when reading.
 * @return TNG_SUCCESS(0) is successful.
 */
static tng_function_status tng_chain_data_read(tng_trajectory_t tng_data,
                                               tng_gen_block_t block,
                                               tng_chain_t chain,
                                               int *offset)
{
    int len;

    memcpy(&chain->id, block->block_contents+*offset,
            sizeof(chain->id));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &chain->id)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(chain->id);

    len = tng_min(strlen(block->block_contents+*offset) + 1,
            TNG_MAX_STR_LEN);
    chain->name = malloc(len);
    strncpy(chain->name,
            block->block_contents+*offset, len);
    *offset += len;

    memcpy(&chain->n_residues, block->block_contents+*offset,
        sizeof(chain->n_residues));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &chain->n_residues)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(chain->n_residues);

    return(TNG_SUCCESS);
}

/** Write the chain data of a molecules block.
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @param chain is the chain data container.
 * @param offset is the offset of the block output and is updated when writing.
 * @return TNG_SUCCESS(0) is successful.
 */
static tng_function_status tng_chain_data_write(tng_trajectory_t tng_data,
                                                tng_gen_block_t block,
                                                tng_chain_t chain,
                                                int *offset)
{
    int len;

    memcpy(block->block_contents+*offset, &chain->id, sizeof(chain->id));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                    (int64_t *)block->header_contents+*offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(chain->id);

    len = tng_min(strlen(chain->name) + 1, TNG_MAX_STR_LEN);
    strncpy(block->block_contents + *offset, chain->name, len);
    *offset += len;

    memcpy(block->block_contents+*offset, &chain->n_residues,
        sizeof(chain->n_residues));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                    (int64_t *)block->header_contents+*offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(chain->n_residues);

    return(TNG_SUCCESS);
}

/** Read the residue data of a molecules block.
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @param residue is the residue data container.
 * @param offset is the offset of the block input and is updated when reading.
 * @return TNG_SUCCESS(0) is successful.
 */
static tng_function_status tng_residue_data_read(tng_trajectory_t tng_data,
                                                 tng_gen_block_t block,
                                                 tng_residue_t residue,
                                                 int *offset)
{
    int len;

    memcpy(&residue->id, block->block_contents+*offset,
        sizeof(residue->id));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &residue->id)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(residue->id);

    len = tng_min(strlen(block->block_contents+*offset) + 1,
            TNG_MAX_STR_LEN);
    residue->name = malloc(len);
    strncpy(residue->name,
            block->block_contents+*offset, len);
    *offset += len;

    memcpy(&residue->n_atoms, block->block_contents+*offset,
            sizeof(residue->n_atoms));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &residue->n_atoms)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(residue->n_atoms);

    return(TNG_SUCCESS);
}

/** Write the residue data of a molecules block.
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @param residue is the residue data container.
 * @param offset is the offset of the block output and is updated when writing.
 * @return TNG_SUCCESS(0) is successful.
 */
static tng_function_status tng_residue_data_write(tng_trajectory_t tng_data,
                                                  tng_gen_block_t block,
                                                  tng_residue_t residue,
                                                  int *offset)
{
    int len;

    memcpy(block->block_contents+*offset, &residue->id, sizeof(residue->id));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                    (int64_t *)block->header_contents+*offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(residue->id);

    len = tng_min(strlen(residue->name) + 1, TNG_MAX_STR_LEN);
    strncpy(block->block_contents + *offset, residue->name, len);
    *offset += len;

    memcpy(block->block_contents+*offset, &residue->n_atoms,
        sizeof(residue->n_atoms));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                    (int64_t *)block->header_contents+*offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(residue->n_atoms);

    return(TNG_SUCCESS);
}

/** Read the atom data of a molecules block.
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @param atom is the atom data container.
 * @param offset is the offset of the block input and is updated when reading.
 * @return TNG_SUCCESS(0) is successful.
 */
static tng_function_status tng_atom_data_read(tng_trajectory_t tng_data,
                                              tng_gen_block_t block,
                                              tng_atom_t atom,
                                              int *offset)
{
    int len;

    memcpy(&atom->id, block->block_contents+*offset,
        sizeof(atom->id));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                    &atom->id)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(atom->id);

    len = tng_min(strlen(block->block_contents+*offset) + 1,
            TNG_MAX_STR_LEN);
    atom->name = malloc(len);
    strncpy(atom->name,
            block->block_contents+*offset, len);
    *offset += len;

    len = tng_min(strlen(block->block_contents+*offset) + 1,
            TNG_MAX_STR_LEN);
    atom->atom_type = malloc(len);
    strncpy(atom->atom_type,
            block->block_contents+*offset, len);
    *offset += len;

    return(TNG_SUCCESS);
}

/** Write the atom data of a molecules block.
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @param atom is the atom data container.
 * @param offset is the offset of the block output and is updated when writing.
 * @return TNG_SUCCESS(0) is successful.
 */
static tng_function_status tng_atom_data_write(tng_trajectory_t tng_data,
                                               tng_gen_block_t block,
                                               tng_atom_t atom,
                                               int *offset)
{
    int len;

    memcpy(block->block_contents+*offset, &atom->id,
            sizeof(atom->id));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                    (int64_t *)block->header_contents+*offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(atom->id);

    len = tng_min(strlen(atom->name) + 1, TNG_MAX_STR_LEN);
    strncpy(block->block_contents + *offset, atom->name, len);
    *offset += len;

    len = tng_min(strlen(atom->atom_type) + 1, TNG_MAX_STR_LEN);
    strncpy(block->block_contents + *offset, atom->atom_type, len);
    *offset += len;

    return(TNG_SUCCESS);
}

/** Read a molecules block. Contains chain, residue and atom data
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_molecules_block_read
                (tng_trajectory_t tng_data,
                 tng_gen_block_t block,
                 const tng_hash_mode hash_mode)
{
    int i, j, k, l, len, offset = 0;
    tng_molecule_t molecule;
    tng_chain_t chain;
    tng_residue_t residue;
    tng_atom_t atom;
    tng_bond_t bond;
    tng_bool same_hash;

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }

    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to disk
     * even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
             tng_data->input_file) == 0)
    {
        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */

    if(hash_mode == TNG_USE_HASH)
    {
        hash_match_verify(block, &same_hash);
        if(same_hash != TNG_TRUE)
        {
            printf("Molecules block contents corrupt. Hashes do not match. "
                "%s: %d\n",
                __FILE__, __LINE__);
        }
    }

    if(tng_data->molecules)
    {
        for(i=tng_data->n_molecules; i--;)
        {
            tng_molecule_destroy(tng_data, &tng_data->molecules[i]);
        }
        free(tng_data->molecules);
        tng_data->molecules = 0;
        tng_data->n_molecules = 0;
    }

    memcpy(&tng_data->n_molecules, block->block_contents,
           sizeof(tng_data->n_molecules));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &tng_data->n_molecules)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->n_molecules);

    if(tng_data->molecules)
    {
        free(tng_data->molecules);
    }

    tng_data->n_particles = 0;

    tng_data->molecules = malloc(tng_data->n_molecules *
                          sizeof(struct tng_molecule));
    if(!tng_data->molecules)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               tng_data->n_molecules * sizeof(struct tng_molecule),
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(!tng_data->var_num_atoms_flag)
    {
        if(tng_data->molecule_cnt_list)
        {
            free(tng_data->molecule_cnt_list);
        }
        tng_data->molecule_cnt_list = malloc(sizeof(int64_t) *
                                      tng_data->n_molecules);
        if(!tng_data->molecule_cnt_list)
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   tng_data->n_molecules * sizeof(struct tng_molecule),
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    /* Read each molecule from file */
    for(i=0; i < tng_data->n_molecules; i++)
    {
        molecule = &tng_data->molecules[i];

        memcpy(&molecule->id, block->block_contents+offset,
               sizeof(molecule->id));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                       &molecule->id)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->id);

//         printf("Read id: %"PRId64" offset: %d\n", molecule->id, offset);
        len = tng_min(strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
        molecule->name = malloc(len);
        strncpy(molecule->name, block->block_contents+offset, len);
        offset += len;

        memcpy(&molecule->quaternary_str, block->block_contents+offset,
               sizeof(molecule->quaternary_str));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                     &molecule->quaternary_str)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->quaternary_str);

        if(!tng_data->var_num_atoms_flag)
        {
            memcpy(&tng_data->molecule_cnt_list[i],
                   block->block_contents+offset,
                   sizeof(int64_t));
            if(tng_data->input_endianness_swap_func_64)
            {
                if(tng_data->input_endianness_swap_func_64(tng_data,
                                               &tng_data->molecule_cnt_list[i])
                    != TNG_SUCCESS)
                {
                    printf("Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }
            offset += sizeof(int64_t);
        }


        memcpy(&molecule->n_chains, block->block_contents+offset,
               sizeof(molecule->n_chains));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                       &molecule->n_chains)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_chains);

        memcpy(&molecule->n_residues, block->block_contents+offset,
               sizeof(molecule->n_residues));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                       &molecule->n_residues)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_residues);

        memcpy(&molecule->n_atoms, block->block_contents+offset,
               sizeof(molecule->n_atoms));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                       &molecule->n_atoms)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_atoms);

        tng_data->n_particles += molecule->n_atoms *
                                 tng_data->molecule_cnt_list[i];

        if(molecule->n_chains > 0)
        {
            molecule->chains = malloc(molecule->n_chains *
                                    sizeof(struct tng_chain));
            if(!molecule->chains)
            {
                printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                    molecule->n_chains * sizeof(struct tng_chain),
                    __FILE__, __LINE__);
                return(TNG_CRITICAL);
            }

            chain = molecule->chains;
        }
        else
        {
            chain = 0;
        }

        if(molecule->n_residues > 0)
        {
            molecule->residues = malloc(molecule->n_residues *
                                sizeof(struct tng_residue));
            if(!molecule->residues)
            {
                printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                    molecule->n_residues * sizeof(struct tng_residue),
                    __FILE__, __LINE__);
                return(TNG_CRITICAL);
            }

            residue = molecule->residues;
        }
        else
        {
            residue = 0;
        }

        molecule->atoms = malloc(molecule->n_atoms *
                                 sizeof(struct tng_atom));
        if(!molecule->atoms)
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   molecule->n_atoms * sizeof(struct tng_atom),
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }

        atom = molecule->atoms;

        if(molecule->n_chains > 0)
        {
            /* Read the chains of the molecule */
            for(j=molecule->n_chains; j--;)
            {
                chain->molecule = molecule;

                tng_chain_data_read(tng_data, block, chain, &offset);

                chain->residues = residue;

                /* Read the residues of the chain */
                for(k=chain->n_residues; k--;)
                {
                    residue->chain = chain;

                    tng_residue_data_read(tng_data, block, residue, &offset);

                    residue->atoms_offset = atom - molecule->atoms;
                    /* Read the atoms of the residue */
                    for(l=residue->n_atoms; l--;)
                    {
                        atom->residue = residue;

                        tng_atom_data_read(tng_data, block, atom, &offset);

                        atom++;
                    }
                    residue++;
                }
                chain++;
            }
        }
        else
        {
            if(molecule->n_residues > 0)
            {
                for(k=molecule->n_residues; k--;)
                {
                    residue->chain = 0;

                    tng_residue_data_read(tng_data, block, residue, &offset);

                    residue->atoms_offset = atom - molecule->atoms;
                    /* Read the atoms of the residue */
                    for(l=residue->n_atoms; l--;)
                    {
                        atom->residue = residue;

                        tng_atom_data_read(tng_data, block, atom, &offset);

                        atom++;
                    }
                    residue++;
                }
            }
            else
            {
                for(l=molecule->n_atoms; l--;)
                {
                    atom->residue = 0;

                    tng_atom_data_read(tng_data, block, atom, &offset);

                    atom++;
                }
            }
        }

        memcpy(&molecule->n_bonds, block->block_contents+offset,
               sizeof(molecule->n_bonds));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                       &molecule->n_bonds)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_bonds);

        molecule->bonds = malloc(molecule->n_bonds *
                                 sizeof(struct tng_bond));
        if(!molecule->bonds)
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   molecule->n_bonds * sizeof(struct tng_bond),
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }

        bond = molecule->bonds;

        for(j=molecule->n_bonds; j--;)
        {
            memcpy(&bond->from_atom_id, block->block_contents+offset,
                sizeof(bond->from_atom_id));
            if(tng_data->input_endianness_swap_func_64)
            {
                if(tng_data->input_endianness_swap_func_64(tng_data,
                                                           &bond->from_atom_id)
                    != TNG_SUCCESS)
                {
                    printf("Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }
            offset += sizeof(bond->from_atom_id);

            memcpy(&bond->to_atom_id, block->block_contents+offset,
                sizeof(bond->to_atom_id));
            if(tng_data->input_endianness_swap_func_64)
            {
                if(tng_data->input_endianness_swap_func_64(tng_data,
                                                           &bond->to_atom_id)
                    != TNG_SUCCESS)
                {
                    printf("Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }
            offset += sizeof(bond->to_atom_id);

            bond++;
        }
    }

    return(TNG_SUCCESS);
}

/** Write a molecules block.
 * @param tng_data is a trajectory data container.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH an md5 hash will be generated and written.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_molecules_block_write
                (tng_trajectory_t tng_data,
                 const tng_hash_mode hash_mode)
{
    int len = 0, name_len;
    int i, j, k, l, offset = 0;
    tng_molecule_t molecule;
    tng_chain_t chain;
    tng_residue_t residue;
    tng_atom_t atom;
    tng_bond_t bond;
    tng_gen_block_t block;

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    /* First predict the size of the block */
    for(i = 0; i < tng_data->n_molecules; i++)
    {
        molecule = &tng_data->molecules[i];
        if(!molecule->name)
        {
            molecule->name = malloc(1);
            if(!molecule->name)
            {
                printf("Cannot allocate memory (1 byte). %s: %d\n",
                       __FILE__, __LINE__);
                return(TNG_CRITICAL);
            }
            molecule->name[0] = 0;
        }
        len += tng_min(strlen(molecule->name) + 1, TNG_MAX_STR_LEN);

        chain = molecule->chains;
        for(j = molecule->n_chains; j--;)
        {
            len += sizeof(chain->id);

            if(!chain->name)
            {
                chain->name = malloc(1);
                if(!chain->name)
                {
                    printf("Cannot allocate memory (1 byte). %s: %d\n",
                           __FILE__, __LINE__);
                    return(TNG_CRITICAL);
                }
                chain->name[0] = 0;
            }
            len += tng_min(strlen(chain->name) + 1, TNG_MAX_STR_LEN);

            len += sizeof(chain->n_residues);

            chain++;
        }

        residue = molecule->residues;
        for(j = molecule->n_residues; j--;)
        {
            len += sizeof(residue->id);

            if(!residue->name)
            {
                residue->name = malloc(1);
                if(!residue->name)
                {
                    printf("Cannot allocate memory (1 byte). %s: %d\n",
                           __FILE__, __LINE__);
                    return(TNG_CRITICAL);
                }
                residue->name[0] = 0;
            }
            len += tng_min(strlen(residue->name) + 1, TNG_MAX_STR_LEN);

            len += sizeof(residue->n_atoms);

            residue++;
        }

        atom = molecule->atoms;
        for(j = molecule->n_atoms; j--;)
        {
            len += sizeof(atom->id);
            if(!atom->name)
            {
                atom->name = malloc(1);
                if(!atom->name)
                {
                    printf("Cannot allocate memory (1 byte). %s: %d\n",
                           __FILE__, __LINE__);
                    return(TNG_CRITICAL);
                }
                atom->name[0] = 0;
            }
            len += tng_min(strlen(atom->name) + 1, TNG_MAX_STR_LEN);

            if(!atom->atom_type)
            {
                atom->atom_type = malloc(1);
                if(!atom->atom_type)
                {
                    printf("Cannot allocate memory (1 byte). %s: %d\n",
                           __FILE__, __LINE__);
                    return(TNG_CRITICAL);
                }
                atom->atom_type[0] = 0;
            }
            len += tng_min(strlen(atom->atom_type) + 1, TNG_MAX_STR_LEN);

            atom++;
        }

        bond = molecule->bonds;
        for(j = molecule->n_bonds; j--;)
        {
            len += sizeof(bond->from_atom_id) + sizeof(bond->to_atom_id);
        }
    }

    tng_block_init(&block);

    name_len = strlen("MOLECULES");

    block->name = malloc(name_len + 1);
    if(!block->name)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n",
                name_len+1, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    strcpy(block->name, "MOLECULES");
    block->id = TNG_MOLECULES;

    block->block_contents_size = sizeof(tng_data->n_molecules) +
                                 (sizeof(molecule->id) +
                                 sizeof(molecule->quaternary_str) +
                                 sizeof(molecule->n_chains) +
                                 sizeof(molecule->n_residues) +
                                 sizeof(molecule->n_atoms) +
                                 sizeof(molecule->n_bonds)) *
                                 tng_data->n_molecules +
                                 len;

    if(!tng_data->var_num_atoms_flag)
    {
        block->block_contents_size += tng_data->n_molecules * sizeof(int64_t);
    }

    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    memcpy(block->block_contents+offset, &tng_data->n_molecules,
           sizeof(tng_data->n_molecules));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->n_molecules);

    for(i = 0; i < tng_data->n_molecules; i++)
    {
        molecule = &tng_data->molecules[i];
//         printf("i=%d\n", i);
        memcpy(block->block_contents+offset, &molecule->id,
               sizeof(molecule->id));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                                        (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->id);

//         printf("Wrote id: %"PRId64" offset: %d\n", molecule->id, offset);
        len = tng_min(strlen(molecule->name) + 1, TNG_MAX_STR_LEN);
        strncpy(block->block_contents + offset, molecule->name, len);
        offset += len;

        memcpy(block->block_contents+offset, &molecule->quaternary_str,
               sizeof(molecule->quaternary_str));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                                        (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->quaternary_str);

        if(!tng_data->var_num_atoms_flag)
        {
            memcpy(block->block_contents+offset,
                   &tng_data->molecule_cnt_list[i], sizeof(int64_t));
            if(tng_data->output_endianness_swap_func_64)
            {
                if(tng_data->output_endianness_swap_func_64(tng_data,
                                            (int64_t *)block->header_contents+offset)
                    != TNG_SUCCESS)
                {
                    printf("Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }
            offset += sizeof(int64_t);
        }

        memcpy(block->block_contents+offset, &molecule->n_chains,
               sizeof(molecule->n_chains));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                                        (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_chains);

        memcpy(block->block_contents+offset, &molecule->n_residues,
               sizeof(molecule->n_residues));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                                        (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_residues);

        memcpy(block->block_contents+offset, &molecule->n_atoms,
               sizeof(molecule->n_atoms));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                                        (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_atoms);

        if(molecule->n_chains > 0)
        {
            chain = molecule->chains;
            for(j = molecule->n_chains; j--;)
            {
                tng_chain_data_write(tng_data, block, chain, &offset);

                residue = chain->residues;
                for(k = chain->n_residues; k--;)
                {
                    tng_residue_data_write(tng_data, block, residue, &offset);

                    atom = molecule->atoms + residue->atoms_offset;
                    for(l = residue->n_atoms; l--;)
                    {
                        tng_atom_data_write(tng_data, block, atom, &offset);

                        atom++;
                    }
                    residue++;
                }
                chain++;
            }
        }
        else
        {
            if(molecule->n_residues > 0)
            {
                for(k = chain->n_residues; k--;)
                {
                    tng_residue_data_write(tng_data, block, residue, &offset);

                    atom = molecule->atoms + residue->atoms_offset;
                    for(l = residue->n_atoms; l--;)
                    {
                        tng_atom_data_write(tng_data, block, atom, &offset);

                        atom++;
                    }
                    residue++;
                }
            }
            else
            {
                atom = molecule->atoms;
                for(l = residue->n_atoms; l--;)
                {
                    tng_atom_data_write(tng_data, block, atom, &offset);

                    atom++;
                }
            }
        }

        memcpy(block->block_contents+offset, &molecule->n_bonds,
               sizeof(molecule->n_bonds));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                                        (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_bonds);

        bond = molecule->bonds;
        for(j = molecule->n_bonds; j--;)
        {
            memcpy(block->block_contents+offset, &bond->from_atom_id,
                   sizeof(bond->from_atom_id));
            if(tng_data->output_endianness_swap_func_64)
            {
                if(tng_data->output_endianness_swap_func_64(tng_data,
                                            (int64_t *)block->header_contents+offset)
                    != TNG_SUCCESS)
                {
                    printf("Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }
            offset += sizeof(bond->from_atom_id);

            memcpy(block->block_contents+offset, &bond->to_atom_id,
                   sizeof(bond->to_atom_id));
            if(tng_data->output_endianness_swap_func_64)
            {
                if(tng_data->output_endianness_swap_func_64(tng_data,
                                            (int64_t *)block->header_contents+offset)
                    != TNG_SUCCESS)
                {
                    printf("Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }
            offset += sizeof(bond->to_atom_id);

            bond++;
        }
    }

    if(tng_block_header_write(tng_data, block, hash_mode) != TNG_SUCCESS)
    {
        printf("Cannot write header of file %s. %s: %d\n",
               tng_data->output_file_path, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
              tng_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n",
               __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

/** Read a frame set block. Update tng_data->current_trajectory_frame_set
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_frame_set_block_read
                (tng_trajectory_t tng_data,
                 tng_gen_block_t block,
                 const tng_hash_mode hash_mode)
{
    int i, file_pos, offset = 0;
    int64_t prev_n_particles;
    tng_bool same_hash;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;
    tng_particle_mapping_t mapping;

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }

    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to
     * disk even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
             tng_data->input_file) == 0)
    {
        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */

    file_pos = ftell(tng_data->input_file) -
               (block->block_contents_size + block->header_contents_size);

    if(hash_mode == TNG_USE_HASH)
    {
        hash_match_verify(block, &same_hash);
        if(same_hash != TNG_TRUE)
        {
            printf("Frame set block contents corrupt. File pos %d Hashes do not match. "
                "%s: %d\n",
                file_pos, __FILE__, __LINE__);
    //         return(TNG_FAILURE);
        }
    }

    tng_data->current_trajectory_frame_set_input_file_pos = file_pos;

    if(frame_set->n_mapping_blocks && frame_set->mappings)
    {
        for(i = frame_set->n_mapping_blocks; i--;)
        {
            mapping = &frame_set->mappings[i];
            if(mapping->real_particle_numbers)
            {
                free(mapping->real_particle_numbers);
                mapping->real_particle_numbers = 0;
            }
        }
        free(frame_set->mappings);
        frame_set->mappings = 0;
        frame_set->n_mapping_blocks = 0;
    }

    if(tng_data->first_trajectory_frame_set_input_file_pos <= 0)
    {
        tng_data->first_trajectory_frame_set_input_file_pos = file_pos;
    }
    /* FIXME: Should check the frame number instead of the file_pos, in case
     * frame sets are not in order */
    if(tng_data->last_trajectory_frame_set_input_file_pos < file_pos)
    {
        tng_data->last_trajectory_frame_set_input_file_pos = file_pos;
    }

    memcpy(&frame_set->first_frame, block->block_contents,
           sizeof(frame_set->first_frame));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &frame_set->first_frame)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->first_frame);

    memcpy(&frame_set->n_frames, block->block_contents + offset,
           sizeof(frame_set->n_frames));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &frame_set->n_frames)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->n_frames);

    if(tng_data->var_num_atoms_flag)
    {
        prev_n_particles = frame_set->n_particles;
        frame_set->n_particles = 0;
        /* If the list of molecule counts has already been created assume that
         * it is of correct size. */
        if(!frame_set->molecule_cnt_list)
        {
                frame_set->molecule_cnt_list =
                malloc(sizeof(int64_t) * tng_data->n_molecules);

                if(!frame_set->molecule_cnt_list)
                {
                    printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                           sizeof(int64_t) * tng_data->n_molecules,
                           __FILE__, __LINE__);
                    return(TNG_CRITICAL);
                }
        }
        for(i = 0; i < tng_data->n_molecules; i++)
        {
            memcpy(&frame_set->molecule_cnt_list[i],
                   block->block_contents + offset,
                   sizeof(int64_t));
            if(tng_data->input_endianness_swap_func_64)
            {
                if(tng_data->input_endianness_swap_func_64(tng_data,
                                              &frame_set->molecule_cnt_list[i])
                    != TNG_SUCCESS)
                {
                    printf("Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }
            offset += sizeof(int64_t);
            frame_set->n_particles += tng_data->molecules[i].n_atoms *
                                      frame_set->molecule_cnt_list[i];
        }
        if(prev_n_particles && frame_set->n_particles != prev_n_particles)
        {
            /* FIXME: Particle dependent data memory management */
        }
    }

    memcpy(&frame_set->next_frame_set_file_pos,
           block->block_contents + offset,
           sizeof(frame_set->next_frame_set_file_pos));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                           &frame_set->next_frame_set_file_pos)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->next_frame_set_file_pos);

    memcpy(&frame_set->prev_frame_set_file_pos,
           block->block_contents + offset,
           sizeof(frame_set->prev_frame_set_file_pos));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                           &frame_set->prev_frame_set_file_pos)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->prev_frame_set_file_pos);

    memcpy(&frame_set->medium_stride_next_frame_set_file_pos,
           block->block_contents + offset,
           sizeof(frame_set->medium_stride_next_frame_set_file_pos));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                             &frame_set->medium_stride_next_frame_set_file_pos)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->medium_stride_next_frame_set_file_pos);

    memcpy(&frame_set->medium_stride_prev_frame_set_file_pos,
           block->block_contents + offset,
           sizeof(frame_set->medium_stride_prev_frame_set_file_pos));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                             &frame_set->medium_stride_prev_frame_set_file_pos)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->medium_stride_prev_frame_set_file_pos);

    memcpy(&frame_set->long_stride_next_frame_set_file_pos,
           block->block_contents + offset,
           sizeof(frame_set->long_stride_next_frame_set_file_pos));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                               &frame_set->long_stride_next_frame_set_file_pos)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->long_stride_next_frame_set_file_pos);

    memcpy(&frame_set->long_stride_prev_frame_set_file_pos,
           block->block_contents + offset,
           sizeof(frame_set->long_stride_prev_frame_set_file_pos));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                               &frame_set->long_stride_prev_frame_set_file_pos)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->long_stride_prev_frame_set_file_pos);

    return(TNG_SUCCESS);
}

/** Write tng_data->current_trajectory_frame_set to file
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH an md5 hash will be generated and written.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_frame_set_block_write
                (tng_trajectory_t tng_data,
                 tng_gen_block_t block,
                 const tng_hash_mode hash_mode)
{
    char *temp_name;
    int64_t i;
    int offset = 0, name_len;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    name_len = strlen("TRAJECTORY FRAME SET");

    if(!block->name || strlen(block->name) < name_len)
    {
        temp_name = realloc(block->name, name_len + 1);
        if(!temp_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                   name_len+1, __FILE__, __LINE__);
            free(block->name);
            return(TNG_CRITICAL);
        }
        block->name = temp_name;
    }
    strcpy(block->name, "TRAJECTORY FRAME SET");
    block->id = TNG_TRAJECTORY_FRAME_SET;

    block->block_contents_size = sizeof(int64_t) * 8;
    if(tng_data->var_num_atoms_flag)
    {
        block->block_contents_size += sizeof(int64_t) * tng_data->n_molecules;
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    memcpy(block->block_contents, &frame_set->first_frame,
           sizeof(frame_set->first_frame));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->first_frame);

    memcpy(block->block_contents+offset, &frame_set->n_frames,
           sizeof(frame_set->n_frames));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->n_frames);

    if(tng_data->var_num_atoms_flag)
    {
        for(i = 0; i < tng_data->n_molecules; i++)
        {
            memcpy(block->block_contents+offset,
                   &frame_set->molecule_cnt_list[i],
                   sizeof(int64_t));
            if(tng_data->output_endianness_swap_func_64)
            {
                if(tng_data->output_endianness_swap_func_64(tng_data,
                                            (int64_t *)block->header_contents+offset)
                    != TNG_SUCCESS)
                {
                    printf("Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }
            offset += sizeof(int64_t);
        }
    }


    memcpy(block->block_contents+offset, &frame_set->next_frame_set_file_pos,
           sizeof(frame_set->next_frame_set_file_pos));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->next_frame_set_file_pos);

    memcpy(block->block_contents+offset, &frame_set->prev_frame_set_file_pos,
           sizeof(frame_set->prev_frame_set_file_pos));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->prev_frame_set_file_pos);

    memcpy(block->block_contents+offset,
           &frame_set->medium_stride_next_frame_set_file_pos,
           sizeof(frame_set->medium_stride_next_frame_set_file_pos));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->medium_stride_next_frame_set_file_pos);

    memcpy(block->block_contents+offset,
           &frame_set->medium_stride_prev_frame_set_file_pos,
           sizeof(frame_set->medium_stride_prev_frame_set_file_pos));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->medium_stride_prev_frame_set_file_pos);

    memcpy(block->block_contents+offset,
           &frame_set->long_stride_next_frame_set_file_pos,
           sizeof(frame_set->long_stride_next_frame_set_file_pos));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->long_stride_next_frame_set_file_pos);

    memcpy(block->block_contents+offset,
           &frame_set->long_stride_prev_frame_set_file_pos,
           sizeof(frame_set->long_stride_prev_frame_set_file_pos));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->long_stride_prev_frame_set_file_pos);

    if(tng_block_header_write(tng_data, block, hash_mode) != TNG_SUCCESS)
    {
        printf("Cannot write header of file %s. %s: %d\n",
               tng_data->output_file_path, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
              tng_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    return(TNG_SUCCESS);
}


/** Read an atom mappings block (translating between real atom indexes and how
 *  the atom info is written in this frame set).
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_trajectory_mapping_block_read
                (tng_trajectory_t tng_data,
                 tng_gen_block_t block,
                 const tng_hash_mode hash_mode)
{
    int64_t i;
    int offset = 0;
    tng_bool same_hash;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;

    tng_particle_mapping_t mapping, mappings;

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }

    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to disk
     *  even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
        tng_data->input_file) == 0)
    {
        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */

    if(hash_mode == TNG_USE_HASH)
    {
        hash_match_verify(block, &same_hash);
        if(same_hash != TNG_TRUE)
        {
            printf("Particle mapping block contents corrupt. Hashes do not match. "
                "%s: %d\n",
                __FILE__, __LINE__);
    //         return(TNG_FAILURE);
        }
    }

    frame_set->n_mapping_blocks++;
    mappings = realloc(frame_set->mappings,
                       sizeof(struct tng_particle_mapping) *
                       frame_set->n_mapping_blocks);
    if(!mappings)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        free(frame_set->mappings);
        return(TNG_CRITICAL);
    }
    frame_set->mappings = mappings;
    mapping = &mappings[frame_set->n_mapping_blocks - 1];


    memcpy(&mapping->num_first_particle, block->block_contents+offset,
           sizeof(mapping->num_first_particle));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &mapping->num_first_particle)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(mapping->num_first_particle);

    memcpy(&mapping->n_particles, block->block_contents+offset,
           sizeof(mapping->n_particles));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &mapping->n_particles)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(mapping->n_particles);

    mapping->real_particle_numbers = malloc(mapping->n_particles *
                                            sizeof(int64_t));
    if(!mapping->real_particle_numbers)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                mapping->n_particles * sizeof(int64_t), __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* If the byte order needs to be swapped the data must be read one value at
     * a time and swapped */
    if(tng_data->input_endianness_swap_func_64)
    {
        for(i = 0; i < mapping->n_particles; i++)
        {
            memcpy(&mapping->real_particle_numbers[i],
                    block->block_contents + offset,
                    sizeof(int64_t));
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                            &mapping->real_particle_numbers[i])
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
            offset += sizeof(int64_t);
        }
    }
    /* Otherwise the data can be read all at once */
    else
    {
        memcpy(mapping->real_particle_numbers, block->block_contents + offset,
               mapping->n_particles * sizeof(int64_t));
    }


    return(TNG_SUCCESS);
}

/** Write the atom mappings of the current trajectory frame set
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @param mapping_block_nr is the index of the mapping block to write.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH an md5 hash will be generated and written.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
static tng_function_status tng_trajectory_mapping_block_write
                (tng_trajectory_t tng_data,
                 tng_gen_block_t block,
                 int mapping_block_nr,
                 const tng_hash_mode hash_mode)
{
    char *temp_name;
    int i, offset = 0, name_len;
    tng_particle_mapping_t mapping =
    &tng_data->current_trajectory_frame_set.mappings[mapping_block_nr];

    if(mapping_block_nr >=
       tng_data->current_trajectory_frame_set.n_mapping_blocks)
    {
        printf("Mapping block index out of bounds. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    name_len = strlen("PARTICLE MAPPING");

    if(!block->name || strlen(block->name) < name_len)
    {
        temp_name = realloc(block->name, name_len + 1);
        if(!temp_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                   name_len+1, __FILE__, __LINE__);
            free(block->name);
            return(TNG_CRITICAL);
        }
        block->name = temp_name;
    }
    strcpy(block->name, "PARTICLE MAPPING");
    block->id = TNG_PARTICLE_MAPPING;

    block->block_contents_size = sizeof(int64_t) * (2 + mapping->n_particles);

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    memcpy(block->block_contents, &mapping->num_first_particle,
           sizeof(mapping->num_first_particle));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(mapping->num_first_particle);

    memcpy(block->block_contents+offset, &mapping->n_particles,
           sizeof(mapping->n_particles));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(mapping->n_particles);

    if(tng_data->output_endianness_swap_func_64)
    {
        for(i = 0; i < mapping->n_particles; i++)
        {
            memcpy(block->block_contents+offset, &mapping->real_particle_numbers[i],
                sizeof(int64_t));
            if(tng_data->output_endianness_swap_func_64(tng_data,
                                        (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
            offset += sizeof(int64_t);
        }
    }
    else
    {
        memcpy(block->block_contents+offset, mapping->real_particle_numbers,
               mapping->n_particles * sizeof(int64_t));
    }


    if(tng_block_header_write(tng_data, block, hash_mode) != TNG_SUCCESS)
    {
        printf("Cannot write header of file %s. %s: %d\n",
               tng_data->output_file_path, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
              tng_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    return(TNG_SUCCESS);
}

/** Prepare a block for storing particle data
 * @param tng_data is a trajectory data container.
 * @param block_type_flag specifies if this is a trajectory block or a
 * non-trajectory block. (TNG_TRAJECTORY_BLOCK or TNG_NON_TRAJECTORY_BLOCK)
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_particle_data_block_create
                (tng_trajectory_t tng_data,
                 const tng_block_type block_type_flag)
{
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;

    tng_particle_data_t data;

    if(block_type_flag == TNG_TRAJECTORY_BLOCK)
    {
        frame_set->n_particle_data_blocks++;
        data = realloc(frame_set->tr_particle_data,
                    sizeof(struct tng_particle_data) *
                    frame_set->n_particle_data_blocks);
        if(!data)
        {
            printf("Cannot allocate memory (%lu bytes). %s: %d\n",
                sizeof(struct tng_particle_data) *
                frame_set->n_particle_data_blocks,
                __FILE__, __LINE__);
            free(frame_set->tr_particle_data);
            return(TNG_CRITICAL);
        }
        frame_set->tr_particle_data = data;
        data = &frame_set->tr_particle_data[frame_set->
                                            n_particle_data_blocks - 1];
    }
    else
    {
        tng_data->n_particle_data_blocks++;
        data = realloc(tng_data->non_tr_particle_data,
                        sizeof(struct tng_particle_data) *
                        tng_data->n_particle_data_blocks);
        if(!data)
        {
            printf("Cannot allocate memory (%lu bytes). %s: %d\n",
                    sizeof(struct tng_particle_data) *
                    tng_data->n_particle_data_blocks,
                    __FILE__, __LINE__);
            free(tng_data->non_tr_particle_data);
            return(TNG_CRITICAL);
        }
        tng_data->non_tr_particle_data = data;
        data = &tng_data->non_tr_particle_data[tng_data->
                                               n_particle_data_blocks - 1];
    }

    return(TNG_SUCCESS);
}

#ifdef USE_ZLIB
static tng_function_status tng_gzip_compress(tng_trajectory_t tng_data,
                                             tng_gen_block_t block,
                                             void *start_pos, int len)
{
    Bytef *dest;
    char *temp;
    uLong max_len, stat;

    max_len = compressBound(len);
    dest = malloc(max_len);

    stat = compress(dest, &max_len, start_pos, len);
    if(stat != Z_OK)
    {
        free(dest);
        if(stat == Z_MEM_ERROR)
        {
            printf("Not enough memory. ");
        }
        else if(stat == Z_BUF_ERROR)
        {
            printf("Destination buffer too small. ");
        }
        printf("Error gzipping data. %s: %d\n", __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    memcpy(start_pos, dest, max_len);

    block->block_contents_size = max_len + (block->block_contents_size - len);

    temp = realloc(block->block_contents, block->block_contents_size);
    if(!temp)
    {
        free(block->block_contents);
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    block->block_contents = temp;

    return(TNG_SUCCESS);
}

static tng_function_status tng_gzip_uncompress(tng_trajectory_t tng_data,
                                               tng_gen_block_t block,
                                               void *start_pos,
                                               unsigned long uncompressed_len)
{
    Bytef *dest;
    char *temp;
    unsigned long stat;
    int offset;

    offset = start_pos - (void *)block->block_contents;

    dest = malloc(uncompressed_len);

    stat = uncompress(dest, &uncompressed_len, (Bytef *) start_pos,
                      block->block_contents_size - offset);

    if(stat != Z_OK)
    {
        free(dest);
        if(stat == Z_MEM_ERROR)
        {
            printf("Not enough memory. ");
        }
        else if(stat == Z_BUF_ERROR)
        {
            printf("Destination buffer too small. ");
        }
        else if(stat == Z_DATA_ERROR)
        {
            printf("Data corrupt. ");
        }
        printf("Error uncompressing gzipped data. %s: %d\n", __FILE__,
               __LINE__);
        return(TNG_FAILURE);
    }


    block->block_contents_size = uncompressed_len + offset;

    temp = realloc(block->block_contents, uncompressed_len + offset);
    if(!temp)
    {
        free(block->block_contents);
        free(dest);
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    memcpy(temp + offset, dest, uncompressed_len);

    block->block_contents = temp;

    free(dest);
    return(TNG_SUCCESS);
}
#endif

/** Allocate memory for storing particle data.
 * The allocated block will be refered to by data->values.
 * @param tng_data is a trajectory data container.
 * @param data is the data struct, which will contain the allocated memory in
 * data->values.
 * @param n_frames is the number of frames of data to store.
 * @param n_particles is the number of particles with data.
 * @param n_values_per_frame is the number of data values per particle and
 * frame.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_allocate_particle_data_mem
                (tng_trajectory_t tng_data,
                 tng_particle_data_t data,
                 int64_t n_frames,
                 int64_t stride_length,
                 const int64_t n_particles,
                 const int64_t n_values_per_frame)
{
    union data_values ***values;
    int64_t i, j, k;

    if(data->values)
    {
        for(i = data->n_frames; i--;)
        {
            for(j = n_particles; j--;)
            {
                free(data->values[i][j]);
            }
            free(data->values[i]);
        }
    }
    data->n_frames = n_frames;
    n_frames = tng_max(1, n_frames);
    data->stride_length = tng_max(1, stride_length);
    data->n_values_per_frame = n_values_per_frame;
    values = realloc(data->values, sizeof(union data_values **) * n_frames);

    if(!values)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
            sizeof(union data_values **) * n_frames,
            __FILE__, __LINE__);
        free(data->values);
        return(TNG_CRITICAL);
    }
    data->values = values;

    for(i = n_frames; i-- ;)
    {
        data->values[i] = malloc(sizeof(union data_values *) *
                                 n_particles);
        if(!data->values[i])
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                sizeof(union data_values *) * n_particles,
                __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        for(j = n_particles; j--;)
        {
            data->values[i][j] = malloc(sizeof(union data_values) *
                                        n_values_per_frame);
            if(!data->values[i][j])
            {
                printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                    sizeof(union data_values) * n_values_per_frame,
                    __FILE__, __LINE__);
                return(TNG_CRITICAL);
            }
            if(data->datatype == TNG_CHAR_DATA)
            {
                for(k = n_values_per_frame; k--;)
                {
                    data->values[i][j][k].c = 0;
                }
            }
        }
    }
    return(TNG_SUCCESS);
}

/** Read the values of a particle data block
 * @param tng_data is a trajectory data container.
 * @param block is the block to store the data (should already contain
 * the block headers and the block contents).
 * @param offset is the reading offset to point at the place where the actual
 * values are stored, starting from the beginning of the block_contents. The
 * offset is changed during the reading.
 * @param datatype is the type of data of the data block (char, int, float or
 * double).
 * @param num_first_particle is the number of the first particle in the data
 * block. This should be the same as in the corresponding particle mapping
 * block.
 * @param n_particles is the number of particles in the data block. This should
 * be the same as in the corresponding particle mapping block.
 * @param first_frame_with_data is the frame number of the first frame with data
 * in this data block.
 * @param stride_length is the number of frames between each data entry.
 * @param n_frames is the number of frames in this data block.
 * @param n_values is the number of values per particle and frame stored in this
 * data block.
 * @param codec_id is the ID of the codec to compress the data.
 * @param multiplier is the multiplication factor applied to each data value
 * before compression. This factor is applied since some compression algorithms
 * work only on integers.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_particle_data_read
                (tng_trajectory_t tng_data,
                 tng_gen_block_t block,
                 int *offset,
                 const char datatype,
                 const int64_t num_first_particle,
                 const int64_t n_particles,
                 const int64_t first_frame_with_data,
                 const int64_t stride_length,
                 int64_t n_frames,
                 const int64_t n_values,
                 const int64_t codec_id,
                 const int64_t multiplier)
{
    int64_t block_index, i, j, k, tot_n_particles;
    int size, len;
    unsigned long data_size;
    union data_values **first_dim_values, *second_dim_values;
    tng_particle_data_t data;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;
    tng_block_type block_type_flag;

    switch(datatype)
    {
    case TNG_CHAR_DATA:
        size = 1;
        break;
    case TNG_INT_DATA:
        size = sizeof(int64_t);
        break;
    case TNG_FLOAT_DATA:
        size = sizeof(float);
        break;
    case TNG_DOUBLE_DATA:
    default:
        size = sizeof(double);
    }

    if(tng_data->current_trajectory_frame_set_input_file_pos > 0)
    {
        block_type_flag = TNG_TRAJECTORY_BLOCK;
    }
    else
    {
        block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
    }

    block_index = -1;
    /* See if there is already a data block of this ID */
    if(block_type_flag == TNG_TRAJECTORY_BLOCK)
    {
        for(i = frame_set->n_particle_data_blocks; i-- ;)
        {
            data = &frame_set->tr_particle_data[i];
            if(data->block_id == block->id)
            {
                block_index = i;
                break;
            }
        }
    }
    else
    {
        for(i = tng_data->n_particle_data_blocks; i-- ;)
        {
            data = &tng_data->non_tr_particle_data[i];
            if(data->block_id == block->id)
            {
                block_index = i;
                break;
            }
        }
    }

    /* Otherwise create a data block */
    if(block_index == -1)
    {
        if(tng_particle_data_block_create(tng_data, block_type_flag) !=
            TNG_SUCCESS)
        {
            printf("Cannot create particle data block. %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        if(block_type_flag == TNG_TRAJECTORY_BLOCK)
        {
            data = &frame_set->tr_particle_data[frame_set->
                                                n_particle_data_blocks - 1];
        }
        else
        {
            data = &tng_data->non_tr_particle_data[tng_data->
                                                   n_particle_data_blocks - 1];
        }
        data->block_id = block->id;

        data->block_name = malloc(strlen(block->name) + 1);
        if(!data->block_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                   (int)strlen(block->name)+1, __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        strcpy(data->block_name, block->name);

        data->datatype = datatype;

        data->values = 0;
        data->n_frames = 0;
        data->codec_id = codec_id;
        data->compression_multiplier = multiplier;
    }

    if(block_type_flag == TNG_TRAJECTORY_BLOCK &&
       tng_data->var_num_atoms_flag)
    {
        tot_n_particles = frame_set->n_particles;
    }
    else
    {
        tot_n_particles = tng_data->n_particles;
    }

    if(codec_id != TNG_UNCOMPRESSED)
    {
        data_size = (n_frames / stride_length) * size * n_particles * n_values;
        switch(codec_id)
        {
        case TNG_XTC_COMPRESSION:
            printf("XTC compression not implemented yet.\n");
            break;
        case TNG_TNG_COMPRESSION:
            printf("TNG compression not implemented yet.\n");
            break;
    #ifdef USE_ZLIB
        case TNG_GZIP_COMPRESSION:
//             printf("Before uncompression: %"PRId64"\n", block->block_contents_size);
            if(tng_gzip_uncompress(tng_data, block,
                                   block->block_contents + *offset,
                                   data_size) != TNG_SUCCESS)
            {
                printf("Could not read gzipped block data. %s: %d\n", __FILE__,
                    __LINE__);
                return(TNG_CRITICAL);
            }
//             printf("After uncompression: %"PRId64"\n", block->block_contents_size);
            break;
    #endif
        }
    }

    /* Allocate memory */
    if(!data->values || data->n_frames != n_frames ||
       data->n_values_per_frame != n_values)
    {
        if(tng_allocate_particle_data_mem(tng_data, data, n_frames,
                                          stride_length,
                                          tot_n_particles, n_values) !=
           TNG_SUCCESS)
        {
            printf("Cannot allocate memory for particle data. %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    data->first_frame_with_data = first_frame_with_data;

    n_frames = tng_max(1, n_frames / stride_length);

    /* FIXME: If not using a union to store data a whole dimension
     * or the whole block can be read at once if byte swapping is not needed */
    switch(datatype)
    {
    case TNG_FLOAT_DATA:
        for(i = 0; i < n_frames; i++)
        {
            first_dim_values = data->values[i];
            for(j = num_first_particle; j < num_first_particle + n_particles;
                j++)
            {
                second_dim_values = first_dim_values[j];
                for(k = 0; k < n_values; k++)
                {
                    memcpy(&second_dim_values[k].f,
                        block->block_contents+*offset,
                        size);
                    if(tng_data->input_endianness_swap_func_32)
                    {
                        if(tng_data->input_endianness_swap_func_32(tng_data,
                           (int32_t *)&second_dim_values[k])
                            != TNG_SUCCESS)
                        {
                            printf("Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }
                    *offset += size;
                }
            }
        }
        break;
    case TNG_INT_DATA:
        for(i = 0; i < n_frames; i++)
        {
            first_dim_values = data->values[i];
            for(j = num_first_particle; j < num_first_particle + n_particles;
                j++)
            {
                second_dim_values = first_dim_values[j];
                for(k = 0; k < n_values; k++)
                {
                    memcpy(&second_dim_values[k].i,
                           block->block_contents+*offset,
                           size);
                    if(tng_data->input_endianness_swap_func_64)
                    {
                        if(tng_data->input_endianness_swap_func_64(tng_data,
                           (int64_t *)&second_dim_values[k])
                            != TNG_SUCCESS)
                        {
                            printf("Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }
                    *offset += size;
                }
            }
        }
        break;
    case TNG_CHAR_DATA:
        for(i = 0; i < n_frames; i++)
        {
            first_dim_values = data->values[i];
            for(j = num_first_particle; j < num_first_particle + n_particles;
                j++)
            {
                second_dim_values = first_dim_values[j];
                for(k = 0; k < n_values; k++)
                {
                    len = tng_min(strlen(block->block_contents+*offset) + 1,
                              TNG_MAX_STR_LEN);
                    if(second_dim_values[k].c)
                    {
                        free(second_dim_values[k].c);
                    }
                    second_dim_values[k].c = malloc(len);
                    if(!second_dim_values[k].c)
                    {
                        printf("Cannot allocate memory (%d bytes). %s: %d\n",
                            len, __FILE__, __LINE__);
                        return(TNG_CRITICAL);
                    }
                    strncpy(second_dim_values[k].c,
                            block->block_contents+*offset, len);
                    *offset += len;
                }
            }
        }
        break;
    case TNG_DOUBLE_DATA:
    default:
        for(i = 0; i < n_frames; i++)
        {
            first_dim_values = data->values[i];
            for(j = num_first_particle; j < num_first_particle + n_particles;
                j++)
            {
                second_dim_values = first_dim_values[j];
                for(k = 0; k < n_values; k++)
                {
                    memcpy(&second_dim_values[k].d,
                           block->block_contents+*offset,
                           size);
                    if(tng_data->input_endianness_swap_func_64)
                    {
                        if(tng_data->input_endianness_swap_func_64(tng_data,
                           (int64_t *)&second_dim_values[k])
                            != TNG_SUCCESS)
                        {
                            printf("Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }
                    *offset += size;
                }
            }
        }
    }
    return(TNG_SUCCESS);
}

/** Write a particle data block
 * @param tng_data is a trajectory data container.
 * @param block is the block to store the data (should already contain
 * the block headers and the block contents).
 * @param block_index is the index number of the data block in the frame set.
 * @param mapping is the particle mapping that is relevant for the data block.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH an md5 hash will be generated and written.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_particle_data_block_write
                (tng_trajectory_t tng_data,
                 tng_gen_block_t block,
                 const int block_index,
                 const tng_particle_mapping_t mapping,
                 const tng_hash_mode hash_mode)
{
    int64_t n_particles, num_first_particle, n_frames, stride_length;
    int i, j, k, offset = 0, size, len, data_start_pos;
    char dependency, temp, *temp_name;
    double multiplier, d_temp;
    float f_temp;
    union data_values **first_dim_values, *second_dim_values;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;

    tng_particle_data_t data;
    tng_block_type block_type_flag;

    /* If we have already started writing frame sets it is too late to write
     * non-trajectory data blocks */
    if(tng_data->current_trajectory_frame_set_output_file_pos > 0)
    {
        block_type_flag = TNG_TRAJECTORY_BLOCK;
    }
    else
    {
        block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
    }

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    if(block_type_flag == TNG_TRAJECTORY_BLOCK)
    {
        data = &frame_set->tr_particle_data[block_index];
        stride_length = tng_max(1, data->stride_length);
    }
    else
    {
        data = &tng_data->non_tr_particle_data[block_index];
        stride_length = 1;
    }

    switch(data->datatype)
    {
    case TNG_CHAR_DATA:
        size = 1;
        break;
    case TNG_INT_DATA:
        size = sizeof(int64_t);
        break;
    case TNG_FLOAT_DATA:
        size = sizeof(float);
        break;
    case TNG_DOUBLE_DATA:
    default:
        size = sizeof(double);
    }

    len = strlen(data->block_name) + 1;

    if(!block->name || strlen(block->name) < len)
    {
        temp_name = realloc(block->name, len);
        if(!temp_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            free(block->name);
            return(TNG_CRITICAL);
        }
        block->name = temp_name;
    }
    strncpy(block->name, data->block_name, len);
    block->id = data->block_id;

    /* If writing frame independent data data->n_frames is 0, but n_frames
       is used for the loop writing the data (and reserving memory) and needs
       to be at least 1 */
    n_frames = tng_max(1, data->n_frames);


    if(mapping && mapping->n_particles != 0)
    {
        n_particles = mapping->n_particles;
        num_first_particle = mapping->num_first_particle;
    }
    else
    {
        num_first_particle = 0;
        if(tng_data->var_num_atoms_flag)
        {
            n_particles = frame_set->n_particles;
        }
        else
        {
            n_particles = tng_data->n_particles;
        }
    }

    block->block_contents_size = sizeof(char) * 2 +
                                 sizeof(data->n_values_per_frame) +
                                 sizeof(data->codec_id) +
                                 sizeof(num_first_particle) +
                                 sizeof(n_particles);

    if(stride_length > 1)
    {
        block->block_contents_size += sizeof(data->first_frame_with_data) +
                                      sizeof(data->stride_length);
    }

    if(data->codec_id != TNG_UNCOMPRESSED)
    {
        block->block_contents_size += sizeof(data->compression_multiplier);
    }

    if(block_type_flag == TNG_TRAJECTORY_BLOCK && data->n_frames > 0)
    {
        dependency = TNG_FRAME_DEPENDENT + TNG_PARTICLE_DEPENDENT;
    }
    else
    {
        dependency = TNG_PARTICLE_DEPENDENT;
    }
    if(dependency & TNG_FRAME_DEPENDENT)
    {
        block->block_contents_size += sizeof(char);
    }

    data_start_pos = block->block_contents_size;

    if(data->datatype == TNG_CHAR_DATA)
    {
        for(i = n_frames; i--;)
        {
            first_dim_values = data->values[i];
            for(j = num_first_particle; j < num_first_particle + n_particles;
                j++)
            {
                second_dim_values = first_dim_values[j];
                for(k = data->n_values_per_frame; k--;)
                {
                    block->block_contents_size +=
                    strlen(second_dim_values[k].c) + 1;
                }
            }
        }
    }
    else
    {
        block->block_contents_size += size * n_frames / stride_length *
                                      n_particles * data->n_values_per_frame;
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }


    memcpy(block->block_contents, &data->datatype, sizeof(char));
    offset += sizeof(char);

    memcpy(block->block_contents+offset, &dependency, sizeof(char));
    offset += sizeof(char);

    if(dependency & TNG_FRAME_DEPENDENT)
    {
        if(stride_length > 1)
        {
            temp = 1;
        }
        else
        {
            temp = 0;
        }
        memcpy(block->block_contents+offset, &temp, sizeof(char));
        offset += sizeof(char);
    }

    memcpy(block->block_contents+offset, &data->n_values_per_frame,
           sizeof(data->n_values_per_frame));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
           (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(data->n_values_per_frame);

    memcpy(block->block_contents+offset, &data->codec_id,
           sizeof(data->codec_id));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
           (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(data->codec_id);

    if(data->codec_id != TNG_UNCOMPRESSED)
    {
        memcpy(block->block_contents+offset, &data->compression_multiplier,
               sizeof(data->compression_multiplier));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
               (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->compression_multiplier);
    }

    if(data->n_frames > 0 && stride_length > 1)
    {
        /* FIXME: Currently only supporting
         * frame_set->n_frames / data->stride_length == 0 */
        data->first_frame_with_data = frame_set->first_frame;
        memcpy(block->block_contents+offset, &data->first_frame_with_data,
               sizeof(data->first_frame_with_data));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
               (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->first_frame_with_data);

        memcpy(block->block_contents+offset, &stride_length,
               sizeof(stride_length));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
               (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(stride_length);
    }


    memcpy(block->block_contents+offset, &num_first_particle,
           sizeof(num_first_particle));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
           (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(num_first_particle);

    memcpy(block->block_contents+offset, &n_particles, sizeof(n_particles));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
           (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(n_particles);


    if(data->values)
    {
        /* FIXME: If not using a union to store data a whole dimension or the
        * whole block can be written at once if byte swapping is not needed */
        switch(data->datatype)
        {
        case TNG_FLOAT_DATA:
            /* For speed reasons the compression multiplier is not used if the data
            * is not compressed. */
            if(data->codec_id == TNG_UNCOMPRESSED)
            {
                for(i = 0; i < data->n_frames / stride_length; i++)
                {
                    first_dim_values = data->values[i];
                    for(j = num_first_particle; j < num_first_particle + n_particles;
                        j++)
                    {
                        second_dim_values = first_dim_values[j];
                        for(k = 0; k < data->n_values_per_frame; k++)
                        {
                            memcpy(block->block_contents+offset,
                                &second_dim_values[k].f,
                                size);
                            if(tng_data->output_endianness_swap_func_32)
                            {
                                if(tng_data->output_endianness_swap_func_32(tng_data,
                                (int32_t *)block->header_contents+offset)
                                    != TNG_SUCCESS)
                                {
                                    printf("Cannot swap byte order. %s: %d\n",
                                            __FILE__, __LINE__);
                                }
                            }
                            offset += size;
                        }
                    }
                }
            }
            else
            {
                multiplier = data->compression_multiplier;
                for(i = 0; i < data->n_frames / stride_length; i++)
                {
                    first_dim_values = data->values[i];
                    for(j = num_first_particle; j < num_first_particle + n_particles;
                        j++)
                    {
                        second_dim_values = first_dim_values[j];
                        for(k = 0; k < data->n_values_per_frame; k++)
                        {
                            f_temp = second_dim_values[k].f * multiplier;
                            memcpy(block->block_contents+offset,
                                &f_temp, size);
                            if(tng_data->output_endianness_swap_func_32)
                            {
                                if(tng_data->output_endianness_swap_func_32(tng_data,
                                (int32_t *)block->header_contents+offset)
                                    != TNG_SUCCESS)
                                {
                                    printf("Cannot swap byte order. %s: %d\n",
                                            __FILE__, __LINE__);
                                }
                            }
                            offset += size;
                        }
                    }
                }
            }
            break;

        case TNG_INT_DATA:
            for(i = 0; i < data->n_frames / stride_length; i++)
            {
                first_dim_values = data->values[i];
                for(j = num_first_particle; j < num_first_particle + n_particles;
                    j++)
                {
                    second_dim_values = first_dim_values[j];
                    for(k = 0; k < data->n_values_per_frame; k++)
                    {
                        memcpy(block->block_contents+offset,
                            &second_dim_values[k].i,
                            size);
                        if(tng_data->output_endianness_swap_func_64)
                        {
                            if(tng_data->output_endianness_swap_func_64(tng_data,
                            (int64_t *)block->header_contents+offset)
                                != TNG_SUCCESS)
                            {
                                printf("Cannot swap byte order. %s: %d\n",
                                        __FILE__, __LINE__);
                            }
                        }
                        offset += size;
                    }
                }
            }
            break;
        case TNG_CHAR_DATA:
            for(i = 0; i < data->n_frames / stride_length; i++)
            {
                first_dim_values = data->values[i];
                for(j = num_first_particle; j < num_first_particle + n_particles;
                    j++)
                {
                    second_dim_values = first_dim_values[j];
                    for(k = 0; k < data->n_values_per_frame; k++)
                    {
                        len = strlen(second_dim_values[k].c) + 1;
                        strncpy(block->block_contents+offset,
                                second_dim_values[k].c, len);
                        offset += len;
                    }
                }
            }
            break;
        case TNG_DOUBLE_DATA:
        default:
            /* For speed reasons the compression multiplier is not used if the data
            * is not compressed.*/
            if(data->codec_id == TNG_UNCOMPRESSED)
            {
                for(i = 0; i < data->n_frames / stride_length; i++)
                {
                    first_dim_values = data->values[i];
                    for(j = num_first_particle; j < num_first_particle + n_particles;
                        j++)
                    {
                        second_dim_values = first_dim_values[j];
                        for(k = 0; k < data->n_values_per_frame; k++)
                        {
                            memcpy(block->block_contents+offset,
                                &second_dim_values[k].d,
                                size);
                            if(tng_data->output_endianness_swap_func_64)
                            {
                                if(tng_data->output_endianness_swap_func_64(tng_data,
                                (int64_t *)block->header_contents+offset)
                                    != TNG_SUCCESS)
                                {
                                    printf("Cannot swap byte order. %s: %d\n",
                                            __FILE__, __LINE__);
                                }
                            }
                            offset += size;
                        }
                    }
                }
            }
            else
            {
                multiplier = data->compression_multiplier;
                for(i = 0; i < data->n_frames / stride_length; i++)
                {
                    first_dim_values = data->values[i];
                    for(j = num_first_particle; j < num_first_particle + n_particles;
                        j++)
                    {
                        second_dim_values = first_dim_values[j];
                        for(k = 0; k < data->n_values_per_frame; k++)
                        {
                            d_temp = second_dim_values[k].d * multiplier;
                            memcpy(block->block_contents+offset,
                                &d_temp,
                                size);
                            if(tng_data->output_endianness_swap_func_64)
                            {
                                if(tng_data->output_endianness_swap_func_64(tng_data,
                                (int64_t *)block->header_contents+offset)
                                    != TNG_SUCCESS)
                                {
                                    printf("Cannot swap byte order. %s: %d\n",
                                            __FILE__, __LINE__);
                                }
                            }
                            offset += size;
                        }
                    }
                }
            }
        }
    }
    else
    {
        memset(block->block_contents+offset, 0, block->block_contents_size - offset);
    }


    switch(data->codec_id)
    {
    case TNG_XTC_COMPRESSION:
        printf("XTC compression not implemented yet.\n");
        break;
    case TNG_TNG_COMPRESSION:
        printf("TNG compression not implemented yet.\n");
        break;
#ifdef USE_ZLIB
    case TNG_GZIP_COMPRESSION:
//         printf("Before compression: %"PRId64"\n", block->block_contents_size);
        if(tng_gzip_compress(tng_data, block,
                             block->block_contents + data_start_pos,
                             block->block_contents_size - data_start_pos) !=
           TNG_SUCCESS)
        {
            printf("Could not write gzipped block data. %s: %d\n", __FILE__,
                   __LINE__);
            return(TNG_CRITICAL);
        }
//         printf("After compression: %"PRId64"\n", block->block_contents_size);
        break;
#endif
    }

    if(tng_block_header_write(tng_data, block, hash_mode) != TNG_SUCCESS)
    {
        printf("Cannot write header of file %s. %s: %d\n",
               tng_data->output_file_path, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
        tng_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n", __FILE__,
                __LINE__);
        return(TNG_CRITICAL);
    }

    return(TNG_SUCCESS);
}

/* TEST: */
/** Create a non-particle data block
 * @param tng_data is a trajectory data container.
 * @param block_type_flag specifies if this is a trajectory block or a
 * non-trajectory block. (TNG_TRAJECTORY_BLOCK or TNG_NON_TRAJECTORY_BLOCK)
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_data_block_create
                (tng_trajectory_t tng_data,
                 const tng_block_type block_type_flag)
{
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;

    tng_non_particle_data_t data;

    if(block_type_flag == TNG_TRAJECTORY_BLOCK)
    {
        frame_set->n_data_blocks++;
        data = realloc(frame_set->tr_data, sizeof(struct tng_non_particle_data) *
                       frame_set->n_data_blocks);
        if(!data)
        {
            printf("Cannot allocate memory (%lu bytes). %s: %d\n",
                sizeof(struct tng_non_particle_data) * frame_set->n_data_blocks,
                __FILE__, __LINE__);
            free(frame_set->tr_data);
            return(TNG_CRITICAL);
        }
        frame_set->tr_data = data;
        data = &frame_set->tr_data[frame_set->n_data_blocks - 1];
    }
    else
    {
        tng_data->n_data_blocks++;
        data = realloc(tng_data->non_tr_data, sizeof(struct tng_non_particle_data) *
                        tng_data->n_data_blocks);
        if(!data)
        {
            printf("Cannot allocate memory (%lu bytes). %s: %d\n",
                sizeof(struct tng_non_particle_data) * tng_data->n_data_blocks,
                __FILE__, __LINE__);
            free(tng_data->non_tr_data);
            return(TNG_CRITICAL);
        }
        tng_data->non_tr_data = data;
        data = &tng_data->non_tr_data[tng_data->n_data_blocks - 1];
    }

    return(TNG_SUCCESS);
}

/* TEST: */
/** Allocate memory for storing non-particle data.
 * The allocated block will be refered to by data->values.
 * @param tng_data is a trajectory data container.
 * @param data is the data struct, which will contain the allocated memory in
 * data->values.
 * @param n_frames is the number of frames of data to store.
 * @param n_values_per_frame is the number of data values per frame.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_allocate_data_mem
                (tng_trajectory_t tng_data,
                 tng_non_particle_data_t data,
                 int64_t n_frames,
                 int64_t stride_length,
                 const int64_t n_values_per_frame)
{
    union data_values **values;
    int64_t i, j;

    if(data->values)
    {
        for(i = data->n_frames; i--;)
        {
            if(data->datatype == TNG_CHAR_DATA)
            {
                for(j = data->n_values_per_frame; j--;)
                {
                    if(data->values[i][j].c)
                    {
                        free(data->values[i][j].c);
                        data->values[i][j].c = 0;
                    }
                }
            }
            free(data->values[i]);
        }
    }
    data->n_frames = n_frames;
    data->stride_length = tng_max(1, stride_length);
    n_frames = tng_max(1, n_frames);
    data->n_values_per_frame = n_values_per_frame;
    values = realloc(data->values,
                     sizeof(union data_values *) *
                     n_frames);
    if(!values)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
            sizeof(union data_values *) * n_frames,
            __FILE__, __LINE__);
        free(data->values);
        return(TNG_CRITICAL);
    }
    data->values = values;

    for(i = n_frames; i-- ;)
    {
        data->values[i] = malloc(sizeof(union data_values) *
                                 n_values_per_frame);
        if(!data->values[i])
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                sizeof(union data_values) * n_values_per_frame,
                __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        if(data->datatype == TNG_CHAR_DATA)
        {
            for(j = n_values_per_frame; j--;)
            {
                data->values[i][j].c = 0;
            }
        }
    }
    return(TNG_SUCCESS);
}

/** Read the values of a non-particle data block
 * @param tng_data is a trajectory data container.
 * @param block is the block to store the data (should already contain
 * the block headers and the block contents).
 * @param offset is the reading offset to point at the place where the actual
 * values are stored, starting from the beginning of the block_contents. The
 * offset is changed during the reading.
 * @param datatype is the type of data of the data block (char, int, float or
 * double).
 * @param first_frame_with_data is the frame number of the first frame with data
 * in this data block.
 * @param stride_length is the number of frames between each data entry.
 * @param n_frames is the number of frames in this data block.
 * @param n_values is the number of values per frame stored in this data block.
 * @param codec_id is the ID of the codec to compress the data.
 * @param multiplier is the multiplication factor applied to each data value
 * before compression. This factor is applied since some compression algorithms
 * work only on integers.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_data_read(tng_trajectory_t tng_data,
                                         tng_gen_block_t block,
                                         int *offset,
                                         const char datatype,
                                         const int64_t first_frame_with_data,
                                         const int64_t stride_length,
                                         int64_t n_frames,
                                         const int64_t n_values,
                                         const int64_t codec_id,
                                         const int64_t multiplier)
{
    int64_t block_index, i, j;
    int size, len;
    unsigned long data_size;
    tng_non_particle_data_t data;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;
    tng_block_type block_type_flag;

//     printf("%s\n", block->name);

    if(tng_data->current_trajectory_frame_set_input_file_pos > 0)
    {
        block_type_flag = TNG_TRAJECTORY_BLOCK;
    }
    else
    {
        block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
    }

    switch(datatype)
    {
    case TNG_CHAR_DATA:
        size = 1;
        break;
    case TNG_INT_DATA:
        size = sizeof(int64_t);
        break;
    case TNG_FLOAT_DATA:
        size = sizeof(float);
        break;
    case TNG_DOUBLE_DATA:
    default:
        size = sizeof(double);
    }

    block_index = -1;
    /* See if there is already a data block of this ID */
    /* FIXME: Do not compare with block->id. Use ID as parameter instead. */
    if(block_type_flag == TNG_TRAJECTORY_BLOCK)
    {
        for(i = frame_set->n_data_blocks; i-- ;)
        {
            data = &frame_set->tr_data[i];
            if(data->block_id == block->id)
            {
                block_index = i;
                break;
            }
        }
    }
    else
    {
        for(i = tng_data->n_data_blocks; i-- ;)
        {
            data = &tng_data->non_tr_data[i];
            if(data->block_id == block->id)
            {
                block_index = i;
                break;
            }
        }
    }

    /* Otherwise create a data block */
    if(block_index == -1)
    {
        if(tng_data_block_create(tng_data, block_type_flag) !=
            TNG_SUCCESS)
        {
            printf("Cannot create particle data block. %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        if(block_type_flag == TNG_TRAJECTORY_BLOCK)
        {
            data = &frame_set->tr_data[frame_set->n_data_blocks - 1];
        }
        else
        {
            data = &tng_data->non_tr_data[tng_data->n_data_blocks - 1];
        }
        data->block_id = block->id;

        data->block_name = malloc(strlen(block->name) + 1);
        if(!data->block_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                   (int)strlen(block->name)+1, __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        strcpy(data->block_name, block->name);

        data->datatype = datatype;

        data->values = 0;
        data->n_frames = 0;
        data->codec_id = codec_id;
        data->compression_multiplier = multiplier;
    }

    if(codec_id != TNG_UNCOMPRESSED)
    {
        data_size = (n_frames / stride_length) * size * n_values;
        switch(codec_id)
        {
        case TNG_XTC_COMPRESSION:
            printf("XTC compression not implemented yet.\n");
            break;
        case TNG_TNG_COMPRESSION:
            printf("TNG compression not implemented yet.\n");
            break;
    #ifdef USE_ZLIB
        case TNG_GZIP_COMPRESSION:
    //         printf("Before compression: %"PRId64"\n", block->block_contents_size);
            if(tng_gzip_uncompress(tng_data, block,
                                   block->block_contents + *offset,
                                   data_size) != TNG_SUCCESS)
            {
                printf("Could not read gzipped block data. %s: %d\n", __FILE__,
                    __LINE__);
                return(TNG_CRITICAL);
            }
    //         printf("After compression: %"PRId64"\n", block->block_contents_size);
            break;
    #endif
        }
    }

    /* Allocate memory */
    if(!data->values || data->n_frames != n_frames ||
       data->n_values_per_frame != n_values)
    {
        if(tng_allocate_data_mem(tng_data, data, n_frames, stride_length,
                                 n_values) !=
           TNG_SUCCESS)
        {
            printf("Cannot allocate memory for data. %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    data->first_frame_with_data = first_frame_with_data;

    n_frames = tng_max(1, n_frames / stride_length);

    /* FIXME: If not using a union to store data a whole dimension
     * or the whole block can be read at once if byte swapping is not needed */
    switch(datatype)
    {
    case TNG_FLOAT_DATA:
        for(i = 0; i < n_frames; i++)
        {
            for(j = 0; j < n_values; j++)
            {
                memcpy(&data->values[i][j].f, block->block_contents+*offset,
                       size);
                if(tng_data->input_endianness_swap_func_32)
                {
                    if(tng_data->input_endianness_swap_func_32(tng_data,
                        (int32_t *)&data->values[i][j])
                        != TNG_SUCCESS)
                    {
                        printf("Cannot swap byte order. %s: %d\n",
                                __FILE__, __LINE__);
                    }
                }
                *offset += size;
            }
        }
        break;
    case TNG_INT_DATA:
        for(i = 0; i < n_frames; i++)
        {
            for(j = 0; j < n_values; j++)
            {
                memcpy(&data->values[i][j].i, block->block_contents+*offset,
                        size);
                if(tng_data->input_endianness_swap_func_64)
                {
                    if(tng_data->input_endianness_swap_func_64(tng_data,
                        (int64_t *)&data->values[i][j])
                        != TNG_SUCCESS)
                    {
                        printf("Cannot swap byte order. %s: %d\n",
                                __FILE__, __LINE__);
                    }
                }
                *offset += size;
            }
        }
        break;
    case TNG_CHAR_DATA:
        for(i = 0; i < n_frames; i++)
        {
            for(j = 0; j < n_values; j++)
            {
                len = tng_min(strlen(block->block_contents+*offset) + 1,
                          TNG_MAX_STR_LEN);
                if(data->values[i][j].c)
                {
                    free(data->values[i][j].c);
                }
                data->values[i][j].c = malloc(len);
                if(!data->values[i][j].c)
                {
                    printf("Cannot allocate memory (%d bytes). %s: %d\n",
                           len, __FILE__, __LINE__);
                    return(TNG_CRITICAL);
                }
                strncpy(data->values[i][j].c, block->block_contents+*offset,
                        len);
                *offset += len;
            }
        }
        break;
    case TNG_DOUBLE_DATA:
    default:
        for(i = 0; i < n_frames; i++)
        {
            for(j = 0; j < n_values; j++)
            {
                memcpy(&data->values[i][j].d, block->block_contents+*offset,
                       size);
                if(tng_data->input_endianness_swap_func_64)
                {
                    if(tng_data->input_endianness_swap_func_64(tng_data,
                        (int64_t *)&data->values[i][j])
                        != TNG_SUCCESS)
                    {
                        printf("Cannot swap byte order. %s: %d\n",
                                __FILE__, __LINE__);
                    }
                }
                *offset += size;
            }
        }
    }
    return(TNG_SUCCESS);
}

/** Write a non-particle data block
 * @param tng_data is a trajectory data container.
 * @param block is the block to store the data (should already contain
 * the block headers and the block contents).
 * @param block_index is the index number of the data block in the frame set.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH an md5 hash will be generated and written.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_data_block_write(tng_trajectory_t tng_data,
                                                tng_gen_block_t block,
                                                const int block_index,
                                                const tng_hash_mode hash_mode)
{
    int64_t n_frames, stride_length;
    int i, j, offset = 0, size, len, data_start_pos;
    char temp, dependency, *temp_name;
    double multiplier, d_temp;
    float f_temp;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;

    tng_non_particle_data_t data;
    tng_block_type block_type_flag;

    /* If we have already started writing frame sets it is too late to write
     * non-trajectory data blocks */
    if(tng_data->current_trajectory_frame_set_output_file_pos > 0)
    {
        block_type_flag = TNG_TRAJECTORY_BLOCK;
    }
    else
    {
        block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
    }

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    if(block_type_flag == TNG_TRAJECTORY_BLOCK)
    {
        data = &frame_set->tr_data[block_index];
        stride_length = tng_max(1, data->stride_length);
    }
    else
    {
        data = &tng_data->non_tr_data[block_index];
        stride_length = 1;
    }

    switch(data->datatype)
    {
    case TNG_CHAR_DATA:
        size = 1;
        break;
    case TNG_INT_DATA:
        size = sizeof(int64_t);
        break;
    case TNG_FLOAT_DATA:
        size = sizeof(float);
        break;
    case TNG_DOUBLE_DATA:
    default:
        size = sizeof(double);
    }

    len = strlen(data->block_name) + 1;

    if(!block->name || strlen(block->name) < len)
    {
        temp_name = realloc(block->name, len);
        if(!temp_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len+1,
                   __FILE__, __LINE__);
            free(block->name);
            return(TNG_CRITICAL);
        }
        block->name = temp_name;
    }
    strncpy(block->name, data->block_name, len);
    block->id = data->block_id;

    /* If writing frame independent data data->n_frames is be 0, but n_frames
       is used for the loop writing the data (and reserving memory) and needs
       to be at least 1 */
    n_frames = tng_max(1, data->n_frames);


    block->block_contents_size = sizeof(char) * 2 +
                                 sizeof(data->n_values_per_frame) +
                                 sizeof(data->codec_id);

    if(stride_length > 1)
    {
        block->block_contents_size += sizeof(data->first_frame_with_data) +
                                      sizeof(data->stride_length);
    }

    if(data->codec_id != TNG_UNCOMPRESSED)
    {
        block->block_contents_size += sizeof(data->compression_multiplier);
    }

    data_start_pos = block->block_contents_size;

    if(data->datatype == TNG_CHAR_DATA)
    {
        for(i = n_frames; i--;)
        {
            for(j = data->n_values_per_frame; j--;)
            {
                block->block_contents_size += strlen(data->values[i][j].c) + 1;
            }
        }
    }
    else
    {
        block->block_contents_size += size * n_frames / stride_length *
        data->n_values_per_frame;
    }

    if(block_type_flag == TNG_TRAJECTORY_BLOCK && data->n_frames > 0)
    {
        dependency = TNG_FRAME_DEPENDENT;
    }
    else
    {
        dependency = 0;
    }
    if(dependency & TNG_FRAME_DEPENDENT)
    {
        block->block_contents_size += sizeof(char);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }


    memcpy(block->block_contents, &data->datatype, sizeof(char));
    offset += sizeof(char);

    memcpy(block->block_contents+offset, &dependency, sizeof(char));
    offset += sizeof(char);

    if(dependency & TNG_FRAME_DEPENDENT)
    {
        if(stride_length > 1)
        {
            temp = 1;
        }
        else
        {
            temp = 0;
        }
        memcpy(block->block_contents+offset, &temp, sizeof(char));
        offset += sizeof(char);
    }

    memcpy(block->block_contents+offset, &data->n_values_per_frame,
           sizeof(data->n_values_per_frame));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
           (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(data->n_values_per_frame);

    memcpy(block->block_contents+offset, &data->codec_id,
           sizeof(data->codec_id));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
           (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(data->codec_id);

    if(data->codec_id != TNG_UNCOMPRESSED)
    {
        memcpy(block->block_contents+offset, &data->compression_multiplier,
               sizeof(data->compression_multiplier));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
            (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->compression_multiplier);
    }

    if(data->n_frames > 0 && stride_length > 1)
    {
        /* FIXME: Currently only supporting
         * frame_set->n_frames / data->stride_length == 0 */
        data->first_frame_with_data = frame_set->first_frame;
        memcpy(block->block_contents+offset, &data->first_frame_with_data,
               sizeof(data->first_frame_with_data));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
            (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->first_frame_with_data);

        memcpy(block->block_contents+offset, &stride_length,
               sizeof(data->stride_length));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
            (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->stride_length);
    }

    if(data->values)
    {
        /* FIXME: If not using a union to store data a whole dimension or the
        * whole block can be written at once if byte swapping is not needed */
        switch(data->datatype)
        {
        case TNG_FLOAT_DATA:
            /* For speed reasons the compression multiplier is not used if the data
            * is not compressed.*/
            if(data->codec_id == TNG_UNCOMPRESSED)
            {
                for(i = 0; i < n_frames / stride_length; i++)
                {
                    for(j = 0; j < data->n_values_per_frame; j++)
                    {
                        memcpy(block->block_contents+offset, &data->values[i][j].f,
                            size);
                        if(tng_data->output_endianness_swap_func_32)
                        {
                            if(tng_data->output_endianness_swap_func_32(tng_data,
                            (int32_t *)block->header_contents+offset)
                                != TNG_SUCCESS)
                            {
                                printf("Cannot swap byte order. %s: %d\n",
                                        __FILE__, __LINE__);
                            }
                        }
                        offset += size;
                    }
                }
            }
            else
            {
                multiplier = data->compression_multiplier;
                for(i = 0; i < n_frames / stride_length; i++)
                {
                    for(j = 0; j < data->n_values_per_frame; j++)
                    {
                        f_temp = data->values[i][j].f * multiplier;
                        memcpy(block->block_contents+offset, &f_temp,
                            size);
                        if(tng_data->output_endianness_swap_func_32)
                        {
                            if(tng_data->output_endianness_swap_func_32(tng_data,
                            (int32_t *)block->header_contents+offset)
                                != TNG_SUCCESS)
                            {
                                printf("Cannot swap byte order. %s: %d\n",
                                        __FILE__, __LINE__);
                            }
                        }
                        offset += size;
                    }
                }
            }
            break;
        case TNG_INT_DATA:
            for(i = 0; i < n_frames / stride_length; i++)
            {
                for(j = 0; j < data->n_values_per_frame; j++)
                {
                    memcpy(block->block_contents+offset, &data->values[i][j].i,
                        size);
                    if(tng_data->output_endianness_swap_func_64)
                    {
                        if(tng_data->output_endianness_swap_func_64(tng_data,
                        (int64_t *)block->header_contents+offset)
                            != TNG_SUCCESS)
                        {
                            printf("Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }
                    offset += size;
                }
            }
            break;
        case TNG_CHAR_DATA:
            for(i = 0; i < n_frames / stride_length; i++)
            {
                for(j = 0; j < data->n_values_per_frame; j++)
                {
                    len = strlen(data->values[i][j].c) + 1;
                    strncpy(block->block_contents+offset, data->values[i][j].c,
                            len);
                    offset += len;
                }
            }
            break;
        case TNG_DOUBLE_DATA:
        default:
            /* For speed reasons the compression multiplier is not used if the data
            * is not compressed.*/
            if(data->codec_id == TNG_UNCOMPRESSED)
            {
                for(i = 0; i < n_frames / stride_length; i++)
                {
                    for(j = 0; j < data->n_values_per_frame; j++)
                    {
                        memcpy(block->block_contents+offset, &data->values[i][j].d,
                            size);
                        if(tng_data->output_endianness_swap_func_64)
                        {
                            if(tng_data->output_endianness_swap_func_64(tng_data,
                            (int64_t *)block->header_contents+offset)
                                != TNG_SUCCESS)
                            {
                                printf("Cannot swap byte order. %s: %d\n",
                                        __FILE__, __LINE__);
                            }
                        }
                        offset += size;
                    }
                }
            }
            else
            {
                multiplier = data->compression_multiplier;
                for(i = 0; i < n_frames / stride_length; i++)
                {
                    for(j = 0; j < data->n_values_per_frame; j++)
                    {
                        d_temp = data->values[i][j].d * multiplier;
                        memcpy(block->block_contents+offset, &d_temp,
                            size);
                        if(tng_data->output_endianness_swap_func_64)
                        {
                            if(tng_data->output_endianness_swap_func_64(tng_data,
                            (int64_t *)block->header_contents+offset)
                                != TNG_SUCCESS)
                            {
                                printf("Cannot swap byte order. %s: %d\n",
                                        __FILE__, __LINE__);
                            }
                        }
                        offset += size;
                    }
                }
            }
        }
    }
    else
    {
        memset(block->block_contents+offset, 0, block->block_contents_size - offset);
    }

    switch(data->codec_id)
    {
    case TNG_XTC_COMPRESSION:
        printf("XTC compression not implemented yet.\n");
        break;
    case TNG_TNG_COMPRESSION:
        printf("TNG compression not implemented yet.\n");
        break;
#ifdef USE_ZLIB
    case TNG_GZIP_COMPRESSION:
//         printf("Before compression: %"PRId64"\n", block->block_contents_size);
        if(tng_gzip_compress(tng_data, block,
                             block->block_contents + data_start_pos,
                             block->block_contents_size - data_start_pos) !=
           TNG_SUCCESS)
        {
            printf("Could not write gzipped block data. %s: %d\n", __FILE__,
                   __LINE__);
            return(TNG_CRITICAL);
        }
//         printf("After compression: %"PRId64"\n", block->block_contents_size);
        break;
#endif
    }

    if(tng_block_header_write(tng_data, block, hash_mode) != TNG_SUCCESS)
    {
        printf("Cannot write header of file %s. %s: %d\n",
               tng_data->output_file_path, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
              tng_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    return(TNG_SUCCESS);
}

/** Read the contents of a data block (particle or non-particle data).
 * @param tng_data is a trajectory data container.
 * @param block is the block to store the data (should already contain
 * the block headers).
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_data_block_contents_read
                (tng_trajectory_t tng_data,
                 tng_gen_block_t block,
                 const tng_hash_mode hash_mode)
{
    int64_t n_values, codec_id, n_frames, first_frame_with_data;
    int64_t stride_length, block_n_particles, num_first_particle;
    double multiplier;
    char datatype, dependency, sparse_data;
    int offset = 0;
    tng_bool same_hash;

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }

    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to
     * disk even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
             tng_data->input_file) == 0)
    {
        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */

    if(hash_mode == TNG_USE_HASH)
    {
        hash_match_verify(block, &same_hash);
        if(same_hash != TNG_TRUE)
        {
            printf("'%s' data block contents corrupt. Hashes do not match. %s: %d\n",
                block->name, __FILE__, __LINE__);
    //         return(TNG_FAILURE);
        }
    }

    memcpy(&datatype, block->block_contents+offset,
           sizeof(datatype));
    offset += sizeof(datatype);

    memcpy(&dependency, block->block_contents+offset,
           sizeof(dependency));
    offset += sizeof(dependency);


    if(dependency & TNG_FRAME_DEPENDENT)
    {
        memcpy(&sparse_data, block->block_contents+offset,
               sizeof(sparse_data));
        offset += sizeof(sparse_data);
    }

    memcpy(&n_values, block->block_contents+offset,
        sizeof(n_values));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &n_values)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(n_values);

    memcpy(&codec_id, block->block_contents+offset,
        sizeof(codec_id));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &codec_id)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(codec_id);

    if(codec_id != TNG_UNCOMPRESSED)
    {
        memcpy(&multiplier, block->block_contents+offset,
            sizeof(multiplier));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                       (int64_t *) &multiplier)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(multiplier);
    }
    else
    {
        multiplier = 1;
    }

    if(dependency & TNG_FRAME_DEPENDENT)
    {
        if(sparse_data)
        {
            memcpy(&first_frame_with_data, block->block_contents+offset,
                sizeof(first_frame_with_data));
            if(tng_data->input_endianness_swap_func_64)
            {
                if(tng_data->input_endianness_swap_func_64(tng_data,
                                                        &first_frame_with_data)
                    != TNG_SUCCESS)
                {
                    printf("Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }
            offset += sizeof(first_frame_with_data);

            memcpy(&stride_length, block->block_contents+offset,
                sizeof(stride_length));
            if(tng_data->input_endianness_swap_func_64)
            {
                if(tng_data->input_endianness_swap_func_64(tng_data,
                                                           &stride_length)
                    != TNG_SUCCESS)
                {
                    printf("Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }
            offset += sizeof(stride_length);
        }
        else
        {
            first_frame_with_data = 0;
            stride_length = 1;
        }
        n_frames = tng_data->current_trajectory_frame_set.n_frames;
    }
    else
    {
        first_frame_with_data = 0;
        stride_length = 1;
        n_frames = 0;
    }

    if (dependency & TNG_PARTICLE_DEPENDENT)
    {
        memcpy(&num_first_particle, block->block_contents+offset,
               sizeof(num_first_particle));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                       &num_first_particle)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(num_first_particle);

        memcpy(&block_n_particles, block->block_contents+offset,
            sizeof(block_n_particles));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                       &block_n_particles)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(block_n_particles);
    }

    if (dependency & TNG_PARTICLE_DEPENDENT)
    {
        return(tng_particle_data_read(tng_data, block,
                                      &offset, datatype,
                                      num_first_particle,
                                      block_n_particles,
                                      first_frame_with_data,
                                      stride_length,
                                      n_frames, n_values,
                                      codec_id, multiplier));
    }
    else
    {
        return(tng_data_read(tng_data, block,
                             &offset, datatype,
                             first_frame_with_data,
                             stride_length,
                             n_frames, n_values,
                             codec_id, multiplier));
    }
}

/** Update the md5 hash of a block already written to the file
 * @param tng_data is a trajectory data container.
 * @param block is the block, of which to update the md5 hash.
 * @param header_start_pos is the file position where the block header starts.
 * @param contents_start_pos is the file position where the block contents
 * start.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_md5_hash_update(tng_trajectory_t tng_data,
                                               tng_gen_block_t block,
                                               const int64_t header_start_pos,
                                               const int64_t contents_start_pos)
{
    if(block->block_contents)
    {
        free(block->block_contents);
    }

    block->block_contents = malloc(block->block_contents_size);
    fseek(tng_data->output_file, contents_start_pos, SEEK_SET);
    if(fread(block->block_contents, block->block_contents_size, 1,
            tng_data->output_file) == 0)
    {
        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    tng_block_hash_generate(block);

    fseek(tng_data->output_file, header_start_pos + 3 * sizeof(int64_t),
          SEEK_SET);
    fwrite(block->hash, TNG_HASH_LEN, 1, tng_data->output_file);

    return(TNG_SUCCESS);
}

/** Update the frame set pointers in the file header (general info block),
 * already written to disk
 * @param tng_data is a trajectory data container.
 * @param hash_mode specifies whether to update the block md5 hash when
 * updating the pointers.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_header_pointers_update
                (tng_trajectory_t tng_data, const tng_hash_mode hash_mode)
{
    tng_gen_block_t block;
    FILE *temp = tng_data->input_file;
    int64_t output_file_pos, pos, contents_start_pos;

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        printf("Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    tng_data->input_file = tng_data->output_file;

    tng_block_init(&block);

    output_file_pos = ftell(tng_data->output_file);
    fseek(tng_data->output_file, 0, SEEK_SET);

    if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
    {
        printf("Cannot read general info header. %s: %d\n",
               __FILE__, __LINE__);
        tng_data->input_file = temp;
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    contents_start_pos = ftell(tng_data->output_file);

    fseek(tng_data->output_file, block->block_contents_size - 4 *
          sizeof(int64_t), SEEK_CUR);

    tng_data->input_file = temp;

    pos = tng_data->first_trajectory_frame_set_output_file_pos;

    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                    &pos)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(fwrite(&pos, sizeof(int64_t), 1, tng_data->output_file) != 1)
    {
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    pos = tng_data->last_trajectory_frame_set_output_file_pos;

    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                    &pos)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(fwrite(&pos,
        sizeof(int64_t), 1, tng_data->output_file) != 1)
    {
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    if(hash_mode == TNG_USE_HASH)
    {
        tng_md5_hash_update(tng_data, block, 0, contents_start_pos);
    }

    tng_block_destroy(&block);

    fseek(tng_data->output_file, output_file_pos, SEEK_SET);

    return(TNG_SUCCESS);
}

/** Update the frame set pointers in the current frame set block, already
 * written to disk
 * @param tng_data is a trajectory data container.
 * @param hash_mode specifies whether to update the block md5 hash when
 * updating the pointers.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_frame_set_pointers_update
                (tng_trajectory_t tng_data, const tng_hash_mode hash_mode)
{
    tng_gen_block_t block;
    tng_trajectory_frame_set_t frame_set;
    FILE *temp = tng_data->input_file;
    int64_t pos, output_file_pos, header_start_pos, contents_start_pos;

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        printf("Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    tng_block_init(&block);
    output_file_pos = ftell(tng_data->output_file);

    tng_data->input_file = tng_data->output_file;

    frame_set = &tng_data->current_trajectory_frame_set;

    /* Update previous frame set */
    if(frame_set->prev_frame_set_file_pos != -1 &&
       frame_set->prev_frame_set_file_pos != 0)
    {
        fseek(tng_data->output_file, frame_set->prev_frame_set_file_pos,
              SEEK_SET);

        header_start_pos = frame_set->prev_frame_set_file_pos;

        if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
        {
            printf("Cannot read frame header. %s: %d\n",
                __FILE__, __LINE__);
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        contents_start_pos = ftell(tng_data->output_file);

        fseek(tng_data->output_file, block->block_contents_size - 6 *
            sizeof(int64_t), SEEK_CUR);

        pos = tng_data->current_trajectory_frame_set_output_file_pos;

        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                        &pos)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }

        if(fwrite(&pos, sizeof(int64_t), 1, tng_data->output_file) != 1)
        {
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        if(hash_mode == TNG_USE_HASH)
        {
            tng_md5_hash_update(tng_data, block, header_start_pos,
                                contents_start_pos);
        }

        fseek(tng_data->output_file, output_file_pos, SEEK_SET);
    }

    /* Update the frame set one medium stride step before */
    if(frame_set->medium_stride_prev_frame_set_file_pos != -1 &&
       frame_set->medium_stride_prev_frame_set_file_pos != 0)
    {
        fseek(tng_data->output_file,
              frame_set->medium_stride_prev_frame_set_file_pos,
              SEEK_SET);

        if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
        {
            printf("Cannot read frame set header. %s: %d\n",
                __FILE__, __LINE__);
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        contents_start_pos = ftell(tng_data->output_file);

        fseek(tng_data->output_file, block->block_contents_size - 4 *
            sizeof(int64_t), SEEK_CUR);

        pos = tng_data->current_trajectory_frame_set_output_file_pos;

        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                        &pos)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }

        if(fwrite(&pos, sizeof(int64_t), 1, tng_data->output_file) != 1)
        {
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        if(hash_mode == TNG_USE_HASH)
        {
            tng_md5_hash_update(tng_data, block,
                                frame_set->medium_stride_prev_frame_set_file_pos,
                                contents_start_pos);
        }
    }

    /* Update the frame set one long stride step before */
    if(frame_set->long_stride_prev_frame_set_file_pos != -1 &&
       frame_set->long_stride_prev_frame_set_file_pos != 0)
    {
        fseek(tng_data->output_file,
              frame_set->long_stride_prev_frame_set_file_pos,
              SEEK_SET);

        if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
        {
            printf("Cannot read frame set header. %s: %d\n",
                __FILE__, __LINE__);
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        contents_start_pos = ftell(tng_data->output_file);

        fseek(tng_data->output_file, block->block_contents_size - 2 *
            sizeof(int64_t), SEEK_CUR);

        pos = tng_data->current_trajectory_frame_set_output_file_pos;

        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                        &pos)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }

        if(fwrite(&pos, sizeof(int64_t), 1, tng_data->output_file) != 1)
        {
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        if(hash_mode == TNG_USE_HASH)
        {
            tng_md5_hash_update(tng_data, block,
                                frame_set->long_stride_prev_frame_set_file_pos,
                                contents_start_pos);
        }
    }

    fseek(tng_data->output_file, output_file_pos, SEEK_SET);

    tng_data->input_file = temp;

    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

/** Move the blocks in a frame set so that there is no unused space between
 * them. This can only be done on the last frame set in the file and should
 * be done e.g. if the last frame set in the file has fewer frames than
 * default or after compressing data blocks in a frame set.
 * @param tng_data is a trajectory data container.
 * @details the current_trajectory_frame_set is the one that will be modified.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the frame set
 * cannot be aligned or TNG_CRITICAL (2) if a major error has occured.
 * FIXME: This function is not finished!!!
 */
static tng_function_status tng_frame_set_align(tng_trajectory_t tng_data)
{
    tng_gen_block_t block;
    tng_trajectory_frame_set_t frame_set;
    FILE *temp = tng_data->input_file;
    int64_t pos, contents_start_pos, output_file_len;

    frame_set = &tng_data->current_trajectory_frame_set;

    if(frame_set->n_written_frames == frame_set->n_frames)
    {
        return(TNG_SUCCESS);
    }

    if(tng_data->current_trajectory_frame_set_output_file_pos !=
       tng_data->last_trajectory_frame_set_output_file_pos)
    {
    }

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        printf("Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    tng_block_init(&block);
//     output_file_pos = ftell(tng_data->output_file);

    tng_data->input_file = tng_data->output_file;

    pos = tng_data->current_trajectory_frame_set_output_file_pos;

    fseek(tng_data->output_file, pos, SEEK_SET);
    if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
    {
        printf("Cannot read frame set header. %s: %d\n",
            __FILE__, __LINE__);
        tng_data->input_file = temp;
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    contents_start_pos = ftell(tng_data->output_file);

    fseek(tng_data->output_file, 0, SEEK_END);
    output_file_len = ftell(tng_data->output_file);
    pos = contents_start_pos + block->block_contents_size;
    fseek(tng_data->output_file, pos,
          SEEK_SET);

    while(pos < output_file_len)
    {
        if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
        {
            printf("Cannot read block header at pos %"PRId64". %s: %d\n", pos,
                   __FILE__, __LINE__);
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
        pos += block->header_contents_size + block->block_contents_size;
        fseek(tng_data->output_file, pos, SEEK_SET);
    }

    return(TNG_SUCCESS);
}

/** Finish writing the current frame set. Update the number of frames
 * and the hashes of the frame set and all its data blocks (if hash_mode
 * == TNG_USE_HASH).
 * @param tng_data is a trajectory data container.
 * @param hash_mode specifies whether to update the block md5 hash when
 * updating the pointers.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_frame_set_finalize
                (tng_trajectory_t tng_data, const tng_hash_mode hash_mode)
{
    tng_gen_block_t block;
    tng_trajectory_frame_set_t frame_set;
    FILE *temp = tng_data->input_file;
    int64_t pos, contents_start_pos, output_file_len;

    frame_set = &tng_data->current_trajectory_frame_set;

    if(frame_set->n_written_frames == frame_set->n_frames)
    {
        return(TNG_SUCCESS);
    }

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        printf("Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    tng_block_init(&block);
//     output_file_pos = ftell(tng_data->output_file);

    tng_data->input_file = tng_data->output_file;

    pos = tng_data->current_trajectory_frame_set_output_file_pos;

    fseek(tng_data->output_file, pos, SEEK_SET);

    if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
    {
        printf("Cannot read frame set header. %s: %d\n",
            __FILE__, __LINE__);
        tng_data->input_file = temp;
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    contents_start_pos = ftell(tng_data->output_file);

    fseek(tng_data->output_file, sizeof(frame_set->first_frame), SEEK_CUR);
    if(fwrite(&frame_set->n_written_frames, sizeof(frame_set->n_frames),
              1, tng_data->output_file) != 1)
    {
        tng_data->input_file = temp;
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }


    if(hash_mode == TNG_USE_HASH)
    {
        tng_md5_hash_update(tng_data, block, pos,
                            pos + block->header_contents_size);
    }

    fseek(tng_data->output_file, 0, SEEK_END);
    output_file_len = ftell(tng_data->output_file);
    pos = contents_start_pos + block->block_contents_size;
    fseek(tng_data->output_file, pos, SEEK_SET);

    while(pos < output_file_len)
    {
        if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
        {
            printf("Cannot read block header at pos %"PRId64". %s: %d\n", pos,
                   __FILE__, __LINE__);
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        if(hash_mode == TNG_USE_HASH)
        {
            tng_md5_hash_update(tng_data, block, pos,
                                pos + block->header_contents_size);
        }
        pos += block->header_contents_size + block->block_contents_size;
        fseek(tng_data->output_file, pos, SEEK_SET);
    }

    tng_data->input_file = temp;
    tng_block_destroy(&block);
    return(TNG_SUCCESS);
}


/** Sets the name of a file contents block
 * @param tng_data is a trajectory data container.
 * @param block is the block, of which to change names.
 * @param new_name is the new name of the block.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_block_name_set(tng_trajectory_t tng_data,
                                              tng_gen_block_t block,
                                              const char *new_name)
{
    int len;

    len = tng_min(strlen(new_name) + 1, TNG_MAX_STR_LEN);

    /* If the currently stored string length is not enough to store the new
     * string it is freed and reallocated. */
    if(block->name && strlen(block->name) < len)
    {
        free(block->name);
        block->name = 0;
    }
    if(!block->name)
    {
        block->name = malloc(len);
        if(!block->name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(block->name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_atom_name_set(tng_trajectory_t tng_data,
                                      tng_atom_t atom,
                                      const char *new_name)
{
    int len;

    len = tng_min(strlen(new_name) + 1, TNG_MAX_STR_LEN);

    /* If the currently stored string length is not enough to store the new
     * string it is freed and reallocated. */
    if(atom->name && strlen(atom->name) < len)
    {
        free(atom->name);
        atom->name = 0;
    }
    if(!atom->name)
    {
        atom->name = malloc(len);
        if(!atom->name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(atom->name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_atom_type_set(tng_trajectory_t tng_data,
                                      tng_atom_t atom,
                                      const char *new_type)
{
    int len;

    len = tng_min(strlen(new_type) + 1, TNG_MAX_STR_LEN);

    /* If the currently stored string length is not enough to store the new
     * string it is freed and reallocated. */
    if(atom->atom_type && strlen(atom->atom_type) < len)
    {
        free(atom->atom_type);
        atom->atom_type = 0;
    }
    if(!atom->atom_type)
    {
        atom->atom_type = malloc(len);
        if(!atom->atom_type)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(atom->atom_type, new_type, len);

    return(TNG_SUCCESS);
}

/** Initialise an atom struct
 * @param atom is the atom to initialise.
 * @return TNG_SUCCESS (0) if successful.
 */
static tng_function_status tng_atom_init(tng_atom_t atom)
{
    atom->name = 0;
    atom->atom_type = 0;

    return(TNG_SUCCESS);
}

/** Free the memory in an atom struct
 * @param atom is the atom to destroy.
 * @return TNG_SUCCESS (0) if successful.
 */
static tng_function_status tng_atom_destroy(tng_atom_t atom)
{
    if(atom->name)
    {
        free(atom->name);
        atom->name = 0;
    }
    if(atom->atom_type)
    {
        free(atom->atom_type);
        atom->atom_type = 0;
    }

    return(TNG_SUCCESS);
}

tng_function_status tng_molecule_add(tng_trajectory_t tng_data,
                                     const char *name,
                                     tng_molecule_t *molecule)
{
    tng_molecule_t new_molecules;
    int64_t *new_molecule_cnt_list;
    int id, i;
    tng_bool found_id = TNG_TRUE;

    new_molecules = realloc(tng_data->molecules,
                            sizeof(struct tng_molecule) *
                            (tng_data->n_molecules + 1));

    if(!new_molecules)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct tng_molecule) * (tng_data->n_molecules + 1),
               __FILE__, __LINE__);
        free(tng_data->molecules);
        return(TNG_CRITICAL);
    }

    new_molecule_cnt_list = realloc(tng_data->molecule_cnt_list,
                                    sizeof(int64_t) *
                                    (tng_data->n_molecules + 1));

    if(!new_molecule_cnt_list)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(int64_t) * (tng_data->n_molecules + 1),
               __FILE__, __LINE__);
        free(tng_data->molecule_cnt_list);
        return(TNG_CRITICAL);
    }

    tng_data->molecules = new_molecules;
    tng_data->molecule_cnt_list = new_molecule_cnt_list;

    *molecule = &new_molecules[tng_data->n_molecules];

    tng_molecule_init(tng_data, *molecule);
    tng_molecule_name_set(tng_data, *molecule, name);

    /* FIXME: Should this be a function argument instead? */
    tng_data->molecule_cnt_list[tng_data->n_molecules] = 0;

    /* Find an unused ID */
    id = 0;
    while(found_id)
    {
        found_id = TNG_FALSE;
        for(i = tng_data->n_molecules; i--;)
        {
            if(tng_data->molecules[i].id == id)
            {
                found_id = TNG_TRUE;
                i = 0;
            }
        }
        if(found_id)
        {
            id++;
        }
    }

    (*molecule)->id = id;

    tng_data->n_molecules++;

    return(TNG_SUCCESS);
}

tng_function_status tng_molecule_name_set(tng_trajectory_t tng_data,
                                          tng_molecule_t molecule,
                                          const char *new_name)
{
    int len;

    len = tng_min(strlen(new_name) + 1, TNG_MAX_STR_LEN);

    /* If the currently stored string length is not enough to store the new
     * string it is freed and reallocated. */
    if(molecule->name && strlen(molecule->name) < len)
    {
        free(molecule->name);
        molecule->name = 0;
    }
    if(!molecule->name)
    {
        molecule->name = malloc(len);
        if(!molecule->name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(molecule->name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_molecule_cnt_get(tng_trajectory_t tng_data,
                                         tng_molecule_t molecule,
                                         int64_t *cnt)
{
    int i, index = -1;

    for(i = tng_data->n_molecules; i--;)
    {
        if(&tng_data->molecules[i] == molecule)
        {
            index = i;
            i = 0;
        }
    }
    if(index == -1)
    {
        return(TNG_FAILURE);
    }
    *cnt = tng_data->molecule_cnt_list[index];

    return(TNG_SUCCESS);
}

tng_function_status tng_molecule_cnt_set(tng_trajectory_t tng_data,
                                         tng_molecule_t molecule,
                                         const int64_t cnt)
{
    int i, index = -1, old_cnt;

    for(i = tng_data->n_molecules; i--;)
    {
        if(&tng_data->molecules[i] == molecule)
        {
            index = i;
            i = 0;
        }
    }
    if(index == -1)
    {
        return(TNG_FAILURE);
    }
    old_cnt = tng_data->molecule_cnt_list[index];
    tng_data->molecule_cnt_list[index] = cnt;

    tng_data->n_particles += (cnt-old_cnt) *
                             tng_data->molecules[index].n_atoms;

    return(TNG_SUCCESS);
}

tng_function_status tng_molecule_chain_find(tng_trajectory_t tng_data,
                                            tng_molecule_t molecule,
                                            const char *name,
                                            int64_t nr,
                                            tng_chain_t *chain)
{
    int i, n_chains;

    n_chains = molecule->n_chains;

    for(i = 0; i < n_chains; i++)
    {
        *chain = &molecule->chains[i];
        if(name[0] != 0 || strcmp(name, (*chain)->name) == 0)
        {
            if(nr == -1 || nr == (*chain)->id)
            {
                return(TNG_SUCCESS);
            }
        }
    }

    *chain = 0;

    return(TNG_FAILURE);
}

/* FIXME: For v. 2 the chain nr should also be possible to specify. */
tng_function_status tng_molecule_chain_add(tng_trajectory_t tng_data,
                                           tng_molecule_t molecule,
                                           const char *name,
                                           tng_chain_t *chain)
{
    tng_chain_t new_chains;

    new_chains = realloc(molecule->chains,
                         sizeof(struct tng_chain) *
                         (molecule->n_chains + 1));

    if(!new_chains)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct tng_chain) * (molecule->n_chains + 1),
               __FILE__, __LINE__);
        free(molecule->chains);
        return(TNG_CRITICAL);
    }

    molecule->chains = new_chains;

    *chain = &new_chains[molecule->n_chains];
    (*chain)->name = 0;

    tng_chain_name_set(tng_data, *chain, name);

    (*chain)->molecule = molecule;
    (*chain)->n_residues = 0;

    molecule->n_chains++;

    (*chain)->id = molecule->n_chains;

    return(TNG_SUCCESS);
}

tng_function_status tng_chain_name_set(tng_trajectory_t tng_data,
                                       tng_chain_t chain,
                                       const char *new_name)
{
    int len;

    len = tng_min(strlen(new_name) + 1, TNG_MAX_STR_LEN);

    /* If the currently stored string length is not enough to store the new
     * string it is freed and reallocated. */
    if(chain->name && strlen(chain->name) < len)
    {
        free(chain->name);
        chain->name = 0;
    }
    if(!chain->name)
    {
        chain->name = malloc(len);
        if(!chain->name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(chain->name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_chain_residue_find(tng_trajectory_t tng_data,
                                           tng_chain_t chain,
                                           const char *name,
                                           int64_t nr,
                                           tng_residue_t *residue)
{
    int i, n_residues;

    n_residues = chain->n_residues;

    for(i = 0; i < n_residues; i++)
    {
        *residue = &chain->residues[i];
        if(name[0] != 0 || strcmp(name, (*residue)->name) == 0)
        {
            if(nr == -1 || nr == (*residue)->id)
            {
                return(TNG_SUCCESS);
            }
        }
    }

    *residue = 0;

    return(TNG_FAILURE);
}

/* FIXME: For v. 2 the residue nr should also be possible to specify. */
tng_function_status tng_chain_residue_add(tng_trajectory_t tng_data,
                                             tng_chain_t chain,
                                             const char *name,
                                             tng_residue_t *residue)
{
    int curr_index;
    tng_residue_t new_residues, temp_residue, last_residue;
    tng_molecule_t molecule = chain->molecule;

    if(chain->n_residues)
    {
        curr_index = chain->residues - molecule->residues;
    }
    else
    {
        curr_index = -1;
    }

    new_residues = realloc(molecule->residues,
                           sizeof(struct tng_residue) *
                           (molecule->n_residues + 1));

    if(!new_residues)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct tng_residue) * (molecule->n_residues + 1),
               __FILE__, __LINE__);
        free(molecule->residues);
        return(TNG_CRITICAL);
    }

    molecule->residues = new_residues;

    if(curr_index != -1)
    {
        chain->residues = new_residues + curr_index;
        if(molecule->n_residues)
        {
            last_residue = &new_residues[molecule->n_residues - 1];

            temp_residue = chain->residues + (chain->n_residues - 1);
            /* Make space in list of residues to add the new residues together with the other
            * residues of this chain */
            if(temp_residue != last_residue)
            {
                temp_residue++;
                memmove(temp_residue + 1, temp_residue,
                        last_residue - temp_residue);
            }
        }
    }
    else
    {
        curr_index = molecule->n_residues;
    }

    *residue = &new_residues[curr_index + chain->n_residues];

    if(!chain->n_residues)
    {
        chain->residues = *residue;
    }

    (*residue)->name = 0;
    tng_residue_name_set(tng_data, *residue, name);

    (*residue)->chain = chain;
    (*residue)->n_atoms = 0;
    (*residue)->atoms_offset = 0;

    chain->n_residues++;
    molecule->n_residues++;

    (*residue)->id = chain->n_residues;

    return(TNG_SUCCESS);
}

tng_function_status tng_residue_name_set(tng_trajectory_t tng_data,
                                         tng_residue_t residue,
                                         const char *new_name)
{
    int len;

    len = tng_min(strlen(new_name) + 1, TNG_MAX_STR_LEN);

    /* If the currently stored string length is not enough to store the new
     * string it is freed and reallocated. */
    if(residue->name && strlen(residue->name) < len)
    {
        free(residue->name);
        residue->name = 0;
    }
    if(!residue->name)
    {
        residue->name = malloc(len);
        if(!residue->name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(residue->name, new_name, len);

    return(TNG_SUCCESS);
}

/* FIXME: For v. 2 the atom nr should also be possible to specify. */
tng_function_status tng_residue_atom_add(tng_trajectory_t tng_data,
                                         tng_residue_t residue,
                                         const char *atom_name,
                                         const char *atom_type,
                                         tng_atom_t *atom)
{
    tng_atom_t new_atoms;
    tng_molecule_t molecule = residue->chain->molecule;

    if(!residue->n_atoms)
    {
        residue->atoms_offset = molecule->n_atoms;
    }

    new_atoms = realloc(molecule->atoms,
                        sizeof(struct tng_atom) *
                        (molecule->n_atoms + 1));

    if(!new_atoms)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct tng_atom) * (molecule->n_atoms + 1),
               __FILE__, __LINE__);
        free(molecule->atoms);
        return(TNG_CRITICAL);
    }

    molecule->atoms = new_atoms;

    *atom = &new_atoms[molecule->n_atoms];

    tng_atom_init(*atom);
    tng_atom_name_set(tng_data, *atom, atom_name);
    tng_atom_type_set(tng_data, *atom, atom_type);

    (*atom)->residue = residue;

    residue->n_atoms++;
    molecule->n_atoms++;

    (*atom)->id = molecule->n_atoms;

    return(TNG_SUCCESS);
}


tng_function_status tng_molecule_init(const tng_trajectory_t tng_data,
                                      tng_molecule_t molecule)
{
    molecule->quaternary_str = 1;
    molecule->name = 0;
    molecule->n_chains = 0;
    molecule->chains = 0;
    molecule->n_residues = 0;
    molecule->residues = 0;
    molecule->n_atoms = 0;
    molecule->atoms = 0;
    molecule->n_bonds = 0;
    molecule->bonds = 0;

    return(TNG_SUCCESS);
}

tng_function_status tng_molecule_destroy(const tng_trajectory_t tng_data,
                                         tng_molecule_t molecule)
{
    int i;

    if(molecule->name)
    {
        free(molecule->name);
        molecule->name = 0;
    }

    if(molecule->chains)
    {
        for(i = molecule->n_chains; i--;)
        {
            if(molecule->chains[i].name)
            {
                free(molecule->chains[i].name);
                molecule->chains[i].name = 0;
            }
        }
        free(molecule->chains);
        molecule->chains = 0;
    }
    molecule->n_chains = 0;

    if(molecule->residues)
    {
        for(i = molecule->n_residues; i--;)
        {
            if(molecule->residues[i].name)
            {
                free(molecule->residues[i].name);
                molecule->residues[i].name = 0;
            }
        }
        free(molecule->residues);
        molecule->residues = 0;
    }
    molecule->n_residues = 0;

    if(molecule->atoms)
    {
        for(i = molecule->n_atoms; i--;)
        {
            tng_atom_destroy(&molecule->atoms[i]);
        }
        free(molecule->atoms);
        molecule->atoms = 0;
    }
    molecule->n_atoms = 0;

    if(molecule->bonds)
    {
        free(molecule->bonds);
        molecule->bonds = 0;
    }
    molecule->n_bonds = 0;

    return(TNG_SUCCESS);
}

tng_function_status tng_molecule_name_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    int64_t cnt = 0, i, *molecule_cnt_list;
    tng_molecule_t mol;
    tng_bool found = TNG_FALSE;

    if(tng_data->var_num_atoms_flag)
    {
        molecule_cnt_list = tng_data->current_trajectory_frame_set.
                            molecule_cnt_list;
    }
    else
    {
        molecule_cnt_list = tng_data->molecule_cnt_list;
    }

    for(i = 0; i < tng_data->n_molecules; i++)
    {
        mol = &tng_data->molecules[i];
        if(cnt + mol->n_atoms * molecule_cnt_list[i] - 1 < nr)
        {
            cnt += mol->n_atoms * molecule_cnt_list[i];
            continue;
        }
        found = TNG_TRUE;
        break;
    }
    if(!found)
    {
        return(TNG_FAILURE);
    }

    strncpy(name, mol->name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(mol->name) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_chain_name_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    int64_t cnt = 0, i, *molecule_cnt_list;
    tng_molecule_t mol;
    tng_atom_t atom;
    tng_bool found = TNG_FALSE;

    if(tng_data->var_num_atoms_flag)
    {
        molecule_cnt_list = tng_data->current_trajectory_frame_set.
                            molecule_cnt_list;
    }
    else
    {
        molecule_cnt_list = tng_data->molecule_cnt_list;
    }

    for(i = 0; i < tng_data->n_molecules; i++)
    {
        mol = &tng_data->molecules[i];
        if(cnt + mol->n_atoms * molecule_cnt_list[i] - 1 < nr)
        {
            cnt += mol->n_atoms * molecule_cnt_list[i];
            continue;
        }
        atom = &mol->atoms[nr % mol->n_atoms];
        found = TNG_TRUE;
        break;
    }
    if(!found)
    {
        return(TNG_FAILURE);
    }
    if(!atom->residue || !atom->residue->chain)
    {
        return(TNG_FAILURE);
    }

    strncpy(name, atom->residue->chain->name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(atom->residue->chain->name) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_residue_name_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    int64_t cnt = 0, i, *molecule_cnt_list;
    tng_molecule_t mol;
    tng_atom_t atom;
    tng_bool found = TNG_FALSE;

    if(tng_data->var_num_atoms_flag)
    {
        molecule_cnt_list = tng_data->current_trajectory_frame_set.
                            molecule_cnt_list;
    }
    else
    {
        molecule_cnt_list = tng_data->molecule_cnt_list;
    }

    for(i = 0; i < tng_data->n_molecules; i++)
    {
        mol = &tng_data->molecules[i];
        if(cnt + mol->n_atoms * molecule_cnt_list[i] - 1 < nr)
        {
            cnt += mol->n_atoms * molecule_cnt_list[i];
            continue;
        }
        atom = &mol->atoms[nr % mol->n_atoms];
        found = TNG_TRUE;
        break;
    }
    if(!found)
    {
        return(TNG_FAILURE);
    }
    if(!atom->residue)
    {
        return(TNG_FAILURE);
    }

    strncpy(name, atom->residue->name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(atom->residue->name) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_atom_name_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    int64_t cnt = 0, i, *molecule_cnt_list;
    tng_molecule_t mol;
    tng_atom_t atom;
    tng_bool found = TNG_FALSE;

    if(tng_data->var_num_atoms_flag)
    {
        molecule_cnt_list = tng_data->current_trajectory_frame_set.
                            molecule_cnt_list;
    }
    else
    {
        molecule_cnt_list = tng_data->molecule_cnt_list;
    }

    for(i = 0; i < tng_data->n_molecules; i++)
    {
        mol = &tng_data->molecules[i];
        if(cnt + mol->n_atoms * molecule_cnt_list[i] - 1 < nr)
        {
            cnt += mol->n_atoms * molecule_cnt_list[i];
            continue;
        }
        atom = &mol->atoms[nr % mol->n_atoms];
        found = TNG_TRUE;
        break;
    }
    if(!found)
    {
        return(TNG_FAILURE);
    }

    strncpy(name, atom->name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(atom->name) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_atom_type_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *type,
                 int max_len)
{
    int64_t cnt = 0, i, *molecule_cnt_list;
    tng_molecule_t mol;
    tng_atom_t atom;
    tng_bool found = TNG_FALSE;

    if(tng_data->var_num_atoms_flag)
    {
        molecule_cnt_list = tng_data->current_trajectory_frame_set.
                            molecule_cnt_list;
    }
    else
    {
        molecule_cnt_list = tng_data->molecule_cnt_list;
    }

    for(i = 0; i < tng_data->n_molecules; i++)
    {
        mol = &tng_data->molecules[i];
        if(cnt + mol->n_atoms * molecule_cnt_list[i] - 1 < nr)
        {
            cnt += mol->n_atoms * molecule_cnt_list[i];
            continue;
        }
        atom = &mol->atoms[nr % mol->n_atoms];
        found = TNG_TRUE;
        break;
    }
    if(!found)
    {
        return(TNG_FAILURE);
    }

    strncpy(type, atom->atom_type, max_len - 1);
    type[max_len - 1] = 0;

    if(strlen(atom->atom_type) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_particle_mapping_add
                (tng_trajectory_t tng_data,
                 const int64_t num_first_particle,
                 const int64_t n_particles,
                 const int64_t *mapping_table)
{
    int64_t i;
    tng_particle_mapping_t mapping;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;

    /* Sanity check of the particle ranges. Split into multiple if
     * statements for improved readability */
    for(i = 0; i < frame_set->n_mapping_blocks; i++)
    {
        mapping = &frame_set->mappings[i];
        if(num_first_particle >= mapping->num_first_particle &&
           num_first_particle < mapping->num_first_particle +
                                   mapping->n_particles)
        {
            printf("Particle mapping overlap. %s: %d\n", __FILE__, __LINE__);
            return(TNG_FAILURE);
        }
        if(num_first_particle + n_particles >=
           mapping->num_first_particle &&
           num_first_particle + n_particles <
           mapping->num_first_particle + mapping->n_particles)
        {
            printf("Particle mapping overlap. %s: %d\n", __FILE__, __LINE__);
            return(TNG_FAILURE);
        }
        if(mapping->num_first_particle >= num_first_particle &&
           mapping->num_first_particle < num_first_particle +
                                            n_particles)
        {
            printf("Particle mapping overlap. %s: %d\n", __FILE__, __LINE__);
            return(TNG_FAILURE);
        }
        if(mapping->num_first_particle + mapping->n_particles >
           num_first_particle &&
           mapping->num_first_particle + mapping->n_particles <
           num_first_particle + n_particles)
        {
            printf("Particle mapping overlap. %s: %d\n", __FILE__, __LINE__);
            return(TNG_FAILURE);
        }
    }

    frame_set->n_mapping_blocks++;

    mapping = realloc(frame_set->mappings, sizeof(struct tng_particle_mapping) *
                      frame_set->n_mapping_blocks);

    if(!mapping)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct tng_particle_mapping)*frame_set->n_mapping_blocks,
               __FILE__, __LINE__);
        free(frame_set->mappings);
        return(TNG_CRITICAL);
    }
    frame_set->mappings = mapping;

    mapping = &frame_set->mappings[frame_set->n_mapping_blocks - 1];

    mapping->num_first_particle = num_first_particle;
    mapping->n_particles = n_particles;

    mapping->real_particle_numbers = malloc(sizeof(int64_t) * n_particles);
    if(!mapping->real_particle_numbers)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(int64_t) * n_particles, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    for(i=0; i<n_particles; i++)
    {
        mapping->real_particle_numbers[i] = mapping_table[i];
    }

    return(TNG_SUCCESS);
}

tng_function_status tng_trajectory_init(tng_trajectory_t *tng_data_p)
{
    time_t seconds;
    tng_trajectory_frame_set_t frame_set;
    tng_trajectory_t tng_data;

    *tng_data_p = malloc(sizeof(struct tng_trajectory));
    if(!tng_data_p)
    {
        printf("Cannot allocate memory (%lu bytes). %s: %d\n",
               sizeof(struct tng_trajectory), __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    tng_data = *tng_data_p;

    frame_set = &tng_data->current_trajectory_frame_set;

    tng_data->input_file_path = 0;
    tng_data->input_file = 0;
    tng_data->input_file_len = 0;
    tng_data->output_file_path = 0;
    tng_data->output_file = 0;

    tng_data->first_program_name = 0;
    tng_data->first_user_name = 0;
    tng_data->first_computer_name = 0;
    tng_data->first_pgp_signature = 0;
    tng_data->last_program_name = 0;
    tng_data->last_user_name = 0;
    tng_data->last_computer_name = 0;
    tng_data->last_pgp_signature = 0;
    tng_data->forcefield_name = 0;

    seconds = time(0);
    if ( seconds == -1)
    {
        printf("Cannot get time. %s: %d\n", __FILE__, __LINE__);
    }
    else
    {
        tng_data->time = seconds;
    }

    tng_data->var_num_atoms_flag = TNG_CONSTANT_N_ATOMS;
    tng_data->first_trajectory_frame_set_input_file_pos = -1;
    tng_data->last_trajectory_frame_set_input_file_pos = -1;
    tng_data->current_trajectory_frame_set_input_file_pos = -1;
    tng_data->first_trajectory_frame_set_output_file_pos = -1;
    tng_data->last_trajectory_frame_set_output_file_pos = -1;
    tng_data->current_trajectory_frame_set_output_file_pos = -1;
    tng_data->frame_set_n_frames = 100;
    tng_data->n_trajectory_frame_sets = 0;
    tng_data->n_trajectory_blocks = 0;
    tng_data->medium_stride_length = 100;
    tng_data->long_stride_length = 10000;

    tng_data->n_particle_data_blocks = 0;
    tng_data->n_data_blocks = 0;

    tng_data->non_tr_particle_data = 0;
    tng_data->non_tr_data = 0;

    frame_set->first_frame = -1;
    frame_set->n_mapping_blocks = 0;
    frame_set->mappings = 0;
    frame_set->molecule_cnt_list = 0;

    frame_set->n_particle_data_blocks = 0;
    frame_set->n_data_blocks = 0;

    frame_set->tr_particle_data = 0;
    frame_set->tr_data = 0;

    frame_set->n_written_frames = 0;

    frame_set->next_frame_set_file_pos = -1;
    frame_set->prev_frame_set_file_pos = -1;
    frame_set->medium_stride_next_frame_set_file_pos = -1;
    frame_set->medium_stride_prev_frame_set_file_pos = -1;
    frame_set->long_stride_next_frame_set_file_pos = -1;
    frame_set->long_stride_prev_frame_set_file_pos = -1;

    tng_data->n_molecules = 0;
    tng_data->molecules = 0;
    tng_data->molecule_cnt_list = 0;
    tng_data->n_particles = 0;

    /* Check the endianness of the computer */
    static int32_t endianness_32 = 0x01234567;
    /* 0x01234567 */
    if ( *(const uint8_t*)&endianness_32 == 0x01 )
    {
        tng_data->endianness_32 = TNG_BIG_ENDIAN_32;
    }

    /* 0x67452301 */
    else if( *(const uint8_t*)&endianness_32 == 0x67 )
    {
        tng_data->endianness_32 = TNG_LITTLE_ENDIAN_32;

    }

    /* 0x45670123 */
    else if ( *(const uint8_t*)&endianness_32 == 0x45 )
    {
        tng_data->endianness_32 = TNG_BYTE_PAIR_SWAP_32;
    }

    static int64_t endianness_64 = 0x0123456789ABCDEF;
    /* 0x0123456789ABCDEF */
    if ( *(const uint8_t*)&endianness_64 == 0x01 )
    {
        tng_data->endianness_64 = TNG_BIG_ENDIAN_64;
    }

    /* 0xEFCDAB8967452301 */
    else if ( *(const uint8_t*)&endianness_64 == 0xEF )
    {
        tng_data->endianness_64 = TNG_LITTLE_ENDIAN_64;
    }

    /* 0x89ABCDEF01234567 */
    else if ( *(const uint8_t*)&endianness_64 == 0x89 )
    {
        tng_data->endianness_64 = TNG_QUAD_SWAP_64;
    }

    /* 0x45670123CDEF89AB */
    else if ( *(const uint8_t*)&endianness_64 == 0x45 )
    {
        tng_data->endianness_64 = TNG_BYTE_PAIR_SWAP_64;
    }

    /* 0x23016745AB89EFCD */
    else if ( *(const uint8_t*)&endianness_64 == 0x23 )
    {
        tng_data->endianness_64 = TNG_BYTE_SWAP_64;
    }

    /* By default do not swap the byte order, i.e. keep the byte order of the
     * architecture. The input file endianness will be set when reading the
     * header. The output endianness can be changed - before the file is
     * written. */
    tng_data->input_endianness_swap_func_32 = 0;
    tng_data->input_endianness_swap_func_64 = 0;
    tng_data->output_endianness_swap_func_32 = 0;
    tng_data->output_endianness_swap_func_64 = 0;

    tng_data->current_trajectory_frame_set.next_frame_set_file_pos = -1;
    tng_data->current_trajectory_frame_set.prev_frame_set_file_pos = -1;
    tng_data->current_trajectory_frame_set.n_frames = 0;

    return(TNG_SUCCESS);
}

tng_function_status tng_trajectory_destroy(tng_trajectory_t *tng_data_p)
{
    int64_t n_particles;
    int i;
    tng_trajectory_t tng_data = *tng_data_p;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;

    tng_particle_mapping_t mapping;

    if(!*tng_data_p)
    {
        return(TNG_SUCCESS);
    }

    if(tng_data->input_file_path)
    {
        free(tng_data->input_file_path);
        tng_data->input_file_path = 0;
    }

    if(tng_data->input_file)
    {
        fclose(tng_data->input_file);
        tng_data->input_file = 0;
    }

    if(tng_data->output_file_path)
    {
        free(tng_data->output_file_path);
        tng_data->output_file_path = 0;
    }

    if(tng_data->output_file)
    {
        /* FIXME: Do not always write the hash */
        tng_frame_set_finalize(tng_data, TNG_USE_HASH);
        fclose(tng_data->output_file);
        tng_data->output_file = 0;
    }

    if(tng_data->first_program_name)
    {
        free(tng_data->first_program_name);
        tng_data->first_program_name = 0;
    }

    if(tng_data->last_program_name)
    {
        free(tng_data->last_program_name);
        tng_data->last_program_name = 0;
    }

    if(tng_data->first_user_name)
    {
        free(tng_data->first_user_name);
        tng_data->first_user_name = 0;
    }

    if(tng_data->last_user_name)
    {
        free(tng_data->last_user_name);
        tng_data->last_user_name = 0;
    }

    if(tng_data->first_computer_name)
    {
        free(tng_data->first_computer_name);
        tng_data->first_computer_name = 0;
    }

    if(tng_data->last_computer_name)
    {
        free(tng_data->last_computer_name);
        tng_data->last_computer_name = 0;
    }

    if(tng_data->first_pgp_signature)
    {
        free(tng_data->first_pgp_signature);
        tng_data->first_pgp_signature = 0;
    }

    if(tng_data->last_pgp_signature)
    {
        free(tng_data->last_pgp_signature);
        tng_data->last_pgp_signature = 0;
    }

    if(tng_data->forcefield_name)
    {
        free(tng_data->forcefield_name);
        tng_data->forcefield_name = 0;
    }

    if(frame_set->mappings)
    {
        for(i = frame_set->n_mapping_blocks; i--;)
        {
            mapping = &frame_set->mappings[i];
            if(mapping->real_particle_numbers)
            {
                free(mapping->real_particle_numbers);
                mapping->real_particle_numbers = 0;
            }
        }
        free(frame_set->mappings);
        frame_set->mappings = 0;
        frame_set->n_mapping_blocks = 0;
    }

    if(frame_set->molecule_cnt_list)
    {
        free(frame_set->molecule_cnt_list);
        frame_set->molecule_cnt_list = 0;
    }

    if(tng_data->var_num_atoms_flag)
    {
        n_particles = tng_data->current_trajectory_frame_set.n_particles;
    }
    else
    {
        n_particles = tng_data->n_particles;
    }

    if(tng_data->non_tr_particle_data)
    {
        for(i = tng_data->n_particle_data_blocks; i--; )
        {
            tng_particle_data_values_free(tng_data,
                                       tng_data->non_tr_particle_data[i].
                                       values, 1,
                                       n_particles,
                                       tng_data->non_tr_particle_data[i].
                                       n_values_per_frame,
                                       tng_data->non_tr_particle_data[i].
                                       datatype);

            if(tng_data->non_tr_particle_data[i].block_name)
            {
                free(tng_data->non_tr_particle_data[i].block_name);
                tng_data->non_tr_particle_data[i].block_name = 0;
            }
        }
        free(tng_data->non_tr_particle_data);
        tng_data->non_tr_particle_data = 0;
    }

    if(tng_data->non_tr_data)
    {
        for(i = tng_data->n_data_blocks; i--;)
        {
            tng_data_values_free(tng_data,
                                 tng_data->non_tr_data[i].values, 1,
                                 tng_data->non_tr_data[i].n_values_per_frame,
                                 tng_data->non_tr_data[i].datatype);

            if(tng_data->non_tr_data[i].block_name)
            {
                free(tng_data->non_tr_data[i].block_name);
                tng_data->non_tr_data[i].block_name = 0;
            }
        }
        free(tng_data->non_tr_data);
        tng_data->non_tr_data = 0;
    }

    tng_data->n_particle_data_blocks = 0;
    tng_data->n_data_blocks = 0;

    if(frame_set->tr_particle_data)
    {
        for(i = frame_set->n_particle_data_blocks; i--; )
        {
            tng_particle_data_values_free(tng_data,
                                          frame_set->tr_particle_data[i].
                                          values,
                                          frame_set->tr_particle_data[i].
                                          n_frames,
                                          n_particles,
                                          frame_set->tr_particle_data[i].
                                          n_values_per_frame,
                                          frame_set->tr_particle_data[i].
                                          datatype);

            if(frame_set->tr_particle_data[i].block_name)
            {
                free(frame_set->tr_particle_data[i].block_name);
                frame_set->tr_particle_data[i].block_name = 0;
            }
        }
        free(frame_set->tr_particle_data);
        frame_set->tr_particle_data = 0;
    }

    if(frame_set->tr_data)
    {
        for(i = frame_set->n_data_blocks; i--;)
        {
            tng_data_values_free(tng_data,
                                 frame_set->tr_data[i].values,
                                 frame_set->tr_data[i].n_frames,
                                 frame_set->tr_data[i].
                                 n_values_per_frame,
                                 frame_set->tr_data[i].datatype);
            if(frame_set->tr_data[i].block_name)
            {
                free(frame_set->tr_data[i].block_name);
                frame_set->tr_data[i].block_name = 0;
            }
        }
        free(frame_set->tr_data);
        frame_set->tr_data = 0;
    }

    frame_set->n_particle_data_blocks = 0;
    frame_set->n_data_blocks = 0;

    if(tng_data->molecules)
    {
        for(i=tng_data->n_molecules; i--;)
        {
            tng_molecule_destroy(tng_data, &tng_data->molecules[i]);
        }
        free(tng_data->molecules);
        tng_data->molecules = 0;
        tng_data->n_molecules = 0;
    }
    if(tng_data->molecule_cnt_list)
    {
        free(tng_data->molecule_cnt_list);
        tng_data->molecule_cnt_list = 0;
    }

    free(*tng_data_p);
    *tng_data_p = 0;

    return(TNG_SUCCESS);
}

tng_function_status tng_trajectory_init_from_src(tng_trajectory_t src,
                                                 tng_trajectory_t *dest_p)
{
    tng_trajectory_frame_set_t frame_set;
    tng_trajectory_t dest;

    *dest_p = malloc(sizeof(struct tng_trajectory));
    if(!dest_p)
    {
        printf("Cannot allocate memory (%lu bytes). %s: %d\n",
               sizeof(struct tng_trajectory), __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    dest = *dest_p;

    frame_set = &dest->current_trajectory_frame_set;

    dest->input_file_path = 0;
    if(src->input_file)
    {
        dest->input_file = fopen(src->input_file_path, "r");
    }
    else
    {
        dest->input_file = 0;
    }

    dest->input_file_len = src->input_file_len;
    dest->output_file_path = 0;
    if(src->output_file)
    {
        dest->output_file = fopen(src->output_file_path, "w+");
    }
    else
    {
        dest->output_file = 0;
    }

    dest->first_program_name = 0;
    dest->first_user_name = 0;
    dest->first_computer_name = 0;
    dest->first_pgp_signature = 0;
    dest->last_program_name = 0;
    dest->last_user_name = 0;
    dest->last_computer_name = 0;
    dest->last_pgp_signature = 0;
    dest->forcefield_name = 0;

    dest->var_num_atoms_flag = src->var_num_atoms_flag;
    dest->first_trajectory_frame_set_input_file_pos =
    src->first_trajectory_frame_set_input_file_pos;
    dest->last_trajectory_frame_set_input_file_pos =
    src->last_trajectory_frame_set_input_file_pos;
    dest->current_trajectory_frame_set_input_file_pos =
    src->current_trajectory_frame_set_input_file_pos;
    dest->first_trajectory_frame_set_output_file_pos =
    src->first_trajectory_frame_set_output_file_pos;
    dest->last_trajectory_frame_set_output_file_pos =
    src->last_trajectory_frame_set_output_file_pos;
    dest->current_trajectory_frame_set_output_file_pos =
    src->current_trajectory_frame_set_output_file_pos;
    dest->frame_set_n_frames = src->frame_set_n_frames;
    dest->n_trajectory_frame_sets = src->n_trajectory_frame_sets;
    dest->n_trajectory_blocks = src->n_trajectory_blocks;
    dest->medium_stride_length = src->medium_stride_length;
    dest->long_stride_length = src->long_stride_length;

    /* Currently the non trajectory data blocks are not copied since it
     * can lead to problems when freeing memory in a parallel block. */
    dest->n_particle_data_blocks = 0;
    dest->n_data_blocks = 0;
    dest->non_tr_particle_data = 0;
    dest->non_tr_data = 0;

    frame_set->first_frame = -1;
    frame_set->n_mapping_blocks = 0;
    frame_set->mappings = 0;
    frame_set->molecule_cnt_list = 0;

    frame_set->n_particle_data_blocks = 0;
    frame_set->n_data_blocks = 0;

    frame_set->tr_particle_data = 0;
    frame_set->tr_data = 0;

    frame_set->next_frame_set_file_pos = -1;
    frame_set->prev_frame_set_file_pos = -1;
    frame_set->medium_stride_next_frame_set_file_pos = -1;
    frame_set->medium_stride_prev_frame_set_file_pos = -1;
    frame_set->long_stride_next_frame_set_file_pos = -1;
    frame_set->long_stride_prev_frame_set_file_pos = -1;

    dest->n_molecules = 0;
    dest->molecules = 0;
    dest->molecule_cnt_list = 0;
    dest->n_particles = src->n_particles;

    dest->endianness_32 = src->endianness_32;
    dest->endianness_64 = src->endianness_64;
    dest->input_endianness_swap_func_32 = src->input_endianness_swap_func_32;
    dest->input_endianness_swap_func_64 = src->input_endianness_swap_func_64;
    dest->output_endianness_swap_func_32 = src->output_endianness_swap_func_32;
    dest->output_endianness_swap_func_64 = src->output_endianness_swap_func_64;

    dest->current_trajectory_frame_set.next_frame_set_file_pos = -1;
    dest->current_trajectory_frame_set.prev_frame_set_file_pos = -1;
    dest->current_trajectory_frame_set.n_frames = 0;

    return(TNG_SUCCESS);
}

tng_function_status tng_input_file_get(const tng_trajectory_t tng_data,
                                       char *file_name, const int max_len)
{
    strncpy(file_name, tng_data->input_file_path, max_len - 1);
    file_name[max_len - 1] = 0;

    if(strlen(tng_data->input_file_path) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_input_file_set(tng_trajectory_t tng_data,
                                       const char *file_name)
{
    int len;
    char *temp;

    if(tng_data->input_file_path && strcmp(tng_data->input_file_path,
                                           file_name) == 0)
    {
        return(TNG_SUCCESS);
    }

    if(tng_data->input_file)
    {
        fclose(tng_data->input_file);
    }

    len = tng_min(strlen(file_name) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->input_file_path, len);
    if(!temp)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->input_file_path);
        return(TNG_CRITICAL);
    }
    tng_data->input_file_path = temp;

    strncpy(tng_data->input_file_path, file_name, len);

    return(tng_input_file_init(tng_data));
}

tng_function_status tng_output_file_get(const tng_trajectory_t tng_data,
                                       char *file_name, const int max_len)
{
    strncpy(file_name, tng_data->output_file_path, max_len - 1);
    file_name[max_len - 1] = 0;

    if(strlen(tng_data->output_file_path) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_output_file_set(tng_trajectory_t tng_data,
                                        const char *file_name)
{
    int len;
    char *temp;

    if(tng_data->output_file_path &&
       strcmp(tng_data->output_file_path, file_name) == 0)
    {
        return(TNG_SUCCESS);
    }

    if(tng_data->output_file)
    {
        fclose(tng_data->output_file);
    }

    len = tng_min(strlen(file_name) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->output_file_path, len);
    if(!temp)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->output_file_path);
        return(TNG_CRITICAL);
    }
    tng_data->output_file_path = temp;

    strncpy(tng_data->output_file_path, file_name, len);

    return(tng_output_file_init(tng_data));
}

tng_function_status tng_output_file_endianness_get
                (tng_trajectory_t tng_data, tng_file_endianness *endianness)
{
    tng_endianness_32 end_32;
    tng_endianness_64 end_64;

    if(tng_data->output_endianness_swap_func_32)
    {
        /* If other endianness variants are added they must be added here as well */
        if(tng_data->output_endianness_swap_func_32 ==
           &tng_swap_byte_order_big_endian_32)
        {
            end_32 = TNG_BIG_ENDIAN_32;
        }
        else if(tng_data->output_endianness_swap_func_32 ==
                &tng_swap_byte_order_little_endian_32)
        {
            end_32 = TNG_LITTLE_ENDIAN_32;
        }
        else
        {
            return(TNG_FAILURE);
        }
    }
    else
    {
        end_32 = tng_data->endianness_32;
    }

    if(tng_data->output_endianness_swap_func_64)
    {
        /* If other endianness variants are added they must be added here as well */
        if(tng_data->output_endianness_swap_func_64 ==
           &tng_swap_byte_order_big_endian_64)
        {
            end_64 = TNG_BIG_ENDIAN_64;
        }
        else if(tng_data->output_endianness_swap_func_64 ==
                &tng_swap_byte_order_little_endian_64)
        {
            end_64 = TNG_LITTLE_ENDIAN_64;
        }
        else
        {
            return(TNG_FAILURE);
        }
    }
    else
    {
        end_64 = tng_data->endianness_64;
    }

    if((int)end_32 != (int)end_64)
    {
        return(TNG_FAILURE);
    }

    if(end_32 == TNG_LITTLE_ENDIAN_32)
    {
        *endianness = TNG_LITTLE_ENDIAN;
    }

    else if(end_32 == TNG_BIG_ENDIAN_32)
    {
        *endianness = TNG_BIG_ENDIAN;
    }
    else
    {
        return(TNG_FAILURE);
    }

    return(TNG_SUCCESS);
}

tng_function_status tng_output_file_endianness_set
                (tng_trajectory_t tng_data,
                 const tng_file_endianness endianness)
{
    /* Tne endianness cannot be changed if the data has already been written
     * to the output file. */
    if(ftell(tng_data->output_file) > 0)
    {
        return(TNG_FAILURE);
    }

    if(endianness == TNG_BIG_ENDIAN)
    {
        if(tng_data->endianness_32 == TNG_BIG_ENDIAN_32)
        {
            tng_data->output_endianness_swap_func_32 = 0;
        }
        else
        {
            tng_data->output_endianness_swap_func_32 =
            &tng_swap_byte_order_big_endian_32;
        }
        if(tng_data->endianness_64 == TNG_BIG_ENDIAN_64)
        {
            tng_data->output_endianness_swap_func_64 = 0;
        }
        else
        {
            tng_data->output_endianness_swap_func_64 =
            &tng_swap_byte_order_big_endian_64;
        }
        return(TNG_SUCCESS);
    }
    else if(endianness == TNG_LITTLE_ENDIAN)
    {
        if(tng_data->endianness_32 == TNG_LITTLE_ENDIAN_32)
        {
            tng_data->output_endianness_swap_func_32 = 0;
        }
        else
        {
            tng_data->output_endianness_swap_func_32 =
            &tng_swap_byte_order_little_endian_32;
        }
        if(tng_data->endianness_64 == TNG_LITTLE_ENDIAN_64)
        {
            tng_data->output_endianness_swap_func_64 = 0;
        }
        else
        {
            tng_data->output_endianness_swap_func_64 =
            &tng_swap_byte_order_little_endian_64;
        }
        return(TNG_SUCCESS);
    }

    /* If the specified endianness is neither big nor little endian return a
     * failure. */
    return(TNG_FAILURE);
}

tng_function_status tng_first_program_name_get(const tng_trajectory_t tng_data,
                                               char *name, const int max_len)
{
    strncpy(name, tng_data->first_program_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->first_program_name) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_first_program_name_set(tng_trajectory_t tng_data,
                                               const char *new_name)
{
    int len;

    len = tng_min(strlen(new_name) + 1, TNG_MAX_STR_LEN);

    if(tng_data->first_program_name && strlen(tng_data->first_program_name) < len)
    {
        free(tng_data->first_program_name);
        tng_data->first_program_name = 0;
    }
    if(!tng_data->first_program_name)
    {
        tng_data->first_program_name = malloc(len);
        if(!tng_data->first_program_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->first_program_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_last_program_name_get(const tng_trajectory_t tng_data,
                                              char *name, const int max_len)
{
    strncpy(name, tng_data->last_program_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->last_program_name) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_last_program_name_set(tng_trajectory_t tng_data,
                                              const char *new_name)
{
    int len;

    len = tng_min(strlen(new_name) + 1, TNG_MAX_STR_LEN);

    if(tng_data->last_program_name && strlen(tng_data->last_program_name) < len)
    {
        free(tng_data->last_program_name);
        tng_data->last_program_name = 0;
    }
    if(!tng_data->last_program_name)
    {
        tng_data->last_program_name = malloc(len);
        if(!tng_data->last_program_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->last_program_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_first_user_name_get(const tng_trajectory_t tng_data,
                                            char *name, const int max_len)
{
    strncpy(name, tng_data->first_user_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->first_user_name) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_first_user_name_set(tng_trajectory_t tng_data,
                                            const char *new_name)
{
    int len;

    len = tng_min(strlen(new_name) + 1, TNG_MAX_STR_LEN);

    /* If the currently stored string length is not enough to store the new
     * string it is freed and reallocated. */
    if(tng_data->first_user_name && strlen(tng_data->first_user_name) < len)
    {
        free(tng_data->first_user_name);
        tng_data->first_user_name = 0;
    }
    if(!tng_data->first_user_name)
    {
        tng_data->first_user_name = malloc(len);
        if(!tng_data->first_user_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->first_user_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_last_user_name_get(const tng_trajectory_t tng_data,
                                           char *name, const int max_len)
{
    strncpy(name, tng_data->last_user_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->last_user_name) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_last_user_name_set(tng_trajectory_t tng_data,
                                           const char *new_name)
{
    int len;

    len = tng_min(strlen(new_name) + 1, TNG_MAX_STR_LEN);

    /* If the currently stored string length is not enough to store the new
     * string it is freed and reallocated. */
    if(tng_data->last_user_name && strlen(tng_data->last_user_name) < len)
    {
        free(tng_data->last_user_name);
        tng_data->last_user_name = 0;
    }
    if(!tng_data->last_user_name)
    {
        tng_data->last_user_name = malloc(len);
        if(!tng_data->last_user_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->last_user_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_first_computer_name_get(const tng_trajectory_t tng_data,
                                                char *name, const int max_len)
{
    strncpy(name, tng_data->first_computer_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->first_computer_name) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_first_computer_name_set(tng_trajectory_t tng_data,
                                                const char *new_name)
{
    int len;

    len = tng_min(strlen(new_name) + 1, TNG_MAX_STR_LEN);

    /* If the currently stored string length is not enough to store the new
     * string it is freed and reallocated. */
    if(tng_data->first_computer_name && strlen(tng_data->first_computer_name) < len)
    {
        free(tng_data->first_computer_name);
        tng_data->first_computer_name = 0;
    }
    if(!tng_data->first_computer_name)
    {
        tng_data->first_computer_name = malloc(len);
        if(!tng_data->first_computer_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->first_computer_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_last_computer_name_get(const tng_trajectory_t tng_data,
                                               char *name, const int max_len)
{
    strncpy(name, tng_data->last_computer_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->last_computer_name) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_last_computer_name_set(tng_trajectory_t tng_data,
                                               const char *new_name)
{
    int len;

    len = tng_min(strlen(new_name) + 1, TNG_MAX_STR_LEN);

    /* If the currently stored string length is not enough to store the new
     * string it is freed and reallocated. */
    if(tng_data->last_computer_name && strlen(tng_data->last_computer_name) <
        len)
    {
        free(tng_data->last_computer_name);
        tng_data->last_computer_name = 0;
    }
    if(!tng_data->last_computer_name)
    {
        tng_data->last_computer_name = malloc(len);
        if(!tng_data->last_computer_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->last_computer_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_first_signature_get(const tng_trajectory_t tng_data,
                                            char *signature, const int max_len)
{
    strncpy(signature, tng_data->first_pgp_signature, max_len - 1);
    signature[max_len - 1] = 0;

    if(strlen(tng_data->first_pgp_signature) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_first_signature_set(tng_trajectory_t tng_data,
                                            const char *signature)
{
    int len;

    len = tng_min(strlen(signature) + 1, TNG_MAX_STR_LEN);

    /* If the currently stored string length is not enough to store the new
     * string it is freed and reallocated. */
    if(tng_data->first_pgp_signature && strlen(tng_data->first_pgp_signature) <
        len)
    {
        free(tng_data->first_pgp_signature);
        tng_data->first_pgp_signature = 0;
    }
    if(!tng_data->first_pgp_signature)
    {
        tng_data->first_pgp_signature = malloc(len);
        if(!tng_data->first_pgp_signature)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->first_pgp_signature, signature, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_last_signature_get(const tng_trajectory_t tng_data,
                                           char *signature, const int max_len)
{
    strncpy(signature, tng_data->last_pgp_signature, max_len - 1);
    signature[max_len - 1] = 0;

    if(strlen(tng_data->last_pgp_signature) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_last_signature_set(tng_trajectory_t tng_data,
                                           const char *signature)
{
    int len;

    len = tng_min(strlen(signature) + 1, TNG_MAX_STR_LEN);

    /* If the currently stored string length is not enough to store the new
     * string it is freed and reallocated. */
    if(tng_data->last_pgp_signature && strlen(tng_data->last_pgp_signature) <
        len)
    {
        free(tng_data->last_pgp_signature);
        tng_data->last_pgp_signature = 0;
    }
    if(!tng_data->last_pgp_signature)
    {
        tng_data->last_pgp_signature = malloc(len);
        if(!tng_data->last_pgp_signature)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->last_pgp_signature, signature, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_forcefield_name_get(const tng_trajectory_t tng_data,
                                            char *name, const int max_len)
{
    strncpy(name, tng_data->forcefield_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->forcefield_name) > max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_forcefield_name_set(tng_trajectory_t tng_data,
                                            const char *new_name)
{
    int len;

    len = tng_min(strlen(new_name) + 1, TNG_MAX_STR_LEN);

    /* If the currently stored string length is not enough to store the new
     * string it is freed and reallocated. */
    if(tng_data->forcefield_name && strlen(tng_data->forcefield_name) < len)
    {
        free(tng_data->forcefield_name);
        tng_data->forcefield_name = 0;
    }
    if(!tng_data->forcefield_name)
    {
        tng_data->forcefield_name = malloc(len);
        if(!tng_data->forcefield_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->forcefield_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_medium_stride_length_get(const tng_trajectory_t tng_data,
                                                 int64_t *len)
{
    *len = tng_data->medium_stride_length;

    return(TNG_SUCCESS);
}

tng_function_status tng_medium_stride_length_set(tng_trajectory_t tng_data,
                                                 const int64_t len)
{
    if(len >= tng_data->long_stride_length)
    {
        return(TNG_FAILURE);
    }
    tng_data->medium_stride_length = len;

    return(TNG_SUCCESS);
}

tng_function_status tng_long_stride_length_get(const tng_trajectory_t tng_data,
                                               int64_t *len)
{
    *len = tng_data->long_stride_length;

    return(TNG_SUCCESS);
}

tng_function_status tng_long_stride_length_set(tng_trajectory_t tng_data,
                                               const int64_t len)
{
    if(len <= tng_data->medium_stride_length)
    {
        return(TNG_FAILURE);
    }
    tng_data->long_stride_length = len;

    return(TNG_SUCCESS);
}

tng_function_status tng_input_file_len_get(const tng_trajectory_t tng_data,
                                           int64_t *len)
{
    *len = tng_data->input_file_len;

    return(TNG_SUCCESS);
}

tng_function_status tng_num_frames_get(const tng_trajectory_t tng_data,
                                       int64_t *n)
{
    tng_gen_block_t block;
    tng_function_status stat;
    int64_t file_pos;

    file_pos = tng_data->last_trajectory_frame_set_input_file_pos;

    if(file_pos <= 0)
    {
        return(TNG_FAILURE);
    }

    tng_block_init(&block);
    fseek(tng_data->input_file,
        file_pos,
        SEEK_SET);
    tng_data->current_trajectory_frame_set_input_file_pos = file_pos;
    /* Read block headers first to see what block is found. */
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        printf("Cannot read block header at pos %"PRId64". %s: %d\n", file_pos,
                __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_FAILURE);
    }

    stat = tng_block_read_next(tng_data, block,
                               TNG_SKIP_HASH);
    tng_block_destroy(&block);

    if(stat != TNG_SUCCESS)
    {
        return(TNG_FAILURE);
    }

    *n = tng_data->current_trajectory_frame_set.first_frame +
         tng_data->current_trajectory_frame_set.n_frames;

    return(TNG_SUCCESS);
}

tng_function_status tng_num_particles_get(const tng_trajectory_t tng_data,
                                          int64_t *n)
{
    if(tng_data->var_num_atoms_flag == TNG_CONSTANT_N_ATOMS)
    {
        *n = tng_data->n_particles;
    }
    else
    {
        *n = tng_data->current_trajectory_frame_set.n_particles;
    }

    return(TNG_SUCCESS);
}

tng_function_status tng_num_molecules_get(const tng_trajectory_t tng_data,
                                          int64_t *n)
{
    int64_t *cnt_list, cnt = 0, i;

    if(tng_data->var_num_atoms_flag == TNG_CONSTANT_N_ATOMS)
    {
        cnt_list = tng_data->molecule_cnt_list;
    }
    else
    {
        cnt_list = tng_data->current_trajectory_frame_set.molecule_cnt_list;
    }

    for(i = tng_data->n_molecules; i --;)
    {
        cnt += cnt_list[i];
    }

    *n = cnt;

    return(TNG_SUCCESS);
}

tng_function_status tng_num_frames_per_frame_set_get
                (const tng_trajectory_t tng_data,
                 int64_t *n)
{
    *n = tng_data->frame_set_n_frames;

    return(TNG_SUCCESS);
}

tng_function_status tng_num_frames_per_frame_set_set
                (const tng_trajectory_t tng_data,
                 const int64_t n)
{
    tng_data->frame_set_n_frames = n;

    return(TNG_SUCCESS);
}

tng_function_status tng_num_frame_sets_get(const tng_trajectory_t tng_data,
                                           int64_t *n)
{
    int64_t long_stride_length, medium_stride_length;
    int64_t file_pos;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;
    tng_gen_block_t block;
    tng_function_status stat;
    int64_t cnt = 0;

    file_pos = tng_data->first_trajectory_frame_set_input_file_pos;

    tng_block_init(&block);
    fseek(tng_data->input_file,
        file_pos,
        SEEK_SET);
    tng_data->current_trajectory_frame_set_input_file_pos = file_pos;
    /* Read block headers first to see what block is found. */
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        printf("Cannot read block header at pos %"PRId64". %s: %d\n", file_pos,
                __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    if(tng_block_read_next(tng_data, block,
                        TNG_SKIP_HASH) != TNG_SUCCESS)
    {
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    ++cnt;

    file_pos = tng_data->current_trajectory_frame_set_input_file_pos;

    long_stride_length = tng_data->long_stride_length;
    medium_stride_length = tng_data->medium_stride_length;

    /* Take long steps forward until a long step forward would be too long or
     * the right frame set is found */
    file_pos = frame_set->long_stride_next_frame_set_file_pos;
    while(file_pos > 0)
    {
        if(file_pos > 0)
        {
            cnt += long_stride_length;
            fseek(tng_data->input_file, file_pos, SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
        }
        file_pos = frame_set->long_stride_next_frame_set_file_pos;
    }

    /* Take medium steps forward until a medium step forward would be too long
     * or the right frame set is found */
    file_pos = frame_set->medium_stride_next_frame_set_file_pos;
    while(file_pos > 0)
    {
        if(file_pos > 0)
        {
            cnt += medium_stride_length;
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
        }
        file_pos = frame_set->medium_stride_next_frame_set_file_pos;
    }

    /* Take one step forward until the right frame set is found */
    file_pos = frame_set->next_frame_set_file_pos;
    while(file_pos > 0)
    {
        if(file_pos > 0)
        {
            ++cnt;
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
        }
        file_pos = frame_set->next_frame_set_file_pos;
    }

    tng_block_destroy(&block);

    *n = tng_data->n_trajectory_frame_sets = cnt;

    return(TNG_SUCCESS);
}

tng_function_status tng_current_frame_set_get
                (tng_trajectory_t tng_data,
                 tng_trajectory_frame_set_t *frame_set_p)
{
    *frame_set_p = &tng_data->current_trajectory_frame_set;

    return(TNG_SUCCESS);
}

tng_function_status tng_frame_set_nr_find(tng_trajectory_t tng_data,
                                       const int64_t nr)
{
    int64_t long_stride_length, medium_stride_length;
    int64_t file_pos, curr_nr = 0, n_frame_sets;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;
    tng_gen_block_t block;
    tng_function_status stat;

    stat = tng_num_frame_sets_get(tng_data, &n_frame_sets);

    if(nr >= n_frame_sets)
    {
        return(TNG_FAILURE);
    }

    long_stride_length = tng_data->long_stride_length;
    medium_stride_length = tng_data->medium_stride_length;

    /* FIXME: The frame set number of the current frame set is not stored */

    if(nr < n_frame_sets - 1 - nr)
    {
        /* Start from the beginning */
        file_pos = tng_data->first_trajectory_frame_set_input_file_pos;
    }
    else
    {
        /* Start from the end */
        file_pos = tng_data->last_trajectory_frame_set_input_file_pos;
        curr_nr = n_frame_sets - 1;
    }
    if(file_pos <= 0)
    {
        return(TNG_FAILURE);
    }

    tng_block_init(&block);
    fseek(tng_data->input_file,
        file_pos,
        SEEK_SET);
    tng_data->current_trajectory_frame_set_input_file_pos = file_pos;
    /* Read block headers first to see what block is found. */
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        printf("Cannot read block header at pos %"PRId64". %s: %d\n", file_pos,
                __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    if(tng_block_read_next(tng_data, block,
                        TNG_SKIP_HASH) != TNG_SUCCESS)
    {
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    if(curr_nr == nr)
    {
        tng_block_destroy(&block);
        return(TNG_SUCCESS);
    }

    file_pos = tng_data->current_trajectory_frame_set_input_file_pos;

    /* Take long steps forward until a long step forward would be too long or
     * the right frame set is found */
    while(file_pos > 0 && curr_nr + long_stride_length <= nr)
    {
        file_pos = frame_set->long_stride_next_frame_set_file_pos;
        if(file_pos > 0)
        {
            curr_nr += long_stride_length;
            fseek(tng_data->input_file, file_pos, SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos,  __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
            if(curr_nr == nr)
            {
                tng_block_destroy(&block);
                return(TNG_SUCCESS);
            }
        }
    }

    /* Take medium steps forward until a medium step forward would be too long
     * or the right frame set is found */
    while(file_pos > 0 && curr_nr + medium_stride_length <= nr)
    {
        file_pos = frame_set->medium_stride_next_frame_set_file_pos;
        if(file_pos > 0)
        {
            curr_nr += medium_stride_length;
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
            if(curr_nr == nr)
            {
                tng_block_destroy(&block);
                return(TNG_SUCCESS);
            }
        }
    }

    /* Take one step forward until the right frame set is found */
    while(file_pos > 0 && curr_nr < nr)
    {
        file_pos = frame_set->next_frame_set_file_pos;

        if(file_pos > 0)
        {
            ++curr_nr;
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
            if(curr_nr == nr)
            {
                tng_block_destroy(&block);
                return(TNG_SUCCESS);
            }
        }
    }

    /* Take long steps backward until a long step backward would be too long
     * or the right frame set is found */
    while(file_pos > 0 && curr_nr - long_stride_length >= nr)
    {
        file_pos = frame_set->long_stride_prev_frame_set_file_pos;
        if(file_pos > 0)
        {
            curr_nr -= long_stride_length;
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
            if(curr_nr == nr)
            {
                tng_block_destroy(&block);
                return(TNG_SUCCESS);
            }
        }
    }

    /* Take medium steps backward until a medium step backward would be too long
     * or the right frame set is found */
    while(file_pos > 0 && curr_nr - medium_stride_length >= nr)
    {
        file_pos = frame_set->medium_stride_prev_frame_set_file_pos;
        if(file_pos > 0)
        {
            curr_nr -= medium_stride_length;
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
            if(curr_nr == nr)
            {
                tng_block_destroy(&block);
                return(TNG_SUCCESS);
            }
        }
    }

    /* Take one step backward until the right frame set is found */
    while(file_pos > 0 && curr_nr > nr)
    {
        file_pos = frame_set->prev_frame_set_file_pos;
        if(file_pos > 0)
        {
            --curr_nr;
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
            if(curr_nr == nr)
            {
                tng_block_destroy(&block);
                return(TNG_SUCCESS);
            }
        }
    }

    /* If for some reason the current frame set is not yet found,
     * take one step forward until the right frame set is found */
    while(file_pos > 0 && curr_nr < nr)
    {
        file_pos = frame_set->next_frame_set_file_pos;
        if(file_pos > 0)
        {
            ++curr_nr;
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
            if(curr_nr == nr)
            {
                tng_block_destroy(&block);
                return(TNG_SUCCESS);
            }
        }
    }

    tng_block_destroy(&block);
    return(TNG_FAILURE);
}

tng_function_status tng_frame_set_of_frame_find(tng_trajectory_t tng_data,
                                       const int64_t frame)
{
    int64_t first_frame, last_frame, n_frames_per_frame_set;
    int64_t long_stride_length, medium_stride_length;
    int64_t file_pos;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;
    tng_gen_block_t block;
    tng_function_status stat;

    first_frame = tng_max(frame_set->first_frame, 0);
    last_frame = first_frame + frame_set->n_frames - 1;
    n_frames_per_frame_set = tng_data->frame_set_n_frames;
    long_stride_length = tng_data->long_stride_length;
    medium_stride_length = tng_data->medium_stride_length;
    tng_block_init(&block);

    /* Is this the right frame set? */
    if(first_frame <= frame && frame <= last_frame)
    {
        tng_block_destroy(&block);
        return(TNG_SUCCESS);
    }

    if(first_frame - frame >= frame ||
       frame - last_frame >
       tng_data->n_trajectory_frame_sets * n_frames_per_frame_set - frame)
    {
        /* Start from the beginning */
        if(first_frame - frame >= frame)
        {
            file_pos = tng_data->first_trajectory_frame_set_input_file_pos;

            if(file_pos <= 0)
            {
                tng_block_destroy(&block);
                return(TNG_FAILURE);
            }
        }
        /* Start from the end */
        else
        {
            file_pos = tng_data->last_trajectory_frame_set_input_file_pos;

            /* If the last frame set position is not set start from the current
             * frame set, since it will be closer than the first frame set. */
        }

        if(file_pos > 0)
        {
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            tng_data->current_trajectory_frame_set_input_file_pos = file_pos;
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
        }
    }

    first_frame = tng_max(frame_set->first_frame, 0);
    last_frame = first_frame + frame_set->n_frames - 1;

    if(frame >= first_frame && frame <= last_frame)
    {
        tng_block_destroy(&block);
        return(TNG_SUCCESS);
    }

    file_pos = tng_data->current_trajectory_frame_set_input_file_pos;

    /* Take long steps forward until a long step forward would be too long or
     * the right frame set is found */
    while(file_pos > 0 && first_frame + long_stride_length *
          n_frames_per_frame_set <= frame)
    {
        file_pos = frame_set->long_stride_next_frame_set_file_pos;
        if(file_pos > 0)
        {
            fseek(tng_data->input_file, file_pos, SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
        }
        first_frame = tng_max(frame_set->first_frame, 0);
        last_frame = first_frame + frame_set->n_frames - 1;
        if(frame >= first_frame && frame <= last_frame)
        {
            tng_block_destroy(&block);
            return(TNG_SUCCESS);
        }
    }

    /* Take medium steps forward until a medium step forward would be too long
     * or the right frame set is found */
    while(file_pos > 0 && first_frame + medium_stride_length *
          n_frames_per_frame_set <= frame)
    {
        file_pos = frame_set->medium_stride_next_frame_set_file_pos;
        if(file_pos > 0)
        {
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
        }
        first_frame = tng_max(frame_set->first_frame, 0);
        last_frame = first_frame + frame_set->n_frames - 1;
        if(frame >= first_frame && frame <= last_frame)
        {
            tng_block_destroy(&block);
            return(TNG_SUCCESS);
        }
    }

    /* Take one step forward until the right frame set is found */
    while(file_pos > 0 && first_frame < frame && last_frame < frame)
    {
        file_pos = frame_set->next_frame_set_file_pos;
        if(file_pos > 0)
        {
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
        }
        first_frame = tng_max(frame_set->first_frame, 0);
        last_frame = first_frame + frame_set->n_frames - 1;
        if(frame >= first_frame && frame <= last_frame)
        {
            tng_block_destroy(&block);
            return(TNG_SUCCESS);
        }
    }

    /* Take long steps backward until a long step backward would be too long
     * or the right frame set is found */
    while(file_pos > 0 && first_frame - long_stride_length *
          n_frames_per_frame_set >= frame)
    {
        file_pos = frame_set->long_stride_prev_frame_set_file_pos;
        if(file_pos > 0)
        {
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
        }
        first_frame = tng_max(frame_set->first_frame, 0);
        last_frame = first_frame + frame_set->n_frames - 1;
        if(frame >= first_frame && frame <= last_frame)
        {
            tng_block_destroy(&block);
            return(TNG_SUCCESS);
        }
    }

    /* Take medium steps backward until a medium step backward would be too long
     * or the right frame set is found */
    while(file_pos > 0 && first_frame - medium_stride_length *
          n_frames_per_frame_set >= frame)
    {
        file_pos = frame_set->medium_stride_prev_frame_set_file_pos;
        if(file_pos > 0)
        {
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
        }
        first_frame = tng_max(frame_set->first_frame, 0);
        last_frame = first_frame + frame_set->n_frames - 1;
        if(frame >= first_frame && frame <= last_frame)
        {
            tng_block_destroy(&block);
            return(TNG_SUCCESS);
        }
    }

    /* Take one step backward until the right frame set is found */
    while(file_pos > 0 && first_frame > frame && last_frame > frame)
    {
        file_pos = frame_set->prev_frame_set_file_pos;
        if(file_pos > 0)
        {
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
        }
        first_frame = tng_max(frame_set->first_frame, 0);
        last_frame = first_frame + frame_set->n_frames - 1;
        if(frame >= first_frame && frame <= last_frame)
        {
            tng_block_destroy(&block);
            return(TNG_SUCCESS);
        }
    }

    /* If for some reason the current frame set is not yet found,
     * take one step forward until the right frame set is found */
    while(file_pos > 0 && first_frame < frame && last_frame < frame)
    {
        file_pos = frame_set->next_frame_set_file_pos;
        if(file_pos > 0)
        {
            fseek(tng_data->input_file,
                file_pos,
                SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_block_read_next(tng_data, block,
                                TNG_SKIP_HASH) != TNG_SUCCESS)
            {
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }
        }
        first_frame = tng_max(frame_set->first_frame, 0);
        last_frame = first_frame + frame_set->n_frames - 1;
        if(frame >= first_frame && frame <= last_frame)
        {
            tng_block_destroy(&block);
            return(TNG_SUCCESS);
        }
    }

    tng_block_destroy(&block);
    return(TNG_FAILURE);
}

tng_function_status tng_frame_set_next_frame_set_file_pos_get
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *pos)
{
    *pos = frame_set->next_frame_set_file_pos;

    return(TNG_SUCCESS);
}

tng_function_status tng_frame_set_prev_frame_set_file_pos_get
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *pos)
{
    *pos = frame_set->prev_frame_set_file_pos;

    return(TNG_SUCCESS);
}

tng_function_status tng_frame_set_frame_range_get
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *first_frame,
                 int64_t *last_frame)
{
    *first_frame = frame_set->first_frame;
    *last_frame = *first_frame + frame_set->n_frames - 1;

    return(TNG_SUCCESS);
}

/** Translate from the particle numbering used in a frame set to the real
 *  particle numbering - used in the molecule description.
 * @param frame_set is the frame_set containing the mappings to use.
 * @param local is the index number of the atom in this frame set
 * @param real is set to the index of the atom in the molecular system.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (1) if the mapping
 * cannot be found.
 */
static inline tng_function_status tng_particle_mapping_get_real_particle
                (const tng_trajectory_frame_set_t frame_set,
                 const int64_t local,
                 int64_t *real)
{
    int64_t i, n_blocks = frame_set->n_mapping_blocks, first;
    tng_particle_mapping_t mapping;
    if(n_blocks <= 0)
    {
        *real = local;
        return(TNG_SUCCESS);
    }
    for(i = 0; i < n_blocks; i++)
    {
        mapping = &frame_set->mappings[i];
        first = mapping->num_first_particle;
        if(local < first ||
           local >= first + mapping->n_particles)
        {
            continue;
        }
        *real = mapping->real_particle_numbers[local-first];
        return(TNG_SUCCESS);
    }
    *real = local;
    return(TNG_FAILURE);
}

/** Translate from the real particle numbering to the particle numbering
 *  used in a frame set.
 * @param frame_set is the frame_set containing the mappings to use.
 * @param real is the index number of the atom in the molecular system.
 * @param local is set to the index of the atom in this frame set.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (1) if the mapping
 * cannot be found.
 */
static inline tng_function_status tng_particle_mapping_get_local_particle
                (const tng_trajectory_frame_set_t frame_set,
                 const int64_t real,
                 int64_t *local)
{
    int64_t i, j, n_blocks = frame_set->n_mapping_blocks;
    tng_particle_mapping_t mapping;
    if(n_blocks <= 0)
    {
        *local = real;
        return(TNG_SUCCESS);
    }
    for(i = 0; i < n_blocks; i++)
    {
        mapping = &frame_set->mappings[i];
        for(j = mapping->n_particles; j--;)
        {
            if(mapping->real_particle_numbers[j] == real)
            {
                *local = j;
                return(TNG_SUCCESS);
            }
        }
    }
    return(TNG_FAILURE);
}


tng_function_status tng_file_headers_read(tng_trajectory_t tng_data,
                                          const tng_hash_mode hash_mode)
{
    int cnt = 0, prev_pos = 0;
    tng_gen_block_t block;

    tng_data->n_trajectory_frame_sets = 0;

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    if(!tng_data->input_file_len)
    {
        fseek(tng_data->input_file, 0, SEEK_END);
        tng_data->input_file_len = ftell(tng_data->input_file);
        fseek(tng_data->input_file, 0, SEEK_SET);
    }

    tng_block_init(&block);
    /* Non trajectory blocks (they come before the trajectory
     * blocks in the file) */
    while (prev_pos < tng_data->input_file_len &&
           tng_block_header_read(tng_data, block) != TNG_CRITICAL &&
           block->id != -1 &&
           block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        if(tng_block_read_next(tng_data, block,
                               hash_mode) == TNG_SUCCESS)
        {
            cnt++;
        }
        prev_pos = ftell(tng_data->input_file);
    }

    /* Go back if a trajectory block was encountered */
    if(block->id == TNG_TRAJECTORY_FRAME_SET)
    {
        fseek(tng_data->input_file, prev_pos, SEEK_SET);
    }

    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

tng_function_status tng_file_headers_write(tng_trajectory_t tng_data,
                                           const tng_hash_mode hash_mode)
{
    int i;
    tng_gen_block_t data_block;

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }


    if(tng_general_info_block_write(tng_data, hash_mode)
       != TNG_SUCCESS)
    {
        printf("Error writing general info block of file %s. %s: %d\n",
                tng_data->input_file_path, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(tng_molecules_block_write(tng_data, hash_mode)
        != TNG_SUCCESS)
    {
        printf("Error writing atom names block of file %s. %s: %d\n",
                tng_data->input_file_path, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* FIXME: Currently writing non-trajectory data blocks here.
     * Should perhaps be moved. */
    tng_block_init(&data_block);
    for(i = 0; i < tng_data->n_data_blocks; i++)
    {
        data_block->id = tng_data->non_tr_data[i].block_id;
        tng_data_block_write(tng_data, data_block,
                             i, hash_mode);
    }

    for(i = 0; i < tng_data->n_particle_data_blocks; i++)
    {
        data_block->id = tng_data->non_tr_particle_data[i].block_id;
        tng_particle_data_block_write(tng_data, data_block,
                                      i, 0, hash_mode);
    }

    tng_block_destroy(&data_block);

    return(TNG_SUCCESS);
}

tng_function_status tng_block_read_next(tng_trajectory_t tng_data,
                                        tng_gen_block_t block,
                                        const tng_hash_mode hash_mode)
{
    switch(block->id)
    {
    case TNG_TRAJECTORY_FRAME_SET:
        return(tng_frame_set_block_read(tng_data, block, hash_mode));
    case TNG_PARTICLE_MAPPING:
        return(tng_trajectory_mapping_block_read(tng_data, block, hash_mode));
    case TNG_GENERAL_INFO:
        return(tng_general_info_block_read(tng_data, block, hash_mode));
    case TNG_MOLECULES:
        return(tng_molecules_block_read(tng_data, block, hash_mode));
    default:
        if(block->id >= TNG_TRAJ_BOX_SHAPE)
        {
            return(tng_data_block_contents_read(tng_data, block, hash_mode));
        }
        else
        {
            /* Skip to the next block */
            fseek(tng_data->input_file, block->block_contents_size, SEEK_CUR);
            return(TNG_FAILURE);
        }
    }
}


tng_function_status tng_frame_set_read_next(tng_trajectory_t tng_data,
                                            const tng_hash_mode hash_mode)
{
    long int file_pos;
    tng_gen_block_t block;
    tng_function_status stat = TNG_SUCCESS;

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    file_pos = tng_data->current_trajectory_frame_set.next_frame_set_file_pos;

    if(file_pos < 0 && tng_data->current_trajectory_frame_set_input_file_pos <= 0)
    {
        file_pos = tng_data->first_trajectory_frame_set_input_file_pos;
    }

    if(file_pos > 0)
    {
        fseek(tng_data->input_file,
              file_pos,
              SEEK_SET);
    }
    else
    {
        return(TNG_FAILURE);
    }

    tng_block_init(&block);

    if(!tng_data->input_file_len)
    {
        fseek(tng_data->input_file, 0, SEEK_END);
        tng_data->input_file_len = ftell(tng_data->input_file);
        fseek(tng_data->input_file, file_pos, SEEK_SET);
    }

    /* Read block headers first to see what block is found. */
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        printf("Cannot read block header at pos %"PRId64". %s: %d\n",
               file_pos, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    tng_data->current_trajectory_frame_set_input_file_pos = file_pos;

    if(tng_block_read_next(tng_data, block,
                           hash_mode) == TNG_SUCCESS)
    {
        tng_data->n_trajectory_frame_sets++;
        file_pos = ftell(tng_data->input_file);
        /* Read all blocks until next frame set block */
        stat = tng_block_header_read(tng_data, block);
        while(file_pos < tng_data->input_file_len &&
              stat != TNG_CRITICAL &&
              block->id != TNG_TRAJECTORY_FRAME_SET)
        {
            stat = tng_block_read_next(tng_data, block,
                                       hash_mode);
            if(stat != TNG_CRITICAL)
            {
                file_pos = ftell(tng_data->input_file);
                if(file_pos < tng_data->input_file_len)
                {
                    stat = tng_block_header_read(tng_data, block);
                }
            }
        }
        if(stat == TNG_CRITICAL)
        {
            printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                   file_pos, __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(stat);
        }

        if(block->id == TNG_TRAJECTORY_FRAME_SET)
        {
            fseek(tng_data->input_file, file_pos, SEEK_SET);
        }
    }

    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

tng_function_status tng_frame_set_write(tng_trajectory_t tng_data,
                                        const tng_hash_mode hash_mode)
{
    int i, j;
    tng_gen_block_t block;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;

    tng_function_status stat;

    tng_data->current_trajectory_frame_set_output_file_pos =
    tng_data->last_trajectory_frame_set_output_file_pos =
    ftell(tng_data->output_file);

    if(tng_data->current_trajectory_frame_set_output_file_pos <= 0)
    {
        return(TNG_FAILURE);
    }

    if(tng_data->first_trajectory_frame_set_output_file_pos == -1)
    {
        tng_data->first_trajectory_frame_set_output_file_pos =
        tng_data->current_trajectory_frame_set_output_file_pos;
    }

    tng_block_init(&block);

    if(tng_frame_set_block_write(tng_data, block, hash_mode) != TNG_SUCCESS)
    {
        tng_block_destroy(&block);
        return(TNG_FAILURE);
    }

    /* Write non-particle data blocks */
    for(i = 0; i<frame_set->n_data_blocks; i++)
    {
        block->id = frame_set->tr_data[i].block_id;
        tng_data_block_write(tng_data, block, i, hash_mode);
    }
    /* Write the mapping blocks and particle data blocks*/
    if(frame_set->n_mapping_blocks)
    {
        for(i = 0; i < frame_set->n_mapping_blocks; i++)
        {
            block->id = TNG_PARTICLE_MAPPING;
            if(frame_set->mappings[i].n_particles > 0)
            {
                tng_trajectory_mapping_block_write(tng_data, block, i, hash_mode);
                for(j = 0; j<frame_set->n_particle_data_blocks; j++)
                {
                    block->id = frame_set->tr_particle_data[j].block_id;
                    tng_particle_data_block_write(tng_data, block,
                                                  j, &frame_set->mappings[i],
                                                  hash_mode);
                }
            }
        }
    }
    else
    {
        for(i = 0; i<frame_set->n_particle_data_blocks; i++)
        {
            block->id = frame_set->tr_particle_data[i].block_id;
            tng_particle_data_block_write(tng_data, block,
                                          i, 0, hash_mode);
        }
    }


    /* Update pointers in the general info block */
    stat = tng_header_pointers_update(tng_data, hash_mode);

    if(stat == TNG_SUCCESS)
    {
        stat = tng_frame_set_pointers_update(tng_data, hash_mode);
    }

    tng_block_destroy(&block);

    return(stat);
}

tng_function_status tng_frame_set_new(tng_trajectory_t tng_data,
                                      const int64_t first_frame,
                                      const int64_t n_frames)
{
    int i;
    tng_gen_block_t block;
    tng_trajectory_frame_set_t frame_set;
    tng_particle_mapping_t mapping;
    FILE *temp = tng_data->input_file;
    int64_t curr_pos;

    frame_set = &tng_data->current_trajectory_frame_set;

    /* Set pointer to previous frame set to the one that was loaded
     * before.
     * FIXME: This is a bit risky. If they are not added in order
     * it will be wrong. */
    if(tng_data->n_trajectory_frame_sets)
    {
        frame_set->prev_frame_set_file_pos =
        tng_data->current_trajectory_frame_set_output_file_pos;
    }

    tng_data->current_trajectory_frame_set_output_file_pos =
    ftell(tng_data->output_file);

    /* Clear mappings if they remain. */
    if(frame_set->n_mapping_blocks && frame_set->mappings)
    {
        for(i = frame_set->n_mapping_blocks; i--;)
        {
            mapping = &frame_set->mappings[i];
            if(mapping->real_particle_numbers)
            {
                free(mapping->real_particle_numbers);
                mapping->real_particle_numbers = 0;
            }
        }
        free(frame_set->mappings);
        frame_set->mappings = 0;
        frame_set->n_mapping_blocks = 0;
    }

    tng_data->n_trajectory_frame_sets++;

    /* Set the medium range pointers */
    if(tng_data->n_trajectory_frame_sets == tng_data->medium_stride_length + 1)
    {
        frame_set->medium_stride_prev_frame_set_file_pos =
        tng_data->first_trajectory_frame_set_output_file_pos;
    }
    else if(tng_data->n_trajectory_frame_sets > tng_data->medium_stride_length + 1)
    {
        /* FIXME: Currently only working if the previous frame set has its
         * medium stride pointer already set. This might need some fixing. */
        if(frame_set->medium_stride_prev_frame_set_file_pos != -1 &&
           frame_set->medium_stride_prev_frame_set_file_pos != 0)
        {
            tng_block_init(&block);
            tng_data->input_file = tng_data->output_file;

            curr_pos = ftell(tng_data->output_file);
            fseek(tng_data->output_file,
                  frame_set->medium_stride_prev_frame_set_file_pos,
                  SEEK_SET);

            if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
            {
                printf("Cannot read frame header. %s: %d\n",
                    __FILE__, __LINE__);
                tng_data->input_file = temp;
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            /* Read the next frame set from the previous frame set and one
             * medium stride step back */
            fseek(tng_data->output_file, block->block_contents_size - 6 *
                sizeof(int64_t), SEEK_CUR);
            if(fread(&frame_set->medium_stride_prev_frame_set_file_pos,
               sizeof(frame_set->medium_stride_prev_frame_set_file_pos),
               1, tng_data->output_file) == 0)
            {
                printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
                tng_data->input_file = temp;
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            if(tng_data->input_endianness_swap_func_64)
            {
                if(tng_data->input_endianness_swap_func_64(tng_data,
                   &frame_set->medium_stride_prev_frame_set_file_pos)
                    != TNG_SUCCESS)
                {
                    printf("Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }

            tng_block_destroy(&block);

            /* Set the long range pointers */
            if(tng_data->n_trajectory_frame_sets == tng_data->long_stride_length + 1)
            {
                frame_set->long_stride_prev_frame_set_file_pos =
                tng_data->first_trajectory_frame_set_output_file_pos;
            }
            else if(tng_data->n_trajectory_frame_sets > tng_data->medium_stride_length + 1)
            {
                /* FIXME: Currently only working if the previous frame set has its
                * long stride pointer already set. This might need some fixing. */
                if(frame_set->long_stride_prev_frame_set_file_pos != -1 &&
                frame_set->long_stride_prev_frame_set_file_pos != 0)
                {
                    tng_block_init(&block);
                    tng_data->input_file = tng_data->output_file;

                    fseek(tng_data->output_file,
                        frame_set->long_stride_prev_frame_set_file_pos,
                        SEEK_SET);

                    if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
                    {
                        printf("Cannot read frame header. %s: %d\n",
                            __FILE__, __LINE__);
                        tng_data->input_file = temp;
                        tng_block_destroy(&block);
                        return(TNG_CRITICAL);
                    }

                    /* Read the next frame set from the previous frame set and one
                    * long stride step back */
                    fseek(tng_data->output_file, block->block_contents_size - 6 *
                        sizeof(int64_t), SEEK_CUR);

                    tng_block_destroy(&block);

                    if(fread(&frame_set->long_stride_prev_frame_set_file_pos,
                    sizeof(frame_set->long_stride_prev_frame_set_file_pos),
                    1, tng_data->output_file) == 0)
                    {
                        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
                        tng_data->input_file = temp;
                        return(TNG_CRITICAL);
                    }

                    if(tng_data->input_endianness_swap_func_64)
                    {
                        if(tng_data->input_endianness_swap_func_64(tng_data,
                           &frame_set->long_stride_prev_frame_set_file_pos)
                            != TNG_SUCCESS)
                        {
                            printf("Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }

                }
            }

            tng_data->input_file = temp;
            fseek(tng_data->output_file, curr_pos, SEEK_SET);
        }
    }

    frame_set->first_frame = first_frame;
    frame_set->n_frames = n_frames;
    frame_set->n_written_frames = 0;

    if(tng_data->first_trajectory_frame_set_output_file_pos == -1 ||
       tng_data->first_trajectory_frame_set_output_file_pos == 0)
    {
        tng_data->first_trajectory_frame_set_output_file_pos =
        tng_data->current_trajectory_frame_set_output_file_pos;
    }
    /* FIXME: Should check the frame number instead of the file_pos,
     * in case frame sets are not in order */
    if(tng_data->last_trajectory_frame_set_output_file_pos == -1 ||
       tng_data->last_trajectory_frame_set_output_file_pos == 0 ||
       tng_data->last_trajectory_frame_set_output_file_pos <
       tng_data->current_trajectory_frame_set_output_file_pos)
    {
        tng_data->last_trajectory_frame_set_output_file_pos =
        tng_data->current_trajectory_frame_set_output_file_pos;
    }

    return(TNG_SUCCESS);
}



tng_function_status tng_data_block_add(tng_trajectory_t tng_data,
                                        const int64_t id,
                                        const char *block_name,
                                        const tng_data_type datatype,
                                        const tng_block_type block_type_flag,
                                        int64_t n_frames,
                                        const int64_t n_values_per_frame,
                                        int64_t stride_length,
                                        const int64_t codec_id,
                                        void *new_data)
{
    int i, j, block_index, size, len;
    tng_trajectory_frame_set_t frame_set;
    tng_non_particle_data_t data;
    union data_values *first_dim_values;
    void *orig;

    frame_set = &tng_data->current_trajectory_frame_set;

    if(stride_length <= 0)
    {
        stride_length = 1;
    }

    block_index = -1;
    /* See if there is already a data block of this ID */
    if(block_type_flag == TNG_TRAJECTORY_BLOCK)
    {
        for(i = frame_set->n_data_blocks; i-- ;)
        {
            data = &frame_set->tr_data[i];
            if(data->block_id == id)
            {
                block_index = i;
                break;
            }
        }
    }
    else
    {
        n_frames = 1;
        for(i = tng_data->n_data_blocks; i-- ;)
        {
            data = &tng_data->non_tr_data[i];
            if(data->block_id == id)
            {
                block_index = i;
                break;
            }
        }
    }

    /* Otherwise create a data block */
    if(block_index == -1)
    {
        if(tng_data_block_create(tng_data, block_type_flag) !=
            TNG_SUCCESS)
        {
            printf("Cannot create data block. %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        if(block_type_flag == TNG_TRAJECTORY_BLOCK)
        {
            data = &frame_set->tr_data[frame_set->n_data_blocks - 1];
        }
        else
        {
            data = &tng_data->non_tr_data[tng_data->n_data_blocks - 1];
        }
        data->block_id = id;

        data->block_name = malloc(strlen(block_name) + 1);
        if(!data->block_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                   (int)strlen(block_name)+1, __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        strncpy(data->block_name, block_name, strlen(block_name) + 1);

        data->datatype = datatype;

        data->stride_length = tng_max(stride_length, 1);
        data->values = 0;
        data->n_values_per_frame = n_values_per_frame;
        data->n_frames = n_frames;
        data->codec_id = codec_id;
        data->compression_multiplier = 1.0;
    }

    switch(datatype)
    {
    case TNG_FLOAT_DATA:
        size = sizeof(float);
        break;
    case TNG_INT_DATA:
        size = sizeof(int64_t);
        break;
    case TNG_CHAR_DATA:
        size = sizeof(char);
    case TNG_DOUBLE_DATA:
    default:
        size = sizeof(double);
        break;
    }

    if(new_data)
    {
        /* Allocate memory */
        if(tng_allocate_data_mem(tng_data, data, n_frames, stride_length,
                                 n_values_per_frame) !=
        TNG_SUCCESS)
        {
            printf("Cannot allocate data memory. %s: %d\n",
                __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }

        orig = new_data;

        if(n_frames > frame_set->n_written_frames)
        {
            frame_set->n_written_frames = n_frames;
        }

        switch(datatype)
        {
        case TNG_FLOAT_DATA:
            for(i = 0; i < n_frames / stride_length; i++)
            {
                first_dim_values = data->values[i];
                for(j = 0; j < n_values_per_frame; j++)
                {
                    memcpy(&first_dim_values[j].f,
                            new_data, size);
                    new_data += size;
                }
            }
            break;
        case TNG_CHAR_DATA:
            for(i = 0; i < n_frames / stride_length; i++)
            {
                first_dim_values = data->values[i];
                for(j = 0; j < n_values_per_frame; j++)
                {
                    len = tng_min(strlen(new_data) + 1,
                                TNG_MAX_STR_LEN);
                    if(first_dim_values[j].c)
                    {
                        free(first_dim_values[j].c);
                    }
                    first_dim_values[j].c = malloc(len);
                    if(!first_dim_values[j].c)
                    {
                        printf("Cannot allocate memory (%d bytes). %s: %d\n",
                            len, __FILE__, __LINE__);
                        return(TNG_CRITICAL);
                    }
                    strncpy(first_dim_values[j].c,
                            new_data, len);
                    new_data += len;
                }
            }
            break;
        case TNG_INT_DATA:
            for(i = 0; i < n_frames / stride_length; i++)
            {
                first_dim_values = data->values[i];
                for(j = 0; j < n_values_per_frame; j++)
                {
                    memcpy(&first_dim_values[j].i,
                            new_data, size);
                    new_data += size;
                }
            }
            break;
        case TNG_DOUBLE_DATA:
        default:
            for(i = 0; i < n_frames / stride_length; i++)
            {
                first_dim_values = data->values[i];
                for(j = 0; j < n_values_per_frame; j++)
                {
                    memcpy(&first_dim_values[j].d,
                            new_data, size);
                    new_data += size;
                }
            }
        }

        new_data = orig;
    }
//     else
//     {
//         for(i = 0; i < n_frames / stride_length; i++)
//         {
//             first_dim_values = data->values[i];
//             for(j = 0; j < n_values_per_frame; j++)
//             {
//                 first_dim_values[j].d = 0;
//             }
//         }
//     }

    return(TNG_SUCCESS);
}

tng_function_status tng_particle_data_block_add(tng_trajectory_t tng_data,
                                        const int64_t id,
                                        const char *block_name,
                                        const tng_data_type datatype,
                                        const tng_block_type block_type_flag,
                                        int64_t n_frames,
                                        const int64_t n_values_per_frame,
                                        int64_t stride_length,
                                        const int64_t num_first_particle,
                                        const int64_t n_particles,
                                        const int64_t codec_id,
                                        void *new_data)
{
    int i, j, k, block_index, size, len;
    int64_t tot_n_particles;
    union data_values **first_dim_values, *second_dim_values;
    tng_trajectory_frame_set_t frame_set;
    tng_particle_data_t data;
    void *orig;

    frame_set = &tng_data->current_trajectory_frame_set;

    if(stride_length <= 0)
    {
        stride_length = 1;
    }

    block_index = -1;
    /* See if there is already a data block of this ID */
    if(block_type_flag == TNG_TRAJECTORY_BLOCK)
    {
        for(i = frame_set->n_particle_data_blocks; i-- ;)
        {
            data = &frame_set->tr_particle_data[i];
            if(data->block_id == id)
            {
                block_index = i;
                break;
            }
        }
    }
    else
    {
        n_frames = 1;
        for(i = tng_data->n_particle_data_blocks; i-- ;)
        {
            data = &tng_data->non_tr_particle_data[i];
            if(data->block_id == id)
            {
                block_index = i;
                break;
            }
        }
    }

    /* Otherwise create a data block */
    if(block_index == -1)
    {
        if(tng_particle_data_block_create(tng_data, block_type_flag) !=
            TNG_SUCCESS)
        {
            printf("Cannot create particle data block. %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        if(block_type_flag == TNG_TRAJECTORY_BLOCK)
        {
            data = &frame_set->tr_particle_data[frame_set->
                                                n_particle_data_blocks - 1];
        }
        else
        {
            data = &tng_data->non_tr_particle_data[tng_data->
                                                   n_particle_data_blocks - 1];
        }
        data->block_id = id;

        data->block_name = malloc(strlen(block_name) + 1);
        if(!data->block_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                   (int)strlen(block_name)+1, __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        strncpy(data->block_name, block_name, strlen(block_name) + 1);

        data->datatype = datatype;

        data->stride_length = tng_max(stride_length, 1);
        data->values = 0;
        data->n_values_per_frame = n_values_per_frame;
        data->n_frames = n_frames;
        data->codec_id = codec_id;
        data->compression_multiplier = 1.0;
    }

    if(block_type_flag == TNG_TRAJECTORY_BLOCK && tng_data->var_num_atoms_flag)
    {
        tot_n_particles = frame_set->n_particles;
    }
    else
    {
        tot_n_particles = tng_data->n_particles;
    }

    /* If data values are supplied add that data to the data block. */
    if(new_data)
    {
        /* Allocate memory */
        if(tng_allocate_particle_data_mem(tng_data, data, n_frames,
                                          stride_length, tot_n_particles,
                                          n_values_per_frame) !=
        TNG_SUCCESS)
        {
            printf("Cannot allocate particle data memory. %s: %d\n",
                __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }

        orig = new_data;

        if(n_frames > frame_set->n_written_frames)
        {
            frame_set->n_written_frames = n_frames;
        }

        switch(datatype)
        {
        case TNG_FLOAT_DATA:
            size = sizeof(float);
            for(i = 0; i < n_frames / stride_length; i++)
            {
                first_dim_values = data->values[i];
                for(j = num_first_particle; j < num_first_particle + n_particles;
                    j++)
                {
                    second_dim_values = first_dim_values[j];
                    for(k = 0; k < n_values_per_frame; k++)
                    {
                        memcpy(&second_dim_values[k].f,
                                new_data, size);
                        new_data += size;
                    }
                }
            }
            break;
        case TNG_INT_DATA:
            size = sizeof(int64_t);
            for(i = 0; i < n_frames / stride_length; i++)
            {
                first_dim_values = data->values[i];
                for(j = num_first_particle; j < num_first_particle + n_particles;
                    j++)
                {
                    second_dim_values = first_dim_values[j];
                    for(k = 0; k < n_values_per_frame; k++)
                    {
                        memcpy(&second_dim_values[k].i,
                                new_data, size);
                        new_data += size;
                    }
                }
            }
            break;
        case TNG_CHAR_DATA:
            for(i = 0; i < n_frames / stride_length; i++)
            {
                first_dim_values = data->values[i];
                for(j = num_first_particle; j < num_first_particle + n_particles;
                    j++)
                {
                    second_dim_values = first_dim_values[j];
                    for(k = 0; k < n_values_per_frame; k++)
                    {
                        len = tng_min(strlen(new_data) + 1,
                                TNG_MAX_STR_LEN);
                        if(second_dim_values[k].c)
                        {
                            free(second_dim_values[k].c);
                        }
                        second_dim_values[k].c = malloc(len);
                        if(!second_dim_values[k].c)
                        {
                            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                                len, __FILE__, __LINE__);
                            return(TNG_CRITICAL);
                        }
                        strncpy(second_dim_values[k].c,
                                new_data, len);
                        new_data += len;
                    }
                }
            }
            break;
        case TNG_DOUBLE_DATA:
        default:
            size = sizeof(double);
            for(i = 0; i < n_frames / stride_length; i++)
            {
                first_dim_values = data->values[i];
                for(j = num_first_particle; j < num_first_particle + n_particles;
                    j++)
                {
                    second_dim_values = first_dim_values[j];
                    for(k = 0; k < n_values_per_frame; k++)
                    {
                        memcpy(&second_dim_values[k].d,
                                new_data, size);
                        new_data += size;
                    }
                }
            }
            break;
        }

        new_data = orig;
    }
    /* Otherwise fill the data block with zeroes */
//     else
//     {
//         for(i = 0; i < n_frames / stride_length; i++)
//         {
//             first_dim_values = data->values[i];
//             for(j = num_first_particle; j < num_first_particle + n_particles;
//                 j++)
//             {
//                 second_dim_values = first_dim_values[j];
//                 for(k = 0; k < n_values_per_frame; k++)
//                 {
//                     second_dim_values[k].d = 0;
//                 }
//             }
//         }
//     }


    return(TNG_SUCCESS);
}

tng_function_status tng_frame_data_write(tng_trajectory_t tng_data,
                                        const int64_t frame_nr,
                                        const int64_t block_id,
                                        const void *values,
                                        const tng_hash_mode hash_mode)
{
    int64_t header_pos, file_pos;
    int64_t output_file_len, n_values_per_frame, size, contents_size;
    int64_t header_size, temp_first, temp_last, temp_current;
    int64_t i, last_frame;
    tng_gen_block_t block;
    tng_trajectory_frame_set_t frame_set;
    FILE *temp = tng_data->input_file;
    struct tng_non_particle_data data;
    tng_function_status stat;
    char dependency, sparse_data, datatype;
    void *copy;

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        printf("Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    temp_first = tng_data->first_trajectory_frame_set_input_file_pos;
    temp_last = tng_data->last_trajectory_frame_set_input_file_pos;
    temp_current = tng_data->current_trajectory_frame_set_input_file_pos;
    tng_data->first_trajectory_frame_set_input_file_pos =
    tng_data->first_trajectory_frame_set_output_file_pos;
    tng_data->last_trajectory_frame_set_input_file_pos =
    tng_data->last_trajectory_frame_set_output_file_pos;
    tng_data->current_trajectory_frame_set_input_file_pos =
    tng_data->current_trajectory_frame_set_output_file_pos;

    tng_data->input_file = tng_data->output_file;

    stat = tng_frame_set_of_frame_find(tng_data, frame_nr);

    frame_set = &tng_data->current_trajectory_frame_set;

    if(stat != TNG_SUCCESS)
    {
        last_frame = frame_set->first_frame +
                     frame_set->n_frames - 1;
        /* If the wanted frame would be in the frame set after the last
            * frame set create a new frame set. */
        if(stat == TNG_FAILURE &&
           (last_frame < frame_nr &&
            tng_data->current_trajectory_frame_set.first_frame +
            tng_data->frame_set_n_frames >= frame_nr))
        {
            tng_frame_set_new(tng_data,
                              last_frame+1,
                              tng_data->frame_set_n_frames);
            file_pos = ftell(tng_data->output_file);
            fseek(tng_data->output_file, 0, SEEK_END);
            output_file_len = ftell(tng_data->output_file);
            fseek(tng_data->output_file, file_pos, SEEK_SET);

            /* Read mapping blocks from the last frame set */
            tng_block_init(&block);

            stat = tng_block_header_read(tng_data, block);
            while(file_pos < output_file_len &&
                  stat != TNG_CRITICAL &&
                  block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                if(block->id == TNG_PARTICLE_MAPPING)
                {
                    tng_trajectory_mapping_block_read(tng_data, block,
                                                      hash_mode);
                }
                else
                {
                    fseek(tng_data->output_file, block->block_contents_size,
                        SEEK_CUR);
                }
                file_pos = ftell(tng_data->output_file);
                if(file_pos < output_file_len)
                {
                    stat = tng_block_header_read(tng_data, block);
                }
            }

            tng_block_destroy(&block);
            /* Write the frame set to disk */
            if(tng_frame_set_write(tng_data, hash_mode) != TNG_SUCCESS)
            {
                printf("Error writing frame set. %s: %d\n", __FILE__, __LINE__);
                return(TNG_CRITICAL);
            }
        }
        else
        {
            tng_data->input_file = temp;
            tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
            tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
            tng_data->current_trajectory_frame_set_input_file_pos = temp_current;
            return(stat);
        }
    }

    tng_block_init(&block);

    file_pos = tng_data->current_trajectory_frame_set_output_file_pos;

    fseek(tng_data->output_file, 0, SEEK_END);
    output_file_len = ftell(tng_data->output_file);
    fseek(tng_data->output_file, file_pos, SEEK_SET);

    /* Read past the frame set block first */
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL)
    {
        printf("Cannot read block header at pos %"PRId64". %s: %d\n",
               file_pos, __FILE__, __LINE__);
        tng_block_destroy(&block);
        tng_data->input_file = temp;

        tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
        tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
        tng_data->current_trajectory_frame_set_input_file_pos = temp_current;
        return(stat);
    }
    fseek(tng_data->output_file, block->block_contents_size,
            SEEK_CUR);

    /* Read all block headers until next frame set block or
     * until the wanted block id is found */
    stat = tng_block_header_read(tng_data, block);
    while(file_pos < output_file_len &&
            stat != TNG_CRITICAL &&
            block->id != block_id &&
            block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        fseek(tng_data->output_file, block->block_contents_size, SEEK_CUR);
        file_pos = ftell(tng_data->output_file);
        if(file_pos < output_file_len)
        {
            stat = tng_block_header_read(tng_data, block);
        }
    }
    if(stat == TNG_CRITICAL)
    {
        printf("Cannot read block header at pos %"PRId64". %s: %d\n",
               file_pos, __FILE__, __LINE__);
        tng_block_destroy(&block);
        tng_data->input_file = temp;
        tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
        tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
        tng_data->current_trajectory_frame_set_input_file_pos = temp_current;
        return(stat);
    }

    contents_size = block->block_contents_size;
    header_size = block->header_contents_size;

    header_pos = ftell(tng_data->output_file) - header_size;
    frame_set = &tng_data->current_trajectory_frame_set;

    if(fread(&datatype, sizeof(datatype), 1, tng_data->input_file) == 0)
    {
        printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(fread(&dependency, sizeof(dependency), 1, tng_data->input_file) == 0)
    {
        printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    data.datatype = datatype;

    if(!(dependency & TNG_FRAME_DEPENDENT) ||
       (dependency & TNG_PARTICLE_DEPENDENT))
    {
        tng_block_destroy(&block);
        tng_data->input_file = temp;

        tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
        tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
        tng_data->current_trajectory_frame_set_input_file_pos = temp_current;
        return(TNG_FAILURE);
    }

    if(fread(&sparse_data, sizeof(sparse_data), 1, tng_data->input_file) == 0)
    {
        printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    if(fread(&data.n_values_per_frame, sizeof(data.n_values_per_frame), 1,
             tng_data->input_file) == 0)
    {
        printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
            &data.n_values_per_frame)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(fread(&data.codec_id, sizeof(data.codec_id), 1,
             tng_data->input_file) == 0)
    {
        printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
            &data.codec_id)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(data.codec_id != TNG_UNCOMPRESSED)
    {
        if(fread(&data.compression_multiplier,
                 sizeof(data.compression_multiplier), 1, tng_data->input_file)
            == 0)
        {
            printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                (int64_t *)&data.compression_multiplier)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
    }
    else
    {
        data.compression_multiplier = 1;
    }

    if(sparse_data)
    {
        if(fread(&data.first_frame_with_data, sizeof(data.first_frame_with_data),
                 1, tng_data->input_file) == 0)
        {
            printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                &data.first_frame_with_data)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }

        if(fread(&data.stride_length, sizeof(data.stride_length),
                 1, tng_data->input_file) == 0)
        {
            printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                &data.stride_length)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
    }
    else
    {
        data.first_frame_with_data = 0;
        data.stride_length = 1;
    }
    data.n_frames = tng_data->current_trajectory_frame_set.n_frames;

    tng_data->input_file = temp;

    tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
    tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
    tng_data->current_trajectory_frame_set_input_file_pos = temp_current;

    switch(data.datatype)
    {
        case(TNG_INT_DATA):
            size = sizeof(int64_t);
            break;
        case(TNG_FLOAT_DATA):
            size = sizeof(float);
            break;
        case(TNG_DOUBLE_DATA):
            size = sizeof(double);
            break;
        default:
            printf("Cannot calculate writing locations. %s: %d.\n", __FILE__,
                   __LINE__);
            tng_block_destroy(&block);
            return(TNG_FAILURE);
    }

    n_values_per_frame = data.n_values_per_frame;

    file_pos = (frame_nr - tng_max(frame_set->first_frame,
                               data.first_frame_with_data)) /
                data.stride_length;
    file_pos *= size * n_values_per_frame;

    if(file_pos > contents_size)
    {
        printf("Attempting to write outside the block. %s: %d\n", __FILE__,
               __LINE__);
        tng_block_destroy(&block);
        return(TNG_FAILURE);
    }

    fseek(tng_data->output_file, file_pos, SEEK_CUR);

    /* If the endianness is not big endian the data needs to be swapped */
    if((data.datatype == TNG_INT_DATA ||
        data.datatype == TNG_DOUBLE_DATA) &&
       tng_data->output_endianness_swap_func_64)
    {
        copy = malloc(n_values_per_frame * size);
        memcpy(copy, values, n_values_per_frame * size);
        for(i = 0; i < n_values_per_frame; i++)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                (int64_t *)copy+i)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        fwrite(copy, n_values_per_frame, size,
               tng_data->output_file);
        free(copy);
    }
    else if(data.datatype == TNG_FLOAT_DATA &&
            tng_data->output_endianness_swap_func_32)
    {
        copy = malloc(n_values_per_frame * size);
        memcpy(copy, values, n_values_per_frame * size);
        for(i = 0; i < n_values_per_frame; i++)
        {
            if(tng_data->output_endianness_swap_func_32(tng_data,
                (int32_t *)copy+i)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        fwrite(copy, n_values_per_frame, size,
               tng_data->output_file);
        free(copy);
    }

    else
    {
        fwrite(values, n_values_per_frame, size, tng_data->output_file);
    }

    fflush(tng_data->output_file);

    /* Update the number of written frames in the frame set. */
    if(frame_nr - frame_set->first_frame + 1 > frame_set->n_written_frames)
    {
        frame_set->n_written_frames = frame_nr - frame_set->first_frame + 1;
    }

    /* If the last frame has been written update the hash */
    if(hash_mode == TNG_USE_HASH && (frame_nr + data.stride_length -
       data.first_frame_with_data) >=
       frame_set->n_frames)
    {
        tng_md5_hash_update(tng_data, block, header_pos, header_pos +
                            header_size);
    }

    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

tng_function_status tng_frame_particle_data_write(tng_trajectory_t tng_data,
                                                  const int64_t frame_nr,
                                                  const int64_t block_id,
                                                  const int64_t val_first_particle,
                                                  const int64_t val_n_particles,
                                                  const void *values,
                                                  const tng_hash_mode hash_mode)
{
    int64_t header_pos, file_pos, tot_n_particles;
    int64_t output_file_len, n_values_per_frame, size, contents_size;
    int64_t header_size, temp_first, temp_last, temp_current;
    int64_t mapping_block_end_pos, num_first_particle, block_n_particles;
    int64_t i, last_frame;
    tng_gen_block_t block;
    tng_trajectory_frame_set_t frame_set;
    FILE *temp = tng_data->input_file;
    struct tng_particle_data data;
    tng_function_status stat;
    tng_particle_mapping_t mapping;
    char dependency, sparse_data, datatype;
    void *copy;

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        printf("Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    temp_first = tng_data->first_trajectory_frame_set_input_file_pos;
    temp_last = tng_data->last_trajectory_frame_set_input_file_pos;
    temp_current = tng_data->current_trajectory_frame_set_input_file_pos;
    tng_data->first_trajectory_frame_set_input_file_pos =
    tng_data->first_trajectory_frame_set_output_file_pos;
    tng_data->last_trajectory_frame_set_input_file_pos =
    tng_data->last_trajectory_frame_set_output_file_pos;
    tng_data->current_trajectory_frame_set_input_file_pos =
    tng_data->current_trajectory_frame_set_output_file_pos;

    tng_data->input_file = tng_data->output_file;

    stat = tng_frame_set_of_frame_find(tng_data, frame_nr);

    frame_set = &tng_data->current_trajectory_frame_set;

    if(stat != TNG_SUCCESS)
    {
        last_frame = frame_set->first_frame +
                     frame_set->n_frames - 1;
//         printf("Frame %"PRId64" not found. Last frame: %"PRId64"\n", frame_nr,
//                last_frame);
        /* If the wanted frame would be in the frame set after the last
         * frame set create a new frame set. */
        if(stat == TNG_FAILURE &&
           (last_frame < frame_nr &&
            last_frame + tng_data->frame_set_n_frames >= frame_nr))
        {
            tng_frame_set_new(tng_data,
                              last_frame+1,
                              tng_data->frame_set_n_frames);

            file_pos = ftell(tng_data->output_file);
            fseek(tng_data->output_file, 0, SEEK_END);
            output_file_len = ftell(tng_data->output_file);
            fseek(tng_data->output_file, file_pos, SEEK_SET);

            /* Read mapping blocks from the last frame set */
            tng_block_init(&block);

            stat = tng_block_header_read(tng_data, block);
            while(file_pos < output_file_len &&
                  stat != TNG_CRITICAL &&
                  block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                if(block->id == TNG_PARTICLE_MAPPING)
                {
                    tng_trajectory_mapping_block_read(tng_data, block,
                                                      hash_mode);
                }
                else
                {
                    fseek(tng_data->output_file, block->block_contents_size,
                        SEEK_CUR);
                }
                file_pos = ftell(tng_data->output_file);
                if(file_pos < output_file_len)
                {
                    stat = tng_block_header_read(tng_data, block);
                }
            }

            tng_block_destroy(&block);
            /* Write the frame set to disk */
            if(tng_frame_set_write(tng_data, hash_mode) != TNG_SUCCESS)
            {
                printf("Error writing frame set. %s: %d\n", __FILE__, __LINE__);
                exit(1);
            }
        }
        else
        {
            tng_data->input_file = temp;
            tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
            tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
            tng_data->current_trajectory_frame_set_input_file_pos = temp_current;
            return(stat);
        }
    }


    tng_block_init(&block);

    file_pos = tng_data->current_trajectory_frame_set_output_file_pos;

    fseek(tng_data->output_file, 0, SEEK_END);
    output_file_len = ftell(tng_data->output_file);
    fseek(tng_data->output_file, file_pos, SEEK_SET);

    /* Read past the frame set block first */
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL)
    {
        printf("Cannot read block header at pos %"PRId64". %s: %d\n",
               file_pos, __FILE__, __LINE__);
        tng_block_destroy(&block);
        tng_data->input_file = temp;

        tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
        tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
        tng_data->current_trajectory_frame_set_input_file_pos = temp_current;
        return(stat);
    }
    fseek(tng_data->output_file, block->block_contents_size,
            SEEK_CUR);

    if(tng_data->var_num_atoms_flag)
    {
        tot_n_particles = frame_set->n_particles;
    }
    else
    {
        tot_n_particles = tng_data->n_particles;
    }

    if(val_n_particles < tot_n_particles)
    {
        mapping_block_end_pos = -1;
        /* Read all mapping blocks to find the right place to put the data */
        stat = tng_block_header_read(tng_data, block);
        while(file_pos < output_file_len &&
                stat != TNG_CRITICAL &&
                block->id != TNG_TRAJECTORY_FRAME_SET)
        {
            if(block->id == TNG_PARTICLE_MAPPING)
            {
                tng_trajectory_mapping_block_read(tng_data, block, hash_mode);
            }
            else
            {
                fseek(tng_data->output_file, block->block_contents_size,
                      SEEK_CUR);
            }
            file_pos = ftell(tng_data->output_file);
            if(block->id == TNG_PARTICLE_MAPPING)
            {
                mapping = &frame_set->mappings[frame_set->n_mapping_blocks - 1];
                if(val_first_particle >= mapping->num_first_particle &&
                   val_first_particle < mapping->num_first_particle +
                   mapping->n_particles &&
                   val_first_particle + val_n_particles <=
                   mapping->num_first_particle + mapping->n_particles)
                {
                    mapping_block_end_pos = file_pos;
                }
            }
            if(file_pos < output_file_len)
            {
                stat = tng_block_header_read(tng_data, block);
            }
        }
        if(stat == TNG_CRITICAL)
        {
            printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                   file_pos, __FILE__, __LINE__);
            tng_block_destroy(&block);
            tng_data->input_file = temp;

            tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
            tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
            tng_data->current_trajectory_frame_set_input_file_pos = temp_current;
            return(stat);
        }
        if(mapping_block_end_pos < 0)
        {
            tng_block_destroy(&block);
            tng_data->input_file = temp;

            tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
            tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
            tng_data->current_trajectory_frame_set_input_file_pos = temp_current;
            return(TNG_FAILURE);
        }
        fseek(tng_data->output_file, mapping_block_end_pos, SEEK_SET);
    }

    /* Read all block headers until next frame set block or
     * until the wanted block id is found */
    stat = tng_block_header_read(tng_data, block);
    while(file_pos < output_file_len &&
            stat != TNG_CRITICAL &&
            block->id != block_id &&
            block->id != TNG_PARTICLE_MAPPING &&
            block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        fseek(tng_data->output_file, block->block_contents_size, SEEK_CUR);
        file_pos = ftell(tng_data->output_file);
        if(file_pos < output_file_len)
        {
            stat = tng_block_header_read(tng_data, block);
        }
    }
    if(stat == TNG_CRITICAL)
    {
        printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                file_pos, __FILE__, __LINE__);
        tng_block_destroy(&block);
        tng_data->input_file = temp;

        tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
        tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
        tng_data->current_trajectory_frame_set_input_file_pos = temp_current;
        return(stat);
    }

    contents_size = block->block_contents_size;
    header_size = block->header_contents_size;

    header_pos = ftell(tng_data->output_file) - header_size;
    frame_set = &tng_data->current_trajectory_frame_set;

    if(fread(&datatype, sizeof(datatype), 1, tng_data->input_file) == 0)
    {
        printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    data.datatype = datatype;

    if(fread(&dependency, sizeof(dependency), 1, tng_data->input_file) == 0)
    {
        printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    if(!(dependency & TNG_FRAME_DEPENDENT) ||
       !(dependency & TNG_PARTICLE_DEPENDENT))
    {
        tng_block_destroy(&block);
        tng_data->input_file = temp;

        tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
        tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
        tng_data->current_trajectory_frame_set_input_file_pos = temp_current;
        return(TNG_FAILURE);
    }

    if(fread(&sparse_data, sizeof(sparse_data), 1, tng_data->input_file) == 0)
    {
        printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    if(fread(&data.n_values_per_frame, sizeof(data.n_values_per_frame), 1,
             tng_data->input_file) == 0)
    {
        printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
            &data.n_values_per_frame)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(fread(&data.codec_id, sizeof(data.codec_id), 1,
             tng_data->input_file) == 0)
    {
        printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
            &data.codec_id)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(data.codec_id != TNG_UNCOMPRESSED)
    {
        if(fread(&data.compression_multiplier,
                 sizeof(data.compression_multiplier), 1, tng_data->input_file)
            == 0)
        {
            printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
               (int64_t *)&data.compression_multiplier)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
    }
    else
    {
        data.compression_multiplier = 1;
    }

    if(sparse_data)
    {
        if(fread(&data.first_frame_with_data,
                 sizeof(data.first_frame_with_data),
                 1, tng_data->input_file) == 0)
        {
            printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                &data.first_frame_with_data)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }

        if(fread(&data.stride_length, sizeof(data.stride_length),
                 1, tng_data->input_file) == 0)
        {
            printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                &data.stride_length)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
    }
    else
    {
        data.first_frame_with_data = 0;
        data.stride_length = 1;
    }
    data.n_frames = tng_data->current_trajectory_frame_set.n_frames;

    if(fread(&num_first_particle, sizeof(num_first_particle), 1,
             tng_data->input_file) == 0)
    {
        printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
            &num_first_particle)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(fread(&block_n_particles, sizeof(block_n_particles), 1,
             tng_data->input_file) == 0)
    {
        printf("Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
            &block_n_particles)
            != TNG_SUCCESS)
        {
            printf("Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }


    tng_data->input_file = temp;

    tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
    tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
    tng_data->current_trajectory_frame_set_input_file_pos = temp_current;


    switch(data.datatype)
    {
        case(TNG_INT_DATA):
            size = sizeof(int64_t);
            break;
        case(TNG_FLOAT_DATA):
            size = sizeof(float);
            break;
        case(TNG_DOUBLE_DATA):
            size = sizeof(double);
            break;
        default:
            printf("Cannot calculate writing locations. %s: %d.\n", __FILE__,
                   __LINE__);
            tng_block_destroy(&block);
            return(TNG_FAILURE);
    }

    n_values_per_frame = data.n_values_per_frame;

    file_pos = (frame_nr - tng_max(frame_set->first_frame,
                               data.first_frame_with_data)) /
                data.stride_length;
    file_pos *= block_n_particles * size * n_values_per_frame;

    if(file_pos > contents_size)
    {
        printf("Attempting to write outside the block. %s: %d\n", __FILE__,
               __LINE__);
        tng_block_destroy(&block);
        return(TNG_FAILURE);
    }

    fseek(tng_data->output_file, file_pos, SEEK_CUR);

    /* If the endianness is not big endian the data needs to be swapped */
    if((data.datatype == TNG_INT_DATA ||
        data.datatype == TNG_DOUBLE_DATA) &&
       tng_data->output_endianness_swap_func_64)
    {
        copy = malloc(val_n_particles * n_values_per_frame * size);
        memcpy(copy, values, val_n_particles * n_values_per_frame * size);
        for(i = 0; i < val_n_particles * n_values_per_frame; i++)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                (int64_t *) copy+i)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        fwrite(copy, val_n_particles * n_values_per_frame, size,
               tng_data->output_file);
        free(copy);
    }
    else if(data.datatype == TNG_FLOAT_DATA &&
       tng_data->output_endianness_swap_func_32)
    {
        copy = malloc(val_n_particles * n_values_per_frame * size);
        memcpy(copy, values, val_n_particles * n_values_per_frame * size);
        for(i = 0; i < val_n_particles * n_values_per_frame; i++)
        {
            if(tng_data->output_endianness_swap_func_32(tng_data,
                (int32_t *) copy+i)
                != TNG_SUCCESS)
            {
                printf("Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        fwrite(copy, val_n_particles * n_values_per_frame, size,
               tng_data->output_file);
        free(copy);
    }

    else
    {
        fwrite(values, val_n_particles * n_values_per_frame, size,
            tng_data->output_file);
    }
    fflush(tng_data->output_file);

    /* Update the number of written frames in the frame set. */
    if(frame_nr - frame_set->first_frame + 1 > frame_set->n_written_frames)
    {
        frame_set->n_written_frames = frame_nr - frame_set->first_frame + 1;
    }

    /* If the last frame has been written update the hash */
    if(hash_mode == TNG_USE_HASH && (frame_nr + data.stride_length -
       data.first_frame_with_data) >=
       frame_set->n_frames)
    {
        tng_md5_hash_update(tng_data, block, header_pos, header_pos +
                            header_size);
    }

    tng_block_destroy(&block);
    return(TNG_SUCCESS);
}

tng_function_status tng_data_values_free(const tng_trajectory_t tng_data,
                                         union data_values **values,
                                         const int64_t n_frames,
                                         const int64_t n_values_per_frame,
                                         const tng_data_type type)
{
    int i, j;

    if(values)
    {
        for(i = 0; i < n_frames; i++)
        {
            if(values[i])
            {
                if(type == TNG_CHAR_DATA)
                {
                    for(j = n_values_per_frame; j--;)
                    {
                        if(values[i][j].c)
                        {
                            free(values[i][j].c);
                            values[i][j].c = 0;
                        }
                    }
                }
                free(values[i]);
                values[i] = 0;
            }
        }
        free(values);
        values = 0;
    }

    return(TNG_SUCCESS);
}


tng_function_status tng_particle_data_values_free
                (const tng_trajectory_t tng_data,
                 union data_values ***values,
                 const int64_t n_frames,
                 const int64_t n_particles,
                 const int64_t n_values_per_frame,
                 const tng_data_type type)
{
    int i, j, k;

    if(values)
    {
        for(i = 0; i < n_frames; i++)
        {
            if(values[i])
            {
                for(j = 0; j < n_particles; j++)
                {
                    if(type == TNG_CHAR_DATA)
                    {
                        for(k = n_values_per_frame; k--;)
                        {
                            if(values[i][j][k].c)
                            {
                                free(values[i][j][k].c);
                                values[i][j][k].c = 0;
                            }
                        }
                    }
                    free(values[i][j]);
                    values[i][j] = 0;
                }
                free(values[i]);
                values[i] = 0;
            }
        }
        free(values);
        values = 0;
    }

    return(TNG_SUCCESS);
}


tng_function_status tng_data_get(tng_trajectory_t tng_data,
                                 const int64_t block_id,
                                 union data_values ***values,
                                 int64_t *n_frames,
                                 int64_t *n_values_per_frame,
                                 tng_data_type *type)
{
    int64_t file_pos;
    int i, j, block_index, len;
    tng_non_particle_data_t data, new_data;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;
    tng_gen_block_t block;
    tng_function_status stat;

    block_index = -1;
    /* See if there is a data block of this ID.
     * Start checking the last read frame set */
    for(i = frame_set->n_data_blocks; i-- ;)
    {
        data = &frame_set->tr_data[i];
        if(data->block_id == block_id)
        {
            block_index = i;
            break;
        }
    }

    if(block_index < 0)
    {
        /* If the data block was not found in the frame set
         * look for it in the non-trajectory data (in tng_data). */
        for(i = tng_data->n_data_blocks; i-- ;)
        {
            data = &tng_data->non_tr_data[i];
            if(data->block_id == block_id)
            {
                block_index = i;
                break;
            }
        }
        if(block_index < 0)
        {
            tng_block_init(&block);
            file_pos = ftell(tng_data->input_file);
            /* Read all blocks until next frame set block */
            stat = tng_block_header_read(tng_data, block);
            while(file_pos < tng_data->input_file_len &&
                    stat != TNG_CRITICAL &&
                    block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                /* Use hash by default */
                stat = tng_block_read_next(tng_data, block,
                                        TNG_USE_HASH);
                if(stat != TNG_CRITICAL)
                {
                    file_pos = ftell(tng_data->input_file);
                    if(file_pos < tng_data->input_file_len)
                    {
                        stat = tng_block_header_read(tng_data, block);
                    }
                }
            }
            tng_block_destroy(&block);
            if(stat == TNG_CRITICAL)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                return(stat);
            }

            for(i = frame_set->n_data_blocks; i-- ;)
            {
                data = &frame_set->tr_data[i];
                if(data->block_id == block_id)
                {
                    block_index = i;
                    break;
                }
            }
            if(block_index < 0)
            {
                return(TNG_FAILURE);
            }
        }
    }

    /* A bit hackish to create a new data struct before returning the data */
    new_data = malloc(sizeof(struct tng_non_particle_data));

    new_data->n_values_per_frame = 0;
    new_data->n_frames = 0;
    new_data->values = 0;
    new_data->datatype = data->datatype;
    *n_values_per_frame = data->n_values_per_frame;
    if(tng_allocate_data_mem(tng_data, new_data, data->n_frames,
                             data->stride_length, *n_values_per_frame)
       != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    *n_frames = tng_max(1, data->n_frames);

    *values = new_data->values;
    *type = data->datatype;
    switch(*type)
    {
    case TNG_CHAR_DATA:
        for(i=*n_frames; i--;)
        {
            for(j=*n_values_per_frame; j--;)
            {
                len = strlen(data->values[i][j].c) + 1;
                (*values)[i][j].c = malloc(len);
                strncpy((*values)[i][j].c, data->values[i][j].c, len);
            }
        }
        break;
    case TNG_INT_DATA:
        for(i=*n_frames; i--;)
        {
            for(j=*n_values_per_frame; j--;)
            {
                (*values)[i][j].i = data->values[i][j].i;
            }
        }
        break;
    case TNG_FLOAT_DATA:
        for(i=*n_frames; i--;)
        {
            for(j=*n_values_per_frame; j--;)
            {
                (*values)[i][j].f = data->values[i][j].f;
            }
        }
        break;
    case TNG_DOUBLE_DATA:
    default:
        for(i=*n_frames; i--;)
        {
            for(j=*n_values_per_frame; j--;)
            {
                (*values)[i][j].d = data->values[i][j].d;
            }
        }
    }

    free(new_data);

    return(TNG_SUCCESS);
}

tng_function_status tng_data_interval_get(tng_trajectory_t tng_data,
                                          const int64_t block_id,
                                          const int64_t start_frame_nr,
                                          const int64_t end_frame_nr,
                                          const tng_hash_mode hash_mode,
                                          union data_values ***values,
                                          int64_t *n_values_per_frame,
                                          tng_data_type *type)
{
    int64_t i, j, n_frames, file_pos, current_frame_pos;
    int block_index, len;
    tng_non_particle_data_t data, new_data;
    tng_trajectory_frame_set_t frame_set;
    tng_gen_block_t block;
    tng_function_status stat;

    block_index = -1;

    stat = tng_frame_set_of_frame_find(tng_data, start_frame_nr);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }


    tng_block_init(&block);
    file_pos = ftell(tng_data->input_file);
    /* Read all blocks until next frame set block */
    stat = tng_block_header_read(tng_data, block);
    while(file_pos < tng_data->input_file_len &&
            stat != TNG_CRITICAL &&
            block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        stat = tng_block_read_next(tng_data, block,
                                   hash_mode);
        if(stat != TNG_CRITICAL)
        {
            file_pos = ftell(tng_data->input_file);
            if(file_pos < tng_data->input_file_len)
            {
                stat = tng_block_header_read(tng_data, block);
            }
        }
    }
    tng_block_destroy(&block);
    if(stat == TNG_CRITICAL)
    {
        printf("Cannot read block header at pos %"PRId64". %s: %d\n",
               file_pos, __FILE__, __LINE__);
        return(stat);
    }

    frame_set = &tng_data->current_trajectory_frame_set;

    /* See if there is a data block of this ID.
     * Start checking the last read frame set */
    for(i = frame_set->n_data_blocks; i-- ;)
    {
        data = &frame_set->tr_data[i];
        if(data->block_id == block_id)
        {
            block_index = i;
            break;
        }
    }

    if(block_index < 0)
    {
        printf("Could not find non-particle data block with id %"PRId64". %s: %d\n",
                block_id, __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    n_frames = end_frame_nr - start_frame_nr + 1;

    /* A bit hackish to create a new data struct before returning the data */

    if(*values == 0)
    {
        new_data = malloc(sizeof(struct tng_non_particle_data));

        new_data->n_values_per_frame = 0;
        new_data->n_frames = 0;
        new_data->values = 0;
        new_data->datatype = data->datatype;
        *n_values_per_frame = data->n_values_per_frame;
        if(tng_allocate_data_mem(tng_data, new_data, n_frames,
                                 data->stride_length,
                                 data->n_values_per_frame) != TNG_SUCCESS)
        {
            free(new_data);
            return(TNG_CRITICAL);
        }

        *values = new_data->values;

        free(new_data);
    }

    *type = data->datatype;
    current_frame_pos = start_frame_nr - frame_set->first_frame;
    /* It's not very elegant to reuse so much of the code in the different case
     * statements, but it's unnecessarily slow to have the switch-case block
     * inside the for loops. */
    switch(*type)
    {
    case TNG_CHAR_DATA:
        for(i=0; i<n_frames; i++)
        {
            if(current_frame_pos == frame_set->n_frames)
            {
                stat = tng_frame_set_read_next(tng_data, hash_mode);
                if(stat != TNG_SUCCESS)
                {
                    return(stat);
                }
                current_frame_pos = 0;
            }
            for(j=*n_values_per_frame; j--;)
            {
                len = strlen(data->values[current_frame_pos][j].c) + 1;
                (*values)[i][j].c = malloc(len);
                strncpy((*values)[i][j].c, data->values[current_frame_pos][j].c, len);
            }
            current_frame_pos++;
        }
        break;
    case TNG_INT_DATA:
        for(i=0; i<n_frames; i++)
        {
            if(current_frame_pos == frame_set->n_frames)
            {
                stat = tng_frame_set_read_next(tng_data, hash_mode);
                if(stat != TNG_SUCCESS)
                {
                    return(stat);
                }
                current_frame_pos = 0;
            }
            for(j=*n_values_per_frame; j--;)
            {
                (*values)[i][j].i = data->values[current_frame_pos][j].i;
            }
            current_frame_pos++;
        }
        break;
    case TNG_FLOAT_DATA:
        for(i=0; i<n_frames; i++)
        {
            if(current_frame_pos == frame_set->n_frames)
            {
                stat = tng_frame_set_read_next(tng_data, hash_mode);
                if(stat != TNG_SUCCESS)
                {
                    return(stat);
                }
                current_frame_pos = 0;
            }
            for(j=*n_values_per_frame; j--;)
            {
                (*values)[i][j].f = data->values[current_frame_pos][j].f;
            }
            current_frame_pos++;
        }
        break;
    case TNG_DOUBLE_DATA:
    default:
        for(i=0; i<n_frames; i++)
        {
            if(current_frame_pos == frame_set->n_frames)
            {
                stat = tng_frame_set_read_next(tng_data, hash_mode);
                if(stat != TNG_SUCCESS)
                {
                    return(stat);
                }
                current_frame_pos = 0;
            }
            for(j=*n_values_per_frame; j--;)
            {
                (*values)[i][j].d = data->values[current_frame_pos][j].d;
            }
            current_frame_pos++;
        }
    }

    return(TNG_SUCCESS);
}


tng_function_status tng_particle_data_get(tng_trajectory_t tng_data,
                                          const int64_t block_id,
                                          union data_values ****values,
                                          int64_t *n_frames,
                                          int64_t *n_particles,
                                          int64_t *n_values_per_frame,
                                          tng_data_type *type)
{
    int64_t i, j, k, mapping, file_pos;
    int block_index, len;
    tng_particle_data_t data, new_data;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;
    tng_gen_block_t block;
    tng_function_status stat;

    tng_block_type block_type_flag;

    block_index = -1;

    /* See if there is already a data block of this ID.
     * Start checking the last read frame set */
    for(i = frame_set->n_particle_data_blocks; i-- ;)
    {
        data = &frame_set->tr_particle_data[i];
        if(data->block_id == block_id)
        {
            block_index = i;
            block_type_flag = TNG_TRAJECTORY_BLOCK;
            break;
        }
    }

    if(block_index < 0)
    {
        /* If the data block was not found in the frame set
         * look for it in the non-trajectory data (in tng_data). */
        for(i = tng_data->n_particle_data_blocks; i-- ;)
        {
            data = &tng_data->non_tr_particle_data[i];
            if(data->block_id == block_id)
            {
                block_index = i;
                block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
                break;
            }
        }
        if(block_index < 0)
        {
            tng_block_init(&block);
            file_pos = ftell(tng_data->input_file);
            /* Read all blocks until next frame set block */
            stat = tng_block_header_read(tng_data, block);
            while(file_pos < tng_data->input_file_len &&
                    stat != TNG_CRITICAL &&
                    block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                /* Use hash by default */
                stat = tng_block_read_next(tng_data, block,
                                        TNG_USE_HASH);
                if(stat != TNG_CRITICAL)
                {
                    file_pos = ftell(tng_data->input_file);
                    if(file_pos < tng_data->input_file_len)
                    {
                        stat = tng_block_header_read(tng_data, block);
                    }
                }
            }
            tng_block_destroy(&block);
            if(stat == TNG_CRITICAL)
            {
                printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                       file_pos, __FILE__, __LINE__);
                return(stat);
            }

            for(i = frame_set->n_particle_data_blocks; i-- ;)
            {
                data = &frame_set->tr_particle_data[i];
                if(data->block_id == block_id)
                {
                    block_index = i;
                    block_type_flag = TNG_TRAJECTORY_BLOCK;
                    break;
                }
            }
            if(block_index < 0)
            {
                return(TNG_FAILURE);
            }
        }
    }

    if(block_type_flag == TNG_TRAJECTORY_BLOCK &&
       tng_data->var_num_atoms_flag)
    {
        *n_particles = frame_set->n_particles;
    }
    else
    {
        *n_particles = tng_data->n_particles;
    }

    /* A bit hackish to create a new data struct before returning the data */

    if(*values == 0)
    {
        new_data = malloc(sizeof(struct tng_particle_data));

        new_data->n_values_per_frame = 0;
        new_data->n_frames = 0;
        new_data->values = 0;
        new_data->datatype = data->datatype;
        *n_values_per_frame = data->n_values_per_frame;
        if(tng_allocate_particle_data_mem(tng_data, new_data, data->n_frames,
                                          data->stride_length,
                                        *n_particles, data->n_values_per_frame)
            != TNG_SUCCESS)
        {
            free(new_data);
            return(TNG_CRITICAL);
        }

        *n_frames = tng_max(1, data->n_frames);

        *values = new_data->values;

        free(new_data);
    }

    *type = data->datatype;
    /* It's not very elegant to reuse so much of the code in the different case
     * statements, but it's unnecessarily slow to have the switch-case block
     * inside the for loops. */
    switch(*type)
    {
    case TNG_CHAR_DATA:
        for(i=*n_frames; i--;)
        {
            for(j=*n_particles; j--;)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k=*n_values_per_frame; k--;)
                {
                    len = strlen(data->values[i][j][k].c) + 1;
                    (*values)[i][j][k].c = malloc(len);
                    strncpy((*values)[i][mapping][k].c, data->values[i][j][k].c, len);
                }
            }
        }
        break;
    case TNG_INT_DATA:
        for(i=*n_frames; i--;)
        {
            for(j=*n_particles; j--;)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k=*n_values_per_frame; k--;)
                {
                    (*values)[i][mapping][k].i = data->values[i][j][k].i;
                }
            }
        }
        break;
    case TNG_FLOAT_DATA:
        for(i=*n_frames; i--;)
        {
            for(j=*n_particles; j--;)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k=*n_values_per_frame; k--;)
                {
                    (*values)[i][mapping][k].f = data->values[i][j][k].f;
                }
            }
        }
        break;
    case TNG_DOUBLE_DATA:
    default:
        for(i=*n_frames; i--;)
        {
            for(j=*n_particles; j--;)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k=*n_values_per_frame; k--;)
                {
                    (*values)[i][mapping][k].d = data->values[i][j][k].d;
                }
            }
        }
    }

    return(TNG_SUCCESS);
}


tng_function_status tng_particle_data_interval_get(tng_trajectory_t tng_data,
                                                  const int64_t block_id,
                                                  const int64_t start_frame_nr,
                                                  const int64_t end_frame_nr,
                                                  const tng_hash_mode hash_mode,
                                                  union data_values ****values,
                                                  int64_t *n_particles,
                                                  int64_t *n_values_per_frame,
                                                  tng_data_type *type)
{
    int64_t i, j, k, mapping, n_frames, file_pos, current_frame_pos;
    int block_index, len;
    tng_particle_data_t data, new_data;
    tng_trajectory_frame_set_t frame_set;
    tng_gen_block_t block;
    tng_function_status stat;
    tng_block_type block_type_flag;

    block_index = -1;

    stat = tng_frame_set_of_frame_find(tng_data, start_frame_nr);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }


    tng_block_init(&block);
    file_pos = ftell(tng_data->input_file);
    /* Read all blocks until next frame set block */
    stat = tng_block_header_read(tng_data, block);
    while(file_pos < tng_data->input_file_len &&
            stat != TNG_CRITICAL &&
            block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        stat = tng_block_read_next(tng_data, block,
                                   hash_mode);
        if(stat != TNG_CRITICAL)
        {
            file_pos = ftell(tng_data->input_file);
            if(file_pos < tng_data->input_file_len)
            {
                stat = tng_block_header_read(tng_data, block);
            }
        }
    }
    tng_block_destroy(&block);
    if(stat == TNG_CRITICAL)
    {
        printf("Cannot read block header at pos %"PRId64". %s: %d\n",
                file_pos, __FILE__, __LINE__);
        return(stat);
    }

    frame_set = &tng_data->current_trajectory_frame_set;

    /* See if there is already a data block of this ID.
     * Start checking the last read frame set */
    for(i = frame_set->n_particle_data_blocks; i-- ;)
    {
        data = &frame_set->tr_particle_data[i];
        if(data->block_id == block_id)
        {
            block_index = i;
            block_type_flag = TNG_TRAJECTORY_BLOCK;
            break;
        }
    }

    if(block_index < 0)
    {
        printf("Could not find particle data block with id %"PRId64". %s: %d\n",
                block_id, __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    if(block_type_flag == TNG_TRAJECTORY_BLOCK &&
       tng_data->var_num_atoms_flag)
    {
        *n_particles = frame_set->n_particles;
    }
    else
    {
        *n_particles = tng_data->n_particles;
    }

    n_frames = end_frame_nr - start_frame_nr + 1;

    /* A bit hackish to create a new data struct before returning the data */

    if(*values == 0)
    {
        new_data = malloc(sizeof(struct tng_particle_data));

        new_data->n_values_per_frame = 0;
        new_data->n_frames = 0;
        new_data->values = 0;
        new_data->datatype = data->datatype;
        *n_values_per_frame = data->n_values_per_frame;
        if(tng_allocate_particle_data_mem(tng_data, new_data, n_frames,
                                          data->stride_length,
                                          *n_particles, data->n_values_per_frame) !=
           TNG_SUCCESS)
        {
            free(new_data);
            return(TNG_CRITICAL);
        }

        *values = new_data->values;

        free(new_data);
    }

    *type = data->datatype;
    current_frame_pos = start_frame_nr - frame_set->first_frame;
    /* It's not very elegant to reuse so much of the code in the different case
     * statements, but it's unnecessarily slow to have the switch-case block
     * inside the for loops. */
    switch(*type)
    {
    case TNG_CHAR_DATA:
        for(i=0; i<n_frames; i++)
        {
            if(current_frame_pos == frame_set->n_frames)
            {
                stat = tng_frame_set_read_next(tng_data, hash_mode);
                if(stat != TNG_SUCCESS)
                {
                    return(stat);
                }
                current_frame_pos = 0;
            }
            for(j=*n_particles; j--;)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k=*n_values_per_frame; k--;)
                {
                    len = strlen(data->values[current_frame_pos][j][k].c) + 1;
                    (*values)[i][j][k].c = malloc(len);
                    strncpy((*values)[i][j][k].c, data->values[current_frame_pos][j][k].c, len);
                }
            }
            current_frame_pos++;
        }
        break;
    case TNG_INT_DATA:
        for(i=0; i<n_frames; i++)
        {
            if(current_frame_pos == frame_set->n_frames)
            {
                stat = tng_frame_set_read_next(tng_data, hash_mode);
                if(stat != TNG_SUCCESS)
                {
                    return(stat);
                }
                current_frame_pos = 0;
            }
            for(j=*n_particles; j--;)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k=*n_values_per_frame; k--;)
                {
                    (*values)[i][mapping][k].i = data->values[current_frame_pos][j][k].i;
                }
            }
            current_frame_pos++;
        }
        break;
    case TNG_FLOAT_DATA:
        for(i=0; i<n_frames; i++)
        {
            if(current_frame_pos == frame_set->n_frames)
            {
                stat = tng_frame_set_read_next(tng_data, hash_mode);
                if(stat != TNG_SUCCESS)
                {
                    return(stat);
                }
                current_frame_pos = 0;
            }
            for(j=*n_particles; j--;)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k=*n_values_per_frame; k--;)
                {
                    (*values)[i][mapping][k].f = data->values[current_frame_pos][j][k].f;
                }
            }
            current_frame_pos++;
        }
        break;
    case TNG_DOUBLE_DATA:
    default:
        for(i=0; i<n_frames; i++)
        {
            if(current_frame_pos == frame_set->n_frames)
            {
                stat = tng_frame_set_read_next(tng_data, hash_mode);
                if(stat != TNG_SUCCESS)
                {
                    return(stat);
                }
                current_frame_pos = 0;
            }
            for(j=*n_particles; j--;)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k=*n_values_per_frame; k--;)
                {
                    (*values)[i][mapping][k].d = data->values[current_frame_pos][j][k].d;
                }
            }
            current_frame_pos++;
        }
    }

    return(TNG_SUCCESS);
}


tng_function_status tng_time_get_str(const tng_trajectory_t tng_data,
                                     char *time)
{
    struct tm *time_data;
    time_t secs;

    secs = tng_data->time;

    time_data = localtime(&secs); // Returns a statically allocated variable.
    snprintf(time, TNG_MAX_DATE_STR_LEN,
             "%4d-%02d-%02d %02d:%02d:%02d",
             time_data->tm_year+1900, time_data->tm_mon+1, time_data->tm_mday,
             time_data->tm_hour, time_data->tm_min, time_data->tm_sec);

    return(TNG_SUCCESS);
}


#ifdef BUILD_FORTRAN
/* The following is for calling the library from fortran */

tng_function_status tng_trajectory_init_(tng_trajectory_t *tng_data_p)
{
    return(tng_trajectory_init(tng_data_p));
}

tng_function_status tng_trajectory_destroy_(tng_trajectory_t *tng_data_p)
{
    return(tng_trajectory_destroy(tng_data_p));
}

tng_function_status tng_trajectory_init_from_src_(tng_trajectory_t src,
                                                   tng_trajectory_t *dest_p)
{
    return(tng_trajectory_init_from_src(src, dest_p));
}

tng_function_status tng_input_file_get_(const tng_trajectory_t tng_data,
                                        char *file_name, const int max_len)
{
    return(tng_input_file_get(tng_data, file_name, max_len));
}

tng_function_status tng_input_file_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_output_file_get_(const tng_trajectory_t tng_data,
                                         char *file_name, const int max_len)
{
    return(tng_output_file_get(tng_data, file_name, max_len));
}

tng_function_status tng_output_file_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_first_program_name_get_(const tng_trajectory_t tng_data,
                                                char *name, const int max_len)
{
    return(tng_first_program_name_get(tng_data, name, max_len));
}

tng_function_status tng_output_file_endianness_get_
                (tng_trajectory_t tng_data, tng_file_endianness *endianness)
{
    return(tng_output_file_endianness_get(tng_data, endianness));
}

tng_function_status tng_output_file_endianness_set_
                (tng_trajectory_t tng_data, const tng_file_endianness *endianness)
{
    return(tng_output_file_endianness_set(tng_data, *endianness));
}

tng_function_status tng_first_program_name_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_last_program_name_get_(const tng_trajectory_t tng_data,
                                                char *name, const int max_len)
{
    return(tng_last_program_name_get(tng_data, name, max_len));
}

tng_function_status tng_last_program_name_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_first_user_name_get_(const tng_trajectory_t tng_data,
                                             char *name, const int max_len)
{
    return(tng_first_user_name_get(tng_data, name, max_len));
}

tng_function_status tng_first_user_name_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_last_user_name_get_(const tng_trajectory_t tng_data,
                                            char *name, const int max_len)
{
    return(tng_last_user_name_get(tng_data, name, max_len));
}

tng_function_status tng_last_user_name_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_first_computer_name_get_(const tng_trajectory_t tng_data,
                                                 char *name, const int max_len)
{
    return(tng_first_computer_name_get(tng_data, name, max_len));
}

tng_function_status tng_first_computer_name_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_last_computer_name_get_(const tng_trajectory_t tng_data,
                                                 char *name, const int max_len)
{
    return(tng_last_computer_name_get(tng_data, name, max_len));
}

tng_function_status tng_last_computer_name_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_first_signature_get_(const tng_trajectory_t tng_data,
                                             char *signature, const int max_len)
{
    return(tng_first_signature_get(tng_data, signature, max_len));
}

tng_function_status tng_first_signature_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_last_signature_get_(const tng_trajectory_t tng_data,
                                             char *signature, const int max_len)
{
    return(tng_last_signature_get(tng_data, signature, max_len));
}

tng_function_status tng_last_signature_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_forcefield_name_get_(const tng_trajectory_t tng_data,
                                             char *name, const int max_len)
{
    return(tng_forcefield_name_get(tng_data, name, max_len));
}

tng_function_status tng_forcefield_name_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_medium_stride_length_get_(const tng_trajectory_t tng_data,
                                                  int64_t *len)
{
    return(tng_medium_stride_length_get(tng_data, len));
}

tng_function_status tng_medium_stride_length_set_(tng_trajectory_t tng_data,
                                                  const int64_t *len)
{
    return(tng_medium_stride_length_set(tng_data, *len));
}

tng_function_status tng_long_stride_length_get_(const tng_trajectory_t tng_data,
                                                int64_t *len)
{
    return(tng_long_stride_length_get(tng_data, len));
}

tng_function_status tng_long_stride_length_set_(tng_trajectory_t tng_data,
                                                const int64_t *len)
{
    return(tng_long_stride_length_set(tng_data, *len));
}

tng_function_status tng_input_file_len_get_(const tng_trajectory_t tng_data,
                                            int64_t *len)
{
    return(tng_input_file_len_get(tng_data, len));
}

tng_function_status tng_num_frames_get_(const tng_trajectory_t tng_data,
                                     int64_t *n)
{
    return(tng_num_frames_get(tng_data, n));
}

tng_function_status tng_num_particles_get_(const tng_trajectory_t tng_data,
                                           int64_t *n)
{
    return(tng_num_particles_get(tng_data, n));
}

tng_function_status tng_num_molecules_get_(const tng_trajectory_t tng_data,
                                           int64_t *n)
{
    return(tng_num_molecules_get(tng_data, n));
}

tng_function_status tng_num_frames_per_frame_set_get_
                (const tng_trajectory_t tng_data,
                 int64_t *n)
{
    return(tng_num_frames_per_frame_set_get(tng_data, n));
}

tng_function_status tng_num_frames_per_frame_set_set_
                (const tng_trajectory_t tng_data,
                 int64_t *n)
{
    return(tng_num_frames_per_frame_set_set(tng_data, *n));
}

tng_function_status tng_num_frame_sets_get_
                (const tng_trajectory_t tng_data,
                 int64_t *n)
{
    return(tng_num_frame_sets_get(tng_data, n));
}

tng_function_status tng_current_frame_set_get_
                (tng_trajectory_t tng_data,
                 tng_trajectory_frame_set_t *frame_set_p)
{
    return(tng_current_frame_set_get(tng_data, frame_set_p));
}

tng_function_status tng_frame_set_nr_find_(tng_trajectory_t tng_data,
                                        const int64_t *nr)
{
    return(tng_frame_set_nr_find(tng_data, *nr));
}

tng_function_status tng_frame_set_of_frame_find_(tng_trajectory_t tng_data,
                                        const int64_t *frame)
{
    return(tng_frame_set_of_frame_find(tng_data, *frame));
}

tng_function_status tng_frame_set_next_frame_set_file_pos_get_
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *pos)
{
    return(tng_frame_set_next_frame_set_file_pos_get(tng_data, frame_set, pos));
}

tng_function_status tng_frame_set_prev_frame_set_file_pos_get_
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *pos)
{
    return(tng_frame_set_prev_frame_set_file_pos_get(tng_data, frame_set, pos));
}

tng_function_status tng_frame_set_frame_range_get_
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *first_frame,
                 int64_t *last_frame)
{
    return(tng_frame_set_frame_range_get(tng_data, frame_set, first_frame,
                                         last_frame));
}

tng_function_status tng_molecule_init_(const tng_trajectory_t tng_data,
                                       tng_molecule_t molecule)
{
    return(tng_molecule_init(tng_data, molecule));
}

tng_function_status tng_molecule_destroy_(const tng_trajectory_t tng_data,
                                          tng_molecule_t molecule)
{
    return(tng_molecule_destroy(tng_data, molecule));
}

tng_function_status tng_molecule_add_(tng_trajectory_t tng_data,
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

tng_function_status tng_molecule_name_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_molecule_cnt_get_(tng_trajectory_t tng_data,
                                          tng_molecule_t molecule,
                                          int64_t *cnt)
{
    return(tng_molecule_cnt_get(tng_data, molecule, cnt));
}

tng_function_status tng_molecule_cnt_set_(tng_trajectory_t tng_data,
                                          tng_molecule_t molecule,
                                          int64_t *cnt)
{
    return(tng_molecule_cnt_set(tng_data, molecule, *cnt));
}

tng_function_status tng_molecule_chain_find_(tng_trajectory_t tng_data,
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

tng_function_status tng_molecule_chain_add_(tng_trajectory_t tng_data,
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

tng_function_status tng_chain_name_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_chain_residue_add_(tng_trajectory_t tng_data,
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

tng_function_status tng_residue_name_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_residue_atom_add_(tng_trajectory_t tng_data,
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

tng_function_status tng_atom_name_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_atom_type_set_(tng_trajectory_t tng_data,
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

tng_function_status tng_molecule_name_of_particle_nr_get_
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    return(tng_molecule_name_of_particle_nr_get(tng_data, nr, name, max_len));
}

tng_function_status tng_chain_name_of_particle_nr_get_
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    return(tng_chain_name_of_particle_nr_get(tng_data, nr, name, max_len));
}

tng_function_status tng_residue_name_of_particle_nr_get_
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    return(tng_residue_name_of_particle_nr_get(tng_data, nr, name, max_len));
}

tng_function_status tng_atom_name_of_particle_nr_get_
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    return(tng_atom_name_of_particle_nr_get(tng_data, nr, name, max_len));
}

tng_function_status tng_atom_type_of_particle_nr_get_
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *type,
                 int max_len)
{
    return(tng_atom_type_of_particle_nr_get(tng_data, nr, type, max_len));
}

tng_function_status tng_particle_mapping_add_
                (tng_trajectory_t tng_data,
                 const int64_t *first_particle_number,
                 const int64_t *n_particles,
                 const int64_t *mapping_table)
{
    return(tng_particle_mapping_add(tng_data, *first_particle_number,
                                    *n_particles, mapping_table));
}

tng_function_status tng_file_headers_read_(tng_trajectory_t tng_data,
                                           const tng_hash_mode *hash_mode)
{
    return(tng_file_headers_read(tng_data, *hash_mode));
}

tng_function_status tng_file_headers_write_(tng_trajectory_t tng_data,
                                            const tng_hash_mode *hash_mode)
{
    return(tng_file_headers_write(tng_data, *hash_mode));
}

tng_function_status tng_block_read_next_(tng_trajectory_t tng_data,
                                         tng_gen_block_t block_data,
                                         const tng_hash_mode *hash_mode)
{
    return(tng_block_read_next(tng_data, block_data, *hash_mode));
}

tng_function_status tng_frame_set_read_next_(tng_trajectory_t tng_data,
                                             const tng_hash_mode *hash_mode)
{
    return(tng_frame_set_read_next(tng_data, *hash_mode));
}

tng_function_status tng_frame_set_write_(tng_trajectory_t tng_data,
                                         const tng_hash_mode *hash_mode)
{
    return(tng_frame_set_write(tng_data, *hash_mode));
}

tng_function_status tng_frame_set_new_(tng_trajectory_t tng_data,
                                       const int64_t *first_frame,
                                       const int64_t *n_frames)
{
    return(tng_frame_set_new(tng_data, *first_frame, *n_frames));
}

tng_function_status tng_data_block_add_(tng_trajectory_t tng_data,
                                        const int64_t *id,
                                        const char *block_name,
                                        const tng_data_type *datatype,
                                        const tng_block_type *block_type_flag,
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

tng_function_status tng_particle_data_block_add_
                (tng_trajectory_t tng_data,
                 const int64_t *id,
                 const char *block_name,
                 const tng_data_type *datatype,
                 const tng_block_type *block_type_flag,
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

tng_function_status tng_frame_data_write_(tng_trajectory_t tng_data,
                                          const int64_t *frame_nr,
                                          const int64_t *block_id,
                                          const void *data,
                                          const tng_hash_mode *hash_mode)
{
    return(tng_frame_data_write(tng_data, *frame_nr, *block_id, data,
                                *hash_mode));
}

tng_function_status tng_frame_particle_data_write_(tng_trajectory_t tng_data,
                                                   const int64_t *frame_nr,
                                                   const int64_t *block_id,
                                                   const int64_t *val_first_particle,
                                                   const int64_t *val_n_particles,
                                                   const void *data,
                                                   const tng_hash_mode *hash_mode)
{
    return(tng_frame_particle_data_write(tng_data, *frame_nr, *block_id,
                                         *val_first_particle, *val_n_particles,
                                         data, *hash_mode));
}

tng_function_status tng_data_values_free_(const tng_trajectory_t tng_data,
                                          union data_values **values,
                                          const int64_t *n_frames,
                                          const int64_t *n_values_per_frame,
                                          const tng_data_type *type)
{
    return(tng_data_values_free(tng_data, values, *n_frames,
                                *n_values_per_frame, *type));
}

tng_function_status tng_particle_data_values_free_(const tng_trajectory_t tng_data,
                                                   union data_values ***values,
                                                   const int64_t *n_frames,
                                                   const int64_t *n_particles,
                                                   const int64_t *n_values_per_frame,
                                                   const tng_data_type *type)
{
    return(tng_particle_data_values_free(tng_data, values, *n_frames, *n_particles,
                                         *n_values_per_frame, *type));
}

tng_function_status tng_data_get_(tng_trajectory_t tng_data,
                                  const int64_t *block_id,
                                  union data_values ***values,
                                  int64_t *n_frames,
                                  int64_t *n_values_per_frame,
                                  tng_data_type *type)
{
    return(tng_data_get(tng_data, *block_id, values, n_frames,
                        n_values_per_frame, type));
}

tng_function_status tng_data_interval_get_(tng_trajectory_t tng_data,
                                           const int64_t *block_id,
                                           const int64_t *start_frame_nr,
                                           const int64_t *end_frame_nr,
                                           const tng_hash_mode *hash_mode,
                                           union data_values ***values,
                                           int64_t *n_values_per_frame,
                                           tng_data_type *type)
{
    return(tng_data_interval_get(tng_data, *block_id, *start_frame_nr,
                                 *end_frame_nr, *hash_mode, values,
                                 n_values_per_frame, type));
}

tng_function_status tng_particle_data_get_(tng_trajectory_t tng_data,
                                           const int64_t *block_id,
                                           union data_values ****values,
                                           int64_t *n_frames,
                                           int64_t *n_particles,
                                           int64_t *n_values_per_frame,
                                           tng_data_type *type)
{
    return(tng_particle_data_get(tng_data, *block_id, values, n_frames,
                                 n_particles, n_values_per_frame, type));
}

tng_function_status tng_particle_data_interval_get_(tng_trajectory_t tng_data,
                                                    const int64_t *block_id,
                                                    const int64_t *start_frame_nr,
                                                    const int64_t *end_frame_nr,
                                                    const tng_hash_mode *hash_mode,
                                                    union data_values ****values,
                                                    int64_t *n_particles,
                                                    int64_t *n_values_per_frame,
                                                    tng_data_type *type)
{
    return(tng_particle_data_interval_get(tng_data, *block_id, *start_frame_nr,
                                          *end_frame_nr, *hash_mode, values,
                                          n_particles, n_values_per_frame,
                                          type));
}

tng_function_status tng_time_get_str_(const tng_trajectory_t tng_data,
                                      char *time, int64_t str_len)
{
    return(tng_time_get_str(tng_data, time));
}
#endif

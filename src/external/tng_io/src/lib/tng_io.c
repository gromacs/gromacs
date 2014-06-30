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

#ifdef USE_STD_INTTYPES_H
#include <inttypes.h>
#endif

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#ifdef USE_ZLIB
#include <zlib.h>
#endif

#include "tng/tng_io.h"
#include "tng/md5.h"
#include "compression/tng_compress.h"
#include "tng/version.h"


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
    char md5_hash[TNG_MD5_HASH_LEN];
    /** The name of the block */
    char *name;
    /** The library version used to write the block */
    int64_t block_version;
    int64_t alt_hash_type;
    int64_t alt_hash_len;
    char *alt_hash;
    int64_t signature_type;
    int64_t signature_len;
    char *signature;
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
     * frame at a time). */
    int64_t n_written_frames;
    /** The number of frames not yet written to file in this frame set
     * (used from the utility functions to finish the writing properly. */
    int64_t n_unwritten_frames;


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
    /** Time stamp (in seconds) of first frame in frame set */
    double first_frame_time;

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
/* FIXME: Make only one data block struct */
struct tng_particle_data {
    /** The block ID of the data block containing this particle data.
     *  This is used to determine the kind of data that is stored */
    int64_t block_id;
    /** The name of the data block. This is used to determine the kind of
     *  data that is stored */
    char *block_name;
    /** The type of data stored. */
    char datatype;
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
    /** If reading one frame at a time this is the last read frame */
    int64_t last_retrieved_frame;
    /** The multiplier used for getting integer values for compression */
    double compression_multiplier;
    /** A 1-dimensional array of values of length
     *  [sizeof (datatype)] * n_frames * n_particles * n_values_per_frame */
    void *values;
    /** If storing character data store it in a 3-dimensional array */
    char ****strings;
};

struct tng_non_particle_data {
    /** The ID of the data block */
    int64_t block_id;
    /** The name of the data block. This is used to determine the kind of
     *  data that is stored */
    char *block_name;
    /** The type of data stored. */
    char datatype;
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
    /** If reading one frame at a time this is the last read frame */
    int64_t last_retrieved_frame;
    /** Compressed data is stored as integers. This compression multiplier is
     *  the multiplication factor to convert from integer to float/double */
    double compression_multiplier;
    /** A 1-dimensional array of values of length
     *  [sizeof (datatype)] * n_frames * n_values_per_frame */
    void *values;
    /** If storing character data store it in a 2-dimensional array */
    char ***strings;
};



struct tng_trajectory {
    /** The path of the input trajectory file */
    char *input_file_path;
    /** A handle to the input file */
    FILE *input_file;
    /** The length of the input file */
    long input_file_len;
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
    char endianness_32;
    /** The endianness of 64 bit values of the current computer */
    char endianness_64;

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
    /** The exponential of the value of the distance unit used. The default
     * distance unit is nm (1e-9), i.e. distance_unit_exponential = -9. If
     * the measurements are in Å the distance_unit_exponential = -10. */
    int64_t distance_unit_exponential;

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
    /** The current (can change from one frame set to another) time length
     *  (in seconds) of one frame */
    double time_per_frame;

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
    long current_trajectory_frame_set_input_file_pos;
    /** The pos in the dest file of the current frame set */
    long current_trajectory_frame_set_output_file_pos;
    /** The number of frame sets in the trajectory N.B. Not saved in file and
     *  cannot be trusted to be up-to-date */
    int64_t n_trajectory_frame_sets;

    /* These data blocks are non-trajectory data blocks */
    /** The number of non-frame dependent particle dependent data blocks */
    int n_particle_data_blocks;
    /** A list of data blocks containing particle dependent data */
    struct tng_particle_data *non_tr_particle_data;

    /** The number of frame and particle independent data blocks */
    int n_data_blocks;
    /** A list of frame and particle indepdendent data blocks */
    struct tng_non_particle_data *non_tr_data;

    /** TNG compression algorithm for compressing positions */
    int *compress_algo_pos;
    /** TNG compression algorithm for compressing velocities */
    int *compress_algo_vel;
    /** The precision used for lossy compression */
    double compression_precision;
};

#ifndef USE_WINDOWS
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#define USE_WINDOWS
#endif /* win32... */
#endif /* not defined USE_WINDOWS */

#ifdef USE_WINDOWS
#define TNG_INLINE __inline
#define TNG_SNPRINTF _snprintf
#else
#define TNG_INLINE inline
#define TNG_SNPRINTF snprintf
#endif

static TNG_INLINE int tng_min_i(int a, int b)
{
    return (a < b ? a : b);
}

/*
static TNG_INLINE int tng_max_i(int a, int b)
{
    return (a > b ? a : b);
}
*/
static TNG_INLINE int64_t tng_min_i64(int64_t a, int64_t b)
{
    return (a < b ? a : b);
}

static TNG_INLINE int64_t tng_max_i64(int64_t a, int64_t b)
{
    return (a > b ? a : b);
}

/*
static TNG_INLINE float tng_min_f(float a, float b)
{
    return (a < b ? a : b);
}

static TNG_INLINE float tng_max_f(float a, float b)
{
    return (a > b ? a : b);
}

static TNG_INLINE double tng_min_d(double a, double b)
{
    return (a < b ? a : b);
}

static TNG_INLINE double tng_max_d(double a, double b)
{
    return (a > b ? a : b);
}
*/

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
    case TNG_LITTLE_ENDIAN_32: /* Byte order is reversed. */
        *v = ((*v & 0xFF000000) >> 24) | /* Move first byte to end */
             ((*v & 0x00FF0000) >> 8) |  /* Move 2nd byte to pos 3 */
             ((*v & 0x0000FF00) << 8) |  /* Move 3rd byte to pos 2 */
             ((*v & 0x000000FF) << 24);  /* Move last byte to first */

        return(TNG_SUCCESS);

    case TNG_BYTE_PAIR_SWAP_32: /* byte pair swap */
        *v = ((*v & 0xFFFF0000) >> 16) |
             ((*v & 0x0000FFFF) << 16);

        return(TNG_SUCCESS);

    case TNG_BIG_ENDIAN_32: /* Already correct */
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
    case TNG_LITTLE_ENDIAN_64: /* Byte order is reversed. */
        *v = ((*v & 0xFF00000000000000LL) >> 56) | /* Move first byte to end */
             ((*v & 0x00FF000000000000LL) >> 40) | /* Move 2nd byte to pos 7 */
             ((*v & 0x0000FF0000000000LL) >> 24) | /* Move 3rd byte to pos 6 */
             ((*v & 0x000000FF00000000LL) >> 8 ) | /* Move 4th byte to pos 5 */
             ((*v & 0x00000000FF000000LL) << 8 ) | /* Move 5th byte to pos 4 */
             ((*v & 0x0000000000FF0000LL) << 24) | /* Move 6th byte to pos 3 */
             ((*v & 0x000000000000FF00LL) << 40) | /* Move 7th byte to pos 2 */
             ((*v & 0x00000000000000FFLL) << 56);  /* Move last byte to first */

        return(TNG_SUCCESS);

    case TNG_QUAD_SWAP_64: /* Byte quad swap */
        *v = ((*v & 0xFFFFFFFF00000000LL) >> 32) |
             ((*v & 0x00000000FFFFFFFFLL) << 32);

        return(TNG_SUCCESS);

    case TNG_BYTE_PAIR_SWAP_64: /* Byte pair swap */
        *v = ((*v & 0xFFFF0000FFFF0000LL) >> 16) |
             ((*v & 0x0000FFFF0000FFFFLL) << 16);

        return(TNG_SUCCESS);

    case TNG_BYTE_SWAP_64: /* Byte swap */
        *v = ((*v & 0xFF00FF00FF00FF00LL) >> 8) |
             ((*v & 0x00FF00FF00FF00FFLL) << 8);

        return(TNG_SUCCESS);

    case TNG_BIG_ENDIAN_64: /* Already correct */
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
    case TNG_LITTLE_ENDIAN_32: /* Already correct */
        return(TNG_SUCCESS);

    case TNG_BYTE_PAIR_SWAP_32: /* byte pair swapped big endian to little endian */
        *v = ((*v & 0xFF00FF00) >> 8) |
             ((*v & 0x00FF00FF) << 8);

        return(TNG_SUCCESS);

    case TNG_BIG_ENDIAN_32: /* Byte order is reversed. */
        *v = ((*v & 0xFF000000) >> 24) | /* Move first byte to end */
             ((*v & 0x00FF0000) >> 8) |  /* Move 2nd byte to pos 3 */
             ((*v & 0x0000FF00) << 8) |  /* Move 3rd byte to pos 2 */
             ((*v & 0x000000FF) << 24);  /* Move last byte to first */

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
    case TNG_LITTLE_ENDIAN_64: /* Already correct */
        return(TNG_SUCCESS);

    case TNG_QUAD_SWAP_64: /* Byte quad swapped big endian to little endian */
        *v = ((*v & 0xFF000000FF000000LL) >> 24) |
             ((*v & 0x00FF000000FF0000LL) >> 8) |
             ((*v & 0x0000FF000000FF00LL) << 8) |
             ((*v & 0x000000FF000000FFLL) << 24);

        return(TNG_SUCCESS);

    case TNG_BYTE_PAIR_SWAP_64: /* Byte pair swapped big endian to little endian */
        *v = ((*v & 0xFF00FF0000000000LL) >> 40) |
             ((*v & 0x00FF00FF00000000LL) >> 24) |
             ((*v & 0x00000000FF00FF00LL) << 24) |
             ((*v & 0x0000000000FF00FFLL) << 40);

        return(TNG_SUCCESS);

    case TNG_BYTE_SWAP_64: /* Byte swapped big endian to little endian */
        *v = ((*v & 0xFFFF000000000000LL) >> 48) |
             ((*v & 0x0000FFFF00000000LL) >> 16) |
             ((*v & 0x00000000FFFF0000LL) << 16) |
             ((*v & 0x000000000000FFFFLL) << 48);

        return(TNG_SUCCESS);

    case TNG_BIG_ENDIAN_64: /* Byte order is reversed. */
        *v = ((*v & 0xFF00000000000000LL) >> 56) | /* Move first byte to end */
             ((*v & 0x00FF000000000000LL) >> 40) | /* Move 2nd byte to pos 7 */
             ((*v & 0x0000FF0000000000LL) >> 24) | /* Move 3rd byte to pos 6 */
             ((*v & 0x000000FF00000000LL) >> 8 ) | /* Move 4th byte to pos 5 */
             ((*v & 0x00000000FF000000LL) << 8 ) | /* Move 5th byte to pos 4 */
             ((*v & 0x0000000000FF0000LL) << 24) | /* Move 6th byte to pos 3 */
             ((*v & 0x000000000000FF00LL) << 40) | /* Move 7th byte to pos 2 */
             ((*v & 0x00000000000000FFLL) << 56);  /* Move last byte to first */

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
static tng_function_status tng_block_md5_hash_generate(tng_gen_block_t block)
{
    md5_state_t md5_state;

    md5_init(&md5_state);
    md5_append(&md5_state, (md5_byte_t *)block->block_contents,
               (int)block->block_contents_size);
    md5_finish(&md5_state, (md5_byte_t *)block->md5_hash);

    return(TNG_SUCCESS);
}

/** Compare the current block md5 hash (e.g. read from file) with the md5 hash
 * calculated from the current contents.
 * If the current md5 hash is not set skip the comparison.
 * @param block is a general block container.
 * @param results If the hashes match results is set to TNG_TRUE, otherwise it is
 * set to TNG_FALSE. If the hash was not set results is set to TNG_TRUE.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (1) if the hash was not
 * set.
 */
static tng_function_status tng_md5_hash_match_verify(tng_gen_block_t block,
                                                     tng_bool *results)
{
    md5_state_t md5_state;
    char hash[TNG_MD5_HASH_LEN];

    TNG_ASSERT(block->block_contents_size > 0, "The block contents size must be > 0");

    *results = TNG_TRUE;
    if(strncmp(block->md5_hash, "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", 16) == 0)
    {
        return(TNG_FAILURE);
    }
    md5_init(&md5_state);
    md5_append(&md5_state, (md5_byte_t *)block->block_contents,
               (int)block->block_contents_size);
    md5_finish(&md5_state, (md5_byte_t *)hash);

    if(strncmp(block->md5_hash, hash, 16) != 0)
    {
        *results = TNG_FALSE;
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
            fprintf(stderr, "TNG library: No file specified for reading. %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        tng_data->input_file = fopen(tng_data->input_file_path, "rb");
        if(!tng_data->input_file)
        {
            fprintf(stderr, "TNG library: Cannot open file %s. %s: %d\n",
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
            fprintf(stderr, "TNG library: No file specified for writing. %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }

        tng_data->output_file = fopen(tng_data->output_file_path, "wb+");

        if(!tng_data->output_file)
        {
            fprintf(stderr, "TNG library: Cannot open file %s. %s: %d\n",
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
    tng_gen_block_t block;

    *block_p = malloc(sizeof(struct tng_gen_block));
    if(!*block_p)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%lu bytes). %s: %d\n",
               sizeof(struct tng_gen_block), __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    block = *block_p;

    block->id = -1;
    /* Reset the md5_hash */
    memcpy(block->md5_hash, "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", TNG_MD5_HASH_LEN);
    block->name = 0;
    block->block_version = TNG_API_VERSION;
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

/*     fprintf(stderr, "TNG library: Destroying block\n"); */
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
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE(1) if a minor
 * error has occured (not able to read the header size, thus skipping
 * the block) or TNG_CRITICAL (2) if a major error has occured.
 */
static tng_function_status tng_block_header_read
                (tng_trajectory_t tng_data, tng_gen_block_t block)
{
    int len, offset = 0;

    TNG_ASSERT(block != 0, "TNG library: Trying to read to uninitialized block (NULL pointer).");

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    /* First read the header size to be able to read the whole header. */
    if(fread(&block->header_contents_size, sizeof(block->header_contents_size),
        1, tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Cannot read header size. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(block->header_contents_size == 0)
    {
        block->id = -1;
        return(TNG_FAILURE);
    }

    /* If this was the size of the general info block check the endianness */
    if(ftell(tng_data->input_file) < 9)
    {
        /* File is little endian */
        if ( *((const char*)&block->header_contents_size) != 0x00 &&
             *((const char*)(&block->header_contents_size) + 7) == 0x00)
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    /* Move the reading position to the beginning of the header. */
    fseek(tng_data->input_file, -(long)sizeof(block->header_contents_size),
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
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->header_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* Read the whole header into header_contents. This way it can be saved
     * even if it cannot be interpreted
     * for one reason or another. */
    if(fread(block->header_contents, block->header_contents_size, 1,
        tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Cannot read header. %s: %d\n", __FILE__, __LINE__);
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    offset += sizeof(block->id);

    memcpy(block->md5_hash, block->header_contents+offset, TNG_MD5_HASH_LEN);
    offset += TNG_MD5_HASH_LEN;

    if(block->name && strcmp(block->name, block->header_contents+offset) != 0)
    {
        free(block->name);
        block->name = 0;
    }
    len = tng_min_i((int)strlen(block->header_contents+offset) + 1, TNG_MAX_STR_LEN);
    if(!block->name)
    {
        block->name = malloc(len);
        if(!block->name)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len,
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    return(TNG_SUCCESS);
}

/** Write a whole block, both header and contents, regardless of it type
 * @param tng_data is a trajectory data container.
 * @param block is a general block container.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
/* Disabled until it is used.*/
/*
// static tng_function_status tng_block_verbatim_write(tng_trajectory_t tng_data,
//                                                     tng_gen_block_t block)
// {
//     if(!block->header_contents)
//     {
//         fprintf(stderr, "TNG library: No contents to write. %s: %d\n", __FILE__, __LINE__);
//         return(TNG_FAILURE);
//     }
//     if(fwrite(block->header_contents, block->header_contents_size, 1,
//                 tng_data->output_file) != 1)
//     {
//         fprintf(stderr, "TNG library: Could not write all header data. %s: %d\n",
//                 __FILE__, __LINE__);
//         return(TNG_CRITICAL);
//     }
//
//     if(!block->block_contents)
//     {
//         fprintf(stderr, "TNG library: No block data to write. %s: %d\n",
//                 __FILE__, __LINE__);
//         return(TNG_FAILURE);
//     }
//     if(fwrite(block->block_contents, block->block_contents_size, 1,
//                 tng_data->output_file) != 1)
//     {
//         fprintf(stderr, "TNG library: Could not write all block data. %s: %d\n",
//                 __FILE__, __LINE__);
//         return(TNG_CRITICAL);
//     }
//     return(TNG_SUCCESS);
// }
*/

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
    if(!block->block_contents)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    fseek(tng_data->output_file, (long)contents_start_pos, SEEK_SET);
    if(fread(block->block_contents, block->block_contents_size, 1,
            tng_data->output_file) == 0)
    {
        fprintf(stderr, "TNG library: Cannot read block. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    tng_block_md5_hash_generate(block);

    fseek(tng_data->output_file, (long)header_start_pos + 3 * sizeof(int64_t),
          SEEK_SET);
    fwrite(block->md5_hash, TNG_MD5_HASH_LEN, 1, tng_data->output_file);

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
                (tng_trajectory_t tng_data, const char hash_mode)
{
    tng_gen_block_t block;
    FILE *temp = tng_data->input_file;
    int64_t output_file_pos, pos, contents_start_pos;

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    tng_data->input_file = tng_data->output_file;

    tng_block_init(&block);

    output_file_pos = ftell(tng_data->output_file);
    fseek(tng_data->output_file, 0, SEEK_SET);

    if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot read general info header. %s: %d\n",
               __FILE__, __LINE__);
        tng_data->input_file = temp;
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    contents_start_pos = ftell(tng_data->output_file);

    fseek(tng_data->output_file, (long)block->block_contents_size - 5 *
          sizeof(int64_t), SEEK_CUR);

    tng_data->input_file = temp;

    pos = tng_data->first_trajectory_frame_set_output_file_pos;

    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                    &pos)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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

    fseek(tng_data->output_file, (long)output_file_pos, SEEK_SET);

    return(TNG_SUCCESS);
}

/** Update the frame set pointers in the current frame set block, already
 * written to disk. It also updates the pointers of the blocks pointing to
 * the current frame set block.
 * @param tng_data is a trajectory data container.
 * @param hash_mode specifies whether to update the block md5 hash when
 * updating the pointers.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_frame_set_pointers_update
                (tng_trajectory_t tng_data, const char hash_mode)
{
    tng_gen_block_t block;
    tng_trajectory_frame_set_t frame_set;
    FILE *temp = tng_data->input_file;
    int64_t pos, output_file_pos, contents_start_pos;

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    tng_block_init(&block);
    output_file_pos = ftell(tng_data->output_file);

    tng_data->input_file = tng_data->output_file;

    frame_set = &tng_data->current_trajectory_frame_set;

    pos = tng_data->current_trajectory_frame_set_output_file_pos;

    /* Update next frame set */
    if(frame_set->next_frame_set_file_pos > 0)
    {
        fseek(tng_data->output_file, (long)frame_set->next_frame_set_file_pos,
              SEEK_SET);

        if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot read frame header. %s: %d\n",
                __FILE__, __LINE__);
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        contents_start_pos = ftell(tng_data->output_file);

        fseek(tng_data->output_file, (long)block->block_contents_size - (5 *
            sizeof(int64_t) + 2 * sizeof(double)), SEEK_CUR);

        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                        &pos)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            tng_md5_hash_update(tng_data, block, frame_set->next_frame_set_file_pos,
                                contents_start_pos);
        }
        fseek(tng_data->output_file, (long)output_file_pos, SEEK_SET);
    }
    /* Update previous frame set */
    if(frame_set->prev_frame_set_file_pos > 0)
    {
        fseek(tng_data->output_file, (long)frame_set->prev_frame_set_file_pos,
              SEEK_SET);

        if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot read frame header. %s: %d\n",
                __FILE__, __LINE__);
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        contents_start_pos = ftell(tng_data->output_file);

        fseek(tng_data->output_file, (long)block->block_contents_size - (6 *
            sizeof(int64_t) + 2 * sizeof(double)), SEEK_CUR);

        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                        &pos)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            tng_md5_hash_update(tng_data, block, frame_set->prev_frame_set_file_pos,
                                contents_start_pos);
        }
        fseek(tng_data->output_file, (long)output_file_pos, SEEK_SET);
    }

    /* Update the frame set one medium stride step after */
    if(frame_set->medium_stride_next_frame_set_file_pos > 0)
    {
        fseek(tng_data->output_file,
              (long)frame_set->medium_stride_next_frame_set_file_pos,
              SEEK_SET);

        if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot read frame set header. %s: %d\n",
                __FILE__, __LINE__);
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        contents_start_pos = ftell(tng_data->output_file);

        fseek(tng_data->output_file, (long)block->block_contents_size - (3 *
            sizeof(int64_t) + 2 * sizeof(double)), SEEK_CUR);

        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                        &pos)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                                frame_set->medium_stride_next_frame_set_file_pos,
                                contents_start_pos);
        }
    }
    /* Update the frame set one medium stride step before */
    if(frame_set->medium_stride_prev_frame_set_file_pos > 0)
    {
        fseek(tng_data->output_file,
              (long)frame_set->medium_stride_prev_frame_set_file_pos,
              SEEK_SET);

        if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot read frame set header. %s: %d\n",
                __FILE__, __LINE__);
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        contents_start_pos = ftell(tng_data->output_file);

        fseek(tng_data->output_file, (long)block->block_contents_size - (4 *
            sizeof(int64_t) + 2 * sizeof(double)), SEEK_CUR);

        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                        &pos)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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

    /* Update the frame set one long stride step after */
    if(frame_set->long_stride_next_frame_set_file_pos > 0)
    {
        fseek(tng_data->output_file,
              (long)frame_set->long_stride_next_frame_set_file_pos,
              SEEK_SET);

        if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot read frame set header. %s: %d\n",
                __FILE__, __LINE__);
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        contents_start_pos = ftell(tng_data->output_file);

        fseek(tng_data->output_file, (long)block->block_contents_size - (1 *
            sizeof(int64_t) + 2 * sizeof(double)), SEEK_CUR);

        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                        &pos)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                                frame_set->long_stride_next_frame_set_file_pos,
                                contents_start_pos);
        }
    }
    /* Update the frame set one long stride step before */
    if(frame_set->long_stride_prev_frame_set_file_pos > 0)
    {
        fseek(tng_data->output_file,
              (long)frame_set->long_stride_prev_frame_set_file_pos,
              SEEK_SET);

        if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot read frame set header. %s: %d\n",
                __FILE__, __LINE__);
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        contents_start_pos = ftell(tng_data->output_file);

        fseek(tng_data->output_file, (long)block->block_contents_size - (2 *
            sizeof(int64_t) + 2 * sizeof(double)), SEEK_CUR);

        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                        &pos)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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

    fseek(tng_data->output_file, (long)output_file_pos, SEEK_SET);

    tng_data->input_file = temp;

    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

static tng_function_status tng_reread_frame_set_at_file_pos
                (tng_trajectory_t tng_data,
                 const int64_t pos)
{
    tng_gen_block_t block;
    tng_function_status stat;

    tng_block_init(&block);

    fseek(tng_data->input_file, pos, SEEK_SET);
    if(pos > 0)
    {
        stat = tng_block_header_read(tng_data, block);
        if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
        {
            fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n", pos,
                    __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_FAILURE);
        }

        if(tng_block_read_next(tng_data, block,
                               TNG_SKIP_HASH) != TNG_SUCCESS)
        {
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
    }

    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

static tng_function_status tng_file_pos_of_subsequent_trajectory_block_get
                (tng_trajectory_t tng_data,
                 int64_t *pos)
{
    int64_t orig_pos, curr_frame_set_pos;
    tng_gen_block_t block;
    tng_function_status stat;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;

    orig_pos = ftell(tng_data->input_file);
    curr_frame_set_pos = tng_data->current_trajectory_frame_set_input_file_pos;

    *pos = tng_data->first_trajectory_frame_set_input_file_pos;

    if(*pos <= 0)
    {
        return(TNG_SUCCESS);
    }

    fseek(tng_data->input_file, *pos, SEEK_SET);

    tng_block_init(&block);
    /* Read block headers first to see that a frame set block is found. */
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n", *pos,
                __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_FAILURE);
    }

    if(tng_block_read_next(tng_data, block,
                           TNG_SKIP_HASH) != TNG_SUCCESS)
    {
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    /* Read all frame set blocks (not the blocks between them) */
    while(frame_set->next_frame_set_file_pos > 0)
    {
        fseek(tng_data->input_file, frame_set->next_frame_set_file_pos, SEEK_SET);
        stat = tng_block_header_read(tng_data, block);
        if(stat == TNG_CRITICAL)
        {
            fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n", *pos,
                    __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
        if(stat != TNG_SUCCESS || block->id != TNG_TRAJECTORY_FRAME_SET)
        {
            return(TNG_FAILURE);
        }

        stat = tng_block_read_next(tng_data, block, TNG_SKIP_HASH);
        if(stat != TNG_SUCCESS)
        {
            tng_block_destroy(&block);
            return(stat);
        }
        /* Update *pos if this is the earliest frame set so far (after orig_pos) */
        if(tng_data->current_trajectory_frame_set_input_file_pos < *pos &&
           tng_data->current_trajectory_frame_set_input_file_pos > orig_pos)
        {
            *pos = tng_data->current_trajectory_frame_set_input_file_pos;
        }
    }

    /* Re-read the frame set that used to be the current one */
    tng_reread_frame_set_at_file_pos(tng_data, curr_frame_set_pos);

    fseek(tng_data->input_file, orig_pos, SEEK_SET);

    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

static tng_function_status tng_frame_set_complete_migrate
                (tng_trajectory_t tng_data,
                 int64_t block_start_pos,
                 int64_t block_len,
                 int64_t new_pos)
{
    int64_t i;
    tng_bool updated = TNG_FALSE;

    char *contents;

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    fseek(tng_data->input_file, block_start_pos, SEEK_SET);

    contents = malloc(block_len);
    if(!contents)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                block_len, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(fread(contents, block_len, 1, tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Cannot read data from file when migrating data. %s: %d\n",
               __FILE__, __LINE__);
        free(contents);
        return(TNG_CRITICAL);
    }
    fseek(tng_data->output_file, new_pos, SEEK_SET);

    if(fwrite(contents, block_len, 1, tng_data->output_file) != 1)
    {
        fprintf(stderr, "TNG library: Could not write data to file when migrating data. %s: %d\n",
                __FILE__, __LINE__);
        free(contents);
        return(TNG_CRITICAL);
    }

    tng_data->current_trajectory_frame_set_output_file_pos = new_pos;

    tng_frame_set_pointers_update(tng_data, TNG_USE_HASH);

    /* Update the general info block if needed */
    if(block_start_pos == tng_data->first_trajectory_frame_set_output_file_pos)
    {
        tng_data->first_trajectory_frame_set_output_file_pos = new_pos;
        updated = TNG_TRUE;
    }
    if(block_start_pos == tng_data->last_trajectory_frame_set_output_file_pos)
    {
        tng_data->last_trajectory_frame_set_output_file_pos = new_pos;
        updated = TNG_TRUE;
    }
    if(updated)
    {
        tng_header_pointers_update(tng_data, TNG_USE_HASH);
    }

    /* Fill the block with NULL to avoid confusion. */
    for(i = 0; i < block_len; i++)
    {
        contents[i] = '\0';
    }
    fseek(tng_data->output_file, block_start_pos, SEEK_SET);

    /* FIXME: casting block_len to size_t is dangerous */
    fwrite(contents, 1, block_len, tng_data->output_file);

    free(contents);

    return(TNG_SUCCESS);
}

static tng_function_status tng_length_of_current_frame_set_contents_get
                (tng_trajectory_t tng_data,
                 int64_t *len)
{
    int64_t orig_pos, pos, curr_frame_set_pos;
    tng_gen_block_t block;
    tng_function_status stat;

    orig_pos = ftell(tng_data->input_file);
    curr_frame_set_pos = pos = tng_data->current_trajectory_frame_set_input_file_pos;

    *len = 0;

    fseek(tng_data->input_file, curr_frame_set_pos, SEEK_SET);

    tng_block_init(&block);
    /* Read block headers first to see that a frame set block is found. */
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
                curr_frame_set_pos, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_FAILURE);
    }

    /* Read the headers of all blocks in the frame set (not the actual contents of them) */
    while(stat == TNG_SUCCESS)
    {
        fseek(tng_data->input_file, block->block_contents_size, SEEK_CUR);
        *len += block->header_contents_size + block->block_contents_size;
        pos += block->header_contents_size + block->block_contents_size;
        if(pos >= tng_data->input_file_len)
        {
            break;
        }
        stat = tng_block_header_read(tng_data, block);
        if(block->id == TNG_TRAJECTORY_FRAME_SET)
        {
            break;
        }
    }

    /* Re-read the frame set that used to be the current one */
    tng_reread_frame_set_at_file_pos(tng_data, curr_frame_set_pos);

    fseek(tng_data->input_file, orig_pos, SEEK_SET);

    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

/** Migrate blocks in the file to make room for new data in a block. This
 * is required e.g. when adding data to a block or extending strings in a
 * block.
 * @param tng_data is a trajectory data container.
 * @param start_pos is the position from which to start moving data, usually
 * the byte after the end of the block to which data was added.
 * @param offset is the number of bytes that were inserted.
 * @details Trajectory blocks (frame sets and their related blocks) are moved
 * to the end of the file (if needed) in order to make room for non-trajectory
 * data.
 */
static tng_function_status tng_migrate_data_in_file
                (tng_trajectory_t tng_data,
                 int64_t start_pos,
                 int64_t offset)
{
    int64_t traj_start_pos, empty_space, orig_file_pos, frame_set_length;
    tng_gen_block_t block;
    tng_function_status stat;
    FILE *temp;

    if(offset <= 0)
    {
        return(TNG_SUCCESS);
    }

    temp = tng_data->input_file;

    stat = tng_file_pos_of_subsequent_trajectory_block_get(tng_data, &traj_start_pos);
    if(stat != TNG_SUCCESS)
    {
        tng_data->input_file = temp;
        return(stat);
    }

    tng_data->current_trajectory_frame_set_input_file_pos = traj_start_pos;

    empty_space = traj_start_pos - (start_pos - 1);

    if(empty_space >= offset)
    {
        return(TNG_SUCCESS);
    }

    orig_file_pos = ftell(tng_data->input_file);
    tng_block_init(&block);

    while(empty_space < offset)
    {
        fseek(tng_data->input_file, traj_start_pos, SEEK_SET);
        stat = tng_block_header_read(tng_data, block);
        if(stat == TNG_CRITICAL)
        {
            fprintf(stderr, "TNG library: Cannot read block header. %s: %d\n",
                    __FILE__, __LINE__);
            tng_block_destroy(&block);
            tng_data->input_file = temp;
            return(TNG_CRITICAL);
        }
        if(stat != TNG_SUCCESS || block->id != TNG_TRAJECTORY_FRAME_SET)
        {
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(TNG_FAILURE);
        }
        stat = tng_length_of_current_frame_set_contents_get(tng_data, &frame_set_length);
        if(stat != TNG_SUCCESS)
        {
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(stat);
        }
        stat = tng_frame_set_complete_migrate(tng_data, traj_start_pos,
                                              frame_set_length, tng_data->input_file_len);
        if(stat != TNG_SUCCESS)
        {
            tng_data->input_file = temp;
            tng_block_destroy(&block);
            return(stat);
        }

        empty_space += frame_set_length;
    }
    fseek(tng_data->input_file, orig_file_pos, SEEK_SET);
    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

static tng_function_status tng_block_header_len_calculate
                (const tng_trajectory_t tng_data,
                 tng_gen_block_t block,
                 int64_t *len)
{
    int name_len;
    (void)tng_data;

    /* If the string is unallocated allocate memory for just string
     * termination */
    if(!block->name)
    {
        block->name = malloc(1);
        if(!block->name)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        block->name[0] = 0;
    }

    name_len = tng_min_i((int)strlen(block->name) + 1, TNG_MAX_STR_LEN);

    /* Calculate the size of the header to write */
    *len = sizeof(block->header_contents_size) +
                  sizeof(block->block_contents_size) +
                  sizeof(block->id) +
                  sizeof(block->block_version) +
                  TNG_MD5_HASH_LEN +
                  name_len;

    return (TNG_SUCCESS);
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
                 const char hash_mode)
{
    int name_len, offset = 0;

    TNG_ASSERT(block != 0, "TNG library: Trying to write uninitialized block (NULL pointer).");

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(tng_block_header_len_calculate(tng_data, block, &block->header_contents_size) !=
        TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot calculate length of block header. %s: %d\n",
                __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(hash_mode == TNG_USE_HASH)
    {
        tng_block_md5_hash_generate(block);
    }

    if(block->header_contents)
    {
        free(block->header_contents);
    }

    block->header_contents = malloc(block->header_contents_size);
    if(!block->header_contents)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->header_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    name_len = tng_min_i((int)strlen(block->name) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(block->id);

    memcpy(block->header_contents+offset, block->md5_hash, TNG_MD5_HASH_LEN);
    offset += TNG_MD5_HASH_LEN;

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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(fwrite(block->header_contents, block->header_contents_size,
       1, tng_data->output_file) != 1)
    {
        fprintf(stderr, "TNG library: Could not write all header data. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }
    return(TNG_SUCCESS);
}

static tng_function_status tng_general_info_block_len_calculate
                (tng_trajectory_t tng_data,
                 int64_t *len)
{
    int first_program_name_len, first_user_name_len;
    int first_computer_name_len, first_pgp_signature_len;
    int last_program_name_len, last_user_name_len;
    int last_computer_name_len, last_pgp_signature_len;
    int forcefield_name_len;

    /* If the strings are unallocated allocate memory for just string
     * termination */
    if(!tng_data->first_program_name)
    {
        tng_data->first_program_name = malloc(1);
        if(!tng_data->first_program_name)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        tng_data->forcefield_name[0] = 0;
    }

    first_program_name_len = tng_min_i((int)strlen(tng_data->first_program_name) + 1,
                           TNG_MAX_STR_LEN);
    last_program_name_len = tng_min_i((int)strlen(tng_data->last_program_name) + 1,
                           TNG_MAX_STR_LEN);
    first_user_name_len = tng_min_i((int)strlen(tng_data->first_user_name) + 1,
                        TNG_MAX_STR_LEN);
    last_user_name_len = tng_min_i((int)strlen(tng_data->last_user_name) + 1,
                        TNG_MAX_STR_LEN);
    first_computer_name_len = tng_min_i((int)strlen(tng_data->first_computer_name) + 1,
                            TNG_MAX_STR_LEN);
    last_computer_name_len = tng_min_i((int)strlen(tng_data->last_computer_name) + 1,
                            TNG_MAX_STR_LEN);
    first_pgp_signature_len = tng_min_i((int)strlen(tng_data->first_pgp_signature) + 1,
                            TNG_MAX_STR_LEN);
    last_pgp_signature_len = tng_min_i((int)strlen(tng_data->last_pgp_signature) + 1,
                            TNG_MAX_STR_LEN);
    forcefield_name_len = tng_min_i((int)strlen(tng_data->forcefield_name) + 1,
                              TNG_MAX_STR_LEN);

    *len = sizeof(tng_data->time) +
                  sizeof(tng_data->var_num_atoms_flag) +
                  sizeof(tng_data->frame_set_n_frames) +
                  sizeof(tng_data->first_trajectory_frame_set_input_file_pos) +
                  sizeof(tng_data->last_trajectory_frame_set_input_file_pos) +
                  sizeof(tng_data->medium_stride_length) +
                  sizeof(tng_data->long_stride_length) +
                  sizeof(tng_data->distance_unit_exponential) +
                  first_program_name_len +
                  last_program_name_len +
                  first_user_name_len +
                  last_user_name_len +
                  first_computer_name_len +
                  last_computer_name_len +
                  first_pgp_signature_len +
                  last_pgp_signature_len +
                  forcefield_name_len;

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
                 const char hash_mode)
{
    int len, offset = 0;
    tng_bool same_hash;

    void *temp;

    TNG_ASSERT(block != 0, "TNG library: Trying to read data to an uninitialized block (NULL pointer)");

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    temp = realloc(block->block_contents, block->block_contents_size);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        free(block->block_contents);
        block->block_contents = 0;
        return(TNG_CRITICAL);
    }
    block->block_contents = temp;

    /* Read the whole block into block_contents to be able to write it to disk
     * even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
             tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Cannot read block. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */

    if(hash_mode == TNG_USE_HASH)
    {
        tng_md5_hash_match_verify(block, &same_hash);
        if(same_hash != TNG_TRUE)
        {
            fprintf(stderr, "TNG library: General info block contents corrupt. Hashes do not match. "
                "%s: %d\n",
                __FILE__, __LINE__);
    /*         return(TNG_FAILURE); */
        }
    }

    len = tng_min_i((int)strlen(block->block_contents) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->first_program_name, len);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->first_program_name);
        tng_data->first_program_name = 0;
        return(TNG_CRITICAL);
    }
    tng_data->first_program_name = temp;
    strncpy(tng_data->first_program_name, block->block_contents, len);
    offset += len;

    len = tng_min_i((int)strlen(block->block_contents + offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->last_program_name, len);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->last_program_name);
        tng_data->last_program_name = 0;
        return(TNG_CRITICAL);
    }
    tng_data->last_program_name = temp;
    strncpy(tng_data->last_program_name, block->block_contents + offset, len);
    offset += len;

    len = tng_min_i((int)strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->first_user_name, len);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->first_user_name);
        tng_data->first_user_name = 0;
        return(TNG_CRITICAL);
    }
    tng_data->first_user_name = temp;
    strncpy(tng_data->first_user_name, block->block_contents+offset, len);
    offset += len;

    len = tng_min_i((int)strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->last_user_name, len);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->last_user_name);
        tng_data->last_user_name = 0;
        return(TNG_CRITICAL);
    }
    tng_data->last_user_name = temp;
    strncpy(tng_data->last_user_name, block->block_contents+offset, len);
    offset += len;

    len = tng_min_i((int)strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->first_computer_name, len);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->first_computer_name);
        tng_data->first_computer_name = 0;
        return(TNG_CRITICAL);
    }
    tng_data->first_computer_name = temp;
    strncpy(tng_data->first_computer_name, block->block_contents+offset, len);
    offset += len;

    len = tng_min_i((int)strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->last_computer_name, len);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->last_computer_name);
        tng_data->last_computer_name = 0;
        return(TNG_CRITICAL);
    }
    tng_data->last_computer_name = temp;
    strncpy(tng_data->last_computer_name, block->block_contents+offset, len);
    offset += len;

    len = tng_min_i((int)strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->first_pgp_signature, len);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->first_pgp_signature);
        tng_data->first_pgp_signature = 0;
        return(TNG_CRITICAL);
    }
    tng_data->first_pgp_signature = temp;
    strncpy(tng_data->first_pgp_signature, block->block_contents+offset, len);
    offset += len;

    len = tng_min_i((int)strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->last_pgp_signature, len);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->last_pgp_signature);
        tng_data->last_pgp_signature = 0;
        return(TNG_CRITICAL);
    }
    tng_data->last_pgp_signature = temp;
    strncpy(tng_data->last_pgp_signature, block->block_contents+offset, len);
    offset += len;

    len = tng_min_i((int)strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->forcefield_name, len);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->forcefield_name);
        tng_data->forcefield_name = 0;
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->long_stride_length);

    if(block->block_version >= 3)
    {
        memcpy(&tng_data->distance_unit_exponential, block->block_contents+offset,
            sizeof(tng_data->distance_unit_exponential));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                          &tng_data->distance_unit_exponential)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
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
                 const char hash_mode)
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

    fseek(tng_data->output_file, 0, SEEK_SET);

    tng_block_init(&block);

    name_len = (int)strlen("GENERAL INFO");

    block->name = malloc(name_len + 1);
    if(!block->name)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
                name_len+1, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    strcpy(block->name, "GENERAL INFO");
    block->id = TNG_GENERAL_INFO;

    if(tng_general_info_block_len_calculate(tng_data, &block->block_contents_size) !=
        TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot calculate length of general info block. %s: %d\n",
                __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    first_program_name_len = tng_min_i((int)strlen(tng_data->first_program_name) + 1,
                           TNG_MAX_STR_LEN);
    last_program_name_len = tng_min_i((int)strlen(tng_data->last_program_name) + 1,
                           TNG_MAX_STR_LEN);
    first_user_name_len = tng_min_i((int)strlen(tng_data->first_user_name) + 1,
                        TNG_MAX_STR_LEN);
    last_user_name_len = tng_min_i((int)strlen(tng_data->last_user_name) + 1,
                        TNG_MAX_STR_LEN);
    first_computer_name_len = tng_min_i((int)strlen(tng_data->first_computer_name) + 1,
                            TNG_MAX_STR_LEN);
    last_computer_name_len = tng_min_i((int)strlen(tng_data->last_computer_name) + 1,
                            TNG_MAX_STR_LEN);
    first_pgp_signature_len = tng_min_i((int)strlen(tng_data->first_pgp_signature) + 1,
                            TNG_MAX_STR_LEN);
    last_pgp_signature_len = tng_min_i((int)strlen(tng_data->last_pgp_signature) + 1,
                            TNG_MAX_STR_LEN);
    forcefield_name_len = tng_min_i((int)strlen(tng_data->forcefield_name) + 1,
                              TNG_MAX_STR_LEN);

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->frame_set_n_frames);

    memcpy(block->block_contents+offset,
           &tng_data->first_trajectory_frame_set_output_file_pos,
           sizeof(tng_data->first_trajectory_frame_set_output_file_pos));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->first_trajectory_frame_set_output_file_pos);

    memcpy(block->block_contents+offset,
           &tng_data->last_trajectory_frame_set_output_file_pos,
           sizeof(tng_data->last_trajectory_frame_set_output_file_pos));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->last_trajectory_frame_set_output_file_pos);

    memcpy(block->block_contents+offset, &tng_data->medium_stride_length,
           sizeof(tng_data->medium_stride_length));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->long_stride_length);

    memcpy(block->block_contents+offset, &tng_data->distance_unit_exponential,
           sizeof(tng_data->distance_unit_exponential));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(tng_block_header_write(tng_data, block, hash_mode) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot write header of file %s. %s: %d\n",
               tng_data->output_file_path, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
        tng_data->output_file) != 1)
    {
        fprintf(stderr, "TNG library: Could not write all block data. %s: %d\n", __FILE__, __LINE__);
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

    TNG_ASSERT(offset != 0, "TNG library: offset must not be a NULL pointer.");

    memcpy(&chain->id, block->block_contents+*offset,
            sizeof(chain->id));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &chain->id)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(chain->id);

    len = tng_min_i((int)strlen(block->block_contents+*offset) + 1,
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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

    TNG_ASSERT(offset != 0, "TNG library: offset must not be a NULL pointer.");

    memcpy(block->block_contents+*offset, &chain->id, sizeof(chain->id));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                    (int64_t *)block->header_contents+*offset)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(chain->id);

    len = tng_min_i((int)strlen(chain->name) + 1, TNG_MAX_STR_LEN);
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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

    TNG_ASSERT(offset != 0, "TNG library: offset must not be a NULL pointer.");

    memcpy(&residue->id, block->block_contents+*offset,
        sizeof(residue->id));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   &residue->id)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(residue->id);

    len = tng_min_i((int)strlen(block->block_contents+*offset) + 1,
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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

    TNG_ASSERT(offset != 0, "TNG library: offset must not be a NULL pointer.");

    memcpy(block->block_contents+*offset, &residue->id, sizeof(residue->id));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                    (int64_t *)block->header_contents+*offset)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(residue->id);

    len = tng_min_i((int)strlen(residue->name) + 1, TNG_MAX_STR_LEN);
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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

    TNG_ASSERT(offset != 0, "TNG library: offset must not be a NULL pointer.");

    memcpy(&atom->id, block->block_contents+*offset,
        sizeof(atom->id));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                    &atom->id)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(atom->id);

    len = tng_min_i((int)strlen(block->block_contents+*offset) + 1,
            TNG_MAX_STR_LEN);
    atom->name = malloc(len);
    strncpy(atom->name,
            block->block_contents+*offset, len);
    *offset += len;

    len = tng_min_i((int)strlen(block->block_contents+*offset) + 1,
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

    TNG_ASSERT(offset != 0, "TNG library: offset must not be a NULL pointer.");

    memcpy(block->block_contents+*offset, &atom->id,
            sizeof(atom->id));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                    (int64_t *)block->header_contents+*offset)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(atom->id);

    len = tng_min_i((int)strlen(atom->name) + 1, TNG_MAX_STR_LEN);
    strncpy(block->block_contents + *offset, atom->name, len);
    *offset += len;

    len = tng_min_i((int)strlen(atom->atom_type) + 1, TNG_MAX_STR_LEN);
    strncpy(block->block_contents + *offset, atom->atom_type, len);
    *offset += len;

    return(TNG_SUCCESS);
}

static tng_function_status tng_molecules_block_len_calculate
                (const tng_trajectory_t tng_data,
                 int64_t *len)
{
    int64_t i, j;
    tng_molecule_t molecule;
    tng_chain_t chain;
    tng_residue_t residue;
    tng_atom_t atom;
    tng_bond_t bond;

    *len = 0;

    for(i = 0; i < tng_data->n_molecules; i++)
    {
        molecule = &tng_data->molecules[i];
        if(!molecule->name)
        {
            molecule->name = malloc(1);
            if(!molecule->name)
            {
                fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
                       __FILE__, __LINE__);
                return(TNG_CRITICAL);
            }
            molecule->name[0] = 0;
        }
        *len += tng_min_i((int)strlen(molecule->name) + 1, TNG_MAX_STR_LEN);

        chain = molecule->chains;
        for(j = 0; j < molecule->n_chains; j++)
        {
            *len += sizeof(chain->id);

            if(!chain->name)
            {
                chain->name = malloc(1);
                if(!chain->name)
                {
                    fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
                           __FILE__, __LINE__);
                    return(TNG_CRITICAL);
                }
                chain->name[0] = 0;
            }
            *len += tng_min_i((int)strlen(chain->name) + 1, TNG_MAX_STR_LEN);

            *len += sizeof(chain->n_residues);

            chain++;
        }

        residue = molecule->residues;
        for(j = 0; j < molecule->n_residues; j++)
        {
            *len += sizeof(residue->id);

            if(!residue->name)
            {
                residue->name = malloc(1);
                if(!residue->name)
                {
                    fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
                           __FILE__, __LINE__);
                    return(TNG_CRITICAL);
                }
                residue->name[0] = 0;
            }
            *len += tng_min_i((int)strlen(residue->name) + 1, TNG_MAX_STR_LEN);

            *len += sizeof(residue->n_atoms);

            residue++;
        }

        atom = molecule->atoms;
        for(j = 0; j < molecule->n_atoms; j++)
        {
            *len += sizeof(atom->id);
            if(!atom->name)
            {
                atom->name = malloc(1);
                if(!atom->name)
                {
                    fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
                           __FILE__, __LINE__);
                    return(TNG_CRITICAL);
                }
                atom->name[0] = 0;
            }
            *len += tng_min_i((int)strlen(atom->name) + 1, TNG_MAX_STR_LEN);

            if(!atom->atom_type)
            {
                atom->atom_type = malloc(1);
                if(!atom->atom_type)
                {
                    fprintf(stderr, "TNG library: Cannot allocate memory (1 byte). %s: %d\n",
                           __FILE__, __LINE__);
                    return(TNG_CRITICAL);
                }
                atom->atom_type[0] = 0;
            }
            *len += tng_min_i((int)strlen(atom->atom_type) + 1, TNG_MAX_STR_LEN);

            atom++;
        }

        for(j = 0; j < molecule->n_bonds; j++)
        {
            *len += sizeof(bond->from_atom_id) + sizeof(bond->to_atom_id);
        }
    }
    *len += sizeof(tng_data->n_molecules) +
            (sizeof(molecule->id) +
            sizeof(molecule->quaternary_str) +
            sizeof(molecule->n_chains) +
            sizeof(molecule->n_residues) +
            sizeof(molecule->n_atoms) +
            sizeof(molecule->n_bonds)) *
            tng_data->n_molecules;

    if(!tng_data->var_num_atoms_flag)
    {
        *len += tng_data->n_molecules * sizeof(int64_t);
    }

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
                 const char hash_mode)
{
    int64_t i, j, k, l;
    int len, offset = 0;
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
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to disk
     * even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
             tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Cannot read block. %s: %d\n", __FILE__, __LINE__);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */

    if(hash_mode == TNG_USE_HASH)
    {
        tng_md5_hash_match_verify(block, &same_hash);
        if(same_hash != TNG_TRUE)
        {
            fprintf(stderr, "TNG library: Molecules block contents corrupt. Hashes do not match. "
                "%s: %d\n",
                __FILE__, __LINE__);
        }
    }

    if(tng_data->molecules)
    {
        for(i=0; i<tng_data->n_molecules; i++)
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->id);

/*         fprintf(stderr, "TNG library: Read id: %"PRId64" offset: %d\n", molecule->id, offset);*/
        len = tng_min_i((int)strlen(block->block_contents+offset) + 1, TNG_MAX_STR_LEN);
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                    fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                    molecule->n_residues * sizeof(struct tng_residue),
                    __FILE__, __LINE__);
                if(molecule->chains)
                {
                    free(molecule->chains);
                    molecule->chains = 0;
                }
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
            fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   molecule->n_atoms * sizeof(struct tng_atom),
                   __FILE__, __LINE__);
            if(molecule->chains)
            {
                free(molecule->chains);
                molecule->chains = 0;
            }
            if(molecule->residues)
            {
                free(molecule->residues);
                molecule->residues = 0;
            }
            return(TNG_CRITICAL);
        }

        atom = molecule->atoms;

        if(molecule->n_chains > 0)
        {
            /* Read the chains of the molecule */
            for(j=0; j<molecule->n_chains; j++)
            {
                chain->molecule = molecule;

                tng_chain_data_read(tng_data, block, chain, &offset);

                chain->residues = molecule->residues;
                residue = chain->residues;

                /* Read the residues of the chain */
                for(k=0; k<chain->n_residues; k++)
                {
                    residue->chain = chain;

                    tng_residue_data_read(tng_data, block, residue, &offset);

                    residue->atoms_offset = atom - molecule->atoms;
                    /* Read the atoms of the residue */
                    for(l=0; l<residue->n_atoms; l++)
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
                for(k=0; k<molecule->n_residues; k++)
                {
                    residue->chain = 0;

                    tng_residue_data_read(tng_data, block, residue, &offset);

                    residue->atoms_offset = atom - molecule->atoms;
                    /* Read the atoms of the residue */
                    for(l=0; l<residue->n_atoms; l++)
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
                for(l=0; l<molecule->n_atoms; l++)
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_bonds);

        if(molecule->n_bonds > 0)
        {
            tng_data->molecules[i].bonds = malloc(molecule->n_bonds *
                                           sizeof(struct tng_bond));
            if(!molecule->bonds)
            {
                fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                       molecule->n_bonds * sizeof(struct tng_bond),
                       __FILE__, __LINE__);
                if(molecule->chains)
                {
                    free(molecule->chains);
                    molecule->chains = 0;
                }
                if(molecule->residues)
                {
                    free(molecule->residues);
                    molecule->residues = 0;
                }
                if(molecule->atoms)
                {
                    free(molecule->atoms);
                    molecule->atoms = 0;
                }
                return(TNG_CRITICAL);
            }

            bond = molecule->bonds;

            for(j=0; j<molecule->n_bonds; j++)
            {
                memcpy(&bond->from_atom_id, block->block_contents+offset,
                    sizeof(bond->from_atom_id));
                if(tng_data->input_endianness_swap_func_64)
                {
                    if(tng_data->input_endianness_swap_func_64(tng_data,
                                                               &bond->from_atom_id)
                        != TNG_SUCCESS)
                    {
                        fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                        fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                __FILE__, __LINE__);
                    }
                }
                offset += sizeof(bond->to_atom_id);

                bond++;
            }
        }
        else
        {
            molecule->bonds = 0;
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
                 const char hash_mode)
{
    int len = 0, name_len, offset = 0;
    int64_t i, j, k, l;
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

    tng_block_init(&block);

    name_len = (int)strlen("MOLECULES");

    block->name = malloc(name_len + 1);
    if(!block->name)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
                name_len+1, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    strcpy(block->name, "MOLECULES");
    block->id = TNG_MOLECULES;

    if(tng_molecules_block_len_calculate(tng_data, &block->block_contents_size) !=
        TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot calculate length of molecules block. %s: %d\n",
                __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(tng_data->n_molecules);

    for(i = 0; i < tng_data->n_molecules; i++)
    {
        molecule = &tng_data->molecules[i];
        memcpy(block->block_contents+offset, &molecule->id,
               sizeof(molecule->id));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                                        (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->id);

/*         fprintf(stderr, "TNG library: Wrote id: %"PRId64" offset: %d\n", molecule->id, offset); */
        len = tng_min_i((int)strlen(molecule->name) + 1, TNG_MAX_STR_LEN);
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                    fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_atoms);

        if(molecule->n_chains > 0)
        {
            chain = molecule->chains;
            for(j = 0; j < molecule->n_chains; j++)
            {
                tng_chain_data_write(tng_data, block, chain, &offset);

                residue = chain->residues;
                for(k = 0; k < chain->n_residues; k++)
                {
                    tng_residue_data_write(tng_data, block, residue, &offset);

                    atom = molecule->atoms + residue->atoms_offset;
                    for(l = 0; l < residue->n_atoms; l++)
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
                residue = molecule->residues;
                for(k = 0; k < molecule->n_residues; k++)
                {
                    tng_residue_data_write(tng_data, block, residue, &offset);

                    atom = molecule->atoms + residue->atoms_offset;
                    for(l = 0; l < residue->n_atoms; l++)
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
                for(l = 0; l < molecule->n_atoms; l++)
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_bonds);

        bond = molecule->bonds;
        for(j = 0; j < molecule->n_bonds; j++)
        {
            memcpy(block->block_contents+offset, &bond->from_atom_id,
                   sizeof(bond->from_atom_id));
            if(tng_data->output_endianness_swap_func_64)
            {
                if(tng_data->output_endianness_swap_func_64(tng_data,
                                            (int64_t *)block->header_contents+offset)
                    != TNG_SUCCESS)
                {
                    fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                    fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }
            offset += sizeof(bond->to_atom_id);

            bond++;
        }
    }

    if(tng_block_header_write(tng_data, block, hash_mode) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot write header of file %s. %s: %d\n",
               tng_data->output_file_path, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
              tng_data->output_file) != 1)
    {
        fprintf(stderr, "TNG library: Could not write all block data. %s: %d\n",
               __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

static tng_function_status tng_frame_set_block_len_calculate
                (const tng_trajectory_t tng_data,
                 int64_t *len)
{
    *len = sizeof(int64_t) * 8;
    *len += sizeof(double) * 2;

    if(tng_data->var_num_atoms_flag)
    {
        *len += sizeof(int64_t) * tng_data->n_molecules;
    }
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
                 const char hash_mode)
{
    long file_pos;
    int offset = 0;
    int64_t i, prev_n_particles;
    tng_bool same_hash;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;

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
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to
     * disk even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
             tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Cannot read block. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */

    file_pos = (int64_t)ftell(tng_data->input_file) -
               (long)(block->block_contents_size + block->header_contents_size);

    if(hash_mode == TNG_USE_HASH)
    {
        tng_md5_hash_match_verify(block, &same_hash);
        if(same_hash != TNG_TRUE)
        {
            fprintf(stderr, "TNG library: Frame set block contents corrupt. File pos %ld Hashes do not match. "
                "%s: %d\n",
                file_pos, __FILE__, __LINE__);
    /*         return(TNG_FAILURE); */
        }
    }

    tng_data->current_trajectory_frame_set_input_file_pos = file_pos;

    tng_frame_set_particle_mapping_free(tng_data);

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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                    fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
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
                    fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->long_stride_prev_frame_set_file_pos);

    if(block->block_version >= 3)
    {
        memcpy(&frame_set->first_frame_time,
            block->block_contents + offset,
            sizeof(frame_set->first_frame_time));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                (int64_t *)&frame_set->first_frame_time)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(frame_set->first_frame_time);

        memcpy(&tng_data->time_per_frame,
            block->block_contents + offset,
            sizeof(tng_data->time_per_frame));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                (int64_t *)&tng_data->time_per_frame)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
    }
    else
    {
        frame_set->first_frame_time = -1;
        tng_data->time_per_frame = -1;
    }

    /* If the output file and the input files are the same the number of
     * frames in the file are the same number as has just been read.
     * This is updated here to later on see if there have been new frames
     * added and thereby the frame set needs to be rewritten. */
    if(tng_data->output_file == tng_data->input_file)
    {
        frame_set->n_written_frames = frame_set->n_frames;
    }

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
                 const char hash_mode)
{
    char *temp_name;
    int64_t i;
    int offset = 0;
    unsigned int name_len;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    name_len = (int)strlen("TRAJECTORY FRAME SET");

    if(!block->name || strlen(block->name) < name_len)
    {
        temp_name = realloc(block->name, name_len + 1);
        if(!temp_name)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
                   name_len+1, __FILE__, __LINE__);
            free(block->name);
            block->name = 0;
            return(TNG_CRITICAL);
        }
        block->name = temp_name;
    }
    strcpy(block->name, "TRAJECTORY FRAME SET");
    block->id = TNG_TRAJECTORY_FRAME_SET;

    if(tng_frame_set_block_len_calculate(tng_data, &block->block_contents_size) !=
        TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot calculate length of frame set block. %s: %d\n",
                __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                    fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->long_stride_prev_frame_set_file_pos);

    memcpy(block->block_contents+offset,
           &frame_set->first_frame_time,
           sizeof(frame_set->first_frame_time));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->first_frame_time);

    memcpy(block->block_contents+offset,
           &tng_data->time_per_frame,
           sizeof(tng_data->time_per_frame));
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
                                      (int64_t *)block->header_contents+offset)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(tng_block_header_write(tng_data, block, hash_mode) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot write header of file %s. %s: %d\n",
               tng_data->output_file_path, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
              tng_data->output_file) != 1)
    {
        fprintf(stderr, "TNG library: Could not write all block data. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    return(TNG_SUCCESS);
}

static tng_function_status tng_trajectory_mapping_block_len_calculate
                (const tng_trajectory_t tng_data,
                 const int64_t n_particles,
                 int64_t *len)
{
    (void)tng_data;
    *len = sizeof(int64_t) * (2 + n_particles);

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
                 const char hash_mode)
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
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to disk
     *  even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
        tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Cannot read block. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */

    if(hash_mode == TNG_USE_HASH)
    {
        tng_md5_hash_match_verify(block, &same_hash);
        if(same_hash != TNG_TRUE)
        {
            fprintf(stderr, "TNG library: Particle mapping block contents corrupt. Hashes do not match. "
                "%s: %d\n",
                __FILE__, __LINE__);
    /*         return(TNG_FAILURE); */
        }
    }

    frame_set->n_mapping_blocks++;
    mappings = realloc(frame_set->mappings,
                       sizeof(struct tng_particle_mapping) *
                       frame_set->n_mapping_blocks);
    if(!mappings)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        free(frame_set->mappings);
        frame_set->mappings = 0;
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(mapping->n_particles);

    mapping->real_particle_numbers = malloc(mapping->n_particles *
                                            sizeof(int64_t));
    if(!mapping->real_particle_numbers)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                 const char hash_mode)
{
    char *temp_name;
    int i, offset = 0;
    unsigned int name_len;
    tng_particle_mapping_t mapping =
    &tng_data->current_trajectory_frame_set.mappings[mapping_block_nr];

    if(mapping_block_nr >=
       tng_data->current_trajectory_frame_set.n_mapping_blocks)
    {
        fprintf(stderr, "TNG library: Mapping block index out of bounds. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    name_len = (int)strlen("PARTICLE MAPPING");

    if(!block->name || strlen(block->name) < name_len)
    {
        temp_name = realloc(block->name, name_len + 1);
        if(!temp_name)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
                   name_len+1, __FILE__, __LINE__);
            free(block->name);
            block->name = 0;
            return(TNG_CRITICAL);
        }
        block->name = temp_name;
    }
    strcpy(block->name, "PARTICLE MAPPING");
    block->id = TNG_PARTICLE_MAPPING;

    if(tng_trajectory_mapping_block_len_calculate(tng_data,
                                                  mapping->n_particles,
                                                  &block->block_contents_size) !=
        TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot calculate length of atom mapping block. %s: %d\n",
                __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
        fprintf(stderr, "TNG library: Cannot write header of file %s. %s: %d\n",
               tng_data->output_file_path, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
              tng_data->output_file) != 1)
    {
        fprintf(stderr, "TNG library: Could not write all block data. %s: %d\n", __FILE__, __LINE__);
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
                 const char block_type_flag)
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
            fprintf(stderr, "TNG library: Cannot allocate memory (%lu bytes). %s: %d\n",
                sizeof(struct tng_particle_data) *
                frame_set->n_particle_data_blocks,
                __FILE__, __LINE__);
            free(frame_set->tr_particle_data);
            frame_set->tr_particle_data = 0;
            return(TNG_CRITICAL);
        }
        frame_set->tr_particle_data = data;
    }
    else
    {
        tng_data->n_particle_data_blocks++;
        data = realloc(tng_data->non_tr_particle_data,
                        sizeof(struct tng_particle_data) *
                        tng_data->n_particle_data_blocks);
        if(!data)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%lu bytes). %s: %d\n",
                    sizeof(struct tng_particle_data) *
                    tng_data->n_particle_data_blocks,
                    __FILE__, __LINE__);
            free(tng_data->non_tr_particle_data);
            tng_data->non_tr_particle_data = 0;
            return(TNG_CRITICAL);
        }
        tng_data->non_tr_particle_data = data;
    }

    return(TNG_SUCCESS);
}

static tng_function_status tng_compress(tng_trajectory_t tng_data,
                                        tng_gen_block_t block,
                                        const int64_t n_frames,
                                        const int64_t n_particles,
                                        const char type,
                                        void *start_pos)
{
    int nalgo;
    int new_len;
    int *alt_algo = 0;
    char *dest, *temp;
    int64_t algo_find_n_frames;
    unsigned long offset;
    float f_precision;
    double d_precision;

    if(block->id != TNG_TRAJ_POSITIONS &&
       block->id != TNG_TRAJ_VELOCITIES)
    {
        fprintf(stderr, "TNG library: Can only compress positions and velocities with the "
               "TNG method. %s: %d\n", __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    if(type != TNG_FLOAT_DATA && type != TNG_DOUBLE_DATA)
    {
        fprintf(stderr, "TNG library: Data type not supported. %s: %d\n", __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    if(n_frames <= 0 || n_particles <= 0)
    {
        fprintf(stderr, "TNG library: Missing frames or particles. Cannot compress data "
               "with the TNG method. %s: %d\n", __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    f_precision = 1/(float)tng_data->compression_precision;
    d_precision = 1/tng_data->compression_precision;

    if(block->id == TNG_TRAJ_POSITIONS)
    {
        /* If there is only one frame in this frame set and there might be more
         * do not store the algorithm as the compression algorithm, but find
         * the best one without storing it */
        if(n_frames == 1 && tng_data->frame_set_n_frames > 1)
        {
            nalgo = tng_compress_nalgo();
            alt_algo=malloc(nalgo * sizeof *tng_data->compress_algo_pos);
            if(type == TNG_FLOAT_DATA)
            {
                dest = tng_compress_pos_float_find_algo(start_pos, (int)n_particles,
                                                        (int)n_frames,
                                                        f_precision,
                                                        0, alt_algo,
                                                        &new_len);

            }
            else
            {
                dest = tng_compress_pos_find_algo(start_pos, (int)n_particles,
                                           (int)n_frames,
                                           d_precision,
                                           0, alt_algo,
                                           &new_len);
            }
        }
        else if(!tng_data->compress_algo_pos)
        {
            if(n_frames > 10)
            {
                algo_find_n_frames = 5;
            }
            else
            {
                algo_find_n_frames = n_frames;
            }

            nalgo = tng_compress_nalgo();
            tng_data->compress_algo_pos=malloc(nalgo *
                                           sizeof *tng_data->compress_algo_pos);
            if(type == TNG_FLOAT_DATA)
            {
                dest = tng_compress_pos_float_find_algo(start_pos, (int)n_particles,
                                                        (int)algo_find_n_frames,
                                                        f_precision,
                                                        0, tng_data->
                                                        compress_algo_pos,
                                                        &new_len);

                if(algo_find_n_frames < n_frames)
                {
                    dest = tng_compress_pos_float(start_pos, (int)n_particles,
                                                  (int)n_frames,
                                                  f_precision,
                                                  0, tng_data->compress_algo_pos,
                                                  &new_len);
                }
            }
            else
            {
                dest = tng_compress_pos_find_algo(start_pos, (int)n_particles,
                                           (int)algo_find_n_frames,
                                           d_precision,
                                           0, tng_data->
                                           compress_algo_pos,
                                           &new_len);

                if(algo_find_n_frames < n_frames)
                {
                    dest = tng_compress_pos(start_pos, (int)n_particles,
                                            (int)n_frames,
                                            d_precision, 0,
                                            tng_data->compress_algo_pos,
                                            &new_len);
                }
            }
        }
        else
        {
            if(type == TNG_FLOAT_DATA)
            {
                dest = tng_compress_pos_float(start_pos, (int)n_particles,
                                              (int)n_frames,
                                              f_precision, 0,
                                              tng_data->compress_algo_pos, &new_len);
            }
            else
            {
                dest = tng_compress_pos(start_pos, (int)n_particles,
                                        (int)n_frames,
                                        d_precision, 0,
                                        tng_data->compress_algo_pos,
                                        &new_len);
            }
        }
    }
    else if(block->id == TNG_TRAJ_VELOCITIES)
    {
        /* If there is only one frame in this frame set and there might be more
         * do not store the algorithm as the compression algorithm, but find
         * the best one without storing it */
        if(n_frames == 1 && tng_data->frame_set_n_frames > 1)
        {
            nalgo = tng_compress_nalgo();
            alt_algo=malloc(nalgo * sizeof *tng_data->compress_algo_pos);
            if(type == TNG_FLOAT_DATA)
            {
                dest = tng_compress_vel_float_find_algo(start_pos, (int)n_particles,
                                                        (int)n_frames,
                                                        f_precision,
                                                        0, alt_algo,
                                                        &new_len);

            }
            else
            {
                dest = tng_compress_vel_find_algo(start_pos, (int)n_particles,
                                                  (int)n_frames,
                                                  d_precision,
                                                  0, alt_algo,
                                                  &new_len);
            }
        }
        else if(!tng_data->compress_algo_vel)
        {
            if(n_frames > 10)
            {
                algo_find_n_frames = 5;
            }
            else
            {
                algo_find_n_frames = n_frames;
            }

            nalgo = tng_compress_nalgo();
            tng_data->compress_algo_vel=malloc(nalgo *
                                           sizeof *tng_data->compress_algo_vel);

            if(type == TNG_FLOAT_DATA)
            {
                dest = tng_compress_vel_float_find_algo(start_pos, (int)n_particles,
                                                        (int)algo_find_n_frames,
                                                        f_precision,
                                                        0, tng_data->
                                                        compress_algo_vel,
                                                        &new_len);
                if(algo_find_n_frames < n_frames)
                {
                    dest = tng_compress_vel_float(start_pos, (int)n_particles,
                                                  (int)n_frames,
                                                  f_precision,
                                                  0, tng_data->compress_algo_vel,
                                                  &new_len);
                }
            }
            else
            {
                dest = tng_compress_vel_find_algo(start_pos, (int)n_particles,
                                                  (int)algo_find_n_frames,
                                                  d_precision,
                                                  0, tng_data->
                                                  compress_algo_vel,
                                                  &new_len);
                if(algo_find_n_frames < n_frames)
                {
                    dest = tng_compress_vel(start_pos, (int)n_particles,
                                            (int)n_frames,
                                            d_precision,
                                            0, tng_data->compress_algo_vel,
                                            &new_len);
                }
            }
        }
        else
        {
            if(type == TNG_FLOAT_DATA)
            {
                dest = tng_compress_vel_float(start_pos, (int)n_particles,
                                              (int)n_frames,
                                              f_precision,
                                              0, tng_data->
                                              compress_algo_vel,
                                              &new_len);
            }
            else
            {
                dest = tng_compress_vel(start_pos, (int)n_particles,
                                        (int)n_frames,
                                        d_precision,
                                        0, tng_data->
                                        compress_algo_vel,
                                        &new_len);
            }
        }
    }
    else
    {
        fprintf(stderr, "TNG library: Can only compress positions and velocities using TNG-MF1 algorithms.\n");
        return(TNG_FAILURE);
    }

    offset = (unsigned long)((char *)start_pos - block->block_contents);

    if(alt_algo)
    {
        free(alt_algo);
    }

    block->block_contents_size = new_len + offset;

    temp = realloc(block->block_contents, block->block_contents_size);
    if(!temp)
    {
        free(block->block_contents);
        block->block_contents = 0;
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }
    block->block_contents = temp;
    if(dest)
    {
        memcpy(temp + offset, dest, new_len);
        free(dest);
    }
    else
    {
        fprintf(stderr, "TNG library: Error during TNG compression. %s: %d\n", __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    return(TNG_SUCCESS);
}

static tng_function_status tng_uncompress(tng_trajectory_t tng_data,
                                          tng_gen_block_t block,
                                          const char type,
                                          void *start_pos,
                                          const unsigned long uncompressed_len)
{
    char *temp;
    double *d_dest = 0;
    float *f_dest = 0;
    unsigned long offset;
    int result;
    (void)tng_data;

    TNG_ASSERT(uncompressed_len, "TNG library: The full length of the uncompressed data must be > 0.");

    if(block->id != TNG_TRAJ_POSITIONS &&
       block->id != TNG_TRAJ_VELOCITIES)
    {
        fprintf(stderr, "TNG library: Can only uncompress positions and velocities with the"
               "TNG method.\n");
        return(TNG_FAILURE);
    }
    if(type != TNG_FLOAT_DATA && type != TNG_DOUBLE_DATA)
    {
        fprintf(stderr, "TNG library: Data type not supported.\n");
        return(TNG_FAILURE);
    }

    if(type == TNG_FLOAT_DATA)
    {
        f_dest = malloc(uncompressed_len);
        if(!f_dest)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%lu bytes). %s: %d\n",
                uncompressed_len, __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        result = tng_compress_uncompress_float(start_pos, f_dest);
    }
    else
    {
        d_dest = malloc(uncompressed_len);
        if(!d_dest)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%lu bytes). %s: %d\n",
                uncompressed_len, __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        result = tng_compress_uncompress(start_pos, d_dest);
    }

    if(result == 1)
    {
        fprintf(stderr, "TNG library: Cannot uncompress TNG compressed block.\n");
        return(TNG_FAILURE);
    }

    offset = (unsigned long)((char *)start_pos - (char *)block->block_contents);

    block->block_contents_size = (int64_t)(uncompressed_len + offset);

    temp = realloc(block->block_contents, uncompressed_len + offset);
    if(!temp)
    {
        free(block->block_contents);
        block->block_contents = 0;
        if(d_dest)
        {
            free(d_dest);
        }
        if(f_dest)
        {
            free(f_dest);
        }
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(type == TNG_FLOAT_DATA)
    {
        memcpy(temp + offset, f_dest, uncompressed_len);
    }
    else
    {
        memcpy(temp + offset, d_dest, uncompressed_len);
    }

    block->block_contents = temp;

    if(d_dest)
    {
        free(d_dest);
    }
    if(f_dest)
    {
        free(f_dest);
    }
    return(TNG_SUCCESS);
}

#ifdef USE_ZLIB
static tng_function_status tng_gzip_compress(tng_trajectory_t tng_data,
                                             tng_gen_block_t block,
                                             void *start_pos, const int len)
{
    Bytef *dest;
    char *temp;
    unsigned long max_len, stat, offset;
    (void)tng_data;

    max_len = compressBound(len);
    dest = malloc(max_len);
    if(!dest)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%ld bytes). %s: %d\n",
               max_len, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    stat = compress(dest, &max_len, start_pos, len);
    if(stat != (unsigned long)Z_OK)
    {
        free(dest);
        if(stat == (unsigned long)Z_MEM_ERROR)
        {
            fprintf(stderr, "TNG library: Not enough memory. ");
        }
        else if(stat == (unsigned long)Z_BUF_ERROR)
        {
            fprintf(stderr, "TNG library: Destination buffer too small. ");
        }
        fprintf(stderr, "TNG library: Error gzipping data. %s: %d\n", __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    offset = (char *)start_pos - block->block_contents;

    block->block_contents_size = max_len + offset;

    temp = realloc(block->block_contents, block->block_contents_size);
    if(!temp)
    {
        free(block->block_contents);
        free(dest);
        block->block_contents = 0;
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    block->block_contents = temp;

    memcpy(temp + offset, dest, max_len);

    free(dest);

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
    (void)tng_data;

    offset = (char *)start_pos - (char *)block->block_contents;

    dest = malloc(uncompressed_len);
    if(!dest)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%lud bytes). %s: %d\n",
               uncompressed_len, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    stat = uncompress(dest, &uncompressed_len, (Bytef *) start_pos,
                      block->block_contents_size - offset);

    if(stat != Z_OK)
    {
        free(dest);
        if(stat == (unsigned long)Z_MEM_ERROR)
        {
            fprintf(stderr, "TNG library: Not enough memory. ");
        }
        else if(stat == (unsigned long)Z_BUF_ERROR)
        {
            fprintf(stderr, "TNG library: Destination buffer too small. ");
        }
        else if(stat == (unsigned long)Z_DATA_ERROR)
        {
            fprintf(stderr, "TNG library: Data corrupt. ");
        }
        fprintf(stderr, "TNG library: Error uncompressing gzipped data. %s: %d\n", __FILE__,
               __LINE__);
        return(TNG_FAILURE);
    }


    block->block_contents_size = uncompressed_len + offset;

    temp = realloc(block->block_contents, uncompressed_len + offset);
    if(!temp)
    {
        free(block->block_contents);
        block->block_contents = 0;
        free(dest);
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
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
    void ***values;
    int64_t i, j, k, size, frame_alloc;
    (void)tng_data;

    if(n_particles == 0 || n_values_per_frame == 0)
    {
        return(TNG_FAILURE);
    }

    if(data->strings && data->datatype == TNG_CHAR_DATA)
    {
        for(i = 0; i < data->n_frames; i++)
        {
            for(j = 0; j < n_particles; j++)
            {
                for(k = 0; k < data->n_values_per_frame; k++)
                {
                    if(data->strings[i][j][k])
                    {
                        free(data->strings[i][j][k]);
                    }
                }
                free(data->strings[i][j]);
            }
            free(data->strings[i]);
        }
        free(data->strings);
    }
    data->n_frames = n_frames;
    n_frames = tng_max_i64(1, n_frames);
    data->stride_length = tng_max_i64(1, stride_length);
    data->n_values_per_frame = n_values_per_frame;
    frame_alloc = (n_frames % stride_length) ? n_frames / stride_length + 1 : n_frames / stride_length;

    if(data->datatype == TNG_CHAR_DATA)
    {
        data->strings = malloc(sizeof(char ***) * frame_alloc);
        for(i = 0; i < frame_alloc; i++)
        {
            data->strings[i] = malloc(sizeof(char **) *
                                    n_particles);
            if(!data->strings[i])
            {
                fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                    sizeof(union data_values *) * n_particles,
                    __FILE__, __LINE__);
                return(TNG_CRITICAL);
            }
            for(j = 0; j < n_particles; j++)
            {
                data->strings[i][j] = malloc(sizeof(char *) *
                                            n_values_per_frame);
                if(!data->strings[i][j])
                {
                    fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                        sizeof(union data_values) * n_values_per_frame,
                        __FILE__, __LINE__);
                    return(TNG_CRITICAL);
                }
                for(k = 0; k < n_values_per_frame; k++)
                {
                    data->strings[i][j][k] = 0;
                }
            }
        }
    }
    else
    {
        switch(data->datatype)
        {
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

        values = realloc(data->values,
                         size * frame_alloc *
                         n_particles * n_values_per_frame);
        if(!values)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   size * frame_alloc *
                   n_particles * n_values_per_frame,
                   __FILE__, __LINE__);
            free(data->values);
            data->values = 0;
            return(TNG_CRITICAL);
        }
        data->values = values;
    }
    return(TNG_SUCCESS);
}

static tng_function_status tng_particle_data_find
                (tng_trajectory_t tng_data,
                 const int64_t id,
                 tng_particle_data_t *data)
{
    int64_t block_index, i;
    tng_trajectory_frame_set_t frame_set = &tng_data->
                                           current_trajectory_frame_set;
    char block_type_flag;

    if(tng_data->current_trajectory_frame_set_input_file_pos > 0 ||
       tng_data->current_trajectory_frame_set_output_file_pos > 0)
    {
        block_type_flag = TNG_TRAJECTORY_BLOCK;
    }
    else
    {
        block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
    }

    block_index = -1;
    if(block_type_flag == TNG_TRAJECTORY_BLOCK)
    {
        for(i = 0; i < frame_set->n_particle_data_blocks; i++)
        {
            *data = &frame_set->tr_particle_data[i];
            if((*data)->block_id == id)
            {
                block_index = i;
                break;
            }
        }
    }
    else
    {
        for(i = 0; i < tng_data->n_particle_data_blocks; i++)
        {
            *data = &tng_data->non_tr_particle_data[i];
            if((*data)->block_id == id)
            {
                block_index = i;
                break;
            }
        }
    }
    if(block_index == -1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

static tng_function_status tng_data_find
                (tng_trajectory_t tng_data,
                 const int64_t id,
                 tng_non_particle_data_t *data)
{
    int64_t block_index, i;
    tng_trajectory_frame_set_t frame_set = &tng_data->
                                           current_trajectory_frame_set;
    char block_type_flag;

    if(tng_data->current_trajectory_frame_set_input_file_pos > 0 ||
       tng_data->current_trajectory_frame_set_output_file_pos > 0)
    {
        block_type_flag = TNG_TRAJECTORY_BLOCK;
    }
    else
    {
        block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
    }

    block_index = -1;
    if(block_type_flag == TNG_TRAJECTORY_BLOCK)
    {
        for(i = 0; i < frame_set->n_data_blocks; i++)
        {
            *data = &frame_set->tr_data[i];
            if((*data)->block_id == id)
            {
                block_index = i;
                break;
            }
        }
        if(block_index == -1)
        {
            for(i = 0; i < tng_data->n_data_blocks; i++)
            {
                *data = &tng_data->non_tr_data[i];
                if((*data)->block_id == id)
                {
                    block_index = i;
                    break;
                }
            }
        }
    }
    else
    {
        for(i = 0; i < tng_data->n_data_blocks; i++)
        {
            *data = &tng_data->non_tr_data[i];
            if((*data)->block_id == id)
            {
                block_index = i;
                break;
            }
        }
    }
    if(block_index == -1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

static tng_function_status tng_data_block_len_calculate
                (const tng_trajectory_t tng_data,
                 const tng_particle_data_t data,
                 const tng_bool is_particle_data,
                 const int64_t n_frames,
                 const int64_t frame_step,
                 const int64_t stride_length,
                 const int64_t num_first_particle,
                 const int64_t n_particles,
                 const char dependency,
                 int64_t *data_start_pos,
                 int64_t *len)
{
    int size;
    int64_t i, j, k;
    char ***first_dim_values, **second_dim_values;
    (void)tng_data;

    if(data == 0)
    {
        return(TNG_SUCCESS);
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

    *len = sizeof(char) * 2 + sizeof(data->n_values_per_frame) +
           sizeof(data->codec_id);
    if(is_particle_data)
    {
        *len += sizeof(num_first_particle) + sizeof(n_particles);
    }

    if(stride_length > 1)
    {
        *len += sizeof(data->first_frame_with_data) +
                sizeof(data->stride_length);
    }

    if(data->codec_id != TNG_UNCOMPRESSED)
    {
        *len += sizeof(data->compression_multiplier);
    }

    if(dependency & TNG_FRAME_DEPENDENT)
    {
        *len += sizeof(char);
    }

    *data_start_pos = *len;

    if(data->datatype == TNG_CHAR_DATA)
    {
        if(is_particle_data)
        {
            for(i = 0; i < n_frames; i++)
            {
                first_dim_values = data->strings[i];
                for(j = num_first_particle; j < num_first_particle + n_particles;
                    j++)
                {
                    second_dim_values = first_dim_values[j];
                    for(k = 0; k < data->n_values_per_frame; k++)
                    {
                        *len += strlen(second_dim_values[k]) + 1;
                    }
                }
            }
        }
        else
        {
            for(i = 0; i < n_frames; i++)
            {
                second_dim_values = ((tng_non_particle_data_t)data)->strings[i];
                for(j = 0; j < data->n_values_per_frame; j++)
                {
                    *len += strlen(second_dim_values[j]) + 1;
                }
            }
        }
    }
    else
    {
        *len += size * frame_step * n_particles * data->n_values_per_frame;
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
                 const double multiplier)
{
    int64_t i, j, k, tot_n_particles, n_frames_div;
    int size, len;
    unsigned long data_size;
    char ***first_dim_values, **second_dim_values;
    tng_particle_data_t data;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;
    char block_type_flag;

    TNG_ASSERT(offset != 0, "TNG library: offset must not be a NULL pointer.");

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

    /* If the block does not exist, create it */
    if(tng_particle_data_find(tng_data, block->id, &data) != TNG_SUCCESS)
    {
        if(tng_data->current_trajectory_frame_set_input_file_pos > 0)
        {
            block_type_flag = TNG_TRAJECTORY_BLOCK;
        }
        else
        {
            block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
        }

        if(tng_particle_data_block_create(tng_data, block_type_flag) !=
           TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot create particle data block. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
                   (int)strlen(block->name)+1, __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        strcpy(data->block_name, block->name);

        data->datatype = datatype;

        data->values = 0;
        /* FIXME: Memory leak from strings. */
        data->strings = 0;
        data->n_frames = 0;
        data->codec_id = codec_id;
        data->compression_multiplier = multiplier;
        data->last_retrieved_frame = -1;
    }

    if(/*block_type_flag == TNG_TRAJECTORY_BLOCK &&*/
       tng_data->current_trajectory_frame_set_input_file_pos > 0 &&
       tng_data->var_num_atoms_flag)
    {
        tot_n_particles = frame_set->n_particles;
    }
    else
    {
        tot_n_particles = tng_data->n_particles;
    }

    n_frames_div = (n_frames % stride_length) ? n_frames / stride_length + 1 : n_frames / stride_length;

    if(codec_id != TNG_UNCOMPRESSED)
    {
        data_size = (unsigned long)(n_frames_div * size * n_particles * n_values);
        switch(codec_id)
        {
        case TNG_XTC_COMPRESSION:
            fprintf(stderr, "TNG library: XTC compression not implemented yet.\n");
            break;
        case TNG_TNG_COMPRESSION:
/*            fprintf(stderr, "TNG library: Before TNG uncompression: %"PRId64"\n", block->block_contents_size);*/
            if(tng_uncompress(tng_data, block, datatype,
                              block->block_contents + *offset,
                              data_size) != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Could not read tng compressed block data. %s: %d\n",
                       __FILE__, __LINE__);
                return(TNG_CRITICAL);
            }
/*            fprintf(stderr, "TNG library: After TNG uncompression: %"PRId64"\n", block->block_contents_size);*/
            break;
#ifdef USE_ZLIB
        case TNG_GZIP_COMPRESSION:
/*            fprintf(stderr, "TNG library: Before GZIP uncompression: %"PRId64"\n", block->block_contents_size);*/
            if(tng_gzip_uncompress(tng_data, block,
                                   block->block_contents + *offset,
                                   data_size) != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Could not read gzipped block data. %s: %d\n", __FILE__,
                    __LINE__);
                return(TNG_CRITICAL);
            }
/*            fprintf(stderr, "TNG library: After GZIP uncompression: %"PRId64"\n", block->block_contents_size);*/
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
            fprintf(stderr, "TNG library: Cannot allocate memory for particle data. %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    data->first_frame_with_data = first_frame_with_data;

    if(datatype == TNG_CHAR_DATA)
    {
        for(i = 0; i < n_frames_div; i++)
        {
            first_dim_values = data->strings[i];
            for(j = num_first_particle; j < num_first_particle + n_particles;
                j++)
            {
                second_dim_values = first_dim_values[j];
                for(k = 0; k < n_values; k++)
                {
                    len = tng_min_i((int)strlen(block->block_contents+*offset) + 1,
                              TNG_MAX_STR_LEN);
                    if(second_dim_values[k])
                    {
                        free(second_dim_values[k]);
                    }
                    second_dim_values[k] = malloc(len);
                    if(!second_dim_values[k])
                    {
                        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
                            len, __FILE__, __LINE__);
                        return(TNG_CRITICAL);
                    }
                    strncpy(second_dim_values[k],
                            block->block_contents+*offset, len);
                    *offset += len;
                }
            }
        }
    }
    else
    {
        memcpy((char *)data->values + n_frames_div * size * n_values *
               num_first_particle,
               block->block_contents + *offset,
               block->block_contents_size - *offset);
        switch(datatype)
        {
        case TNG_FLOAT_DATA:
            if(tng_data->input_endianness_swap_func_32)
            {
                for(i = 0; i < (block->block_contents_size - *offset); i+=size)
                {
                    if(tng_data->input_endianness_swap_func_32(tng_data,
                        (int32_t *)((char *)data->values + i))
                        != TNG_SUCCESS)
                    {
                        fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                __FILE__, __LINE__);
                    }
                }
            }
            break;
        case TNG_INT_DATA:
        case TNG_DOUBLE_DATA:
            if(tng_data->input_endianness_swap_func_64)
            {
                for(i = 0; i < (block->block_contents_size - *offset); i+=size)
                {
                    if(tng_data->input_endianness_swap_func_64(tng_data,
                        (int64_t *)((char *)data->values + i))
                        != TNG_SUCCESS)
                    {
                        fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                __FILE__, __LINE__);
                    }
                }
            }
            break;
        case TNG_CHAR_DATA:
            break;
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
                 const int64_t block_index,
                 const tng_particle_mapping_t mapping,
                 const char hash_mode)
{
    int64_t n_particles, num_first_particle, n_frames, stride_length;
    int64_t frame_step, data_start_pos;
    int64_t i, j, k;
    int size;
    size_t len, offset = 0;
    char dependency, temp, *temp_name;
    double multiplier;
    char ***first_dim_values, **second_dim_values;
    tng_trajectory_frame_set_t frame_set;
    tng_function_status stat;

    tng_particle_data_t data;
    char block_type_flag;

    frame_set = &tng_data->current_trajectory_frame_set;

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

        /* If this data block has not had any data added in this frame set
         * do not write it. */
        if(data->first_frame_with_data < frame_set->first_frame)
        {
            return(TNG_SUCCESS);
        }

        stride_length = tng_max_i64(1, data->stride_length);
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
            fprintf(stderr, "TNG library: Cannot allocate memory (%lud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            free(block->name);
            block->name = 0;
            return(TNG_CRITICAL);
        }
        block->name = temp_name;
    }
    strncpy(block->name, data->block_name, len);
    block->id = data->block_id;

    /* If writing frame independent data data->n_frames is 0, but n_frames
       is used for the loop writing the data (and reserving memory) and needs
       to be at least 1 */
    n_frames = tng_max_i64(1, data->n_frames);

    if(block_type_flag == TNG_TRAJECTORY_BLOCK)
    {
        /* If the frame set is finished before writing the full number of frames
           make sure the data block is not longer than the frame set. */
        n_frames = tng_min_i64(n_frames, frame_set->n_frames);

        n_frames -= (data->first_frame_with_data - frame_set->first_frame);
    }

    frame_step = (n_frames % stride_length) ? n_frames / stride_length + 1:
                 n_frames / stride_length;

    /* TNG compression will use compression precision to get integers from
     * floating point data. The compression multiplier stores that information
     * to be able to return the precision of the compressed data. */
    if(data->codec_id == TNG_TNG_COMPRESSION)
    {
        data->compression_multiplier = tng_data->compression_precision;
    }
    /* Uncompressed data blocks do not use compression multipliers at all.
     * GZip compression does not need it either. */
    else if(data->codec_id == TNG_UNCOMPRESSED || data->codec_id == TNG_GZIP_COMPRESSION)
    {
        data->compression_multiplier = 1.0;
    }

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

    if(block_type_flag == TNG_TRAJECTORY_BLOCK && data->n_frames > 0)
    {
        dependency = TNG_FRAME_DEPENDENT + TNG_PARTICLE_DEPENDENT;
    }
    else
    {
        dependency = TNG_PARTICLE_DEPENDENT;
    }

    if(tng_data_block_len_calculate(tng_data, data, TNG_TRUE, n_frames,
                                    frame_step, stride_length, num_first_particle,
                                    n_particles, dependency, &data_start_pos,
                                    &block->block_contents_size) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot calculate length of particle data block. %s: %d\n",
                __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->compression_multiplier);
    }

    if(data->n_frames > 0 && stride_length > 1)
    {
        /* FIXME: first_frame_with_data is not reliably set */
        if(data->first_frame_with_data == 0)
        {
            data->first_frame_with_data = frame_set->first_frame;
        }
        memcpy(block->block_contents+offset, &data->first_frame_with_data,
               sizeof(data->first_frame_with_data));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
               (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    offset += sizeof(n_particles);

    if(data->datatype == TNG_CHAR_DATA)
    {
        if(data->strings)
        {
            for(i = 0; i < frame_step; i++)
            {
                first_dim_values = data->strings[i];
                for(j = num_first_particle; j < num_first_particle + n_particles;
                    j++)
                {
                    second_dim_values = first_dim_values[j];
                    for(k = 0; k < data->n_values_per_frame; k++)
                    {
                        len = (unsigned int)strlen(second_dim_values[k]) + 1;
                        strncpy(block->block_contents+offset,
                                second_dim_values[k], len);
                        offset += len;
                    }
                }
            }
        }
    }
    else if(data->values)
    {
        memcpy(block->block_contents + offset, data->values,
               block->block_contents_size - offset);

        switch(data->datatype)
        {
        case TNG_FLOAT_DATA:
            if(data->codec_id == TNG_UNCOMPRESSED || data-> codec_id == TNG_GZIP_COMPRESSION ||
               data->codec_id == TNG_TNG_COMPRESSION)
            {
                if(tng_data->input_endianness_swap_func_32)
                {
                    for(i = offset; i < block->block_contents_size; i+=size)
                    {
                        if(tng_data->input_endianness_swap_func_32(tng_data,
                           (int32_t *)(block->block_contents + i))
                           != TNG_SUCCESS)
                        {
                            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }
                }
            }
            else
            {
                multiplier = data->compression_multiplier;
                if(fabs(multiplier - 1.0) > 0.00001 ||
                   tng_data->input_endianness_swap_func_32)
                {
                    for(i = offset; i < block->block_contents_size; i+=size)
                    {
                        *(float *)(block->block_contents + i) *= (float)multiplier;
                        if(tng_data->input_endianness_swap_func_32 &&
                        tng_data->input_endianness_swap_func_32(tng_data,
                        (int32_t *)(block->block_contents + i))
                        != TNG_SUCCESS)
                        {
                            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }
                }
            }
            break;
        case TNG_INT_DATA:
            if(tng_data->input_endianness_swap_func_64)
            {
                for(i = offset; i < block->block_contents_size; i+=size)
                {
                    if(tng_data->input_endianness_swap_func_64(tng_data,
                       (int64_t *)(block->block_contents + i))
                       != TNG_SUCCESS)
                    {
                        fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                __FILE__, __LINE__);
                    }
                }
            }
            break;
        case TNG_DOUBLE_DATA:
            if(data->codec_id == TNG_UNCOMPRESSED || data-> codec_id == TNG_GZIP_COMPRESSION ||
               data->codec_id == TNG_TNG_COMPRESSION)
            {
                if(tng_data->input_endianness_swap_func_64)
                {
                    for(i = offset; i < block->block_contents_size; i+=size)
                    {
                        if(tng_data->input_endianness_swap_func_64(tng_data,
                           (int64_t *)(block->block_contents + i))
                           != TNG_SUCCESS)
                        {
                            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }
                }
            }
            else
            {
                multiplier = data->compression_multiplier;
                if(fabs(multiplier - 1.0) > 0.00001 ||
                   tng_data->input_endianness_swap_func_64)
                {
                    for(i = offset; i < block->block_contents_size; i+=size)
                    {
                        *(double *)(block->block_contents + i) *= multiplier;
                        if(tng_data->input_endianness_swap_func_64 &&
                        tng_data->input_endianness_swap_func_64(tng_data,
                        (int64_t *)(block->block_contents + i))
                        != TNG_SUCCESS)
                        {
                            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }
                }
            }
            break;
        case TNG_CHAR_DATA:
            break;
        }
    }
    else
    {
        memset(block->block_contents+offset, 0, block->block_contents_size - offset);
    }

    frame_set->n_written_frames += frame_set->n_unwritten_frames;
    frame_set->n_unwritten_frames = 0;

    if(block_type_flag == TNG_NON_TRAJECTORY_BLOCK || frame_set->n_written_frames > 0)
    {
        switch(data->codec_id)
        {
        case TNG_XTC_COMPRESSION:
            fprintf(stderr, "TNG library: XTC compression not implemented yet.\n");
            data->codec_id = TNG_UNCOMPRESSED;
            break;
        case TNG_TNG_COMPRESSION:
            stat = tng_compress(tng_data, block, frame_step,
                                n_particles, data->datatype,
                                block->block_contents + data_start_pos);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Could not write tng compressed block data. %s: %d\n",
                    __FILE__, __LINE__);
                if(stat == TNG_CRITICAL)
                {
                    return(TNG_CRITICAL);
                }
                /* Set the data again, but with no compression (to write only
                 * the relevant data) */
                data->codec_id = TNG_UNCOMPRESSED;
                stat = tng_particle_data_block_write(tng_data, block,
                                                     block_index, mapping,
                                                     hash_mode);
                return(stat);
            }
            break;
#ifdef USE_ZLIB
        case TNG_GZIP_COMPRESSION:
    /*         fprintf(stderr, "TNG library: Before compression: %"PRId64"\n", block->block_contents_size);*/
            stat = tng_gzip_compress(tng_data, block,
                                     block->block_contents + data_start_pos,
                                     block->block_contents_size - data_start_pos);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Could not write gzipped block data. %s: %d\n", __FILE__,
                    __LINE__);
                if(stat == TNG_CRITICAL)
                {
                    return(TNG_CRITICAL);
                }
                /* Set the data again, but with no compression (to write only
                 * the relevant data) */
                data->codec_id = TNG_UNCOMPRESSED;
                stat = tng_particle_data_block_write(tng_data, block,
                                                     block_index, mapping,
                                                     hash_mode);
                return(stat);
            }
    /*         fprintf(stderr, "TNG library: After compression: %"PRId64"\n", block->block_contents_size);*/
            break;
#endif
        }
    }

    if(tng_block_header_write(tng_data, block, hash_mode) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot write header of file %s. %s: %d\n",
               tng_data->output_file_path, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
        tng_data->output_file) != 1)
    {
        fprintf(stderr, "TNG library: Could not write all block data. %s: %d\n", __FILE__,
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
                 const char block_type_flag)
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
            fprintf(stderr, "TNG library: Cannot allocate memory (%lu bytes). %s: %d\n",
                sizeof(struct tng_non_particle_data) * frame_set->n_data_blocks,
                __FILE__, __LINE__);
            free(frame_set->tr_data);
            frame_set->tr_data = 0;
            return(TNG_CRITICAL);
        }
        frame_set->tr_data = data;
    }
    else
    {
        tng_data->n_data_blocks++;
        data = realloc(tng_data->non_tr_data, sizeof(struct tng_non_particle_data) *
                        tng_data->n_data_blocks);
        if(!data)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%lu bytes). %s: %d\n",
                sizeof(struct tng_non_particle_data) * tng_data->n_data_blocks,
                __FILE__, __LINE__);
            free(tng_data->non_tr_data);
            tng_data->non_tr_data = 0;
            return(TNG_CRITICAL);
        }
        tng_data->non_tr_data = data;
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
    void **values;
    int64_t i, j, size, frame_alloc;
    (void)tng_data;

    if(n_values_per_frame == 0)
    {
        return(TNG_FAILURE);
    }

    if(data->strings && data->datatype == TNG_CHAR_DATA)
    {
        for(i = 0; i < data->n_frames; i++)
        {
            for(j = 0; j < data->n_values_per_frame; j++)
            {
                if(data->strings[i][j])
                {
                    free(data->strings[i][j]);
                    data->strings[i][j] = 0;
                }
            }
            free(data->strings[i]);
            data->strings[i] = 0;
        }
        free(data->strings);
    }
    data->n_frames = n_frames;
    data->stride_length = tng_max_i64(1, stride_length);
    n_frames = tng_max_i64(1, n_frames);
    data->n_values_per_frame = n_values_per_frame;
    frame_alloc = (n_frames % stride_length) ? n_frames / stride_length + 1 : n_frames / stride_length;

    if(data->datatype == TNG_CHAR_DATA)
    {
        data->strings = malloc(sizeof(char **) * frame_alloc);
        for(i = 0; i < frame_alloc; i++)
        {
            data->strings[i] = malloc(sizeof(char *) * n_values_per_frame);
            if(!data->strings[i])
            {
                fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                       n_values_per_frame,
                       __FILE__, __LINE__);
                return(TNG_CRITICAL);
            }
            for(j = 0; j < n_values_per_frame; j++)
            {
                data->strings[i][j] = 0;
            }
        }
    }
    else
    {
        switch(data->datatype)
        {
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

        values = realloc(data->values,
                         size * frame_alloc *
                         n_values_per_frame);
        if(!values)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   size * frame_alloc *
                   n_values_per_frame,
                   __FILE__, __LINE__);
            free(data->values);
            data->values = 0;
            return(TNG_CRITICAL);
        }
        data->values = values;
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
                                         const double multiplier)
{
    int64_t i, j, n_frames_div;
    int size, len;
#ifdef USE_ZLIB
    unsigned long data_size;
#endif
    tng_non_particle_data_t data;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;
    char block_type_flag;

    TNG_ASSERT(offset != 0, "TNG library: offset must not be a NULL pointer.");

/*     fprintf(stderr, "TNG library: %s\n", block->name);*/

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

    /* If the block does not exist, create it */
    if(tng_data_find(tng_data, block->id, &data) != TNG_SUCCESS)
    {
        if(tng_data->current_trajectory_frame_set_input_file_pos > 0)
        {
            block_type_flag = TNG_TRAJECTORY_BLOCK;
        }
        else
        {
            block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
        }

        if(tng_data_block_create(tng_data, block_type_flag) !=
            TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot create particle data block. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
                   (int)strlen(block->name)+1, __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        strcpy(data->block_name, block->name);

        data->datatype = datatype;

        data->values = 0;
        /* FIXME: Memory leak from strings. */
        data->strings = 0;
        data->n_frames = 0;
        data->codec_id = codec_id;
        data->compression_multiplier = multiplier;
        data->last_retrieved_frame = -1;
    }

    n_frames_div = (n_frames % stride_length) ? n_frames / stride_length + 1 : n_frames / stride_length;

    if(codec_id != TNG_UNCOMPRESSED)
    {
        switch(codec_id)
        {
#ifdef USE_ZLIB
        case TNG_GZIP_COMPRESSION:
            data_size = n_frames_div * size * n_values;
    /*         fprintf(stderr, "TNG library: Before compression: %"PRId64"\n", block->block_contents_size); */
            if(tng_gzip_uncompress(tng_data, block,
                                   block->block_contents + *offset,
                                   data_size) != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Could not read gzipped block data. %s: %d\n", __FILE__,
                    __LINE__);
                return(TNG_CRITICAL);
            }
    /*         fprintf(stderr, "TNG library: After compression: %"PRId64"\n", block->block_contents_size); */
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
            fprintf(stderr, "TNG library: Cannot allocate memory for data. %s: %d\n",
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    data->first_frame_with_data = first_frame_with_data;

    if(datatype == TNG_CHAR_DATA)
    {
        for(i = 0; i < n_frames_div; i++)
        {
            for(j = 0; j < n_values; j++)
            {
                len = tng_min_i((int)strlen(block->block_contents+*offset) + 1,
                              TNG_MAX_STR_LEN);
                if(data->strings[i][j])
                {
                    free(data->strings[i][j]);
                }
                data->strings[i][j] = malloc(len);
                if(!data->strings[i][j])
                {
                    fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
                           len, __FILE__, __LINE__);
                    return(TNG_CRITICAL);
                }
                strncpy(data->strings[i][j], block->block_contents+*offset,
                        len);
                *offset += len;
            }
        }
    }
    else
    {
        memcpy(data->values, block->block_contents + *offset,
               block->block_contents_size - *offset);
        switch(datatype)
        {
        case TNG_FLOAT_DATA:
            if(tng_data->input_endianness_swap_func_32)
            {
                for(i = 0; i < (block->block_contents_size - *offset); i+=size)
                {
                    if(tng_data->input_endianness_swap_func_32(tng_data,
                        (int32_t *)((char *)data->values + i))
                        != TNG_SUCCESS)
                    {
                        fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                __FILE__, __LINE__);
                    }
                }
            }
            break;
        case TNG_INT_DATA:
        case TNG_DOUBLE_DATA:
            if(tng_data->input_endianness_swap_func_64)
            {
                for(i = 0; i < (block->block_contents_size - *offset); i+=size)
                {
                    if(tng_data->input_endianness_swap_func_64(tng_data,
                        (int64_t *)((char *)data->values + i))
                        != TNG_SUCCESS)
                    {
                        fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                __FILE__, __LINE__);
                    }
                }
            }
            break;
        case TNG_CHAR_DATA:
            break;
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
                                                const int64_t block_index,
                                                const char hash_mode)
{
    int64_t n_frames, stride_length, frame_step, data_start_pos;
    int64_t i, j;
    int offset = 0, size;
    unsigned int len;
#ifdef USE_ZLIB
    tng_function_status stat;
#endif
    char temp, dependency, *temp_name;
    double multiplier;
    tng_trajectory_frame_set_t frame_set =
    &tng_data->current_trajectory_frame_set;

    tng_non_particle_data_t data;
    char block_type_flag;

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

        /* If this data block has not had any data added in this frame set
         * do not write it. */
        if(data->first_frame_with_data < frame_set->first_frame)
        {
            return(TNG_SUCCESS);
        }

        stride_length = tng_max_i64(1, data->stride_length);
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

    len = (unsigned int)strlen(data->block_name) + 1;

    if(!block->name || strlen(block->name) < len)
    {
        temp_name = realloc(block->name, len);
        if(!temp_name)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len+1,
                   __FILE__, __LINE__);
            free(block->name);
            block->name = 0;
            return(TNG_CRITICAL);
        }
        block->name = temp_name;
    }
    strncpy(block->name, data->block_name, len);
    block->id = data->block_id;

    /* If writing frame independent data data->n_frames is 0, but n_frames
       is used for the loop writing the data (and reserving memory) and needs
       to be at least 1 */
    n_frames = tng_max_i64(1, data->n_frames);

    if(block_type_flag == TNG_TRAJECTORY_BLOCK)
    {
        /* If the frame set is finished before writing the full number of frames
           make sure the data block is not longer than the frame set. */
        n_frames = tng_min_i64(n_frames, frame_set->n_frames);

        n_frames -= (data->first_frame_with_data - frame_set->first_frame);
    }

    frame_step = (n_frames % stride_length) ? n_frames / stride_length + 1:
                 n_frames / stride_length;

    /* TNG compression will use compression precision to get integers from
     * floating point data. The compression multiplier stores that information
     * to be able to return the precision of the compressed data. */
    if(data->codec_id == TNG_TNG_COMPRESSION)
    {
        data->compression_multiplier = tng_data->compression_precision;
    }
    /* Uncompressed data blocks do not use compression multipliers at all.
     * GZip compression does not need it either. */
    else if(data->codec_id == TNG_UNCOMPRESSED || data->codec_id == TNG_GZIP_COMPRESSION)
    {
        data->compression_multiplier = 1.0;
    }

    if(block_type_flag == TNG_TRAJECTORY_BLOCK && data->n_frames > 0)
    {
        dependency = TNG_FRAME_DEPENDENT;
    }
    else
    {
        dependency = 0;
    }

    if(tng_data_block_len_calculate(tng_data, (tng_particle_data_t)data, TNG_FALSE, n_frames,
                                    frame_step, stride_length, 0,
                                    1, dependency, &data_start_pos,
                                    &block->block_contents_size) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot calculate length of non-particle data block. %s: %d\n",
                __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->compression_multiplier);
    }

    if(data->n_frames > 0 && stride_length > 1)
    {
        /* FIXME: first_frame_with_data is not reliably set */
        if(data->first_frame_with_data == 0)
        {
            data->first_frame_with_data = frame_set->first_frame;
        }
        memcpy(block->block_contents+offset, &data->first_frame_with_data,
               sizeof(data->first_frame_with_data));
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
            (int64_t *)block->header_contents+offset)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->stride_length);
    }

    if(data->datatype == TNG_CHAR_DATA)
    {
        if(data->strings)
        {
            for(i = 0; i < frame_step; i++)
            {
                for(j = 0; j < data->n_values_per_frame; j++)
                {
                    len = (unsigned int)strlen(data->strings[i][j]) + 1;
                    strncpy(block->block_contents+offset, data->strings[i][j],
                            len);
                    offset += len;
                }
            }
        }
    }
    else if(data->values)
    {
        memcpy(block->block_contents + offset, data->values,
               block->block_contents_size - offset);
        switch(data->datatype)
        {
        case TNG_FLOAT_DATA:
            if(data->codec_id == TNG_UNCOMPRESSED || data-> codec_id == TNG_GZIP_COMPRESSION ||
               data->codec_id == TNG_TNG_COMPRESSION)
            {
                if(tng_data->input_endianness_swap_func_32)
                {
                    for(i = offset; i < block->block_contents_size; i+=size)
                    {
                        if(tng_data->input_endianness_swap_func_32(tng_data,
                           (int32_t *)(block->block_contents + i))
                           != TNG_SUCCESS)
                        {
                            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }
                }
            }
            else
            {
                multiplier = data->compression_multiplier;
                if(fabs(multiplier - 1.0) > 0.00001 ||
                   tng_data->input_endianness_swap_func_32)
                {
                    for(i = offset; block->block_contents_size; i+=size)
                    {
                        *(float *)(block->block_contents + i) *= (float)multiplier;
                        if(tng_data->input_endianness_swap_func_32 &&
                        tng_data->input_endianness_swap_func_32(tng_data,
                        (int32_t *)(block->block_contents + i))
                        != TNG_SUCCESS)
                        {
                            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }
                }
            }
            break;
        case TNG_INT_DATA:
            if(tng_data->input_endianness_swap_func_64)
            {
                for(i = offset; i < block->block_contents_size; i+=size)
                {
                    if(tng_data->input_endianness_swap_func_64(tng_data,
                       (int64_t *)(block->block_contents + i))
                       != TNG_SUCCESS)
                    {
                        fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                __FILE__, __LINE__);
                    }
                }
            }
            break;
        case TNG_DOUBLE_DATA:
            if(data->codec_id == TNG_UNCOMPRESSED || data-> codec_id == TNG_GZIP_COMPRESSION ||
               data->codec_id == TNG_TNG_COMPRESSION)
            {
                if(tng_data->input_endianness_swap_func_64)
                {
                    for(i = offset; i < block->block_contents_size; i+=size)
                    {
                        if(tng_data->input_endianness_swap_func_64(tng_data,
                           (int64_t *)(block->block_contents + i))
                           != TNG_SUCCESS)
                        {
                            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }
                }
            }
            else
            {
                multiplier = data->compression_multiplier;
                if(fabs(multiplier - 1.0) > 0.00001 ||
                   tng_data->input_endianness_swap_func_64)
                {
                    for(i = offset; i < block->block_contents_size; i+=size)
                    {
                        *(double *)(block->block_contents + i) *= multiplier;
                        if(tng_data->input_endianness_swap_func_64 &&
                        tng_data->input_endianness_swap_func_64(tng_data,
                        (int64_t *)(block->block_contents + i))
                        != TNG_SUCCESS)
                        {
                            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }
                }
            }
            break;
        case TNG_CHAR_DATA:
            break;
        }
    }
    else
    {
        memset(block->block_contents+offset, 0, block->block_contents_size - offset);
    }

    frame_set->n_written_frames += frame_set->n_unwritten_frames;
    frame_set->n_unwritten_frames = 0;

    if(block_type_flag == TNG_NON_TRAJECTORY_BLOCK || frame_set->n_written_frames > 0)
    {
        switch(data->codec_id)
        {
#ifdef USE_ZLIB
        case TNG_GZIP_COMPRESSION:
    /*         fprintf(stderr, "TNG library: Before compression: %"PRId64"\n", block->block_contents_size); */
            stat = tng_gzip_compress(tng_data, block,
                                     block->block_contents + data_start_pos,
                                     block->block_contents_size - data_start_pos);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Could not write gzipped block data. %s: %d\n", __FILE__,
                    __LINE__);
                if(stat == TNG_CRITICAL)
                {
                    return(TNG_CRITICAL);
                }
                data->codec_id = TNG_UNCOMPRESSED;
            }
    /*         fprintf(stderr, "TNG library: After compression: %"PRId64"\n", block->block_contents_size); */
            break;
#endif
        }
    }

    if(tng_block_header_write(tng_data, block, hash_mode) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot write header of file %s. %s: %d\n",
               tng_data->output_file_path, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
              tng_data->output_file) != 1)
    {
        fprintf(stderr, "TNG library: Could not write all block data. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    return(TNG_SUCCESS);
}

/** Read the meta information of a data block (particle or non-particle data).
 * @param tng_data is a trajectory data container.
 * @param block is the block to store the data (should already contain
 * the block headers).
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
static tng_function_status tng_data_block_meta_information_read
                (tng_trajectory_t tng_data,
                 tng_gen_block_t block,
                 int *offset,
                 char *datatype,
                 char *dependency,
                 char *sparse_data,
                 int64_t *n_values,
                 int64_t *codec_id,
                 int64_t *first_frame_with_data,
                 int64_t *stride_length,
                 int64_t *n_frames,
                 int64_t *num_first_particle,
                 int64_t *block_n_particles,
                 double *multiplier)
{
    int meta_size;
    char *contents;

    if(block->block_contents)
    {
        contents = block->block_contents;
    }
    else
    {
        meta_size = 3 * sizeof(char) + sizeof(double) + 6 * sizeof(int64_t);
        contents = malloc(meta_size);
        if(!contents)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
               meta_size, __FILE__, __LINE__);
        }

        if(fread(contents, meta_size, 1, tng_data->input_file) == 0)
        {
            fprintf(stderr, "TNG library: Cannot read data block meta information. %s: %d\n", __FILE__, __LINE__);
            free(contents);
            return(TNG_CRITICAL);
        }
    }

    memcpy(datatype, contents+*offset,
           sizeof(*datatype));
    *offset += sizeof(*datatype);

    memcpy(dependency, contents+*offset,
           sizeof(*dependency));
    *offset += sizeof(*dependency);

    if(*dependency & TNG_FRAME_DEPENDENT)
    {
        memcpy(sparse_data, contents+*offset,
               sizeof(*sparse_data));
        *offset += sizeof(*sparse_data);
    }

    memcpy(n_values, contents+*offset,
        sizeof(*n_values));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   n_values)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(*n_values);

    memcpy(codec_id, contents+*offset,
        sizeof(*codec_id));
    if(tng_data->input_endianness_swap_func_64)
    {
        if(tng_data->input_endianness_swap_func_64(tng_data,
                                                   codec_id)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
    *offset += sizeof(*codec_id);

    if(*codec_id != TNG_UNCOMPRESSED)
    {
        memcpy(multiplier, contents+*offset,
            sizeof(*multiplier));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                       (int64_t *) multiplier)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        *offset += sizeof(*multiplier);
    }
    else
    {
        *multiplier = 1;
    }

    if(*dependency & TNG_FRAME_DEPENDENT)
    {
        if(*sparse_data)
        {
            memcpy(first_frame_with_data, contents+*offset,
                sizeof(*first_frame_with_data));
            if(tng_data->input_endianness_swap_func_64)
            {
                if(tng_data->input_endianness_swap_func_64(tng_data,
                                                           first_frame_with_data)
                    != TNG_SUCCESS)
                {
                    fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }
            *offset += sizeof(*first_frame_with_data);

            memcpy(stride_length, contents+*offset,
                sizeof(*stride_length));
            if(tng_data->input_endianness_swap_func_64)
            {
                if(tng_data->input_endianness_swap_func_64(tng_data,
                                                           stride_length)
                    != TNG_SUCCESS)
                {
                    fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }
            *offset += sizeof(*stride_length);
            *n_frames = tng_data->current_trajectory_frame_set.n_frames -
                        (*first_frame_with_data -
                        tng_data->current_trajectory_frame_set.first_frame);
        }
        else
        {
            *first_frame_with_data = 0;
            *stride_length = 1;
            *n_frames = tng_data->current_trajectory_frame_set.n_frames;
        }
    }
    else
    {
        *first_frame_with_data = 0;
        *stride_length = 1;
        *n_frames = 1;
    }

    if (*dependency & TNG_PARTICLE_DEPENDENT)
    {
        memcpy(num_first_particle, contents+*offset,
               sizeof(*num_first_particle));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                       num_first_particle)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        *offset += sizeof(*num_first_particle);

        memcpy(block_n_particles, contents+*offset,
            sizeof(*block_n_particles));
        if(tng_data->input_endianness_swap_func_64)
        {
            if(tng_data->input_endianness_swap_func_64(tng_data,
                                                       block_n_particles)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }
        *offset += sizeof(*block_n_particles);
    }

    if(!block->block_contents)
    {
        free(contents);
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
                 const char hash_mode)
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
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to
     * disk even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
             tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Cannot read block. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */

    if(hash_mode == TNG_USE_HASH)
    {
        tng_md5_hash_match_verify(block, &same_hash);
        if(same_hash != TNG_TRUE)
        {
            fprintf(stderr, "TNG library: '%s' data block contents corrupt. Hashes do not match. %s: %d\n",
                block->name, __FILE__, __LINE__);
    /*         return(TNG_FAILURE); */
        }
    }

    if(tng_data_block_meta_information_read(tng_data, block,
                                            &offset, &datatype,
                                            &dependency, &sparse_data,
                                            &n_values, &codec_id,
                                            &first_frame_with_data,
                                            &stride_length, &n_frames,
                                            &num_first_particle,
                                            &block_n_particles,
                                            &multiplier) == TNG_CRITICAL)
    {
        fprintf(stderr, "TNG library: Cannot read data block (%s) meta information. %s: %d\n",
            block->name, __FILE__, __LINE__);
        return(TNG_CRITICAL);
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

/*
// ** Move the blocks in a frame set so that there is no unused space between
//  * them. This can only be done on the last frame set in the file and should
//  * be done e.g. if the last frame set in the file has fewer frames than
//  * default or after compressing data blocks in a frame set.
//  * @param tng_data is a trajectory data container.
//  * @details the current_trajectory_frame_set is the one that will be modified.
//  * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the frame set
//  * cannot be aligned or TNG_CRITICAL (2) if a major error has occured.
//  * FIXME: This function is not finished!!!
//  *
// static tng_function_status tng_frame_set_align(tng_trajectory_t tng_data)
// {
//     tng_gen_block_t block;
//     tng_trajectory_frame_set_t frame_set;
//     FILE *temp = tng_data->input_file;
//     int64_t pos, contents_start_pos, output_file_len;
//
//     frame_set = &tng_data->current_trajectory_frame_set;
//
//     if(frame_set->n_written_frames == frame_set->n_frames)
//     {
//         return(TNG_SUCCESS);
//     }
//
//     if(tng_data->current_trajectory_frame_set_output_file_pos !=
//        tng_data->last_trajectory_frame_set_output_file_pos)
//     {
//     }
//
//     if(tng_output_file_init(tng_data) != TNG_SUCCESS)
//     {
//         fprintf(stderr, "TNG library: Cannot initialise destination file. %s: %d\n",
//                __FILE__, __LINE__);
//         return(TNG_CRITICAL);
//     }
//
//     tng_block_init(&block);
// //     output_file_pos = ftell(tng_data->output_file);
//
//     tng_data->input_file = tng_data->output_file;
//
//     pos = tng_data->current_trajectory_frame_set_output_file_pos;
//
//     fseek(tng_data->output_file, pos, SEEK_SET);
//     if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
//     {
//         fprintf(stderr, "TNG library: Cannot read frame set header. %s: %d\n",
//             __FILE__, __LINE__);
//         tng_data->input_file = temp;
//         tng_block_destroy(&block);
//         return(TNG_CRITICAL);
//     }
//
//     contents_start_pos = ftell(tng_data->output_file);
//
//     fseek(tng_data->output_file, 0, SEEK_END);
//     output_file_len = ftell(tng_data->output_file);
//     pos = contents_start_pos + block->block_contents_size;
//     fseek(tng_data->output_file, pos,
//           SEEK_SET);
//
//     while(pos < output_file_len)
//     {
//         if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
//         {
//             fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n", pos,
//                    __FILE__, __LINE__);
//             tng_data->input_file = temp;
//             tng_block_destroy(&block);
//             return(TNG_CRITICAL);
//         }
//         pos += block->header_contents_size + block->block_contents_size;
//         fseek(tng_data->output_file, pos, SEEK_SET);
//     }
//
//     return(TNG_SUCCESS);
// }
*/
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
                (tng_trajectory_t tng_data, const char hash_mode)
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
        fprintf(stderr, "TNG library: Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    tng_block_init(&block);
/*     output_file_pos = ftell(tng_data->output_file); */

    tng_data->input_file = tng_data->output_file;

    pos = tng_data->current_trajectory_frame_set_output_file_pos;

    fseek(tng_data->output_file, (long)pos, SEEK_SET);

    if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot read frame set header. %s: %d\n",
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
    fseek(tng_data->output_file, (long)pos, SEEK_SET);

    while(pos < output_file_len)
    {
        if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n", pos,
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
        fseek(tng_data->output_file, (long)pos, SEEK_SET);
    }

    tng_data->input_file = temp;
    tng_block_destroy(&block);
    return(TNG_SUCCESS);
}

/*
// ** Sets the name of a file contents block
//  * @param tng_data is a trajectory data container.
//  * @param block is the block, of which to change names.
//  * @param new_name is the new name of the block.
//  * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
//  * error has occured.
//
// static tng_function_status tng_block_name_set(tng_trajectory_t tng_data,
//                                               tng_gen_block_t block,
//                                               const char *new_name)
// {
//     int len;
//
//     len = tng_min_i((int)strlen(new_name) + 1, TNG_MAX_STR_LEN);
//
//      * If the currently stored string length is not enough to store the new
//      * string it is freed and reallocated. *
//     if(block->name && strlen(block->name) < len)
//     {
//         free(block->name);
//         block->name = 0;
//     }
//     if(!block->name)
//     {
//         block->name = malloc(len);
//         if(!block->name)
//         {
//             fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len,
//                    __FILE__, __LINE__);
//             return(TNG_CRITICAL);
//         }
//     }
//
//     strncpy(block->name, new_name, len);
//
//     return(TNG_SUCCESS);
// }
*/

tng_function_status tng_atom_residue_get(const tng_trajectory_t tng_data,
                                         const tng_atom_t atom,
                                         tng_residue_t *residue)
{
    (void) tng_data;

    TNG_ASSERT(atom, "TNG library: atom must not be NULL");

    *residue = atom->residue;

    return(TNG_SUCCESS);
}

tng_function_status tng_atom_name_get(const tng_trajectory_t tng_data,
                                      const tng_atom_t atom,
                                      char *name,
                                      const int max_len)
{
    (void) tng_data;
    TNG_ASSERT(atom, "TNG library: atom must not be NULL");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer");

    strncpy(name, atom->name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(atom->name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_atom_name_set(tng_trajectory_t tng_data,
                                      tng_atom_t atom,
                                      const char *new_name)
{
    unsigned int len;
    (void)tng_data;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(new_name, "TNG library: new_name must not be a NULL pointer.");

    len = tng_min_i((int)strlen(new_name) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(atom->name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status tng_atom_type_get(const tng_trajectory_t tng_data,
                                      const tng_atom_t atom,
                                      char *type,
                                      const int max_len)
{
    (void) tng_data;
    TNG_ASSERT(atom, "TNG library: atom must not be NULL");
    TNG_ASSERT(type, "TNG library: type must not be a NULL pointer");

    strncpy(type, atom->atom_type, max_len - 1);
    type[max_len - 1] = 0;

    if(strlen(atom->atom_type) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status tng_atom_type_set(tng_trajectory_t tng_data,
                                      tng_atom_t atom,
                                      const char *new_type)
{
    unsigned int len;
    (void)tng_data;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(new_type, "TNG library: new_type must not be a NULL pointer.");

    len = tng_min_i((int)strlen(new_type) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
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

tng_function_status DECLSPECDLLEXPORT tng_version_major
                (const tng_trajectory_t tng_data,
                 int *version)
{
    (void)tng_data;

    *version = TNG_VERSION_MAJOR;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_version_minor
                (const tng_trajectory_t tng_data,
                 int *version)
{
    (void)tng_data;

    *version = TNG_VERSION_MINOR;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_version_patchlevel
                (const tng_trajectory_t tng_data,
                 int *patch_level)
{
    (void)tng_data;

    *patch_level = TNG_VERSION_PATCHLEVEL;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_version
                (const tng_trajectory_t tng_data,
                 char *version,
                 const int max_len)
{
    (void)tng_data;
    TNG_ASSERT(version, "TNG library: version must not be a NULL pointer");

    TNG_SNPRINTF(version, max_len, "%s", TNG_VERSION);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_add
                (tng_trajectory_t tng_data,
                 const char *name,
                 tng_molecule_t *molecule)
{
    int64_t id;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

    /* Set ID to the ID of the last molecule + 1 */
    if(tng_data->n_molecules)
    {
        id = tng_data->molecules[tng_data->n_molecules-1].id + 1;
    }
    else
    {
        id = 1;
    }

    return(tng_molecule_w_id_add(tng_data, name, id, molecule));
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_w_id_add
                (tng_trajectory_t tng_data,
                 const char *name,
                 const int64_t id,
                 tng_molecule_t *molecule)
{
    tng_molecule_t new_molecules;
    int64_t *new_molecule_cnt_list;
    tng_function_status stat = TNG_SUCCESS;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

    new_molecules = realloc(tng_data->molecules,
                            sizeof(struct tng_molecule) *
                            (tng_data->n_molecules + 1));

    if(!new_molecules)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct tng_molecule) * (tng_data->n_molecules + 1),
               __FILE__, __LINE__);
        free(tng_data->molecules);
        tng_data->molecules = 0;
        return(TNG_CRITICAL);
    }

    new_molecule_cnt_list = realloc(tng_data->molecule_cnt_list,
                                    sizeof(int64_t) *
                                    (tng_data->n_molecules + 1));

    if(!new_molecule_cnt_list)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(int64_t) * (tng_data->n_molecules + 1),
               __FILE__, __LINE__);
        free(tng_data->molecule_cnt_list);
        tng_data->molecule_cnt_list = 0;
        free(new_molecules);
        return(TNG_CRITICAL);
    }

    tng_data->molecules = new_molecules;
    tng_data->molecule_cnt_list = new_molecule_cnt_list;

    *molecule = &new_molecules[tng_data->n_molecules];

    tng_molecule_init(tng_data, *molecule);
    tng_molecule_name_set(tng_data, *molecule, name);

    /* FIXME: Should this be a function argument instead? */
    tng_data->molecule_cnt_list[tng_data->n_molecules] = 0;

    (*molecule)->id = id;

    tng_data->n_molecules++;

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_existing_add
                (tng_trajectory_t tng_data,
                 tng_molecule_t *molecule_p)
{
    int64_t *new_molecule_cnt_list, id;
    tng_molecule_t new_molecules, molecule;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    /* Set ID to the ID of the last molecule + 1 */
    if(tng_data->n_molecules)
    {
        id = tng_data->molecules[tng_data->n_molecules-1].id + 1;
    }
    else
    {
        id = 1;
    }

    new_molecules = realloc(tng_data->molecules,
                            sizeof(struct tng_molecule) *
                            (tng_data->n_molecules + 1));

    if(!new_molecules)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct tng_molecule) * (tng_data->n_molecules + 1),
               __FILE__, __LINE__);
        free(tng_data->molecules);
        tng_data->molecules = 0;
        return(TNG_CRITICAL);
    }

    new_molecule_cnt_list = realloc(tng_data->molecule_cnt_list,
                                    sizeof(int64_t) *
                                    (tng_data->n_molecules + 1));

    if(!new_molecule_cnt_list)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(int64_t) * (tng_data->n_molecules + 1),
               __FILE__, __LINE__);
        free(tng_data->molecule_cnt_list);
        tng_data->molecule_cnt_list = 0;
        free(new_molecules);
        return(TNG_CRITICAL);
    }

    molecule = *molecule_p;

    tng_data->molecules = new_molecules;
    tng_data->molecule_cnt_list = new_molecule_cnt_list;

    new_molecules[tng_data->n_molecules] = *molecule;

    tng_data->molecule_cnt_list[tng_data->n_molecules] = 0;

    free(*molecule_p);

    molecule = &new_molecules[tng_data->n_molecules];

    *molecule_p = molecule;

    molecule->id = id;

    tng_data->n_molecules++;

    return(TNG_SUCCESS);
}

tng_function_status tng_molecule_name_get(const tng_trajectory_t tng_data,
                                          const tng_molecule_t molecule,
                                          char *name,
                                          const int max_len)
{
    (void) tng_data;
    TNG_ASSERT(molecule, "TNG library: molecule must not be NULL");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer");

    strncpy(name, molecule->name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(molecule->name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_name_set
                (tng_trajectory_t tng_data,
                 tng_molecule_t molecule,
                 const char *new_name)
{
    unsigned int len;
    (void)tng_data;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(new_name, "TNG library: new_name must not be a NULL pointer.");

    len = tng_min_i((int)strlen(new_name) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(molecule->name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_cnt_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 int64_t *cnt)
{
    int64_t i, index = -1;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(cnt, "TNG library: cnt must not be a NULL pointer.");

    for(i = 0; i < tng_data->n_molecules; i++)
    {
        if(&tng_data->molecules[i] == molecule)
        {
            index = i;
            break;
        }
    }
    if(index == -1)
    {
        return(TNG_FAILURE);
    }
    *cnt = tng_data->molecule_cnt_list[index];

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_cnt_set
                (tng_trajectory_t tng_data,
                 tng_molecule_t molecule,
                 const int64_t cnt)
{
    int64_t i, old_cnt, index = -1;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    for(i = 0; i < tng_data->n_molecules; i++)
    {
        if(&tng_data->molecules[i] == molecule)
        {
            index = i;
            break;
        }
    }
    if(index == -1)
    {
        fprintf(stderr, "TNG library: Could not find molecule in TNG trajectory. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }
    if(tng_data->var_num_atoms_flag == TNG_CONSTANT_N_ATOMS)
    {
        old_cnt = tng_data->molecule_cnt_list[index];
        tng_data->molecule_cnt_list[index] = cnt;

        tng_data->n_particles += (cnt-old_cnt) *
                                 tng_data->molecules[index].n_atoms;
    }
    else
    {
        old_cnt = tng_data->current_trajectory_frame_set.molecule_cnt_list[index];
        tng_data->current_trajectory_frame_set.molecule_cnt_list[index] = cnt;

        tng_data->current_trajectory_frame_set.n_particles += (cnt-old_cnt) *
                tng_data->molecules[index].n_atoms;
    }

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_find
                (tng_trajectory_t tng_data,
                 const char *name,
                 int64_t nr,
                 tng_molecule_t *molecule)
{
    int64_t i, n_molecules;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");
    TNG_ASSERT(molecule, "TNG library: molecule must not be a NULL pointer.");

    n_molecules = tng_data->n_molecules;

    for(i = n_molecules - 1; i >= 0; i--)
    {
        *molecule = &tng_data->molecules[i];
        if(name[0] == 0 || strcmp(name, (*molecule)->name) == 0)
        {
            if(nr == -1 || nr == (*molecule)->id)
            {
                return(TNG_SUCCESS);
            }
        }
    }

    *molecule = 0;

    return(TNG_FAILURE);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_of_index_get
                (tng_trajectory_t tng_data,
                 int64_t index,
                 tng_molecule_t *molecule)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(molecule, "TNG library: molecule must not be a NULL pointer.");

    if(index >= tng_data->n_molecules)
    {
        *molecule = 0;
        return(TNG_FAILURE);
    }
    *molecule = &tng_data->molecules[index];
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_system_copy(tng_trajectory_t tng_data_src,
                                                               tng_trajectory_t tng_data_dest)
{
    tng_molecule_t molecule, molecule_temp;
    tng_chain_t chain, chain_temp;
    tng_residue_t residue, residue_temp;
    tng_atom_t atom, atom_temp;
    tng_bond_t bond_temp;
    tng_function_status stat;
    int64_t i, j, k, l, *list_temp;

    TNG_ASSERT(tng_data_src, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(tng_data_dest, "TNG library: Trajectory container not properly setup.");

    for(i = 0; i < tng_data_dest->n_molecules; i++)
    {
        molecule = &tng_data_dest->molecules[i];
        tng_molecule_destroy(tng_data_dest, molecule);
    }

    tng_data_dest->n_molecules = 0;
    tng_data_dest->n_particles = 0;

    molecule_temp = realloc(tng_data_dest->molecules,
                    sizeof(struct tng_molecule) * tng_data_src->n_molecules);
    if(!molecule_temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct tng_molecule) * tng_data_src->n_molecules,
               __FILE__, __LINE__);
        free(tng_data_dest->molecules);
        tng_data_dest->molecules = 0;
        return(TNG_CRITICAL);
    }
    list_temp = realloc(tng_data_dest->molecule_cnt_list,
                                     sizeof(int64_t) * tng_data_src->n_molecules);
    if(!list_temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(int64_t) * tng_data_src->n_molecules,
               __FILE__, __LINE__);
        free(tng_data_dest->molecule_cnt_list);
        tng_data_dest->molecule_cnt_list = 0;
        free(molecule_temp);
        return(TNG_CRITICAL);
    }

    tng_data_dest->molecules = molecule_temp;
    tng_data_dest->molecule_cnt_list = list_temp;

    for(i = 0; i < tng_data_src->n_molecules; i++)
    {
        molecule = &tng_data_src->molecules[i];
        stat = tng_molecule_w_id_add(tng_data_dest, molecule->name, molecule->id,
                                     &molecule_temp);
        if(stat != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot create new molecule to make a copy. %s: %d\n",
                   __FILE__, __LINE__);
            return(stat);
        }
        molecule_temp->quaternary_str = molecule->quaternary_str;
        for(j = 0; j < molecule->n_chains; j++)
        {
            chain = &molecule->chains[j];
            stat = tng_molecule_chain_w_id_add(tng_data_dest, molecule_temp,
                                               chain->name, chain->id,
                                               &chain_temp);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot create new chain to make a copy. %s: %d\n",
                       __FILE__, __LINE__);
                return(stat);
            }
            for(k = 0; k < chain->n_residues; k++)
            {
                residue = &chain->residues[k];
                stat = tng_chain_residue_w_id_add(tng_data_dest, chain_temp,
                                                  residue->name, residue->id,
                                                  &residue_temp);
                if(stat != TNG_SUCCESS)
                {
                    fprintf(stderr, "TNG library: Cannot create new residue to make a copy. %s: %d\n",
                           __FILE__, __LINE__);
                    return(stat);
                }
                for(l = 0; l < residue->n_atoms; l++)
                {
                    atom = &molecule->atoms[residue->atoms_offset + l];
                    stat = tng_residue_atom_w_id_add(tng_data_dest, residue_temp,
                                                     atom->name, atom->atom_type,
                                                     atom->id, &atom_temp);
                    if(stat != TNG_SUCCESS)
                    {
                    fprintf(stderr, "TNG library: Cannot create new atom to make a copy. %s: %d\n",
                           __FILE__, __LINE__);
                        return(stat);
                    }
                }
            }
        }
        molecule_temp->n_bonds = molecule->n_bonds;
        bond_temp = realloc(molecule_temp->bonds, sizeof(struct tng_bond) *
                            molecule->n_bonds);
        if(!bond_temp)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   sizeof(struct tng_bond) * molecule->n_bonds,
                   __FILE__, __LINE__);
            free(molecule_temp->bonds);
            molecule_temp->n_bonds = 0;
            return(TNG_CRITICAL);
        }
        molecule_temp->bonds = bond_temp;
        for(j = 0; j < molecule->n_bonds; j++)
        {
            molecule_temp->bonds[j] = molecule->bonds[j];
        }
        stat = tng_molecule_cnt_set(tng_data_dest, molecule_temp,
                                    tng_data_src->molecule_cnt_list[i]);
        if(stat != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot set molecule count. %s: %d.\n",
                   __FILE__, __LINE__);
            return(stat);
        }
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_num_chains_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 int64_t *n)
{
    (void) tng_data;
    TNG_ASSERT(molecule, "TNG library: molecule must not be NULL");
    TNG_ASSERT(n, "TNG library: n must not be a NULL pointer");

    *n = molecule->n_chains;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_chain_of_index_get
                (tng_trajectory_t tng_data,
                 tng_molecule_t molecule,
                 int64_t index,
                 tng_chain_t *chain)
{
    (void) tng_data;
    TNG_ASSERT(molecule, "TNG library: molecule must not be a NULL pointer.");
    TNG_ASSERT(chain, "TNG library: chain must not be a NULL pointer.");

    if(index >= molecule->n_chains)
    {
        *chain = 0;
        return(TNG_FAILURE);
    }
    *chain = &molecule->chains[index];
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_num_residues_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 int64_t *n)
{
    (void) tng_data;
    TNG_ASSERT(molecule, "TNG library: molecule must not be NULL");
    TNG_ASSERT(n, "TNG library: n must not be a NULL pointer");

    *n = molecule->n_residues;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_residue_of_index_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 const int64_t index,
                 tng_residue_t *residue)
{
    (void) tng_data;
    TNG_ASSERT(molecule, "TNG library: molecule must not be a NULL pointer.");
    TNG_ASSERT(residue, "TNG library: residue must not be a NULL pointer.");

    if(index >= molecule->n_residues)
    {
        *residue = 0;
        return(TNG_FAILURE);
    }
    *residue = &molecule->residues[index];
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_num_atoms_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 int64_t *n)
{
    (void) tng_data;
    TNG_ASSERT(molecule, "TNG library: molecule must not be NULL");
    TNG_ASSERT(n, "TNG library: n must not be a NULL pointer");

    *n = molecule->n_atoms;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_atom_of_index_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 const int64_t index,
                 tng_atom_t *atom)
{
    (void) tng_data;
    TNG_ASSERT(molecule, "TNG library: molecule must not be a NULL pointer.");
    TNG_ASSERT(atom, "TNG library: atom must not be a NULL pointer.");

    if(index >= molecule->n_atoms)
    {
        *atom = 0;
        return(TNG_FAILURE);
    }
    *atom = &molecule->atoms[index];
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_chain_find
                (tng_trajectory_t tng_data,
                 tng_molecule_t molecule,
                 const char *name,
                 int64_t nr,
                 tng_chain_t *chain)
{
    int64_t i, n_chains;
    (void)tng_data;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

    n_chains = molecule->n_chains;

    for(i = n_chains - 1; i >= 0; i--)
    {
        *chain = &molecule->chains[i];
        if(name[0] == 0 || strcmp(name, (*chain)->name) == 0)
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

tng_function_status DECLSPECDLLEXPORT tng_molecule_chain_add
                (tng_trajectory_t tng_data,
                 tng_molecule_t molecule,
                 const char *name,
                 tng_chain_t *chain)
{
    int64_t id;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

    /* Set ID to the ID of the last chain + 1 */
    if(molecule->n_chains)
    {
        id = molecule->chains[molecule->n_chains-1].id + 1;
    }
    else
    {
        id = 1;
    }

    return(tng_molecule_chain_w_id_add(tng_data, molecule, name,
                                       id, chain));
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_chain_w_id_add
                (tng_trajectory_t tng_data,
                 tng_molecule_t molecule,
                 const char *name,
                 const int64_t id,
                 tng_chain_t *chain)
{
    tng_chain_t new_chains;
    tng_function_status stat = TNG_SUCCESS;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

    new_chains = realloc(molecule->chains,
                         sizeof(struct tng_chain) *
                         (molecule->n_chains + 1));

    if(!new_chains)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct tng_chain) * (molecule->n_chains + 1),
               __FILE__, __LINE__);
        free(molecule->chains);
        molecule->chains = 0;
        return(TNG_CRITICAL);
    }

    molecule->chains = new_chains;

    *chain = &new_chains[molecule->n_chains];
    (*chain)->name = 0;

    tng_chain_name_set(tng_data, *chain, name);

    (*chain)->molecule = molecule;
    (*chain)->n_residues = 0;

    molecule->n_chains++;

    (*chain)->id = id;

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_bond_add
                (const tng_trajectory_t tng_data,
                 tng_molecule_t molecule,
                 const int64_t from_atom_id,
                 const int64_t to_atom_id,
                 tng_bond_t *bond)
{
    tng_bond_t new_bonds;
    (void)tng_data;

    new_bonds = realloc(molecule->bonds,
                        sizeof(struct tng_bond) *
                        (molecule->n_bonds + 1));

    if(!new_bonds)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct tng_bond) * (molecule->n_bonds + 1),
               __FILE__, __LINE__);
        *bond = 0;
        free(molecule->bonds);
        molecule->bonds = 0;
        return(TNG_CRITICAL);
    }

    molecule->bonds = new_bonds;

    *bond = &new_bonds[molecule->n_bonds];

    (*bond)->from_atom_id = from_atom_id;
    (*bond)->to_atom_id = to_atom_id;

    molecule->n_bonds++;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_atom_find
                (tng_trajectory_t tng_data,
                 tng_molecule_t molecule,
                 const char *name,
                 int64_t id,
                 tng_atom_t *atom)
{
    int64_t i, n_atoms;
    (void)tng_data;

    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

    n_atoms = molecule->n_atoms;

    for(i = n_atoms - 1; i >= 0; i--)
    {
        *atom = &molecule->atoms[i];
        if(name[0] == 0 || strcmp(name, (*atom)->name) == 0)
        {
            if(id == -1 || id == (*atom)->id)
            {
                return(TNG_SUCCESS);
            }
        }
    }

    *atom = 0;

    return(TNG_FAILURE);
}

tng_function_status tng_chain_name_get(const tng_trajectory_t tng_data,
                                       const tng_chain_t chain,
                                       char *name,
                                       const int max_len)
{
    (void) tng_data;
    TNG_ASSERT(chain, "TNG library: chain must not be NULL");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer");

    strncpy(name, chain->name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(chain->name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_chain_name_set
                (tng_trajectory_t tng_data,
                 tng_chain_t chain,
                 const char *new_name)
{
    unsigned int len;
    (void)tng_data;

    TNG_ASSERT(new_name, "TNG library: new_name must not be a NULL pointer.");

    len = tng_min_i((int)strlen(new_name) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(chain->name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_chain_num_residues_get
                (const tng_trajectory_t tng_data,
                 const tng_chain_t chain,
                 int64_t *n)
{
    (void) tng_data;
    TNG_ASSERT(chain, "TNG library: chain must not be NULL");
    TNG_ASSERT(n, "TNG library: n must not be a NULL pointer");

    *n = chain->n_residues;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_chain_residue_of_index_get
                (const tng_trajectory_t tng_data,
                 const tng_chain_t chain,
                 const int64_t index,
                 tng_residue_t *residue)
{
    (void) tng_data;
    TNG_ASSERT(chain, "TNG library: chain must not be a NULL pointer.");
    TNG_ASSERT(residue, "TNG library: residue must not be a NULL pointer.");

    if(index >= chain->n_residues)
    {
        *residue = 0;
        return(TNG_FAILURE);
    }
    *residue = &chain->residues[index];
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_chain_residue_find
                (tng_trajectory_t tng_data,
                 tng_chain_t chain,
                 const char *name,
                 int64_t id,
                 tng_residue_t *residue)
{
    int64_t i, n_residues;
    (void)tng_data;

    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

    n_residues = chain->n_residues;

    for(i = n_residues - 1; i >= 0; i--)
    {
        *residue = &chain->residues[i];
        if(name[0] == 0 || strcmp(name, (*residue)->name) == 0)
        {
            if(id == -1 || id == (*residue)->id)
            {
                return(TNG_SUCCESS);
            }
        }
    }

    *residue = 0;

    return(TNG_FAILURE);
}

tng_function_status DECLSPECDLLEXPORT tng_chain_residue_add
                (tng_trajectory_t tng_data,
                 tng_chain_t chain,
                 const char *name,
                 tng_residue_t *residue)
{
    int64_t id;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

    /* Set ID to the ID of the last residue + 1 */
    if(chain->n_residues)
    {
        id = chain->residues[chain->n_residues-1].id + 1;
    }
    else
    {
        id = 0;
    }

    return(tng_chain_residue_w_id_add(tng_data, chain, name,
                                      id, residue));
}

tng_function_status DECLSPECDLLEXPORT tng_chain_residue_w_id_add
                (tng_trajectory_t tng_data,
                 tng_chain_t chain,
                 const char *name,
                 const int64_t id,
                 tng_residue_t *residue)
{
    int64_t curr_index;
    tng_residue_t new_residues, temp_residue, last_residue;
    tng_molecule_t molecule = chain->molecule;
    tng_function_status stat = TNG_SUCCESS;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

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
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct tng_residue) * (molecule->n_residues + 1),
               __FILE__, __LINE__);
        free(molecule->residues);
        molecule->residues = 0;
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
                ++temp_residue;
                memmove(temp_residue + 1, temp_residue,
                        last_residue - temp_residue);
            }
        }
    }
    else
    {
        curr_index = molecule->n_residues;
    }

    *residue = &molecule->residues[curr_index + chain->n_residues];

    if(!chain->n_residues)
    {
        chain->residues = *residue;
    }
    else
    {
        chain->residues = &molecule->residues[curr_index];
    }

    (*residue)->name = 0;
    tng_residue_name_set(tng_data, *residue, name);

    (*residue)->chain = chain;
    (*residue)->n_atoms = 0;
    (*residue)->atoms_offset = 0;

    chain->n_residues++;
    molecule->n_residues++;

    (*residue)->id = id;

    return(stat);
}

tng_function_status tng_residue_name_get(const tng_trajectory_t tng_data,
                                         const tng_residue_t residue,
                                         char *name,
                                         const int max_len)
{
    (void) tng_data;
    TNG_ASSERT(residue, "TNG library: residue must not be NULL");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer");

    strncpy(name, residue->name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(residue->name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_residue_name_set(tng_trajectory_t tng_data,
                                                           tng_residue_t residue,
                                                           const char *new_name)
{
    unsigned int len;
    (void)tng_data;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(new_name, "TNG library: new_name must not be a NULL pointer");

    len = tng_min_i((int)strlen(new_name) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(residue->name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_residue_num_atoms_get
                (const tng_trajectory_t tng_data,
                 const tng_residue_t residue,
                 int64_t *n)
{
    (void) tng_data;
    TNG_ASSERT(residue, "TNG library: residue must not be NULL");
    TNG_ASSERT(n, "TNG library: n must not be a NULL pointer");

    *n = residue->n_atoms;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_residue_atom_of_index_get
                (const tng_trajectory_t tng_data,
                 const tng_residue_t residue,
                 const int64_t index,
                 tng_atom_t *atom)
{
    tng_chain_t chain;
    tng_molecule_t molecule;

    (void) tng_data;
    TNG_ASSERT(residue, "TNG library: residue must not be a NULL pointer.");
    TNG_ASSERT(atom, "TNG library: atom must not be a NULL pointer.");

    if(index >= residue->n_atoms)
    {
        *atom = 0;
        return(TNG_FAILURE);
    }
    chain = residue->chain;
    molecule = chain->molecule;

    if(index + residue->atoms_offset >= molecule->n_atoms)
    {
        *atom = 0;
        return(TNG_FAILURE);
    }

    *atom = &molecule->atoms[residue->atoms_offset + index];
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_residue_atom_add
                (tng_trajectory_t tng_data,
                 tng_residue_t residue,
                 const char *atom_name,
                 const char *atom_type,
                 tng_atom_t *atom)
{
    int64_t id;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(atom_name, "TNG library: atom_name must not be a NULL pointer.");
    TNG_ASSERT(atom_type, "TNG library: atom_type must not be a NULL pointer.");

    /* Set ID to the ID of the last atom + 1 */
    if(residue->chain->molecule->n_atoms)
    {
        id = residue->chain->molecule->atoms[residue->chain->molecule->n_atoms-1].id + 1;
    }
    else
    {
        id = 0;
    }

    return(tng_residue_atom_w_id_add(tng_data, residue, atom_name, atom_type,
                                     id, atom));
}

tng_function_status DECLSPECDLLEXPORT tng_residue_atom_w_id_add
                (tng_trajectory_t tng_data,
                 tng_residue_t residue,
                 const char *atom_name,
                 const char *atom_type,
                 const int64_t id,
                 tng_atom_t *atom)
{
    tng_atom_t new_atoms;
    tng_molecule_t molecule = residue->chain->molecule;
    tng_function_status stat = TNG_SUCCESS;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(atom_name, "TNG library: atom_name must not be a NULL pointer.");
    TNG_ASSERT(atom_type, "TNG library: atom_type must not be a NULL pointer.");

    if(!residue->n_atoms)
    {
        residue->atoms_offset = molecule->n_atoms;
    }

    new_atoms = realloc(molecule->atoms,
                        sizeof(struct tng_atom) *
                        (molecule->n_atoms + 1));

    if(!new_atoms)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct tng_atom) * (molecule->n_atoms + 1),
               __FILE__, __LINE__);
        free(molecule->atoms);
        molecule->atoms = 0;
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

    (*atom)->id = id;

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_alloc(const tng_trajectory_t tng_data,
                                                         tng_molecule_t *molecule_p)
{
    *molecule_p = malloc(sizeof(struct tng_molecule));
    if(!*molecule_p)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%lu bytes). %s: %d\n",
               sizeof(struct tng_molecule), __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    tng_molecule_init(tng_data, *molecule_p);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_free(const tng_trajectory_t tng_data,
                                                        tng_molecule_t *molecule_p)
{
    if(!*molecule_p)
    {
        return(TNG_SUCCESS);
    }

    tng_molecule_destroy(tng_data, *molecule_p);

    free(*molecule_p);
    *molecule_p = 0;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_init(const tng_trajectory_t tng_data,
                                                        tng_molecule_t molecule)
{
    (void)tng_data;
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

tng_function_status DECLSPECDLLEXPORT tng_molecule_destroy(const tng_trajectory_t tng_data,
                                                           tng_molecule_t molecule)
{
    int64_t i;
    (void)tng_data;

    if(molecule->name)
    {
        free(molecule->name);
        molecule->name = 0;
    }

    if(molecule->chains)
    {
        for(i = 0; i < molecule->n_chains; i++)
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
        for(i = 0; i < molecule->n_residues; i++)
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
        for(i = 0; i < molecule->n_atoms; i++)
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

tng_function_status DECLSPECDLLEXPORT tng_molecule_name_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    int64_t cnt = 0, i, *molecule_cnt_list = 0;
    tng_molecule_t mol;
    tng_bool found = TNG_FALSE;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

    tng_molecule_cnt_list_get(tng_data, &molecule_cnt_list);

    if(!molecule_cnt_list)
    {
        return(TNG_FAILURE);
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

    if(strlen(mol->name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_id_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 int64_t *id)
{
    int64_t cnt = 0, i, *molecule_cnt_list = 0;
    tng_molecule_t mol;
    tng_bool found = TNG_FALSE;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(id, "TNG library: id must not be a NULL pointer.");

    tng_molecule_cnt_list_get(tng_data, &molecule_cnt_list);

    if(!molecule_cnt_list)
    {
        return(TNG_FAILURE);
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

    *id = mol->id;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molsystem_bonds_get
                (const tng_trajectory_t tng_data,
                 int64_t *n_bonds,
                 int64_t **from_atoms,
                 int64_t **to_atoms)
{
    int64_t atom_cnt = 0, cnt, mol_cnt, i, j, k;
    int64_t from_atom, to_atom, *molecule_cnt_list = 0;
    tng_molecule_t mol;
    tng_bond_t bond;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(n_bonds, "TNG library: n_bonds must not be a NULL pointer.");
    TNG_ASSERT(from_atoms, "TNG library: from_atoms must not be a NULL pointer.");
    TNG_ASSERT(to_atoms, "TNG library: to_atoms must not be a NULL pointer.");

    tng_molecule_cnt_list_get(tng_data, &molecule_cnt_list);

    if(!molecule_cnt_list)
    {
        return(TNG_FAILURE);
    }

    *n_bonds = 0;
    /* First count the total number of bonds to allocate memory */
    for(i = 0; i < tng_data->n_molecules; i++)
    {
        mol = &tng_data->molecules[i];
        mol_cnt = molecule_cnt_list[i];
        *n_bonds += mol_cnt * mol->n_bonds;
    }
    if(*n_bonds == 0)
    {
        return(TNG_SUCCESS);
    }

    *from_atoms = malloc(sizeof(int64_t) * (*n_bonds));
    if(!*from_atoms)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(int64_t) * (*n_bonds), __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }
    *to_atoms = malloc(sizeof(int64_t) * (*n_bonds));
    if(!*to_atoms)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(int64_t) * (*n_bonds), __FILE__, __LINE__);
        free(*from_atoms);
        *from_atoms = 0;
        return(TNG_CRITICAL);
    }

    cnt = 0;
    for(i = 0; i < tng_data->n_molecules; i++)
    {
        mol = &tng_data->molecules[i];
        mol_cnt = molecule_cnt_list[i];
        for(j = 0; j < mol_cnt; j++)
        {
            for(k = 0; k < mol->n_bonds; k++)
            {
                bond = &mol->bonds[k];
                from_atom = atom_cnt + bond->from_atom_id;
                to_atom = atom_cnt + bond->to_atom_id;
                (*from_atoms)[cnt] = from_atom;
                (*to_atoms)[cnt++] = to_atom;
            }
            atom_cnt += mol->n_atoms;
        }
    }

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_chain_name_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    int64_t cnt = 0, i, *molecule_cnt_list = 0;
    tng_molecule_t mol;
    tng_atom_t atom;
    tng_bool found = TNG_FALSE;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

    tng_molecule_cnt_list_get(tng_data, &molecule_cnt_list);

    if(!molecule_cnt_list)
    {
        return(TNG_FAILURE);
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

    if(strlen(atom->residue->chain->name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_residue_name_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    int64_t cnt = 0, i, *molecule_cnt_list = 0;
    tng_molecule_t mol;
    tng_atom_t atom;
    tng_bool found = TNG_FALSE;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

    tng_molecule_cnt_list_get(tng_data, &molecule_cnt_list);

    if(!molecule_cnt_list)
    {
        return(TNG_FAILURE);
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

    if(strlen(atom->residue->name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_residue_id_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 int64_t *id)
{
    int64_t cnt = 0, i, *molecule_cnt_list = 0;
    tng_molecule_t mol;
    tng_atom_t atom;
    tng_bool found = TNG_FALSE;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(id, "TNG library: id must not be a NULL pointer.");

    tng_molecule_cnt_list_get(tng_data, &molecule_cnt_list);

    if(!molecule_cnt_list)
    {
        return(TNG_FAILURE);
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

    *id = atom->residue->id;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_global_residue_id_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 int64_t *id)
{
    int64_t cnt = 0, i, offset = 0, *molecule_cnt_list = 0;
    tng_molecule_t mol;
    tng_atom_t atom;
    tng_bool found = TNG_FALSE;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(id, "TNG library: id must not be a NULL pointer.");

    tng_molecule_cnt_list_get(tng_data, &molecule_cnt_list);

    if(!molecule_cnt_list)
    {
        return(TNG_FAILURE);
    }

    for(i = 0; i < tng_data->n_molecules; i++)
    {
        mol = &tng_data->molecules[i];
        if(cnt + mol->n_atoms * molecule_cnt_list[i] - 1 < nr)
        {
            cnt += mol->n_atoms * molecule_cnt_list[i];
            offset += mol->n_residues * molecule_cnt_list[i];
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

    offset += mol->n_residues * ((nr - cnt) / mol->n_atoms);

    *id = atom->residue->id + offset;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_atom_name_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 int max_len)
{
    int64_t cnt = 0, i, *molecule_cnt_list = 0;
    tng_molecule_t mol;
    tng_atom_t atom;
    tng_bool found = TNG_FALSE;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

    tng_molecule_cnt_list_get(tng_data, &molecule_cnt_list);

    if(!molecule_cnt_list)
    {
        return(TNG_FAILURE);
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

    if(strlen(atom->name) > (unsigned int)max_len - 1)
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
    int64_t cnt = 0, i, *molecule_cnt_list = 0;
    tng_molecule_t mol;
    tng_atom_t atom;
    tng_bool found = TNG_FALSE;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(type, "TNG library: type must not be a NULL pointer.");

    tng_molecule_cnt_list_get(tng_data, &molecule_cnt_list);

    if(!molecule_cnt_list)
    {
        return(TNG_FAILURE);
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

    if(strlen(atom->atom_type) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_particle_mapping_add
                (tng_trajectory_t tng_data,
                 const int64_t num_first_particle,
                 const int64_t n_particles,
                 const int64_t *mapping_table)
{
    int64_t i;
    tng_particle_mapping_t mapping;
    tng_trajectory_frame_set_t frame_set;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    frame_set = &tng_data->current_trajectory_frame_set;

    /* Sanity check of the particle ranges. Split into multiple if
     * statements for improved readability */
    for(i = 0; i < frame_set->n_mapping_blocks; i++)
    {
        mapping = &frame_set->mappings[i];
        if(num_first_particle >= mapping->num_first_particle &&
           num_first_particle < mapping->num_first_particle +
                                   mapping->n_particles)
        {
            fprintf(stderr, "TNG library: Particle mapping overlap. %s: %d\n", __FILE__, __LINE__);
            return(TNG_FAILURE);
        }
        if(num_first_particle + n_particles >=
           mapping->num_first_particle &&
           num_first_particle + n_particles <
           mapping->num_first_particle + mapping->n_particles)
        {
            fprintf(stderr, "TNG library: Particle mapping overlap. %s: %d\n", __FILE__, __LINE__);
            return(TNG_FAILURE);
        }
        if(mapping->num_first_particle >= num_first_particle &&
           mapping->num_first_particle < num_first_particle +
                                            n_particles)
        {
            fprintf(stderr, "TNG library: Particle mapping overlap. %s: %d\n", __FILE__, __LINE__);
            return(TNG_FAILURE);
        }
        if(mapping->num_first_particle + mapping->n_particles >
           num_first_particle &&
           mapping->num_first_particle + mapping->n_particles <
           num_first_particle + n_particles)
        {
            fprintf(stderr, "TNG library: Particle mapping overlap. %s: %d\n", __FILE__, __LINE__);
            return(TNG_FAILURE);
        }
    }

    frame_set->n_mapping_blocks++;

    mapping = realloc(frame_set->mappings, sizeof(struct tng_particle_mapping) *
                      frame_set->n_mapping_blocks);

    if(!mapping)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct tng_particle_mapping)*frame_set->n_mapping_blocks,
               __FILE__, __LINE__);
        free(frame_set->mappings);
        frame_set->mappings = 0;
        return(TNG_CRITICAL);
    }
    frame_set->mappings = mapping;

    frame_set->mappings[frame_set->n_mapping_blocks - 1].num_first_particle = num_first_particle;
    frame_set->mappings[frame_set->n_mapping_blocks - 1].n_particles = n_particles;

    frame_set->mappings[frame_set->n_mapping_blocks - 1].real_particle_numbers = malloc(sizeof(int64_t) * n_particles);
    if(!frame_set->mappings[frame_set->n_mapping_blocks - 1].real_particle_numbers)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(int64_t) * n_particles, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    for(i=0; i<n_particles; i++)
    {
        frame_set->mappings[frame_set->n_mapping_blocks - 1].real_particle_numbers[i] = mapping_table[i];
    }

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_particle_mapping_free(tng_trajectory_t tng_data)
{
    tng_trajectory_frame_set_t frame_set;
    tng_particle_mapping_t mapping;
    int64_t i;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    frame_set = &tng_data->current_trajectory_frame_set;

    if(frame_set->n_mapping_blocks && frame_set->mappings)
    {
        for(i = 0; i < frame_set->n_mapping_blocks; i++)
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

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_trajectory_init(tng_trajectory_t *tng_data_p)
{
    time_t seconds;
    tng_trajectory_frame_set_t frame_set;
    tng_trajectory_t tng_data;

    *tng_data_p = malloc(sizeof(struct tng_trajectory));
    if(!*tng_data_p)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%lu bytes). %s: %d\n",
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
        fprintf(stderr, "TNG library: Cannot get time. %s: %d\n", __FILE__, __LINE__);
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
    tng_data->medium_stride_length = 100;
    tng_data->long_stride_length = 10000;

    tng_data->time_per_frame = -1;

    tng_data->n_particle_data_blocks = 0;
    tng_data->n_data_blocks = 0;

    tng_data->non_tr_particle_data = 0;
    tng_data->non_tr_data = 0;

    tng_data->compress_algo_pos = 0;
    tng_data->compress_algo_vel = 0;
    tng_data->compression_precision = 1000;
    tng_data->distance_unit_exponential = -9;

    frame_set->first_frame = -1;
    frame_set->n_mapping_blocks = 0;
    frame_set->mappings = 0;
    frame_set->molecule_cnt_list = 0;

    frame_set->n_particle_data_blocks = 0;
    frame_set->n_data_blocks = 0;

    frame_set->tr_particle_data = 0;
    frame_set->tr_data = 0;

    frame_set->n_written_frames = 0;
    frame_set->n_unwritten_frames = 0;

    frame_set->next_frame_set_file_pos = -1;
    frame_set->prev_frame_set_file_pos = -1;
    frame_set->medium_stride_next_frame_set_file_pos = -1;
    frame_set->medium_stride_prev_frame_set_file_pos = -1;
    frame_set->long_stride_next_frame_set_file_pos = -1;
    frame_set->long_stride_prev_frame_set_file_pos = -1;

    frame_set->first_frame_time = -1;

    tng_data->n_molecules = 0;
    tng_data->molecules = 0;
    tng_data->molecule_cnt_list = 0;
    tng_data->n_particles = 0;

    {
      /* Check the endianness of the computer */
      static int32_t endianness_32 = 0x01234567;
      /* 0x01234567 */
      if ( *(const unsigned char*)&endianness_32 == 0x01 )
        {
          tng_data->endianness_32 = TNG_BIG_ENDIAN_32;
        }

      /* 0x67452301 */
      else if( *(const unsigned char*)&endianness_32 == 0x67 )
        {
          tng_data->endianness_32 = TNG_LITTLE_ENDIAN_32;

        }

      /* 0x45670123 */
      else if ( *(const unsigned char*)&endianness_32 == 0x45 )
        {
          tng_data->endianness_32 = TNG_BYTE_PAIR_SWAP_32;
        }
    }
    {
      static int64_t endianness_64 = 0x0123456789ABCDEFLL;
      /* 0x0123456789ABCDEF */
      if ( *(const unsigned char*)&endianness_64 == 0x01 )
        {
          tng_data->endianness_64 = TNG_BIG_ENDIAN_64;
        }

      /* 0xEFCDAB8967452301 */
      else if ( *(const unsigned char*)&endianness_64 == 0xEF )
        {
          tng_data->endianness_64 = TNG_LITTLE_ENDIAN_64;
        }

      /* 0x89ABCDEF01234567 */
      else if ( *(const unsigned char*)&endianness_64 == 0x89 )
        {
          tng_data->endianness_64 = TNG_QUAD_SWAP_64;
        }

      /* 0x45670123CDEF89AB */
      else if ( *(const unsigned char*)&endianness_64 == 0x45 )
        {
          tng_data->endianness_64 = TNG_BYTE_PAIR_SWAP_64;
        }

      /* 0x23016745AB89EFCD */
      else if ( *(const unsigned char*)&endianness_64 == 0x23 )
        {
          tng_data->endianness_64 = TNG_BYTE_SWAP_64;
        }
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

tng_function_status DECLSPECDLLEXPORT tng_trajectory_destroy(tng_trajectory_t *tng_data_p)
{
    int64_t i, j, k, l;
    int64_t n_particles, n_values_per_frame;
    tng_trajectory_t tng_data = *tng_data_p;
    tng_trajectory_frame_set_t frame_set;

    if(!*tng_data_p)
    {
        return(TNG_SUCCESS);
    }

    frame_set = &tng_data->current_trajectory_frame_set;

    if(tng_data->input_file_path)
    {
        free(tng_data->input_file_path);
        tng_data->input_file_path = 0;
    }

    if(tng_data->input_file)
    {
        if(tng_data->output_file == tng_data->input_file)
        {
            tng_frame_set_finalize(tng_data, TNG_USE_HASH);
            tng_data->output_file = 0;
        }
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

    tng_frame_set_particle_mapping_free(tng_data);

    if(frame_set->molecule_cnt_list)
    {
        free(frame_set->molecule_cnt_list);
        frame_set->molecule_cnt_list = 0;
    }

    if(tng_data->var_num_atoms_flag)
    {
        n_particles = frame_set->n_particles;
    }
    else
    {
        n_particles = tng_data->n_particles;
    }

    if(tng_data->non_tr_particle_data)
    {
        for(i = 0; i < tng_data->n_particle_data_blocks; i++)
        {
            if(tng_data->non_tr_particle_data[i].values)
            {
                free(tng_data->non_tr_particle_data[i].values);
                tng_data->non_tr_particle_data[i].values = 0;
            }

            if(tng_data->non_tr_particle_data[i].strings)
            {
                n_values_per_frame = tng_data->non_tr_particle_data[i].
                                     n_values_per_frame;
                if(tng_data->non_tr_particle_data[i].strings[0])
                {
                    for(j = 0; j < n_particles; j++)
                    {
                        if(tng_data->non_tr_particle_data[i].strings[0][j])
                        {
                            for(k = 0; k < n_values_per_frame; k++)
                            {
                                if(tng_data->non_tr_particle_data[i].
                                   strings[0][j][k])
                                {
                                    free(tng_data->non_tr_particle_data[i].
                                         strings[0][j][k]);
                                    tng_data->non_tr_particle_data[i].
                                    strings[0][j][k] = 0;
                                }
                            }
                            free(tng_data->non_tr_particle_data[i].
                                 strings[0][j]);
                            tng_data->non_tr_particle_data[i].strings[0][j] = 0;
                        }
                    }
                    free(tng_data->non_tr_particle_data[i].strings[0]);
                    tng_data->non_tr_particle_data[i].strings[0] = 0;
                }
                free(tng_data->non_tr_particle_data[i].strings);
                tng_data->non_tr_particle_data[i].strings = 0;
            }

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
        for(i = 0; i < tng_data->n_data_blocks; i++)
        {
            if(tng_data->non_tr_data[i].values)
            {
                free(tng_data->non_tr_data[i].values);
                tng_data->non_tr_data[i].values = 0;
            }

            if(tng_data->non_tr_data[i].strings)
            {
                n_values_per_frame = tng_data->non_tr_data[i].
                                     n_values_per_frame;
                if(tng_data->non_tr_data[i].strings[0])
                {
                    for(j = 0; j < n_values_per_frame; j++)
                    {
                        if(tng_data->non_tr_data[i].strings[0][j])
                        {
                            free(tng_data->non_tr_data[i].strings[0][j]);
                            tng_data->non_tr_data[i].strings[0][j] = 0;
                        }
                    }
                    free(tng_data->non_tr_data[i].strings[0]);
                    tng_data->non_tr_data[i].strings[0] = 0;
                }
                free(tng_data->non_tr_data[i].strings);
                tng_data->non_tr_data[i].strings = 0;
            }

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

    if(tng_data->compress_algo_pos)
    {
        free(tng_data->compress_algo_pos);
        tng_data->compress_algo_pos = 0;
    }
    if(tng_data->compress_algo_vel)
    {
        free(tng_data->compress_algo_vel);
        tng_data->compress_algo_vel = 0;
    }

    if(frame_set->tr_particle_data)
    {
        for(i = 0; i < frame_set->n_particle_data_blocks; i++)
        {
            if(frame_set->tr_particle_data[i].values)
            {
                free(frame_set->tr_particle_data[i].values);
                frame_set->tr_particle_data[i].values = 0;
            }

            if(frame_set->tr_particle_data[i].strings)
            {
                n_values_per_frame = frame_set->tr_particle_data[i].
                                     n_values_per_frame;
                for(j = 0; j < frame_set->tr_particle_data[i].n_frames; j++)
                {
                    if(frame_set->tr_particle_data[i].strings[j])
                    {
                        for(k = 0; k < n_particles; k++)
                        {
                            if(frame_set->tr_particle_data[i].
                                strings[j][k])
                            {
                                for(l = 0; l < n_values_per_frame; l++)
                                {
                                    if(frame_set->tr_particle_data[i].
                                        strings[j][k][l])
                                    {
                                        free(frame_set->tr_particle_data[i].
                                                strings[j][k][l]);
                                        frame_set->tr_particle_data[i].
                                        strings[j][k][l] = 0;
                                    }
                                }
                                free(frame_set->tr_particle_data[i].
                                        strings[j][k]);
                                frame_set->tr_particle_data[i].
                                strings[j][k] = 0;
                            }
                        }
                        free(frame_set->tr_particle_data[i].strings[j]);
                        frame_set->tr_particle_data[i].strings[j] = 0;
                    }
                }
                free(frame_set->tr_particle_data[i].strings);
                frame_set->tr_particle_data[i].strings = 0;
            }

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
        for(i = 0; i < frame_set->n_data_blocks; i++)
        {
            if(frame_set->tr_data[i].values)
            {
                free(frame_set->tr_data[i].values);
                frame_set->tr_data[i].values = 0;
            }

            if(frame_set->tr_data[i].strings)
            {
                n_values_per_frame = frame_set->tr_data[i].
                                     n_values_per_frame;
                for(j = 0; j < frame_set->tr_data[i].n_frames; j++)
                {
                    if(frame_set->tr_data[i].strings[j])
                    {
                        for(k = 0; k < n_values_per_frame; k++)
                        {
                            if(frame_set->tr_data[i].strings[j][k])
                            {
                                free(frame_set->tr_data[i].strings[j][k]);
                                frame_set->tr_data[i].strings[j][k] = 0;
                            }
                        }
                        free(frame_set->tr_data[i].strings[j]);
                        frame_set->tr_data[i].strings[j] = 0;
                    }
                }
                free(frame_set->tr_data[i].strings);
                frame_set->tr_data[i].strings = 0;
            }

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
        for(i = 0; i < tng_data->n_molecules; i++)
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

tng_function_status DECLSPECDLLEXPORT tng_trajectory_init_from_src(tng_trajectory_t src,
                                                 tng_trajectory_t *dest_p)
{
    tng_trajectory_frame_set_t frame_set;
    tng_trajectory_t dest;

    TNG_ASSERT(src != 0, "TNG library: Source trajectory must not be NULL.");

    *dest_p = malloc(sizeof(struct tng_trajectory));
    if(!*dest_p)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%lu bytes). %s: %d\n",
               sizeof(struct tng_trajectory), __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    dest = *dest_p;

    frame_set = &dest->current_trajectory_frame_set;

    dest->input_file_path = malloc(strlen(src->input_file_path) + 1);
    if(!dest->input_file_path)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
               (int)strlen(src->input_file_path) + 1, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }
    strcpy(dest->input_file_path, src->input_file_path);
    dest->input_file = 0;
    dest->input_file_len = src->input_file_len;
    dest->output_file_path = malloc(strlen(src->output_file_path) + 1);
    if(!dest->output_file_path)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
               (int)strlen(src->output_file_path) + 1, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }
    strcpy(dest->output_file_path, src->output_file_path);
    dest->output_file = 0;

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
    dest->medium_stride_length = src->medium_stride_length;
    dest->long_stride_length = src->long_stride_length;

    dest->time_per_frame = src->time_per_frame;

    /* Currently the non trajectory data blocks are not copied since it
     * can lead to problems when freeing memory in a parallel block. */
    dest->n_particle_data_blocks = 0;
    dest->n_data_blocks = 0;
    dest->non_tr_particle_data = 0;
    dest->non_tr_data = 0;

    dest->compress_algo_pos = 0;
    dest->compress_algo_vel = 0;
    dest->distance_unit_exponential = -9;
    dest->compression_precision = 1000;

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
    frame_set->first_frame = -1;

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

tng_function_status DECLSPECDLLEXPORT tng_input_file_get(const tng_trajectory_t tng_data,
                                       char *file_name, const int max_len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(file_name, "TNG library: file_name must not be a NULL pointer");

    strncpy(file_name, tng_data->input_file_path, max_len - 1);
    file_name[max_len - 1] = 0;

    if(strlen(tng_data->input_file_path) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_input_file_set(tng_trajectory_t tng_data,
                                                         const char *file_name)
{
    unsigned int len;
    char *temp;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(file_name, "TNG library: file_name must not be a NULL pointer");


    if(tng_data->input_file_path && strcmp(tng_data->input_file_path,
                                           file_name) == 0)
    {
        return(TNG_SUCCESS);
    }

    if(tng_data->input_file)
    {
        fclose(tng_data->input_file);
    }

    len = tng_min_i((int)strlen(file_name) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->input_file_path, len);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->input_file_path);
        tng_data->input_file_path = 0;
        return(TNG_CRITICAL);
    }
    tng_data->input_file_path = temp;

    strncpy(tng_data->input_file_path, file_name, len);

    return(tng_input_file_init(tng_data));
}

tng_function_status tng_output_file_get(const tng_trajectory_t tng_data,
                                       char *file_name, const int max_len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(file_name, "TNG library: file_name must not be a NULL pointer");

    strncpy(file_name, tng_data->output_file_path, max_len - 1);
    file_name[max_len - 1] = 0;

    if(strlen(tng_data->output_file_path) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_output_file_set(tng_trajectory_t tng_data,
                                                          const char *file_name)
{
    int len;
    char *temp;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(file_name, "TNG library: file_name must not be a NULL pointer");

    if(tng_data->output_file_path &&
       strcmp(tng_data->output_file_path, file_name) == 0)
    {
        return(TNG_SUCCESS);
    }

    if(tng_data->output_file)
    {
        fclose(tng_data->output_file);
    }

    len = tng_min_i((int)strlen(file_name) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->output_file_path, len);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->output_file_path);
        tng_data->output_file_path = 0;
        return(TNG_CRITICAL);
    }
    tng_data->output_file_path = temp;

    strncpy(tng_data->output_file_path, file_name, len);

    return(tng_output_file_init(tng_data));
}

tng_function_status DECLSPECDLLEXPORT tng_output_append_file_set
                (tng_trajectory_t tng_data,
                 const char *file_name)
{
    int len;
    char *temp;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(file_name, "TNG library: file_name must not be a NULL pointer");

    if(tng_data->output_file_path &&
       strcmp(tng_data->output_file_path, file_name) == 0)
    {
        return(TNG_SUCCESS);
    }

    if(tng_data->output_file)
    {
        fclose(tng_data->output_file);
    }

    len = tng_min_i((int)strlen(file_name) + 1, TNG_MAX_STR_LEN);
    temp = realloc(tng_data->output_file_path, len);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        free(tng_data->output_file_path);
        tng_data->output_file_path = 0;
        return(TNG_CRITICAL);
    }
    tng_data->output_file_path = temp;

    strncpy(tng_data->output_file_path, file_name, len);

    tng_data->output_file = fopen(tng_data->output_file_path, "rb+");
    if(!tng_data->output_file)
    {
        fprintf(stderr, "TNG library: Cannot open file %s. %s: %d\n",
                tng_data->output_file_path, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }
    tng_data->input_file = tng_data->output_file;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_output_file_endianness_get
                (const tng_trajectory_t tng_data, tng_file_endianness *endianness)
{
    tng_endianness_32 end_32;
    tng_endianness_64 end_64;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(endianness, "TNG library: endianness must not be a NULL pointer");

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
        end_32 = (tng_endianness_32)tng_data->endianness_32;
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
        end_64 = (tng_endianness_64)tng_data->endianness_64;
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

tng_function_status DECLSPECDLLEXPORT tng_output_file_endianness_set
                (tng_trajectory_t tng_data,
                 const tng_file_endianness endianness)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

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

tng_function_status DECLSPECDLLEXPORT tng_first_program_name_get
                    (const tng_trajectory_t tng_data,
                     char *name, const int max_len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer");

    strncpy(name, tng_data->first_program_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->first_program_name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_first_program_name_set(tng_trajectory_t tng_data,
                                                                 const char *new_name)
{
    unsigned int len;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(new_name, "TNG library: new_name must not be a NULL pointer");

    len = tng_min_i((int)strlen(new_name) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->first_program_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_last_program_name_get
                    (const tng_trajectory_t tng_data,
                     char *name, const int max_len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer");

    strncpy(name, tng_data->last_program_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->last_program_name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_last_program_name_set
                    (tng_trajectory_t tng_data,
                     const char *new_name)
{
    unsigned int len;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(new_name, "TNG library: new_name must not be a NULL pointer");

    len = tng_min_i((int)strlen(new_name) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->last_program_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_first_user_name_get
                    (const tng_trajectory_t tng_data,
                     char *name, const int max_len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer");

    strncpy(name, tng_data->first_user_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->first_user_name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_first_user_name_set
                    (tng_trajectory_t tng_data,
                     const char *new_name)
{
    unsigned int len;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(new_name, "TNG library: new_name must not be a NULL pointer");

    len = tng_min_i((int)strlen(new_name) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->first_user_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_last_user_name_get
                    (const tng_trajectory_t tng_data,
                     char *name, const int max_len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer");

    strncpy(name, tng_data->last_user_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->last_user_name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_last_user_name_set
                    (tng_trajectory_t tng_data,
                     const char *new_name)
{
    unsigned int len;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(new_name, "TNG library: new_name must not be a NULL pointer");

    len = tng_min_i((int)strlen(new_name) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->last_user_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_first_computer_name_get
                    (const tng_trajectory_t tng_data,
                     char *name, const int max_len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer");

    strncpy(name, tng_data->first_computer_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->first_computer_name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_first_computer_name_set
                    (tng_trajectory_t tng_data,
                     const char *new_name)
{
    unsigned int len;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(new_name, "TNG library: new_name must not be a NULL pointer");

    len = tng_min_i((int)strlen(new_name) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->first_computer_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_last_computer_name_get
                    (const tng_trajectory_t tng_data,
                     char *name, const int max_len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer");

    strncpy(name, tng_data->last_computer_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->last_computer_name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_last_computer_name_set
                    (tng_trajectory_t tng_data,
                     const char *new_name)
{
    unsigned int len;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(new_name, "TNG library: new_name must not be a NULL pointer");

    len = tng_min_i((int)strlen(new_name) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->last_computer_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_first_signature_get
                    (const tng_trajectory_t tng_data,
                     char *signature, const int max_len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(signature, "TNG library: signature must not be a NULL pointer");

    strncpy(signature, tng_data->first_pgp_signature, max_len - 1);
    signature[max_len - 1] = 0;

    if(strlen(tng_data->first_pgp_signature) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_first_signature_set
                    (tng_trajectory_t tng_data,
                     const char *signature)
{
    unsigned int len;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(signature, "TNG library: signature must not be a NULL pointer");

    len = tng_min_i((int)strlen(signature) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->first_pgp_signature, signature, len);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_last_signature_get
                    (const tng_trajectory_t tng_data,
                     char *signature, const int max_len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(signature, "TNG library: signature must not be a NULL pointer");

    strncpy(signature, tng_data->last_pgp_signature, max_len - 1);
    signature[max_len - 1] = 0;

    if(strlen(tng_data->last_pgp_signature) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_last_signature_set
                    (tng_trajectory_t tng_data,
                     const char *signature)
{
    unsigned int len;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(signature, "TNG library: signature must not be a NULL pointer");

    len = tng_min_i((int)strlen(signature) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->last_pgp_signature, signature, len);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_forcefield_name_get
                    (const tng_trajectory_t tng_data,
                     char *name, const int max_len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer");

    strncpy(name, tng_data->forcefield_name, max_len - 1);
    name[max_len - 1] = 0;

    if(strlen(tng_data->forcefield_name) > (unsigned int)max_len - 1)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_forcefield_name_set
                    (tng_trajectory_t tng_data,
                     const char *new_name)
{
    unsigned int len;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(new_name, "TNG library: new_name must not be a NULL pointer");

    len = tng_min_i((int)strlen(new_name) + 1, TNG_MAX_STR_LEN);

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
            fprintf(stderr, "TNG library: Cannot allocate memory (%ud bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
    }

    strncpy(tng_data->forcefield_name, new_name, len);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_medium_stride_length_get
                    (const tng_trajectory_t tng_data,
                     int64_t *len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(len, "TNG library: len must not be a NULL pointer");

    *len = tng_data->medium_stride_length;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_medium_stride_length_set
                    (tng_trajectory_t tng_data,
                     const int64_t len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    if(len >= tng_data->long_stride_length)
    {
        return(TNG_FAILURE);
    }
    tng_data->medium_stride_length = len;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_long_stride_length_get
                (const tng_trajectory_t tng_data,
                 int64_t *len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(len, "TNG library: len must not be a NULL pointer");

    *len = tng_data->long_stride_length;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_long_stride_length_set
                (tng_trajectory_t tng_data,
                 const int64_t len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    if(len <= tng_data->medium_stride_length)
    {
        return(TNG_FAILURE);
    }
    tng_data->long_stride_length = len;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_time_per_frame_get
                (const tng_trajectory_t tng_data,
                 double *time)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(time, "TNG library: time must not be a NULL pointer");

    *time = tng_data->time_per_frame;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_time_per_frame_set
                (tng_trajectory_t tng_data,
                 const double time)
{
    tng_trajectory_frame_set_t frame_set;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(time >= 0, "TNG library: The time per frame must be >= 0.");

    if(fabs(time - tng_data->time_per_frame) < 0.00001)
    {
        return(TNG_SUCCESS);
    }

    frame_set = &tng_data->current_trajectory_frame_set;

    /* If the current frame set is not finished write it to disk before
       changing time per frame. */
    if(tng_data->time_per_frame > 0 && frame_set->n_unwritten_frames > 0)
    {
        frame_set->n_frames = frame_set->n_unwritten_frames;
        tng_frame_set_write(tng_data, TNG_USE_HASH);
    }
    tng_data->time_per_frame = time;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_input_file_len_get
                    (const tng_trajectory_t tng_data,
                     int64_t *len)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(len, "TNG library: len must not be a NULL pointer");

    *len = tng_data->input_file_len;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_num_frames_get
                    (const tng_trajectory_t tng_data,
                     int64_t *n)
{
    tng_gen_block_t block;
    tng_function_status stat;
    long file_pos;
    int64_t last_file_pos, first_frame, n_frames;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(tng_data->input_file, "TNG library: An input file must be open to find the next frame set");
    TNG_ASSERT(n, "TNG library: n must not be a NULL pointer");

    file_pos = ftell(tng_data->input_file);
    last_file_pos = tng_data->last_trajectory_frame_set_input_file_pos;

    if(last_file_pos <= 0)
    {
        return(TNG_FAILURE);
    }

    tng_block_init(&block);
    fseek(tng_data->input_file,
          (long)last_file_pos,
          SEEK_SET);
    /* Read block headers first to see that a frame set block is found. */
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n", last_file_pos,
                __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_FAILURE);
    }
    tng_block_destroy(&block);

    if(fread(&first_frame, sizeof(int64_t), 1, tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Cannot read first frame of frame set. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }
    if(fread(&n_frames, sizeof(int64_t), 1, tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Cannot read n frames of frame set. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }
    fseek(tng_data->input_file, file_pos, SEEK_SET);

    *n = first_frame + n_frames;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_compression_precision_get
                (const tng_trajectory_t tng_data,
                 double *precision)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    *precision = tng_data->compression_precision;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_compression_precision_set
                (tng_trajectory_t tng_data,
                 const double precision)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    tng_data->compression_precision = precision;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_implicit_num_particles_set
                (tng_trajectory_t tng_data,
                 const int64_t n)
{
    tng_molecule_t mol;
    tng_chain_t chain;
    tng_residue_t res;
    tng_atom_t atom;
    tng_function_status stat;
    int64_t diff, n_mod, n_impl;

    TNG_ASSERT(n >= 0, "TNG library: The number of molecules must be >= 0");

    diff = n - tng_data->n_particles;

    stat = tng_molecule_find(tng_data, "TNG_IMPLICIT_MOL", -1, &mol);
    if(stat == TNG_SUCCESS)
    {
        if(tng_molecule_cnt_get(tng_data, mol, &n_impl) != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot get the number of implicit molecules. %s: %d\n",
                    __FILE__, __LINE__);
            return(TNG_FAILURE);
        }
        diff -= n_impl * mol->n_atoms;
    }

    if(diff == 0)
    {
        if(stat == TNG_SUCCESS)
        {
            stat = tng_molecule_cnt_set(tng_data, mol, 0);
            return(stat);
        }
        return(TNG_SUCCESS);
    }
    else if(diff < 0)
    {
        fprintf(stderr, "TNG library: Already more actual particles than requested implicit ");
        fprintf(stderr, "particle count.\n");
        fprintf(stderr, "TNG library: Cannot set implicit particle count. %s: %d\n",
                __FILE__, __LINE__);
        /* FIXME: Should we set the count of all other molecules to 0 and add
         * implicit molecules? */
        return(TNG_FAILURE);
    }
    if(stat != TNG_SUCCESS)
    {
        stat = tng_molecule_add(tng_data,
                                "TNG_IMPLICIT_MOL",
                                &mol);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
        stat = tng_molecule_chain_add(tng_data, mol, "", &chain);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
        stat = tng_chain_residue_add(tng_data, chain, "", &res);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
        stat = tng_residue_atom_add(tng_data, res, "", "", &atom);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
    }
    else
    {
        if(mol->n_atoms > 1)
        {
            n_mod = diff % mol->n_atoms;
            if(n_mod != 0)
            {
                fprintf(stderr, "TNG library: Number of atoms in implicit molecule ");
                fprintf(stderr, "not compatible with requested implicit particle cnt.\n");
                fprintf(stderr, "TNG library: Cannot set implicit particle count. %s: %d\n",
                        __FILE__, __LINE__);
                return(TNG_FAILURE);
            }
            diff /= mol->n_atoms;
        }
    }
    stat = tng_molecule_cnt_set(tng_data, mol, diff);

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_num_particles_get
                (const tng_trajectory_t tng_data,
                 int64_t *n)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(n, "TNG library: n must not be a NULL pointer");

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

tng_function_status DECLSPECDLLEXPORT tng_num_particles_variable_get
                (const tng_trajectory_t tng_data,
                 char *variable)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(variable, "TNG library: variable must not be a NULL pointer");

    *variable = tng_data->var_num_atoms_flag;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_num_molecule_types_get
                    (const tng_trajectory_t tng_data,
                     int64_t *n)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(n, "TNG library: n must not be a NULL pointer");

    *n = tng_data->n_molecules;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_num_molecules_get
                    (const tng_trajectory_t tng_data,
                     int64_t *n)
{
    int64_t *cnt_list = 0, cnt = 0, i;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(n, "TNG library: n must not be a NULL pointer");

    tng_molecule_cnt_list_get(tng_data, &cnt_list);

    if(!cnt_list)
    {
        return(TNG_FAILURE);
    }

    for(i = 0; i < tng_data->n_molecules; i++)
    {
        cnt += cnt_list[i];
    }

    *n = cnt;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_molecule_cnt_list_get
                (const tng_trajectory_t tng_data,
                 int64_t **mol_cnt_list)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    if(tng_data->var_num_atoms_flag)
    {
        *mol_cnt_list = tng_data->current_trajectory_frame_set.
                       molecule_cnt_list;
    }
    else
    {
        *mol_cnt_list = tng_data->molecule_cnt_list;
    }
    if(*mol_cnt_list == 0)
    {
        return(TNG_FAILURE);
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_distance_unit_exponential_get
                (const tng_trajectory_t tng_data,
                 int64_t *exp)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(exp, "TNG library: exp must not be a NULL pointer");

    *exp = tng_data->distance_unit_exponential;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_distance_unit_exponential_set
                (const tng_trajectory_t tng_data,
                 const int64_t exp)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    tng_data->distance_unit_exponential = exp;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_num_frames_per_frame_set_get
                (const tng_trajectory_t tng_data,
                 int64_t *n)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(n, "TNG library: n must not be a NULL pointer");

    *n = tng_data->frame_set_n_frames;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_num_frames_per_frame_set_set
                (const tng_trajectory_t tng_data,
                 const int64_t n)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    tng_data->frame_set_n_frames = n;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_num_frame_sets_get
                (const tng_trajectory_t tng_data,
                 int64_t *n)
{
    int64_t long_stride_length, medium_stride_length;
    long file_pos, orig_frame_set_file_pos;
    tng_trajectory_frame_set_t frame_set;
    struct tng_trajectory_frame_set orig_frame_set;
    tng_gen_block_t block;
    tng_function_status stat;
    int64_t cnt = 0;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(n, "TNG library: n must not be a NULL pointer");

    orig_frame_set = tng_data->current_trajectory_frame_set;

    frame_set = &tng_data->current_trajectory_frame_set;

    orig_frame_set_file_pos = tng_data->current_trajectory_frame_set_input_file_pos;
    file_pos = (long)tng_data->first_trajectory_frame_set_input_file_pos;

    if(file_pos < 0)
    {
        *n = tng_data->n_trajectory_frame_sets = cnt;
        return(TNG_SUCCESS);
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
        fprintf(stderr, "TNG library: Cannot read block header at pos %ld. %s: %d\n", file_pos,
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

    long_stride_length = tng_data->long_stride_length;
    medium_stride_length = tng_data->medium_stride_length;

    /* Take long steps forward until a long step forward would be too long or
     * the last frame set is found */
    file_pos = (long)frame_set->long_stride_next_frame_set_file_pos;
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
                fprintf(stderr, "TNG library: Cannot read block header at pos %ld. %s: %d\n",
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
        file_pos = (long)frame_set->long_stride_next_frame_set_file_pos;
    }

    /* Take medium steps forward until a medium step forward would be too long
     * or the last frame set is found */
    file_pos = (long)frame_set->medium_stride_next_frame_set_file_pos;
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
                fprintf(stderr, "TNG library: Cannot read block header at pos %ld. %s: %d\n",
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
        file_pos = (long)frame_set->medium_stride_next_frame_set_file_pos;
    }

    /* Take one step forward until the last frame set is found */
    file_pos = (long)frame_set->next_frame_set_file_pos;
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
                fprintf(stderr, "TNG library: Cannot read block header at pos %ld. %s: %d\n",
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
        file_pos = (long)frame_set->next_frame_set_file_pos;
    }

    tng_block_destroy(&block);

    *n = tng_data->n_trajectory_frame_sets = cnt;

    *frame_set = orig_frame_set;

    fseek(tng_data->input_file,
          (long)tng_data->first_trajectory_frame_set_input_file_pos,
          SEEK_SET);

    tng_data->current_trajectory_frame_set_input_file_pos = orig_frame_set_file_pos;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_current_frame_set_get
                (const tng_trajectory_t tng_data,
                 tng_trajectory_frame_set_t *frame_set_p)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    *frame_set_p = &tng_data->current_trajectory_frame_set;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_nr_find
                (tng_trajectory_t tng_data,
                 const int64_t nr)
{
    int64_t long_stride_length, medium_stride_length;
    int64_t file_pos, curr_nr = 0, n_frame_sets;
    tng_trajectory_frame_set_t frame_set;
    tng_gen_block_t block;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(nr >= 0, "The frame set number (nr) must be >= 0");

    frame_set = &tng_data->current_trajectory_frame_set;

    stat = tng_num_frame_sets_get(tng_data, &n_frame_sets);

    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

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
          (long)file_pos,
          SEEK_SET);
    tng_data->current_trajectory_frame_set_input_file_pos = (long)file_pos;
    /* Read block headers first to see what block is found. */
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n", file_pos,
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
            fseek(tng_data->input_file, (long)file_pos, SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
                  (long)file_pos,
                  SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
                  (long)file_pos,
                  SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
                  (long)file_pos,
                  SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
                  (long)file_pos,
                  SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
                  (long)file_pos,
                  SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
                  (long)file_pos,
                  SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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

tng_function_status DECLSPECDLLEXPORT tng_frame_set_of_frame_find
                (tng_trajectory_t tng_data,
                 const int64_t frame)
{
    int64_t first_frame, last_frame, n_frames_per_frame_set;
    int64_t long_stride_length, medium_stride_length;
    int64_t file_pos, temp_frame, n_frames;
    tng_trajectory_frame_set_t frame_set;
    tng_gen_block_t block;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame >= 0, "TNG library: frame must be >= 0.");

    frame_set = &tng_data->current_trajectory_frame_set;

    tng_block_init(&block);

    if(tng_data->current_trajectory_frame_set_input_file_pos < 0)
    {
        file_pos = tng_data->first_trajectory_frame_set_input_file_pos;
        fseek(tng_data->input_file,
                (long)file_pos,
                SEEK_SET);
        tng_data->current_trajectory_frame_set_input_file_pos = (long)file_pos;
        /* Read block headers first to see what block is found. */
        stat = tng_block_header_read(tng_data, block);
        if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
        {
            fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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

    first_frame = tng_max_i64(frame_set->first_frame, 0);
    last_frame = first_frame + frame_set->n_frames - 1;
    /* Is this the right frame set? */
    if(first_frame <= frame && frame <= last_frame)
    {
        tng_block_destroy(&block);
        return(TNG_SUCCESS);
    }

    n_frames_per_frame_set = tng_data->frame_set_n_frames;
    long_stride_length = tng_data->long_stride_length;
    medium_stride_length = tng_data->medium_stride_length;

    if(tng_first_frame_nr_of_next_frame_set_get(tng_data, &temp_frame) ==
       TNG_SUCCESS)
    {
        if(temp_frame - first_frame > n_frames_per_frame_set)
        {
            n_frames_per_frame_set = temp_frame - first_frame;
        }
    }

    tng_num_frames_get(tng_data, &n_frames);

    if(frame >= n_frames)
    {
        tng_block_destroy(&block);
        return(TNG_FAILURE);
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
        else if(frame - first_frame > (n_frames - 1) - frame)
        {
            file_pos = tng_data->last_trajectory_frame_set_input_file_pos;

            /* If the last frame set position is not set start from the current
             * frame set, since it will be closer than the first frame set. */
        }
        /* Start from current */
        else
        {
            file_pos = tng_data->current_trajectory_frame_set_input_file_pos;
        }

        if(file_pos > 0)
        {
            fseek(tng_data->input_file,
                  (long)file_pos,
                  SEEK_SET);
            tng_data->current_trajectory_frame_set_input_file_pos = (long)file_pos;
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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

    first_frame = tng_max_i64(frame_set->first_frame, 0);
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
            fseek(tng_data->input_file, (long)file_pos, SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
        first_frame = tng_max_i64(frame_set->first_frame, 0);
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
                  (long)file_pos,
                  SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
        first_frame = tng_max_i64(frame_set->first_frame, 0);
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
                  (long)file_pos,
                  SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
        first_frame = tng_max_i64(frame_set->first_frame, 0);
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
                  (long)file_pos,
                  SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
        first_frame = tng_max_i64(frame_set->first_frame, 0);
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
                  (long)file_pos,
                  SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
        first_frame = tng_max_i64(frame_set->first_frame, 0);
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
                  (long)file_pos,
                  SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
        first_frame = tng_max_i64(frame_set->first_frame, 0);
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
                  (long)file_pos,
                  SEEK_SET);
            /* Read block headers first to see what block is found. */
            stat = tng_block_header_read(tng_data, block);
            if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
            {
                fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
        first_frame = tng_max_i64(frame_set->first_frame, 0);
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

tng_function_status DECLSPECDLLEXPORT tng_frame_set_next_frame_set_file_pos_get
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *pos)
{
    (void)tng_data;

    TNG_ASSERT(frame_set, "TNG library: frame_set not initialised before accessing data.");
    TNG_ASSERT(pos, "TNG library: pos must not be a NULL pointer");

    *pos = frame_set->next_frame_set_file_pos;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_prev_frame_set_file_pos_get
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *pos)
{
    (void)tng_data;

    TNG_ASSERT(frame_set, "TNG library: frame_set not initialised before accessing data.");
    TNG_ASSERT(pos, "TNG library: pos must not be a NULL pointer");

    *pos = frame_set->prev_frame_set_file_pos;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_frame_range_get
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *first_frame,
                 int64_t *last_frame)
{
    (void)tng_data;

    TNG_ASSERT(first_frame, "TNG library: first_frame must not be a NULL pointer");
    TNG_ASSERT(last_frame, "TNG library: last_frame must not be a NULL pointer");
    TNG_ASSERT(frame_set, "TNG library: frame_set must not be a NULL pointer");

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
static TNG_INLINE tng_function_status tng_particle_mapping_get_real_particle
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
/*static TNG_INLINE tng_function_status tng_particle_mapping_get_local_particle
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
*/

static tng_function_status tng_file_headers_len_get
                (tng_trajectory_t tng_data,
                 int64_t *len)
{
    int64_t orig_pos;
    tng_gen_block_t block;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    *len = 0;

    orig_pos = ftell(tng_data->input_file);

    if(!tng_data->input_file_len)
    {
        fseek(tng_data->input_file, 0, SEEK_END);
        tng_data->input_file_len = ftell(tng_data->input_file);
    }
    fseek(tng_data->input_file, 0, SEEK_SET);

    tng_block_init(&block);
    /* Read through the headers of non-trajectory blocks (they come before the
     * trajectory blocks in the file) */
    while (*len < tng_data->input_file_len &&
           tng_block_header_read(tng_data, block) != TNG_CRITICAL &&
           block->id != -1 &&
           block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        *len += block->header_contents_size + block->block_contents_size;
        fseek(tng_data->input_file, block->block_contents_size, SEEK_CUR);
    }

    fseek(tng_data->input_file, orig_pos, SEEK_SET);

    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_file_headers_read
                (tng_trajectory_t tng_data,
                 const char hash_mode)
{
    int64_t prev_pos = 0;
    tng_gen_block_t block;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    tng_data->n_trajectory_frame_sets = 0;

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    if(!tng_data->input_file_len)
    {
        fseek(tng_data->input_file, 0, SEEK_END);
        tng_data->input_file_len = ftell(tng_data->input_file);
    }
    fseek(tng_data->input_file, 0, SEEK_SET);

    tng_block_init(&block);
    /* Non trajectory blocks (they come before the trajectory
     * blocks in the file) */
    while (prev_pos < tng_data->input_file_len &&
           tng_block_header_read(tng_data, block) != TNG_CRITICAL &&
           block->id != -1 &&
           block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        tng_block_read_next(tng_data, block, hash_mode);
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

tng_function_status DECLSPECDLLEXPORT tng_file_headers_write
                (tng_trajectory_t tng_data,
                 const char hash_mode)
{
    int i;
    int64_t len, orig_len, tot_len = 0, data_start_pos;
    tng_function_status stat;
    tng_gen_block_t block;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    if(tng_data->n_trajectory_frame_sets > 0)
    {
        stat = tng_file_headers_len_get(tng_data, &orig_len);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }

        tng_block_init(&block);
        block->name = malloc(TNG_MAX_STR_LEN);
        if(!block->name)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
                    TNG_MAX_STR_LEN, __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
        strcpy(block->name, "GENERAL INFO");
        tng_block_header_len_calculate(tng_data, block, &len);
        tot_len += len;
        tng_general_info_block_len_calculate(tng_data, &len);
        tot_len += len;
        strcpy(block->name, "MOLECULES");
        tng_block_header_len_calculate(tng_data, block, &len);
        tot_len += len;
        tng_molecules_block_len_calculate(tng_data, &len);
        tot_len += len;

        for(i = 0; i < tng_data->n_data_blocks; i++)
        {
            strcpy(block->name, tng_data->non_tr_data[i].block_name);
            tng_block_header_len_calculate(tng_data, block, &len);
            tot_len += len;
            tng_data_block_len_calculate(tng_data,
                                        (tng_particle_data_t)&tng_data->non_tr_data[i],
                                        TNG_FALSE, 1, 1, 1, 0,
                                        1, 0, &data_start_pos,
                                        &len);
            tot_len += len;
        }
        for(i = 0; i < tng_data->n_particle_data_blocks; i++)
        {
            strcpy(block->name, tng_data->non_tr_particle_data[i].block_name);
            tng_block_header_len_calculate(tng_data, block, &len);
            tot_len += len;
            tng_data_block_len_calculate(tng_data,
                                        &tng_data->non_tr_particle_data[i],
                                        TNG_TRUE, 1, 1, 1, 0,
                                        tng_data->n_particles, TNG_PARTICLE_DEPENDENT,
                                        &data_start_pos,
                                        &len);
            tot_len += len;
        }
        tng_block_destroy(&block);

        if(tot_len > orig_len)
        {
            tng_migrate_data_in_file(tng_data, orig_len+1, tot_len - orig_len);
        }

        tng_data->current_trajectory_frame_set_output_file_pos = -1;
    }

    /* TODO: If there is already frame set data written to this file (e.g. when
     * appending to an already existing file we might need to move frame sets to
     * the end of the file. */

    if(tng_general_info_block_write(tng_data, hash_mode)
       != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Error writing general info block of file %s. %s: %d\n",
                tng_data->input_file_path, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    if(tng_molecules_block_write(tng_data, hash_mode)
        != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Error writing atom names block of file %s. %s: %d\n",
                tng_data->input_file_path, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    /* FIXME: Currently writing non-trajectory data blocks here.
     * Should perhaps be moved. */
    tng_block_init(&block);
    for(i = 0; i < tng_data->n_data_blocks; i++)
    {
        block->id = tng_data->non_tr_data[i].block_id;
        tng_data_block_write(tng_data, block,
                             i, hash_mode);
    }

    for(i = 0; i < tng_data->n_particle_data_blocks; i++)
    {
        block->id = tng_data->non_tr_particle_data[i].block_id;
        tng_particle_data_block_write(tng_data, block,
                                      i, 0, hash_mode);
    }

    tng_block_destroy(&block);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_block_read_next(tng_trajectory_t tng_data,
                                        tng_gen_block_t block,
                                        const char hash_mode)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(block, "TNG library: block must be initialised and must not be a NULL pointer.");

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
            fseek(tng_data->input_file, (long)block->block_contents_size, SEEK_CUR);
            return(TNG_FAILURE);
        }
    }
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_read
                (tng_trajectory_t tng_data,
                 const char hash_mode)
{
    long file_pos;
    tng_gen_block_t block;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    file_pos = ftell(tng_data->input_file);

    tng_block_init(&block);

    if(!tng_data->input_file_len)
    {
        fseek(tng_data->input_file, 0, SEEK_END);
        tng_data->input_file_len = ftell(tng_data->input_file);
        fseek(tng_data->input_file, file_pos, SEEK_SET);
    }

    /* Read block headers first to see what block is found. */
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET ||
       block->id == -1)
    {
        fprintf(stderr, "TNG library: Cannot read block header at pos %ld. %s: %d\n",
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
              block->id != TNG_TRAJECTORY_FRAME_SET &&
              block->id != -1)
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
            fprintf(stderr, "TNG library: Cannot read block header at pos %ld. %s: %d\n",
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


tng_function_status DECLSPECDLLEXPORT tng_frame_set_read_current_only_data_from_block_id
                (tng_trajectory_t tng_data,
                 const char hash_mode,
                 const int64_t block_id)
{
    long file_pos;
    tng_gen_block_t block;
    tng_function_status stat;
    int found_flag = 1;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    file_pos = (long)tng_data->current_trajectory_frame_set_input_file_pos;

    if(file_pos < 0)
    {
        /* No current frame set. This means that the first frame set must be
         * read */
        found_flag = 0;
        file_pos = (long)tng_data->first_trajectory_frame_set_input_file_pos;
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
        fprintf(stderr, "TNG library: Cannot read block header at pos %ld. %s: %d\n",
               file_pos, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    /* If the current frame set had already been read skip its block contents */
    if(found_flag)
    {
        fseek(tng_data->input_file, (long)block->block_contents_size, SEEK_CUR);
    }
    /* Otherwiese read the frame set block */
    else
    {
        stat = tng_block_read_next(tng_data, block,
                                   hash_mode);
        if(stat != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot read frame set block. %s: %d\n", __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(stat);
        }
    }
    file_pos = ftell(tng_data->input_file);

    found_flag = 0;

    /* Read only blocks of the requested ID
        * until next frame set block */
    stat = tng_block_header_read(tng_data, block);
    while(file_pos < tng_data->input_file_len &&
            stat != TNG_CRITICAL &&
            block->id != TNG_TRAJECTORY_FRAME_SET &&
            block->id != -1)
    {
        if(block->id == block_id)
        {
            stat = tng_block_read_next(tng_data, block,
                                       hash_mode);
            if(stat != TNG_CRITICAL)
            {
                file_pos = ftell(tng_data->input_file);
                found_flag = 1;
                if(file_pos < tng_data->input_file_len)
                {
                    stat = tng_block_header_read(tng_data, block);
                }
            }
        }
        else
        {
            file_pos += (long)(block->block_contents_size + block->header_contents_size);
            fseek(tng_data->input_file, (long)block->block_contents_size, SEEK_CUR);
            if(file_pos < tng_data->input_file_len)
            {
                stat = tng_block_header_read(tng_data, block);
            }
        }
    }
    if(stat == TNG_CRITICAL)
    {
        fprintf(stderr, "TNG library: Cannot read block header at pos %ld. %s: %d\n",
                file_pos, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(stat);
    }

    if(block->id == TNG_TRAJECTORY_FRAME_SET)
    {
        fseek(tng_data->input_file, file_pos, SEEK_SET);
    }

    tng_block_destroy(&block);

    if(found_flag)
    {
        return(TNG_SUCCESS);
    }
    else
    {
        return(TNG_FAILURE);
    }
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_read_next
                (tng_trajectory_t tng_data,
                 const char hash_mode)
{
    long file_pos;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    file_pos = (long)tng_data->current_trajectory_frame_set.next_frame_set_file_pos;

    if(file_pos < 0 && tng_data->current_trajectory_frame_set_input_file_pos <= 0)
    {
        file_pos = (long)tng_data->first_trajectory_frame_set_input_file_pos;
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

    return(tng_frame_set_read(tng_data, hash_mode));
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_read_next_only_data_from_block_id
                (tng_trajectory_t tng_data,
                 const char hash_mode,
                 const int64_t block_id)
{
    long file_pos;
    tng_gen_block_t block;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    if(tng_input_file_init(tng_data) != TNG_SUCCESS)
    {
        return(TNG_CRITICAL);
    }

    file_pos = (long)tng_data->current_trajectory_frame_set.next_frame_set_file_pos;

    if(file_pos < 0 && tng_data->current_trajectory_frame_set_input_file_pos <= 0)
    {
        file_pos = (long)tng_data->first_trajectory_frame_set_input_file_pos;
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
        fprintf(stderr, "TNG library: Cannot read block header at pos %ld. %s: %d\n",
               file_pos, __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    tng_data->current_trajectory_frame_set_input_file_pos = file_pos;

    if(tng_block_read_next(tng_data, block,
                           hash_mode) == TNG_SUCCESS)
    {
        stat = tng_frame_set_read_current_only_data_from_block_id(tng_data, hash_mode, block_id);
    }

    tng_block_destroy(&block);

    return(stat);
}

tng_function_status tng_frame_set_write(tng_trajectory_t tng_data,
                                        const char hash_mode)
{
    int i, j;
    tng_gen_block_t block;
    tng_trajectory_frame_set_t frame_set;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    frame_set = &tng_data->current_trajectory_frame_set;

    if(frame_set->n_written_frames == frame_set->n_frames)
    {
        return(TNG_SUCCESS);
    }

    tng_data->current_trajectory_frame_set_output_file_pos =
    ftell(tng_data->output_file);
    tng_data->last_trajectory_frame_set_output_file_pos =
    tng_data->current_trajectory_frame_set_output_file_pos;

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

    frame_set->n_unwritten_frames = 0;

    fflush(tng_data->output_file);

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_premature_write
                (tng_trajectory_t tng_data,
                 const char hash_mode)
{
    tng_trajectory_frame_set_t frame_set;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");

    frame_set = &tng_data->current_trajectory_frame_set;

    if(frame_set->n_unwritten_frames == 0)
    {
        return(TNG_SUCCESS);
    }
    frame_set->n_frames = frame_set->n_unwritten_frames;

    return(tng_frame_set_write(tng_data, hash_mode));
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_new
                (tng_trajectory_t tng_data,
                 const int64_t first_frame,
                 const int64_t n_frames)
{
    tng_gen_block_t block;
    tng_trajectory_frame_set_t frame_set;
    FILE *temp = tng_data->input_file;
    int64_t curr_pos;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(first_frame >= 0, "TNG library: first_frame must be >= 0.");
    TNG_ASSERT(n_frames >= 0, "TNG library: n_frames must be >= 0.");

    frame_set = &tng_data->current_trajectory_frame_set;

    curr_pos = ftell(tng_data->output_file);

    if(curr_pos <= 10)
    {
        tng_file_headers_write(tng_data, TNG_USE_HASH);
    }

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
                  (long)frame_set->medium_stride_prev_frame_set_file_pos,
                  SEEK_SET);

            if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot read frame set header. %s: %d\n",
                    __FILE__, __LINE__);
                tng_data->input_file = temp;
                tng_block_destroy(&block);
                return(TNG_CRITICAL);
            }

            /* Read the next frame set from the previous frame set and one
             * medium stride step back */
            fseek(tng_data->output_file, (long)block->block_contents_size - (6 *
            sizeof(int64_t) + 2 * sizeof(double)), SEEK_CUR);
            if(fread(&frame_set->medium_stride_prev_frame_set_file_pos,
               sizeof(frame_set->medium_stride_prev_frame_set_file_pos),
               1, tng_data->output_file) == 0)
            {
                fprintf(stderr, "TNG library: Cannot read block. %s: %d\n", __FILE__, __LINE__);
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
                    fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                          (long)frame_set->long_stride_prev_frame_set_file_pos,
                          SEEK_SET);

                    if(tng_block_header_read(tng_data, block) != TNG_SUCCESS)
                    {
                        fprintf(stderr, "TNG library: Cannot read frame set header. %s: %d\n",
                            __FILE__, __LINE__);
                        tng_data->input_file = temp;
                        tng_block_destroy(&block);
                        return(TNG_CRITICAL);
                    }

                    /* Read the next frame set from the previous frame set and one
                    * long stride step back */
                    fseek(tng_data->output_file, (long)block->block_contents_size - (6 *
                          sizeof(int64_t) + 2 * sizeof(double)), SEEK_CUR);

                    tng_block_destroy(&block);

                    if(fread(&frame_set->long_stride_prev_frame_set_file_pos,
                    sizeof(frame_set->long_stride_prev_frame_set_file_pos),
                    1, tng_data->output_file) == 0)
                    {
                        fprintf(stderr, "TNG library: Cannot read block. %s: %d\n", __FILE__, __LINE__);
                        tng_data->input_file = temp;
                        return(TNG_CRITICAL);
                    }

                    if(tng_data->input_endianness_swap_func_64)
                    {
                        if(tng_data->input_endianness_swap_func_64(tng_data,
                           &frame_set->long_stride_prev_frame_set_file_pos)
                            != TNG_SUCCESS)
                        {
                            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                                    __FILE__, __LINE__);
                        }
                    }

                }
            }

            tng_data->input_file = temp;
            fseek(tng_data->output_file, (long)curr_pos, SEEK_SET);
        }
    }

    frame_set->first_frame = first_frame;
    frame_set->n_frames = n_frames;
    frame_set->n_written_frames = 0;
    frame_set->n_unwritten_frames = 0;
    frame_set->first_frame_time = -1;

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

tng_function_status DECLSPECDLLEXPORT tng_frame_set_with_time_new
                (tng_trajectory_t tng_data,
                 const int64_t first_frame,
                 const int64_t n_frames,
                 const double first_frame_time)
{
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(first_frame >= 0, "TNG library: first_frame must be >= 0.");
    TNG_ASSERT(n_frames >= 0, "TNG library: n_frames must be >= 0.");
    TNG_ASSERT(first_frame_time >= 0, "TNG library: first_frame_time must be >= 0.");


    stat = tng_frame_set_new(tng_data, first_frame, n_frames);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }
    stat = tng_frame_set_first_frame_time_set(tng_data, first_frame_time);

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_frame_set_first_frame_time_set
                (tng_trajectory_t tng_data,
                 const double first_frame_time)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(first_frame_time >= 0, "TNG library: first_frame_time must be >= 0.");

    tng_data->current_trajectory_frame_set.first_frame_time = first_frame_time;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_first_frame_nr_of_next_frame_set_get
                (const tng_trajectory_t tng_data,
                 int64_t *frame)
{
    long file_pos, next_frame_set_file_pos;
    tng_gen_block_t block;
    tng_function_status stat;

    tng_trajectory_frame_set_t frame_set;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(tng_data->input_file, "TNG library: An input file must be open to find the next frame set");
    TNG_ASSERT(frame, "TNG library: frame must not be a NULL pointer");

    file_pos = ftell(tng_data->input_file);

    if(tng_data->current_trajectory_frame_set_input_file_pos <= 0)
    {
        next_frame_set_file_pos = (long)tng_data->first_trajectory_frame_set_input_file_pos;
    }
    else
    {
        frame_set = &tng_data->current_trajectory_frame_set;
        next_frame_set_file_pos = (long)frame_set->next_frame_set_file_pos;
    }

    if(next_frame_set_file_pos <= 0)
    {
        return(TNG_FAILURE);
    }

    fseek(tng_data->input_file, (long)next_frame_set_file_pos, SEEK_SET);
    /* Read block headers first to see that a frame set block is found. */
    tng_block_init(&block);
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL || block->id != TNG_TRAJECTORY_FRAME_SET)
    {
        fprintf(stderr, "TNG library: Cannot read block header at pos %ld. %s: %d\n",
               file_pos, __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }
/*    if(tng_data->current_trajectory_frame_set_input_file_pos <= 0)
    {
        tng_block_read_next(tng_data, block, TNG_USE_HASH);
    }*/
    tng_block_destroy(&block);

    if(fread(frame, sizeof(int64_t), 1, tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Cannot read first frame of next frame set. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }
    fseek(tng_data->input_file, file_pos, SEEK_SET);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_data_block_add
                (tng_trajectory_t tng_data,
                 const int64_t id,
                 const char *block_name,
                 const char datatype,
                 const char block_type_flag,
                 int64_t n_frames,
                 const int64_t n_values_per_frame,
                 int64_t stride_length,
                 const int64_t codec_id,
                 void *new_data)
{
    int i, j, size, len;
    tng_trajectory_frame_set_t frame_set;
    tng_non_particle_data_t data;
    char **first_dim_values;
    char *new_data_c=new_data;
    int64_t n_frames_div;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(block_name, "TNG library: block_name must not be a NULL pointer.");
    TNG_ASSERT(n_values_per_frame > 0, "TNG library: n_values_per_frame must be a positive integer.");

    frame_set = &tng_data->current_trajectory_frame_set;

    if(stride_length <= 0)
    {
        stride_length = 1;
    }

    /* If the block does not exist, create it */
    if(tng_data_find(tng_data, id, &data) != TNG_SUCCESS)
    {
        if(tng_data_block_create(tng_data, block_type_flag) !=
            TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot create data block. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
                   (int)strlen(block_name)+1, __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        strncpy(data->block_name, block_name, strlen(block_name) + 1);

        data->values = 0;
        /* FIXME: Memory leak from strings. */
        data->strings = 0;
        data->last_retrieved_frame = -1;
    }

    data->datatype = datatype;
    data->stride_length = tng_max_i64(stride_length, 1);
    data->n_values_per_frame = n_values_per_frame;
    data->n_frames = n_frames;
    data->codec_id = codec_id;
    data->compression_multiplier = 1.0;
    /* FIXME: This can cause problems. */
    data->first_frame_with_data = frame_set->first_frame;

    switch(datatype)
    {
    case TNG_FLOAT_DATA:
        size = sizeof(float);
        break;
    case TNG_INT_DATA:
        size = sizeof(int64_t);
        break;
    case TNG_DOUBLE_DATA:
    default:
        size = sizeof(double);
        break;
    }

    if(new_data_c)
    {
        /* Allocate memory */
        if(tng_allocate_data_mem(tng_data, data, n_frames, stride_length,
                                 n_values_per_frame) !=
        TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot allocate data memory. %s: %d\n",
                __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }

        if(n_frames > frame_set->n_unwritten_frames)
        {
            frame_set->n_unwritten_frames = n_frames;
        }

        n_frames_div = (n_frames % stride_length) ?
                     n_frames / stride_length + 1:
                     n_frames / stride_length;

        if(datatype == TNG_CHAR_DATA)
        {
            for(i = 0; i < n_frames_div; i++)
            {
                first_dim_values = data->strings[i];
                for(j = 0; j < n_values_per_frame; j++)
                {
                    len = tng_min_i((int)strlen(new_data_c) + 1,
                                TNG_MAX_STR_LEN);
                    if(first_dim_values[j])
                    {
                        free(first_dim_values[j]);
                    }
                    first_dim_values[j] = malloc(len);
                    if(!first_dim_values[j])
                    {
                        fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
                            len, __FILE__, __LINE__);
                        return(TNG_CRITICAL);
                    }
                    strncpy(first_dim_values[j],
                            new_data_c, len);
                    new_data_c += len;
                }
            }
        }
        else
        {
            memcpy(data->values, new_data, size * n_frames_div *
                   n_values_per_frame);
        }
    }

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_particle_data_block_add
                (tng_trajectory_t tng_data,
                 const int64_t id,
                 const char *block_name,
                 const char datatype,
                 const char block_type_flag,
                 int64_t n_frames,
                 const int64_t n_values_per_frame,
                 int64_t stride_length,
                 const int64_t num_first_particle,
                 const int64_t n_particles,
                 const int64_t codec_id,
                 void *new_data)
{
    int i, size, len;
    int64_t j, k;
    int64_t tot_n_particles, n_frames_div;
    char ***first_dim_values, **second_dim_values;
    tng_trajectory_frame_set_t frame_set;
    tng_particle_data_t data;
    char *new_data_c=new_data;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(block_name, "TNG library: block_name mustnot be a NULL pointer.");
    TNG_ASSERT(n_values_per_frame > 0, "TNG library: n_values_per_frame must be a positive integer.");
    TNG_ASSERT(num_first_particle >= 0, "TNG library: num_first_particle must be >= 0.");
    TNG_ASSERT(n_particles >= 0, "TNG library: n_particles must be >= 0.");


    frame_set = &tng_data->current_trajectory_frame_set;

    if(stride_length <= 0)
    {
        stride_length = 1;
    }

    /* If the block does not exist, create it */
    if(tng_particle_data_find(tng_data, id, &data) != TNG_SUCCESS)
    {
        if(tng_particle_data_block_create(tng_data, block_type_flag) !=
            TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot create particle data block. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
                   (int)strlen(block_name)+1, __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }
        strncpy(data->block_name, block_name, strlen(block_name) + 1);

        data->datatype = datatype;

        data->values = 0;
        /* FIXME: Memory leak from strings. */
        data->strings = 0;
        data->last_retrieved_frame = -1;
    }

    data->stride_length = tng_max_i64(stride_length, 1);
    data->n_values_per_frame = n_values_per_frame;
    data->n_frames = n_frames;
    data->codec_id = codec_id;
    data->compression_multiplier = 1.0;
    /* FIXME: This can cause problems. */
    data->first_frame_with_data = frame_set->first_frame;

    if(block_type_flag == TNG_TRAJECTORY_BLOCK && tng_data->var_num_atoms_flag)
    {
        tot_n_particles = frame_set->n_particles;
    }
    else
    {
        tot_n_particles = tng_data->n_particles;
    }

    /* If data values are supplied add that data to the data block. */
    if(new_data_c)
    {
        /* Allocate memory */
        if(tng_allocate_particle_data_mem(tng_data, data, n_frames,
                                          stride_length, tot_n_particles,
                                          n_values_per_frame) !=
        TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot allocate particle data memory. %s: %d\n",
                __FILE__, __LINE__);
            return(TNG_CRITICAL);
        }

        if(n_frames > frame_set->n_unwritten_frames)
        {
            frame_set->n_unwritten_frames = n_frames;
        }

        n_frames_div = (n_frames % stride_length) ?
                     n_frames / stride_length + 1:
                     n_frames / stride_length;

        if(datatype == TNG_CHAR_DATA)
        {
            for(i = 0; i < n_frames_div; i++)
            {
                first_dim_values = data->strings[i];
                for(j = num_first_particle; j < num_first_particle + n_particles;
                    j++)
                {
                    second_dim_values = first_dim_values[j];
                    for(k = 0; k < n_values_per_frame; k++)
                    {
                        len = tng_min_i((int)strlen(new_data_c) + 1,
                                TNG_MAX_STR_LEN);
                        if(second_dim_values[k])
                        {
                            free(second_dim_values[k]);
                        }
                        second_dim_values[k] = malloc(len);
                        if(!second_dim_values[k])
                        {
                            fprintf(stderr, "TNG library: Cannot allocate memory (%d bytes). %s: %d\n",
                                len, __FILE__, __LINE__);
                            return(TNG_CRITICAL);
                        }
                        strncpy(second_dim_values[k],
                                new_data_c, len);
                        new_data_c += len;
                    }
                }
            }
        }
        else
        {
            switch(datatype)
            {
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

            memcpy(data->values, new_data, size * n_frames_div *
                   n_particles * n_values_per_frame);
        }
    }

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_data_block_name_get
                (tng_trajectory_t tng_data,
                 int64_t block_id,
                 char *name,
                 int max_len)
{
    int64_t i;
    tng_trajectory_frame_set_t frame_set;
    tng_function_status stat;
    tng_particle_data_t p_data;
    tng_non_particle_data_t np_data;
    int block_type = -1;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer.");

    for(i = 0; i < tng_data->n_particle_data_blocks; i++)
    {
        p_data = &tng_data->non_tr_particle_data[i];
        if(p_data->block_id == block_id)
        {
            strncpy(name, p_data->block_name, max_len);
            name[max_len - 1] = '\0';
            return(TNG_SUCCESS);
        }
    }
    for(i = 0; i < tng_data->n_data_blocks; i++)
    {
        np_data = &tng_data->non_tr_data[i];
        if(np_data->block_id == block_id)
        {
            strncpy(name, np_data->block_name, max_len);
            name[max_len - 1] = '\0';
            return(TNG_SUCCESS);
        }
    }

    frame_set = &tng_data->current_trajectory_frame_set;

    stat = tng_particle_data_find(tng_data, block_id, &p_data);
    if(stat == TNG_SUCCESS)
    {
        block_type = TNG_PARTICLE_BLOCK_DATA;
    }
    else
    {
        stat = tng_data_find(tng_data, block_id, &np_data);
        if(stat == TNG_SUCCESS)
        {
            block_type = TNG_NON_PARTICLE_BLOCK_DATA;
        }
        else
        {
            stat = tng_frame_set_read_current_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
            if(stat != TNG_SUCCESS)
            {
                return(stat);
            }
            stat = tng_particle_data_find(tng_data, block_id, &p_data);
            if(stat == TNG_SUCCESS)
            {
                block_type = TNG_PARTICLE_BLOCK_DATA;
            }
            else
            {
                stat = tng_data_find(tng_data, block_id, &np_data);
                if(stat == TNG_SUCCESS)
                {
                    block_type = TNG_NON_PARTICLE_BLOCK_DATA;
                }
            }
        }
    }
    if(block_type == TNG_PARTICLE_BLOCK_DATA)
    {
        for(i = 0; i < frame_set->n_particle_data_blocks; i++)
        {
            p_data = &frame_set->tr_particle_data[i];
            if(p_data->block_id == block_id)
            {
                strncpy(name, p_data->block_name, max_len);
                name[max_len - 1] = '\0';
                return(TNG_SUCCESS);
            }
        }
    }
    else if(block_type == TNG_NON_PARTICLE_BLOCK_DATA)
    {
        for(i = 0; i < frame_set->n_data_blocks; i++)
        {
            np_data = &frame_set->tr_data[i];
            if(np_data->block_id == block_id)
            {
                strncpy(name, np_data->block_name, max_len);
                name[max_len - 1] = '\0';
                return(TNG_SUCCESS);
            }
        }
    }

    return(TNG_FAILURE);
}

tng_function_status DECLSPECDLLEXPORT tng_data_block_dependency_get
                (const tng_trajectory_t tng_data,
                 int64_t block_id,
                 int *block_dependency)
{
    int64_t i;
    tng_function_status stat;
    tng_particle_data_t p_data;
    tng_non_particle_data_t np_data;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(block_dependency, "TNG library: block_dependency must not be a NULL pointer.");

    for(i = 0; i < tng_data->n_particle_data_blocks; i++)
    {
        p_data = &tng_data->non_tr_particle_data[i];
        if(p_data->block_id == block_id)
        {
            *block_dependency = TNG_PARTICLE_DEPENDENT;
            return(TNG_SUCCESS);
        }
    }
    for(i = 0; i < tng_data->n_data_blocks; i++)
    {
        np_data = &tng_data->non_tr_data[i];
        if(np_data->block_id == block_id)
        {
            *block_dependency = 0;
            return(TNG_SUCCESS);
        }
    }

    stat = tng_particle_data_find(tng_data, block_id, &p_data);
    if(stat == TNG_SUCCESS)
    {
        *block_dependency = TNG_PARTICLE_DEPENDENT + TNG_FRAME_DEPENDENT;
        return(TNG_SUCCESS);
    }
    else
    {
        stat = tng_data_find(tng_data, block_id, &np_data);
        if(stat == TNG_SUCCESS)
        {
            *block_dependency = TNG_FRAME_DEPENDENT;
            return(TNG_SUCCESS);
        }
        else
        {
            stat = tng_frame_set_read_current_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
            if(stat != TNG_SUCCESS)
            {
                return(stat);
            }
            stat = tng_particle_data_find(tng_data, block_id, &p_data);
            if(stat == TNG_SUCCESS)
            {
                *block_dependency = TNG_PARTICLE_DEPENDENT + TNG_FRAME_DEPENDENT;
                return(TNG_SUCCESS);
            }
            else
            {
                stat = tng_data_find(tng_data, block_id, &np_data);
                if(stat == TNG_SUCCESS)
                {
                    *block_dependency = TNG_FRAME_DEPENDENT;
                    return(TNG_SUCCESS);
                }
            }
        }
    }

    return(TNG_FAILURE);
}

tng_function_status DECLSPECDLLEXPORT tng_data_block_num_values_per_frame_get
                (const tng_trajectory_t tng_data,
                 int64_t block_id,
                 int64_t *n_values_per_frame)
{
    int64_t i;
    tng_function_status stat;
    tng_particle_data_t p_data;
    tng_non_particle_data_t np_data;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(n_values_per_frame, "TNG library: n_values_per_frame must not be a NULL pointer.");

    for(i = 0; i < tng_data->n_particle_data_blocks; i++)
    {
        p_data = &tng_data->non_tr_particle_data[i];
        if(p_data->block_id == block_id)
        {
            *n_values_per_frame = p_data->n_values_per_frame;
            return(TNG_SUCCESS);
        }
    }
    for(i = 0; i < tng_data->n_data_blocks; i++)
    {
        np_data = &tng_data->non_tr_data[i];
        if(np_data->block_id == block_id)
        {
            *n_values_per_frame = np_data->n_values_per_frame;
            return(TNG_SUCCESS);
        }
    }

    stat = tng_particle_data_find(tng_data, block_id, &p_data);
    if(stat == TNG_SUCCESS)
    {
        *n_values_per_frame = p_data->n_values_per_frame;
        return(TNG_SUCCESS);
    }
    else
    {
        stat = tng_data_find(tng_data, block_id, &np_data);
        if(stat == TNG_SUCCESS)
        {
            *n_values_per_frame = np_data->n_values_per_frame;
            return(TNG_SUCCESS);
        }
        else
        {
            stat = tng_frame_set_read_current_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
            if(stat != TNG_SUCCESS)
            {
                return(stat);
            }
            stat = tng_particle_data_find(tng_data, block_id, &p_data);
            if(stat == TNG_SUCCESS)
            {
                *n_values_per_frame = p_data->n_values_per_frame;
                return(TNG_SUCCESS);
            }
            else
            {
                stat = tng_data_find(tng_data, block_id, &np_data);
                if(stat == TNG_SUCCESS)
                {
                    *n_values_per_frame = np_data->n_values_per_frame;
                    return(TNG_SUCCESS);
                }
            }
        }
    }

    return(TNG_FAILURE);
}

tng_function_status DECLSPECDLLEXPORT tng_frame_data_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const int64_t block_id,
                 const void *values,
                 const char hash_mode)
{
    int64_t header_pos, file_pos;
    int64_t output_file_len, n_values_per_frame, size, contents_size;
    int64_t header_size, temp_first, temp_last;
    int64_t i, last_frame;
    long temp_current;
    tng_gen_block_t block;
    tng_trajectory_frame_set_t frame_set;
    FILE *temp = tng_data->input_file;
    struct tng_non_particle_data data;
    tng_function_status stat;
    char dependency, sparse_data, datatype;
    void *copy;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(values, "TNG library: values must not be a NULL pointer.");

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot initialise destination file. %s: %d\n",
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
            last_frame < frame_nr)
/*           (last_frame < frame_nr &&
            tng_data->current_trajectory_frame_set.first_frame +
            tng_data->frame_set_n_frames >= frame_nr))*/
        {
            if(last_frame + tng_data->frame_set_n_frames < frame_nr)
            {
                last_frame = frame_nr - 1;
            }
            tng_frame_set_new(tng_data,
                              last_frame+1,
                              tng_data->frame_set_n_frames);
            file_pos = ftell(tng_data->output_file);
            fseek(tng_data->output_file, 0, SEEK_END);
            output_file_len = ftell(tng_data->output_file);
            fseek(tng_data->output_file, (long)file_pos, SEEK_SET);

            /* Read mapping blocks from the last frame set */
            tng_block_init(&block);

            stat = tng_block_header_read(tng_data, block);
            while(file_pos < output_file_len &&
                  stat != TNG_CRITICAL &&
                  block->id != TNG_TRAJECTORY_FRAME_SET &&
                  block->id != -1)
            {
                if(block->id == TNG_PARTICLE_MAPPING)
                {
                    tng_trajectory_mapping_block_read(tng_data, block,
                                                      hash_mode);
                }
                else
                {
                    fseek(tng_data->output_file, (long)block->block_contents_size,
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
                fprintf(stderr, "TNG library: Error writing frame set. %s: %d\n", __FILE__, __LINE__);
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
    fseek(tng_data->output_file, (long)file_pos, SEEK_SET);

    /* Read past the frame set block first */
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL)
    {
        fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
               file_pos, __FILE__, __LINE__);
        tng_block_destroy(&block);
        tng_data->input_file = temp;

        tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
        tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
        tng_data->current_trajectory_frame_set_input_file_pos = temp_current;
        return(stat);
    }
    fseek(tng_data->output_file, (long)block->block_contents_size,
            SEEK_CUR);

    /* Read all block headers until next frame set block or
     * until the wanted block id is found */
    stat = tng_block_header_read(tng_data, block);
    while(file_pos < output_file_len &&
            stat != TNG_CRITICAL &&
            block->id != block_id &&
            block->id != TNG_TRAJECTORY_FRAME_SET &&
            block->id != -1)
    {
        fseek(tng_data->output_file, (long)block->block_contents_size, SEEK_CUR);
        file_pos = ftell(tng_data->output_file);
        if(file_pos < output_file_len)
        {
            stat = tng_block_header_read(tng_data, block);
        }
    }
    if(stat == TNG_CRITICAL)
    {
        fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
        fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(fread(&dependency, sizeof(dependency), 1, tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
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
        fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    if(fread(&data.n_values_per_frame, sizeof(data.n_values_per_frame), 1,
             tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
            &data.n_values_per_frame)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(fread(&data.codec_id, sizeof(data.codec_id), 1,
             tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
            &data.codec_id)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(data.codec_id != TNG_UNCOMPRESSED)
    {
        if(fread(&data.compression_multiplier,
                 sizeof(data.compression_multiplier), 1, tng_data->input_file)
            == 0)
        {
            fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                (int64_t *)&data.compression_multiplier)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                &data.first_frame_with_data)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }

        if(fread(&data.stride_length, sizeof(data.stride_length),
                 1, tng_data->input_file) == 0)
        {
            fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                &data.stride_length)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot calculate writing locations. %s: %d.\n", __FILE__,
                   __LINE__);
            tng_block_destroy(&block);
            return(TNG_FAILURE);
    }

    n_values_per_frame = data.n_values_per_frame;

    file_pos = (frame_nr - tng_max_i64(frame_set->first_frame,
                               data.first_frame_with_data)) /
                data.stride_length;
    file_pos *= size * n_values_per_frame;

    if(file_pos > contents_size)
    {
        fprintf(stderr, "TNG library: Attempting to write outside the block. %s: %d\n", __FILE__,
               __LINE__);
        tng_block_destroy(&block);
        return(TNG_FAILURE);
    }

    fseek(tng_data->output_file, (long)file_pos, SEEK_CUR);

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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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

tng_function_status DECLSPECDLLEXPORT tng_frame_particle_data_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const int64_t block_id,
                 const int64_t val_first_particle,
                 const int64_t val_n_particles,
                 const void *values,
                 const char hash_mode)
{
    int64_t header_pos, file_pos, tot_n_particles;
    int64_t output_file_len, n_values_per_frame, size, contents_size;
    int64_t header_size, temp_first, temp_last;
    int64_t mapping_block_end_pos, num_first_particle, block_n_particles;
    int64_t i, last_frame;
    long temp_current;
    tng_gen_block_t block;
    tng_trajectory_frame_set_t frame_set;
    FILE *temp = tng_data->input_file;
    struct tng_particle_data data;
    tng_function_status stat;
    tng_particle_mapping_t mapping;
    char dependency, sparse_data, datatype;
    void *copy;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(values, "TNG library: values must not be a NULL pointer.");
    TNG_ASSERT(val_first_particle >= 0, "TNG library: val_first_particle must be >= 0.");
    TNG_ASSERT(val_n_particles >= 0, "TNG library: val_n_particles must be >= 0.");

    if(tng_output_file_init(tng_data) != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot initialise destination file. %s: %d\n",
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
/*         fprintf(stderr, "TNG library: Frame %"PRId64" not found. Last frame: %"PRId64"\n", frame_nr,
                  last_frame); */
        /* If the wanted frame would be in the frame set after the last
         * frame set create a new frame set. */
        if(stat == TNG_FAILURE &&
           (last_frame < frame_nr &&
            last_frame + tng_data->frame_set_n_frames >= frame_nr))
        {
            if(last_frame + tng_data->frame_set_n_frames < frame_nr)
            {
                last_frame = frame_nr - 1;
            }
            tng_frame_set_new(tng_data,
                              last_frame+1,
                              tng_data->frame_set_n_frames);

            file_pos = ftell(tng_data->output_file);
            fseek(tng_data->output_file, 0, SEEK_END);
            output_file_len = ftell(tng_data->output_file);
            fseek(tng_data->output_file, (long)file_pos, SEEK_SET);

            /* Read mapping blocks from the last frame set */
            tng_block_init(&block);

            stat = tng_block_header_read(tng_data, block);
            while(file_pos < output_file_len &&
                  stat != TNG_CRITICAL &&
                  block->id != TNG_TRAJECTORY_FRAME_SET &&
                  block->id != -1)
            {
                if(block->id == TNG_PARTICLE_MAPPING)
                {
                    tng_trajectory_mapping_block_read(tng_data, block,
                                                      hash_mode);
                }
                else
                {
                    fseek(tng_data->output_file, (long)block->block_contents_size,
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
                fprintf(stderr, "TNG library: Error writing frame set. %s: %d\n", __FILE__, __LINE__);
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
    fseek(tng_data->output_file, (long)file_pos, SEEK_SET);

    /* Read past the frame set block first */
    stat = tng_block_header_read(tng_data, block);
    if(stat == TNG_CRITICAL)
    {
        fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
               file_pos, __FILE__, __LINE__);
        tng_block_destroy(&block);
        tng_data->input_file = temp;

        tng_data->first_trajectory_frame_set_input_file_pos = temp_first;
        tng_data->last_trajectory_frame_set_input_file_pos = temp_last;
        tng_data->current_trajectory_frame_set_input_file_pos = temp_current;
        return(stat);
    }
    fseek(tng_data->output_file, (long)block->block_contents_size,
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
                block->id != TNG_TRAJECTORY_FRAME_SET &&
                block->id != -1)
        {
            if(block->id == TNG_PARTICLE_MAPPING)
            {
                tng_trajectory_mapping_block_read(tng_data, block, hash_mode);
            }
            else
            {
                fseek(tng_data->output_file, (long)block->block_contents_size,
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
            fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
        fseek(tng_data->output_file, (long)mapping_block_end_pos, SEEK_SET);
    }

    /* Read all block headers until next frame set block or
     * until the wanted block id is found */
    stat = tng_block_header_read(tng_data, block);
    while(file_pos < output_file_len &&
            stat != TNG_CRITICAL &&
            block->id != block_id &&
            block->id != TNG_PARTICLE_MAPPING &&
            block->id != TNG_TRAJECTORY_FRAME_SET &&
            block->id != -1)
    {
        fseek(tng_data->output_file, (long)block->block_contents_size, SEEK_CUR);
        file_pos = ftell(tng_data->output_file);
        if(file_pos < output_file_len)
        {
            stat = tng_block_header_read(tng_data, block);
        }
    }
    if(stat == TNG_CRITICAL)
    {
        fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
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
        fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    data.datatype = datatype;

    if(fread(&dependency, sizeof(dependency), 1, tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
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
        fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }

    if(fread(&data.n_values_per_frame, sizeof(data.n_values_per_frame), 1,
             tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
            &data.n_values_per_frame)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(fread(&data.codec_id, sizeof(data.codec_id), 1,
             tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
            &data.codec_id)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(data.codec_id != TNG_UNCOMPRESSED)
    {
        if(fread(&data.compression_multiplier,
                 sizeof(data.compression_multiplier), 1, tng_data->input_file)
            == 0)
        {
            fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }

        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
               (int64_t *)&data.compression_multiplier)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                &data.first_frame_with_data)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }

        if(fread(&data.stride_length, sizeof(data.stride_length),
                 1, tng_data->input_file) == 0)
        {
            fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
            tng_block_destroy(&block);
            return(TNG_CRITICAL);
        }
        if(tng_data->output_endianness_swap_func_64)
        {
            if(tng_data->output_endianness_swap_func_64(tng_data,
                &data.stride_length)
                != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
        fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
            &num_first_particle)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(fread(&block_n_particles, sizeof(block_n_particles), 1,
             tng_data->input_file) == 0)
    {
        fprintf(stderr, "TNG library: Error reading file. %s: %d\n", __FILE__, __LINE__);
        tng_block_destroy(&block);
        return(TNG_CRITICAL);
    }
    if(tng_data->output_endianness_swap_func_64)
    {
        if(tng_data->output_endianness_swap_func_64(tng_data,
            &block_n_particles)
            != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
            fprintf(stderr, "TNG library: Cannot calculate writing locations. %s: %d.\n", __FILE__,
                   __LINE__);
            tng_block_destroy(&block);
            return(TNG_FAILURE);
    }

    n_values_per_frame = data.n_values_per_frame;

    file_pos = (frame_nr - tng_max_i64(frame_set->first_frame,
                               data.first_frame_with_data)) /
                data.stride_length;
    file_pos *= block_n_particles * size * n_values_per_frame;

    if(file_pos > contents_size)
    {
        fprintf(stderr, "TNG library: Attempting to write outside the block. %s: %d\n", __FILE__,
               __LINE__);
        tng_block_destroy(&block);
        return(TNG_FAILURE);
    }

    fseek(tng_data->output_file, (long)file_pos, SEEK_CUR);

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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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
                fprintf(stderr, "TNG library: Cannot swap byte order. %s: %d\n",
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

static tng_function_status tng_data_values_alloc
                (const tng_trajectory_t tng_data,
                 union data_values ***values,
                 const int64_t n_frames,
                 const int64_t n_values_per_frame,
                 const char type)
{
    int64_t i;
    tng_function_status stat;

    if(n_frames <= 0 || n_values_per_frame <= 0)
    {
        return(TNG_FAILURE);
    }

    if(*values)
    {
        stat = tng_data_values_free(tng_data, *values, n_frames,
                                    n_values_per_frame,
                                    type);
        if(stat != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot free particle data values. %s: %d\n",
                   __FILE__, __LINE__);
            return(stat);
        }
    }
    *values = malloc(sizeof(union data_values *) * n_frames);
    if(!*values)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(union data_values **) * n_frames,
               __FILE__, __LINE__);
        return(TNG_CRITICAL);

    }

    for(i = 0; i < n_frames; i++)
    {
        (*values)[i] = malloc(sizeof(union data_values) *
                           n_values_per_frame);
        if(!(*values)[i])
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   sizeof(union data_values) * n_values_per_frame,
                   __FILE__, __LINE__);
            free(values);
            values = 0;
            return(TNG_CRITICAL);
        }
    }
    return(TNG_SUCCESS);
}

/* FIXME: This needs ***values */
tng_function_status DECLSPECDLLEXPORT tng_data_values_free
                (const tng_trajectory_t tng_data,
                 union data_values **values,
                 const int64_t n_frames,
                 const int64_t n_values_per_frame,
                 const char type)
{
    int64_t i, j;
    (void)tng_data;

    if(values)
    {
        for(i = 0; i < n_frames; i++)
        {
            if(values[i])
            {
                if(type == TNG_CHAR_DATA)
                {
                    for(j = 0; j < n_values_per_frame; j++)
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

static tng_function_status tng_particle_data_values_alloc
                (const tng_trajectory_t tng_data,
                 union data_values ****values,
                 const int64_t n_frames,
                 const int64_t n_particles,
                 const int64_t n_values_per_frame,
                 const char type)
{
    int64_t i, j;
    tng_function_status stat;

    if(n_particles == 0 || n_values_per_frame == 0)
    {
        return(TNG_FAILURE);
    }

    if(*values)
    {
        stat = tng_particle_data_values_free(tng_data, *values, n_frames,
                                             n_particles, n_values_per_frame,
                                             type);
        if(stat != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot free particle data values. %s: %d\n",
                   __FILE__, __LINE__);
            return(stat);
        }
    }
    *values = malloc(sizeof(union data_values **) * n_frames);
    if(!*values)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(union data_values **) * n_frames,
               __FILE__, __LINE__);
        return(TNG_CRITICAL);

    }

    for(i = 0; i < n_frames; i++)
    {
        (*values)[i] = malloc(sizeof(union data_values *) *
                           n_particles);
        if(!(*values)[i])
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   sizeof(union data_values *) * n_particles,
                   __FILE__, __LINE__);
            free(*values);
            *values = 0;
            return(TNG_CRITICAL);
        }
        for(j = 0; j < n_particles; j++)
        {
            (*values)[i][j] = malloc(sizeof(union data_values) *
                                  n_values_per_frame);
            if(!(*values)[i][j])
            {
                fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                    sizeof(union data_values *) * n_particles,
                    __FILE__, __LINE__);
                tng_particle_data_values_free(tng_data, *values, n_frames,
                                              n_particles, n_values_per_frame,
                                              type);
                *values = 0;
                return(TNG_CRITICAL);
            }
        }
    }
    return(TNG_SUCCESS);
}

/* FIXME: This needs ****values */
tng_function_status DECLSPECDLLEXPORT tng_particle_data_values_free
                (const tng_trajectory_t tng_data,
                 union data_values ***values,
                 const int64_t n_frames,
                 const int64_t n_particles,
                 const int64_t n_values_per_frame,
                 const char type)
{
    int64_t i, j, k;
    (void)tng_data;

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
                        for(k = 0; k < n_values_per_frame; k++)
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


tng_function_status DECLSPECDLLEXPORT tng_data_get
                (tng_trajectory_t tng_data,
                 const int64_t block_id,
                 union data_values ***values,
                 int64_t *n_frames,
                 int64_t *n_values_per_frame,
                 char *type)
{
    int64_t i, j, file_pos, block_index;
    int size;
    size_t len;
    tng_non_particle_data_t data;
    tng_trajectory_frame_set_t frame_set;
    tng_gen_block_t block;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(n_frames, "TNG library: n_frames must not be a NULL pointer.");
    TNG_ASSERT(n_values_per_frame, "TNG library: n_values_per_frame must not be a NULL pointer.");
    TNG_ASSERT(type, "TNG library: type must not be a NULL pointer.");

    frame_set = &tng_data->current_trajectory_frame_set;

    block_index = -1;
    data = 0;

    if(tng_data_find(tng_data, block_id, &data) != TNG_SUCCESS)
    {
        tng_block_init(&block);
        file_pos = ftell(tng_data->input_file);
        /* Read all blocks until next frame set block */
        stat = tng_block_header_read(tng_data, block);
        while(file_pos < tng_data->input_file_len &&
                stat != TNG_CRITICAL &&
                block->id != TNG_TRAJECTORY_FRAME_SET &&
                block->id != -1)
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
            fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
                    file_pos, __FILE__, __LINE__);
            return(stat);
        }

        for(i = 0; i < frame_set->n_data_blocks; i++)
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

    *n_frames = tng_max_i64(1, data->n_frames);
    *n_values_per_frame = data->n_values_per_frame;
    *type = data->datatype;

    if(*values == 0)
    {
        if(tng_data_values_alloc(tng_data, values, *n_frames,
                                 *n_values_per_frame,
                                 *type)
        != TNG_SUCCESS)
        {
            return(TNG_CRITICAL);
        }
    }

    switch(*type)
    {
    case TNG_CHAR_DATA:
        for(i = 0; i < *n_frames; i++)
        {
            for(j = 0; j < *n_values_per_frame; j++)
            {
                len = strlen(data->strings[i][j]) + 1;
                (*values)[i][j].c = malloc(len);
                strncpy((*values)[i][j].c, data->strings[i][j], len);
            }
        }
        break;
    case TNG_INT_DATA:
        size = sizeof(int);
        for(i = 0; i < *n_frames; i++)
        {
            for(j = 0; j < *n_values_per_frame; j++)
            {
                (*values)[i][j].i = *(int *)((char *)data->values + size *
                                             (i*(*n_values_per_frame) + j));
            }
        }
        break;
    case TNG_FLOAT_DATA:
        size = sizeof(float);
        for(i = 0; i < *n_frames; i++)
        {
            for(j = 0; j < *n_values_per_frame; j++)
            {
                (*values)[i][j].f = *(float *)((char *)data->values + size *
                                               (i*(*n_values_per_frame) + j));
            }
        }
        break;
    case TNG_DOUBLE_DATA:
    default:
        size = sizeof(double);
        for(i = 0; i < *n_frames; i++)
        {
            for(j = 0; j < *n_values_per_frame; j++)
            {
                (*values)[i][j].d = *(double *)((char *)data->values + size *
                                                (i*(*n_values_per_frame) + j));
            }
        }
    }

    data->last_retrieved_frame = frame_set->first_frame + data->n_frames - 1;

    return(TNG_SUCCESS);
}

tng_function_status tng_data_vector_get(tng_trajectory_t tng_data,
                                        const int64_t block_id,
                                        void **values,
                                        int64_t *n_frames,
                                        int64_t *stride_length,
                                        int64_t *n_values_per_frame,
                                        char *type)
{
    int64_t file_pos, data_size, n_frames_div, block_index;
    int i, size;
    tng_non_particle_data_t data;
    tng_trajectory_frame_set_t frame_set;
    tng_gen_block_t block;
    void *temp;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(n_frames, "TNG library: n_frames must not be a NULL pointer.");
    TNG_ASSERT(stride_length, "TNG library: stride_length must not be a NULL pointer.");
    TNG_ASSERT(n_values_per_frame, "TNG library: n_values_per_frame must not be a NULL pointer.");
    TNG_ASSERT(type, "TNG library: type must not be a NULL pointer.");

    frame_set = &tng_data->current_trajectory_frame_set;

    block_index = -1;
    data = 0;

    if(tng_data_find(tng_data, block_id, &data) != TNG_SUCCESS)
    {
        tng_block_init(&block);
        file_pos = ftell(tng_data->input_file);
        /* Read all blocks until next frame set block */
        stat = tng_block_header_read(tng_data, block);
        while(file_pos < tng_data->input_file_len &&
                stat != TNG_CRITICAL &&
                block->id != TNG_TRAJECTORY_FRAME_SET &&
                block->id != -1)
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
            fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
                    file_pos, __FILE__, __LINE__);
            return(stat);
        }

        for(i = 0; i < frame_set->n_data_blocks; i++)
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

    *type = data->datatype;

    switch(*type)
    {
    case TNG_CHAR_DATA:
        return(TNG_FAILURE);
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

    *n_frames = data->n_frames;
    *n_values_per_frame = data->n_values_per_frame;
    *stride_length = data->stride_length;
    n_frames_div = (*n_frames % *stride_length) ? *n_frames / *stride_length + 1:
                   *n_frames / *stride_length;

    data_size = n_frames_div * size *
                *n_values_per_frame;

    temp = realloc(*values, data_size);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               data_size, __FILE__, __LINE__);
        free(*values);
        *values = 0;
        return(TNG_CRITICAL);
    }

    *values = temp;

    memcpy(*values, data->values, data_size);

    data->last_retrieved_frame = frame_set->first_frame + data->n_frames - 1;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_data_interval_get
                (tng_trajectory_t tng_data,
                 const int64_t block_id,
                 const int64_t start_frame_nr,
                 const int64_t end_frame_nr,
                 const char hash_mode,
                 union data_values ***values,
                 int64_t *n_values_per_frame,
                 char *type)
{
    int64_t i, j, n_frames, file_pos, current_frame_pos, first_frame;
    int64_t block_index;
    int size;
    size_t len;
    tng_non_particle_data_t data;
    tng_trajectory_frame_set_t frame_set;
    tng_gen_block_t block;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(start_frame_nr <= end_frame_nr, "TNG library: start_frame_nr must not be higher than tne end_frame_nr.");
    TNG_ASSERT(n_values_per_frame, "TNG library: n_values_per_frame must not be a NULL pointer.");
    TNG_ASSERT(type, "TNG library: type must not be a NULL pointer.");

    block_index = -1;

    frame_set = &tng_data->current_trajectory_frame_set;
    first_frame = frame_set->first_frame;

    stat = tng_frame_set_of_frame_find(tng_data, start_frame_nr);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }


    /* Do not re-read the frame set. */
    if(first_frame != frame_set->first_frame ||
       frame_set->n_data_blocks <= 0)
    {
        tng_block_init(&block);
        file_pos = ftell(tng_data->input_file);
        /* Read all blocks until next frame set block */
        stat = tng_block_header_read(tng_data, block);
        while(file_pos < tng_data->input_file_len &&
            stat != TNG_CRITICAL &&
            block->id != TNG_TRAJECTORY_FRAME_SET &&
            block->id != -1)
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
            fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
                    file_pos, __FILE__, __LINE__);
            return(stat);
        }
    }


    /* See if there is a data block of this ID.
     * Start checking the last read frame set */
    for(i = 0; i < frame_set->n_data_blocks; i++)
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
        fprintf(stderr, "TNG library: Could not find non-particle data block with id %"PRId64". %s: %d\n",
                block_id, __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    n_frames = end_frame_nr - start_frame_nr + 1;
    *n_values_per_frame = data->n_values_per_frame;
    *type = data->datatype;

    if(*values == 0)
    {
        if(tng_data_values_alloc(tng_data, values, n_frames,
                                 *n_values_per_frame,
                                 *type) != TNG_SUCCESS)
        {
            return(TNG_CRITICAL);
        }
    }

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
            for(j = 0; j < *n_values_per_frame; j++)
            {
                len = strlen(data->strings[current_frame_pos][j]) + 1;
                (*values)[i][j].c = malloc(len);
                strncpy((*values)[i][j].c, data->strings[current_frame_pos][j], len);
            }
            current_frame_pos++;
        }
        break;
    case TNG_INT_DATA:
        size = sizeof(int);
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
            for(j = 0; j < *n_values_per_frame; j++)
            {
                (*values)[i][j].i = *(int *)((char *)data->values + size *
                                            (current_frame_pos *
                                             (*n_values_per_frame) + j));
            }
            current_frame_pos++;
        }
        break;
    case TNG_FLOAT_DATA:
        size = sizeof(float);
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
            for(j = 0; j < *n_values_per_frame; j++)
            {
                (*values)[i][j].f = *(float *)((char *)data->values + size *
                                               (current_frame_pos *
                                                (*n_values_per_frame) + j));
            }
            current_frame_pos++;
        }
        break;
    case TNG_DOUBLE_DATA:
    default:
        size = sizeof(double);
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
            for(j = 0; j < *n_values_per_frame; j++)
            {
                (*values)[i][j].d = *(double *)((char *)data->values + size *
                                                (current_frame_pos *
                                                 (*n_values_per_frame) + j));
            }
            current_frame_pos++;
        }
    }

    data->last_retrieved_frame = end_frame_nr;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_data_vector_interval_get
                (tng_trajectory_t tng_data,
                 const int64_t block_id,
                 const int64_t start_frame_nr,
                 const int64_t end_frame_nr,
                 const char hash_mode,
                 void **values,
                 int64_t *stride_length,
                 int64_t *n_values_per_frame,
                 char *type)
{
    int64_t n_frames, tot_n_frames, n_frames_div, n_frames_div_2, first_frame;
    int64_t file_pos, current_frame_pos, data_size, frame_size;
    int64_t last_frame_pos;
    int size;
    tng_trajectory_frame_set_t frame_set;
    tng_non_particle_data_t np_data;
    tng_gen_block_t block;
    void *current_values = 0, *temp;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(start_frame_nr <= end_frame_nr, "TNG library: start_frame_nr must not be higher than the end_frame_nr.");
    TNG_ASSERT(stride_length, "TNG library: stride_length must not be a NULL pointer.");
    TNG_ASSERT(n_values_per_frame, "TNG library: n_values_per_frame must not be a NULL pointer.");
    TNG_ASSERT(type, "TNG library: type must not be a NULL pointer.");

    frame_set = &tng_data->current_trajectory_frame_set;
    first_frame = frame_set->first_frame;

    stat = tng_frame_set_of_frame_find(tng_data, start_frame_nr);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    /* Do not re-read the frame set and only need the requested block. */
    /* TODO: Test that blocks are read correctly now that not all of them are read at the same time. */
    stat = tng_data_find(tng_data, block_id, &np_data);
    if(first_frame != frame_set->first_frame ||
       stat != TNG_SUCCESS)
    {
        tng_block_init(&block);
        if(stat != TNG_SUCCESS)
        {
            fseek(tng_data->input_file,
                  (long)tng_data->current_trajectory_frame_set_input_file_pos,
                  SEEK_SET);
            stat = tng_block_header_read(tng_data, block);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot read block header. %s: %d\n",
                        __FILE__, __LINE__);
                return(stat);
            }

            fseek(tng_data->input_file, (long)block->block_contents_size, SEEK_CUR);
        }
        file_pos = ftell(tng_data->input_file);
        /* Read until next frame set block */
        stat = tng_block_header_read(tng_data, block);
        while(file_pos < tng_data->input_file_len &&
            stat != TNG_CRITICAL &&
            block->id != TNG_TRAJECTORY_FRAME_SET &&
            block->id != -1)
        {
            if(block->id == block_id)
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
            else
            {
                file_pos += block->block_contents_size + block->header_contents_size;
                fseek(tng_data->input_file, (long)block->block_contents_size, SEEK_CUR);
                if(file_pos < tng_data->input_file_len)
                {
                    stat = tng_block_header_read(tng_data, block);
                }
            }
        }
        tng_block_destroy(&block);
        if(stat == TNG_CRITICAL)
        {
            fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
                    file_pos, __FILE__, __LINE__);
            return(stat);
        }
    }

    stat = tng_data_find(tng_data, block_id, &np_data);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    stat = tng_data_vector_get(tng_data, block_id, &current_values,
                               &n_frames, stride_length,
                               n_values_per_frame, type);

    if(stat != TNG_SUCCESS)
    {
        if(current_values)
        {
            free(current_values);
        }
        return(stat);
    }

    if(n_frames == 1 && n_frames < frame_set->n_frames)
    {
        tot_n_frames = 1;
    }
    else
    {
        tot_n_frames = end_frame_nr - start_frame_nr + 1;
    }

    switch(*type)
    {
    case TNG_CHAR_DATA:
        return(TNG_FAILURE);
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

    n_frames_div = (tot_n_frames % *stride_length) ?
                 tot_n_frames / *stride_length + 1:
                 tot_n_frames / *stride_length;
    data_size = n_frames_div * size * (*n_values_per_frame);

/*     fprintf(stderr, "TNG library: size: %d, n_frames_div: %"PRId64", data_size: %"PRId64"\n",
              size, n_frames_div, data_size);
*/
    temp = realloc(*values, data_size);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               data_size, __FILE__, __LINE__);
        free(*values);
        *values = 0;
        return(TNG_CRITICAL);
    }

    *values = temp;

    if( n_frames == 1 && n_frames < frame_set->n_frames)
    {
        memcpy(*values, current_values, size * (*n_values_per_frame));
    }
    else
    {
        current_frame_pos = start_frame_nr - frame_set->first_frame;

        frame_size = size * (*n_values_per_frame);

        last_frame_pos = tng_min_i64(n_frames,
                                     end_frame_nr - start_frame_nr);

        n_frames_div = current_frame_pos / *stride_length;
        n_frames_div_2 = (last_frame_pos % *stride_length) ?
                       last_frame_pos / *stride_length + 1:
                       last_frame_pos / *stride_length;
        n_frames_div_2 = tng_max_i64(1, n_frames_div_2);

        memcpy(*values, (char *)current_values + n_frames_div * frame_size,
               n_frames_div_2 * frame_size);

        current_frame_pos += n_frames - current_frame_pos;

        while(current_frame_pos <= end_frame_nr - start_frame_nr)
        {
            stat = tng_frame_set_read_next(tng_data, hash_mode);
            if(stat != TNG_SUCCESS)
            {
                if(current_values)
                {
                    free(current_values);
                }
                free(*values);
                *values = 0;
                return(stat);
            }

            stat = tng_data_vector_get(tng_data, block_id, &current_values,
                                    &n_frames, stride_length,
                                    n_values_per_frame, type);

            if(stat != TNG_SUCCESS)
            {
                if(current_values)
                {
                    free(current_values);
                }
                free(*values);
                *values = 0;
                return(stat);
            }

            last_frame_pos = tng_min_i64(n_frames,
                                         end_frame_nr - current_frame_pos);

            n_frames_div = current_frame_pos / *stride_length;
            n_frames_div_2 = (last_frame_pos % *stride_length) ?
                           last_frame_pos / *stride_length + 1:
                           last_frame_pos / *stride_length;
            n_frames_div_2 = tng_max_i64(1, n_frames_div_2);

            memcpy(((char *)*values) + n_frames_div * frame_size,
                   current_values,
                   n_frames_div_2 * frame_size);

            current_frame_pos += n_frames;
        }
    }

    if(current_values)
    {
        free(current_values);
    }

    np_data->last_retrieved_frame = end_frame_nr;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_particle_data_get
                (tng_trajectory_t tng_data,
                 const int64_t block_id,
                 union data_values ****values,
                 int64_t *n_frames,
                 int64_t *n_particles,
                 int64_t *n_values_per_frame,
                 char *type)
{
    int64_t i, j, k, mapping, file_pos, i_step, block_index;
    int size;
    size_t len;
    tng_particle_data_t data;
    tng_trajectory_frame_set_t frame_set;
    tng_gen_block_t block;
    char block_type_flag;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(n_frames, "TNG library: n_frames must not be a NULL pointer.");
    TNG_ASSERT(n_particles, "TNG library: n_particles must not be a NULL pointer.");
    TNG_ASSERT(n_values_per_frame, "TNG library: n_values_per_frame must not be a NULL pointer.");
    TNG_ASSERT(type, "TNG library: type must not be a NULL pointer.");

    frame_set = &tng_data->current_trajectory_frame_set;

    block_index = -1;
    data = 0;

    if(tng_particle_data_find(tng_data, block_id, &data) != TNG_SUCCESS)
    {
        if(tng_data->current_trajectory_frame_set_input_file_pos > 0)
        {
            block_type_flag = TNG_TRAJECTORY_BLOCK;
        }
        else
        {
            block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
        }

        tng_block_init(&block);
        file_pos = ftell(tng_data->input_file);
        /* Read all blocks until next frame set block */
        stat = tng_block_header_read(tng_data, block);
        while(file_pos < tng_data->input_file_len &&
                stat != TNG_CRITICAL &&
                block->id != TNG_TRAJECTORY_FRAME_SET &&
                block->id != -1)
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
            fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
                    file_pos, __FILE__, __LINE__);
            return(stat);
        }

        for(i = 0; i < frame_set->n_particle_data_blocks; i++)
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
    else
    {
        if(tng_data->current_trajectory_frame_set_input_file_pos > 0)
        {
            block_type_flag = TNG_TRAJECTORY_BLOCK;
        }
        else
        {
            block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
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

    *n_frames = tng_max_i64(1, data->n_frames);
    *n_values_per_frame = data->n_values_per_frame;
    *type = data->datatype;

    if(*values == 0)
    {
        if(tng_particle_data_values_alloc(tng_data, values, *n_frames,
                                         *n_particles, *n_values_per_frame,
                                         *type)
            != TNG_SUCCESS)
        {
            return(TNG_CRITICAL);
        }
    }

    /* It's not very elegant to reuse so much of the code in the different case
     * statements, but it's unnecessarily slow to have the switch-case block
     * inside the for loops. */
    switch(*type)
    {
    case TNG_CHAR_DATA:
        for(i = 0; i < *n_frames; i++)
        {
            for(j = 0; j < *n_particles; j++)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k = 0; k < *n_values_per_frame; k++)
                {
                    len = strlen(data->strings[i][j][k]) + 1;
                    (*values)[i][mapping][k].c = malloc(len);
                    strncpy((*values)[i][mapping][k].c,
                            data->strings[i][j][k], len);
                }
            }
        }
        break;
    case TNG_INT_DATA:
        size = sizeof(int);
        i_step = (*n_particles) * (*n_values_per_frame);
        for(i = 0; i < *n_frames; i++)
        {
            for(j = 0; j < *n_particles; j++)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k = 0; k < *n_values_per_frame; k++)
                {
                    (*values)[i][mapping][k].i = *(int *)
                                                 ((char *)data->values + size *
                                                 (i * i_step + j *
                                                  (*n_values_per_frame) + k));
                }
            }
        }
        break;
    case TNG_FLOAT_DATA:
        size = sizeof(float);
        i_step = (*n_particles) * (*n_values_per_frame);
        for(i = 0; i < *n_frames; i++)
        {
            for(j = 0; j < *n_particles; j++)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k = 0; k < *n_values_per_frame; k++)
                {
                    (*values)[i][mapping][k].f = *(float *)
                                                 ((char *)data->values + size *
                                                 (i * i_step + j *
                                                  (*n_values_per_frame) + k));
                }
            }
        }
        break;
    case TNG_DOUBLE_DATA:
    default:
        size = sizeof(double);
        i_step = (*n_particles) * (*n_values_per_frame);
        for(i = 0; i < *n_frames; i++)
        {
            for(j = 0; j < *n_particles; j++)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k = 0; k < *n_values_per_frame; k++)
                {
                    (*values)[i][mapping][k].d = *(double *)
                                                 ((char *)data->values + size *
                                                 (i * i_step + j *
                                                  (*n_values_per_frame) + k));
                }
            }
        }
    }

    data->last_retrieved_frame = frame_set->first_frame + data->n_frames - 1;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_particle_data_vector_get
                (tng_trajectory_t tng_data,
                 const int64_t block_id,
                 void **values,
                 int64_t *n_frames,
                 int64_t *stride_length,
                 int64_t *n_particles,
                 int64_t *n_values_per_frame,
                 char *type)
{
    int64_t i, j, mapping, file_pos, i_step, data_size, n_frames_div;
    int64_t block_index;
    int size;
    tng_particle_data_t data;
    tng_trajectory_frame_set_t frame_set;
    tng_gen_block_t block;
    void *temp;
    char block_type_flag;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(n_particles, "TNG library: n_particles must not be a NULL pointer.");
    TNG_ASSERT(stride_length, "TNG library: stride_length must not be a NULL pointer.");
    TNG_ASSERT(n_values_per_frame, "TNG library: n_values_per_frame must not be a NULL pointer.");
    TNG_ASSERT(type, "TNG library: type must not be a NULL pointer.");

    frame_set = &tng_data->current_trajectory_frame_set;

    block_index = -1;
    data = 0;

    if(tng_particle_data_find(tng_data, block_id, &data) != TNG_SUCCESS)
    {
        tng_block_init(&block);
        file_pos = ftell(tng_data->input_file);
        /* Read all blocks until next frame set block */
        stat = tng_block_header_read(tng_data, block);
        while(file_pos < tng_data->input_file_len &&
                stat != TNG_CRITICAL &&
                block->id != TNG_TRAJECTORY_FRAME_SET &&
                block->id != -1)
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
            fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
                    file_pos, __FILE__, __LINE__);
            return(stat);
        }

        for(i = 0; i < frame_set->n_particle_data_blocks; i++)
        {
            data = &frame_set->tr_particle_data[i];
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

    if(tng_data->current_trajectory_frame_set_input_file_pos > 0)
    {
        block_type_flag = TNG_TRAJECTORY_BLOCK;
    }
    else
    {
        block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
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

    *type = data->datatype;

    switch(*type)
    {
    case TNG_CHAR_DATA:
        return(TNG_FAILURE);
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

    *n_frames = tng_max_i64(1, data->n_frames);
    *n_values_per_frame = data->n_values_per_frame;
    *stride_length = data->stride_length;

    n_frames_div = (*n_frames % *stride_length) ?
                   *n_frames / *stride_length + 1:
                   *n_frames / *stride_length;

    data_size = n_frames_div * size * (*n_particles) *
                (*n_values_per_frame);

    temp = realloc(*values, data_size);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               data_size, __FILE__, __LINE__);
        free(*values);
        *values = 0;
        return(TNG_CRITICAL);
    }

    *values = temp;

    if(frame_set->n_mapping_blocks <= 0)
    {
        memcpy(*values, data->values, data_size);
    }
    else
    {
        i_step = (*n_particles) * (*n_values_per_frame);
        for(i = 0; i < *n_frames; i++)
        {
            for(j = 0; j < *n_particles; j++)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                memcpy(((char *)*values) + size * (i * i_step + mapping *
                       (*n_values_per_frame)),
                       (char *)data->values + size *
                       (i * i_step + j * (*n_values_per_frame)),
                       size * (*n_values_per_frame));
            }
        }
    }

    data->last_retrieved_frame = frame_set->first_frame + data->n_frames - 1;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_particle_data_interval_get
                (tng_trajectory_t tng_data,
                 const int64_t block_id,
                 const int64_t start_frame_nr,
                 const int64_t end_frame_nr,
                 const char hash_mode,
                 union data_values ****values,
                 int64_t *n_particles,
                 int64_t *n_values_per_frame,
                 char *type)
{
    int64_t i, j, k, mapping, n_frames, file_pos, current_frame_pos, i_step;
    int64_t first_frame, block_index;
    int size;
    size_t len;
    tng_particle_data_t data;
    tng_trajectory_frame_set_t frame_set;
    tng_gen_block_t block;
    char block_type_flag;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(start_frame_nr <= end_frame_nr, "TNG library: start_frame_nr must not be higher than tne end_frame_nr.");
    TNG_ASSERT(n_particles, "TNG library: n_particles must not be a NULL pointer.");
    TNG_ASSERT(n_values_per_frame, "TNG library: n_values_per_frame must not be a NULL pointer.");
    TNG_ASSERT(type, "TNG library: type must not be a NULL pointer.");

    block_index = -1;

    frame_set = &tng_data->current_trajectory_frame_set;
    first_frame = frame_set->first_frame;

    stat = tng_frame_set_of_frame_find(tng_data, start_frame_nr);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    /* Do not re-read the frame set. */
    if(first_frame != frame_set->first_frame ||
       frame_set->n_particle_data_blocks <= 0)
    {
        tng_block_init(&block);
        file_pos = ftell(tng_data->input_file);
        /* Read all blocks until next frame set block */
        stat = tng_block_header_read(tng_data, block);
        while(file_pos < tng_data->input_file_len &&
                stat != TNG_CRITICAL &&
                block->id != TNG_TRAJECTORY_FRAME_SET &&
                block->id != -1)
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
            fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
                    file_pos, __FILE__, __LINE__);
            return(stat);
        }
    }

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
        fprintf(stderr, "TNG library: Could not find particle data block with id %"PRId64". %s: %d\n",
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
    *n_values_per_frame = data->n_values_per_frame;
    *type = data->datatype;

    if(*values == 0)
    {
        if(tng_particle_data_values_alloc(tng_data, values, n_frames,
                                         *n_particles, *n_values_per_frame,
                                         *type)
            != TNG_SUCCESS)
        {
            return(TNG_CRITICAL);
        }
    }

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
            for(j = 0; j < *n_particles; j++)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k = 0; k < *n_values_per_frame; k++)
                {
                    len = strlen(data->strings[current_frame_pos][j][k]) + 1;
                    (*values)[i][mapping][k].c = malloc(len);
                    strncpy((*values)[i][mapping][k].c, data->strings[current_frame_pos][j][k], len);
                }
            }
            current_frame_pos++;
        }
        break;
    case TNG_INT_DATA:
        size = sizeof(int);
        i_step = (*n_particles) * (*n_values_per_frame);
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
            for(j = 0; j < *n_particles; j++)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k = 0; k < *n_values_per_frame; k++)
                {
                    (*values)[i][mapping][k].i = *(int *)
                                                 ((char *)data->values + size *
                                                  (current_frame_pos *
                                                   i_step + j *
                                                   (*n_values_per_frame) + k));
                }
            }
            current_frame_pos++;
        }
        break;
    case TNG_FLOAT_DATA:
        size = sizeof(float);
        i_step = (*n_particles) * (*n_values_per_frame);
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
            for(j=0; j<*n_particles; j++)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k=0; k<*n_values_per_frame; k++)
                {
                    (*values)[i][mapping][k].f = *(float *)
                                                 ((char *)data->values + size *
                                                  (current_frame_pos *
                                                   i_step + j *
                                                   (*n_values_per_frame) + k));
                }
            }
            current_frame_pos++;
        }
        break;
    case TNG_DOUBLE_DATA:
    default:
        size = sizeof(double);
        i_step = (*n_particles) * (*n_values_per_frame);
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
            for(j=0; j<*n_particles; j++)
            {
                tng_particle_mapping_get_real_particle(frame_set, j, &mapping);
                for(k=0; k<*n_values_per_frame; k++)
                {
                    (*values)[i][mapping][k].d = *(double *)
                                                 ((char *)data->values + size *
                                                  (current_frame_pos *
                                                   i_step + j *
                                                   (*n_values_per_frame) + k));
                }
            }
            current_frame_pos++;
        }
    }

    data->last_retrieved_frame = end_frame_nr;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_particle_data_vector_interval_get
                (tng_trajectory_t tng_data,
                 const int64_t block_id,
                 const int64_t start_frame_nr,
                 const int64_t end_frame_nr,
                 const char hash_mode,
                 void **values,
                 int64_t *n_particles,
                 int64_t *stride_length,
                 int64_t *n_values_per_frame,
                 char *type)
{
    int64_t n_frames, tot_n_frames, n_frames_div, n_frames_div_2, first_frame;
    int64_t file_pos, current_frame_pos, last_frame_pos, data_size, frame_size;
    int size;
    tng_trajectory_frame_set_t frame_set;
    tng_particle_data_t p_data;
    tng_gen_block_t block;
    void *current_values = 0, *temp;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(start_frame_nr <= end_frame_nr, "TNG library: start_frame_nr must not be higher than tne end_frame_nr.");
    TNG_ASSERT(n_particles, "TNG library: n_particles must not be a NULL pointer.");
    TNG_ASSERT(stride_length, "TNG library: stride_length must not be a NULL pointer.");
    TNG_ASSERT(n_values_per_frame, "TNG library: n_values_per_frame must not be a NULL pointer.");
    TNG_ASSERT(type, "TNG library: type must not be a NULL pointer.");

    frame_set = &tng_data->current_trajectory_frame_set;
    first_frame = frame_set->first_frame;

    stat = tng_frame_set_of_frame_find(tng_data, start_frame_nr);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    /* Do not re-read the frame set and only need the requested block + particle mapping blocks. */
    /* TODO: Test that blocks are read correctly now that now all of them are read at the same time. */
    stat = tng_particle_data_find(tng_data, block_id, &p_data);
    if(first_frame != frame_set->first_frame ||
       stat != TNG_SUCCESS)
    {
        tng_block_init(&block);
        if(stat != TNG_SUCCESS)
        {
            fseek(tng_data->input_file,
                  (long)tng_data->current_trajectory_frame_set_input_file_pos,
                  SEEK_SET);
            stat = tng_block_header_read(tng_data, block);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot read block header. %s: %d\n",
                        __FILE__, __LINE__);
                return(stat);
            }

            fseek(tng_data->input_file, (long)block->block_contents_size, SEEK_CUR);
        }
        file_pos = ftell(tng_data->input_file);
        /* Read until next frame set block */
        stat = tng_block_header_read(tng_data, block);
        while(file_pos < tng_data->input_file_len &&
            stat != TNG_CRITICAL &&
            block->id != TNG_TRAJECTORY_FRAME_SET &&
            block->id != -1)
        {
            if(block->id == block_id || block->id == TNG_PARTICLE_MAPPING)
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
            else
            {
                file_pos += block->block_contents_size + block->header_contents_size;
                fseek(tng_data->input_file, (long)block->block_contents_size, SEEK_CUR);
                if(file_pos < tng_data->input_file_len)
                {
                    stat = tng_block_header_read(tng_data, block);
                }
            }
        }
        tng_block_destroy(&block);
        if(stat == TNG_CRITICAL)
        {
            fprintf(stderr, "TNG library: Cannot read block header at pos %"PRId64". %s: %d\n",
                    file_pos, __FILE__, __LINE__);
            return(stat);
        }
    }
    stat = tng_particle_data_find(tng_data, block_id, &p_data);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    stat = tng_particle_data_vector_get(tng_data, block_id, &current_values,
                                        &n_frames, stride_length, n_particles,
                                        n_values_per_frame, type);

    if(stat != TNG_SUCCESS || *n_particles == 0)
    {
        if(current_values)
        {
            free(current_values);
        }
        return(stat);
    }

    if(n_frames == 1 && n_frames < frame_set->n_frames)
    {
        tot_n_frames = 1;
    }
    else
    {
        tot_n_frames = end_frame_nr - start_frame_nr + 1;
    }

    switch(*type)
    {
    case TNG_CHAR_DATA:
        return(TNG_FAILURE);
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

    n_frames_div = (tot_n_frames % *stride_length) ?
                 tot_n_frames / *stride_length + 1:
                 tot_n_frames / *stride_length;

    data_size = n_frames_div * size * (*n_particles) *
                (*n_values_per_frame);

    temp = realloc(*values, data_size);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               data_size, __FILE__, __LINE__);
        free(*values);
        *values = 0;
        return(TNG_CRITICAL);
    }

    *values = temp;

    if( n_frames == 1 && n_frames < frame_set->n_frames)
    {
        memcpy(*values, current_values, size * (*n_particles) *
               (*n_values_per_frame));
    }
    else
    {
        current_frame_pos = start_frame_nr - frame_set->first_frame;

        frame_size = size * (*n_particles) * (*n_values_per_frame);

        last_frame_pos = tng_min_i64(n_frames,
                                     end_frame_nr - start_frame_nr);

        n_frames_div = current_frame_pos / *stride_length;
        n_frames_div_2 = (last_frame_pos % *stride_length) ?
                       last_frame_pos / *stride_length + 1:
                       last_frame_pos / *stride_length;
        n_frames_div_2 = tng_max_i64(1, n_frames_div_2 + 1);

        memcpy(*values, (char *)current_values + n_frames_div * frame_size,
               n_frames_div_2 * frame_size);

        current_frame_pos += n_frames - current_frame_pos;

        while(current_frame_pos <= end_frame_nr - start_frame_nr)
        {
            stat = tng_frame_set_read_next(tng_data, hash_mode);
            if(stat != TNG_SUCCESS)
            {
                if(current_values)
                {
                    free(current_values);
                }
                free(*values);
                *values = 0;
                return(stat);
            }

            stat = tng_particle_data_vector_get(tng_data, block_id, &current_values,
                                                &n_frames, stride_length, n_particles,
                                                n_values_per_frame, type);

            if(stat != TNG_SUCCESS)
            {
                if(current_values)
                {
                    free(current_values);
                }
                free(*values);
                *values = 0;
                return(stat);
            }

            last_frame_pos = tng_min_i64(n_frames,
                                         end_frame_nr - current_frame_pos);

            n_frames_div = current_frame_pos / *stride_length;
            n_frames_div_2 = (last_frame_pos % *stride_length) ?
                           last_frame_pos / *stride_length + 1:
                           last_frame_pos / *stride_length;
            n_frames_div_2 = tng_max_i64(1, n_frames_div_2);

            memcpy(((char *)*values) + n_frames_div * frame_size,
                   current_values,
                   n_frames_div_2 * frame_size);

            current_frame_pos += n_frames;
        }
    }

    if(current_values)
    {
        free(current_values);
    }

    p_data->last_retrieved_frame = end_frame_nr;

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_data_get_stride_length
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 int64_t frame,
                 int64_t *stride_length)
{
    tng_function_status stat;
    tng_non_particle_data_t np_data;
    tng_particle_data_t p_data;
    long orig_file_pos, file_pos;
    int is_particle_data;

    if(tng_data->current_trajectory_frame_set_input_file_pos <= 0)
    {
        frame = 0;
    }

    if(frame >= 0)
    {
        stat = tng_frame_set_of_frame_find(tng_data, frame);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
    }
    orig_file_pos = tng_data->current_trajectory_frame_set_input_file_pos;
    stat = tng_data_find(tng_data, block_id, &np_data);
    if(stat != TNG_SUCCESS)
    {
        stat = tng_particle_data_find(tng_data, block_id, &p_data);
        if(stat != TNG_SUCCESS)
        {
            stat = tng_frame_set_read_current_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
            /* If no specific frame was required read until this data block is found */
            if(frame < 0)
            {
                file_pos = ftell(tng_data->input_file);
                while(stat != TNG_SUCCESS && file_pos < tng_data->input_file_len)
                {
                    stat = tng_frame_set_read_next_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
                    file_pos = ftell(tng_data->input_file);
                }
            }
            if(stat != TNG_SUCCESS)
            {
                tng_reread_frame_set_at_file_pos(tng_data, orig_file_pos);

                return(stat);
            }
            stat = tng_data_find(tng_data, block_id, &np_data);
            if(stat != TNG_SUCCESS)
            {
                stat = tng_particle_data_find(tng_data, block_id, &p_data);
                if(stat != TNG_SUCCESS)
                {
                    tng_reread_frame_set_at_file_pos(tng_data, orig_file_pos);

                    return(stat);
                }
                else
                {
                    is_particle_data = 1;
                }
            }
            else
            {
                is_particle_data = 0;
            }
        }
        else
        {
            is_particle_data = 1;
        }
    }
    else
    {
        is_particle_data = 0;
    }
    if(is_particle_data)
    {
        *stride_length = p_data->stride_length;
    }
    else
    {
        *stride_length = np_data->stride_length;
    }
    tng_reread_frame_set_at_file_pos(tng_data, orig_file_pos);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_time_get_str
                (const tng_trajectory_t tng_data,
                 char *time)
{
    struct tm *time_data;
    time_t secs;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(time, "TNG library: time must not be a NULL pointer");

    secs = tng_data->time;

    time_data = localtime(&secs); /* Returns a statically allocated variable. */
    TNG_SNPRINTF(time, TNG_MAX_DATE_STR_LEN,
             "%4d-%02d-%02d %02d:%02d:%02d",
             time_data->tm_year+1900, time_data->tm_mon+1, time_data->tm_mday,
             time_data->tm_hour, time_data->tm_min, time_data->tm_sec);

    return(TNG_SUCCESS);
}


tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_open
                (const char *filename,
                 const char mode,
                 tng_trajectory_t *tng_data_p)
{
    tng_function_status stat;

    TNG_ASSERT(filename, "TNG library: filename must not be a NULL pointer.");

    if(mode != 'r' && mode != 'w' && mode != 'a')
    {
        return(TNG_FAILURE);
    }

    if(tng_trajectory_init(tng_data_p) != TNG_SUCCESS)
    {
        tng_trajectory_destroy(tng_data_p);
        return(TNG_CRITICAL);
    }

    if(mode == 'r' || mode == 'a')
    {
        tng_input_file_set(*tng_data_p, filename);

        /* Read the file headers */
        tng_file_headers_read(*tng_data_p, TNG_USE_HASH);

        stat = tng_num_frame_sets_get(*tng_data_p, &(*tng_data_p)->n_trajectory_frame_sets);

        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
    }

    if(mode == 'w')
    {
        tng_output_file_set(*tng_data_p, filename);
    }
    else if(mode == 'a')
    {
        if((*tng_data_p)->output_file)
        {
            fclose((*tng_data_p)->output_file);
        }
        (*tng_data_p)->output_file = (*tng_data_p)->input_file;
        fseek((*tng_data_p)->input_file,
                (long)(*tng_data_p)->last_trajectory_frame_set_input_file_pos,
                SEEK_SET);

        stat = tng_frame_set_read(*tng_data_p, TNG_USE_HASH);
        if(stat != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot read frame set and related blocks. %s: %d\n",
                   __FILE__, __LINE__);
        }
        (*tng_data_p)->output_file = 0;

        (*tng_data_p)->first_trajectory_frame_set_output_file_pos =
        (*tng_data_p)->first_trajectory_frame_set_input_file_pos;
        (*tng_data_p)->last_trajectory_frame_set_output_file_pos =
        (*tng_data_p)->last_trajectory_frame_set_input_file_pos;
        (*tng_data_p)->current_trajectory_frame_set_output_file_pos =
        (*tng_data_p)->current_trajectory_frame_set_input_file_pos;
        if((*tng_data_p)->input_file)
        {
            fclose((*tng_data_p)->input_file);
            (*tng_data_p)->input_file = 0;
        }
        if((*tng_data_p)->input_file_path)
        {
            free((*tng_data_p)->input_file_path);
            (*tng_data_p)->input_file_path = 0;
        }
        tng_output_append_file_set(*tng_data_p, filename);

        fseek((*tng_data_p)->output_file, 0, SEEK_END);
    }

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_close
                (tng_trajectory_t *tng_data_p)
{
    tng_trajectory_frame_set_t frame_set;

    if(tng_data_p == 0)
    {
        fprintf(stderr, "TNG library: Empty pointer to trajectory when attempting to close. %s: %d\n",
               __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    if(*tng_data_p == 0)
    {
        return(TNG_SUCCESS);
    }

    frame_set = &(*tng_data_p)->current_trajectory_frame_set;

    if(frame_set->n_unwritten_frames > 0)
    {
        frame_set->n_frames = frame_set->n_unwritten_frames;
        tng_frame_set_write(*tng_data_p, TNG_USE_HASH);
    }

    return(tng_trajectory_destroy(tng_data_p));
}

tng_function_status DECLSPECDLLEXPORT tng_util_time_of_frame_get
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 double *time)
{
    int64_t first_frame;
    tng_trajectory_frame_set_t frame_set;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(time, "TNG library: time must not be a NULL pointer");

    stat = tng_frame_set_of_frame_find(tng_data, frame_nr);
    if(stat != TNG_SUCCESS)
    {
        fprintf(stderr, "TNG library: Cannot find frame nr %"PRId64". %s: %d\n",
               frame_nr, __FILE__, __LINE__);
        return(stat);
    }

    frame_set = &tng_data->current_trajectory_frame_set;
    first_frame = frame_set->first_frame;

    if(tng_data->time_per_frame <= 0)
    {
        return(TNG_FAILURE);
    }

    *time = frame_set->first_frame_time + (tng_data->time_per_frame * (frame_nr - first_frame));

    return(TNG_SUCCESS);
}

/*
tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_molecules_get
                (tng_trajectory_t tng_data,
                 int64_t *n_mols,
                 int64_t **molecule_cnt_list,
                 tng_molecule_t *mols)
{
    tng_trajectory_frame_set_t frame_set;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(n_mols, "TNG library: n_mols must not be a NULL pointer.");

    *n_mols = tng_data->n_molecules;

    frame_set = &tng_data->current_trajectory_frame_set;
    if(tng_data->var_num_atoms_flag && frame_set && frame_set->molecule_cnt_list)
    {
        *molecule_cnt_list = frame_set->molecule_cnt_list;
    }
    else
    {
        *molecule_cnt_list = tng_data->molecule_cnt_list;
    }

    *mols = tng_data->molecules;

    return(TNG_SUCCESS);
}
*/
/*
tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_molecule_add
                (tng_trajectory_t tng_data,
                 const char *name,
                 const int64_t cnt,
                 tng_molecule_t *mol)
{
    tng_function_status stat;

    TNG_ASSERT(name, "TNG library: name must not be a NULL pointer");
    TNG_ASSERT(cnt>=0, "TNG library: cnt must be >= 0");

    stat = tng_molecule_add(tng_data, name, mol);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }
    stat = tng_molecule_cnt_set(tng_data, *mol, cnt);

    return(stat);
}
*/
tng_function_status DECLSPECDLLEXPORT tng_util_molecule_particles_get
                (tng_trajectory_t tng_data,
                 const tng_molecule_t mol,
                 int64_t *n_particles,
                 char ***names,
                 char ***types,
                 char ***res_names,
                 int64_t **res_ids,
                 char ***chain_names,
                 int64_t **chain_ids)
{
    tng_atom_t atom;
    tng_residue_t res;
    tng_chain_t chain;
    int64_t i;
    (void)tng_data;

    *n_particles = mol->n_atoms;

    *names = malloc(sizeof(char *) * *n_particles);
    *types = malloc(sizeof(char *) * *n_particles);
    *res_names = malloc(sizeof(char *) * *n_particles);
    *chain_names = malloc(sizeof(char *) * *n_particles);
    *res_ids = malloc(sizeof(int64_t) * *n_particles);
    *chain_ids = malloc(sizeof(int64_t) * *n_particles);

    for(i = 0; i < *n_particles; i++)
    {
        atom = &mol->atoms[i];
        res = atom->residue;
        chain = res->chain;
        (*names)[i] = malloc(strlen(atom->name));
        strcpy(*names[i], atom->name);
        (*types)[i] = malloc(strlen(atom->atom_type));
        strcpy(*types[i], atom->atom_type);
        (*res_names)[i] = malloc(strlen(res->name));
        strcpy(*res_names[i], res->name);
        (*chain_names)[i] = malloc(strlen(chain->name));
        strcpy(*chain_names[i], chain->name);
        (*res_ids)[i] = res->id;
        (*chain_ids)[i] = chain->id;
    }

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_util_molecule_particles_set
                (tng_trajectory_t tng_data,
                 tng_molecule_t mol,
                 const int64_t n_particles,
                 const char **names,
                 const char **types,
                 const char **res_names,
                 const int64_t *res_ids,
                 const char **chain_names,
                 const int64_t *chain_ids)
{
    int64_t i;
    tng_chain_t chain;
    tng_residue_t residue;
    tng_atom_t atom;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(names, "TNG library: names must not be a NULL pointer");
    TNG_ASSERT(types, "TNG library: types must not be a NULL pointer");
    TNG_ASSERT(res_names, "TNG library: res_names must not be a NULL pointer");
    TNG_ASSERT(res_ids, "TNG library: res_ids must not be a NULL pointer");
    TNG_ASSERT(chain_names, "TNG library: chain_names must not be a NULL pointer");
    TNG_ASSERT(chain_ids, "TNG library: chain_ids must not be a NULL pointer");

    for(i = 0; i < n_particles; i++)
    {
        if(tng_molecule_chain_find(tng_data, mol, chain_names[i], chain_ids[i],
           &chain) == TNG_FAILURE)
        {
            stat = tng_molecule_chain_add(tng_data, mol, chain_names[i],
                                          &chain);
            if(stat != TNG_SUCCESS)
            {
                return(stat);
            }
        }
        if(tng_chain_residue_find(tng_data, chain, res_names[i], res_ids[i],
           &residue) == TNG_FAILURE)
        {
            stat = tng_chain_residue_add(tng_data, chain, res_names[i],
                                         &residue);
            if(stat != TNG_SUCCESS)
            {
                return(stat);
            }
        }
        stat = tng_residue_atom_add(tng_data, residue, names[i], types[i], &atom);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_util_pos_read
                (tng_trajectory_t tng_data,
                 float **positions, int64_t *stride_length)
{
    int64_t n_frames, n_particles, n_values_per_frame;
    char type;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(positions, "TNG library: positions must not be a NULL pointer");
    TNG_ASSERT(stride_length, "TNG library: stride_length must not be a NULL pointer");

    stat = tng_num_frames_get(tng_data, &n_frames);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    stat = tng_particle_data_vector_interval_get(tng_data, TNG_TRAJ_POSITIONS,
                                                 0, n_frames - 1, TNG_USE_HASH,
                                                 (void **)positions,
                                                 &n_particles,
                                                 stride_length,
                                                 &n_values_per_frame,
                                                 &type);

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_read
                (tng_trajectory_t tng_data,
                 float **velocities, int64_t *stride_length)
{
    int64_t n_frames, n_particles, n_values_per_frame;
    char type;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(velocities, "TNG library: velocities must not be a NULL pointer");
    TNG_ASSERT(stride_length, "TNG library: stride_length must not be a NULL pointer");

    stat = tng_num_frames_get(tng_data, &n_frames);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    stat = tng_particle_data_vector_interval_get(tng_data, TNG_TRAJ_VELOCITIES,
                                                 0, n_frames - 1, TNG_USE_HASH,
                                                 (void **)velocities,
                                                 &n_particles,
                                                 stride_length,
                                                 &n_values_per_frame,
                                                 &type);

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_read
                (tng_trajectory_t tng_data,
                 float **forces, int64_t *stride_length)
{
    int64_t n_frames, n_particles, n_values_per_frame;
    char type;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(forces, "TNG library: forces must not be a NULL pointer");
    TNG_ASSERT(stride_length, "TNG library: stride_length must not be a NULL pointer");

    stat = tng_num_frames_get(tng_data, &n_frames);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    stat = tng_particle_data_vector_interval_get(tng_data, TNG_TRAJ_FORCES,
                                                 0, n_frames - 1, TNG_USE_HASH,
                                                 (void **)forces,
                                                 &n_particles,
                                                 stride_length,
                                                 &n_values_per_frame,
                                                 &type);

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_read
                (tng_trajectory_t tng_data,
                 float **box_shape,
                 int64_t *stride_length)
{
    int64_t n_frames, n_values_per_frame;
    char type;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(box_shape, "TNG library: box_shape must not be a NULL pointer");
    TNG_ASSERT(stride_length, "TNG library: stride_length must not be a NULL pointer");

    stat = tng_num_frames_get(tng_data, &n_frames);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    stat = tng_data_vector_interval_get(tng_data, TNG_TRAJ_BOX_SHAPE,
                                        0, n_frames - 1, TNG_USE_HASH,
                                        (void **)box_shape,
                                        stride_length,
                                        &n_values_per_frame,
                                        &type);

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_util_particle_data_next_frame_read
                (tng_trajectory_t tng_data,
                 const int64_t block_id,
                 void **values,
                 char *data_type,
                 int64_t *retrieved_frame_number,
                 double *retrieved_time)
{
    tng_trajectory_frame_set_t frame_set;
    tng_particle_data_t data = 0;
    tng_function_status stat;
    int size;
    int64_t i, data_size, n_particles;
    void *temp;
    long file_pos;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(values, "TNG library: The pointer to the values array must not be a NULL pointer");
    TNG_ASSERT(data_type, "TNG library: The pointer to the data type of the returned data must not be a NULL pointer");
    TNG_ASSERT(retrieved_frame_number, "TNG library: The pointer to the frame number of the returned data must not be a NULL pointer");
    TNG_ASSERT(retrieved_time, "TNG library: The pointer to the time of the returned data must not be a NULL pointer");

    frame_set = &tng_data->current_trajectory_frame_set;

    stat = tng_particle_data_find(tng_data, block_id, &data);
    if(stat != TNG_SUCCESS)
    {
        stat = tng_frame_set_read_current_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
        file_pos = ftell(tng_data->input_file);
        while(stat != TNG_SUCCESS && file_pos < tng_data->input_file_len)
        {
            stat = tng_frame_set_read_next_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
            file_pos = ftell(tng_data->input_file);
        }
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
        stat = tng_particle_data_find(tng_data, block_id, &data);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
    }
    if(data->last_retrieved_frame < 0)
    {
        fseek(tng_data->input_file,
              (long)tng_data->first_trajectory_frame_set_input_file_pos,
              SEEK_SET);
        stat = tng_frame_set_read(tng_data, TNG_USE_HASH);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
        stat = tng_frame_set_read_current_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }

        i = data->first_frame_with_data;
    }
    else
    {
        if(data->n_frames == 1)
        {
            i = data->last_retrieved_frame + 1;
        }
        else
        {
            i = data->last_retrieved_frame + data->stride_length;
        }
        if(i < frame_set->first_frame || i >= frame_set->first_frame + frame_set->n_frames)
        {
            stat = tng_frame_set_of_frame_find(tng_data, i);
            if(stat != TNG_SUCCESS)
            {
                /* If the frame set search found the frame set after the starting
                 * frame set there is a gap in the frame sets. So, even if the frame
                 * was not found the next frame with data is still in the found
                 * frame set. */
                if(stat == TNG_CRITICAL)
                {
                    return(stat);
                }
                i = frame_set->first_frame;
            }
        }
        if(data->last_retrieved_frame < frame_set->first_frame)
        {
            stat = tng_frame_set_read_current_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
            if(stat != TNG_SUCCESS)
            {
                return(stat);
            }
        }
    }
    data->last_retrieved_frame = i;
    *retrieved_frame_number = i;
    if(frame_set->first_frame_time >= 0 && tng_data->time_per_frame >= 0)
    {
        *retrieved_time = frame_set->first_frame_time +
                        (i - frame_set->first_frame) *
                        tng_data->time_per_frame;
    }
    else
    {
        *retrieved_time = 0;
    }

    if(data->stride_length > 1)
    {
        i = (i - data->first_frame_with_data) / data->stride_length;
    }
    else
    {
        i = (i - frame_set->first_frame);
    }

    tng_num_particles_get(tng_data, &n_particles);

    *data_type = data->datatype;

    switch(*data_type)
    {
    case TNG_CHAR_DATA:
        return(TNG_FAILURE);
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

    data_size = size * n_particles * data->n_values_per_frame;

//     fprintf(stderr, "TNG library: TEMP: i = %"PRId64", data_size = %"PRId64", size = %d, n_particles = %"PRId64", n_values_per_frame = %"PRId64"\n",
//            i, data_size, size, n_particles, data->n_values_per_frame);

    temp = realloc(*values, data_size);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               data_size, __FILE__, __LINE__);
        free(*values);
        *values = 0;
        return(TNG_CRITICAL);
    }

    *values = temp;

    memcpy(*values, (char *)data->values + i * data_size, data_size);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_util_non_particle_data_next_frame_read
                (tng_trajectory_t tng_data,
                 const int64_t block_id,
                 void **values,
                 char *data_type,
                 int64_t *retrieved_frame_number,
                 double *retrieved_time)
{
    tng_trajectory_frame_set_t frame_set;
    tng_non_particle_data_t data = 0;
    tng_function_status stat;
    int size;
    int64_t i, data_size;
    void *temp;
    long file_pos;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(values, "TNG library: The pointer to the values array must not be a NULL pointer");
    TNG_ASSERT(data_type, "TNG library: The pointer to the data type of the returned data must not be a NULL pointer");
    TNG_ASSERT(retrieved_frame_number, "TNG library: The pointer to the frame number of the returned data must not be a NULL pointer");
    TNG_ASSERT(retrieved_time, "TNG library: The pointer to the time of the returned data must not be a NULL pointer");

    frame_set = &tng_data->current_trajectory_frame_set;

    stat = tng_data_find(tng_data, block_id, &data);
    if(stat != TNG_SUCCESS)
    {
        stat = tng_frame_set_read_current_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
        file_pos = ftell(tng_data->input_file);
        while(stat != TNG_SUCCESS && file_pos < tng_data->input_file_len)
        {
            stat = tng_frame_set_read_next_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
            file_pos = ftell(tng_data->input_file);
        }
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
        stat = tng_data_find(tng_data, block_id, &data);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
    }
    if(data->last_retrieved_frame < 0)
    {
        fseek(tng_data->input_file,
                (long)tng_data->first_trajectory_frame_set_input_file_pos,
                SEEK_SET);
        stat = tng_frame_set_read(tng_data, TNG_USE_HASH);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
        stat = tng_frame_set_read_current_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }

        i = data->first_frame_with_data;
    }
    else
    {
        if(data->n_frames == 1)
        {
            i = data->last_retrieved_frame + 1;
        }
        else
        {
            i = data->last_retrieved_frame + data->stride_length;
        }
        if(i < frame_set->first_frame || i >= frame_set->first_frame + frame_set->n_frames)
        {
            stat = tng_frame_set_of_frame_find(tng_data, i);
            if(stat != TNG_SUCCESS)
            {
                /* If the frame set search found the frame set after the starting
                 * frame set there is a gap in the frame sets. So, even if the frame
                 * was not found the next frame with data is still in the found
                 * frame set. */
                if(stat == TNG_CRITICAL)
                {
                    return(stat);
                }
                i = frame_set->first_frame;
            }
        }
        if(data->last_retrieved_frame < frame_set->first_frame)
        {
            stat = tng_frame_set_read_current_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
            if(stat != TNG_SUCCESS)
            {
                return(stat);
            }
        }
    }
    data->last_retrieved_frame = i;
    *retrieved_frame_number = i;
    if(frame_set->first_frame_time >= 0 && tng_data->time_per_frame >= 0)
    {
        *retrieved_time = frame_set->first_frame_time +
                        (i - frame_set->first_frame) *
                        tng_data->time_per_frame;
    }
    else
    {
        *retrieved_time = 0;
    }

    if(data->stride_length > 1)
    {
        i = (i - data->first_frame_with_data) / data->stride_length;
    }
    else
    {
        i = (i - frame_set->first_frame);
    }

    *data_type = data->datatype;

    switch(*data_type)
    {
    case TNG_CHAR_DATA:
        return(TNG_FAILURE);
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

    data_size = size * data->n_values_per_frame;

    temp = realloc(*values, data_size);
    if(!temp)
    {
        fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               data_size, __FILE__, __LINE__);
        free(*values);
        *values = 0;
        return(TNG_CRITICAL);
    }

    *values = temp;

    memcpy(*values, (char *)data->values + i * data_size, data_size);

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_util_pos_read_range
                (tng_trajectory_t tng_data,
                 const int64_t first_frame,
                 const int64_t last_frame,
                 float **positions,
                 int64_t *stride_length)
{
    int64_t n_particles, n_values_per_frame;
    char type;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(positions, "TNG library: positions must not be a NULL pointer");
    TNG_ASSERT(first_frame <= last_frame, "TNG library: first_frame must be lower or equal to last_frame.");
    TNG_ASSERT(stride_length, "TNG library: stride_length must not be a NULL pointer");

    stat = tng_particle_data_vector_interval_get(tng_data, TNG_TRAJ_POSITIONS,
                                                 first_frame, last_frame,
                                                 TNG_USE_HASH,
                                                 (void **)positions,
                                                 &n_particles,
                                                 stride_length,
                                                 &n_values_per_frame,
                                                 &type);

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_read_range
                (tng_trajectory_t tng_data,
                 const int64_t first_frame,
                 const int64_t last_frame,
                 float **velocities,
                 int64_t *stride_length)
{
    int64_t n_particles, n_values_per_frame;
    char type;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(velocities, "TNG library: velocities must not be a NULL pointer");
    TNG_ASSERT(first_frame <= last_frame, "TNG library: first_frame must be lower or equal to last_frame.");
    TNG_ASSERT(stride_length, "TNG library: stride_length must not be a NULL pointer");

    stat = tng_particle_data_vector_interval_get(tng_data, TNG_TRAJ_VELOCITIES,
                                                 first_frame, last_frame,
                                                 TNG_USE_HASH,
                                                 (void **)velocities,
                                                 &n_particles,
                                                 stride_length,
                                                 &n_values_per_frame,
                                                 &type);

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_read_range
                (tng_trajectory_t tng_data,
                 const int64_t first_frame,
                 const int64_t last_frame,
                 float **forces,
                 int64_t *stride_length)
{
    int64_t n_particles, n_values_per_frame;
    char type;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(forces, "TNG library: forces must not be a NULL pointer");
    TNG_ASSERT(first_frame <= last_frame, "TNG library: first_frame must be lower or equal to last_frame.");
    TNG_ASSERT(stride_length, "TNG library: stride_length must not be a NULL pointer");

    stat = tng_particle_data_vector_interval_get(tng_data, TNG_TRAJ_FORCES,
                                                 first_frame, last_frame,
                                                 TNG_USE_HASH,
                                                 (void **)forces,
                                                 &n_particles,
                                                 stride_length,
                                                 &n_values_per_frame,
                                                 &type);

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_read_range
                (tng_trajectory_t tng_data,
                 const int64_t first_frame,
                 const int64_t last_frame,
                 float **box_shape,
                 int64_t *stride_length)
{
    int64_t n_values_per_frame;
    char type;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(box_shape, "TNG library: box_shape must not be a NULL pointer");
    TNG_ASSERT(first_frame <= last_frame, "TNG library: first_frame must be lower or equal to last_frame.");
    TNG_ASSERT(stride_length, "TNG library: stride_length must not be a NULL pointer");

    stat = tng_data_vector_interval_get(tng_data, TNG_TRAJ_BOX_SHAPE,
                                        first_frame, last_frame,
                                        TNG_USE_HASH,
                                        (void **)box_shape,
                                        stride_length,
                                        &n_values_per_frame,
                                        &type);

    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_util_generic_write_interval_set
                (tng_trajectory_t tng_data,
                 const int64_t i,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression)
{
    tng_trajectory_frame_set_t frame_set;
    tng_particle_data_t p_data;
    tng_non_particle_data_t np_data;
    int64_t n_particles, n_frames;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(i >= 0, "TNG library: i (writing interval) must be >= 0.");

    if(i <= 0)
    {
        fprintf(stderr, "TNG library: Cannot set writing frequency to %"PRId64". %s: %d\n",
               i, __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    frame_set = &tng_data->current_trajectory_frame_set;

    if(!frame_set || tng_data->n_trajectory_frame_sets <= 0)
    {
        n_frames = tng_data->frame_set_n_frames;

        stat = tng_frame_set_new(tng_data, 0, n_frames);
        if(stat != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot create frame set.  %s: %d\n", __FILE__,
                __LINE__);
            return(stat);
        }
    }
    else
    {
        n_frames = frame_set->n_frames;
    }

    if(particle_dependency == TNG_PARTICLE_BLOCK_DATA)
    {
        tng_num_particles_get(tng_data, &n_particles);
        if(n_particles <= 0)
        {
            return(TNG_FAILURE);
        }

        if(tng_particle_data_find(tng_data, block_id, &p_data)
        != TNG_SUCCESS)
        {
            stat = tng_particle_data_block_add(tng_data, block_id,
                                               block_name,
                                               TNG_FLOAT_DATA,
                                               TNG_TRAJECTORY_BLOCK,
                                               n_frames, n_values_per_frame, i,
                                               0, n_particles,
                                               compression, 0);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error %s adding data block. %s: %d\n", block_name,
                       __FILE__, __LINE__);
                return(stat);
            }
            p_data = &frame_set->tr_particle_data[frame_set->
                                                  n_particle_data_blocks - 1];
            stat = tng_allocate_particle_data_mem(tng_data, p_data, n_frames,
                                                  i, n_particles,
                                                  n_values_per_frame);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error allocating particle data memory. %s: %d\n",
                       __FILE__, __LINE__);
                return(stat);
            }
        }
        else
        {
            if(p_data->stride_length != i)
            {
                p_data->stride_length = i;
                stat = tng_allocate_particle_data_mem(tng_data, p_data, n_frames,
                                                      i, n_particles,
                                                      n_values_per_frame);
                if(stat != TNG_SUCCESS)
                {
                    fprintf(stderr, "TNG library: Error allocating particle data memory. %s: %d\n",
                           __FILE__, __LINE__);
                    return(stat);
                }
            }
        }
    }
    else
    {
        if(tng_data_find(tng_data, block_id, &np_data) != TNG_SUCCESS)
        {
            stat = tng_data_block_add(tng_data, block_id, block_name,
                                      TNG_FLOAT_DATA, TNG_TRAJECTORY_BLOCK,
                                      n_frames, n_values_per_frame,
                                      i, compression, 0);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error %s adding data block. %s: %d\n", block_name,
                       __FILE__, __LINE__);
                return(stat);
            }
            np_data = &frame_set->tr_data[frame_set->
                                          n_data_blocks - 1];
            stat = tng_allocate_data_mem(tng_data, np_data, n_frames,
                                         i, n_values_per_frame);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error allocating particle data memory. %s: %d\n",
                       __FILE__, __LINE__);
                return(stat);
            }
        }
        else
        {
            if(np_data->stride_length != i)
            {
                np_data->stride_length = i;
                stat = tng_allocate_data_mem(tng_data, np_data, n_frames,
                                             i, n_values_per_frame);
                if(stat != TNG_SUCCESS)
                {
                    fprintf(stderr, "TNG library: Error allocating particle data memory. %s: %d\n",
                           __FILE__, __LINE__);
                    return(stat);
                }
            }
        }
    }

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_util_generic_write_interval_double_set
                (tng_trajectory_t tng_data,
                 const int64_t i,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression)
{
    tng_trajectory_frame_set_t frame_set;
    tng_particle_data_t p_data;
    tng_non_particle_data_t np_data;
    int64_t n_particles, n_frames;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(i >= 0, "TNG library: i (writing interval) must be >= 0.");

    if(i <= 0)
    {
        fprintf(stderr, "TNG library: Cannot set writing frequency to %"PRId64". %s: %d\n",
               i, __FILE__, __LINE__);
        return(TNG_FAILURE);
    }

    frame_set = &tng_data->current_trajectory_frame_set;

    if(!frame_set || tng_data->n_trajectory_frame_sets <= 0)
    {
        n_frames = tng_data->frame_set_n_frames;

        stat = tng_frame_set_new(tng_data, 0, n_frames);
        if(stat != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot create frame set.  %s: %d\n", __FILE__,
                __LINE__);
            return(stat);
        }
    }
    else
    {
        n_frames = frame_set->n_frames;
    }

    if(particle_dependency == TNG_PARTICLE_BLOCK_DATA)
    {
        tng_num_particles_get(tng_data, &n_particles);

        if(n_particles <= 0)
        {
            return(TNG_FAILURE);
        }

        if(tng_particle_data_find(tng_data, block_id, &p_data)
        != TNG_SUCCESS)
        {
            stat = tng_particle_data_block_add(tng_data, block_id,
                                            block_name,
                                            TNG_DOUBLE_DATA,
                                            TNG_TRAJECTORY_BLOCK,
                                            n_frames, n_values_per_frame, i,
                                            0, n_particles,
                                            compression, 0);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error %s adding data block. %s: %d\n", block_name,
                       __FILE__, __LINE__);
                return(stat);
            }
            p_data = &frame_set->tr_particle_data[frame_set->
                                                  n_particle_data_blocks - 1];
            stat = tng_allocate_particle_data_mem(tng_data, p_data, n_frames,
                                                  i, n_particles,
                                                  n_values_per_frame);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error allocating particle data memory. %s: %d\n",
                       __FILE__, __LINE__);
                return(stat);
            }
        }
        else
        {
            p_data->stride_length = i;
        }
    }
    else
    {
        if(tng_data_find(tng_data, block_id, &np_data) != TNG_SUCCESS)
        {
            stat = tng_data_block_add(tng_data, block_id, block_name,
                                      TNG_DOUBLE_DATA, TNG_TRAJECTORY_BLOCK,
                                      n_frames, n_values_per_frame,
                                      i, compression, 0);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error %s adding data block. %s: %d\n", block_name,
                       __FILE__, __LINE__);
                return(stat);
            }
            np_data = &frame_set->tr_data[frame_set->
                                          n_data_blocks - 1];
            stat = tng_allocate_data_mem(tng_data, np_data, n_frames,
                                         i, n_values_per_frame);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error allocating particle data memory. %s: %d\n",
                       __FILE__, __LINE__);
                return(stat);
            }
        }
        else
        {
            np_data->stride_length = i;
        }
    }

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_util_generic_write_frequency_set
                (tng_trajectory_t tng_data,
                 const int64_t i,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression)
{
    fprintf(stderr, "TNG library: Using obsolete function tng_util_generic_write_frequency_set(). "
           "See documentation. %s: %d", __FILE__, __LINE__);
    return(tng_util_generic_write_interval_set(tng_data, i, n_values_per_frame,
                                               block_id, block_name,
                                               particle_dependency,
                                               compression));
}
tng_function_status DECLSPECDLLEXPORT tng_util_pos_write_interval_set
                (tng_trajectory_t tng_data,
                 const int64_t i)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(i > 0, "TNG library: i (writing interval) must be >= 0.");

    return(tng_util_generic_write_interval_set(tng_data, i, 3,
                                               TNG_TRAJ_POSITIONS,
                                               "POSITIONS",
                                               TNG_PARTICLE_BLOCK_DATA,
                                               TNG_TNG_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_pos_write_interval_double_set
                (tng_trajectory_t tng_data,
                 const int64_t i)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(i > 0, "TNG library: i (writing interval) must be >= 0.");

    return(tng_util_generic_write_interval_double_set(tng_data, i, 3,
                                                      TNG_TRAJ_POSITIONS,
                                                      "POSITIONS",
                                                      TNG_PARTICLE_BLOCK_DATA,
                                                      TNG_TNG_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_pos_write_frequency_set
                (tng_trajectory_t tng_data,
                 const int64_t i)
{
    fprintf(stderr, "TNG library: Using obsolete function tng_util_pos_write_frequency_set(). "
           "See documentation. %s: %d", __FILE__, __LINE__);
    return(tng_util_pos_write_interval_set(tng_data, i));
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_write_interval_set
                (tng_trajectory_t tng_data,
                 const int64_t i)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(i > 0, "TNG library: i (writing interval) must be >= 0.");

    return(tng_util_generic_write_interval_set(tng_data, i, 3,
                                               TNG_TRAJ_VELOCITIES,
                                               "VELOCITIES",
                                               TNG_PARTICLE_BLOCK_DATA,
                                               TNG_TNG_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_write_interval_double_set
                (tng_trajectory_t tng_data,
                 const int64_t i)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(i > 0, "TNG library: i (writing interval) must be >= 0.");

    return(tng_util_generic_write_interval_double_set(tng_data, i, 3,
                                                      TNG_TRAJ_VELOCITIES,
                                                      "VELOCITIES",
                                                      TNG_PARTICLE_BLOCK_DATA,
                                                      TNG_TNG_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_write_frequency_set
                (tng_trajectory_t tng_data,
                 const int64_t i)
{
    fprintf(stderr, "TNG library: Using obsolete function tng_util_vel_write_frequency_set(). "
           "See documentation. %s: %d", __FILE__, __LINE__);
    return(tng_util_vel_write_interval_set(tng_data, i));
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_write_interval_set
                (tng_trajectory_t tng_data,
                 const int64_t i)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(i > 0, "TNG library: i (writing interval) must be >= 0.");

    return(tng_util_generic_write_interval_set(tng_data, i, 3,
                                               TNG_TRAJ_FORCES,
                                               "FORCES",
                                               TNG_PARTICLE_BLOCK_DATA,
                                               TNG_GZIP_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_write_interval_double_set
                (tng_trajectory_t tng_data,
                 const int64_t i)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(i > 0, "TNG library: i (writing interval) must be >= 0.");

    return(tng_util_generic_write_interval_double_set(tng_data, i, 3,
                                                      TNG_TRAJ_FORCES,
                                                      "FORCES",
                                                      TNG_PARTICLE_BLOCK_DATA,
                                                      TNG_GZIP_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_write_frequency_set
                (tng_trajectory_t tng_data,
                 const int64_t i)
{
    fprintf(stderr, "TNG library: Using obsolete function tng_util_force_write_frequency_set(). "
           "See documentation. %s: %d", __FILE__, __LINE__);
    return(tng_util_force_write_interval_set(tng_data, i));
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_write_interval_set
                (tng_trajectory_t tng_data,
                 const int64_t i)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(i > 0, "TNG library: i (writing interval) must be >= 0.");

    return(tng_util_generic_write_interval_set(tng_data, i, 9,
                                               TNG_TRAJ_BOX_SHAPE,
                                               "BOX SHAPE",
                                               TNG_NON_PARTICLE_BLOCK_DATA,
                                               TNG_GZIP_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_write_interval_double_set
                (tng_trajectory_t tng_data,
                 const int64_t i)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(i > 0, "TNG library: i (writing interval) must be >= 0.");

    return(tng_util_generic_write_interval_double_set(tng_data, i, 9,
                                                      TNG_TRAJ_BOX_SHAPE,
                                                      "BOX SHAPE",
                                                      TNG_NON_PARTICLE_BLOCK_DATA,
                                                      TNG_GZIP_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_write_frequency_set
                (tng_trajectory_t tng_data,
                 const int64_t i)
{
    fprintf(stderr, "TNG library: Using obsolete function tng_util_box_shape_write_frequency_set(). "
           "See documentation. %s: %d", __FILE__, __LINE__);
    return(tng_util_box_shape_write_interval_set(tng_data, i));
}

tng_function_status DECLSPECDLLEXPORT tng_util_generic_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const float *values,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression)
{
    tng_trajectory_frame_set_t frame_set;
    tng_particle_data_t p_data;
    tng_non_particle_data_t np_data;
    int64_t n_particles = 0, n_frames, stride_length = 100, frame_pos;
    int64_t last_frame;
    int is_first_frame_flag = 0;
    char block_type_flag;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(values, "TNG library: values must not be a NULL pointer");

    if(particle_dependency == TNG_PARTICLE_BLOCK_DATA)
    {
        tng_num_particles_get(tng_data, &n_particles);
        TNG_ASSERT(n_particles > 0, "TNG library: There must be particles in the system to write particle data.");
    }

    if(values == 0)
    {
        return(TNG_FAILURE);
    }

    frame_set = &tng_data->current_trajectory_frame_set;

    if(frame_nr < 0)
    {
        block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
        n_frames = stride_length = 1;
    }
    else
    {
        block_type_flag = TNG_TRAJECTORY_BLOCK;

        if(!frame_set || tng_data->n_trajectory_frame_sets <= 0)
        {
            stat = tng_frame_set_new(tng_data, 0, tng_data->frame_set_n_frames);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot create frame set.  %s: %d\n", __FILE__,
                    __LINE__);
                return(stat);
            }
        }
        last_frame = frame_set->first_frame +
                     frame_set->n_frames - 1;
        if(frame_nr > last_frame)
        {
            stat = tng_frame_set_write(tng_data, TNG_USE_HASH);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot write frame set.  %s: %d\n", __FILE__,
                    __LINE__);
                return(stat);
            }
            if(last_frame + tng_data->frame_set_n_frames < frame_nr)
            {
                last_frame = frame_nr - 1;
            }
            stat = tng_frame_set_new(tng_data, last_frame + 1,
                                     tng_data->frame_set_n_frames);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot create frame set.  %s: %d\n", __FILE__,
                    __LINE__);
                return(stat);
            }
        }
        if(frame_set->n_unwritten_frames == 0)
        {
            is_first_frame_flag = 1;
        }
        frame_set->n_unwritten_frames = frame_nr -
                                        frame_set->first_frame + 1;

        n_frames = frame_set->n_frames;
    }

    if(particle_dependency == TNG_PARTICLE_BLOCK_DATA)
    {
        if(tng_particle_data_find(tng_data, block_id, &p_data)
        != TNG_SUCCESS)
        {
            stat = tng_particle_data_block_add(tng_data, block_id,
                                               block_name,
                                               TNG_FLOAT_DATA,
                                               block_type_flag,
                                               n_frames, n_values_per_frame,
                                               stride_length,
                                               0, n_particles,
                                               compression, 0);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error %s adding data block. %s: %d\n", block_name,
                       __FILE__, __LINE__);
                return(stat);
            }
            if(block_type_flag == TNG_TRAJECTORY_BLOCK)
            {
                p_data = &frame_set->tr_particle_data[frame_set->
                                                    n_particle_data_blocks - 1];
            }
            else
            {
                p_data = &tng_data->non_tr_particle_data[tng_data->
                                                    n_particle_data_blocks - 1];
            }
            stat = tng_allocate_particle_data_mem(tng_data, p_data, n_frames,
                                                  stride_length, n_particles,
                                                  n_values_per_frame);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error allocating particle data memory. %s: %d\n",
                       __FILE__, __LINE__);
                return(stat);
            }
        }

        if(block_type_flag == TNG_TRAJECTORY_BLOCK)
        {
            stride_length = p_data->stride_length;

            if(is_first_frame_flag || p_data->first_frame_with_data < frame_set->first_frame)
            {
                p_data->first_frame_with_data = frame_nr;
                frame_pos = 0;
            }
            else
            {
                frame_pos = (frame_nr - frame_set->first_frame) / stride_length;
            }

            memcpy((char *)p_data->values + sizeof(float) * frame_pos * n_particles *
                   n_values_per_frame, values, sizeof(float) *
                   n_particles * n_values_per_frame);
        }
        else
        {
            memcpy(p_data->values, values, sizeof(float) * n_particles *
                   n_values_per_frame);
        }
    }
    else
    {
        if(tng_data_find(tng_data, block_id, &np_data) != TNG_SUCCESS)
        {
            stat = tng_data_block_add(tng_data, block_id, block_name,
                                      TNG_FLOAT_DATA, block_type_flag,
                                      n_frames, n_values_per_frame,
                                      stride_length, compression, 0);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error %s adding data block. %s: %d\n", block_name,
                       __FILE__, __LINE__);
                return(stat);
            }
            if(block_type_flag == TNG_TRAJECTORY_BLOCK)
            {
                np_data = &frame_set->tr_data[frame_set->
                                              n_data_blocks - 1];
            }
            else
            {
                np_data = &tng_data->non_tr_data[tng_data->
                                                 n_data_blocks - 1];
            }
            stat = tng_allocate_data_mem(tng_data, np_data, n_frames,
                                         stride_length, n_values_per_frame);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error allocating particle data memory. %s: %d\n",
                       __FILE__, __LINE__);
                return(stat);
            }
        }

        if(block_type_flag == TNG_TRAJECTORY_BLOCK)
        {
            stride_length = np_data->stride_length;

            if(is_first_frame_flag || np_data->first_frame_with_data < frame_set->first_frame)
            {
                np_data->first_frame_with_data = frame_nr;
                frame_pos = 0;
            }
            else
            {
                frame_pos = (frame_nr - frame_set->first_frame) / stride_length;
            }

            memcpy((char *)np_data->values + sizeof(float) * frame_pos *
                   n_values_per_frame, values, sizeof(float) *
                   n_values_per_frame);
        }
        else
        {
            memcpy(np_data->values, values, sizeof(float) * n_values_per_frame);
        }
    }

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_util_generic_double_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double *values,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression)
{
    tng_trajectory_frame_set_t frame_set;
    tng_particle_data_t p_data;
    tng_non_particle_data_t np_data;
    int64_t n_particles = 0, n_frames, stride_length = 100, frame_pos;
    int64_t last_frame;
    int is_first_frame_flag = 0;
    char block_type_flag;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(values, "TNG library: values must not be a NULL pointer");

    if(particle_dependency == TNG_PARTICLE_BLOCK_DATA)
    {
        tng_num_particles_get(tng_data, &n_particles);
        TNG_ASSERT(n_particles > 0, "TNG library: There must be particles in the system to write particle data.");
    }

    if(values == 0)
    {
        return(TNG_FAILURE);
    }

    frame_set = &tng_data->current_trajectory_frame_set;

    if(frame_nr < 0)
    {
        block_type_flag = TNG_NON_TRAJECTORY_BLOCK;
        n_frames = stride_length = 1;
    }
    else
    {
        block_type_flag = TNG_TRAJECTORY_BLOCK;

        n_frames = tng_data->frame_set_n_frames;

        if(!frame_set || tng_data->n_trajectory_frame_sets <= 0)
        {
            stat = tng_frame_set_new(tng_data, 0, n_frames);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot create frame set.  %s: %d\n", __FILE__,
                    __LINE__);
                return(stat);
            }
        }
        else
        {
            n_frames = frame_set->n_frames;
        }
        last_frame = frame_set->first_frame +
                     frame_set->n_frames - 1;
        if(frame_nr > last_frame)
        {
            stat = tng_frame_set_write(tng_data, TNG_USE_HASH);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot write frame set.  %s: %d\n", __FILE__,
                    __LINE__);
                return(stat);
            }
            if(last_frame + tng_data->frame_set_n_frames < frame_nr)
            {
                last_frame = frame_nr - 1;
            }
            stat = tng_frame_set_new(tng_data, last_frame + 1, n_frames);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Cannot create frame set.  %s: %d\n", __FILE__,
                    __LINE__);
                return(stat);
            }
        }
        if(frame_set->n_unwritten_frames == 0)
        {
            is_first_frame_flag = 1;
        }
        frame_set->n_unwritten_frames = frame_nr -
                                        frame_set->first_frame + 1;
    }

    if(particle_dependency == TNG_PARTICLE_BLOCK_DATA)
    {
        if(tng_particle_data_find(tng_data, block_id, &p_data)
        != TNG_SUCCESS)
        {
            stat = tng_particle_data_block_add(tng_data, block_id,
                                            block_name,
                                            TNG_DOUBLE_DATA,
                                            block_type_flag,
                                            n_frames, n_values_per_frame,
                                            stride_length,
                                            0, n_particles,
                                            compression, 0);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error %s adding data block. %s: %d\n", block_name,
                       __FILE__, __LINE__);
                return(stat);
            }
            if(block_type_flag == TNG_TRAJECTORY_BLOCK)
            {
                p_data = &frame_set->tr_particle_data[frame_set->
                                                    n_particle_data_blocks - 1];
            }
            else
            {
                p_data = &tng_data->non_tr_particle_data[tng_data->
                                                    n_particle_data_blocks - 1];
            }
            stat = tng_allocate_particle_data_mem(tng_data, p_data, n_frames,
                                                  stride_length, n_particles,
                                                  n_values_per_frame);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error allocating particle data memory. %s: %d\n",
                       __FILE__, __LINE__);
                return(stat);
            }
        }

        if(block_type_flag == TNG_TRAJECTORY_BLOCK)
        {
            stride_length = p_data->stride_length;

            if(is_first_frame_flag)
            {
                p_data->first_frame_with_data = frame_nr;
                frame_pos = 0;
            }
            else
            {
                frame_pos = (frame_nr - frame_set->first_frame) / stride_length;
            }

            memcpy((char *)p_data->values + sizeof(double) * frame_pos * n_particles *
                   n_values_per_frame, values, sizeof(double) *
                   n_particles * n_values_per_frame);
        }
        else
        {
            memcpy(p_data->values, values, sizeof(double) * n_particles *
                   n_values_per_frame);
        }
    }
    else
    {
        if(tng_data_find(tng_data, block_id, &np_data) != TNG_SUCCESS)
        {
            stat = tng_data_block_add(tng_data, block_id, block_name,
                                      TNG_DOUBLE_DATA, block_type_flag,
                                      n_frames, n_values_per_frame,
                                      stride_length, compression, 0);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error %s adding data block. %s: %d\n", block_name,
                       __FILE__, __LINE__);
                return(stat);
            }
            if(block_type_flag == TNG_TRAJECTORY_BLOCK)
            {
                np_data = &frame_set->tr_data[frame_set->
                                              n_data_blocks - 1];
            }
            else
            {
                np_data = &tng_data->non_tr_data[tng_data->
                                                 n_data_blocks - 1];
            }
            stat = tng_allocate_data_mem(tng_data, np_data, n_frames,
                                         stride_length, n_values_per_frame);
            if(stat != TNG_SUCCESS)
            {
                fprintf(stderr, "TNG library: Error allocating particle data memory. %s: %d\n",
                       __FILE__, __LINE__);
                return(stat);
            }
        }

        if(block_type_flag == TNG_TRAJECTORY_BLOCK)
        {
            stride_length = np_data->stride_length;

            if(is_first_frame_flag)
            {
                np_data->first_frame_with_data = frame_nr;
                frame_pos = 0;
            }
            else
            {
                frame_pos = (frame_nr - frame_set->first_frame) / stride_length;
            }

            memcpy((char *)np_data->values + sizeof(double) * frame_pos *
                   n_values_per_frame, values, sizeof(double) *
                   n_values_per_frame);
        }
        else
        {
            memcpy(np_data->values, values, sizeof(double) * n_values_per_frame);
        }
    }

    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_util_pos_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const float *positions)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(positions, "TNG library: positions must not be a NULL pointer");

    return(tng_util_generic_write(tng_data, frame_nr, positions, 3,
                                  TNG_TRAJ_POSITIONS, "POSITIONS",
                                  TNG_PARTICLE_BLOCK_DATA,
                                  TNG_TNG_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_pos_double_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double *positions)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(positions, "TNG library: positions must not be a NULL pointer");

    return(tng_util_generic_double_write(tng_data, frame_nr, positions, 3,
                                         TNG_TRAJ_POSITIONS, "POSITIONS",
                                         TNG_PARTICLE_BLOCK_DATA,
                                         TNG_TNG_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const float *velocities)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(velocities, "TNG library: velocities must not be a NULL pointer");

    return(tng_util_generic_write(tng_data, frame_nr, velocities, 3,
                                  TNG_TRAJ_VELOCITIES, "VELOCITIES",
                                  TNG_PARTICLE_BLOCK_DATA,
                                  TNG_TNG_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_double_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double *velocities)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(velocities, "TNG library: velocities must not be a NULL pointer");

    return(tng_util_generic_double_write(tng_data, frame_nr, velocities, 3,
                                         TNG_TRAJ_VELOCITIES, "VELOCITIES",
                                         TNG_PARTICLE_BLOCK_DATA,
                                         TNG_TNG_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const float *forces)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(forces, "TNG library: forces must not be a NULL pointer");

    return(tng_util_generic_write(tng_data, frame_nr, forces, 3,
                                  TNG_TRAJ_FORCES, "FORCES",
                                  TNG_PARTICLE_BLOCK_DATA,
                                  TNG_GZIP_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_double_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double *forces)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(forces, "TNG library: forces must not be a NULL pointer");

    return(tng_util_generic_double_write(tng_data, frame_nr, forces, 3,
                                         TNG_TRAJ_FORCES, "FORCES",
                                         TNG_PARTICLE_BLOCK_DATA,
                                         TNG_GZIP_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const float *box_shape)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(box_shape, "TNG library: box_shape must not be a NULL pointer");

    return(tng_util_generic_write(tng_data, frame_nr, box_shape, 9,
                                  TNG_TRAJ_BOX_SHAPE, "BOX SHAPE",
                                  TNG_NON_PARTICLE_BLOCK_DATA,
                                  TNG_GZIP_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_double_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double *box_shape)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(box_shape, "TNG library: box_shape must not be a NULL pointer");

    return(tng_util_generic_double_write(tng_data, frame_nr, box_shape, 9,
                                         TNG_TRAJ_BOX_SHAPE, "BOX SHAPE",
                                         TNG_NON_PARTICLE_BLOCK_DATA,
                                         TNG_GZIP_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_generic_with_time_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const float *values,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression)
{
    tng_trajectory_frame_set_t frame_set;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(time >= 0, "TNG library: time must be >= 0.");
    TNG_ASSERT(values, "TNG library: values must not be a NULL pointer");

    stat = tng_util_generic_write(tng_data, frame_nr, values, n_values_per_frame,
                                  block_id, block_name,
                                  particle_dependency,
                                  compression);

    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    frame_set = &tng_data->current_trajectory_frame_set;

    /* first_frame_time is -1 when it is not yet set. */
    if(frame_set->first_frame_time < -0.1)
    {
        if(frame_nr > frame_set->first_frame)
        {
            stat = tng_frame_set_first_frame_time_set(tng_data,
                                                      time -
                                                      (frame_nr -
                                                       frame_set->first_frame) *
                                                      tng_data->time_per_frame);
        }
        else
        {
            stat = tng_frame_set_first_frame_time_set(tng_data, time);
        }
    }
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_util_generic_with_time_double_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const double *values,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression)
{
    tng_trajectory_frame_set_t frame_set;
    tng_function_status stat;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(time >= 0, "TNG library: time must be >= 0.");
    TNG_ASSERT(values, "TNG library: values must not be a NULL pointer");

    stat = tng_util_generic_double_write(tng_data, frame_nr, values, n_values_per_frame,
                                         block_id, block_name,
                                         particle_dependency,
                                         compression);

    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    frame_set = &tng_data->current_trajectory_frame_set;

    /* first_frame_time is -1 when it is not yet set. */
    if(frame_set->first_frame_time < -0.1)
    {
        if(frame_nr > frame_set->first_frame)
        {
            stat = tng_frame_set_first_frame_time_set(tng_data,
                                                      time -
                                                      (frame_nr -
                                                       frame_set->first_frame) *
                                                      tng_data->time_per_frame);
        }
        else
        {
            stat = tng_frame_set_first_frame_time_set(tng_data, time);
        }
    }
    return(stat);
}

tng_function_status DECLSPECDLLEXPORT tng_util_pos_with_time_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const float *positions)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(time >= 0, "TNG library: time must be >= 0.");
    TNG_ASSERT(positions, "TNG library: positions must not be a NULL pointer");

    return(tng_util_generic_with_time_write(tng_data, frame_nr, time, positions,
                                            3, TNG_TRAJ_POSITIONS, "POSITIONS",
                                            TNG_PARTICLE_BLOCK_DATA,
                                            TNG_TNG_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_pos_with_time_double_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const double *positions)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(time >= 0, "TNG library: time must be >= 0.");
    TNG_ASSERT(positions, "TNG library: positions must not be a NULL pointer");

    return(tng_util_generic_with_time_double_write(tng_data, frame_nr, time,
                                                   positions, 3,
                                                   TNG_TRAJ_POSITIONS,
                                                   "POSITIONS",
                                                   TNG_PARTICLE_BLOCK_DATA,
                                                   TNG_TNG_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_with_time_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const float *velocities)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(time >= 0, "TNG library: time must be >= 0.");
    TNG_ASSERT(velocities, "TNG library: velocities must not be a NULL pointer");

    return(tng_util_generic_with_time_write(tng_data, frame_nr, time,
                                            velocities, 3,
                                            TNG_TRAJ_VELOCITIES,
                                            "VELOCITIES",
                                            TNG_PARTICLE_BLOCK_DATA,
                                            TNG_TNG_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_vel_with_time_double_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const double *velocities)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(time >= 0, "TNG library: time must be >= 0.");
    TNG_ASSERT(velocities, "TNG library: velocities must not be a NULL pointer");

    return(tng_util_generic_with_time_double_write(tng_data, frame_nr, time,
                                                   velocities, 3,
                                                   TNG_TRAJ_VELOCITIES,
                                                   "VELOCITIES",
                                                   TNG_PARTICLE_BLOCK_DATA,
                                                   TNG_TNG_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_with_time_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const float *forces)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(time >= 0, "TNG library: time must be >= 0.");
    TNG_ASSERT(forces, "TNG library: forces must not be a NULL pointer");

    return(tng_util_generic_with_time_write(tng_data, frame_nr, time, forces,
                                            3, TNG_TRAJ_FORCES, "FORCES",
                                            TNG_PARTICLE_BLOCK_DATA,
                                            TNG_GZIP_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_force_with_time_double_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const double *forces)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(time >= 0, "TNG library: time must be >= 0.");
    TNG_ASSERT(forces, "TNG library: forces must not be a NULL pointer");

    return(tng_util_generic_with_time_double_write(tng_data, frame_nr, time,
                                                   forces, 3,
                                                   TNG_TRAJ_FORCES, "FORCES",
                                                   TNG_PARTICLE_BLOCK_DATA,
                                                   TNG_GZIP_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_with_time_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const float *box_shape)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(time >= 0, "TNG library: time must be >= 0.");
    TNG_ASSERT(box_shape, "TNG library: box_shape must not be a NULL pointer");

    return(tng_util_generic_with_time_write(tng_data, frame_nr, time, box_shape,
                                            9, TNG_TRAJ_BOX_SHAPE, "BOX SHAPE",
                                            TNG_NON_PARTICLE_BLOCK_DATA,
                                            TNG_GZIP_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_with_time_double_write
                (tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const double *box_shape)
{
    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(frame_nr >= 0, "TNG library: frame_nr must be >= 0.");
    TNG_ASSERT(time >= 0, "TNG library: time must be >= 0.");
    TNG_ASSERT(box_shape, "TNG library: box_shape must not be a NULL pointer");

    return(tng_util_generic_with_time_double_write(tng_data, frame_nr,
                                                   time, box_shape, 9,
                                                   TNG_TRAJ_BOX_SHAPE,
                                                   "BOX SHAPE",
                                                   TNG_NON_PARTICLE_BLOCK_DATA,
                                                   TNG_GZIP_COMPRESSION));
}

tng_function_status DECLSPECDLLEXPORT tng_util_frame_current_compression_get
                (tng_trajectory_t tng_data,
                 const int64_t block_id,
                 int64_t *codec_id,
                 double *factor)
{
    tng_trajectory_frame_set_t frame_set;
    tng_particle_data_t p_data = 0;
    tng_non_particle_data_t np_data = 0;
    tng_function_status stat;
    int64_t i;
    int block_type = -1;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(codec_id, "TNG library: The pointer to the returned codec id must not be a NULL pointer.");
    TNG_ASSERT(factor, "TNG library: The pointer to the returned multiplication factor must not be a NULL pointer.");

    frame_set = &tng_data->current_trajectory_frame_set;

    stat = tng_particle_data_find(tng_data, block_id, &p_data);
    if(stat == TNG_SUCCESS)
    {
        block_type = TNG_PARTICLE_BLOCK_DATA;
    }
    else
    {
        stat = tng_data_find(tng_data, block_id, &np_data);
        if(stat == TNG_SUCCESS)
        {
            block_type = TNG_NON_PARTICLE_BLOCK_DATA;
        }
        else
        {
            stat = tng_frame_set_read_current_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
            if(stat != TNG_SUCCESS)
            {
                return(stat);
            }
            stat = tng_particle_data_find(tng_data, block_id, &p_data);
            if(stat == TNG_SUCCESS)
            {
                block_type = TNG_PARTICLE_BLOCK_DATA;
            }
            else
            {
                stat = tng_data_find(tng_data, block_id, &np_data);
                if(stat == TNG_SUCCESS)
                {
                    block_type = TNG_NON_PARTICLE_BLOCK_DATA;
                }
                else
                {
                    return(stat);
                }
            }
        }
    }
    if(block_type == TNG_PARTICLE_BLOCK_DATA)
    {
        if(p_data->last_retrieved_frame < 0)
        {
            i = p_data->first_frame_with_data;
        }
        else
        {
            i = p_data->last_retrieved_frame;
        }
    }
    else if(block_type == TNG_NON_PARTICLE_BLOCK_DATA)
    {
        if(np_data->last_retrieved_frame < 0)
        {
            i = np_data->first_frame_with_data;
        }
        else
        {
            i = np_data->last_retrieved_frame;
        }
    }
    else
    {
        return(TNG_FAILURE);
    }
    if(i < frame_set->first_frame || i >= frame_set->first_frame + frame_set->n_frames)
    {
        stat = tng_frame_set_of_frame_find(tng_data, i);
        if(stat != TNG_SUCCESS)
        {
            return(stat);
        }
        stat = tng_frame_set_read_current_only_data_from_block_id(tng_data, TNG_USE_HASH, block_id);
        if(stat != TNG_SUCCESS)
        {
            fprintf(stderr, "TNG library: Cannot read data block of frame set. %s: %d\n",
                __FILE__, __LINE__);
            return(stat);
        }
    }
    if(block_type == TNG_PARTICLE_BLOCK_DATA)
    {
        *codec_id = p_data->codec_id;
        *factor   = p_data->compression_multiplier;
    }
    else if(block_type == TNG_NON_PARTICLE_BLOCK_DATA)
    {
        *codec_id = np_data->codec_id;
        *factor   = np_data->compression_multiplier;
    }
    return(TNG_SUCCESS);
}

tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_next_frame_present_data_blocks_find
                (tng_trajectory_t tng_data,
                 int64_t current_frame,
                 const int64_t n_requested_data_block_ids,
                 const int64_t *requested_data_block_ids,
                 int64_t *next_frame,
                 int64_t *n_data_blocks_in_next_frame,
                 int64_t **data_block_ids_in_next_frame)
{
    tng_trajectory_frame_set_t frame_set;
    tng_function_status stat;
    tng_particle_data_t p_data;
    tng_non_particle_data_t np_data;
    tng_gen_block_t block;
    int64_t i, j, block_id, *temp;
    int64_t data_frame, frame_diff, min_diff;
    int64_t size, frame_set_file_pos;
    int found, read_all = 0;
    long file_pos;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(next_frame, "TNG library: The pointer to the next frame must not be NULL.");
    TNG_ASSERT(n_data_blocks_in_next_frame, "TNG library: The pointer to n_data_blocks_in_next_frame must not be NULL.");
    TNG_ASSERT(data_block_ids_in_next_frame, "TNG library: The pointer to the list of data block IDs must not be NULL.");

    if(n_requested_data_block_ids)
    {
        TNG_ASSERT(requested_data_block_ids, "TNG library: If the number of requested data blocks is > 0 then the array of data block IDs must not be NULL.");
        size = sizeof(int64_t) * n_requested_data_block_ids;
        temp = realloc(*data_block_ids_in_next_frame, size);
        if(!temp)
        {
            fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                    sizeof(int64_t) * (*n_data_blocks_in_next_frame),
                    __FILE__, __LINE__);
            free(*data_block_ids_in_next_frame);
            *data_block_ids_in_next_frame = 0;
            return(TNG_CRITICAL);
        }
        *data_block_ids_in_next_frame = temp;
    }

    frame_set = &tng_data->current_trajectory_frame_set;

    current_frame += 1;

    if(current_frame < frame_set->first_frame ||
       current_frame >= frame_set->first_frame + frame_set->n_frames)
    {
        frame_set_file_pos = tng_data->current_trajectory_frame_set_input_file_pos;
        stat = tng_frame_set_of_frame_find(tng_data, current_frame);
        if(stat != TNG_SUCCESS)
        {
            /* If the frame set search found the frame set after the starting
             * frame set there is a gap in the frame sets. So, even if the frame
             * was not found the next frame with data is still in the found
             * frame set. */
            if(stat == TNG_CRITICAL || frame_set->prev_frame_set_file_pos !=
               frame_set_file_pos)
            {
                return(stat);
            }
            current_frame = frame_set->first_frame;
        }
    }

    /* Check for data blocks only if they have not already been found. */
    if(frame_set->n_particle_data_blocks <= 0 && frame_set->n_data_blocks <= 0)
    {
        file_pos = ftell(tng_data->input_file);
        if(file_pos < tng_data->input_file_len)
        {
            tng_block_init(&block);
            stat = tng_block_header_read(tng_data, block);
            while(file_pos < tng_data->input_file_len &&
                stat != TNG_CRITICAL &&
                block->id != TNG_TRAJECTORY_FRAME_SET &&
                block->id != -1)
            {
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
                fprintf(stderr, "TNG library: Cannot read block header at pos %ld. %s: %d\n",
                        file_pos, __FILE__, __LINE__);
                return(stat);
            }
        }
        read_all = 1;
    }

    min_diff = -1;

    *n_data_blocks_in_next_frame = 0;

    for(i = 0; i < frame_set->n_particle_data_blocks; i++)
    {
        p_data = &frame_set->tr_particle_data[i];
        block_id = p_data->block_id;

        if(n_requested_data_block_ids > 0)
        {
            found = 0;
            for(j = 0; j < n_requested_data_block_ids; j++)
            {
                if(block_id == requested_data_block_ids[j])
                {
                    found = 1;
                    break;
                }
            }
            if(!found)
            {
                continue;
            }
        }

        if(!read_all && (p_data->last_retrieved_frame < frame_set->first_frame ||
           p_data->last_retrieved_frame >=
           frame_set->first_frame + frame_set->n_frames))
        {
            stat = tng_frame_set_read_current_only_data_from_block_id(tng_data,
                                                                      TNG_USE_HASH, block_id);
            if(stat == TNG_CRITICAL)
            {
                fprintf(stderr, "TNG library: Cannot read data block of frame set. %s: %d\n",
                    __FILE__, __LINE__);
                return(stat);
            }
            if(stat == TNG_FAILURE)
            {
                continue;
            }
        }
        if(frame_set->first_frame != current_frame &&
           p_data->last_retrieved_frame >= 0)
        {
            data_frame = p_data->last_retrieved_frame + p_data->stride_length;
        }
        else
        {
            data_frame = p_data->first_frame_with_data;
        }
        frame_diff = data_frame - current_frame;
        if(frame_diff < 0)
        {
            continue;
        }
        if(min_diff == -1 || frame_diff <= min_diff)
        {
            if(frame_diff < min_diff)
            {
                *n_data_blocks_in_next_frame = 1;
            }
            else
            {
                *n_data_blocks_in_next_frame += 1;
            }
            if(n_requested_data_block_ids <= 0)
            {
                size = sizeof(int64_t) * (*n_data_blocks_in_next_frame);
                temp = realloc(*data_block_ids_in_next_frame, size);
                if(!temp)
                {
                    fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                           sizeof(int64_t) * (*n_data_blocks_in_next_frame),
                           __FILE__, __LINE__);
                    free(*data_block_ids_in_next_frame);
                    *data_block_ids_in_next_frame = 0;
                    return(TNG_CRITICAL);
                }
                *data_block_ids_in_next_frame = temp;
            }
            else
            {
                TNG_ASSERT(*n_data_blocks_in_next_frame <= n_requested_data_block_ids, "TNG library: Array of data block IDs out of bounds");
            }
            (*data_block_ids_in_next_frame)[(*n_data_blocks_in_next_frame) - 1] = block_id;

            min_diff = frame_diff;
        }
    }
    for(i = 0; i < frame_set->n_data_blocks; i++)
    {
        np_data = &frame_set->tr_data[i];
        block_id = np_data->block_id;

        if(n_requested_data_block_ids > 0)
        {
            found = 0;
            for(j = 0; j < n_requested_data_block_ids; j++)
            {
                if(block_id == requested_data_block_ids[j])
                {
                    found = 1;
                    break;
                }
            }
            if(!found)
            {
                continue;
            }
        }

        if(!read_all && (np_data->last_retrieved_frame < frame_set->first_frame ||
           np_data->last_retrieved_frame >=
           frame_set->first_frame + frame_set->n_frames))
        {
            stat = tng_frame_set_read_current_only_data_from_block_id(tng_data,
                                                                      TNG_USE_HASH, block_id);
            if(stat == TNG_CRITICAL)
            {
                fprintf(stderr, "TNG library: Cannot read data block of frame set. %s: %d\n",
                    __FILE__, __LINE__);
                return(stat);
            }
            if(stat == TNG_FAILURE)
            {
                continue;
            }
        }
        if(frame_set->first_frame != current_frame &&
           np_data->last_retrieved_frame >= 0)
        {
            data_frame = np_data->last_retrieved_frame + np_data->stride_length;
        }
        else
        {
            data_frame = np_data->first_frame_with_data;
        }
        frame_diff = data_frame - current_frame;
        if(frame_diff < 0)
        {
            continue;
        }
        if(min_diff == -1 || frame_diff <= min_diff)
        {
            if(frame_diff < min_diff)
            {
                *n_data_blocks_in_next_frame = 1;
            }
            else
            {
                *n_data_blocks_in_next_frame += 1;
            }
            if(n_requested_data_block_ids <= 0)
            {
                size = sizeof(int64_t) * (*n_data_blocks_in_next_frame);
                temp = realloc(*data_block_ids_in_next_frame, size);
                if(!temp)
                {
                    fprintf(stderr, "TNG library: Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                           sizeof(int64_t) * (*n_data_blocks_in_next_frame),
                           __FILE__, __LINE__);
                    free(*data_block_ids_in_next_frame);
                    *data_block_ids_in_next_frame = 0;
                    return(TNG_CRITICAL);
                }
                *data_block_ids_in_next_frame = temp;
            }
            else
            {
                TNG_ASSERT(*n_data_blocks_in_next_frame <= n_requested_data_block_ids, "TNG library: Array of data block IDs out of bounds");
            }
            (*data_block_ids_in_next_frame)[(*n_data_blocks_in_next_frame) - 1] = block_id;

            min_diff = frame_diff;
        }
    }
    if(min_diff < 0)
    {
        return(TNG_FAILURE);
    }
    *next_frame = current_frame + min_diff;

    return(TNG_SUCCESS);
}

/*
tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_all_data_block_types_get
                (tng_trajectory_t tng_data,
                 int64_t *n_data_blocks,
                 int64_t **data_block_ids,
                 char ***data_block_names,
                 int64_t **stride_lengths,
                 int64_t **n_values_per_frame,
                 char **block_types,
                 char **dependencies,
                 char **compressions)
{
    tng_gen_block_t block;
    long orig_file_pos, file_pos;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(n_data_blocks, "TNG library: The pointer to n_data_blocks must not be NULL.");
    TNG_ASSERT(data_block_ids, "TNG library: The pointer to the list of data block IDs must not be NULL.");
    TNG_ASSERT(data_block_names, "TNG library: The pointer to the list of data block names must not be NULL.");
    TNG_ASSERT(stride_lengths, "TNG library: The pointer to the list of stride lengths must not be NULL.");

    orig_file_pos = ftell(tng_data->input_file);

    if(!tng_data->input_file_len)
    {
        fseek(tng_data->input_file, 0, SEEK_END);
        tng_data->input_file_len = ftell(tng_data->input_file);
    }

    fseek(tng_data->input_file, 0, SEEK_SET);
    file_pos = 0;

    *n_data_blocks = 0;

    tng_block_init(&block);

    while(file_pos < tng_data->input_file_len &&
          tng_block_header_read(tng_data, block) != TNG_CRITICAL)
    {
        if(block->id > TNG_TRAJECTORY_FRAME_SET)
        {

        }
        file_pos += (long)(block->block_contents_size + block->header_contents_size);
        fseek(tng_data->input_file, (long)block->block_contents_size, SEEK_CUR);
    }

    fseek(tng_data->input_file, orig_file_pos, SEEK_SET);

    return(TNG_SUCCESS);
}
*/
tng_function_status DECLSPECDLLEXPORT tng_util_prepare_append_after_frame
                (tng_trajectory_t tng_data,
                 const int64_t prev_frame)
{
    tng_function_status stat;
    FILE *temp = tng_data->input_file;

    TNG_ASSERT(tng_data, "TNG library: Trajectory container not properly setup.");
    TNG_ASSERT(prev_frame >= 0, "TNG library: The previous frame must not be negative.");

    tng_data->input_file = tng_data->output_file;

    stat = tng_frame_set_of_frame_find(tng_data, prev_frame);
    if(stat != TNG_SUCCESS)
    {
        return(stat);
    }

    tng_data->current_trajectory_frame_set_output_file_pos =
    tng_data->current_trajectory_frame_set_input_file_pos;

    tng_data->input_file = temp;

    return(TNG_SUCCESS);
}

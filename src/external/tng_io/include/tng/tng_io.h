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

/** @file tng_io.h
 *  @brief API for input and output of tng trajectory files
 *  @mainpage TNG: A flexible binary trajectory format
 *  @section intro_sec Introduction
 *
 * The TNG format is developed as part of the ScalaLife EU project.
 * It is flexible by design to allow parallel writing, custom data blocks,
 * different output frequencies and different compression algorithms.
 *
 * Each block can contain MD5 hashes to verify data integrity and the file
 * can be signed by the user to ensure that the origin is correct.
 *
 * The intention is that the API and ABI should be stable, but it is
 * still possible that future changes might make that impossible, in which
 * case that will be clarified.
 *
 * The API and all examples are released without any warranties. Use them at
 * your own risk.
 *
 * @section authors_sec Authors
 *
 * The TNG trajectory format is developed by:
 *
 * Magnus Lundborg magnus.lundborg@scilifelab.se
 *
 * Daniel Spångberg daniels@mkem.uu.se
 *
 * Rossen Apostolov rossen@kth.se
 *
 * The API is implemented mainly by:
 *
 * Magnus Lundborg
 *
 * @section License
 *
 * Copyright (c) 2012, The GROMACS development team.
 * check out http://www.gromacs.org for more information.
 *
 * The TNG API is released under the Revised BSD License and is free to
 * redistribute according to that license.
 *
 * A license file (named COPYING) should be included with each copy of the API.
 *
 * @section install_sec Installation
 *
 * \code
 * mkdir build
 *
 * cd build
 *
 * cmake ..
 *
 * make
 *
 * make install
 * \endcode
 * Test by running:
 * \code
 * bin/tests/tng_testing
 * \endcode
 *
 * @section change_sec Change Log
 *
 * See git log for full revision history.
 *
 * Revisions
 *
 * v. 1.7 - Fifth stable release of the API
 *
 *        - Added function tng_util_num_frames_with_data_of_block_id_get().
 *        - Merged some functions and data structures
 *          to make less difference between data blocks.
 *        - Bugs fixed
 *
 * v. 1.6 - Fourth stable release of the API.
 *
 *        - Removed OpenMP option when building.
 *        - Functionality for migrating data blocks.
 *        - Improved handling of molecules.
 *        - Improved installation of TNG documentation.
 *        - Enhancements to CMake usage.
 *        - Required CMake version raised to 2.8.8.
 *        - Bugs fixed.
 *
 * v. 1.5 - Third stable release of the API.
 *
 *        - Fortran wrapper split into separate file
 *        - Added more block IDs.
 *        - Some new functions and utility functions added.
 *        - Improved compression precision settings.
 *        - Improved tests.
 *        - Make appending to file work better.
 *        - Modified CMake settings
 *        - Bugs fixed
 *
 * v. 1.4 - Changed from LGPL to the Revised BSD License.
 *
 *        - More flexible support for digital signatures in header.
 *        - Block ID numbers changed.
 *
 * v. 1.3 - Second stable release of the API.
 *
 *      - Added multiplication factor for coordinate units to general info.
 *      - Added time stamps and time per frame in frame sets.
 *      - High-level API functions added (not for managing molecules yet)
 *      - Added functions for reading data blocks into 1D arrays.
 *      - TNG compression added.
 *      - C++ interface added.
 *      - Avoid memory allocation if no data is submitted when adding data
 *        blocks.
 *      - Added function tng_num_frames_per_frame_set_set
 *      - Added data block IDs for charges, b-factors and occupancy.
 *      - GZIP compression added.
 *      - Fixed bug when updating MD5 hashes of data blocks.
 *      - Fixed bug in chain_name_of_particle_get(...)
 *      - Update frame set pointers properly.
 *      - Moved fortran wrapper from header file to source file.
 *      - Write sparse data in mdrun examples.
 *      - Fixed bugs related to reading and writing sparse data.
 *      - Fixed memory leak for non-trajectory particle data blocks.
 *      - Fixed bug when writing data blocks.
 *      - Fixed wrong values in dependency constants
 *      - Write box shape, partial charges and annotation data in tng_testing
 *      - Bug fixes in tng_testing (frame sets not written before)
 *
 * v. 1.0 - First stable release of the API.
 *
 *
 * @section examples_sec Examples
 *
 * There are some examples of how to use the library located in src/tests/
 *
 * @subsection tng_subsec TNG files
 *
 * The build directory contains an example_files directory, which in turn
 * contains a very short example of a TNG file containing a few water molecules,
 * a box shape description and positions in 10 frames.
 *
 * It is also possible to run the bin/examples/md_openmp_util
 * (see src/tests/md_openmp_util.c)
 * testing program, which will save MD simulations output to a new file
 * (saved in the example_files directory).
 *
 * These files can be read using the bin/examples/tng_io_read_pos_util
 * program.
 *
 * @subsection c_subsec C
 *
 * Example writing data to a TNG file (just an excerpt):
 * \code
 *     for ( step = 1; step < step_num; step++ )
 *     {
 *         compute ( np, nd, pos, vel, mass, force, &potential, &kinetic );
 *
 *         if(step % step_save == 0)
 *         {
 *             // Write positions, velocities and forces
 *             if(tng_util_pos_write(traj, step, pos) != TNG_SUCCESS)
 *             {
 *                 printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
 *                 break;
 *             }
 *             if(tng_util_vel_write(traj, step, vel) != TNG_SUCCESS)
 *             {
 *                 printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
 *                 break;
 *             }
 *             if(tng_util_force_write(traj, step, force) != TNG_SUCCESS)
 *             {
 *                 printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
 *                 break;
 *             }
 *         }
 *         update ( np, nd, pos, vel, force, acc, mass, dt );
 *     }
 * \endcode
 *
 * Example reading positions from a TNG file:
 * \code
 * #include <stdlib.h>
 * #include <stdio.h>
 * #include "tng/tng_io.h"
 *
 * int main(int argc, char **argv)
 * {
 *     tng_trajectory_t traj;
 *     // Assume that the data is stored as floats. The data is placed in 1-D
 *     // arrays
 *     float *positions = 0, *box_shape = 0;
 *     int64_t n_particles, n_frames, tot_n_frames, stride_length, i, j;
 *     // Set a default frame range
 *     int64_t first_frame = 0, last_frame = 5000;
 *     int k;
 *
 *     // A reference must be passed to allocate memory
 *     tng_util_trajectory_open(argv[1], 'r', &traj);
 *
 *     if(tng_num_frames_get(traj, &tot_n_frames) != TNG_SUCCESS)
 *     {
 *         printf("Cannot determine the number of frames in the file\n");
 *         tng_util_trajectory_close(&traj);
 *         exit(1);
 *     }
 *
 *     if(tng_num_particles_get(traj, &n_particles) != TNG_SUCCESS)
 *     {
 *         printf("Cannot determine the number of particles in the file\n");
 *         tng_util_trajectory_close(&traj);
 *         exit(1);
 *     }
 *
 *     printf("%"PRId64" frames in file\n", tot_n_frames);
 *
 *     if(last_frame > tot_n_frames - 1)
 *     {
 *         last_frame = tot_n_frames - 1;
 *     }
 *
 *     if(tng_util_box_shape_read(traj, &box_shape, &stride_length) ==
 *         TNG_SUCCESS)
 *     {
 *         printf("Simulation box shape: ");
 *         for(i=0; i < 9; i++)
 *         {
 *             printf("%f ", box_shape[i]);
 *         }
 *         printf("\n");
 *     }
 *     else
 *     {
 *         printf("Simulation box shape not set in the file (or could not be read)\n");
 *     }
 *
 *     n_frames = last_frame - first_frame + 1;
 *
 *
 *     // Get the positions of all particles in the requested frame range.
 *     // The positions are stored in the positions array.
 *     // N.B. No proper error checks.
 *     if(tng_util_pos_read_range(traj, 0, last_frame, &positions, &stride_length)
 *        == TNG_SUCCESS)
 *     {
 *         // Print the positions of the wanted particle (zero based)
 *         for(i=0; i < n_frames; i += stride_length)
 *         {
 *             printf("\nFrame %"PRId64":\n", first_frame + i);
 *             for(j=0; j < n_particles; j++)
 *             {
 *                 printf("Atom nr: %"PRId64"", j);
 *                 for(k=0; k < 3; k++)
 *                 {
 *                     printf("\t%f", positions[i/stride_length*n_particles*
 *                                              3+j*3+k]);
 *                 }
 *                 printf("\n");
 *             }
 *         }
 *     }
 *     else
 *     {
 *         printf("Cannot read positions\n");
 *     }
 *
 *     // Free memory
 *     if(positions)
 *     {
 *         free(positions);
 *     }
 *     tng_util_trajectory_close(&traj);
 *
 *     return(0);
 * }
 *
 * \endcode
 *
 * @subsection fortran_subsec Fortran
 *
 * The TNG library can be used from Fortran. It requires cray pointers, which
 * are not part of the Fortran 77 standard, but available in most compilers.
 *
 * To compile the fortran example -DTNG_BUILD_FORTRAN=ON needs to be specified when
 * running cmake.
 *
 */

#ifndef TNG_IO_H
#define TNG_IO_H     1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "tng_io_fwd.h"

#ifdef USE_STD_INTTYPES_H
#include <inttypes.h>
#else
/* Visual Studio does not contain inttypes.h and stdint.h. Some defines and
 * typedefs are used from the GNU C Library */
#ifdef _MSC_VER

typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;

#else
#include <stdint.h>
#endif /* _MSC_VER */

/* This is from inttypes.h  (GNU C Library) */
/* The ISO C99 standard specifies that these macros must only be
   defined if explicitly requested.  */
#if !defined __cplusplus || defined __STDC_FORMAT_MACROS

# if __WORDSIZE == 64
#  define __PRI64_PREFIX        "l"
#  define __PRIPTR_PREFIX       "l"
# else
#  define __PRI64_PREFIX        "ll"
#  define __PRIPTR_PREFIX
# endif

/* From stdint.h (GNU C Library) */
/* Macros for printing format specifiers. */
/* Decimal notation.  */
#ifndef PRId64
# define PRId64         __PRI64_PREFIX "d"
#endif

#ifndef PRIu64
# define PRIu64         __PRI64_PREFIX "u"
#endif

#ifndef PRIuPTR
# define PRIuPTR         __PRIPTR_PREFIX "u"
#endif

#endif

#endif /* USE_STD_INTTYPES_H */

#ifndef USE_WINDOWS
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#define USE_WINDOWS
#endif /* win32... */
#endif /* not defined USE_WINDOWS */

#ifndef DECLSPECDLLEXPORT
#ifdef USE_WINDOWS
#define DECLSPECDLLEXPORT __declspec(dllexport)
#else /* USE_WINDOWS */
#define DECLSPECDLLEXPORT
#endif /* USE_WINDOWS */
#endif /* DECLSPECDLLEXPORT */

/** Flag to indicate frame dependent data. */
#define TNG_FRAME_DEPENDENT 1
/** Flag to indicate particle dependent data. */
#define TNG_PARTICLE_DEPENDENT 2

/** The maximum length of a date string */
#define TNG_MAX_DATE_STR_LEN 24
/** The length of an MD5 hash */
#define TNG_MD5_HASH_LEN 16
/** The maximum allowed length of a string */
#define TNG_MAX_STR_LEN 1024

#ifndef NDEBUG
#define TNG_ASSERT(cnd, msg) if(!(cnd)) {printf("%s\n", msg); assert(cnd);}
#else
#define TNG_ASSERT(cnd, msg) (void)0;
#endif

/** Flag to specify the endianness of a TNG file */
typedef enum {TNG_BIG_ENDIAN,
              TNG_LITTLE_ENDIAN} tng_file_endianness;

/** Flag to specify the endianness of 32 bit values of the current architecture. */
typedef enum {TNG_BIG_ENDIAN_32,
              TNG_LITTLE_ENDIAN_32,
              TNG_BYTE_PAIR_SWAP_32} tng_endianness_32;

/** Flag to specify the endianness of 64 bit values of the current architecture. */
typedef enum {TNG_BIG_ENDIAN_64,
              TNG_LITTLE_ENDIAN_64,
              TNG_QUAD_SWAP_64,
              TNG_BYTE_PAIR_SWAP_64,
              TNG_BYTE_SWAP_64} tng_endianness_64;

/** Compression mode is specified in each data block */
typedef enum {TNG_UNCOMPRESSED,
              TNG_XTC_COMPRESSION,
              TNG_TNG_COMPRESSION,
              TNG_GZIP_COMPRESSION} tng_compression;

/** Hash types */
typedef enum {TNG_NO_HASH,
              TNG_MD5,
              TNG_SHA256} tng_hash_type;

/** Non trajectory blocks come before the first frame set block */
typedef enum {TNG_NON_TRAJECTORY_BLOCK, TNG_TRAJECTORY_BLOCK} tng_block_type;

/** @defgroup def1 Standard non-trajectory blocks
 *  Block IDs of standard non-trajectory blocks.
 * @{
 */
#define TNG_GENERAL_INFO                0x0000000000000000LL
#define TNG_MOLECULES                   0x0000000000000001LL
#define TNG_TRAJECTORY_FRAME_SET        0x0000000000000002LL
#define TNG_PARTICLE_MAPPING            0x0000000000000003LL
/** @} */

/** @defgroup def2 Standard trajectory blocks
 * Block IDs of standard trajectory blocks. Box shape and partial charges can
 * be either trajectory blocks or non-trajectory blocks
 * @{
 */
#define TNG_TRAJ_BOX_SHAPE              0x0000000010000000LL
#define TNG_TRAJ_POSITIONS              0x0000000010000001LL
#define TNG_TRAJ_VELOCITIES             0x0000000010000002LL
#define TNG_TRAJ_FORCES                 0x0000000010000003LL
#define TNG_TRAJ_PARTIAL_CHARGES        0x0000000010000004LL
#define TNG_TRAJ_FORMAL_CHARGES         0x0000000010000005LL
#define TNG_TRAJ_B_FACTORS              0x0000000010000006LL
#define TNG_TRAJ_ANISOTROPIC_B_FACTORS  0x0000000010000007LL
#define TNG_TRAJ_OCCUPANCY              0x0000000010000008LL
#define TNG_TRAJ_GENERAL_COMMENTS       0x0000000010000009LL
/** @} */


/** @defgroup def3 GROMACS data block IDs
 *  Block IDs of data blocks specific to GROMACS.
 * @{
 */
#define TNG_GMX_LAMBDA                  0x1000000010000000LL
#define TNG_GMX_ENERGY_ANGLE            0x1000000010000001LL
#define TNG_GMX_ENERGY_RYCKAERT_BELL    0x1000000010000002LL
#define TNG_GMX_ENERGY_LJ_14            0x1000000010000003LL
#define TNG_GMX_ENERGY_COULOMB_14       0x1000000010000004LL
#define TNG_GMX_ENERGY_LJ_(SR)          0x1000000010000005LL
#define TNG_GMX_ENERGY_COULOMB_(SR)     0x1000000010000006LL
#define TNG_GMX_ENERGY_COUL_RECIP       0x1000000010000007LL
#define TNG_GMX_ENERGY_POTENTIAL        0x1000000010000008LL
#define TNG_GMX_ENERGY_KINETIC_EN       0x1000000010000009LL
#define TNG_GMX_ENERGY_TOTAL_ENERGY     0x1000000010000010LL
#define TNG_GMX_ENERGY_TEMPERATURE      0x1000000010000011LL
#define TNG_GMX_ENERGY_PRESSURE         0x1000000010000012LL
#define TNG_GMX_ENERGY_CONSTR_RMSD      0x1000000010000013LL
#define TNG_GMX_ENERGY_BOX_X            0x1000000010000014LL
#define TNG_GMX_ENERGY_BOX_Y            0x1000000010000015LL
#define TNG_GMX_ENERGY_BOX_Z            0x1000000010000016LL
#define TNG_GMX_ENERGY_VOLUME           0x1000000010000017LL
#define TNG_GMX_ENERGY_DENSITY          0x1000000010000018LL
#define TNG_GMX_ENERGY_PV               0x1000000010000019LL
#define TNG_GMX_ENERGY_ENTHALPY         0x1000000010000020LL
#define TNG_GMX_ENERGY_VIR_XX           0x1000000010000021LL
#define TNG_GMX_ENERGY_VIR_XY           0x1000000010000022LL
#define TNG_GMX_ENERGY_VIR_XZ           0x1000000010000023LL
#define TNG_GMX_ENERGY_VIR_YX           0x1000000010000024LL
#define TNG_GMX_ENERGY_VIR_YY           0x1000000010000025LL
#define TNG_GMX_ENERGY_VIR_YZ           0x1000000010000026LL
#define TNG_GMX_ENERGY_VIR_ZX           0x1000000010000027LL
#define TNG_GMX_ENERGY_VIR_ZY           0x1000000010000028LL
#define TNG_GMX_ENERGY_VIR_ZZ           0x1000000010000029LL
#define TNG_GMX_ENERGY_PRES_XX          0x1000000010000030LL
#define TNG_GMX_ENERGY_PRES_XY          0x1000000010000031LL
#define TNG_GMX_ENERGY_PRES_XZ          0x1000000010000032LL
#define TNG_GMX_ENERGY_PRES_YX          0x1000000010000033LL
#define TNG_GMX_ENERGY_PRES_YY          0x1000000010000034LL
#define TNG_GMX_ENERGY_PRES_YZ          0x1000000010000035LL
#define TNG_GMX_ENERGY_PRES_ZX          0x1000000010000036LL
#define TNG_GMX_ENERGY_PRES_ZY          0x1000000010000037LL
#define TNG_GMX_ENERGY_PRES_ZZ          0x1000000010000038LL
#define TNG_GMX_ENERGY_SURFXSURFTEN     0x1000000010000039LL
#define TNG_GMX_ENERGY_T_SYSTEM         0x1000000010000040LL
#define TNG_GMX_ENERGY_LAMB_SYSTEM      0x1000000010000041LL
#define TNG_GMX_SELECTION_GROUP_NAMES   0x1000000010000042LL
#define TNG_GMX_ATOM_SELECTION_GROUP    0x1000000010000043LL
/** @} */

/** Flag to specify if a data block contains data related to particles or not.*/
typedef enum {TNG_NON_PARTICLE_BLOCK_DATA,
              TNG_PARTICLE_BLOCK_DATA} tng_particle_dependency;


typedef enum {TNG_FALSE, TNG_TRUE} tng_bool;

/** Flag to specify if the number of atoms change throughout the trajectory or
 *  if it is constant. */
typedef enum {TNG_CONSTANT_N_ATOMS, TNG_VARIABLE_N_ATOMS}
             tng_variable_n_atoms_flag;

/** Return values of API functions. TNG_SUCCESS means that the operation
 *  was successful. TNG_FAILURE means that the operation failed for some
 *  reason, but it is possible to try to continue anyhow. TNG_CRITICAL
 *  means that the error is irrecoverable. */
typedef enum {TNG_SUCCESS, TNG_FAILURE, TNG_CRITICAL} tng_function_status;

/** If tng_hash_mode == TNG_USE_HASH md5 hashes will be written to output files
 *  and when reading a file the md5 hashes of the contents will be compared to
 *  those in the file (for each block) in order to ensure data integrity */
typedef enum {TNG_SKIP_HASH, TNG_USE_HASH} tng_hash_mode;

/** Possible formats of data block contents */
typedef enum {TNG_CHAR_DATA,
              TNG_INT_DATA,
              TNG_FLOAT_DATA,
              TNG_DOUBLE_DATA} tng_data_type;


struct tng_trajectory;
struct tng_molecule;
struct tng_chain;
struct tng_residue;
struct tng_atom;
struct tng_bond;
struct tng_gen_block;
struct tng_particle_mapping;
struct tng_trajectory_frame_set;
struct tng_particle_data;
struct tng_non_particle_data;

/** Data can be either double, float, int or a string */
union data_values {
    double d;
    float f;
    int64_t i;
    char *c;
};


#ifdef __cplusplus
extern "C"
{
#endif

/** @defgroup group1 Low-level API
 *  These functions give detailed control of the TNG data management. Most
 *  things can be done using the more convenient high-level API functions
 *  instead.
 *  @{
 */

/**
 * @brief Get the major version of the TNG library.
 * @param tng_data is a trajectory data container, it does not have
 * to be initialized beforehand.
 * @param version is pointing to a value set to the major version of
 * the library.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_version_major
                (const tng_trajectory_t tng_data,
                 int *version);

/**
 * @brief Get the minor version of the TNG library.
 * @param tng_data is a trajectory data container, it does not have
 * to be initialized beforehand.
 * @param version is pointing to a value set to the minor version of
 * the library.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_version_minor
                (const tng_trajectory_t tng_data,
                 int *version);

/**
 * @brief Get the patch level of the TNG library.
 * @param tng_data is a trajectory data container, it does not have
 * to be initialized beforehand.
 * @param patch_level is the string to fill with the full version,
 * memory must be allocated before.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_version_patchlevel
                (const tng_trajectory_t tng_data,
                 int *patch_level);

/**
 * @brief Get the full version string of the TNG library.
 * @param tng_data is a trajectory data container, it does not have
 * to be initialized beforehand.
 * @param version is pointing to a value set to the major version of
 * the library.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for version. This includes \0 terminating character.
 * @pre \code version != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_version
                (const tng_trajectory_t tng_data,
                 char *version,
                 const int max_len);

/**
 * @brief Setup a trajectory data container.
 * @param tng_data_p a pointer to memory to initialise as a trajectory.
 * @pre tng_data_p must not be pointing at a reserved memory block.
 * @details Memory is allocated during initialisation.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_trajectory_init
                (tng_trajectory_t *tng_data_p);

/**
 * @brief Clean up a trajectory data container.
 * @param tng_data_p a pointer to the trajectory data to destroy.
 * @details All allocated memory in the data structure is freed, as well as
 * tng_data_p itself.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_trajectory_destroy
                (tng_trajectory_t *tng_data_p);

/**
 * @brief Copy a trajectory data container (dest is setup as well).
 * @details This initialises dest and copies only what is absolute necessary for
 * parallel i/o. This can be used inside pragma omp for setting up a thread
 * local copy of src. It can be freed (using tng_trajectory_destroy) at the
 * end of the parallel block.
 * @param src the original trajectory.
 * @param dest_p a pointer to memory to initialise as a trajectory.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre tng_data_p must not be pointing at a reserved memory block.
 * @details Memory is allocated during initialisation.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_trajectory_init_from_src
                (const tng_trajectory_t src, tng_trajectory_t *dest_p);

/**
 * @brief Get the name of the input file.
 * @param tng_data the trajectory of which to get the input file name.
 * @param file_name the string to fill with the name of the input file,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for file_name. This includes \0 terminating character.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code file_name != 0 \endcode The pointer to the file name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_input_file_get
                (const tng_trajectory_t tng_data,
                 char *file_name, const int max_len);

/**
 * @brief Set the name of the input file.
 * @param tng_data the trajectory of which to set the input file name.
 * @param file_name the name of the input file.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code file_name != 0 \endcode The pointer to the file name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_input_file_set
                (const tng_trajectory_t tng_data,
                 const char *file_name);

/**
 * @brief Get the name of the output file.
 * @param tng_data the trajectory of which to get the input file name.
 * @param file_name the string to fill with the name of the output file,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for file_name. This includes \0 terminating character.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code file_name != 0 \endcode The pointer to the file name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_output_file_get
                (const tng_trajectory_t tng_data,
                 char *file_name, const int max_len);

/**
 * @brief Set the name of the output file.
 * @param tng_data the trajectory of which to set the output file name.
 * @param file_name the name of the output file.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code file_name != 0 \endcode The pointer to the file name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_output_file_set
                (const tng_trajectory_t tng_data,
                 const char *file_name);

/**
 * @brief Set the name of the output file for appending. The output file
 * will not be overwritten.
 * @param tng_data the trajectory of which to set the output file name.
 * @param file_name the name of the output file to append to.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code file_name != 0 \endcode The pointer to the file name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_output_append_file_set
                (const tng_trajectory_t tng_data,
                 const char *file_name);

/**
 * @brief Get the endianness of the output file.
 * @param tng_data the trajectory of which to get the endianness of the current
 * output file.
 * @param endianness will contain the enumeration of the endianness.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code endianness != 0 \endcode The pointer to the endianness container
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (1) if the endianness
 * could not be retrieved.
 */
tng_function_status DECLSPECDLLEXPORT tng_output_file_endianness_get
                (const tng_trajectory_t tng_data, tng_file_endianness *endianness);

/**
 * @brief Set the endianness of the output file.
 * @param tng_data the trajectory of which to set the endianness of the current
 * output file.
 * @param endianness the enumeration of the endianness, can be either
 * TNG_BIG_ENDIAN (0) or TNG_LITTLE_ENDIAN (1).
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @details The endianness cannot be changed after file output has started.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (1) if the endianness
 * could not be set.
 */
tng_function_status DECLSPECDLLEXPORT tng_output_file_endianness_set
                (const tng_trajectory_t tng_data,
                 const tng_file_endianness endianness);

/**
 * @brief Get the name of the program used when creating the trajectory.
 * @param tng_data the trajectory of which to get the program name.
 * @param name the string to fill with the name of the program,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for name. This includes \0 terminating character.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_first_program_name_get
                (const tng_trajectory_t tng_data,
                 char *name, const int max_len);

/**
 * @brief Set the name of the program used when creating the trajectory.
 * @param tng_data the trajectory of which to set the program name.
 * @param new_name is a string containing the wanted name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code new_name != 0 \endcode The pointer to the new_name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_first_program_name_set
                (const tng_trajectory_t tng_data,
                 const char *new_name);

/**
 * @brief Get the name of the program used when last modifying the trajectory.
 * @param tng_data the trajectory of which to get the program name.
 * @param name the string to fill with the name of the program,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for name. This includes \0 terminating character.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_last_program_name_get
                (const tng_trajectory_t tng_data,
                 char *name, const int max_len);

/**
 * @brief Set the name of the program used when last modifying the trajectory.
 * @param tng_data the trajectory of which to set the program name.
 * @param new_name is a string containing the wanted name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code new_name != 0 \endcode The pointer to the new_name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_last_program_name_set
                (const tng_trajectory_t tng_data,
                 const char *new_name);

/**
 * @brief Get the name of the user who created the trajectory.
 * @param tng_data the trajectory of which to get the user name.
 * @param name the string to fill with the name of the user,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for name. This includes \0 terminating character.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_first_user_name_get
                (const tng_trajectory_t tng_data,
                 char *name, const int max_len);

/**
 * @brief Set the name of the user who created the trajectory.
 * @param tng_data the trajectory of which to set the user name.
 * @param new_name is a string containing the wanted name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code new_name != 0 \endcode The pointer to the new_name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_first_user_name_set
                (const tng_trajectory_t tng_data,
                 const char *new_name);

/**
 * @brief Get the name of the user who last modified the trajectory.
 * @param tng_data the trajectory of which to get the user name.
 * @param name the string to fill with the name of the user,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for name. This includes \0 terminating character.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_last_user_name_get
                (const tng_trajectory_t tng_data,
                 char *name, const int max_len);

/**
 * @brief Set the name of the user who last modified the trajectory.
 * @param tng_data the trajectory of which to set the user name.
 * @param new_name is a string containing the wanted name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code new_name != 0 \endcode The pointer to the new_name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_last_user_name_set
                (const tng_trajectory_t tng_data,
                 const char *new_name);

/**
 * @brief Get the name of the computer used when creating the trajectory.
 * @param tng_data the trajectory of which to get the computer name.
 * @param name the string to fill with the name of the computer,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for name. This includes \0 terminating character.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_first_computer_name_get
                (const tng_trajectory_t tng_data,
                 char *name, const int max_len);

/**
 * @brief Set the name of the computer used when creating the trajectory.
 * @param tng_data the trajectory of which to set the computer name.
 * @param new_name is a string containing the wanted name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code new_name != 0 \endcode The pointer to the new_name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_first_computer_name_set
                (const tng_trajectory_t tng_data,
                 const char *new_name);

/**
 * @brief Get the name of the computer used when last modifying the trajectory.
 * @param tng_data the trajectory of which to get the computer name.
 * @param name the string to fill with the name of the computer,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for name. This includes \0 terminating character.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_last_computer_name_get
                (const tng_trajectory_t tng_data,
                 char *name, const int max_len);

/**
 * @brief Set the name of the computer used when last modifying the trajectory.
 * @param tng_data the trajectory of which to set the computer name.
 * @param new_name is a string containing the wanted name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code new_name != 0 \endcode The pointer to the new_name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_last_computer_name_set
                (const tng_trajectory_t tng_data,
                 const char *new_name);

/**
 * @brief Get the pgp_signature of the user creating the trajectory.
 * @param tng_data the trajectory of which to get the computer name.
 * @param signature the string to fill with the signature,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for name. This includes \0 terminating character.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code signature != 0 \endcode The pointer to the signature
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_first_signature_get
                (const tng_trajectory_t tng_data,
                 char *signature, const int max_len);

/**
 * @brief Set the pgp_signature of the user creating the trajectory.
 * @param tng_data the trajectory of which to set the computer name.
 * @param signature is a string containing the pgp_signature.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code signature != 0 \endcode The pointer to the signature
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_first_signature_set
                (const tng_trajectory_t tng_data,
                 const char *signature);

/**
 * @brief Get the pgp_signature of the user last modifying the trajectory.
 * @param tng_data the trajectory of which to get the computer name.
 * @param signature the string to fill with the signature,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for name. This includes \0 terminating character.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code signature != 0 \endcode The pointer to the signature
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_last_signature_get
                (const tng_trajectory_t tng_data,
                 char *signature, const int max_len);

/**
 * @brief Set the pgp_signature of the user last modifying the trajectory.
 * @param tng_data the trajectory of which to set the computer name.
 * @param signature is a string containing the pgp_signature.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code signature != 0 \endcode The pointer to the signature
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_last_signature_set
                (const tng_trajectory_t tng_data,
                 const char *signature);

/**
 * @brief Get the name of the forcefield used in the trajectory.
 * @param tng_data the trajectory of which to get the forcefield name.
 * @param name the string to fill with the name of the forcefield,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for name. This includes \0 terminating character.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_forcefield_name_get
                (const tng_trajectory_t tng_data,
                 char *name, const int max_len);

/**
 * @brief Set the name of the forcefield used in the trajectory.
 * @param tng_data the trajectory of which to set the forcefield name.
 * @param new_name is a string containing the wanted name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code new_name != 0 \endcode The pointer to the new_name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_forcefield_name_set
                (const tng_trajectory_t tng_data,
                 const char *new_name);

/**
 * @brief Get the medium stride length of the trajectory.
 * @param tng_data is the trajectory from which to get the stride length.
 * @param len is pointing to a value set to the stride length.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code len != 0 \endcode The pointer to len must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_medium_stride_length_get
                (const tng_trajectory_t tng_data,
                 int64_t *len);

/**
 * @brief Set the medium stride length of the trajectory.
 * @param tng_data is the trajectory of which to set the stride length.
 * @param len is the wanted medium stride length.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred.
 */
tng_function_status DECLSPECDLLEXPORT tng_medium_stride_length_set
                (const tng_trajectory_t tng_data,
                 const int64_t len);

/**
 * @brief Get the long stride length of the trajectory.
 * @param tng_data is the trajectory from which to get the stride length.
 * @param len is pointing to a value set to the stride length.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code len != 0 \endcode The pointer to len must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_long_stride_length_get
                (const tng_trajectory_t tng_data,
                 int64_t *len);

/**
 * @brief Set the long stride length of the trajectory.
 * @param tng_data is the trajectory of which to set the stride length.
 * @param len is the wanted long stride length.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred.
 */
tng_function_status DECLSPECDLLEXPORT tng_long_stride_length_set
                (const tng_trajectory_t tng_data,
                 const int64_t len);

/**
 * @brief Get the current time per frame of the trajectory.
 * @param tng_data is the trajectory from which to get the time per frame.
 * @param time is pointing to a value set to the time per frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code time != 0 \endcode The pointer to time must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_time_per_frame_get
                (const tng_trajectory_t tng_data,
                 double *time);

/**
 * @brief Set the time per frame of the trajectory.
 * @param tng_data is the trajectory of which to set the time per frame.
 * @param time is the new time per frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code time > 0 \endcode The time per frame must be >= 0.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred.
 */
tng_function_status DECLSPECDLLEXPORT tng_time_per_frame_set
                (const tng_trajectory_t tng_data,
                 const double time);

/**
 * @brief Get the length of the input file.
 * @param tng_data is the trajectory from which to get the input file length.
 * @param len is pointing to a value set to the file length.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code len != 0 \endcode The pointer to len must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_input_file_len_get
                (const tng_trajectory_t tng_data,
                 int64_t *len);

/**
 * @brief Get the number of frames in the trajectory
 * @param tng_data is the trajectory of which to get the number of frames.
 * @param n is pointing to a value set to the number of frames.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code tng_data->input_file != 0 \endcode An input file must be open
 * to find the next frame set.
 * @pre \code n != 0 \endcode The pointer to n must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (could not find last frame set).
 */
tng_function_status DECLSPECDLLEXPORT tng_num_frames_get
                (const tng_trajectory_t tng_data,
                 int64_t *n);

/**
 * @brief Get the precision of lossy compression.
 * @param tng_data is the trajectory of which to get the compression precision.
 * @param precision will be pointing to the retrieved compression precision.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @details A compression precision of 0.001 (the default) means that the
 * compressed values are accurate to the third decimal. This function does
 * not check actual precision of compressed data, but just returns what has
 * previously been set using tng_compression_precision_set().
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_compression_precision_get
                (const tng_trajectory_t tng_data,
                 double *precision);

/**
 * @brief Set the precision of lossy compression.
 * @param tng_data is the trajectory of which to set the compression precision.
 * @param precision is the new compression precision.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @details A compression precision of 0.001 (the default) means that the
 * compressed values are accurate to the third decimal.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_compression_precision_set
                (const tng_trajectory_t tng_data,
                 const double precision);

/**
 * @brief Set the number of particles, in the case no molecular system is used.
 * @param tng_data is the trajectory of which to get the number of particles.
 * @param n is the number of particles to use.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @details When creating a molecular system the number of particles are set
 * automatically. This should only be used when there is no molecular system
 * specified or if the number of atoms needs to be overridden for some reason.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_implicit_num_particles_set
                (const tng_trajectory_t tng_data,
                 const int64_t n);

/**
 * @brief Get the current number of particles.
 * @param tng_data is the trajectory from which to get the number of particles.
 * @param n is pointing to a value set to the number of particles.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n != 0 \endcode The pointer to n must not be a NULL pointer.
 * @details If variable number of particles are used this function will return
 * the number of particles in the current frame set.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_num_particles_get
                (const tng_trajectory_t tng_data,
                 int64_t *n);

/**
 * @brief Get if the number of particle can be varied during the simulation.
 * @param tng_data is the trajectory from which to get the number of particles.
 * @param variable is pointing to a value set to TNG_CONSTANT_N_ATOMS if the
 * number of particles cannot change or TNG_VARIABLE_N_ATOMS if the number of
 * particles can change.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code variable != 0 \endcode The pointer to variable must not be
 * a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_num_particles_variable_get
                (const tng_trajectory_t tng_data,
                 char *variable);

/**
 * @brief Get the number of molecule types (length of tng_data->molecules).
 * @param tng_data is the trajectory from which to get the number of molecules.
 * @param n is pointing to a value set to the number of molecule types.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n != 0 \endcode The pointer to n must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_num_molecule_types_get
                (const tng_trajectory_t tng_data,
                 int64_t *n);

/**
 * @brief Get the current total number of molecules.
 * @param tng_data is the trajectory from which to get the number of molecules.
 * @param n is pointing to a value set to the number of molecules.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n != 0 \endcode The pointer to n must not be a NULL pointer.
 * @details If variable number of particles are used this function will return
 * the total number of molecules in the current frame set.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_num_molecules_get
                (const tng_trajectory_t tng_data,
                 int64_t *n);

/** @brief Get the list of the count of each molecule.
 * @param tng_data is the trajectory from which to get the molecule count list.
 * @param mol_cnt_list is a list of the count of each molecule in the
 * mol system. This is a pointer to the list in the TNG container, which
 * means that it should be handled carefully, e.g. not freed.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE(1) if the list of
 * molecule counts was not valid.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_cnt_list_get
                (const tng_trajectory_t tng_data,
                 int64_t **mol_cnt_list);

/**
 * @brief Get the exponent used for distances in the trajectory.
 * @param tng_data is the trajectory from which to get the information.
 * @param exp is pointing to a value set to the distance unit exponent.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code exp != 0 \endcode The pointer to exp must not be a NULL pointer.
 * @details Example: If the distances are specified in nm (default) exp is -9.
 * If the distances are specified in Å exp is -10.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_distance_unit_exponential_get
                (const tng_trajectory_t tng_data,
                 int64_t *exp);

/**
 * @brief Set the exponent used for distances in the trajectory.
 * @param tng_data is the trajectory of which to set the unit exponent.
 * @param exp is the distance unit exponent to use.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @details Example: If the distances are specified in nm (default) exp is -9.
 * If the distances are specified in Å exp is -10.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_distance_unit_exponential_set
                (const tng_trajectory_t tng_data,
                 const int64_t exp);

/**
 * @brief Get the number of frames per frame set.
 * @param tng_data is the trajectory from which to get the number of frames
 * per frame set.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n != 0 \endcode The pointer to n must not be a NULL pointer.
 * @param n is pointing to a value set to the number of frames per frame set.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_num_frames_per_frame_set_get
                (const tng_trajectory_t tng_data,
                 int64_t *n);

/**
 * @brief Set the number of frames per frame set.
 * @param tng_data is the trajectory of which to set the number of frames
 * per frame set.
 * @param n is the number of frames per frame set.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @details This does not affect already existing frame sets. For
 * consistency the number of frames per frame set should be set
 * betfore creating any frame sets.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_num_frames_per_frame_set_set
                (const tng_trajectory_t tng_data,
                 const int64_t n);

/**
 * @brief Get the number of frame sets.
 * @details This updates tng_data->n_trajectory_frame_sets before returning it.
 * @param tng_data is the trajectory from which to get the number of frame sets.
 * @param n is pointing to a value set to the number of frame sets.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n != 0 \endcode The pointer to n must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_num_frame_sets_get
                (const tng_trajectory_t tng_data,
                 int64_t *n);

/**
 * @brief Get the current trajectory frame set.
 * @param tng_data is the trajectory from which to get the frame set.
 * @param frame_set_p will be set to point at the memory position of
 * the found frame set.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_current_frame_set_get
                (const tng_trajectory_t tng_data,
                 tng_trajectory_frame_set_t *frame_set_p);

/**
 * @brief Find the requested frame set number.
 * @param tng_data is the trajectory from which to get the frame set.
 * @param nr is the frame set number to search for.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code nr >= 0 \endcode The frame set number (nr) must be >= 0.
 * @details tng_data->current_trajectory_frame_set will contain the
 * found trajectory if successful.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_nr_find
                (const tng_trajectory_t tng_data,
                 const int64_t nr);

/**
 * @brief Find the frame set containing a specific frame.
 * @param tng_data is the trajectory from which to get the frame set.
 * @param frame is the frame number to search for.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code frame >= 0 \endcode The frame number must be >= 0.
 * @details tng_data->current_trajectory_frame_set will contain the
 * found trajectory if successful.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_of_frame_find
                (const tng_trajectory_t tng_data,
                 const int64_t frame);

/**
 * @brief Get the file position of the next frame set in the input file.
 * @param tng_data is a trajectory data container.
 * @param frame_set is the frame set of which to get the position of the
 * following frame set.
 * @param pos is pointing to a value set to the file position.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code pos != 0 \endcode The pointer to pos must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_next_frame_set_file_pos_get
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *pos);

/**
 * @brief Get the file position of the previous frame set in the input file.
 * @param tng_data is a trajectory data container.
 * @param frame_set is the frame set of which to get the position of the
 * previous frame set.
 * @param pos is pointing to a value set to the file position.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code pos != 0 \endcode The pointer to pos must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_prev_frame_set_file_pos_get
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *pos);

/**
 * @brief Get the first and last frames of the frame set.
 * @param tng_data is a trajectory data container.
 * @param frame_set is the frame set of which to get the frame range.
 * @param first_frame is set to the first frame of the frame set.
 * @param last_frame is set to the last frame of the frame set.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code first_frame != 0 \endcode The pointer to first_frame must
 * not be a NULL pointer.
 * @pre \code last_frame != 0 \endcode The pointer to last_frame must
 * not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_frame_range_get
                (const tng_trajectory_t tng_data,
                 const tng_trajectory_frame_set_t frame_set,
                 int64_t *first_frame,
                 int64_t *last_frame);

/**
 * @brief Allocate memory for and setup a molecule container.
 * @param tng_data is a trajectory data container.
 * @param molecule_p is a pointer to molecule to allocate and initialise.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_alloc(const tng_trajectory_t tng_data,
                                                         tng_molecule_t *molecule_p);

/**
 * @brief Clean up a molecule container and free its allocated memory.
 * @param tng_data is a trajectory data container.
 * @param molecule_p is the molecule to destroy.
 * @details All allocated memory in the data structure is freed and also the memory
 * of the molecule itself.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_free(const tng_trajectory_t tng_data,
                                                        tng_molecule_t *molecule_p);

/**
 * @brief Setup a molecule container.
 * @param tng_data is a trajectory data container.
 * @param molecule is the molecule to initialise. Memory must be preallocated.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_init
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule);

/**
 * @brief Clean up a molecule container.
 * @param tng_data is a trajectory data container.
 * @param molecule is the molecule to destroy.
 * @details All allocated memory in the data structure is freed, but not the
 * memory of molecule itself.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_destroy
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule);

/**
 * @brief Add a molecule to the trajectory.
 * @param tng_data is the trajectory data container containing the block..
 * @param name is a pointer to the string containing the name of the new molecule.
 * @param molecule is a pointer to the newly created molecule.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the ID could
 * not be set properly or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_add
                (const tng_trajectory_t tng_data,
                 const char *name,
                 tng_molecule_t *molecule);

/**
 * @brief Add a molecule with a specific ID to the trajectory.
 * @param tng_data is the trajectory data container containing the block..
 * @param name is a pointer to the string containing the name of the new molecule.
 * @param id is the ID of the created molecule.
 * @param molecule is a pointer to the newly created molecule.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the ID could
 * not be set properly or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_w_id_add
                (const tng_trajectory_t tng_data,
                 const char *name,
                 const int64_t id,
                 tng_molecule_t *molecule);

/**
 * @brief Add an existing molecule (from a molecule container) to the trajectory.
 * @param tng_data is the trajectory data container containing the block..
 * @param molecule is a pointer to the molecule to add to the trajectory and will
 * afterwards point to the molecule in the trajectory.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_existing_add
                (const tng_trajectory_t tng_data,
                 tng_molecule_t *molecule);

/**
 * @brief Get the name of a molecule.
 * @param tng_data the trajectory containing the molecule.
 * @param molecule the molecule of which to get the name.
 * @param name the string to fill with the name of the molecule,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for name. This includes \0 terminating character.
 * @pre \code molecule != 0 \endcode The molecule must not be NULL.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_name_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 char *name,
                 const int max_len);

/**
 * @brief Set the name of a molecule.
 * @param tng_data is the trajectory data container containing the molecule..
 * @param molecule is the molecule to rename.
 * @param new_name is a string containing the wanted name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code new_name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_name_set
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 const char *new_name);

/**
 * @brief Get the count of a molecule.
 * @param tng_data is the trajectory data container containing the molecule..
 * @param molecule is the molecule of which to get the count.
 * @param cnt is a pointer to the variable to be populated with the count.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code cnt != 0 \endcode The pointer to the molecule count
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_cnt_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 int64_t *cnt);

/**
 * @brief Set the count of a molecule.
 * @param tng_data is the trajectory data container containing the molecule..
 * @param molecule is the molecule of which to set the count.
 * @param cnt is the number of instances of this molecule.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_cnt_set
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 const int64_t cnt);

/**
 * @brief Find a molecule.
 * @param tng_data is the trajectory data container containing the molecule.
 * @param name is a string containing the name of the molecule. If name is empty
 * only id will be used for finding the molecule.
 * @param id is the id of the molecule to look for. If id is -1 only the name of
 * the molecule will be used for finding the molecule.
 * @param molecule is a pointer to the molecule if it was found - otherwise 0.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if the molecule is found or TNG_FAILURE (1) if the
 * molecule is not found.
 * @details If name is an empty string and id == -1 the first residue will
 * be found.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_find
                (const tng_trajectory_t tng_data,
                 const char *name,
                 const int64_t id,
                 tng_molecule_t *molecule);

/**
 * @brief Retrieve the molecule with specified index in the list of molecules.
 * @param tng_data is the trajectory data container containing the molecule.
 * @param index is the index (in tng_data->molecules) of the molecule to return
 * @param molecule is a pointer to the molecule if it was found - otherwise 0.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code molecule != 0 \endcode molecule must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if the molecule is found or TNG_FAILURE (1) if the
 * molecule is not found.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_of_index_get
                (const tng_trajectory_t tng_data,
                 const int64_t index,
                 tng_molecule_t *molecule);

/**
 * @brief Copy all molecules and the molecule counts from one TNG trajectory
 * to another.
 * @param tng_data_src is the source trajectory containing the molecular
 * system to copy.
 * @param tng_data_dest is the destination trajectory.
 * @pre \code tng_data_src != 0 \endcode The trajectory container (tng_data_src)
 * must be initialised before using it.
 * @pre \code tng_data_dest != 0 \endcode The trajectory container (tng_data_dest)
 * must be initialised before using it.
 * @details The molecular system in tng_data_dest will be overwritten.
 * @return TNG_SUCCESS(0) if the copying is successful, TNG_FAILURE if a minor
 * error has occured or TNG_CRITICAL(2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_system_copy(const tng_trajectory_t tng_data_src,
                                                               const tng_trajectory_t tng_data_dest);

/**
 * @brief Get the number of chains in a molecule.
 * @param tng_data is the trajectory containing the molecule.
 * @param molecule is the molecule of which to get the number of chains.
 * @param n is pointing to a value set to the number of chains.
 * @pre \code molecule != 0 \endcode The molecule must not be NULL.
 * @pre \code n != 0 \endcode The pointer to n must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_num_chains_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 int64_t *n);

/**
 * @brief Retrieve the chain of a molecule with specified index in the list
 * of chains.
 * @param tng_data is the trajectory data container containing the molecule.
 * @param index is the index (in molecule->chains) of the chain to return
 * @param molecule is the molecule from which to get the chain.
 * @param chain is a pointer to the chain if it was found - otherwise 0.
 * @pre \code molecule != 0 \endcode molecule must not be a NULL pointer.
 * @pre \code chain != 0 \endcode chain must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if the chain is found or TNG_FAILURE (1) if the
 * chain is not found.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_chain_of_index_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 const int64_t index,
                 tng_chain_t *chain);

/**
 * @brief Get the number of residues in a molecule.
 * @param tng_data is the trajectory containing the molecule.
 * @param molecule is the molecule of which to get the number residues.
 * @param n is pointing to a value set to the number of residues.
 * @pre \code molecule != 0 \endcode The molecule must not be NULL.
 * @pre \code n != 0 \endcode The pointer to n must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_num_residues_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 int64_t *n);

/**
 * @brief Retrieve the residue of a molecule with specified index in the list
 * of chains.
 * @param tng_data is the trajectory data container containing the molecule.
 * @param index is the index (in molecule->residues) of the residue to return
 * @param molecule is the molecule from which to get the residue.
 * @param residue is a pointer to the residue if it was found - otherwise 0.
 * @pre \code molecule != 0 \endcode molecule must not be a NULL pointer.
 * @pre \code residue != 0 \endcode residue must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if the residue is found or TNG_FAILURE (1) if the
 * residue is not found.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_residue_of_index_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 const int64_t index,
                 tng_residue_t *residue);

/**
 * @brief Get the number of atoms in a molecule.
 * @param tng_data is the trajectory containing the molecule.
 * @param molecule is the molecule of which to get the number of atoms.
 * @param n is pointing to a value set to the number of atoms.
 * @pre \code molecule != 0 \endcode The molecule must not be NULL.
 * @pre \code n != 0 \endcode The pointer to n must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_num_atoms_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 int64_t *n);

/**
 * @brief Retrieve the atom of a molecule with specified index in the list
 * of atoms.
 * @param tng_data is the trajectory data container containing the molecule.
 * @param index is the index (in molecule->atoms) of the atom to return
 * @param molecule is the molecule from which to get the atom.
 * @param atom is a pointer to the atom if it was found - otherwise 0.
 * @pre \code molecule != 0 \endcode molecule must not be a NULL pointer.
 * @pre \code atom != 0 \endcode atom must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if the atom is found or TNG_FAILURE (1) if the
 * atom is not found.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_atom_of_index_get
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 const int64_t index,
                 tng_atom_t *atom);

/**
 * @brief Find a chain in a molecule.
 * @param tng_data is the trajectory data container containing the molecule.
 * @param molecule is the molecule in which to search for the chain.
 * @param name is a string containing the name of the chain. If name is empty
 * only id will be used for finding the chain.
 * @param id is the id of the chain to look for. If id is -1 only the name of
 * the chain will be used for finding the chain.
 * @param chain is a pointer to the chain if it was found - otherwise 0.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if the chain is found or TNG_FAILURE (1) if the
 * chain is not found.
 * @details If name is an empty string and id == -1 the first residue will
 * be found.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_chain_find
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 const char *name,
                 const int64_t id,
                 tng_chain_t *chain);

/**
 * @brief Add a chain to a molecule.
 * @param tng_data is the trajectory data container containing the molecule..
 * @param molecule is the molecule to add a chain to.
 * @param name is a string containing the name of the chain.
 * @param chain is a pointer to the newly created chain.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the ID could
 * not be set properly or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_chain_add
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 const char *name,
                 tng_chain_t *chain);

/**
 * @brief Add a chain with a specific id to a molecule.
 * @param tng_data is the trajectory data container containing the molecule..
 * @param molecule is the molecule to add a chain to.
 * @param name is a string containing the name of the chain.
 * @param id is the ID of the created chain.
 * @param chain is a pointer to the newly created chain.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the ID could
 * not be set properly or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_chain_w_id_add
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 const char *name,
                 const int64_t id,
                 tng_chain_t *chain);

/**
 * @brief Add a bond between two atoms to a molecule.
 * @param tng_data is the trajectory data container containing the molecule.
 * @param molecule is the molecule containing the atoms to connect.
 * @param from_atom_id is the id of one of the two atoms in the bond.
 * @param to_atom_id is the id of the other atom in the bond.
 * @param bond is a pointer to the newly created bond.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (!) if a minor error
 * has occured or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_bond_add
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 const int64_t from_atom_id,
                 const int64_t to_atom_id,
                 tng_bond_t *bond);

/**
 * @brief Find an atom in a molecule.
 * @param tng_data is the trajectory data container containing the molecule.
 * @param molecule is the molecule in which to search for the atom.
 * @param name is a string containing the name of the atom. If name is an
 * empty string only id will be used for searching.
 * @param id is the id of the atom to find. If id == -1 the first atom
 * that matches the specified name will be found.
 * @param atom is a pointer to the atom if it was found - otherwise 0.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if the atom is found or TNG_FAILURE (1) if the
 * atom is not found.
 * @details If name is an empty string and id == -1 the first residue will
 * be found.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_atom_find
                (const tng_trajectory_t tng_data,
                 const tng_molecule_t molecule,
                 const char *name,
                 const int64_t id,
                 tng_atom_t *atom);

/**
 * @brief Get the name of a chain.
 * @param tng_data the trajectory containing the chain.
 * @param chain the chain of which to get the name.
 * @param name the string to fill with the name of the chain,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for name. This includes \0 terminating character.
 * @pre \code chain != 0 \endcode The chain must not be NULL.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_chain_name_get
                (const tng_trajectory_t tng_data,
                 const tng_chain_t chain,
                 char *name,
                 const int max_len);

/**
 * @brief Set the name of a chain.
 * @param tng_data is the trajectory data container containing the atom..
 * @param chain is the chain to rename.
 * @param new_name is a string containing the wanted name.
 * @pre \code new_name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_chain_name_set
                (const tng_trajectory_t tng_data,
                 const tng_chain_t chain,
                 const char *new_name);

/**
 * @brief Get the number of residues in a molecule chain.
 * @param tng_data is the trajectory containing the chain.
 * @param chain is the chain of which to get the number of residues.
 * @param n is pointing to a value set to the number of residues.
 * @pre \code chain != 0 \endcode The chain must not be NULL.
 * @pre \code n != 0 \endcode The pointer to n must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_chain_num_residues_get
                (const tng_trajectory_t tng_data,
                 const tng_chain_t chain,
                 int64_t *n);

/**
 * @brief Retrieve the residue of a chain with specified index in the list
 * of residues.
 * @param tng_data is the trajectory data container containing the chain.
 * @param index is the index (in chain->residues) of the residue to return
 * @param chain is the chain from which to get the residue.
 * @param residue is a pointer to the residue if it was found - otherwise 0.
 * @pre \code chain != 0 \endcode chain must not be a NULL pointer.
 * @pre \code residue != 0 \endcode residue must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if the residue is found or TNG_FAILURE (1) if the
 * residue is not found.
 */
tng_function_status DECLSPECDLLEXPORT tng_chain_residue_of_index_get
                (const tng_trajectory_t tng_data,
                 const tng_chain_t chain,
                 const int64_t index,
                 tng_residue_t *residue);

/**
 * @brief Find a residue in a chain.
 * @param tng_data is the trajectory data container containing the chain.
 * @param chain is the chain in which to search for the residue.
 * @param name is a string containing the name of the residue.  If name is an
 * empty string only id will be used for searching.
 * @param id is the id of the residue to find. If id == -1 the first residue
 * that matches the specified name will be found.
 * @param residue is a pointer to the residue if it was found - otherwise 0.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if the residue is found or TNG_FAILURE (1) if the
 * residue is not found.
 * @details If name is an empty string and id == -1 the first residue will
 * be found.
 */
tng_function_status DECLSPECDLLEXPORT tng_chain_residue_find
                (const tng_trajectory_t tng_data,
                 const tng_chain_t chain,
                 const char *name,
                 const int64_t id,
                 tng_residue_t *residue);

/**
 * @brief Add a residue to a chain.
 * @param tng_data is the trajectory data container containing the chain..
 * @param chain is the chain to add a residue to.
 * @param name is a string containing the name of the residue.
 * @param residue is a pointer to the newly created residue.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the ID could
 * not be set properly or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_chain_residue_add
                (const tng_trajectory_t tng_data,
                 const tng_chain_t chain,
                 const char *name,
                 tng_residue_t *residue);

/**
 * @brief Add a residue with a specific ID to a chain.
 * @param tng_data is the trajectory data container containing the chain..
 * @param chain is the chain to add a residue to.
 * @param name is a string containing the name of the residue.
 * @param id is the ID of the created residue.
 * @param residue is a pointer to the newly created residue.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the ID could
 * not be set properly or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_chain_residue_w_id_add
                (const tng_trajectory_t tng_data,
                 const tng_chain_t chain,
                 const char *name,
                 const int64_t id,
                 tng_residue_t *residue);

/**
 * @brief Get the name of a residue.
 * @param tng_data the trajectory containing the residue.
 * @param residue the residue of which to get the name.
 * @param name the string to fill with the name of the residue,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for name. This includes \0 terminating character.
 * @pre \code residue != 0 \endcode The residue must not be NULL.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_residue_name_get
                (const tng_trajectory_t tng_data,
                 const tng_residue_t residue,
                 char *name,
                 const int max_len);

/**
 * @brief Set the name of a residue.
 * @param tng_data is the trajectory data container containing the residue.
 * @param residue is the residue to rename.
 * @param new_name is a string containing the wanted name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code new_name != 0 \endcode The new name to set (new_name) must
 * not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_residue_name_set
                (const tng_trajectory_t tng_data,
                 const tng_residue_t residue,
                 const char *new_name);

/**
 * @brief Get the number of atoms in a residue.
 * @param tng_data is the trajectory containing the residue.
 * @param residue is the residue of which to get the number atoms.
 * @param n is pointing to a value set to the number of atoms.
 * @pre \code residue != 0 \endcode The residue must not be NULL.
 * @pre \code n != 0 \endcode The pointer to n must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_residue_num_atoms_get
                (const tng_trajectory_t tng_data,
                 const tng_residue_t residue,
                 int64_t *n);

/**
 * @brief Retrieve the atom of a residue with specified index in the list
 * of atoms.
 * @param tng_data is the trajectory data container containing the residue.
 * @param index is the index (in residue->atoms) of the atom to return
 * @param residue is the residue from which to get the atom.
 * @param atom is a pointer to the atom if it was found - otherwise 0.
 * @pre \code residue != 0 \endcode residue must not be a NULL pointer.
 * @pre \code atom != 0 \endcode atom must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if the atom is found or TNG_FAILURE (1) if the
 * atom is not found.
 */
tng_function_status DECLSPECDLLEXPORT tng_residue_atom_of_index_get
                (const tng_trajectory_t tng_data,
                 const tng_residue_t residue,
                 const int64_t index,
                 tng_atom_t *atom);

/**
 * @brief Add an atom to a residue.
 * @param tng_data is the trajectory containing the residue.
 * @param residue is the residue to add an atom to.
 * @param atom_name is a string containing the name of the atom.
 * @param atom_type is a string containing the atom type of the atom.
 * @param atom is a pointer to the newly created atom.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code atom_name != 0 \endcode The pointer to the atom name string
 * must not be a NULL pointer.
 * @pre \code atom_type != 0 \endcode The pointer to the atom_type string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the ID could
 * not be set properly or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_residue_atom_add
                (const tng_trajectory_t tng_data,
                 const tng_residue_t residue,
                 const char *atom_name,
                 const char *atom_type,
                 tng_atom_t *atom);

/**
 * @brief Add an atom with a specific ID to a residue.
 * @param tng_data is the trajectory containing the residue.
 * @param residue is the residue to add an atom to.
 * @param atom_name is a string containing the name of the atom.
 * @param atom_type is a string containing the atom type of the atom.
 * @param id is the ID of the created atom.
 * @param atom is a pointer to the newly created atom.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code atom_name != 0 \endcode The pointer to the atom name string
 * must not be a NULL pointer.
 * @pre \code atom_type != 0 \endcode The pointer to the atom_type string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if the ID could
 * not be set properly or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_residue_atom_w_id_add
                (const tng_trajectory_t tng_data,
                 const tng_residue_t residue,
                 const char *atom_name,
                 const char *atom_type,
                 const int64_t id,
                 tng_atom_t *atom);

/**
 * @brief Get the residue of an atom.
 * @param tng_data the trajectory containing the atom.
 * @param atom the atom of which to get the name.
 * @param residue is set to the residue of the atom.
 * @pre \code atom != 0 \endcode The atom must not be NULL.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_atom_residue_get
                (const tng_trajectory_t tng_data,
                 const tng_atom_t atom,
                 tng_residue_t *residue);

/**
 * @brief Get the name of an atom.
 * @param tng_data the trajectory containing the atom.
 * @param atom the atom of which to get the name.
 * @param name the string to fill with the name of the atom,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for name. This includes \0 terminating character.
 * @pre \code atom != 0 \endcode The atom must not be NULL.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_atom_name_get
                (const tng_trajectory_t tng_data,
                 const tng_atom_t atom,
                 char *name,
                 const int max_len);

/**
 * @brief Set the name of an atom.
 * @param tng_data is the trajectory data container containing the atom.
 * @param atom is the atom to rename.
 * @param new_name is a string containing the wanted name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code new_name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_atom_name_set
                (const tng_trajectory_t tng_data,
                 const tng_atom_t atom,
                 const char *new_name);

/**
 * @brief Get the type of an atom.
 * @param tng_data the trajectory containing the atom.
 * @param atom the atom of which to get the type.
 * @param type the string to fill with the type of the atom,
 * memory must be allocated before.
 * @param max_len maximum char length of the string, i.e. how much memory has
 * been reserved for type. This includes \0 terminating character.
 * @pre \code atom != 0 \endcode The atom must not be NULL.
 * @pre \code type != 0 \endcode The pointer to the type string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred (source string longer than destination string).
 */
tng_function_status DECLSPECDLLEXPORT tng_atom_type_get
                (const tng_trajectory_t tng_data,
                 const tng_atom_t atom,
                 char *type,
                 const int max_len);

/**
 * @brief Set the atom type of an atom.
 * @param tng_data is the trajectory data container containing the atom.
 * @param atom is the atom to change.
 * @param new_type is a string containing the atom type.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code new_type != 0 \endcode The pointer to the atom type string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_atom_type_set
                (const tng_trajectory_t tng_data,
                 const tng_atom_t atom,
                 const char *new_type);

/**
 * @brief Get the molecule name of real particle number (number in mol system).
 * @param tng_data is the trajectory data container containing the atom.
 * @param nr is the real number of the particle in the molecular system.
 * @param name is a string, which is set to the name of the molecule. Memory
 * must be reserved beforehand.
 * @param max_len is the maximum length of name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (!) if a minor error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_name_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 const int max_len);

/**
 * @brief Get the molecule id of real particle number (number in mol system).
 * @param tng_data is the trajectory data container containing the atom.
 * @param nr is the real number of the particle in the molecular system.
 * @param id is will be set to the id of the molecule.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code id != 0 \endcode The pointer to id must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (!) if a minor error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molecule_id_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 int64_t *id);

/**
 * @brief Get the bonds of the current molecular system.
 * @param tng_data is the trajectory data container containing the molecular
 * system.
 * @param n_bonds is set to the number of bonds in the molecular system and
 * thereby also the lengths of the two lists: from_atoms and to_atoms.
 * @param from_atoms is a list (memory reserved by this function) of atoms
 * (number of atom in mol system) in bonds.
 * @param to_atoms is a list (memory reserved by this function) of atoms
 * (number of atom in mol system) in bonds.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n_bonds != 0 \endcode The pointer to n_bonds must not be a
 * NULL pointer.
 * @pre \code from_atoms != 0 \endcode The pointer to from_atoms must not
 * be a NULL pointer.
 * @pre \code to_atoms != 0 \endcode The pointer to to_atoms must not
 * be a NULL pointer.
 * @details The two lists of atoms use the same index, i.e. from_atoms[0]
 * and to_atoms[0] are linked with a bond. Since memory is reserved in
 * this function it must be freed afterwards.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_molsystem_bonds_get
                (const tng_trajectory_t tng_data,
                 int64_t *n_bonds,
                 int64_t **from_atoms,
                 int64_t **to_atoms);

/**
 * @brief Get the chain name of real particle number (number in mol system).
 * @param tng_data is the trajectory data container containing the atom.
 * @param nr is the real number of the particle in the molecular system.
 * @param name is a string, which is set to the name of the chain. Memory
 * must be reserved beforehand.
 * @param max_len is the maximum length of name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (!) if a minor error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_chain_name_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 const int max_len);

/**
 * @brief Get the residue name of real particle number (number in mol system).
 * @param tng_data is the trajectory data container containing the atom.
 * @param nr is the real number of the particle in the molecular system.
 * @param name is a string, which is set to the name of the residue. Memory
 * must be reserved beforehand.
 * @param max_len is the maximum length of name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (!) if a minor error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_residue_name_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 const int max_len);

/**
 * @brief Get the residue id (local to molecule) of real particle number
 * (number in mol system).
 * @param tng_data is the trajectory data container containing the atom.
 * @param nr is the real number of the particle in the molecular system.
 * @param id is a pointer to the variable, which will be set to the ID.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code id != 0 \endcode The pointer to id must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (!) if a minor error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_residue_id_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 int64_t *id);

/**
 * @brief Get the residue id (based on other molecules and molecule counts)
 * of real particle number (number in mol system).
 * @param tng_data is the trajectory data container containing the atom.
 * @param nr is the real number of the particle in the molecular system.
 * @param id is a pointer to the variable, which will be set to the ID.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code id != 0 \endcode The pointer to id must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (!) if a minor error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_global_residue_id_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 int64_t *id);

/**
 * @brief Get the atom name of real particle number (number in mol system).
 * @param tng_data is the trajectory data container containing the atom.
 * @param nr is the real number of the particle in the molecular system.
 * @param name is a string, which is set to the name of the atom. Memory
 * must be reserved beforehand.
 * @param max_len is the maximum length of name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (!) if a minor error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_atom_name_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *name,
                 const int max_len);

/**
 * @brief Get the atom type of real particle number (number in mol system).
 * @param tng_data is the trajectory data container containing the atom.
 * @param nr is the real number of the particle in the molecular system.
 * @param type is a string, which is set to the type of the atom. Memory
 * must be reserved beforehand.
 * @param max_len is the maximum length of type.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code type != 0 \endcode The pointer to the type string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (!) if a minor error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_atom_type_of_particle_nr_get
                (const tng_trajectory_t tng_data,
                 const int64_t nr,
                 char *type,
                 const int max_len);

/**
 * @brief Add a particle mapping table.
 * @details Each particle mapping table will be written as a separate block,
 * followed by the data blocks for the corresponding particles. In most cases
 * there is one particle mapping block for each thread writing the trajectory.
 * @param tng_data is the trajectory, with the frame set to which to add
 * the mapping block.
 * @details The mapping information is added to the currently active frame set
 * of tng_data
 * @param num_first_particle is the first particle number of this mapping
 * block.
 * @param n_particles is the number of particles in this mapping block.
 * @param mapping_table is a list of the real particle numbers (i.e. the numbers
 * used in the molecular system). The list is n_particles long.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @details mapping_table[0] is the real particle number of the first particle
 * in the following data blocks.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_particle_mapping_add
                (const tng_trajectory_t tng_data,
                 const int64_t num_first_particle,
                 const int64_t n_particles,
                 const int64_t *mapping_table);

/**
 * @brief Remove all particle mappings (in memory) from the current frame set.
 * @details Clears the currently setup particle mappings of the current frame
 * set.
 * @param tng_data is the trajectory, with the frame set of which to clear
 * all particle mappings.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_particle_mapping_free
                (const tng_trajectory_t tng_data);

/**
 * @brief Read the header blocks from the input_file of tng_data.
 * @details The trajectory blocks must be read separately and iteratively in chunks
 * to fit in memory.
 * @param tng_data is a trajectory data container.
 * @details tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_file_headers_read
                (const tng_trajectory_t tng_data,
                 const char hash_mode);

/**
 * @brief Write the header blocks to the output_file of tng_data.
 * @details The trajectory blocks must be written separately and iteratively in chunks
 * to fit in memory.
 * @param tng_data is a trajectory data container.
 * @details tng_data->output_file_path
 * specifies which file to write to. If the file (output_file) is not open it
 * will be opened.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH an md5 hash for each header block will be generated.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_file_headers_write
                (const tng_trajectory_t tng_data,
                 const char hash_mode);

/**
 * @brief Read one (the next) block (of any kind) from the input_file of tng_data.
 * @param tng_data is a trajectory data container.
 * @details tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param block_data is a pointer to the struct which will be populated with the
 * data.
 * @details If block_data->input_file_pos > 0 it is the position from where the
 * reading starts otherwise it starts from the current position.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code block != 0 \endcode The block container (block) must be
 * initialised before using it.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_block_read_next
                (const tng_trajectory_t tng_data,
                 const tng_gen_block_t block_data,
                 const char hash_mode);

/**
 * @brief Read one frame set, including all particle mapping blocks and data
 * blocks, starting from the current file position.
 * @param tng_data is a trajectory data container.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_read
                (const tng_trajectory_t tng_data,
                 const char hash_mode);

/**
 * @brief Read data from the current frame set from the input_file. Only read
 * particle mapping and data blocks matching the specified block_id.
 * @param tng_data is a trajectory data container.
 * @details  tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @param block_id is the ID of the data block to read from file.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_read_current_only_data_from_block_id
                (const tng_trajectory_t tng_data,
                 const char hash_mode,
                 const int64_t block_id);

/**
 * @brief Read one (the next) frame set, including particle mapping and related data blocks
 * from the input_file of tng_data.
 * @param tng_data is a trajectory data container.
 * @details  tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_read_next
                (const tng_trajectory_t tng_data,
                 const char hash_mode);

/**
 * @brief Read one (the next) frame set, including particle mapping and data blocks with a
 * specific block id from the input_file of tng_data.
 * @param tng_data is a trajectory data container.
 * @details  tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @param block_id is the ID number of the blocks that should be read from file.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_read_next_only_data_from_block_id
                (const tng_trajectory_t tng_data,
                 const char hash_mode,
                 const int64_t block_id);

/**
 * @brief Write one frame set, including mapping and related data blocks
 * to the output_file of tng_data.
 * @param tng_data is a trajectory data container.
 * @details  tng_data->output_file_path specifies
 * which file to write to. If the file (output_file) is not open it will be
 * opened.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH an md5 hash for each header block will be generated.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_write
                (const tng_trajectory_t tng_data,
                 const char hash_mode);

/**
 * @brief Write one frame set even if it does not have as many frames as
 * expected. The function also writes mapping and related data blocks
 * to the output_file of tng_data.
 * @param tng_data is a trajectory data container.
 * @details  tng_data->output_file_path specifies
 * which file to write to. If the file (output_file) is not open it will be
 * opened.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH an md5 hash for each header block will be generated.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @details The number of frames in the frame set is set to the number of
 * frames of the data blocks before writing it to disk.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_premature_write
                (const tng_trajectory_t tng_data,
                 const char hash_mode);

/**
 * @brief Create and initialise a frame set.
 * @details Particle mappings are retained from previous frame set (if any).
 * To explicitly clear particle mappings use tng_frame_set_particle_mapping_free().
 * @param tng_data is the trajectory data container in which to add the frame
 * set.
 * @param first_frame is the first frame of the frame set.
 * @param n_frames is the number of frames in the frame set.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code first_frame >= 0 \endcode The first frame must not be negative.
 * @pre \code n_frames >= 0 \endcode The number of frames must not be negative.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_new
                (const tng_trajectory_t tng_data,
                 const int64_t first_frame,
                 const int64_t n_frames);

/**
 * @brief Create and initialise a frame set with the time of the first frame
 * specified.
 * @param tng_data is the trajectory data container in which to add the frame
 * set.
 * @param first_frame is the first frame of the frame set.
 * @param n_frames is the number of frames in the frame set.
 * @param first_frame_time is the time stamp of the first frame (in seconds).
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code first_frame >= 0 \endcode The first frame must not be negative.
 * @pre \code n_frames >= 0 \endcode The number of frames must not be negative.
 * @pre \code first_frame_time >= 0 \endcode The time stamp of the first frame
 * must not be negative.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_with_time_new
                (const tng_trajectory_t tng_data,
                 const int64_t first_frame,
                 const int64_t n_frames,
                 const double first_frame_time);

/**
 * @brief Set the time stamp of the first frame of the current frame set.
 * @param tng_data is the trajectory containing the frame set.
 * @param first_frame_time is the time stamp of the first frame in the
 * frame set.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code first_frame_time >= 0 \endcode The time stamp of the first frame
 * must not be negative.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_set_first_frame_time_set
                (const tng_trajectory_t tng_data,
                 const double first_frame_time);

/**
 * @brief Read the number of the first frame of the next frame set.
 * @param tng_data is the trajectory containing the frame set.
 * @param frame is set to the frame number of the first frame in the
 * next frame set.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code tng_data->input_file != 0 \endcode An input file must be open
 * to find the next frame set.
 * @pre \code frame != 0 \endcode The pointer to the frame must not be a NULL
 * pointer.
 * @return TNG_SUCCESS(0) if successful, TNG_FAILURE(1) if there is no next
 * frame set or TNG_CRITICAL(2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_first_frame_nr_of_next_frame_set_get
                (const tng_trajectory_t tng_data,
                 int64_t *frame);

/**
 * @brief Add a non-particle dependent data block.
 * @param tng_data is the trajectory data container in which to add the data
 * block
 * @param id is the block ID of the block to add.
 * @param block_name is a descriptive name of the block to add
 * @param datatype is the datatype of the data in the block (e.g. int/float)
 * @param block_type_flag indicates if this is a non-trajectory block (added
 * directly to tng_data) or if it is a trajectory block (added to the
 * frame set)
 * @param n_frames is the number of frames of the data block (automatically
 * set to 1 if adding a non-trajectory data block)
 * @param n_values_per_frame is how many values a stored each frame (e.g. 9
 * for a box shape block)
 * @param stride_length is how many frames are between each entry in the
 * data block
 * @param codec_id is the ID of the codec to compress the data.
 * @param new_data is an array of data values to add.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code block_name != 0 \endcode The pointer to the block name must
 * not be a NULL pointer.
 * @pre \code n_values_per_frame > 0 \endcode n_values_per_frame must be
 * a positive integer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_data_block_add
                (const tng_trajectory_t tng_data,
                 const int64_t id,
                 const char *block_name,
                 const char datatype,
                 const char block_type_flag,
                 int64_t n_frames,
                 const int64_t n_values_per_frame,
                 int64_t stride_length,
                 const int64_t codec_id,
                 void *new_data);

/**
 * @brief Add a particle dependent data block.
 * @param tng_data is the trajectory data container in which to add the data
 * block
 * @param id is the block ID of the block to add.
 * @param block_name is a descriptive name of the block to add
 * @param datatype is the datatype of the data in the block (e.g. int/float)
 * @param block_type_flag indicates if this is a non-trajectory block (added
 * directly to tng_data) or if it is a trajectory block (added to the
 * frame set)
 * @param n_frames is the number of frames of the data block (automatically
 * set to 1 if adding a non-trajectory data block)
 * @param n_values_per_frame is how many values a stored each frame (e.g. 9
 * for a box shape block)
 * @param stride_length is how many frames are between each entry in the
 * data block
 * @param num_first_particle is the number of the first particle stored
 * in this data block
 * @param n_particles is the number of particles stored in this data block
 * @param codec_id is the ID of the codec to compress the data.
 * @param new_data is an array of data values to add.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code block_name != 0 \endcode The pointer to the block name must
 * not be a NULL pointer.
 * @pre \code n_values_per_frame > 0 \endcode n_values_per_frame must be
 * a positive integer.
 * @pre \code num_first_particle >= 0 \endcode The number of the
 * first particle must be >= 0.
 * @pre \code n_particles >= 0 \endcode n_particles must be >= 0.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_particle_data_block_add
                (const tng_trajectory_t tng_data,
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
                 void *new_data);

/** @brief Get the name of a data block of a specific ID.
 * @param tng_data is the trajectory data container.
 * @param block_id is the ID of the data block of which to get the name.
 * @param name is a string, which is set to the name of the data block.
 * Memory must be reserved beforehand.
 * @param max_len is the maximum length of name.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code name != 0 \endcode The pointer to the name string
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if the data block is found, TNG_FAILURE (1)
 * if a minor error has occured or the data block is not found or
 * TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_data_block_name_get
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 char *name,
                 const int max_len);

/** @brief Get the dependency of a data block of a specific ID.
 * @param tng_data is the trajectory data container.
 * @param block_id is the ID of the data block of which to get the name.
 * @param block_dependency is a pointer to the dependency of the data block.
 * If the block is frame dependent it will be set to TNG_FRAME_DEPENDENT,
 * if it is particle dependent it will be set to TNG_PARTICLE_DEPENDENT and
 * if it is both it will be set to TNG_FRAME_DEPENDENT & TNG_PARTICLE_DEPENDENT.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code block_dependency != 0 \endcode The pointer to the block dependency
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if the data block is found, TNG_FAILURE (1)
 * if a minor error has occured or the data block is not found or
 * TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_data_block_dependency_get
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 int *block_dependency);

/** @brief Get the number of values per frame of a data block of a specific ID.
 * @param tng_data is the trajectory data container.
 * @param block_id is the ID of the data block of which to get the name.
 * @param n_values_per_frame is a pointer set to the number of values per frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n_values_per_frame != 0 \endcode The pointer to the number of values
 * per frame must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if the data block is found, TNG_FAILURE (1)
 * if a minor error has occured or the data block is not found or
 * TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_data_block_num_values_per_frame_get
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 int64_t *n_values_per_frame);

/**
 * @brief Write data of one trajectory frame to the output_file of tng_data.
 * @param tng_data is a trajectory data container. tng_data->output_file_path
 * specifies which file to write to. If the file (output_file) is not open it
 * will be opened.
 * @param frame_nr is the index number of the frame to write.
 * @param block_id is the ID of the data block to write the data to.
 * @param values is an array of data to write. The length of the array should
 * equal n_values_per_frame.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code frame_nr >= 0 \endcode The frame number to write must be >= 0.
 * @pre \code values != 0 \endcode The pointer to the values must not be a NULL
 * pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_data_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const int64_t block_id,
                 const void *values,
                 const char hash_mode);

/**
 * @brief Write particle data of one trajectory frame to the output_file of
 * tng_data.
 * @param tng_data is a trajectory data container. tng_data->output_file_path
 * specifies which file to write to. If the file (output_file) is not open it
 * will be opened.
 * @param frame_nr is the index number of the frame to write.
 * @param block_id is the ID of the data block to write the data to.
 * @param val_first_particle is the number of the first particle in the data
 * array.
 * @param val_n_particles is the number of particles in the data array.
 * @param values is a 1D-array of data to write. The length of the array should
 * equal n_particles * n_values_per_frame.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code frame_nr >= 0 \endcode The frame number to write must be >= 0.
 * @pre \code val_first_particle >= 0 \endcode The number of the
 * first particle must be >= 0.
 * @pre \code val_n_particles >= 0 \endcode The number of particles must be >= 0.
 * @pre \code values != 0 \endcode The pointer to the values must not be a NULL
 * pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_frame_particle_data_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const int64_t block_id,
                 const int64_t val_first_particle,
                 const int64_t val_n_particles,
                 const void *values,
                 const char hash_mode);

/**
 * @brief Free data of an array of values (2D).
 * @param tng_data is a trajectory data container.
 * @param values is the 2D array to free and will be set to 0 afterwards.
 * @param n_frames is the number of frames in the data array.
 * @param n_values_per_frame is the number of values per frame in the data array.
 * @param type is the data type of the data in the array (e.g. int/float/char).
 * @details This function should not be used. The data_values union is obsolete.
 * This function also causes memory leaks, but its signature cannot be changed
 * without disturbing the API.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_data_values_free
                (const tng_trajectory_t tng_data,
                 union data_values **values,
                 const int64_t n_frames,
                 const int64_t n_values_per_frame,
                 const char type);

/**
 * @brief Free data of an array of values (3D).
 * @param tng_data is a trajectory data container.
 * @param values is the array to free and will be set to 0 afterwards.
 * @param n_frames is the number of frames in the data array.
 * @param n_particles is the number of particles in the data array.
 * @param n_values_per_frame is the number of values per frame in the data array.
 * @param type is the data type of the data in the array (e.g. int/float/char).
 * @details This function should not be used. The data_values union is obsolete.
 * This function also causes memory leaks, but its signature cannot be changed
 * without disturbing the API.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_particle_data_values_free
                (const tng_trajectory_t tng_data,
                 union data_values ***values,
                 const int64_t n_frames,
                 const int64_t n_particles,
                 const int64_t n_values_per_frame,
                 const char type);

/**
 * @brief Retrieve non-particle data, from the last read frame set. Obsolete!
 * @param tng_data is a trajectory data container. tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param block_id is the id number of the particle data block to read.
 * @param values is a pointer to a 2-dimensional array (memory unallocated), which
 * will be filled with data. The array will be sized
 * (n_frames * n_values_per_frame).
 * Since ***values is allocated in this function it is the callers
 * responsibility to free the memory.
 * @param n_frames is set to the number of frames in the returned data. This is
 * needed to properly reach and/or free the data afterwards.
 * @param n_values_per_frame is set to the number of values per frame in the data.
 * This is needed to properly reach and/or free the data afterwards.
 * @param type is set to the data type of the data in the array.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n_frames != 0 \endcode The pointer to the number of frames
 * must not be a NULL pointer.
 * @pre \code n_values_per_frame != 0 \endcode The pointer to the number of
 * values per frame must not be a NULL pointer.
 * @pre \code type != 0 \endcode The pointer to the data type must not
 * be a NULL pointer.
 * @details This function is obsolete and only retained for compatibility. Use
 * tng_data_vector_get() instead.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_data_get(const tng_trajectory_t tng_data,
                                                   const int64_t block_id,
                                                   union data_values ***values,
                                                   int64_t *n_frames,
                                                   int64_t *n_values_per_frame,
                                                   char *type);

/**
 * @brief Retrieve a vector (1D array) of non-particle data, from the last read frame set.
 * @param tng_data is a trajectory data container. tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param block_id is the id number of the particle data block to read.
 * @param values is a pointer to a 1-dimensional array (memory unallocated), which
 * will be filled with data. The length of the array will be (n_frames * n_values_per_frame).
 * Since **values is allocated in this function it is the callers
 * responsibility to free the memory.
 * @param n_frames is set to the number of particles in the returned data. This is
 * needed to properly reach and/or free the data afterwards.
 * @param stride_length is set to the stride length of the returned data.
 * @param n_values_per_frame is set to the number of values per frame in the data.
 * This is needed to properly reach and/or free the data afterwards.
 * @param type is set to the data type of the data in the array.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n_frames != 0 \endcode The pointer to the number of frames
 * must not be a NULL pointer.
 * @pre \code stride_length != 0 \endcode The pointer to the stride length
 * must not be a NULL pointer.
 * @pre \code n_values_per_frame != 0 \endcode The pointer to the number of
 * values per frame must not be a NULL pointer.
 * @pre \code type != 0 \endcode The pointer to the data type must not
 * be a NULL pointer.
 * @details This does only work for numerical (int, float, double) data.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_data_vector_get
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 void **values,
                 int64_t *n_frames,
                 int64_t *stride_length,
                 int64_t *n_values_per_frame,
                 char *type);

/**
 * @brief Read and retrieve non-particle data, in a specific interval. Obsolete!
 * @param tng_data is a trajectory data container. tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param block_id is the id number of the particle data block to read.
 * @param start_frame_nr is the index number of the first frame to read.
 * @param end_frame_nr is the index number of the last frame to read.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @param values is a pointer to a 2-dimensional array (memory unallocated), which
 * will be filled with data. The array will be sized
 * (n_frames * n_values_per_frame).
 * Since ***values is allocated in this function it is the callers
 * responsibility to free the memory.
 * @param n_values_per_frame is set to the number of values per frame in the data.
 * This is needed to properly reach and/or free the data afterwards.
 * @param type is set to the data type of the data in the array.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code start_frame_nr <= end_frame_nr \endcode The first frame must be before
 * the last frame.
 * @pre \code n_values_per_frame != 0 \endcode The pointer to the number of
 * values per frame must not be a NULL pointer.
 * @pre \code type != 0 \endcode The pointer to the data type must not
 * be a NULL pointer.
 * @details This function is obsolete and only retained for compatibility. Use
 * tng_data_vector_interval_get() instead.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_data_interval_get
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 const int64_t start_frame_nr,
                 const int64_t end_frame_nr,
                 const char hash_mode,
                 union data_values ***values,
                 int64_t *n_values_per_frame,
                 char *type);

/**
 * @brief Read and retrieve a vector (1D array) of non-particle data,
 * in a specific interval.
 * @param tng_data is a trajectory data container. tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param block_id is the id number of the particle data block to read.
 * @param start_frame_nr is the index number of the first frame to read.
 * @param end_frame_nr is the index number of the last frame to read.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @param values is a pointer to a 1-dimensional array (memory unallocated), which
 * will be filled with data. The length of the array will be (n_frames * n_values_per_frame).
 * Since **values is allocated in this function it is the callers
 * responsibility to free the memory.
 * @param stride_length is set to the stride length (writing interval) of
 * the data.
 * @param n_values_per_frame is set to the number of values per frame in the data.
 * This is needed to properly reach and/or free the data afterwards.
 * @param type is set to the data type of the data in the array.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code start_frame_nr <= end_frame_nr \endcode The first frame must be before
 * the last frame.
 * @pre \code stride_length != 0 \endcode The pointer to the stride length
 * must not be a NULL pointer.
 * @pre \code n_values_per_frame != 0 \endcode The pointer to the number of
 * values per frame must not be a NULL pointer.
 * @pre \code type != 0 \endcode The pointer to the data type must not
 * be a NULL pointer.
 * @details This does only work for numerical (int, float, double) data.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_data_vector_interval_get
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 const int64_t start_frame_nr,
                 const int64_t end_frame_nr,
                 const char hash_mode,
                 void **values,
                 int64_t *stride_length,
                 int64_t *n_values_per_frame,
                 char *type);

/**
 * @brief Retrieve particle data, from the last read frame set. Obsolete!
 * @details The particle dimension of the returned values array is translated
 * to real particle numbering, i.e. the numbering of the actual molecular
 * system.
 * @param tng_data is a trajectory data container. tng_data->input_file_path
 * specifies which file to read from. If the file (input_file) is not open it
 * will be opened.
 * @param block_id is the id number of the particle data block to read.
 * @param values is a pointer to a 3-dimensional array (memory unallocated), which
 * will be filled with data. The array will be sized
 * (n_frames * n_particles * n_values_per_frame).
 * Since ****values is allocated in this function it is the callers
 * responsibility to free the memory.
 * @param n_frames is set to the number of frames in the returned data. This is
 * needed to properly reach and/or free the data afterwards.
 * @param n_particles is set to the number of particles in the returned data. This is
 * needed to properly reach and/or free the data afterwards.
 * @param n_values_per_frame is set to the number of values per frame in the data.
 * This is needed to properly reach and/or free the data afterwards.
 * @param type is set to the data type of the data in the array.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n_frames != 0 \endcode The pointer to the number of frames
 * must not be a NULL pointer.
 * @pre \code n_particles != 0 \endcode The pointer to the number of particles must
 * not be a NULL pointer.
 * @pre \code n_values_per_frame != 0 \endcode The pointer to the number of
 * values per frame must not be a NULL pointer.
 * @pre \code type != 0 \endcode The pointer to the data type must not
 * be a NULL pointer.
 * @details This function is obsolete and only retained for compatibility. Use
 * tng_particle_data_vector_get() instead.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_particle_data_get
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 union data_values ****values,
                 int64_t *n_frames,
                 int64_t *n_particles,
                 int64_t *n_values_per_frame,
                 char *type);

/**
 * @brief Retrieve a vector (1D array) of particle data, from the last read frame set.
 * @details The particle dimension of the returned values array is translated
 * to real particle numbering, i.e. the numbering of the actual molecular
 * system.
 * @param tng_data is a trajectory data container. tng_data->input_file_path
 * specifies which file to read from. If the file (input_file) is not open it
 * will be opened.
 * @param block_id is the id number of the particle data block to read.
 * @param values is a pointer to a 1-dimensional array (memory unallocated), which
 * will be filled with data. The length of the array will be
 * (n_frames * n_particles * n_values_per_frame).
 * Since **values is allocated in this function it is the callers
 * responsibility to free the memory.
 * @param n_frames is set to the number of frames in the returned data. This is
 * needed to properly reach and/or free the data afterwards.
 * @param stride_length is set to the stride length of the returned data.
 * @param n_particles is set to the number of particles in the returned data. This is
 * needed to properly reach and/or free the data afterwards.
 * @param n_values_per_frame is set to the number of values per frame in the data.
 * This is needed to properly reach and/or free the data afterwards.
 * @param type is set to the data type of the data in the array.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n_particles != 0 \endcode The pointer to the number of particles must
 * not be a NULL pointer.
 * @pre \code stride_length != 0 \endcode The pointer to the stride length
 * must not be a NULL pointer.
 * @pre \code n_values_per_frame != 0 \endcode The pointer to the number of
 * values per frame must not be a NULL pointer.
 * @pre \code type != 0 \endcode The pointer to the data type must not
 * be a NULL pointer.
 * @details This does only work for numerical (int, float, double) data.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_particle_data_vector_get
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 void **values,
                 int64_t *n_frames,
                 int64_t *stride_length,
                 int64_t *n_particles,
                 int64_t *n_values_per_frame,
                 char *type);

/**
 * @brief Read and retrieve particle data, in a specific interval. Obsolete!
 * @details The particle dimension of the returned values array is translated
 * to real particle numbering, i.e. the numbering of the actual molecular
 * system.
 * @param tng_data is a trajectory data container. tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param block_id is the id number of the particle data block to read.
 * @param start_frame_nr is the index number of the first frame to read.
 * @param end_frame_nr is the index number of the last frame to read.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @param values is a pointer to a 3-dimensional array (memory unallocated), which
 * will be filled with data. The array will be sized
 * (n_frames * n_particles * n_values_per_frame).
 * Since ****values is allocated in this function it is the callers
 * responsibility to free the memory.
 * @param n_particles is set to the number of particles in the returned data. This is
 * needed to properly reach and/or free the data afterwards.
 * @param n_values_per_frame is set to the number of values per frame in the data.
 * This is needed to properly reach and/or free the data afterwards.
 * @param type is set to the data type of the data in the array.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n_frames != 0 \endcode The pointer to the number of frames
 * must not be a NULL pointer.
 * @pre \code start_frame_nr <= end_frame_nr \endcode The first frame must be before
 * the last frame.
 * @pre \code n_particles != 0 \endcode The pointer to the number of particles must
 * not be a NULL pointer.
 * @pre \code n_values_per_frame != 0 \endcode The pointer to the number of
 * values per frame must not be a NULL pointer.
 * @pre \code type != 0 \endcode The pointer to the data type must not
 * be a NULL pointer.
 * @details This function is obsolete and only retained for compatibility. Use
 * tng_particle_data_vector_interval_get() instead.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_particle_data_interval_get
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 const int64_t start_frame_nr,
                 const int64_t end_frame_nr,
                 const char hash_mode,
                 union data_values ****values,
                 int64_t *n_particles,
                 int64_t *n_values_per_frame,
                 char *type);

/**
 * @brief Read and retrieve a vector (1D array) particle data, in a
 * specific interval.
 * @details The particle dimension of the returned values array is translated
 * to real particle numbering, i.e. the numbering of the actual molecular
 * system.
 * @param tng_data is a trajectory data container. tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param block_id is the id number of the particle data block to read.
 * @param start_frame_nr is the index number of the first frame to read.
 * @param end_frame_nr is the index number of the last frame to read.
 * @param hash_mode is an option to decide whether to use the md5 hash or not.
 * If hash_mode == TNG_USE_HASH the md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @param values is a pointer to a 1-dimensional array (memory unallocated), which
 * will be filled with data. The length of the array will be
 * (n_frames * n_particles * n_values_per_frame).
 * Since **values is allocated in this function it is the callers
 * responsibility to free the memory.
 * @param stride_length is set to the stride length (writing interval) of
 * the data.
 * @param n_particles is set to the number of particles in the returned data. This is
 * needed to properly reach and/or free the data afterwards.
 * @param n_values_per_frame is set to the number of values per frame in the data.
 * This is needed to properly reach and/or free the data afterwards.
 * @param type is set to the data type of the data in the array.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code start_frame_nr <= end_frame_nr \endcode The first frame must be before
 * the last frame.
 * @pre \code n_particles != 0 \endcode The pointer to the number of particles must
 * not be a NULL pointer.
 * @pre \code stride_length != 0 \endcode The pointer to the stride length
 * must not be a NULL pointer.
 * @pre \code n_values_per_frame != 0 \endcode The pointer to the number of
 * values per frame must not be a NULL pointer.
 * @pre \code type != 0 \endcode The pointer to the data type must not
 * be a NULL pointer.
 * @details This does only work for numerical (int, float, double) data.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_particle_data_vector_interval_get
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 const int64_t start_frame_nr,
                 const int64_t end_frame_nr,
                 const char hash_mode,
                 void **values,
                 int64_t *n_particles,
                 int64_t *stride_length,
                 int64_t *n_values_per_frame,
                 char *type);

/**
 * @brief Get the stride length of a specific data (particle dependency does not matter)
 * block, either in the current frame set or of a specific frame.
 * @param tng_data is the trajectory data container.
 * @param block_id is the block ID of the data block, of which to retrieve the
 * stride length of the data.
 * @param frame is the frame from which to get the stride length. If frame is set to -1
 * no specific frame will be used, but instead the first frame, starting from the last read
 * frame set, containing the data block will be used.
 * @param stride_length is set to the value of the stride length of the data block.
 * @return  TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_data_get_stride_length
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 int64_t frame,
                 int64_t *stride_length);

/**
 * @brief Get the date and time of initial file creation in ISO format (string).
 * @param tng_data is a trajectory data container.
 * @param time is a pointer to the string in which the date will be stored. Memory
 * must be reserved beforehand.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code time != 0 \endcode The pointer to the time must not be a NULL
 * pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_time_get_str
                (const tng_trajectory_t tng_data,
                 char *time);
/** @} */ /* end of group1 */

/** @defgroup group2 High-level API
 *  These functions make it easier to access and output TNG data. They
 *  are recommended unless there is a special reason to use the more
 *  detailed functions available in the low-level API.
 *  @{
 */

/**
 * @brief High-level function for opening and initializing a TNG trajectory.
 * @param filename is a string containing the name of the trajectory to open.
 * @param mode specifies the file mode of the trajectory. Can be set to 'r',
 * 'w' or 'a' for reading, writing or appending respectively.
 * @param tng_data_p is a pointer to the opened trajectory. This will be
 * allocated by the TNG library. The trajectory must be
 * closed by the user, whereby memory is freed.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code filename != 0 \endcode The pointer to the filename must not be a
 * NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_open
                (const char *filename,
                 const char mode,
                 tng_trajectory_t *tng_data_p);

/**
 * @brief High-level function for closing a TNG trajectory.
 * @param tng_data_p is a pointer to the trajectory to close. The memory
 * will be freed after finalising the writing.
 * @return TNG_SUCCESS (0) if successful.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_close
                (tng_trajectory_t *tng_data_p);

/**
 * @brief High-level function for getting the time (in seconds) of a frame.
 * @param tng_data is the trajectory containing the frame.
 * @param frame_nr is the frame number of which to get the time.
 * @param time is set to the time (in seconds) of the specified frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code time != 0 \endcode The pointer to the time must not be a
 * NULL pointer.
 * @return TNG_SUCCESS (0) if successful or TNG_FAILURE (1) if a
 * minor error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_time_of_frame_get
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 double *time);

/*
 * @brief High-level function for getting the molecules in the mol system.
 * @param tng_data is the trajectory containing the mol system.
 * @param n_mols is set to the number of molecules in the system.
 * @param molecule_cnt_list will be pointing to the list of counts of each molecule
 * in the mol system.
 * @param mols pointing to the list of molecules in the mol system.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n_mols != 0 \endcode The pointer to the number of molecules must
 * not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful.
 */
/*tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_molecules_get
                (const tng_trajectory_t tng_data,
                 int64_t *n_mols,
                 int64_t **molecule_cnt_list,
                 tng_molecule_t *mols);
*/
/*
 * @brief High-level function for adding a molecule to the mol system.
 * @param tng_data is the trajectory containing the mol system.
 * @param name is the name of the molecule to add.
 * @param cnt is the count of the molecule.
 * @param mol is set to point to the newly created molecule.
 * @pre \code name != 0 \endcode The pointer to the name must not be a
 * NULL pointer.
 * @pre \code cnt >= 0 \endcode The requested count must be >= 0.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured or TNG_CRITICAL (2) if a major error has occured.
 */
/*tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_molecule_add
                (const tng_trajectory_t tng_data,
                 const char *name,
                 const int64_t cnt,
                 tng_molecule_t *mol);
*/
/*
// tng_function_status DECLSPECDLLEXPORT tng_util_molecule_particles_get
//                 (const tng_trajectory_t tng_data,
//                  const tng_molecule_t mol,
//                  int64_t *n_particles,
//                  char ***names,
//                  char ***types,
//                  char ***res_names,
//                  int64_t **res_ids,
//                  char ***chain_names,
//                  int64_t **chain_ids);
//
// tng_function_status DECLSPECDLLEXPORT tng_util_molecule_particles_set
//                 (const tng_trajectory_t tng_data,
//                  tng_molecule_t mol,
//                  const int64_t n_particles,
//                  const char **names,
//                  const char **types,
//                  const char **res_names,
//                  const int64_t *res_ids,
//                  const char **chain_names,
//                  const int64_t *chain_ids);
*/
/**
 * @brief High-level function for reading the positions of all particles
 * from all frames.
 * @param tng_data is the trajectory to read from.
 * @param positions will be set to point at a 1-dimensional array of floats,
 * which will contain the positions. The data is stored sequentially in order
 * of frames. For each frame the positions (x, y and z coordinates) are stored.
 * The variable may point at already allocated memory or be a NULL pointer.
 * The memory must be freed afterwards.
 * @param stride_length will be set to the writing interval of the stored data.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code positions != 0 \endcode The pointer to the positions array
 * must not be a NULL pointer.
 * @pre \code stride_length != 0 \endcode The pointer to the stride length
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_pos_read
                (const tng_trajectory_t tng_data,
                 float **positions,
                 int64_t *stride_length);

/**
 * @brief High-level function for reading the velocities of all particles
 * from all frames.
 * @param tng_data is the trajectory to read from.
 * @param velocities will be set to point at a 1-dimensional array of floats,
 * which will contain the velocities. The data is stored sequentially in order
 * of frames. For each frame the velocities (in x, y and z) are stored. The
 * variable may point at already allocated memory or be a NULL pointer.
 * The memory must be freed afterwards.
 * @param stride_length will be set to the writing interval of the stored data.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code velocities != 0 \endcode The pointer to the velocities array
 * must not be a NULL pointer.
 * @pre \code stride_length != 0 \endcode The pointer to the stride length
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_vel_read
                (const tng_trajectory_t tng_data,
                 float **velocities,
                 int64_t *stride_length);

/**
 * @brief High-level function for reading the forces of all particles
 * from all frames.
 * @param tng_data is the trajectory to read from.
 * @param forces will be set to point at a 1-dimensional array of floats,
 * which will contain the forces. The data is stored sequentially in order
 * of frames. For each frame the forces (in x, y and z) are stored. The
 * variable may point at already allocated memory or be a NULL pointer.
 * The memory must be freed afterwards.
 * @param stride_length will be set to the writing interval of the stored data.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code forces != 0 \endcode The pointer to the forces array
 * must not be a NULL pointer.
 * @pre \code stride_length != 0 \endcode The pointer to the stride length
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_force_read
                (const tng_trajectory_t tng_data,
                 float **forces,
                 int64_t *stride_length);

/**
 * @brief High-level function for reading the box shape from all frames.
 * @param tng_data is the trajectory to read from.
 * @param box_shape will be set to point at a 1-dimensional array of floats,
 * which will contain the box shape. The data is stored sequentially in order
 * of frames. The variable may point at already allocated memory or be a NULL pointer.
 * If the box shape is not modified during the trajectory, but as general data,
 * that will be returned instead.
 * @param stride_length will be set to the writing interval of the stored data.
 * @details This function should only be used if number of values used to specify
 * the box shape is known (by default TNG uses 9 values) since it does not
 * return the number of values in the array. It is recommended to use
 * tng_data_vector_interval_get() instead.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code box_shape != 0 \endcode The pointer to the box_shape array
 * must not be a NULL pointer.
 * @pre \code stride_length != 0 \endcode The pointer to the stride length
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_read
                (const tng_trajectory_t tng_data,
                 float **box_shape,
                 int64_t *stride_length);

/**
 * @brief High-level function for reading the next frame of particle-dependent
 * data of a specific type.
 * @param tng_data is the trajectory to read from.
 * @param block_id is the ID number of the block containing the data of interest.
 * @param values will be set to point at a 1-dimensional array containing the
 * requested data. The variable may point at already allocated memory or be a
 * NULL pointer. The memory must be freed afterwards.
 * @param data_type will be pointing to a character indicating the size of the
 * data of the returned values, e.g. TNG_INT_DATA, TNG_FLOAT_DATA or TNG_DOUBLE_DATA.
 * @param retrieved_frame_number will be pointing at the frame number of the
 * returned frame.
 * @param retrieved_time will be pointing at the time stamp of the returned
 * frame.
 * @details If no frame has been read before the first frame of the trajectory
 * is read.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code values != 0 \endcode The pointer to the values array
 * must not be a NULL pointer.
 * @pre \code data_type != 0 \endcode The pointer to the data type of the
 * returned data must not be a NULL pointer.
 * @pre \code retrieved_frame_number != 0 \endcode The pointer to the frame
 * number of the returned data must not be a NULL pointer.
 * @pre \code retrieved_time != 0 \endcode The pointer to the time of the
 * returned data must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_particle_data_next_frame_read
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 void **values,
                 char *data_type,
                 int64_t *retrieved_frame_number,
                 double *retrieved_time);

/**
 * @brief High-level function for reading the next frame of non-particle-dependent
 * data of a specific type.
 * @param tng_data is the trajectory to read from.
 * @param block_id is the ID number of the block containing the data of interest.
 * @param values will be set to point at a 1-dimensional array containing the
 * requested data. The variable may point at already allocated memory or be a
 * NULL pointer. The memory must be freed afterwards.
 * @param data_type will be pointing to a character indicating the size of the
 * data of the returned values, e.g. TNG_INT_DATA, TNG_FLOAT_DATA or TNG_DOUBLE_DATA.
 * @param retrieved_frame_number will be pointing at the frame number of the
 * returned frame.
 * @param retrieved_time will be pointing at the time stamp of the returned
 * frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code values != 0 \endcode The pointer to the values array
 * must not be a NULL pointer.
 * @pre \code data_type != 0 \endcode The pointer to the data type of the
 * returned data must not be a NULL pointer.
 * @pre \code retrieved_frame_number != 0 \endcode The pointer to the frame
 * number of the returned data must not be a NULL pointer.
 * @pre \code retrieved_time != 0 \endcode The pointer to the time of the
 * returned data must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_non_particle_data_next_frame_read
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 void **values,
                 char *data_type,
                 int64_t *retrieved_frame_number,
                 double *retrieved_time);

/**
 * @brief High-level function for reading the positions of all particles
 * from a specific range of frames.
 * @param tng_data is the trajectory to read from.
 * @param first_frame is the first frame to return position data from.
 * @param last_frame is the last frame to return position data from.
 * @param positions will be set to point at a 1-dimensional array of floats,
 * which will contain the positions. The data is stored sequentially in order
 * of frames. For each frame the positions (x, y and z coordinates) are stored.
 * The variable may point at already allocated memory or be a NULL pointer.
 * The memory must be freed afterwards.
 * @param stride_length will be set to the writing interval of the stored data.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code start_frame_nr <= end_frame_nr \endcode The first frame must be before
 * the last frame.
 * @pre \code positions != 0 \endcode The pointer to the positions array
 * must not be a NULL pointer.
 * @pre \code stride_length != 0 \endcode The pointer to the stride length
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_pos_read_range
                (const tng_trajectory_t tng_data,
                 const int64_t first_frame,
                 const int64_t last_frame,
                 float **positions,
                 int64_t *stride_length);

/**
 * @brief High-level function for reading the velocities of all particles
 * from a specific range of frames.
 * @param tng_data is the trajectory to read from.
 * @param first_frame is the first frame to return position data from.
 * @param last_frame is the last frame to return position data from.
 * @param velocities will be set to point at a 1-dimensional array of floats,
 * which will contain the velocities. The data is stored sequentially in order
 * of frames. For each frame the velocities (in x, y and z) are stored. The
 * variable may point at already allocated memory or be a NULL pointer.
 * The memory must be freed afterwards.
 * @param stride_length will be set to the writing interval of the stored data.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code start_frame_nr <= end_frame_nr \endcode The first frame must be before
 * the last frame.
 * @pre \code velocities != 0 \endcode The pointer to the velocities array
 * must not be a NULL pointer.
 * @pre \code stride_length != 0 \endcode The pointer to the stride length
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_vel_read_range
                (const tng_trajectory_t tng_data,
                 const int64_t first_frame,
                 const int64_t last_frame,
                 float **velocities,
                 int64_t *stride_length);

/**
 * @brief High-level function for reading the forces of all particles
 * from a specific range of frames.
 * @param tng_data is the trajectory to read from.
 * @param first_frame is the first frame to return position data from.
 * @param last_frame is the last frame to return position data from.
 * @param forces will be set to point at a 1-dimensional array of floats,
 * which will contain the forces. The data is stored sequentially in order
 * of frames. For each frame the forces (in x, y and z) are stored. The
 * variable may point at already allocated memory or be a NULL pointer.
 * The memory must be freed afterwards.
 * @param stride_length will be set to the writing interval of the stored data.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code start_frame_nr <= end_frame_nr \endcode The first frame must be before
 * the last frame.
 * @pre \code forces != 0 \endcode The pointer to the forces array
 * must not be a NULL pointer.
 * @pre \code stride_length != 0 \endcode The pointer to the stride length
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_force_read_range
                (const tng_trajectory_t tng_data,
                 const int64_t first_frame,
                 const int64_t last_frame,
                 float **forces,
                 int64_t *stride_length);

/**
 * @brief High-level function for reading the box shape
 * from a specific range of frames.
 * @param tng_data is the trajectory to read from.
 * @param first_frame is the first frame to return position data from.
 * @param last_frame is the last frame to return position data from.
 * @param box_shape will be set to point at a 1-dimensional array of floats,
 * which will contain the box shape. The data is stored sequentially in order
 * of frames.
 * If the box shape is not modified during the trajectory, but as general data,
 * that will be returned instead. The
 * variable may point at already allocated memory or be a NULL pointer.
 * The memory must be freed afterwards.
 * @param stride_length will be set to the writing interval of the stored data.
 * @details This function should only be used if number of values used to specify
 * the box shape is known (by default TNG uses 9 values) since it does not
 * return the number of values in the array. It is recommended to use
 * tng_data_vector_interval_get() instead.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code start_frame_nr <= end_frame_nr \endcode The first frame must be before
 * the last frame.
 * @pre \code box_shape != 0 \endcode The pointer to the box_shape array
 * must not be a NULL pointer.
 * @pre \code stride_length != 0 \endcode The pointer to the stride length
 * must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_read_range
                (const tng_trajectory_t tng_data,
                 const int64_t first_frame,
                 const int64_t last_frame,
                 float **box_shape,
                 int64_t *stride_length);

/**
 * @brief High-level function for setting the writing interval of data blocks.
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @param n_values_per_frame is the number of values to store per frame. If the
 * data is particle dependent there will be n_values_per_frame stored per
 * particle each frame.
 * @param block_id is the ID of the block, of which to set the output interval.
 * @param block_name is a string that will be used as name of the block. Only
 * required if the block did not exist, i.e. a new block is created.
 * @param particle_dependency should be TNG_NON_PARTICLE_BLOCK_DATA (0) if the
 * data is not related to specific particles (e.g. box shape) or
 * TNG_PARTICLE_BLOCK_DATA (1) is it is related to specific particles (e.g.
 * positions). Only required if the block did not exist, i.e. a new block is
 * created.
 * @param compression is the compression routine to use when writing the data.
 * Only required if the block did not exist, i.e. a new block is created.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code i >= 0 \endcode The writing interval must be >= 0.
 * @details n_values_per_frame, block_name, particle_dependency and
 * compression are only used if the data block did not exist before calling
 * this function, in which case it is created.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_generic_write_interval_set
                (const tng_trajectory_t tng_data,
                 const int64_t i,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression);

/**
 * @brief High-level function for setting the writing interval of data blocks
 * containing double precision data.
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @param n_values_per_frame is the number of values to store per frame. If the
 * data is particle dependent there will be n_values_per_frame stored per
 * particle each frame.
 * @param block_id is the ID of the block, of which to set the output interval.
 * @param block_name is a string that will be used as name of the block. Only
 * required if the block did not exist, i.e. a new block is created.
 * @param particle_dependency should be TNG_NON_PARTICLE_BLOCK_DATA (0) if the
 * data is not related to specific particles (e.g. box shape) or
 * TNG_PARTICLE_BLOCK_DATA (1) is it is related to specific particles (e.g.
 * positions). Only required if the block did not exist, i.e. a new block is
 * created.
 * @param compression is the compression routine to use when writing the data.
 * Only required if the block did not exist, i.e. a new block is created.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code i >= 0 \endcode The writing interval must be >= 0.
 * @details n_values_per_frame, block_name, particle_dependency and
 * compression are only used if the data block did not exist before calling
 * this function, in which case it is created.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_generic_write_interval_double_set
                (const tng_trajectory_t tng_data,
                 const int64_t i,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression);

/**
 * @brief High-level function for setting the writing interval of data blocks.
 * Obsolete! Use tng_util_generic_write_interval_set()
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @param n_values_per_frame is the number of values to store per frame. If the
 * data is particle dependent there will be n_values_per_frame stored per
 * particle each frame.
 * @param block_id is the ID of the block, of which to set the output interval.
 * @param block_name is a string that will be used as name of the block. Only
 * required if the block did not exist, i.e. a new block is created.
 * @param particle_dependency should be TNG_NON_PARTICLE_BLOCK_DATA (0) if the
 * data is not related to specific particles (e.g. box shape) or
 * TNG_PARTICLE_BLOCK_DATA (1) is it is related to specific particles (e.g.
 * positions). Only required if the block did not exist, i.e. a new block is
 * created.
 * @param compression is the compression routine to use when writing the data.
 * Only required if the block did not exist, i.e. a new block is created.
 * @details n_values_per_frame, block_name, particle_dependency and
 * compression are only used if the data block did not exist before calling
 * this function, in which case it is created.
 * This function is replaced by the more correcly named
 * tng_util_generic_write_interval_set(), but is kept for compatibility.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_generic_write_frequency_set
                (const tng_trajectory_t tng_data,
                 const int64_t i,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression);

/**
 * @brief High-level function for setting the writing interval of position
 * data blocks.
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code i >= 0 \endcode The writing interval must be >= 0.
 * @details This function uses tng_util_generic_write_interval_set() and will
 * create a positions data block if none exists.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_pos_write_interval_set
                (const tng_trajectory_t tng_data,
                 const int64_t i);

/**
 * @brief High-level function for setting the writing interval of position
 * data blocks containing double precision data.
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code i >= 0 \endcode The writing interval must be >= 0.
 * @details This function uses tng_util_generic_write_interval_set() and will
 * create a positions data block if none exists.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_pos_write_interval_double_set
                (const tng_trajectory_t tng_data,
                 const int64_t i);

/**
 * @brief High-level function for setting the writing interval of position
 * data blocks. Obsolete! Use tng_util_pos_write_interval_set()
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code i >= 0 \endcode The writing interval must be >= 0.
 * @details This function uses tng_util_generic_write_interval_set() and will
 * create a positions data block if none exists.
 * This function is replaced by the more correcly named
 * tng_util_pos_write_interval_set(), but is kept for compatibility.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_pos_write_frequency_set
                (const tng_trajectory_t tng_data,
                 const int64_t i);

/**
 * @brief High-level function for setting the writing interval of velocity
 * data blocks.
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code i >= 0 \endcode The writing interval must be >= 0.
 * @details This function uses tng_util_generic_write_interval_set() and will
 * create a velocities data block if none exists.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_vel_write_interval_set
                (const tng_trajectory_t tng_data,
                 const int64_t i);

/**
 * @brief High-level function for setting the writing interval of velocity
 * data blocks containing double precision data.
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code i >= 0 \endcode The writing interval must be >= 0.
 * @details This function uses tng_util_generic_write_interval_set() and will
 * create a velocities data block if none exists.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_vel_write_interval_double_set
                (const tng_trajectory_t tng_data,
                 const int64_t i);

/**
 * @brief High-level function for setting the writing interval of velocity
 * data blocks. Obsolete! Use tng_util_vel_write_interval_set()
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @details This function uses tng_util_generic_write_interval_set() and will
 * create a velocities data block if none exists.
 * This function is replaced by the more correcly named
 * tng_util_vel_write_interval_set(), but is kept for compatibility.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_vel_write_frequency_set
                (const tng_trajectory_t tng_data,
                 const int64_t i);

/**
 * @brief High-level function for setting the writing interval of force
 * data blocks.
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code i >= 0 \endcode The writing interval must be >= 0.
 * @details This function uses tng_util_generic_write_interval_set() and will
 * create a forces data block if none exists.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_force_write_interval_set
                (const tng_trajectory_t tng_data,
                 const int64_t i);

/**
 * @brief High-level function for setting the writing interval of force
 * data blocks containing double precision data.
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code i >= 0 \endcode The writing interval must be >= 0.
 * @details This function uses tng_util_generic_write_interval_set() and will
 * create a forces data block if none exists.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_force_write_interval_double_set
                (const tng_trajectory_t tng_data,
                 const int64_t i);

/**
 * @brief High-level function for setting the writing interval of force
 * data blocks. Obsolete! Use tng_util_force_write_interval_set()
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @details This function uses tng_util_generic_write_interval_set() and will
 * create a forces data block if none exists.
 * This function is replaced by the more correcly named
 * tng_util_force_write_interval_set(), but is kept for compatibility.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_force_write_frequency_set
                (const tng_trajectory_t tng_data,
                 const int64_t i);

/**
 * @brief High-level function for setting the writing interval of box shape
 * data blocks.
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code i >= 0 \endcode The writing interval must be >= 0.
 * @details This function uses tng_util_generic_write_interval_set() and will
 * create a box shape data block if none exists.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_write_interval_set
                (const tng_trajectory_t tng_data,
                 const int64_t i);

/**
 * @brief High-level function for setting the writing interval of box shape
 * data blocks containing double precision data.
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code i >= 0 \endcode The writing interval must be >= 0.
 * @details This function uses tng_util_generic_write_interval_set() and will
 * create a box shape data block if none exists.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_write_interval_double_set
                (const tng_trajectory_t tng_data,
                 const int64_t i);

/**
 * @brief High-level function for setting the writing interval of velocity
 * data blocks. Obsolete! Use tng_util_box_shape_write_interval_set()
 * @param tng_data is the trajectory to use.
 * @param i is the output interval, i.e. i == 10 means data written every 10th
 * frame.
 * @details This function uses tng_util_generic_write_interval_set() and will
 * create a box shape data block if none exists.
 * This function is replaced by the more correcly named
 * tng_util_box_shape_write_interval_set(), but is kept for compatibility.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_write_frequency_set
                (const tng_trajectory_t tng_data,
                 const int64_t i);

/**
 * @brief High-level function for writing data of one frame to a data block.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data. If frame_nr < 0 the
 * data is written as non-trajectory data.
 * @param values is a 1D array of data to add. The array should be of length
 * n_particles * n_values_per_frame if writing particle related data, otherwise
 * it should be n_values_per_frame.
 * @param n_values_per_frame is the number of values to store per frame. If the
 * data is particle dependent there will be n_values_per_frame stored per
 * particle each frame.
 * @param block_id is the ID of the block, of which to set the output interval.
 * @param block_name is a string that will be used as name of the block. Only
 * required if the block did not exist, i.e. a new block is created.
 * @param particle_dependency should be TNG_NON_PARTICLE_BLOCK_DATA (0) if the
 * data is not related to specific particles (e.g. box shape) or
 * TNG_PARTICLE_BLOCK_DATA (1) is it is related to specific particles (e.g.
 * positions). Only required if the block did not exist, i.e. a new block is
 * created.
 * @param compression is the compression routine to use when writing the data.
 * Only required if the block did not exist, i.e. a new block is created.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code values != 0 \endcode The pointer to the values array must not
 * be a NULL pointer.
 * @details n_values_per_frame, block_name, particle_dependency and
 * compression are only used if the data block did not exist before calling
 * this function, in which case it is created.
 * N.b. Data is written a whole block at a time. The data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_generic_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const float *values,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression);

/**
 * @brief High-level function for writing data of one frame to a double precision
 * data block.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data. If frame_nr < 0 the
 * data is written as non-trajectory data.
 * @param values is a 1D array of data to add. The array should be of length
 * n_particles * n_values_per_frame if writing particle related data, otherwise
 * it should be n_values_per_frame.
 * @param n_values_per_frame is the number of values to store per frame. If the
 * data is particle dependent there will be n_values_per_frame stored per
 * particle each frame.
 * @param block_id is the ID of the block, of which to set the output interval.
 * @param block_name is a string that will be used as name of the block. Only
 * required if the block did not exist, i.e. a new block is created.
 * @param particle_dependency should be TNG_NON_PARTICLE_BLOCK_DATA (0) if the
 * data is not related to specific particles (e.g. box shape) or
 * TNG_PARTICLE_BLOCK_DATA (1) is it is related to specific particles (e.g.
 * positions). Only required if the block did not exist, i.e. a new block is
 * created.
 * @param compression is the compression routine to use when writing the data.
 * Only required if the block did not exist, i.e. a new block is created.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code values != 0 \endcode The pointer to the values array must not
 * be a NULL pointer.
 * @details n_values_per_frame, block_name, particle_dependency and
 * compression are only used if the data block did not exist before calling
 * this function, in which case it is created.
 * N.b. Data is written a whole block at a time. The data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_generic_double_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double *values,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression);

/**
 * @brief High-level function for adding data to positions data blocks.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data. If frame_nr < 0 the
 * data is written as non-trajectory data.
 * @param positions is a 1D array of data to add. The array should be of length
 * n_particles * 3.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code positions != 0 \endcode The pointer to the positions array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_write() and will
 * create a positions data block if none exists. Positions are stored as three
 * values per frame and compressed using TNG compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_pos_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const float *positions);

/**
 * @brief High-level function for adding data to positions data blocks at double
 * precision.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data. If frame_nr < 0 the
 * data is written as non-trajectory data.
 * @param positions is a 1D array of data to add. The array should be of length
 * n_particles * 3.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code positions != 0 \endcode The pointer to the positions array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_write() and will
 * create a positions data block if none exists. Positions are stored as three
 * values per frame and compressed using TNG compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_pos_double_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double *positions);

/**
 * @brief High-level function for adding data to velocities data blocks.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data. If frame_nr < 0 the
 * data is written as non-trajectory data.
 * @param velocities is a 1D array of data to add. The array should be of length
 * n_particles * 3.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code velocities != 0 \endcode The pointer to the velocities array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_write() and will
 * create a velocities data block if none exists. Velocities are stored as three
 * values per frame and compressed using TNG compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_vel_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const float *velocities);

/**
 * @brief High-level function for adding data to velocities data blocks at double
 * precision.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data. If frame_nr < 0 the
 * data is written as non-trajectory data.
 * @param velocities is a 1D array of data to add. The array should be of length
 * n_particles * 3.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code velocities != 0 \endcode The pointer to the velocities array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_write() and will
 * create a velocities data block if none exists. Velocities are stored as three
 * values per frame and compressed using TNG compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_vel_double_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double *velocities);

/**
 * @brief High-level function for adding data to forces data blocks.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data. If frame_nr < 0 the
 * data is written as non-trajectory data.
 * @param forces is a 1D array of data to add. The array should be of length
 * n_particles * 3.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code forces != 0 \endcode The pointer to the forces array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_write() and will
 * create a forces data block if none exists. Forces are stored as three
 * values per frame and compressed using gzip compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_force_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const float *forces);

/**
 * @brief High-level function for adding data to forces data blocks at double
 * precision.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data. If frame_nr < 0 the
 * data is written as non-trajectory data.
 * @param forces is a 1D array of data to add. The array should be of length
 * n_particles * 3.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code forces != 0 \endcode The pointer to the forces array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_write() and will
 * create a forces data block if none exists. Forces are stored as three
 * values per frame and compressed using gzip compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_force_double_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double *forces);

/**
 * @brief High-level function for adding data to box shape data blocks.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data. If frame_nr < 0 the
 * data is written as non-trajectory data.
 * @param box_shape is a 1D array of data to add. The array should be of length 9.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code box_shape != 0 \endcode The pointer to the box_shape array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_write() and will
 * create a box shape data block if none exists. Box shapes are stored as 9
 * values per frame and compressed using TNG compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const float *box_shape);

/**
 * @brief High-level function for adding data to box shape data blocks at double
 * precision.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data. If frame_nr < 0 the
 * data is written as non-trajectory data.
 * @param box_shape is a 1D array of data to add. The array should be of length 9.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code box_shape != 0 \endcode The pointer to the box_shape array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_write() and will
 * create a box shape data block if none exists. Box shapes are stored as 9
 * values per frame and compressed using TNG compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_double_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double *box_shape);

/**
 * @brief High-level function for writing data of one frame to a data block.
 * If the frame is at the beginning of a frame set the time stamp of the frame
 * set is set.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data.
 * @param time is the time stamp of the frame (in seconds).
 * @param values is a 1D array of data to add. The array should be of length
 * n_particles * n_values_per_frame if writing particle related data, otherwise
 * it should be n_values_per_frame.
 * @param n_values_per_frame is the number of values to store per frame. If the
 * data is particle dependent there will be n_values_per_frame stored per
 * particle each frame.
 * @param block_id is the ID of the block, of which to set the output interval.
 * @param block_name is a string that will be used as name of the block. Only
 * required if the block did not exist, i.e. a new block is created.
 * @param particle_dependency should be TNG_NON_PARTICLE_BLOCK_DATA (0) if the
 * data is not related to specific particles (e.g. box shape) or
 * TNG_PARTICLE_BLOCK_DATA (1) is it is related to specific particles (e.g.
 * positions). Only required if the block did not exist, i.e. a new block is
 * created.
 * @param compression is the compression routine to use when writing the data.
 * Only required if the block did not exist, i.e. a new block is created.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code frame_nr >= 0 \endcode The frame number to write must be >= 0.
 * @pre \code time >= 0 \endcode The time stamp must be >= 0.
 * @pre \code values != 0 \endcode The pointer to the values array must not
 * be a NULL pointer.
 * @details n_values_per_frame, block_name, particle_dependency and
 * compression are only used if the data block did not exist before calling
 * this function, in which case it is created.
 * N.b. Data is written a whole block at a time. The data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_generic_with_time_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const float *values,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression);

/**
 * @brief High-level function for writing data of one frame to a double precision
 * data block. If the frame is at the beginning of a frame set the time stamp of
 * the frame set is set.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data.
 * @param time is the time stamp of the frame (in seconds).
 * @param values is a 1D array of data to add. The array should be of length
 * n_particles * n_values_per_frame if writing particle related data, otherwise
 * it should be n_values_per_frame.
 * @param n_values_per_frame is the number of values to store per frame. If the
 * data is particle dependent there will be n_values_per_frame stored per
 * particle each frame.
 * @param block_id is the ID of the block, of which to set the output interval.
 * @param block_name is a string that will be used as name of the block. Only
 * required if the block did not exist, i.e. a new block is created.
 * @param particle_dependency should be TNG_NON_PARTICLE_BLOCK_DATA (0) if the
 * data is not related to specific particles (e.g. box shape) or
 * TNG_PARTICLE_BLOCK_DATA (1) is it is related to specific particles (e.g.
 * positions). Only required if the block did not exist, i.e. a new block is
 * created.
 * @param compression is the compression routine to use when writing the data.
 * Only required if the block did not exist, i.e. a new block is created.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code frame_nr >= 0 \endcode The frame number to write must be >= 0.
 * @pre \code time >= 0 \endcode The time stamp must be >= 0.
 * @pre \code values != 0 \endcode The pointer to the values array must not
 * be a NULL pointer.
 * @details n_values_per_frame, block_name, particle_dependency and
 * compression are only used if the data block did not exist before calling
 * this function, in which case it is created.
 * N.b. Data is written a whole block at a time. The data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_generic_with_time_double_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const double *values,
                 const int64_t n_values_per_frame,
                 const int64_t block_id,
                 const char *block_name,
                 const char particle_dependency,
                 const char compression);

/**
 * @brief High-level function for adding data to positions data blocks. If the
 * frame is at the beginning of a frame set the time stamp of the frame set
 * is set.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data.
 * @param time is the time stamp of the frame (in seconds).
 * @param positions is a 1D array of data to add. The array should be of length
 * n_particles * 3.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code frame_nr >= 0 \endcode The frame number to write must be >= 0.
 * @pre \code time >= 0 \endcode The time stamp must be >= 0.
 * @pre \code positions != 0 \endcode The pointer to the positions array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_with_time_write() and will
 * create a positions data block if none exists. Positions are stored as three
 * values per frame and compressed using TNG compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_pos_with_time_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const float *positions);

/**
 * @brief High-level function for adding data to positions data blocks at double
 * precision. If the frame is at the beginning of a frame set the time stamp of
 * the frame set is set.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data.
 * @param time is the time stamp of the frame (in seconds).
 * @param positions is a 1D array of data to add. The array should be of length
 * n_particles * 3.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code frame_nr >= 0 \endcode The frame number to write must be >= 0.
 * @pre \code time >= 0 \endcode The time stamp must be >= 0.
 * @pre \code positions != 0 \endcode The pointer to the positions array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_with_time_double_write() and will
 * create a positions data block if none exists. Positions are stored as three
 * values per frame and compressed using TNG compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_pos_with_time_double_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const double *positions);

/**
 * @brief High-level function for adding data to velocities data blocks. If the
 * frame is at the beginning of a frame set the time stamp of the frame set
 * is set.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data.
 * @param time is the time stamp of the frame (in seconds).
 * @param velocities is a 1D array of data to add. The array should be of length
 * n_particles * 3.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code frame_nr >= 0 \endcode The frame number to write must be >= 0.
 * @pre \code time >= 0 \endcode The time stamp must be >= 0.
 * @pre \code velocities != 0 \endcode The pointer to the velocities array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_with_time_write() and will
 * create a velocities data block if none exists. Velocities are stored as three
 * values per frame and compressed using TNG compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_vel_with_time_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const float *velocities);

/**
 * @brief High-level function for adding data to velocities data blocks at
 * double precision. If the frame is at the beginning of a frame set the
 * time stamp of the frame set is set.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data.
 * @param time is the time stamp of the frame (in seconds).
 * @param velocities is a 1D array of data to add. The array should be of length
 * n_particles * 3.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code frame_nr >= 0 \endcode The frame number to write must be >= 0.
 * @pre \code time >= 0 \endcode The time stamp must be >= 0.
 * @pre \code velocities != 0 \endcode The pointer to the velocities array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_with_time_double_write() and will
 * create a velocities data block if none exists. Velocities are stored as three
 * values per frame and compressed using TNG compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_vel_with_time_double_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const double *velocities);

/**
 * @brief High-level function for adding data to forces data blocks. If the
 * frame is at the beginning of a frame set the time stamp of the frame set
 * is set.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data.
 * @param time is the time stamp of the frame (in seconds).
 * @param forces is a 1D array of data to add. The array should be of length
 * n_particles * 3.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code frame_nr >= 0 \endcode The frame number to write must be >= 0.
 * @pre \code time >= 0 \endcode The time stamp must be >= 0.
 * @pre \code forces != 0 \endcode The pointer to the forces array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_with_time_write() and will
 * create a forces data block if none exists. Forces are stored as three
 * values per frame and compressed using gzip compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_force_with_time_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const float *forces);

/**
 * @brief High-level function for adding data to forces data blocks at
 * double precision. If the frame is at the beginning of a frame set
 * the time stamp of the frame set is set.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data.
 * @param time is the time stamp of the frame (in seconds).
 * @param forces is a 1D array of data to add. The array should be of length
 * n_particles * 3.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code frame_nr >= 0 \endcode The frame number to write must be >= 0.
 * @pre \code time >= 0 \endcode The time stamp must be >= 0.
 * @pre \code forces != 0 \endcode The pointer to the forces array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_with_time_double_write() and will
 * create a forces data block if none exists. Forces are stored as three
 * values per frame and compressed using gzip compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_force_with_time_double_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const double *forces);

/**
 * @brief High-level function for adding data to box shape data blocks. If the
 * frame is at the beginning of a frame set the time stamp of the frame set
 * is set.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data.
 * @param time is the time stamp of the frame (in seconds).
 * @param box_shape is a 1D array of data to add. The array should be of length 9.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code frame_nr >= 0 \endcode The frame number to write must be >= 0.
 * @pre \code time >= 0 \endcode The time stamp must be >= 0.
 * @pre \code box_shape != 0 \endcode The pointer to the box_shape array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_with_time_write() and will
 * create a box shape data block if none exists. Box shapes are stored as 9
 * values per frame and compressed using TNG compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_with_time_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const float *box_shape);

/**
 * @brief High-level function for adding data to box shape data blocks at
 * double precision. If the frame is at the beginning of a frame set the
 * time stamp of the frame set is set.
 * @param tng_data is the trajectory to use.
 * @param frame_nr is the frame number of the data.
 * @param time is the time stamp of the frame (in seconds).
 * @param box_shape is a 1D array of data to add. The array should be of length 9.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code frame_nr >= 0 \endcode The frame number to write must be >= 0.
 * @pre \code time >= 0 \endcode The time stamp must be >= 0.
 * @pre \code box_shape != 0 \endcode The pointer to the box_shape array must not
 * be a NULL pointer.
 * @details This function uses tng_util_generic_with_time_double_write() and will
 * create a box shape data block if none exists. Box shapes are stored as 9
 * values per frame and compressed using TNG compression.
 * N.b. Since compressed data is written a whole block at a time the data is not
 * actually written to disk until the frame set is finished or the TNG
 * trajectory is closed.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_box_shape_with_time_double_write
                (const tng_trajectory_t tng_data,
                 const int64_t frame_nr,
                 const double time,
                 const double *box_shape);

/**
 * @brief High-level function for getting the compression method and
 * multiplication factor of the last read frame of a specific data block.
 * @param tng_data is the trajectory to use.
 * @param block_id is the ID number of the block containing the data of
 * interest.
 * @param codec_id will be set to the value of the codec_id of the
 * compression of the data block. See tng_compression for more details.
 * @param factor will be set to the multiplication factor applied to
 * the values before compression, in order to get integers from them.
 * factor is 1/precision.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code codec_id != 0 \endcode  The pointer to the returned codec id
 * must not be a NULL pointer.
 * @pre \code factor != 0 \endcode The pointer to the returned multiplication
 * factor must not be a NULL pointer.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as invalid mode) or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_frame_current_compression_get
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 int64_t *codec_id,
                 double  *factor);

/** @brief High-level function for determining the next frame with data and what
 * data blocks have data for that frame. The search can be limited to certain
 * data blocks.
 * @param tng_data is the trajectory to use.
 * @param current_frame is the frame that was last read, from where to start
 * looking for data.
 * @param n_requested_data_block_ids is the number of data blocks listed in
 * requested_data_block_ids. If this is 0 all data blocks will be taken into
 * account.
 * @param requested_data_block_ids is an array of data blocks to look for.
 * @param next_frame will be set to the next frame with data.
 * @param n_data_blocks_in_next_frame is set to the number of data blocks with
 * data for next_frame.
 * @param data_block_ids_in_next_frame is set to an array (of length
 * n_data_blocks_in_next_frame) that lists the data block IDs with data for
 * next_frame. It must be pointing at NULL or previously allocated memory.
 * Memory for the array is allocated by this function.
 * The memory must be freed by the client afterwards or
 * there will be a memory leak.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code next_frame != 0 \endcode The pointer to the next frame must not
 * be NULL.
 * @pre \code n_data_blocks_in_next_frame != 0 \endcode The pointer to
 * n_data_blocks_in_next_frame must not be NULL.
 * @pre \code *data_block_ids_in_next_frame != 0 \endcode The pointer to the
 * list of data block IDs must not be NULL.
 * @pre \code n_requested_data_block_ids == 0 || requested_data_block_ids != 0 \endcode
 * If the number of requested data blocks != 0 then the array of data block IDs must not be NULL.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured or TNG_CRITICAL (2) if a major error
 * has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_next_frame_present_data_blocks_find
                (const tng_trajectory_t tng_data,
                 int64_t current_frame,
                 const int64_t n_requested_data_block_ids,
                 const int64_t *requested_data_block_ids,
                 int64_t *next_frame,
                 int64_t *n_data_blocks_in_next_frame,
                 int64_t **data_block_ids_in_next_frame);

/* @brief High-level function for getting all data block ids and their names
 * and stride lengths.
 * @param tng_data is the trajectory to use.
 * @param n_data_blocks is set to the number of data blocks in the trajectory.
 * @param data_block_ids is set to an array (of length
 * n_data_blocks) that lists the data block IDs in the trajectory.
 * It must be pointing at NULL or previously allocated memory.
 * Memory for the array is allocated by this function.
 * The memory must be freed by the client afterwards or
 * there will be a memory leak.
 * @param data_block_names is set to an array (of length
 * n_data_blocks) that contains the names of the data blocks.
 * It must be pointing at NULL or previously allocated memory.
 * Memory for the array is allocated by this function.
 * The memory must be freed by the client afterwards or
 * there will be a memory leak.
 * @param stride_lengths is set to an array (of length
 * n_data_blocks) that lists the stride lengths of the data blocks.
 * It must be pointing at NULL or previously allocated memory.
 * Memory for the array is allocated by this function.
 * The memory must be freed by the client afterwards or
 * there will be a memory leak.
 * @param n_values_per_frame is set to an array (of length
 * n_data_blocks) that lists the number of values per frame of the data blocks.
 * It must be pointing at NULL or previously allocated memory.
 * Memory for the array is allocated by this function.
 * The memory must be freed by the client afterwards or
 * there will be a memory leak.
 * @param block_types is set to an array (of length
 * n_data_blocks) that lists the block types of the data blocks.
 * It must be pointing at NULL or previously allocated memory.
 * Memory for the array is allocated by this function.
 * The memory must be freed by the client afterwards or
 * there will be a memory leak.
 * @param dependencies is set to an array (of length
 * n_data_blocks) that lists the dependencies of the data blocks.
 * It must be pointing at NULL or previously allocated memory.
 * Memory for the array is allocated by this function.
 * The memory must be freed by the client afterwards or
 * there will be a memory leak.
 * @param compressions is set to an array (of length
 * n_data_blocks) that lists the compressions of the data blocks.
 * It must be pointing at NULL or previously allocated memory.
 * Memory for the array is allocated by this function.
 * The memory must be freed by the client afterwards or
 * there will be a memory leak.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code n_data_blocks != 0 \endcode The pointer to
 * n_data_blocks must not be NULL.
 * @pre \code data_block_ids != 0 \endcode The pointer to the
 * list of data block IDs must not be NULL.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured or TNG_CRITICAL (2) if a major error
 * has occured.
 */
/*
tng_function_status DECLSPECDLLEXPORT tng_util_trajectory_all_data_block_types_get
                (const tng_trajectory_t tng_data,
                 int64_t *n_data_blocks,
                 int64_t **data_block_ids,
                 char ***data_block_names,
                 int64_t **stride_lengths,
                 int64_t **n_values_per_frame,
                 char **block_types,
                 char **dependencies,
                 char **compressions);
*/

/** @brief Finds the frame set of the specified frame in order to prepare for writing
 * after it.
 * @param tng_data is the trajectory to use.
 * @param prev_frame is the frame after which to start appending.
 * @pre \code tng_data != 0 \endcode The trajectory container (tng_data)
 * must be initialised before using it.
 * @pre \code prev_frame >= 0 \endcode The previous frame must not be negative.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occured (such as not finding the requested frame) or TNG_CRITICAL (2)
 * if a major error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_prepare_append_after_frame
                (const tng_trajectory_t tng_data,
                 const int64_t prev_frame);


/** @brief Get the number of frames containing data of a specific type.
 * @param tng_data is the trajectory to use.
 * @param block_id is the id of the block of the data type.
 * @param n_frames is set to the number of frames containing data of
 * the requested data type.
 * @return TNG_SUCCESS (0) if successful or TNG_CRITICAL (2) if a major
 * error has occured.
 */
tng_function_status DECLSPECDLLEXPORT tng_util_num_frames_with_data_of_block_id_get
                (const tng_trajectory_t tng_data,
                 const int64_t block_id,
                 int64_t *n_frames);
/** @} */ /* end of group2 */


#ifdef __cplusplus
}  /* end extern "C" */
#endif

#endif /* TNG_IO_H */

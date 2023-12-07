/*
 * This code is part of the tng binary trajectory format.
 *
 * Copyright (c) 2012-2013, The GROMACS development team.
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

#ifndef TNG_IO_FWD_H
#define TNG_IO_FWD_H 1

/** A pointer to the main trajectory data storage. */
typedef struct tng_trajectory* tng_trajectory_t;
/** A pointer to a molecule description. */
typedef struct tng_molecule* tng_molecule_t;
/** A pointer to a molecular chain description. */
typedef struct tng_chain* tng_chain_t;
/** A pointer to a molecular residue description. */
typedef struct tng_residue* tng_residue_t;
/** A pointer to a molecular atom description. */
typedef struct tng_atom* tng_atom_t;
/** A pointer to a bond between two atoms. */
typedef struct tng_bond* tng_bond_t;
/** A pointer to a structure containing data common to all trajectory blocks,
 *  such as header and contents. */
typedef struct tng_gen_block* tng_gen_block_t;
/** A pointer to particle mapping information. */
typedef struct tng_particle_mapping* tng_particle_mapping_t;
/** A pointer to a structure containing frame set information. */
typedef struct tng_trajectory_frame_set* tng_trajectory_frame_set_t;
/** A pointer to a data container. */
typedef struct tng_data* tng_data_t;

#endif

/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Should be called after generating or reading mtop,
 * to set some compute intesive variables to avoid
 * N^2 operations later on.
 */
void
gmx_mtop_finalize(gmx_mtop_t *mtop);


/* Returns the total number of charge groups in mtop */
int
ncg_mtop(const gmx_mtop_t *mtop);


/* Returns a pointer to the t_atom struct belonging to atnr_global.
 * This can be an expensive operation, so if possible use
 * one of the atom loop constructs below.
 */
void
gmx_mtop_atomnr_to_atom(const gmx_mtop_t *mtop,int atnr_global,
			t_atom **atom);


/* Returns a pointer to the molecule interaction array ilist_mol[F_NRE]
 * and the local atom number in the molecule belonging to atnr_global.
 */
void
gmx_mtop_atomnr_to_ilist(const gmx_mtop_t *mtop,int atnr_global,
			 t_ilist **ilist_mol,int *atnr_offset);


/* Returns the molecule block index
 * and the molecule number in the block
 * and the atom number offset for the atom indices in moltype
 * belonging to atnr_global.
 */
void
gmx_mtop_atomnr_to_molblock_ind(const gmx_mtop_t *mtop,int atnr_global,
				int *molb,int *molnr,int *atnr_mol);


/* Returns atom name, global resnr and residue name  of atom atnr_global */
void
gmx_mtop_atominfo_global(const gmx_mtop_t *mtop,int atnr_global,
			 char **atomname,int *resnr,char **resname);


/* Abstract type for atom loop over all atoms */
typedef struct gmx_mtop_atomloop_all *gmx_mtop_atomloop_all_t;

/* Initialize an atom loop over all atoms in the system.
 * The order of the atoms will be as in the state struct.
 * Only use this when you really need to loop over all atoms,
 * i.e. when you use groups which might differ per molecule,
 * otherwise use gmx_mtop_atomloop_block.
 */
gmx_mtop_atomloop_all_t
gmx_mtop_atomloop_all_init(const gmx_mtop_t *mtop);

/* Loop to the next atom.
 * When not at the end:
 *   returns TRUE and at_global,
 *   writes the global atom number in *at_global
 *   and sets the pointer atom to the t_atom struct of that atom.
 * When at the end, destroys aloop and returns FALSE.
 * Use as:
 * gmx_mtop_atomloop_all_t aloop;
 * aloop = gmx_mtop_atomloop_all_init(mtop)
 * while (gmx_mtop_atomloop_all_next(aloop,&at_global,&atom)) {
 *     ...
 * }
 */
gmx_bool
gmx_mtop_atomloop_all_next(gmx_mtop_atomloop_all_t aloop,
			   int *at_global,t_atom **atom);

/* Return the atomname, the residue number and residue name
 * of the current atom in the loop.
 */
void
gmx_mtop_atomloop_all_names(gmx_mtop_atomloop_all_t aloop,
			    char **atomname,int *resnr,char **resname);

/* Return the a pointer to the moltype struct of the current atom
 * in the loop and the atom number in the molecule.
 */
void
gmx_mtop_atomloop_all_moltype(gmx_mtop_atomloop_all_t aloop,
			      gmx_moltype_t **moltype,int *at_mol);


/* Abstract type for atom loop over atoms in all molecule blocks */
typedef struct gmx_mtop_atomloop_block *gmx_mtop_atomloop_block_t;

/* Initialize an atom loop over atoms in all molecule blocks the system.
 */
gmx_mtop_atomloop_block_t
gmx_mtop_atomloop_block_init(const gmx_mtop_t *mtop);

/* Loop to the next atom.
 * When not at the end:
 *   returns TRUE
 *   sets the pointer atom to the t_atom struct of that atom
 *   and return the number of molecules corresponding to this atom.
 * When at the end, destroys aloop and returns FALSE.
 * Use as:
 * gmx_mtop_atomloop_block_t aloop;
 * aloop = gmx_mtop_atomloop_block_init(mtop)
 * while (gmx_mtop_atomloop_block_next(aloop,&atom,&nmol)) {
 *     ...
 * }
 */
gmx_bool
gmx_mtop_atomloop_block_next(gmx_mtop_atomloop_block_t aloop,
			     t_atom **atom,int *nmol);


/* Abstract type for ilist loop over all ilists */
typedef struct gmx_mtop_ilistloop *gmx_mtop_ilistloop_t;

/* Initialize an ilist loop over all molecule types in the system. */
gmx_mtop_ilistloop_t
gmx_mtop_ilistloop_init(const gmx_mtop_t *mtop);


/* Loop to the next molecule,
 * When not at the end:
 *   returns TRUE and a pointer to the next array ilist_mol[F_NRE],
 *   writes the number of molecules for this ilist in *nmol.
 * When at the end, destroys iloop and returns FALSE.
 */
gmx_bool
gmx_mtop_ilistloop_next(gmx_mtop_ilistloop_t iloop,
			t_ilist **ilist_mol,int *nmol);


/* Abstract type for ilist loop over all ilists of all molecules */
typedef struct gmx_mtop_ilistloop_all *gmx_mtop_ilistloop_all_t;

/* Initialize an ilist loop over all molecule types in the system.
 * Only use this when you really need to loop over all molecules,
 * i.e. when you use groups which might differ per molecule,
 * otherwise use gmx_mtop_ilistloop.
 */
gmx_mtop_ilistloop_all_t
gmx_mtop_ilistloop_all_init(const gmx_mtop_t *mtop);

/* Loop to the next molecule,
 * When not at the end:
 *   returns TRUE and a pointer to the next array ilist_mol[F_NRE],
 *   writes the atom offset which should be added to iatoms in atnr_offset.
 * When at the end, destroys iloop and returns FALSE.
 */
gmx_bool
gmx_mtop_ilistloop_all_next(gmx_mtop_ilistloop_all_t iloop,
			    t_ilist **ilist_mol,int *atnr_offset);


/* Returns the total number of interactions in the system of type ftype */
int
gmx_mtop_ftype_count(const gmx_mtop_t *mtop,int ftype);


/* Returns a charge group index for the whole system */
t_block
gmx_mtop_global_cgs(const gmx_mtop_t *mtop);


/* Returns a single t_atoms struct for the whole system */ 
t_atoms
gmx_mtop_global_atoms(const gmx_mtop_t *mtop);


/* Make all charge groups the size of one atom.
 * When bKeepSingleMolCG==TRUE keep charge groups for molecules
 * that consist of a single charge group.
 */
void
gmx_mtop_make_atomic_charge_groups(gmx_mtop_t *mtop,gmx_bool bKeepSingleMolCG);


/* Generate a 'local' topology for the whole system.
 * When ir!=NULL the free energy interactions will be sorted to the end.
 */
gmx_localtop_t *
gmx_mtop_generate_local_top(const gmx_mtop_t *mtop,const t_inputrec *ir);


/* Converts a gmx_mtop_t struct to t_topology.
 * All memory relating only to mtop will be freed.
 */
t_topology
gmx_mtop_t_to_t_topology(gmx_mtop_t *mtop);

#ifdef __cplusplus
}
#endif


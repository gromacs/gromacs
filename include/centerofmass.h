/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
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
 */
/*! \file
 * \brief API for calculation of centers of mass/geometry.
 *
 * This header defines a few functions that can be used to calculate
 * centers of mass/geometry for a group of atoms.
 * These routines can be used independently of the other parts of the
 * library, but they are also used internally by the selection engine.
 * In most cases, it should not be necessary to call these functions
 * directly.
 * Instead, one should write an analysis tool such that it gets all
 * positions through selections.
 *
 * The functions in the header can be divided into a few groups based on the
 * parameters they take. The simplest group of functions calculates the center
 * of a single group of atoms:
 *  - gmx_calc_cog(): Calculates the center of geometry (COG) of a given
 *    group of atoms.
 *  - gmx_calc_com(): Calculates the center of mass (COM) of a given group
 *    of atoms.
 *  - gmx_calc_comg(): Calculates either the COM or COG, based on a
 *    gmx_boolean flag.
 *
 * A second set of routines is provided for calculating the centers for groups
 * that wrap over periodic boundaries (gmx_calc_cog_pbc(), gmx_calc_com_pbc(),
 * gmx_calc_comg_pbc()). These functions are slower, because they need to
 * adjust the center iteratively.
 *
 * It is also possible to calculate centers for several groups of atoms in
 * one call. The functions gmx_calc_cog_block(), gmx_calc_com_block() and
 * gmx_calc_comg_block() take an index group and a partitioning of that index
 * group (as a \c t_block structure), and calculate the centers for
 * each group defined by the \c t_block structure separately.
 *
 * Finally, there is a function gmx_calc_comg_blocka() that takes both the
 * index group and the partitioning as a single \c t_blocka structure.
 */
#ifndef CENTEROFMASS_H
#define CENTEROFMASS_H

#include "typedefs.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** Calculate a single center of geometry. */
int
gmx_calc_cog(t_topology *top, rvec x[], int nrefat, atom_id index[], rvec xout);
/** Calculate a single center of mass. */
int
gmx_calc_com(t_topology *top, rvec x[], int nrefat, atom_id index[], rvec xout);
/** Calculate force on a single center of geometry. */
int
gmx_calc_cog_f(t_topology *top, rvec f[], int nrefat, atom_id index[], rvec fout);
/** Calculate a single center of mass/geometry. */
int
gmx_calc_comg(t_topology *top, rvec x[], int nrefat, atom_id index[],
              gmx_bool bMass, rvec xout);
/** Calculate force on a single center of mass/geometry. */
int
gmx_calc_comg_f(t_topology *top, rvec f[], int nrefat, atom_id index[],
                gmx_bool bMass, rvec fout);

/** Calculate a single center of geometry iteratively, taking PBC into account. */
int
gmx_calc_cog_pbc(t_topology *top, rvec x[], t_pbc *pbc,
                 int nrefat, atom_id index[], rvec xout);
/** Calculate a single center of mass iteratively, taking PBC into account. */
int
gmx_calc_com_pbc(t_topology *top, rvec x[], t_pbc *pbc,
                 int nrefat, atom_id index[], rvec xout);
/** Calculate a single center of mass/geometry iteratively with PBC. */
int
gmx_calc_comg_pbc(t_topology *top, rvec x[], t_pbc *pbc,
                  int nrefat, atom_id index[], gmx_bool bMass, rvec xout);

/** Calculate centers of geometry for a blocked index. */
int
gmx_calc_cog_block(t_topology *top, rvec x[], t_block *block,
                   atom_id index[], rvec xout[]);
/** Calculate centers of mass for a blocked index. */
int
gmx_calc_com_block(t_topology *top, rvec x[], t_block *block,
                   atom_id index[], rvec xout[]);
/** Calculate forces on centers of geometry for a blocked index. */
int
gmx_calc_cog_f_block(t_topology *top, rvec f[], t_block *block,
                     atom_id index[], rvec fout[]);
/** Calculate centers of mass/geometry for a blocked index. */
int
gmx_calc_comg_block(t_topology *top, rvec x[], t_block *block,
                    atom_id index[], gmx_bool bMass, rvec xout[]);
/** Calculate forces on centers of mass/geometry for a blocked index. */
int
gmx_calc_comg_f_block(t_topology *top, rvec f[], t_block *block,
                      atom_id index[], gmx_bool bMass, rvec fout[]);
/** Calculate centers of mass/geometry for a set of blocks; */
int
gmx_calc_comg_blocka(t_topology *top, rvec x[], t_blocka *block,
                     gmx_bool bMass, rvec xout[]);
/** Calculate forces on centers of mass/geometry for a set of blocks; */
int
gmx_calc_comg_f_blocka(t_topology *top, rvec x[], t_blocka *block,
                       gmx_bool bMass, rvec xout[]);

#ifdef __cplusplus
}
#endif

#endif

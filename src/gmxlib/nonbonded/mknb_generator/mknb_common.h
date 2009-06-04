/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifndef _MKNB_H_
#define _MKNB_H_



/*! \file  mknb_common.h
 *  \brief Kernel generator (only for compile): Typedefs/variables
 *
 *  \internal
 *
 *  \note This file is only used to generate the inner loop kernels
 *        at compile time, which in turn are included in Gromacs. 
 *        This code itself is NOT linked into the Gromacs library, so
 *        it does not need to be threadsafe.
 *
 *  mknb_common.h contains definitions of enumerated types and structures
 *  used in the nonbonded kernel generator itself, i.e. not in
 *  the generated code. Declarations for the generated code are
 *  in the file mknb_declarations.c
 */



/*! \brief Kernel generator (only for compile): Options from command line.
 * 
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  The options in this structure control the overall behavior of 
 *  the kernel generator, for instance whether threads and/or software
 *  inverse square root is used.
 *
 *  For each alternative, 1 means enabled, and 0 disabled.
 */
struct mknb_options
{
	int     threads;          /*!< Create kernels with thread support  */
	int     software_invsqrt; /*!< Use software optimization of 1/sqrt */
	int     ppc_invsqrt;      /*!< Use IBM PowerPC intrinsics of 1/sqrt */
	int     prefetch_forces;  /*!< Prefetch forces in innermost loop   */
};




/*! \brief Kernel generator (only for compile): Coulomb interaction alternatives
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  Type of Coulomb interaction to generate in a 
 *  nonbonded kernel. There are different kernels for different
 *  alternatives, so the generator code will be looping over these.
 *
 *  Make sure you update the text descriptions in mknb_coul_names
 *  if you change this.
 */ 
enum mknb_coul {
	MKNB_COUL_NO,        /*!< No Coulomb interaction           */
	MKNB_COUL_NORMAL,    /*!< Standard 1/r Coulomb interaction */
	MKNB_COUL_RF,        /*!< Reaction-Field Coulomb           */
	MKNB_COUL_TAB,       /*!< Tabulated Coulomb                */
	MKNB_COUL_GB,        /*!< Generalized Born Coulomb         */
	MKNB_COUL_NR         /*!< Number of choices for Coulomb    */
};




/*! \brief Kernel generator (only for compile): Van der Waals interaction alternatives
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  Type of Van der Waals interaction to generate in a 
 *  nonbonded kernel. There are different kernels for different
 *  alternatives, so the generator code will be looping over these.
 *
 *  Make sure you update the text descriptions in mknb_vdw_names
 *  if you change this.
 */ 
enum mknb_vdw {
	MKNB_VDW_NO,         /*!< No Van der Waals interaction     */
	MKNB_VDW_LJ,         /*!< Lennard-Jones 6-12 interactions  */
	MKNB_VDW_BHAM,       /*!< Buckingham 6-exp interactions    */
	MKNB_VDW_TAB,        /*!< Tabulated VdW interactions       */
	MKNB_VDW_NR          /*!< Number of choices for Vdw        */
};



/*! \brief Kernel generator (only for compile): Water optimization alternatives
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  Since water is extremely common in biological systems, Gromacs 
 *  includes special nonbonded kernels optimized for interactions between 
 *  water molecules and other atoms, as well as interactions between
 *  pairs of water molecules to improve performance.
 *
 *  Make sure you update the text descriptions in mknb_water_names
 *  if you change this.
 */ 
enum mknb_water {
	MKNB_WATER_NO,           /*!< No water optimization                    */
	MKNB_WATER_SPC_SINGLE,   /*!< 3-atom water - other atom (SPC,TIP3P)    */
	MKNB_WATER_SPC_PAIR,     /*!< 3-atom water - 3-atom water  (SPC,TIP3P) */
	MKNB_WATER_TIP4P_SINGLE, /*!< 4-atom water - other atom (TIP4P)        */
	MKNB_WATER_TIP4P_PAIR,   /*!< 4-atom water - 4-atom water (TIP4P)      */
	MKNB_WATER_NR            /*!< Number of choices for water optimization */
};




/*! \brief Kernel generator (only for compile): Text for enum mknb_coul
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  Make sure you update these when you change the enumerated
 *  alternatives in enum mknb_coul. 
 */
extern const char *
mknb_coul_names[MKNB_COUL_NR];




/*! \brief Kernel generator (only for compile): Text for enum mknb_vdw 
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  Make sure you update these when you change the enumerated
 *  alternatives in enum mknb_vdw.
 */
extern const char *
mknb_vdw_names[MKNB_VDW_NR];




/*! \brief Kernel generator (only for compile): Text for enum mknb_water
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  Make sure you update these when you change the enumerated
 *  alternatives in enum mknb_water.
 */
extern const char *
mknb_water_names[MKNB_WATER_NR];




/*! \brief Kernel generator (only for compile): Options for a nonbonded kernel
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  Some fields deserve a comment:
 * 
 *  ni,nj: *  For a normal (non-water) kernel we treat one atom
 *  per iteration in both the outer and inner loop. Then ni=nj=1.
 *  In water-optimized kernels we use 3 (SPC/TIP3P) or 4 (TIP4P)
 *  atoms in the outer loop (ni=3/4, nj=1). Finally, in kernels
 *  optimized for water-water interactions you have multiple 
 *  atoms in the inner loop too, i.e. ni=nj=3 (or 4).
 *
 *  table_element_size: Each interaction uses four floating-point
 *  variables per table point. To optimize, we merge the Coulomb
 *  and VdW tables into a single structure if both are tabulated,
 *  i.e. depending on the type of interaction the kernel is called
 *  either with a table containing only Coulomb, only VdW, or both.
 *  I.e., table_element_size is 4 if only Coulomb is tabulated,
 *  8 if only VdW is tabulated, and 12 if both are tabulated.
 * 
 *  nvdw_parameters: Lennard-Jones interactions use 2 (c6 & c12),
 *  while Buckingham interactions need 3 (c6 & a & b).
 */
struct mknb_func {
	char             name[32];          /*!< Name for this kernel         */
	enum mknb_coul   coul;              /*!< Type of Coulomb interaction  */
	enum mknb_vdw    vdw;               /*!< Type of VdW interaction      */
	enum mknb_water  water;             /*!< Water optimization type      */
	int              ni;                /*!< Number of atoms in outerloop */
    int              nj;                /*!< Number of atoms in innerloop */
	int              do_force;          /*!< Calculate Forces in loop?    */
	int              table_element_size;/*!< Number of floats/table point */
	int              nvdw_parameters;   /*!< LJ needs 2, Buckingham 3     */
};




/*! \brief Kernel generator (only for compile): Global variable for general mknb options.
 *
 *  \internal
 *
 *  Variable is defined in mknb_common.c
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  The fields in this structure are set when the program is 
 *  started, and specifies general options like thread support, software
 *  inverse square root, force prefetching.
 */
extern struct mknb_options
mknb_options;    




/*! \brief Kernel generator (only for compile):  Global variable for current kernel options.
 *
 *  \internal
 *
 *  Variable is defined in mknb_common.c
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  The fields in this structure are updated before each new kernel is
 *  generated, and determine the type of interactions, water optimization,
 *  etc.
 */
extern struct mknb_func     
mknb_func;






#endif /* _MKNB_H_ */

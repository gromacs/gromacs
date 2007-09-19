/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.2
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2007, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#include "typedefs.h"

#define FLAG_DOTS       01
#define FLAG_VOLUME     02
#define FLAG_ATOM_AREA  04

#define NSC nsc_dclm

extern int NSC(
  real * ,	/* atom coordinates xyz0, xyz1, ... */
  real * ,	/* atom radii r0, r1, r2, ... */
  int ,         /* number of atoms */
  int ,		/* number of dots per fully accessible sphere */
  int ,		/* flag : dots, volume and/or area per atom */
  real * ,	/* 1 output: overall area */
  real ** ,     /* 2 output: pointer to list of areas per atom */
  real * ,	/* 3 output: overall volume */
  real ** ,	/* 4 output: pointer to list of surface dots x0, y0, z0, ... */
  int *	        /* 5 output: number of surface dots */
  );

extern int nsc_dclm2(rvec *coords, real *radius, int nat, atom_id index[],
		     int  densit, int mode,
		     real *value_of_area, real **at_area,
		     real *value_of_vol,
		     real **lidots, int *nu_dots,
		     matrix box);

extern int nsc_dclm_pbc(rvec *coords, real *radius, int nat,
			int  densit, int mode,
			real *value_of_area, real **at_area,
			real *value_of_vol,
			real **lidots, int *nu_dots,
			atom_id index[],matrix box);

/* 
    User notes :
The input requirements :
  The arrays with atom coordinates and radii are thought to start
  with index 0, i.e., places 0, 1, and 2 are the x-, y-, and z-
  coordinates of the zero-th atom and place 0 in the other array
  is its radius.

  PLEASE TAKE INTO ACCOUNT THAT THE RADII GIVEN HERE ARE DIRECTLY
  USED FOR SURFACE CALCULATION. NSC does not increment with a probe 
  radius.

  The user can define any number of dots. The program selects a
  dot density that is the lowest possible with at least the required
  number of dots. The points are distributed in accordance with the
  icosahedron-based or the dodecahedron-based method as described in
  ref. 1.

The output requirements are :
  1 and 3 :  pointer to an existing real
  2 and 4 :  pointer to an existing pointer to real
             NSC allocates memory for an array
  5       :  pointer to an existing integer

The subroutine NSC makes use of variant 2 described in reference 1.
By selecting the necessary output via flags, the requirements for
cpu-time and computer memory can be adapted to the actual needs.

Example : flag = FLAG_VOLUME | FLAG_ATOM_AREA | FLAG_DOTS
          The routine calculates the area, volume and the dot surface. The
          program allocates arrays for the atomwise areas and for the surface
          dots. The addresses are returned in the pointers to pointers to
          real. 
          This variant is not recommended because normally the dot surface
          is needed for low point density (e.g.42) at which area and volume 
          are inaccurate. The sign "|" is used as binary AND !

          flag = FLAG_VOLUME | FLAG_ATOM_AREA
          In this case the large arrays for storing the surface dots
          are not allocated. A large point number of the fully accessible
          sphere can be selected. Good accuracy is already achieved with
          600-700 points per sphere (accuracy of about 1.5 square Angstrem
          per atomic sphere).
          Output pointers 4 and 5 may be NULL.

          flag = FLAG_DOTS
          Only the dot surface is produced.
          Output pointers 2 and 3 may be NULL.

The output pointer 1 cannot be set to NULL in any circumstances. The
overall area value is returned in every mode.

All files calling NSC should include nsc.h !!


Example for calling NSC (contents of user file): 

  ...
#include "nsc.h"

int routine_calling_NSC(int n_atom, real *coordinates, real *radii) {
  real area, volume, *atomwise_area, *surface_dots;
  int    i, density = 300, n_dots;

  ...

  for (i=0; i<n_atom; i++) {
   radii[i]  += 1.4      /# add the probe radius if necessary #/

  if (NSC(coordinates, radii, n_atom, density,
          FLAG_AREA | FLAG_VOLUME | FLAG_DOTS,
          &area, &atomwise_area, &volume, &surface_dots, &n_dots))
    printf("error occured\n");
    return 1;
    }

  ...

/# do something with areas, volume and surface dots #/

  ...

  return 0;
  }

*/

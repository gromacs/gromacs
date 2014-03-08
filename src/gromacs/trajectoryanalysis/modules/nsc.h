/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "typedefs.h"

#define FLAG_DOTS       01
#define FLAG_VOLUME     02
#define FLAG_ATOM_AREA  04

#ifdef __cplusplus
extern "C"
{
#endif

int nsc_dclm_pbc(const rvec *coords, real *radius, int nat,
                 int  densit, int mode,
                 real *value_of_area, real **at_area,
                 real *value_of_vol,
                 real **lidots, int *nu_dots,
                 atom_id index[], int ePBC, matrix box);

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

#ifdef __cplusplus
}
#endif

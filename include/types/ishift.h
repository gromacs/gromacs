/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#define D_BOX_Z 1
#ifdef ALLOW_OFFDIAG_LT_HALFDIAG
#define D_BOX_Y 2
#define D_BOX_X 2
#else
#define D_BOX_Y 1
#define D_BOX_X 1
#endif 
#define N_BOX_Z (2*D_BOX_Z+1)
#define N_BOX_Y (2*D_BOX_Y+1)
#define N_BOX_X (2*D_BOX_X+1)
#define N_IVEC  (N_BOX_Z*N_BOX_Y*N_BOX_X)
#define CENTRAL (N_IVEC/2)
#define SHIFTS  N_IVEC

#define XYZ2IS(x,y,z) (N_BOX_X*(N_BOX_Y*((z)+D_BOX_Z)+(y)+D_BOX_Y)+(x)+D_BOX_X)
#define IVEC2IS(iv)   (XYZ2IS((iv)[XX],(iv)[YY],(iv)[ZZ]))

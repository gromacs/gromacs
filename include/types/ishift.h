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

#define D_BOX   1
#define N_BOX   (2*D_BOX+1)
#define N_BOX2  (N_BOX*N_BOX)
#define N_IVEC  (N_BOX*N_BOX*N_BOX)
#define CENTRAL (N_IVEC/2)
#define SHIFTS  N_IVEC

#define XYZ2IS(x,y,z) (N_BOX2*((x)+D_BOX)+N_BOX*((y)+D_BOX)+(z+D_BOX))
#define IVEC2IS(iv)   (XYZ2IS((iv)[XX],(iv)[YY],(iv)[ZZ]))

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

#ifndef _fftw_wrapper_h
#define _fftw_wrapper_h

static char *SRCID_fftw_wrapper_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef FFTW_NAME_FFTW
#  include<fftw.h>
#  include<rfftw.h>
#  ifdef USE_MPI
#    include<fftw_mpi.h>
#    include<rfftw_mpi.h>
#  endif
#elif defined FFTW_NAME_SFFTW
#  include<sfftw.h>
#  include<srfftw.h>
#  ifdef USE_MPI
#    include<sfftw_mpi.h>
#    include<srfftw_mpi.h>
#  endif
#elif defined FFTW_NAME_DFFTW
#  include<dfftw.h>
#  include<drfftw.h>
#  ifdef USE_MPI
#    include<dfftw_mpi.h>
#    include<drfftw_mpi.h>
#  endif
#endif

#endif			

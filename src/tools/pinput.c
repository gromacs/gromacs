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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_pinput_c = "$Id$";

#include "paramio.h"
#include "pinput.h"

#define PPIN \
  ITYPE("nframes",              p->nframes,        1) \
  ETYPE("optimize",             p->nSel)              \
        "radius", "twist", "rise", "len", "nhx", "dipole", "rms", "cphi", NULL, \
  ETYPE("funct",                p->funct)              \
        "MC", "RECOMB", "PROPTRJ", NULL, \
  ITYPE("nsteps",               p->nsteps,         100) \
  ITYPE("nev",                  p->nev,            10) \
  ITYPE("nskip",                p->nskip,          0) \
  STYPE("projection",           p->base,           "WEDPRJVEC10.DAT") \
  STYPE("recomb",               p->recomb,         "WEDRECOMB10.DAT") \
  STYPE("gamma",                p->gamma,          "WEDGAMMA10.DAT") \
  RTYPE("stepsize",             p->step,           0.1) \
  RTYPE("tolerance",            p->tol,            1e-6) \
  RTYPE("ref-fluc",             p->v0,             1e-3) \
  NULL

void read_inp(char *fnin,char *fnout,t_pinp *p)
{
  read_params(fnin,PPIN);
  write_params(fnout,PPIN);
}

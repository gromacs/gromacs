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
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_general_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef __GENERAL__
#define __GENERAL__

/*------------------------------------------------------------------------
 *  edition history
 *------------------------------------------------------------------------
 *
 *  version  date    author  description
 *  ----------------------------------------------------------------------
 *  1.1      07Dec93 mkr     creation
 */

/*------------------------------------------------------------------------
 *  type definitions
 *------------------------------------------------------------------------
 */
#define LONG	unsigned int			/* 32 bits						*/
#define WORD	unsigned short			/* 16 bits						*/
#define BYTE	unsigned char			/*  8 bits						*/

/*------------------------------------------------------------------------
 *  constant definitions
 *------------------------------------------------------------------------
 */
#define	FALSE	0
#define	TRUE	1

#endif  /* __GENERAL__ */

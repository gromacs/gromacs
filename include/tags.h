/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * S  C  A  M  O  R  G
 */

#ifndef	_tags_h
#define	_tags_h

#ifdef HAVE_IDENT
#ident	"@(#) tags.h 1.3 11/23/92"
#endif /* HAVE_IDENT */

#define SYSCALL_TAG	0x11		/* Tag for server system calls      */
#define SEMGET_TAG	0x10		/* Server subcommand i860_semget()  */
#define SEMCTL_TAG	0x11		/* Server subcommand i860_semctl()  */
#define SEMOP_TAG	0x12		/* Server subcommand i860_semop()   */
#define SYNCALL_TAG	0x13		/* Server subcommand i860_syncall() */

#endif	/* _tags_h */

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

#ifndef	_systest_h
#define	_systest_h

#ifdef HAVE_IDENT
#ident	"@(#) systest.h 1.3 11/23/92"
#endif /* HAVE_IDENT */

#ifndef MASTER_ID
#define MASTER_ID	2
#endif
#ifndef BUFSIZE
#define BUFSIZE		25000
#endif
#define MEMMARGIN	10000

/*
 * Led usage:
 */
 
#define PARM_LED		0
#define	ACTIVE_LED		1
#define	HOP_LED			2
#define	ERR_LED			3
#define RANDOM_TEST_LED		0
#define WALKBIT_TEST_LED	1
#define WORD_TEST_LED		2
#define RANDOM_SWITCH		0xffff
#define WALKBIT_SWITCH		0x7ffff
#define WORD_SWITCH		0x7ffff

#endif	/* _systest_h */

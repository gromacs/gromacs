/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_index_h = "$Id$";

extern t_block *new_block(void);
/* allocate new block */

extern void write_index(char *outf, t_block *b,char **gnames);
/* Writes index blocks to outf (writes an indexfile) */

void add_grp(t_block *b,char ***gnames,int nra,atom_id a[],char *name);
/* Ads group a with name name to block b and namelist gnames */ 

extern void analyse(t_atoms *atoms,t_block *gb,char ***gn,
                    bool bASK,bool bVerb);
/* Makes index groups gb with names gn for atoms in atoms.
 * bASK=FALSE gives default groups.
 */


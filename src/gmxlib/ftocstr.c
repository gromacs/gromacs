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
static char *SRCID_ftocstr_c = "$Id$";

int ftocstr(char *ds, int dl, char *ss, int sl)
    /* dst, src ptrs */
    /* dst max len */
    /* src len */
{
    char *p;

    p = ss + sl;
    while ( --p >= ss && *p == ' ' );
    sl = p - ss + 1;
    dl--;
    ds[0] = 0;
    if (sl > dl)
      return 1;
    while (sl--)
      (*ds++ = *ss++);
    *ds = '\0';
    return 0;
}


int ctofstr(char *ds, int dl, char *ss)
     /* dest space */
     /* max dest length */
     /* src string (0-term) */
{
    while (dl && *ss) {
	*ds++ = *ss++;
	dl--;
    }
    while (dl--)
	*ds++ = ' ';
    return 0;
}

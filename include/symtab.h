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
 * GRowing Old MAkes el Chrono Sweat
 */

#ifndef _symtab_h
#define _symtab_h

static char *SRCID_symtab_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) symtab.h 1.6 12/16/92"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "typedefs.h"

/*
 * This module handles symbol table manipulation. All text strings 
 * needed by an application are allocated only once. All references
 * to these text strings use handles returned from the put_symtab()
 * routine. These handles can easily be converted to address independent
 * values by invoking lookup_symtab(). So when writing structures to
 * a file which contains text strings, this value can be written in stead
 * of the text string or its address. This value can easily be converted
 * back to a text string handle by get_symtab_handle().
 */

extern void open_symtab(t_symtab *symtab);
     /* Initialises the symbol table symtab.
      */

extern void close_symtab(t_symtab *symtab);
     /* Undoes the effect of open_symtab(), after invoking this function, 
      * no value can be added to the symbol table, only values can be 
      * retrieved using get_symtab().
      */

extern void free_symtab(t_symtab *symtab);
     /* Frees the space allocated by the symbol table itself */

extern void done_symtab(t_symtab *symtab);
     /* Frees the space allocated by the symbol table, including all
      * entries in it */

extern char **put_symtab(t_symtab *symtab,char *name);
     /* Enters a string into the symbol table symtab, if it was not
      * available, a reference to a copy is returned else a reference 
      * to the earlier entered value is returned. Strings are trimmed
      * of spaces.
      */

extern int lookup_symtab(t_symtab *symtab,char **name);
     /* Returns a unique handle for **name, without a memory reference.
      * It is a failure when name cannot be found in the symbol table,
      * it should be entered before with put_symtab().
      */

extern char **get_symtab_handle(t_symtab *symtab,int name);
     /* Returns a text string handle for name. Name should be a value
      * returned from lookup_symtab(). So get_symtab_handle() and 
      * lookup_symtab() are inverse functions.
      */

extern long wr_symtab(FILE *fp,t_symtab *symtab);
     /* Writes the symbol table symtab to the file, specified by fp.
      * The function returns the number of bytes written.
      */

extern long rd_symtab(FILE *fp,t_symtab *symtab);
     /* Reads the symbol table symtab from the file, specified by fp.
      * This will include allocating the needed space. The function 
      * returns the number of bytes read. The symtab is in the closed
      * state afterwards, so no strings can be added to it.
      */

extern void pr_symtab(FILE *fp,int indent,char *title,t_symtab *symtab);
     /* This routine prints out a (human) readable representation of 
      * the symbol table symtab to the file fp. Ident specifies the
      * number of spaces the text should be indented. Title is used
      * to print a header text.
      */

#endif	/* _symtab_h */

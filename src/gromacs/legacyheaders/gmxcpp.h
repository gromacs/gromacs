
/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */


typedef struct gmx_cpp *gmx_cpp_t; 
	
	/* The possible return codes for these functions */
enum { eCPP_OK, eCPP_FILE_NOT_FOUND, eCPP_EOF, eCPP_SYNTAX, eCPP_INTERRUPT,
       eCPP_INVALID_HANDLE,
       eCPP_FILE_NOT_OPEN, eCPP_UNKNOWN, eCPP_NR };

/* THESE FUNCTIONS ARE NOT THREAD SAFE!! */

/* Open the file to be processed. The handle variable holds internal
   info for the cpp emulator. The cppopt variable (null terminated)
   can hold cpp options like -IXXX and -DXXX. Return integer status. 
   
   NOT THREAD SAFE 
   */
int cpp_open_file(const char *filenm,gmx_cpp_t *handlep, char **cppopts);

/* Return one whole line from the file into buf which holds at most n
   characters, for subsequent processing. Returns integer status. 
   
   NOT THREAD SAFE 
   */
int cpp_read_line(gmx_cpp_t *handlep,int n,char buf[]);

/* Return the file currently being read. 
   
   NOT THREAD SAFE 
   */
char *cpp_cur_file(const gmx_cpp_t *handlep);

/* Return the current line number. 
   
   NOT THREAD SAFE 
   */
int cpp_cur_linenr(const gmx_cpp_t *handlep);

/* Close the file! Return integer status. 
   
   NOT THREAD SAFE 
   */
int cpp_close_file(gmx_cpp_t *handlep);

/* Return a string containing the error message coresponding to status
   variable. 
   
   NOT THREAD SAFE 
   */
char *cpp_error(gmx_cpp_t *handlep,int status);

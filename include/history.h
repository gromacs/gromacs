/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 4.6.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */

#ifndef HISTORY_H_
#define HISTORY_H_


/*! \brief Initialize program for history functionality
 *  Pass in command line arguments for recording.
 *
 *  \param[in] argc     number of arguments
 *  \param[in] argv     command line arguments, including command name.
 */
int history_init(int argc, char** argv);

/*! \brief Write history information
 * Should not be called directly. Is called from thanx.  Has to be called after all files are closed.
 */
int history_write();

/*! \brief Record file usage in history
 *  Should normally not be called directly. Is called by ffopen. Has to be called after the file has
 *  already be opened.
 *
 *  \param[in] file    file pointer to the opene file
 *  \param[in] fn      file name of the file
 *  \param[in] mode    file opening mode string
 */
int history_fopen(FILE* file, const char* fn, const char* mode);
int history_fclose(FILE** file);

int history_addinput(char* str);

#endif /* HISTORY_H_ */

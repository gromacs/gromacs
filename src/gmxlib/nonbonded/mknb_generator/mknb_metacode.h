/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifndef _MKNB_METACODE_H_
#define _MKNB_METACODE_H_

#include <stdio.h>

/*! \file   mknb_metacode.h
 *  \brief Kernel generator (only for compile): C/Fortran statements
 *
 *  \internal
 *
 *  <b>This file is only used to generate the inner loop kernels
 *  at compile time, which in turn are included in Gromacs. 
 *  This code itself is NOT linked into the Gromacs library, so
 *  it does not need to be threadsafe.</b>
 *
 * mknb_metacode.h contains utility routines and variables for
 * writing code in either C or Fortran using our very primitive
 * and limited meta-language. It is a separate file since it
 * is essentially independent of the actual loops we write
 * for the nonbonded routines.
 *
 * It tries to split fortran lines correctly over multiple lines
 * if they are more than 72 characters, but check the generated
 * code to make sure...
 */





/*! \brief Kernel generator (only for compile): Current Language
 *  
 *  \internal
 *
 * <b>Only defined/used in the nonbonded kernel generator 
 * program mknb. This program is run once at compile
 * time to create the inner loops, and then discarded.
 * This source is NOT linked into any Gromacs library.</b>
 *
 *  This code is only being executed at compile time,
 *  so we allow the use of some ugly global variables.
 *
 *  If this variable has the value 1, Fortran code is being
 *  generated, otherwise C code.
 */
extern int
mknb_fortran;     




/*! \brief Kernel generator (only for compile): Current fp precision 
 *  
 *  \internal
 *
 * <b>Only defined/used in the nonbonded kernel generator 
 * program mknb. This program is run once at compile
 * time to create the inner loops, and then discarded.
 * This source is NOT linked into any Gromacs library.</b>
 *
 *  Global variable, 1 if using double, 0 for single.
 */
extern int                      
mknb_double; 




/*! \brief Kernel generator (only for compile): Generate comments or not
 *  
 *  \internal
 *
 * <b>Only defined/used in the nonbonded kernel generator 
 * program mknb. This program is run once at compile
 * time to create the inner loops, and then discarded.
 * This source is NOT linked into any Gromacs library.</b>
 *
 *  Global variable, comments are thrown away unless this is 1.
 */
extern int                 
mknb_keep_comments;    




/*! \brief Kernel generator (only for compile): Current indentation level
 * 
 * \internal
 *
 * <b>Only defined/used in the nonbonded kernel generator 
 * program mknb. This program is run once at compile
 * time to create the inner loops, and then discarded.
 * This source is NOT linked into any Gromacs library.</b>
 *
 * For each level (integer), lines will be indented 4(C) or 2(F77) spaces. 
 * This will be changed when you open/close new loops, or sometimes 
 * manually. Since this code is only run single-threaded at
 * compile time that is ok (although a bit ugly).
 */
extern int                 
mknb_indent_level;     

#define MKNB_C_INDENT_STEP         4
#define MKNB_FORTRAN_INDENT_STEP   2



/*! \brief Kernel generator (only for compile): Name of current output file
 *
 *  \internal
 *
 * <b>Only defined/used in the nonbonded kernel generator 
 * program mknb. This program is run once at compile
 * time to create the inner loops, and then discarded.
 * This source is NOT linked into any Gromacs library.</b>
 *
 *  All code generated will be written to this file.
 */
extern FILE *
mknb_output;          




/*! \brief Kernel generator (only for compile): Add comment to output
 *
 *  \internal
 *
 * <b>Only defined/used in the nonbonded kernel generator 
 * program mknb. This program is run once at compile
 * time to create the inner loops, and then discarded.
 * This source is NOT linked into any Gromacs library.</b>
 *
 *  Note that the comment will only be written to the C/Fortran file
 *  if keep_comments is 1, since this makes them compile faster in the 
 *  distribution.
 *
 *  \param text     Comment text
 */
void
mknb_comment         (char *     text);




/*! \brief Kernel generator (only for compile): Declare floating-point variable
 *
 *  \internal
 *
 * <b>Only defined/used in the nonbonded kernel generator 
 * program mknb. This program is run once at compile
 * time to create the inner loops, and then discarded.
 * This source is NOT linked into any Gromacs library.</b>
 *
 *  Depending on the current precision this will be either single or double.
 *
 *  \param name    Name of the variable
 */
void
mknb_declare_real         (char *     name);





/*! \brief Kernel generator (only for compile): Declare 4-byte fp variable
 *
 *  \internal
 *
 * <b>Only defined/used in the nonbonded kernel generator 
 * program mknb. This program is run once at compile
 * time to create the inner loops, and then discarded.
 * This source is NOT linked into any Gromacs library.</b>
 *
 *  This variable will always be exactly 4 bytes, even
 *  when the current precision is double.
 *  This is necessary to move between integer and floating-point
 *  registers to perform the single precision table lookup used for
 *  our fast inverse square root algorithm.
 * 
 *  \param name    Name of the variable
 */
void
mknb_declare_real4        (char *     name);




/*! \brief Kernel generator (only for compile): Declare fp constant
 *
 *  \internal
 *
 * <b>Only defined/used in the nonbonded kernel generator 
 * program mknb. This program is run once at compile
 * time to create the inner loops, and then discarded.
 * This source is NOT linked into any Gromacs library.</b>
 *
 *  This will be a 'const' variable in C, and a parameter in Fortran.
 *
 *  \param name    Name of the constant
 *  \param val     Constant value
 */
void
mknb_declare_const_real  (char *      name,
						  double      val);




/*! \brief Kernel generator (only for compile): Declare an integer variable
 *
 *  \internal
 *
 * <b>Only defined/used in the nonbonded kernel generator 
 * program mknb. This program is run once at compile
 * time to create the inner loops, and then discarded.
 * This source is NOT linked into any Gromacs library.</b>
 *
 *  \param name    Name of the variable
 */
void
mknb_declare_int         (char *      name);




/*! \brief Kernel generator (only for compile): Declare 4-byte integer variable
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  This variable will always be exactly 4 bytes. 
 *  This is necessary to move between integer and floating-point
 *  registers to perform the single precision table lookup used for
 *  our fast inverse square root algorithm.
 * 
 *  \param name    Name of the variable
 */
void
mknb_declare_int4         (char *      name);




/*! \brief Kernel generator (only for compile): Declare an integer constant
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  This will be a 'const' variable in C, and a parameter in Fortran.
 *
 *  \param name    Name of the constant
 *  \param val     Constant value
 */
void
mknb_declare_const_int     (char *     name,
							int        val);



/*! \brief Kernel generator (only for compile): Generic type declaration
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  Arbitrary variable declaration where you provide your own type.
 *  Note that you will have to take care of C/Fortran differences
 *  manually. In the generated code it will appear as 
 *  "type_name name;" in C, and without the semicolon in Fortran.
 *
 *  \param type_name  Name of the variable type
 *  \param name      Name of the variable
 */
void
mknb_declare_other          (char *     type_name,
							 char *     name);





/*! \brief Kernel generator (only for compile): Reference element in a vector
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  This statement will generate the code vector[idx] in C,
 *  and vector(idx) in Fortran. 
 *  Note that the code does not take into account that C start
 *  counting on 0 and Fortran on 1 - you will have to account
 *  for that manually when necessary.
 *
 *  \param vector   Name of array
 *  \param idx      Index of element to access
 *
 *  \return         Text string with the code
 *
 */   
char *
mknb_array                   (char *    vector,
			      char *    idx);




/*! \brief Kernel generator (only for compile): Formatted output to file
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 * Print a line of code to the output file.
 * This routine does proper indentation based on the global
 * indentation variable in mknb_metacode.h, and also supports the
 * same type of variable-argument lists as printf, with the
 * exception of field widths.
 *
 * This is meant as a low-level routine for raw output, so
 * we do not append any semicolon when the language is C.
 *
 * \param format   printf-like format string
 */
void
mknb_code                    (char *    format, 
							  ... );




/*! \brief Kernel generator (only for compile): Generate an assignment a=b
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 * This routine does proper indentation, and also supports the
 * same type of variable-argument lists as printf both in
 * the left and right-hand side buffers, apart from field widths.
 *
 * a statement like <tt>assign("atom%d","data%d",3,5)</tt> will produce
 *
 * <tt>atom3 = data5;</tt>     (in C)
 *
 * In contrast to code(), this routine appends a semicolon when 
 * the language is set to C.
 */
void
mknb_assign                  (char *     left,
							  char *     right, 
							  ...);




/*! \brief Kernel generator (only for compile): Start for loop block
 *
 *  \internal 
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  In C, this will generate the statement: 
 *
 *  <tt>for(var=from;var<from;var++) {</tt>
 *
 *  i.e., the loop will be executed to-from times,
 *  with the first value being lvar=from, and the
 *  last lvar=to-1.
 *
 *  In Fortran, you will get the statement:
 *
 *  <tt>FOR var=from,to DO</tt>
 *
 *  i.e., the loop will be executed to-from+1 times,
 *  with the first value being lvar=from, and the
 *  last lvar=to.
 *
 *  Take care - you will probably need to provide
 *  different limits depending on language!
 * 
 *  \param   var    Loop variable 
 *  \param   from   First value for variable
 *  \param   to     Last value for variable when using
 *                  Fortran, one more than last value in C.
 */
void
mknb_start_loop              (char *      var,
							  char *      from,
							  char *      to);




/*! \brief Kernel generator (only for compile): End an open for loop block
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 */
void
mknb_end_loop                (void);




/*! \brief Kernel generator (only for compile): Open a new if-block 
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  Starts a new if-block, which will be executed
 *  when the provided condition is true.
 *
 *  \param   cond   Condition to determine if
 *                  block should be executed.
 *                  This will be written literally to
 *                  the output, so you will need to
 *                  account for C/Fortran differences
 *                  yourself.    
 */
void
mknb_start_if                (char *      cond);




/*! \brief Kernel generator (only for compile): Close if-, open else-block
 *
 * \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 */
void
mknb_do_else                 (void);




/*! \brief Kernel generator (only for compile): End an open if/else block
 *
 *  \internal
 *
 *  <b>Only defined/used in the nonbonded kernel generator 
 *  program mknb. This program is run once at compile
 *  time to create the inner loops, and then discarded.
 *  This source is NOT linked into any Gromacs library.</b>
 *
 *  Each start_if() statement should be followed by this,
 *  possible with a do_else() somewhere in between.
 */
void
mknb_end_if                  (void);



#endif /* _MKNB_METACODE_H_ */

 

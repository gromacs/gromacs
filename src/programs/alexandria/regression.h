#ifndef ALEAXNDRIA_REGRESSION_H
#define ALEAXNDRIA_REGRESSION_H

/*! \brief Solve a matrix equation A x = b for x
 *
 * Solve a matrix equation A x = b where A is a rectangular matrix.
 * Both a and b are modified in the routine.
 * Routine may throw if unsuccessful.
 *
 * \param[in]  nrow  Number of rows in the matrix
 * \param[in]  b     Right hand side of the equation, length nrow
 * \param[in]  ncol  Number of columns in the matrix
 * \param[in]  a     The matrix
 * \param[out] x     The values to be solved, length ncol
 */
extern void multi_regression2(int      nrow, 
                              double   b[], 
                              int      ncol,
                              double **a,
                              double   x[]);

#endif

real calc_grid(matrix box,real gr_sp,int *nx,int *ny,int *nz,int nprocs);
/* Sets the number of grid points for each zero n* to the first reasonable
 * number which gives a spacing equal to or smaller than gr_sp.
 * nx and ny should be divisible by nprocs, an error is generated when this
 * can not be achieved by calc_grid.
 * Returns the maximum grid spacing.
 */

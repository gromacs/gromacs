/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_FILEIO_GMXFIO_XDR_H
#define GMX_FILEIO_GMXFIO_XDR_H

#include "gromacs/fileio/xdrf.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

struct t_fileio;

void gmx_fio_setprecision(struct t_fileio *fio, gmx_bool bDouble);
/* Select the floating point precision for reading and writing files */

XDR *gmx_fio_getxdr(struct t_fileio *fio);
/* Return the file pointer itself */

gmx_bool gmx_fio_writee_string(struct t_fileio *fio, const char *item,
                               const char *desc, const char *srcfile, int line);

/* reading or writing, depending on the file's opening mode string */
gmx_bool gmx_fio_doe_real(struct t_fileio *fio, real *item,
                          const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_float(struct t_fileio *fio, float *item,
                           const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_double(struct t_fileio *fio, double *item,
                            const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_gmx_bool(struct t_fileio *fio, gmx_bool *item,
                              const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_int(struct t_fileio *fio, int *item,
                         const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_int64(struct t_fileio *fio, gmx_int64_t *item,
                           const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_uchar(struct t_fileio *fio, unsigned char *item,
                           const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_ushort(struct t_fileio *fio, unsigned short *item,
                            const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_rvec(struct t_fileio *fio, rvec *item,
                          const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_ivec(struct t_fileio *fio, ivec *item,
                          const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_string(struct t_fileio *fio, char *item,
                            const char *desc, const char *srcfile, int line);

/* array reading & writing */
gmx_bool gmx_fio_ndoe_real(struct t_fileio *fio, real *item, int n,
                           const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_float(struct t_fileio *fio, float *item, int n,
                            const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_double(struct t_fileio *fio, double *item, int n,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_gmx_bool(struct t_fileio *fio, gmx_bool *item, int n,
                               const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_int(struct t_fileio *fio, int *item, int n,
                          const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_int64(struct t_fileio *fio, gmx_int64_t *item, int n,
                            const char *desc, const char *srcfile,
                            int line);
gmx_bool gmx_fio_ndoe_uchar(struct t_fileio *fio, unsigned char *item, int n,
                            const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_ushort(struct t_fileio *fio, unsigned short *item, int n,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_rvec(struct t_fileio *fio, rvec *item, int n,
                           const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_ivec(struct t_fileio *fio, ivec *item, int n,
                           const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_string(struct t_fileio *fio, char *item[], int n,
                             const char *desc, const char *srcfile, int line);



/* convenience macros */
#define gmx_fio_write_string(fio, item)         gmx_fio_writee_string(fio, item, (#item), __FILE__, __LINE__)

#define gmx_fio_do_real(fio, item)              gmx_fio_doe_real(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_float(fio, item)             gmx_fio_doe_float(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_double(fio, item)            gmx_fio_doe_double(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_gmx_bool(fio, item)          gmx_fio_doe_gmx_bool(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_int(fio, item)               gmx_fio_doe_int(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_int64(fio, item)             gmx_fio_doe_int64(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_uchar(fio, item)             gmx_fio_doe_uchar(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_ushort(fio, item)            gmx_fio_doe_ushort(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_rvec(fio, item)              gmx_fio_doe_rvec(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_ivec(fio, item)              gmx_fio_doe_ivec(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_string(fio, item)            gmx_fio_doe_string(fio, item, (#item), __FILE__, __LINE__)


#define gmx_fio_ndo_real(fio, item, n)              gmx_fio_ndoe_real(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_float(fio, item, n)             gmx_fio_ndoe_float(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_double(fio, item, n)            gmx_fio_ndoe_double(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_gmx_bool(fio, item, n)          gmx_fio_ndoe_gmx_bool(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_int(fio, item, n)               gmx_fio_ndoe_int(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_int64(fio, item, n)             gmx_fio_ndoe_int64(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_uchar(fio, item, n)             gmx_fio_ndoe_uchar(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_ushort(fio, item, n)            gmx_fio_ndoe_ushort(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_rvec(fio, item, n)              gmx_fio_ndoe_rvec(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_ivec(fio, item, n)              gmx_fio_ndoe_ivec(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_string(fio, item, n)            gmx_fio_ndoe_string(fio, item, n, (#item), __FILE__, __LINE__)

#ifdef __cplusplus
}
#endif

#endif

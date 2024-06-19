/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_FILEIO_GMXFIO_XDR_H
#define GMX_FILEIO_GMXFIO_XDR_H

#include <cstddef>
#include <cstdint>

#include <string>

#include "gromacs/fileio/xdrf.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/real.h"

struct t_fileio;

void gmx_fio_setprecision(struct t_fileio* fio, gmx_bool bDouble);
/* Select the floating point precision for reading and writing files */

bool gmx_fio_is_double(struct t_fileio* fio);

XDR* gmx_fio_getxdr(struct t_fileio* fio);
/* Return the file pointer itself */

gmx_bool gmx_fio_writee_string(struct t_fileio* fio, const char* item, const char* desc, const char* srcfile, int line);

/* reading or writing, depending on the file's opening mode string */
gmx_bool gmx_fio_doe_real(struct t_fileio* fio, real* item, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_doe_float(struct t_fileio* fio, float* item, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_doe_double(struct t_fileio* fio, double* item, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_doe_gmx_bool(struct t_fileio* fio, gmx_bool* item, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_doe_int(struct t_fileio* fio, int* item, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_doe_int32(struct t_fileio* fio, int32_t* item, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_doe_int64(struct t_fileio* fio, int64_t* item, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_doe_uchar(struct t_fileio* fio, unsigned char* item, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_doe_char(struct t_fileio* fio, char* item, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_doe_ushort(struct t_fileio* fio,
                            unsigned short*  item,
                            const char*      desc,
                            const char*      srcfile,
                            int              line);
gmx_bool gmx_fio_doe_rvec(struct t_fileio* fio, rvec* item, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_doe_ivec(struct t_fileio* fio, ivec* item, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_doe_string(struct t_fileio* fio, char* item, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_doe_opaque(struct t_fileio* fio,
                            char*            item,
                            std::size_t      size,
                            const char*      desc,
                            const char*      srcfile,
                            int              line);

/* array reading & writing */
gmx_bool gmx_fio_ndoe_real(struct t_fileio* fio, real* item, int n, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_ndoe_float(struct t_fileio* fio, float* item, int n, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_ndoe_double(struct t_fileio* fio, double* item, int n, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_ndoe_gmx_bool(struct t_fileio* fio,
                               gmx_bool*        item,
                               int              n,
                               const char*      desc,
                               const char*      srcfile,
                               int              line);
gmx_bool gmx_fio_ndoe_int(struct t_fileio* fio, int* item, int n, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_ndoe_int32(struct t_fileio* fio, int32_t* item, int n, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_ndoe_int64(struct t_fileio* fio, int64_t* item, int n, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_ndoe_uchar(struct t_fileio* fio,
                            unsigned char*   item,
                            int              n,
                            const char*      desc,
                            const char*      srcfile,
                            int              line);
gmx_bool gmx_fio_ndoe_char(struct t_fileio* fio, char* item, int n, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_ndoe_ushort(struct t_fileio* fio,
                             unsigned short*  item,
                             int              n,
                             const char*      desc,
                             const char*      srcfile,
                             int              line);
gmx_bool gmx_fio_ndoe_rvec(struct t_fileio* fio, rvec* item, int n, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_ndoe_ivec(struct t_fileio* fio, ivec* item, int n, const char* desc, const char* srcfile, int line);
gmx_bool gmx_fio_ndoe_string(struct t_fileio* fio, char* item[], int n, const char* desc, const char* srcfile, int line);


/* convenience macros */
#define gmx_fio_write_string(fio, item) \
    gmx_fio_writee_string(fio, item, (#item), __FILE__, __LINE__)

#define gmx_fio_do_real(fio, item) gmx_fio_doe_real(fio, &(item), (#item), __FILE__, __LINE__)
#define gmx_fio_do_float(fio, item) gmx_fio_doe_float(fio, &(item), (#item), __FILE__, __LINE__)
#define gmx_fio_do_double(fio, item) gmx_fio_doe_double(fio, &(item), (#item), __FILE__, __LINE__)
#define gmx_fio_do_gmx_bool(fio, item) \
    gmx_fio_doe_gmx_bool(fio, &(item), (#item), __FILE__, __LINE__)
#define gmx_fio_do_int(fio, item) gmx_fio_doe_int(fio, &(item), (#item), __FILE__, __LINE__)
#define gmx_fio_do_int32(fio, item) gmx_fio_doe_int32(fio, &(item), (#item), __FILE__, __LINE__)
#define gmx_fio_do_int64(fio, item) gmx_fio_doe_int64(fio, &(item), (#item), __FILE__, __LINE__)
#define gmx_fio_do_uchar(fio, item) gmx_fio_doe_uchar(fio, &(item), (#item), __FILE__, __LINE__)
#define gmx_fio_do_char(fio, item) gmx_fio_doe_char(fio, &(item), (#item), __FILE__, __LINE__)
#define gmx_fio_do_ushort(fio, item) gmx_fio_doe_ushort(fio, &(item), (#item), __FILE__, __LINE__)
#define gmx_fio_do_rvec(fio, item) gmx_fio_doe_rvec(fio, &(item), (#item), __FILE__, __LINE__)
#define gmx_fio_do_ivec(fio, item) gmx_fio_doe_ivec(fio, &(item), (#item), __FILE__, __LINE__)
#define gmx_fio_do_string(fio, item) gmx_fio_doe_string(fio, item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_opaque(fio, item, size) \
    gmx_fio_doe_opaque(fio, item, size, (#item), __FILE__, __LINE__)


#define gmx_fio_ndo_real(fio, item, n) gmx_fio_ndoe_real(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_float(fio, item, n) \
    gmx_fio_ndoe_float(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_double(fio, item, n) \
    gmx_fio_ndoe_double(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_gmx_bool(fio, item, n) \
    gmx_fio_ndoe_gmx_bool(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_int(fio, item, n) gmx_fio_ndoe_int(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_int32(fio, item, n) \
    gmx_fio_ndoe_int32(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_int64(fio, item, n) \
    gmx_fio_ndoe_int64(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_uchar(fio, item, n) \
    gmx_fio_ndoe_uchar(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_char(fio, item, n) gmx_fio_ndoe_char(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_ushort(fio, item, n) \
    gmx_fio_ndoe_ushort(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_rvec(fio, item, n) gmx_fio_ndoe_rvec(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_ivec(fio, item, n) gmx_fio_ndoe_ivec(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_string(fio, item, n) \
    gmx_fio_ndoe_string(fio, item, n, (#item), __FILE__, __LINE__)

namespace gmx
{
/*!\internal \brief
 * Serializer to read/write XDR data.
 */
class FileIOXdrSerializer : public ISerializer
{
public:
    //! Only create with valid file I/O handle.
    explicit FileIOXdrSerializer(t_fileio* fio);

    //! If file is open in reading mode.
    bool reading() const override;
    //! Handle bool I/O.
    void doBool(bool* value) override;
    //! Handle unsigned char I/O.
    void doUChar(unsigned char* value) override;
    //! Handle char I/O.
    void doChar(char* value) override;
    //! Handle unsigned short I/O.
    void doUShort(unsigned short* value) override;
    //! Handle default integer I/O.
    void doInt(int* value) override;
    //! Handle int32 I/O.
    void doInt32(int32_t* value) override;
    //! Handle int64 I/O.
    void doInt64(int64_t* value) override;
    //! Handle single precision float I/O.
    void doFloat(float* value) override;
    //! Handle double precision float I/O.
    void doDouble(double* value) override;
    //! Handle GROMACS floating point number I/O.
    void doReal(real* value) override;
    //! Handle I/O of integer vector of size DIM.
    void doIvec(ivec* value) override;
    //! Handle I/O of GROMACS real vector of size DIM.
    void doRvec(rvec* value) override;
    //! Handle I/O if string.
    void doString(std::string* value) override;
    //! Handle opaque data.
    void doOpaque(char* data, std::size_t size) override;
    //! Special case for handling I/O of a vector of characters.
    void doCharArray(char* values, int elements) override;
    //! Special case for handling I/O of a vector of unsigned characters.
    void doUCharArray(unsigned char* values, int elements) override;
    //! Special case for handling I/O of a vector of rvecs.
    void doRvecArray(rvec* values, int elements) override;

private:
    //! File I/O handle.
    t_fileio* fio_;
};

} // namespace gmx

#endif

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements methods from voxels.h
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_fileio
 * \ingroup module_griddata
 */
#include "gmxpre.h"

#include "griddataio.h"
#include "mrcmetadata.h"

#include <cstdio>

#include <algorithm>
#include <complex>
#include <string>
#include <vector>
#include <set>
#include <type_traits>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/math/units.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/griddata/operations/realfieldmeasure.h"

namespace gmx
{

/*******************************************************************************
 * MrcFile::Impl
 */
class MrcFile::Impl
{
    public:
        Impl();
        ~Impl();

        /*! \brief Opens file stream in file_.
         * \param[in] filename typically *.mrc, *.map, or *.ccp4
         * \param[in] bRead Open the file for reading if true, otherwise for writing
         */
        void open_file(std::string filename, bool bRead);
        void close_file();
        std::string filetype(std::string filename);
        bool known_extension(std::string filename);
        void read_file_size();

        void do_mrc_data_(GridDataReal3D &grid_data, bool bRead);
        void do_mrc_header_(bool bRead);

        // std::array<int, 3> to_crs_order(const std::array<int, 3> &xyz_order);

        /*! \brief Guess, whether endianess differs between input file and reading architecture .
         *
         * If the number of columns in the density file is negative or larger than 65534,
         * assume endianess missmatch between input file and reading machine architecture.*/
        void check_swap_bytes();

        bool has_skew_matrix();
        void set_skew_matrix(matrix skew);


        template <typename T> void do_(T * result, bool bRead)
        {
            if (bRead)
            {
                if (fread(result, sizeof(T), 1, file_) != 1)
                {
                    GMX_THROW(gmx::FileIOError("Cannot read from volume data at " + std::to_string(ftell(file_)) + "."));
                }
            }
            // swap bytes for correct endianness
            if (meta_.swap_bytes)
            {
                // byte swap for complex numbers, swap real and imaginary part seperately
                if (std::is_same<T, t_complex>())
                {
                    GMX_THROW(NotImplementedError("Cannot read complex numbers for now."));
                }
                // byte swap for real numbers
                if (std::is_same<T, double>())
                {
                    gmx_int64_t int_tmp = gmx_int64_t(*result);
                    *result = (int_tmp & 0xFF00000000000000) >> 7*8 | (int_tmp & 0x00FF000000000000) >> 5*8 |
                        (int_tmp & 0x0000FF0000000000) >> 3*8 | (int_tmp & 0x000000FF00000000) >> 1*8 |
                        (int_tmp & 0x00000000FF000000) << 1*8 | (int_tmp & 0x0000000000FF0000) << 3*8 |
                        (int_tmp & 0x000000000000FF00) << 5*8 | (int_tmp & 0x00000000000000FF) << 7*8;
                }

                // byte swap for real numbers
                if (std::is_same<T, float>())
                {
                    gmx_int32_t int_tmp = gmx_int32_t(*result);
                    *result = (int_tmp & 0xFF000000) >> 24 | (int_tmp & 0x00FF0000) >> 8 | (int_tmp & 0x0000FF00) << 8 | (int_tmp & 0x000000FF) << 24;
                }
                if (std::is_same<T, int16_t>())
                {
                    int16_t int_tmp = int16_t(*result);
                    *result = (int_tmp & 0xFF00) >> 8 | (int_tmp & 0x00FF) << 8;
                }
            }

            if (!bRead)
            {
                fwrite(result, sizeof(T), 1, file_);
            }
        }

        void do_float32_rvec_(RVec * result, bool bRead);
        void do_int32_ivec_(std::array<int, 3> * result, bool bRead);

        bool colummn_row_section_order_valid_(std::array<int, 3> crs_to_xyz);

        FILE     *file_;
        long      file_size_;
        const int number_labels;
        const int label_size;
        const int header_bytes;
        const std::vector<std::string> filetypes;

        MrcMetaData                    meta_;
};


void MrcFile::Impl::do_int32_ivec_(std::array<int, 3> * result, bool bRead)
{
    do_<gmx_int32_t>(&(*result)[XX], bRead);
    do_<gmx_int32_t>(&(*result)[YY], bRead);
    do_<gmx_int32_t>(&(*result)[ZZ], bRead);
}

void MrcFile::Impl::do_float32_rvec_(RVec * result, bool bRead)
{
    do_<float>(&((*result)[XX]), bRead);
    do_<float>(&((*result)[YY]), bRead);
    do_<float>(&((*result)[ZZ]), bRead);
}


std::string MrcFile::Impl::filetype(std::string filename)
{
    std::string result = "";
    size_t      pos    = filename.rfind(".");
    if (pos == 0 || pos == std::string::npos)
    {
        return result;
    }
    else
    {
        return filename.substr(pos+1);
    }
}

bool MrcFile::Impl::colummn_row_section_order_valid_(std::array<int, 3> crs_to_xyz)
{
    const std::set<int> valid_crs_set {
        0, 1, 2
    };
    std::set<int> crs_set {
        crs_to_xyz[XX], crs_to_xyz[YY], crs_to_xyz[ZZ]
    };
    return valid_crs_set == crs_set;
};

void MrcFile::Impl::do_mrc_header_(bool bRead)
{
    if (bRead)
    {
        check_swap_bytes();
        read_file_size();
        if (file_size_ < header_bytes)
        {
            GMX_THROW(gmx::FileIOError("Density format file is smaller than expected header."));
        }
    }

    /* Supports reading according to
       ftp://ftp.wwpdb.org/pub/emdb/doc/Map-format/current/EMDB_map_format.pdf
       note, that
       http://www.ccpem.ac.uk/mrc_format/mrc2014.php
       differs slightly in definition */

    /* 1-3 | NC, NR, NS | signed int >0
     * # of columns (fastest changing),rows, sections (slowest changing)
     * emdb convention: NC=NR=NS                     */

    do_int32_ivec_(&meta_.num_crs, bRead);

    /* 4   | MODE | signed int | 0,1,2,3,4
     * voxel datatype
     * emdb convention: 2       */
    do_<gmx_int32_t>(&meta_.mrc_data_mode, bRead);

    /* MODE = 0: 8 bits, density stored as a signed byte (range -128 to 127, ISO/IEC 10967)
     * MODE = 1: 16 bits, density stored as a signed integer (range -32768 to 32767, ISO/IEC 10967)
     * MODE = 2: 32 bits, density stored as a floating point number (IEEE 754)
     * MODE = 3: 32 bits, Fourier transform stored as complex signed integers (ISO/IEC 10967)
     * MODE = 4: 64 bits, Fourier transform stored as complex floating point numbers (IEEE 754)     */
    if (meta_.mrc_data_mode < 0 || meta_.mrc_data_mode > 4)
    {
        GMX_THROW(gmx::FileIOError("Read invalid mrc/cpp4/imod data mode. Mode " + std::to_string(meta_.mrc_data_mode) + " not in [0..4] ."));
    }
    if (meta_.mrc_data_mode != 2)
    {
        GMX_THROW(gmx::NotImplementedError("Other mrc/ccp4/imod data modes than 32 bit float not currently implemented. Upgrading to the newest gromacs verison might possibly help."));
    }

    /* 5-7 | NCSTART, NRSTART, NSSTART | signed int
     * position of first column, first row, and first section (voxel grid units)
     *
     * The position of the first voxel is defined in grid units by NCSTART, NRSTART, and NSSTART.
     * The center of the voxel with grid position (0,0,0) corresponds to the Cartesian coordinate origin.*/
    do_int32_ivec_(&meta_.crs_start, bRead);

    /* 8-10 | NX, NY, NZ | signed int >0 |
     * intervals per unit cell repeat along X,Y Z
     * intervals per map length along X,Y,Z;
     * emdb convention: same as NC, NR, NS
     *
     * Lengths in Aangstroms for a single voxel are as follows:
     * Xvoxel = X_LENGTH/NX Yvoxel = Y_LENGTH/NY Zvoxel = Z_LENGTH/NZ */

    do_int32_ivec_(&meta_.extend, bRead);

    /* 11-13 | X_LENGTH, Y_LENGTH, Z_LENGTH | floating pt >0
     * Unit Cell repeats along X, Y, Z In Aangstrom
     * emdb Map lengths along X,Y,Z in Aangstrom
     *
     * Lengths in Aangstroms for a single voxel are as follows:
     * Xvoxel = X_LENGTH/NX Yvoxel = Y_LENGTH/NY Zvoxel = Z_LENGTH/NZ */
    do_float32_rvec_(&meta_.cell_length, bRead);

    /* 14-16 | ALPHA,BETA,GAMMA | floating pt >0, <180
     * Unit Cell angles (degrees)
     * emdb convention: 90, 90, 90
     *
     * By convention, cell angles (ALPHA, BETA, GAMMA)
     * are 90 degrees for single particle or tomogram EM maps;
     * they follow IUCr space group conventions for crystals.*/
    do_float32_rvec_(&meta_.cell_angles, bRead);
    // By convention, unset cell angles (all 0) are interpreted as 90 deg.
    if (meta_.cell_angles[XX]*meta_.cell_angles[YY]*meta_.cell_angles[ZZ] < 1e-5)
    {
        meta_.cell_angles = {90, 90, 90};
    }

    /* 17-19 | MAPC, MAPR, MAPS | signed int | 1 (=X) 2 (=Y) 3 (=Z)
     * relationship of X,Y,Z axes to columns, rows, sections
     * emdb convention: 1, 2, 3 */
    std::array<int, 3> crs_to_xyz {{
                                       meta_.crs_to_xyz[XX]+1, meta_.crs_to_xyz[YY]+1, meta_.crs_to_xyz[ZZ]+1
                                   }};
    do_int32_ivec_(&crs_to_xyz, bRead);

    if (bRead)
    {
        meta_.crs_to_xyz[XX] = crs_to_xyz[XX]-1;
        meta_.crs_to_xyz[YY] = crs_to_xyz[YY]-1;
        meta_.crs_to_xyz[ZZ] = crs_to_xyz[ZZ]-1;
        if (!colummn_row_section_order_valid_(meta_.crs_to_xyz))
        {
            meta_.crs_to_xyz = {{0, 1, 2}};
        }
        meta_.xyz_to_crs = meta_.to_xyz_order(meta_.crs_to_xyz);
    }


    /* 20-22 | AMIN, AMAX, AMEAN | floating pt
     * Minimum, maximum, average density */
    do_<float>(&(meta_.min_value ), bRead);
    do_<float>(&(meta_.max_value ), bRead);
    do_<float>(&(meta_.mean_value), bRead);

    /* 23 | ISPG | signed int 1-230 |
     * space group #
     * emdb convention 1
     *
     * Space Group Numbers are defined by IUCr conventions
     * (Table 12.3.4.1 Standard space-group symbols”, pages 824-831,
     * International Tables for Crystallography, Volume A, fifth edition).
     *
     * For 3D volumes of single particle or tomogram entries, ISPG=1 and NSYMBT=0.
     * For image stacks ISPG = 0 */
    do_<gmx_int32_t>(&meta_.space_group, bRead);

    /* 24 | NSYMBT | signed int | 80n
     * # of bytes in symmetry table (multiple of 80)
     * emdb convention 0 */
    do_<gmx_int32_t>(&meta_.num_bytes_extened_header, bRead);
    if (meta_.num_bytes_extened_header%80 != 0)
    {
        GMX_THROW(gmx::FileIOError("Read invalid number of bytes in symbol table from mrc/cpp4/imod file. Should be 80, but is " + std::to_string(meta_.num_bytes_extened_header) + "instead."));
    }

    if (meta_.is_crystallographic)
    {
        /* 25 | LSKFLG | signed int | 0,1
         * flag for skew matrix
         * emdb convention 0 */
        gmx_int32_t hasSkewMatrix = meta_.has_skew_matrix ? 1 : 0;
        do_<gmx_int32_t>(&hasSkewMatrix, bRead);
        if (bRead)
        {
            if (!(hasSkewMatrix == 0 || hasSkewMatrix == 1))
            {
                GMX_THROW(gmx::FileIOError("Skew matrix flag set to invalid value in mrc/cpp4/imod file. Should be 0 or 1 but is " + std::to_string(hasSkewMatrix) + "instead."));
            }
            meta_.has_skew_matrix = (hasSkewMatrix == 1) ? true : false;
        }

        if (meta_.has_skew_matrix)
        {
            /* TODO: A2NM conversion for skew matrix if necessary */
            /* 26-34 | SKWMAT | floating pt
             * skew matrix-S11, S12, S13, S21, S22, S23, S31, S32, S33
             * emdb convention: not set
             *
             * SKWMAT, SKWTRN, and EXTRA fields are not currently used by EMDB. */

            /* 35-37 | SKWTRN | floating pt
             * skew translation-T1, T2, T3
             * emdb convention: not set
             *
             * SKWMAT, SKWTRN, and EXTRA fields are not currently used by EMDB. */
            for (auto && i : meta_.skew_matrix)
            {
                do_<float>(&i, bRead);
            }
            do_float32_rvec_(&meta_.skew_translation, bRead);
        }
    }
    else
    {
        /* 25-37 not used in EMDB */
        for (auto && i : meta_.extraskew)
        {
            do_<float>(&i, bRead);
        }
    }

    /* 38-52 | EXTRA | 32 bit binary
     * user-defined metadata
     *
     * SKWMAT, SKWTRN, and EXTRA fields are not currently used by EMDB.
     * EMDB might use fields 50,51 and 52 for setting the coordinate system origin */
    for (auto && i : meta_.extra)
    {
        do_<float>(&i, bRead);
    }

    /* 53 | MAP | ASCII char
     * "MAP "
     * MRC/CCP4 MAP format identifier */
    meta_.format_identifier.resize(4);
    for (size_t i = 0; i < 4; i++)
    {
        do_<char>(&meta_.format_identifier[i], bRead);
    }
    if (!(meta_.format_identifier.compare("MAP ") == 0))
    {
        fprintf(stderr, "\nWARNING: Expected " "MAP " " as format identifier.\n");
    }
    /* 54 | MACHST | 32 bit
     * binary machine stamp
     *
     * MACHST is (written/read as 4 hex byte sequence)
     * 0x44,0x41,0x00,0x00  for little endian machines
     * 0x11,0x11,0x00,0x00  for big endian machines
     */
    do_<gmx_int32_t>(&(meta_.machine_stamp), bRead);

    /* 55 | RMS | floating pt
     * Density root-mean-square deviation */
    do_<float>(&(meta_.rms_value), bRead);

    /* 56 | NLABL | signed int | 0-10
     * # of labels
     *
     * Following the 2010 remediation, maps distributed by EMDB
     * now have a single label of form “::::EMDataBank.org::::EMD-1234::::”.  */
    do_<gmx_int32_t>(&meta_.num_labels, bRead);

    /* 57-256 | LABEL_N | ASCII char
     * Up to 10 user-defined labels */
    if (bRead)
    {
        std::string mrc_label(label_size, ' ');
        meta_.labels.clear();
        for (int i = 0; i < number_labels; i++)
        {
            int read_size = fread(&mrc_label[0], 1, label_size, file_);
            if (read_size != label_size)
            {
                GMX_THROW(gmx::FileIOError("Could not read label from file."));
            }
            meta_.labels.push_back(std::string(mrc_label));
        }
    }
    else
    {
        for (auto label : meta_.labels)
        {
            fprintf(file_, "%.*s", label_size, label.c_str());
        }
    }

    /* 257-257+NSYMBT | anything
     */
    if (bRead)
    {
        if (file_size_ < meta_.num_bytes_extened_header + header_bytes)
        {
            GMX_THROW(gmx::FileIOError("Density format file is smaller than expected extended header."));
        }
        meta_.extended_header.resize(meta_.num_bytes_extened_header);
        if (fgets(&((meta_.extended_header)[0]), meta_.num_bytes_extened_header, file_) != nullptr)
        {
            GMX_THROW(gmx::FileIOError("Cannot read extended header from file."));
        }
    }
    else
    {
        fwrite(&((meta_.extended_header)[0]), sizeof(char), meta_.extended_header.size(), file_);
    }

};

void MrcFile::Impl::do_mrc_data_(GridDataReal3D &grid_data, bool bRead)
{
    const auto &lattice = grid_data.getGrid().lattice();
    if (bRead)
    {
        if (file_size_ < header_bytes + meta_.num_bytes_extened_header + long(lattice.getNumLatticePoints()*sizeof(float)) )
        {
            GMX_THROW(gmx::FileIOError("Density format file is smaller than indicated in its header."));
        }
        if (file_size_ > header_bytes + meta_.num_bytes_extened_header +  long(lattice.getNumLatticePoints()*sizeof(float)) )
        {
            fprintf(stderr, "WARNING : Density format file size is %ld, however header (%d) + symbol table (%d) + data  (%ld) is larger than indicated in its header. Reading anyway.. \n", file_size_, header_bytes, meta_.num_bytes_extened_header, lattice.getNumLatticePoints()*sizeof(float));
        }
    }

    auto num_crs = meta_.to_crs_order(lattice.extend());
    for (int section = 0; section  < num_crs[ZZ]; section++)
    {
        for (int row  = 0; row  < num_crs[YY]; row++)
        {
            for (int column  = 0; column  < num_crs[XX]; column++)
            {
                const auto &currentPosition =  grid_data.iteratorAtMultiIndex(meta_.to_xyz_order({{column, row, section}}));
                do_<float>(&(*currentPosition), bRead);
            }
        }
    }
}

void MrcFile::Impl::read_file_size()
{
    fpos_t current;
    fgetpos(file_, &current);
    fseek(file_, 0, SEEK_END);
    file_size_ = ftell(file_);
    fsetpos(file_, &current);
}

void MrcFile::Impl::check_swap_bytes()
{

    fpos_t      current;
    meta_.swap_bytes = false;
    gmx_int32_t number_columns;
    fgetpos(file_, &current);

    fseek(file_, 0, SEEK_SET);
    do_<gmx_int32_t>(&number_columns, true);
    if (number_columns <= 0 || number_columns >= 65536)
    {
        meta_.swap_bytes = true;
    }

    // rewind the file
    fsetpos(file_, &current);

}

bool MrcFile::Impl::known_extension(std::string filename)
{
    return std::find(filetypes.begin(), filetypes.end(), filetype(filename)) != filetypes.end();
}

void MrcFile::Impl::open_file(std::string filename, bool bRead)
{

    if (filename.empty())
    {
        GMX_THROW(gmx::FileIOError("Filename empty."));
    }

    filename.erase(std::remove(filename.begin(), filename.end(), '\n'), filename.end());

    if (bRead)
    {
        if (!known_extension(filename) )
        {
            GMX_THROW(gmx::FileIOError("Cannot read filetype " + filetype(filename) + "."));
        }
        file_ = gmx_fio_fopen(filename.c_str(), "r");
    }
    else
    {
        file_ = gmx_fio_fopen(filename.c_str(), "w");
    }
    if (file_ == nullptr)
    {
        GMX_THROW(gmx::FileIOError("Cannot open " + filename + ". "));
    }

}

void MrcFile::Impl::close_file()
{
    gmx_fio_fclose(file_);
}

MrcFile::Impl::Impl() : file_(nullptr), file_size_(0),
                        number_labels(10), label_size(80), header_bytes(1024),
                        filetypes({"mrc", "ccp4", "imod", "map"}
                                  )
{
    meta_.setEMDBDefaults();
};

MrcFile::Impl::~Impl()
{
    if (file_ != nullptr)
    {
        close_file();
    }
};



/*******************************************************************************
 * MrcFile
 */
MrcFile::MrcFile() : impl_(new MrcFile::Impl)
{
}

MrcFile::~MrcFile()
{
    impl_->close_file();
}

void MrcFile::write_with_own_meta(std::string filename, GridDataReal3D &grid_data, const MrcMetaData &meta, bool bOwnGridStats)
{
    bool bRead = false;

    impl_->meta_ = meta;

    if (bOwnGridStats)
    {
        impl_->meta_.set_grid_stats(grid_data);
    }

    impl_->open_file(filename, bRead);
    impl_->do_mrc_header_(bRead);
    impl_->do_mrc_data_(*(const_cast<GridDataReal3D*>(&grid_data)), bRead);
    impl_->close_file();
}

void MrcFile::write(std::string filename, const GridDataReal3D &grid_data)
{
    bool bRead = false;
    impl_->meta_.set_grid_stats(grid_data);
    impl_->open_file(filename, bRead);
    impl_->meta_.fromGrid(grid_data.getGrid());
    impl_->do_mrc_header_(bRead);
    impl_->do_mrc_data_(*(const_cast<GridDataReal3D*>(&grid_data)), bRead);
    impl_->close_file();
}


void MrcFile::read_meta(std::string filename, MrcMetaData *meta)
{
    bool        bRead = true;
    impl_->open_file(filename, bRead);
    impl_->do_mrc_header_(bRead);
    impl_->close_file();
    *meta = impl_->meta_;
}

GridDataReal3D MrcFile::read_with_meta(std::string filename, MrcMetaData *meta)
{
    auto result = read(filename);
    *meta = impl_->meta_;
    return result;
}

GridDataReal3D MrcFile::read(std::string filename)
{
    bool bRead = true;
    impl_->open_file(filename, bRead);
    impl_->do_mrc_header_(bRead);
    auto result = GridDataReal3D(impl_->meta_.toGrid());
    impl_->do_mrc_data_(result, bRead);
    impl_->close_file();
    return result;
}

} // namespace gmx

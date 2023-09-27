/*
 * Copyright © 2014-2015 Klaus Reuter
 * Copyright © 2013-2014 Felix Höfling
 * Copyright © 2013      Manuel Dibak
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */


#ifndef H5XX_FILE_HPP
#define H5XX_FILE_HPP

#include <h5xx/hdf5_compat.hpp>
#include <h5xx/error.hpp>
#include <h5xx/utility.hpp>

#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>

#ifdef H5XX_USE_MPI
#include <mpi.h>
#endif



namespace h5xx {

inline htri_t is_hdf5_file(std::string const& filename)
{
    H5E_BEGIN_TRY {
        return H5Fis_hdf5(filename.c_str());
    } H5E_END_TRY
}

/**
 *  Represent HDF5 file. Instances of the class cannot be copied.
 *
 *  The following flags are defined to specify the opening mode:
 *
 *      file::in       read access
 *      file::out      write access, append to existing file
 *      file::trunc    write access, truncate existing file
 *      file::excl     write access, file must not exist
 *
 *  The flags may be combined by bitwise OR. Read access is always granted by
 *  the HDF5 library, so file::in may be omitted. The flags file::trunc and
 *  file::excl are mutually exclusive and imply file::out.
 */
class file
{
public:
    enum mode
    {
        in = 0x0000u        /* H5F_ACC_RDONLY absence of rdwr => rd-only */
      , out = 0x0001u       /* H5F_ACC_RDWR open for read and write    */
      , trunc = 0x0002u     /* H5F_ACC_TRUNC overwrite existing files   */
      , excl = 0x0004u      /* H5F_ACC_EXCL fail if file already exists*/
    };

    /** default constructor */
    file() : hid_(-1), plid_(H5P_DEFAULT) {}

    /** open file upon construction */
    explicit file(std::string const& filename, unsigned mode = in | out);

#ifdef H5XX_USE_MPI
    // --- default arguments require mode to be the last argument
    explicit file(std::string const& filename, MPI_Comm comm, MPI_Info info, unsigned mode = in | out);
#endif

    /**
     * deleted copy constructor
     *
     * Calling the constructor throws an exception. Copying must be elided by
     * return value optimisation. See also "file h5xx::move(file&)".
     */
    file(file const& other);

    /**
     * assignment operator
     *
     * Uses the copy-and-swap idiom. Move semantics is achieved in conjunction
     * with "file h5xx::move(file&)", i.e., the group object on the right
     * hand side is empty after move assignment.
     */
    file & operator=(file other);

    /** close file upon destruction */
    ~file();

    /** return HDF5 object ID */
    hid_t hid() const
    {
        return hid_;
    }

    /** returns true if associated to a valid HDF5 object */
    bool valid() const
    {
        return hid_ >= 0;
    }

    /** open HDF5 file in specified mode */
    void open(std::string const& filename, unsigned mode = in | out);

    /** close HDF5 file
     *
     *  If strict is true, throw h5xx::error if open HDF5 objects are
     *  associated with the file
     */
    void close(bool strict=false);

    /** flush all buffers associated with the file to disk */
    void flush() const;

    /** return filename on disk */
    std::string name() const;

private:
    /** HDF5 object ID */
    hid_t hid_;

    /** file access property list ID, required for MPI parallel functionality */
    hid_t plid_;

    template <typename h5xxObject>
    friend void swap(h5xxObject& left, h5xxObject& right);
};

inline file::file(std::string const& filename, unsigned mode)
  : hid_(-1),plid_(H5P_DEFAULT)
{
    open(filename, mode);
}

#ifdef H5XX_USE_MPI
inline file::file(std::string const& filename, MPI_Comm comm, MPI_Info info, unsigned mode)
  : hid_(-1),plid_(H5P_DEFAULT)
{
    plid_ = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plid_, comm, info);
    open(filename, mode);
}
#endif

inline file::file(file const& other)
{
    // copying would be safe if the exception were disabled.
    throw error("h5xx::file can not be copied. Copying must be elided by return value optimisation.");
    hid_ = H5Freopen(other.hid_);
    plid_ = H5Fget_access_plist(hid_);
}

inline file & file::operator=(file other)
{
    swap(*this, other);
    return *this;
}

inline file::~file()
{
    close();
}

inline void file::open(std::string const& filename, unsigned mode)
{
    // check that object is not yet in use
    if (hid_ >= 0) {
        throw error("h5xx::file object is already open");
    }

    // check for conflicting combination of opening flags
    if ((mode & trunc) && (mode & excl)) {
        throw error("h5xx::file: conflicting opening mode: " + boost::lexical_cast<std::string>(mode));
    }

    htri_t is_hdf5 = is_hdf5_file(filename);
    if (is_hdf5 >= 0 && !(mode & trunc)) {
        // file exists and may be valid HDF5, but shall not be truncated
        if (mode & excl) {
            throw error("refuse to overwrite existing HDF5 file: " + filename);
        }
        else { // open file, either to append or read-only
            if (is_hdf5 == 0) {
                throw error("not a valid HDF5 file: " + filename);
            }
            // use that "in" and "out" are equal to H5F_ACC_RDONLY and H5F_ACC_RDWR, resp.
            hid_ = H5Fopen(filename.c_str(), mode & (in | out), plid_);
        }
    }
    else {
        // file does not exist (or other error), or it exists, but shall be truncated
        if (mode == in) { // read-only
            throw error("read-only access to non-existing HDF5 file: " + filename);
        }
        // create new file
        hid_ = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plid_);
    }
    if (hid_ < 0) {
        throw error("opening or creation of HDF5 file \"" + filename + "\" failed");
    }
}

inline void file::flush() const
{
    if (hid_ < 0) {
        return;
    }
    if (H5Fflush(hid_, H5F_SCOPE_LOCAL) < 0) {
        throw error("flushing HDF5 file: " + name());
    }
}

inline void file::close(bool strict)
{
    if (hid_ < 0) {
        return;
    }
    if (strict) {
        ssize_t count = H5Fget_obj_count(hid_, H5F_OBJ_ALL | H5F_OBJ_LOCAL) - 1; // don't count the file itself
        if (count > 0) {
            throw error("closing HDF5 file would leave " + boost::lexical_cast<std::string>(count) + " open objects behind");
        }
    }
    if (plid_ >= 0) {
        H5Pclose(plid_);
    }
    if (H5Fclose(hid_) < 0) {
        throw error("closing HDF5 file: " + name() +
                    ", file ID: " + boost::lexical_cast<std::string>(hid_));
    }
    plid_ = -1;
    hid_ = -1;
}

inline std::string file::name() const
{
    if (hid_ < 0) {
        throw error("no HDF5 file associated to h5xx::file object");
    }
    ssize_t size = H5Fget_name(hid_, NULL, 0);        // determine string length
    if (size < 0) {
        throw error("retrieving name of HDF5 file with ID " + boost::lexical_cast<std::string>(hid_));
    }
    std::vector<char> buffer(size + 1);
    H5Fget_name(hid_, &*buffer.begin(), buffer.size()); // get string data
    return &*buffer.begin();
}

} // namespace h5xx

#endif // ! H5XX_FILE_HPP

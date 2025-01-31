/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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

/*! \brief Declares the i/o interface to H5MD HDF5 files.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 */

#ifndef GMX_FILEIO_H5MD_H
#define GMX_FILEIO_H5MD_H

#include "config.h" // To define GMX_USE_HDF5

#include <filesystem>
#include <optional>
#include <string>

enum class PbcType : int;

namespace gmx
{

// \brief ID for HDF5 files, groups, sets, etc., returned by the library.
typedef int64_t hid_t;
// \brief Status return values from some HDF5 library functions.
typedef int herr_t;

enum class H5mdFileMode : char
{
    Read  = 'r', //! Only read from the file.
    Write = 'w', //! Write to the file, replaces it if it exists. Also allows reading from the file.
    Append = 'a' //! Write to the file without replacing it if it exists. Also allows reading from the file.
};

/*! \brief Manager of an H5MD filehandle.
 * The class is designed to read/write data according to de Buyl et al., 2014
 * (https://doi.org/10.1016/j.cpc.2014.01.018) and https://www.nongnu.org/h5md/h5md.html
 */
class H5md
{
#if GMX_USE_HDF5
private:
    hid_t file_;            //!< The HDF5 identifier of the file. This is the H5MD root.
    H5mdFileMode filemode_; //!< Whether the file is open for reading ('r'), writing ('w') or appending ('a')
#endif

public:
    /*! \brief Open an H5MD file and manage its filehandle.
     *
     * \param[in] fileName    Name of the file to open. The same as the file path.
     * \param[in] mode        The mode to open the file.
     * \throws FileIOError if fileName is specified and the file cannot be opened.
     */
    H5md(const std::filesystem::path& fileName, const H5mdFileMode mode);

    ~H5md();

    H5md(const H5md&)            = delete;
    H5md& operator=(const H5md&) = delete;
    H5md(H5md&&)                 = delete;
    H5md& operator=(H5md&&)      = delete;

    /*! \brief Write all unwritten data to the file.
     * \param[in] throwExceptionUponError Whether to throw an exception if an error occurs.
     * Assumes a valid file_ identifier.
     * \throws FileIOError If there were errors during flushing (and throwExceptionUponError is true).
     */
    void flush(bool throwExceptionUponError = true);

    /*! \brief Set the author name attribute in the H5MD file
     * \param[in] authorName The author name.
     * \throws FileIOError If the author name attribute could not be set.
     */
    void setAuthor(const std::string& authorName);

    /*! \brief Get the author name attribute from the H5MD file
     * \returns the author name if the attribute was set.
     */
    std::optional<std::string> author();

    /*! \brief Set the name of the creating program as an attribute in the H5MD file.
     * \param[in] creatorName The creator name, i.e. the name of the program that created the file.
     * \throws FileIOError If the creator name attribute could not be set.
     */
    void setCreatorProgramName(const std::string& creatorName);

    /*! \brief Get the name of the creating program attribute from the H5MD file.
     * \returns The creator name, i.e. the name of the program that created the file,
     * if the attribute was set.
     */
    std::optional<std::string> creatorProgramName();

    /*! \brief Set the version of the creating program as an attribute in the H5MD file.
     * \param[in] version The version.
     * \throws FileIOError If the version attribute could not be set.
     */
    void setCreatorProgramVersion(const std::string& version);

    /*! \brief Get the version of the creating program attribute from the H5MD file.
     * \returns the version if the attribute was set.
     */
    std::optional<std::string> creatorProgramVersion();
};

} // namespace gmx
#endif // GMX_FILEIO_H5MD_H

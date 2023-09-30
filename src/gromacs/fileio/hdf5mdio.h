/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
/* This file was inspired by ch5md by Pierre de Buyl (BSD license). */

#ifndef GMX_FILEIO_HDF5MDIO_H
#define GMX_FILEIO_HDF5MDIO_H

// FIXME: TEMPORARY FOR EASIER EDITIING:
#define GMX_USE_HDF5 1

#include <string>

namespace h5xx
{
    // class attribute;
    // class datatype;
    class file;
}

class GmxHdf5MdElement
{
#ifdef GMX_USE_HDF5
private:
    // h5xx::attribute *group_;
    // h5xx::attribute *step_;
    // h5xx::attribute *time_;
    // h5xx::attribute *value_;
    // h5xx::datatype *datatype_;
    int isTime_;
    int currentStep_;
    // GmxHdf5MdParticlesGroup *particlesGroup_;
#endif
public:
    GmxHdf5MdElement();
    ~GmxHdf5MdElement();
    void append(void *data, int step, double time);
};

class GmxHdf5MdParticlesGroup
{
#ifdef GMX_USE_HDF5
private:
    // h5xx::attribute *group_;
    // GmxHdf5MdElement *position_;
    // h5xx::attribute *box_;
    // GmxHdf5MdElement *boxEdges_;
    // GmxHdf5MdElement *image_;
    // GmxHdf5MdElement *velocity_;
    // GmxHdf5MdElement *force_;
/*    GmxHdf5MdElement mass_;
    GmxHdf5MdElement species_;
    GmxHdf5MdElement id_;
    GmxHdf5MdElement charge_;*/
    // int localSizeMax;
#endif
public:
    GmxHdf5MdParticlesGroup();
    ~GmxHdf5MdParticlesGroup();
    void createBox(int dim, char *boundary[], bool isTime, double value[], GmxHdf5MdElement *link);
};

class GmxHdf5MdIo
{
#ifdef GMX_USE_HDF5
private:
    h5xx::file      *file_;
    // h5xx::attribute *particles_;
    // h5xx::attribute *observables_;
    // h5xx::attribute *parameters_;
#endif
public:
    GmxHdf5MdIo();
    /*! Construct a GmxHdf5MdIo object and open a GmxHdf5 file.
     *
     * \param[in] fileName    Name of the file to open. The same as the file path.
     * \param[in] modeString  The mode to open the file, described by a case-insensitive string of
     *                        letters, up to three characters long. Reading is always assumed.
     *                        'w' means writing,
     *                        't' means truncate, i.e., that existing files will be overwritten
     *                        'e' results in a failure if the file already exists.
     *                        All these modes can be combined.
     */
    GmxHdf5MdIo(const std::string &fileName, const std::string &modeString);
    ~GmxHdf5MdIo();

    /*! Open an GmxHdf5 file.
     *
     * \param[in] fileName    Name of the file to open. The same as the file path.
     * \param[in] modeString  The mode to open the file, described by a case-insensitive string of
     *                        letters, up to three characters long. Reading is always assumed.
     *                        'w' means writing,
     *                        't' means truncate, i.e., that existing files will be overwritten
     *                        'e' results in a failure if the file already exists.
     *                        All these modes can be combined.
     */
    void openFile(const std::string &fileName, const std::string &modeString);
    void closeFile();
    // GmxHdf5MdParticlesGroup createParticlesGroup(std::string name);
//     GmxHdf5MdElement createTimeData(hid_t loc, std::string name, int rank, int dims[], hid_t datatype, GmxHdf5MdElement *link);
//     GmxHdf5MdElement createFixedDataSimple(hid_t loc, std::string name, int rank, int dims[], hid_t datatype, void *data);
//     GmxHdf5MdElement createFixedDataScalar(hid_t loc, std::string name, hid_t datatype, void *data);
//     int writeStringAttribute(hid_t loc, std::string objectName, std::string attributeName, std::string value);
};

#endif // GMX_FILEIO_HDF5MDIO_H

/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
/*! \file
 * \brief
 * Declares gmx::ProgramContextInterface and related methods.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_PROGRAMCONTEXT_H
#define GMX_UTILITY_PROGRAMCONTEXT_H

namespace gmx
{

//! \addtogroup module_utility
//! \{

/*! \brief
 * Provides information about installation prefix (see
 * ProgramContextInterface::installationPrefix()).
 *
 * \inpublicapi
 */
struct InstallationPrefixInfo
{
    //! Initializes the structure with given values.
    InstallationPrefixInfo(const char *path, bool bSource)
        : path(path), bSourceLayout(bSource)
    {
    }

    /*! \brief
     * Path to the installation prefix of the current \Gromacs instance.
     *
     * If this is `NULL` or empty, data files cannot be looked up from the
     * install tree and \Gromacs functions that access such files may fail.
     * This can also contain a path to the source tree (see \a bSourceLayout).
     */
    const char *const path;
    /*! \brief
     * Whether \a path points to a source tree -like layout.
     *
     * For testing, it is useful to read data files from the source tree.
     * For such cases, the program context can return the source tree root path
     * in \a path, and set this to `true` to indicate that the data files
     * should be searched using the layout of the source tree instead of the
     * installation.
     */
    const bool        bSourceLayout;
};


/*! \brief
 * Provides context information about the program that is calling the library.
 *
 * This class provides access to information about the program that is
 * currently running.  This information is used to provide better information
 * to the user, and to locate the \Gromacs data files.
 *
 * setProgramContext() should be called before any other \Gromacs calls in
 * a program (except for gmx::init()).  This avoids thread safety issues, and
 * allows nicely formatted error messages.
 *
 * For thread safety, the implementations of the interface should ensure that
 * the returned strings remain valid until the end of the program (see
 * getProgramContext() for discussion on deinitialization).  Callers of the
 * interface within \Gromacs are prepared for exceptions, but generally
 * terminate the program on any exception.  Implementations should avoid
 * exception except for truly fatal conditions.
 *
 * The destructor is protected to signify that the context is never destroyed
 * through the interface.
 *
 * \see getProgramContext()
 * \see setProgramContext()
 * \inpublicapi
 */
class ProgramContextInterface
{
    public:
        /*! \brief
         * Returns the name of the binary as it was invoked without any path.
         *
         * This is typically `argv[0]` with any leading directory stripped.
         * Currently, this should be a valid file name.
         */
        virtual const char *programName() const = 0;
        /*! \brief
         * Returns a display name for the program.
         *
         * For simple programs, this can equal programName().  For the \Gromacs
         * `gmx` wrapper binary, this includes the name of the module (e.g.,
         * `gmx angle`).  This is used only for informational purposes, and
         * there are no constraints on contents, except that it should not be
         * `NULL`.
         */
        virtual const char *displayName() const = 0;
        /*! \brief
         * Returns the full path of the running binary.
         *
         * This is mainly used for informational purposes.  There are no
         * constraints on contents, except that it should not be `NULL`.
         * Currently, this is also used for sanity checks in checkpointing.
         *
         * The implementation can provide an empty string if the path to the
         * binary is not available.  In such a case, the information is not
         * shown.
         */
        virtual const char *fullBinaryPath() const = 0;
        /*! \brief
         * Returns the installation prefix for \Gromacs.
         *
         * This path is used to locate the data files that are in `share/top/`
         * in the source directory.
         * The implementation can provide an empty string if the path is not
         * available; in such a case, functions that require data files may
         * fail.
         *
         * The returned structure also contains a flag to indicate whether the
         * prefix actually points to the source tree.  This is used for tests
         * and to support running binaries directly from the build tree.
         */
        virtual InstallationPrefixInfo installationPrefix() const = 0;
        /*! \brief
         * Returns the full command line used to invoke the binary.
         *
         * The implementation can provide an empty string if no command line is
         * available.
         */
        virtual const char *commandLine() const = 0;

    protected:
        virtual ~ProgramContextInterface() {}
};

/*! \brief
 * Returns the global ProgramContextInterface instance.
 *
 * \returns The context set with setProgramContext().
 *
 * If nothing has been set with setProgramContext(), returns a default
 * implementation that returns `"GROMACS"` for the program and display names,
 * and empty strings for other values.
 * The default implementation never throws.
 *
 * Does not throw.
 *
 * See setProgramContext() for thread safety notes.  You should not call this
 * method in global deinitialization methods (e.g., destructors of global
 * variables), since it is very difficult to clean up the state correctly in
 * the presence of such calls.  For example, initForCommandLine() assumes that
 * such calls do not exist to be able to free the context before exiting.
 *
 * \see ProgramContextInterface
 */
const ProgramContextInterface &getProgramContext();
/*! \brief
 * Sets the global ProgramContextInterface instance.
 *
 * \param[in] context  Program context to set
 *     (can be NULL to restore the default context).
 *
 * The library does not take ownership of \p context.
 * The provided object must remain valid until the global instance is changed
 * by another call to setProgramContext().
 *
 * This method is not thread-safe.  It must be the first call to the library
 * after gmx::init(), and multi-threaded access is only supported after the
 * call completes.  If \Gromacs is getting called from multiple threads, or
 * uses multiple threads simultaneously, changing the program context is not
 * supported once it is set.
 * If the context is cleared at the end of the program, the caller must ensure
 * that all other threads have been terminated at this point.
 * These constraints simplify the implementation significantly.
 *
 * Does not throw.
 *
 * \see ProgramContextInterface
 */
void setProgramContext(const ProgramContextInterface *context);

//! \}

} // namespace gmx

#endif

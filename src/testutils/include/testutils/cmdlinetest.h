/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
/*! \libinternal \file
 * \brief
 * Declares utilities testing command-line programs.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_CMDLINETEST_H
#define GMX_TESTUTILS_CMDLINETEST_H

#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>

#include <gtest/gtest.h>

// arrayref.h is not strictly necessary for this header, but nearly all
// callers will need it to use the constructor that takes ArrayRef.
#include "gromacs/commandline/cmdlinemodulesettings.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

class ICommandLineModule;
class ICommandLineOptionsModule;

namespace test
{

class FloatingPointTolerance;
class IFileMatcherSettings;
class ITextBlockMatcherSettings;
class TestFileManager;
class TestReferenceChecker;

/*! \libinternal \brief
 * Helper class for tests that need to construct command lines.
 *
 * This class helps in writing tests for command-line handling.
 * The constructor method takes an array of const char pointers, specifying the
 * command-line arguments, each as one array element.  It is also possible to
 * construct the command line by adding individual arguments with append() and
 * addOption().
 * The argc() and argv() methods can then be used to obtain `argc` and `argv`
 * (non-const char pointers) arrays for passing into methods that expect these.
 *
 * Note that although the interface allows passing the argc and argv pointers
 * to methods that modify them (typically as \p f(&argc(), argv())), currently
 * the CommandLine object is not in a consistent state internally if the
 * parameters are actually modified.  Reading the command line is possible
 * afterwards, but modification is not.
 *
 * If you need to construct command lines that refer to files on the file
 * system, see CommandLineTestHelper and CommandLineTestBase for additional
 * convenience utilities.
 *
 * All constructors and methods that modify this class may throw an
 * std::bad_alloc.  Const methods and accessors do not throw.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class CommandLine
{
public:
    //! Initializes an empty command-line object.
    CommandLine();
    /*! \brief
     * Initializes a command-line object from an array.
     *
     * \param[in] cmdline  Array of command-line arguments.
     *
     * \p cmdline should include the binary name as the first element if
     * that is desired in the output.
     *
     * This constructor is not explicit to make it possible to create a
     * CommandLine object directly from a C array.
     */
    CommandLine(const ArrayRef<const char* const>& cmdline);
    //! \copydoc CommandLine(const ArrayRef<const char *const> &)
    CommandLine(const ArrayRef<const std::string>& cmdline);
    //! Creates a deep copy of a command-line object.
    CommandLine(const CommandLine& other);
    ~CommandLine();

    /*! \brief
     * Initializes a command-line object in-place from an array.
     *
     * \param[in] cmdline  Array of command-line arguments.
     *
     * \p cmdline should include the binary name as the first element if
     * that is desired in the output.
     *
     * This function does the same as the constructor that takes a
     * ArrayRef.  Any earlier contents of the object are discarded.
     *
     * Strong exception safety.
     */
    void initFromArray(const ArrayRef<const char* const>& cmdline);

    /*! \brief
     * Appends an argument to the command line.
     *
     * \param[in] arg  Argument to append.
     *
     * Strong exception safety.
     */
    void append(const char* arg);
    //! Convenience overload taking a std::string.
    void append(const std::string& arg) { append(arg.c_str()); }
    /*! \brief
     * Adds an option to the command line, typically a boolean.
     *
     * \param[in] name   Name of the option to append, which
     *                   should start with "-".
     */
    void addOption(const char* name);
    /*! \brief
     * Adds an option-value pair to the command line.
     *
     * \param[in] name   Name of the option to append, which
     *                   should start with "-".
     * \param[in] value  Value of the argument to append.
     */
    void addOption(const char* name, const char* value);
    //! Convenience overload taking a std::string.
    void addOption(const char* name, const std::string& value);
    //! Convenience overload taking a std::filesystem::path.
    void addOption(const char* name, const std::filesystem::path& value);
    //! Overload taking an int.
    void addOption(const char* name, int value);
    //! Overload taking a double.
    void addOption(const char* name, double value);
    /*! \brief
     * Appends all arguments from \p args to the command line.
     *
     * If the first argument of \p args does not start with a `-`, it is
     * skipped, assuming it is a gmx module name and thus useless.
     */
    void merge(const CommandLine& args);

    //! Returns argc for passing into C-style command-line handling.
    int& argc();
    //! Returns argv for passing into C-style command-line handling.
    char** argv();
    //! Returns argc for passing into C-style command-line handling.
    int argc() const;
    //! Returns argv for passing into C-style command-line handling.
    const char* const* argv() const;
    //! Returns a single argument.
    const char* arg(int i) const;

    //! Returns the command line formatted as a single string.
    std::string toString() const;

    //! Whether the command line contains the given option.
    bool contains(std::string_view name) const;

    //! The value of the argument following the given option, if found
    std::optional<std::string_view> argumentOf(std::string_view name) const;

private:
    class Impl;

    std::unique_ptr<Impl> impl_;
};

/*! \libinternal \brief
 * Helper class for tests that construct command lines that need to reference
 * existing files.
 *
 * This class provides helper methods for:
 *
 *   1. Adding input files to a CommandLine instance by generating them from a
 *      string provided in the test (setInputFileContents()).
 *   2. Adding output files to a CommandLine instance (setOutputFile()).
 *   3. Checking the contents of some of the output files using
 *      TestReferenceData (setOutputFile() and checkOutputFiles()).
 *   4. Static methods for easily executing command-line modules
 *      (various overloads of runModule()).
 *
 * All files created during the test are cleaned up at the end of the test.
 *
 * All methods can throw std::bad_alloc.
 *
 * \see TestFileManager
 * \inlibraryapi
 * \ingroup module_testutils
 */
class CommandLineTestHelper
{
public:
    /*! \brief
     * Runs a command-line program that implements ICommandLineModule.
     *
     * \param[in,out] module       Module to run.
     *     The function does not take ownership.
     * \param[in,out] commandLine  Command line parameters to pass.
     *     This is only modified if \p module modifies it.
     * \param[in]     settings     Settings for the command-line module
     * \returns The return value of the module.
     * \throws  unspecified  Any exception thrown by the module.
     */
    static int runModuleDirect(ICommandLineModule*         module,
                               CommandLine*                commandLine,
                               CommandLineModuleSettings&& settings = {});
    /*! \brief
     * Runs a command-line program that implements
     * ICommandLineOptionsModule.
     *
     * \param[in,out] module       Module to run.
     * \param[in,out] commandLine  Command line parameters to pass.
     *     This is only modified if \p module modifies it.
     * \param[in]     settings     Settings for the command-line module
     * \returns The return value of the module.
     * \throws  unspecified  Any exception thrown by the module.
     */
    static int runModuleDirect(std::unique_ptr<ICommandLineOptionsModule> module,
                               CommandLine*                               commandLine,
                               CommandLineModuleSettings&&                settings = {});
    /*! \brief
     * Runs a command-line program that implements
     * ICommandLineOptionsModule.
     *
     * \param[in] factory          Factory method for the module to run.
     * \param[in,out] commandLine  Command line parameters to pass.
     *     This is only modified if the module modifies it.
     * \param[in]     settings     Settings for the command-line module
     * \returns The return value of the module.
     * \throws  unspecified  Any exception thrown by the factory or the
     *     module.
     */
    static int runModuleFactory(const std::function<std::unique_ptr<ICommandLineOptionsModule>()>& factory,
                                CommandLine*                commandLine,
                                CommandLineModuleSettings&& settings = {});

    /*! \brief
     * Initializes an instance.
     *
     * \param  fileManager  File manager to use for generating temporary
     *     file names and to track temporary files.
     */
    explicit CommandLineTestHelper(TestFileManager* fileManager);
    ~CommandLineTestHelper();

    /*! \brief
     * Generates and sets an input file.
     *
     * \param[in,out] args      CommandLine to which to add the option.
     * \param[in]     option    Option to set.
     * \param[in]     extension Extension for the file to create.
     * \param[in]     contents  Text to write to the input file.
     *
     * Creates a temporary file with contents from \p contents, and adds
     * \p option to \p args with a value that points to the generated file.
     */
    void setInputFileContents(CommandLine*       args,
                              const char*        option,
                              const char*        extension,
                              const std::string& contents);
    /*! \brief
     * Generates and sets an input file.
     *
     * \param[in,out] args      CommandLine to which to add the option.
     * \param[in]     option    Option to set.
     * \param[in]     extension Extension for the file to create.
     * \param[in]     contents  Text to write to the input file.
     *
     * Creates a temporary file with contents from \p contents (each array
     * entry on its own line), and adds \p option to \p args with a value
     * that points to the generated file.
     */
    void setInputFileContents(CommandLine*                       args,
                              const char*                        option,
                              const char*                        extension,
                              const ArrayRef<const char* const>& contents);
    /*! \brief
     * Sets an output file parameter and adds it to the set of tested files.
     *
     * \param[in,out] args      CommandLine to which to add the option.
     * \param[in]     option    Option to set.
     * \param[in]     filename  Name of the output file.
     * \param[in]     matcher   Specifies how the contents of the file are
     *     tested.
     * \return        The name of the output file
     *
     * This method does the following:
     *  - Adds \p option to \p args to point a temporary file name
     *    constructed from \p filename.
     *  - Makes checkOutputFiles() to check the contents of the file
     *    against reference data, using \p matcher.
     *  - Marks the temporary file for removal at test teardown.
     *
     * \p filename is given to TestTemporaryFileManager to make a unique
     * filename for the temporary file.
     * If \p filename starts with a dot, a unique number is prefixed (such
     * that it is possible to create multiple files with the same extension
     * by just specifying the extension for every call of setOutputFile()).
     *
     * If the output file is needed to trigger some computation, or is
     * unconditionally produced by the code under test, but the contents
     * are not interesting for the test, use NoContentsMatch as the matcher.
     * Note that the existence of the output file is still verified.
     */
    std::string setOutputFile(CommandLine*                     args,
                              const char*                      option,
                              const char*                      filename,
                              const ITextBlockMatcherSettings& matcher);
    //! \copydoc setOutputFile(CommandLine *, const char *, const char *, const ITextBlockMatcherSettings &)
    std::string setOutputFile(CommandLine*                args,
                              const char*                 option,
                              const char*                 filename,
                              const IFileMatcherSettings& matcher);
    /*! \brief As for \c setOutputFile() but does not create an option
     *
     * Sometimes \c label and \c filename can be the same, depending
     * how the tool generates the filename. When the tool has
     * generated the filename from the name of another output file,
     * rather than as its own string, it can be necessary to prefix \c
     * filename with the temporary directory name so it can be found.
     *
     * This method is useful when a tool generates a series of output
     * files by modifying a common base name. Call this method with
     * every file name that is expected to be generated by the tool.
     * It returns the name of the output file */
    std::string setOutputFileWithGeneratedName(const std::string_view       label,
                                               const std::filesystem::path& filename,
                                               const IFileMatcherSettings&  matcher);

    /*! \brief
     * Checks output files added with setOutputFile() against reference
     * data.
     *
     * \param     checker  Reference data root location where the reference
     *     data is stored.
     *
     * The file contents are tested verbatim, using direct string
     * comparison.  The text can be found verbatim in the reference data
     * XML files for manual inspection.
     *
     * Generates non-fatal test failures if some output file contents do
     * not match the reference data.
     */
    void checkOutputFiles(TestReferenceChecker checker) const;

private:
    class Impl;

    std::unique_ptr<Impl> impl_;
};

/*! \libinternal \brief
 * Test fixture for tests that call a single command-line program with
 * input/output files.
 *
 * This class provides a convenient package for using CommandLineTestHelper in
 * a test that do not need special customization.  It takes care of creating
 * the other necessary objects (like TestFileManager, TestReferenceData, and
 * CommandLine) and wrapping the methods from CommandLineTestHelper such that
 * extra parameters are not needed.  Additionally, it provides setInputFile()
 * as a convenience function for adding a fixed input file, pointing to a file
 * that resides in the source tree.
 *
 * \see CommandLineTestHelper
 * \inlibraryapi
 * \ingroup module_testutils
 */
class CommandLineTestBase : public ::testing::Test
{
public:
    CommandLineTestBase();
    ~CommandLineTestBase() override;

    /*! \brief
     * Sets an input file.
     *
     * \param[in]     option    Option to set.
     * \param[in]     filename  Name of the input file.
     *
     * \see TestFileManager::getInputFilePath()
     */
    void setInputFile(const char* option, const char* filename);
    //! \copydoc setInputFile(const char *, const char *);
    void setInputFile(const char* option, const std::string& filename);
    /*! \brief
     * Sets an input file that may be modified. The file is copied to a
     * temporary file, which is used as the test input
     *
     * \param[in]     option    Option to set.
     * \param[in]     filename  Name of the input file.
     *
     */
    void setModifiableInputFile(const char* option, const char* filename);
    //! \copydoc setModifiableInputFile(const char *, const char *);
    void setModifiableInputFile(const char* option, const std::string& filename);
    /*! \brief
     * Generates and sets an input file.
     *
     * \see CommandLineTestHelper::setInputFileContents()
     */
    void setInputFileContents(const char* option, const char* extension, const std::string& contents);
    /*! \brief
     * Generates and sets an input file.
     *
     * \see CommandLineTestHelper::setInputFileContents()
     */
    void setInputFileContents(const char*                        option,
                              const char*                        extension,
                              const ArrayRef<const char* const>& contents);
    /*! \brief
     * Sets an output file whose name is passed via an option and adds it to the set of tested files.
     *
     * \see CommandLineTestHelper::setOutputFile()
     */
    std::string setOutputFile(const char* option, const char* filename, const ITextBlockMatcherSettings& matcher);
    //! \copydoc setOutputFile(const char *, const ITextBlockMatcherSettings &)
    std::string setOutputFile(const char* option, const char* filename, const IFileMatcherSettings& matcher);
    /*! \brief
     * Sets an output file whose name is generated by the tool and adds it to the set of tested files.
     *
     * Sometimes \c label and \c filename can be the same, depending
     * how the tool generates the filename. When the tool has
     * generated the filename from the name of another output file,
     * rather than as its own string, it can be necessary to prefix \c
     * filename with the temporary directory name so it can be found.
     *
     * This method is useful when a tool generates a series of output
     * files by modifying a common base name. Call this method with
     * every file name that is expected to be generated by the tool.
     *
     * \see CommandLineTestHelper::setOutputFileWithGeneratedName() */
    std::string setOutputFileWithGeneratedName(const std::string_view           label,
                                               const std::filesystem::path&     filename,
                                               const ITextBlockMatcherSettings& matcher);
    //! \copydoc setOutputFileWithGeneratedName(const std::string_view, const std::string_view, const ITextBlockMatcherSettings &)
    std::string setOutputFileWithGeneratedName(const std::string_view       label,
                                               const std::filesystem::path& filename,
                                               const IFileMatcherSettings&  matcher);
    /*! \brief
     * Sets a file parameter that is used for input and modified as output. The input file
     * is copied to a temporary file that is used as input and can be modified.
     */
    std::string setInputAndOutputFile(const char*                      option,
                                      const char*                      filename,
                                      const ITextBlockMatcherSettings& matcher);
    //! \copydoc setInputAndOutputFile(const char *, const char *, const ITextBlockMatcherSettings&);
    std::string setInputAndOutputFile(const char*                 option,
                                      const char*                 filename,
                                      const IFileMatcherSettings& matcher);

    /*! \brief
     * Returns the internal CommandLine object used to construct the
     * command line for the test.
     *
     * Derived test fixtures can use this to add additional options, and
     * to access the final command line to do the actual call that is being
     * tested.
     *
     * Does not throw.
     */
    CommandLine& commandLine();
    /*! \brief
     * Returns the internal TestFileManager object used to manage the
     * files.
     *
     * Derived test fixtures can use this to manage files in cases the
     * canned methods are not sufficient.
     *
     * Does not throw.
     */
    TestFileManager& fileManager();
    /*! \brief
     * Returns the root reference data checker.
     *
     * Derived test fixtures can use this to check other things than output
     * file contents.
     */
    TestReferenceChecker rootChecker();
    /*! \brief
     * Sets the tolerance for floating-point comparisons.
     *
     * All following floating-point comparisons using the checker will use
     * the new tolerance.
     *
     * Does not throw.
     */
    void setDefaultTolerance(const FloatingPointTolerance& tolerance);
    /*! \brief
     * Checks the output of writeHelp() against reference data.
     */
    void testWriteHelp(ICommandLineModule* module);
    /*! \brief
     * Checks output files added with setOutputFile() against reference
     * data.
     *
     * \see CommandLineTestHelper::checkOutputFiles()
     */
    void checkOutputFiles();

private:
    class Impl;

    std::unique_ptr<Impl> impl_;
};

} // namespace test
} // namespace gmx

#endif

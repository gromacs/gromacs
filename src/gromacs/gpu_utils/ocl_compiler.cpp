/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
 *  \brief Define infrastructure for OpenCL JIT compilation for Gromacs
 *
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 *
 * TODO Currently this file handles compilation of NBNXN kernels,
 * but e.g. organizing the defines for various physics models
 * is leaking in here a bit.
 */

#include "gmxpre.h"

#include "ocl_compiler.h"

#include "config.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <string>

#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

/*! \brief Path separator
 */
#define SEPARATOR '/'

/*! \brief Compiler options index
 */
typedef enum {
    b_invalid_option          = 0,
    b_amd_cpp,
    b_nvidia_verbose,
    b_generic_cl11,
    b_generic_cl12,
    b_generic_fast_relaxed_math,
    b_generic_noopt_compilation,
    b_generic_debug_symbols,
    b_amd_dump_temp_files,
    b_include_install_opencl_dir,
    b_include_source_opencl_dirs,
    b_num_build_options
} build_options_index_t;

/*! \brief List of available OpenCL compiler options
 */
static const char* build_options_list[] = {
    "",
    "-x clc++",                         /**< AMD C++ extension */
    "-cl-nv-verbose",                   /**< Nvidia verbose build log */
    "-cl-std=CL1.1",                    /**< Force CL 1.1  */
    "-cl-std=CL1.2",                    /**< Force CL 1.2  */
    "-cl-fast-relaxed-math",            /**< Fast math */
    "-cl-opt-disable",                  /**< Disable optimisations */
    "-g",                               /**< Debug symbols */
    "-save-temps"                       /**< AMD option to dump intermediate temporary
                                             files such as IL or ISA code */
};

/*! \brief Available sources
 */
static const char * kernel_filenames[] = {"nbnxn_ocl_kernels.cl"};

/*! \brief Defines to enable specific kernels based on vendor
 */
static const char * kernel_vendor_spec_definitions[] = {
    "-D_WARPLESS_SOURCE_",     /**< nbnxn_ocl_kernel_nowarp.clh  */
    "-D_NVIDIA_SOURCE_",       /**< nbnxn_ocl_kernel_nvidia.clh  */
    "-D_AMD_SOURCE_"           /**< nbnxn_ocl_kernel_amd.clh     */
};


/*! \brief Get the string of a build option of the specific id
 * \param  build_option_id  The option id as defines in the header
 * \return String containing the actual build option string for the compiler
 */
static const char* get_ocl_build_option(build_options_index_t build_option_id)
{
    if (build_option_id < b_num_build_options)
    {
        return build_options_list[build_option_id];
    }
    else
    {
        return build_options_list[b_invalid_option];
    }
}

/*! \brief Get the size of the string (without null termination) required
 *  for the build option of the specific id
 * \param  build_option_id  The option id as defines in the header
 * \return size_t containing the size in bytes of the build option string
 */
static size_t get_ocl_build_option_length(build_options_index_t build_option_id)
{

    if (build_option_id < b_num_build_options)
    {
        return strlen(build_options_list[build_option_id]);
    }
    else
    {
        return strlen(build_options_list[b_invalid_option]);
    }
}

/*! \brief Get the size of final composed build options literal
 *
 * \param build_device_vendor_id  Device vendor id. Used to
 *          automatically enable some vendor specific options
 * \param custom_build_options_prepend Prepend options string
 * \param custom_build_options_append  Append  options string
 * \return size_t containing the size in bytes of the composed
 *             build options string including null termination
 */
static size_t
create_ocl_build_options_length(
        ocl_vendor_id_t build_device_vendor_id,
        const char *    custom_build_options_prepend,
        const char *    custom_build_options_append)
{
    size_t build_options_length = 0;
    size_t whitespace           = 1;

    assert(build_device_vendor_id <= OCL_VENDOR_UNKNOWN);

    if (custom_build_options_prepend)
    {
        build_options_length +=
            strlen(custom_build_options_prepend)+whitespace;
    }

    if ( (build_device_vendor_id == OCL_VENDOR_AMD) && getenv("GMX_OCL_DEBUG") && getenv("GMX_OCL_FORCE_CPU") )
    {
        build_options_length += get_ocl_build_option_length(b_generic_debug_symbols)+whitespace;
    }

    if (getenv("GMX_OCL_NOOPT"))
    {
        build_options_length +=
            get_ocl_build_option_length(b_generic_noopt_compilation)+whitespace;
    }

    if (getenv("GMX_OCL_FASTMATH"))
    {
        build_options_length +=
            get_ocl_build_option_length(b_generic_fast_relaxed_math)+whitespace;
    }

    if ((build_device_vendor_id == OCL_VENDOR_NVIDIA) && getenv("GMX_OCL_VERBOSE"))
    {
        build_options_length +=
            get_ocl_build_option_length(b_nvidia_verbose) + whitespace;
    }

    if ((build_device_vendor_id == OCL_VENDOR_AMD) && getenv("GMX_OCL_DUMP_INTERM_FILES"))
    {
        /* To dump OpenCL build intermediate files, caching must be off */
        if (NULL != getenv("GMX_OCL_NOGENCACHE"))
        {
            build_options_length +=
                get_ocl_build_option_length(b_amd_dump_temp_files) + whitespace;
        }
    }

    if (custom_build_options_append)
    {
        build_options_length +=
            strlen(custom_build_options_append)+whitespace;
    }

    return build_options_length+1;
}

/*! \brief Get the size of final composed build options literal
 *
 * \param build_options_string The string where to save the
 *                                  resulting build options in
 * \param build_options_length The size of the build options
 * \param build_device_vendor_id  Device vendor id. Used to
 *          automatically enable some vendor specific options
 * \param custom_build_options_prepend Prepend options string
 * \param custom_build_options_append  Append  options string
 * \return The string build_options_string with the build options
 */
static char *
create_ocl_build_options(
        char *             build_options_string,
        size_t gmx_unused  build_options_length,
        ocl_vendor_id_t    build_device_vendor_id,
        const char *       custom_build_options_prepend,
        const char *       custom_build_options_append)
{
    size_t char_added = 0;

    if (custom_build_options_prepend)
    {
        strncpy( build_options_string+char_added,
                 custom_build_options_prepend,
                 strlen(custom_build_options_prepend));

        char_added += strlen(custom_build_options_prepend);
        build_options_string[char_added++] = ' ';
    }

    if (getenv("GMX_OCL_NOOPT") )
    {
        strncpy( build_options_string+char_added,
                 get_ocl_build_option(b_generic_noopt_compilation),
                 get_ocl_build_option_length(b_generic_noopt_compilation) );

        char_added += get_ocl_build_option_length(b_generic_noopt_compilation);
        build_options_string[char_added++] = ' ';

    }

    if (getenv("GMX_OCL_FASTMATH") )
    {
        strncpy( build_options_string+char_added,
                 get_ocl_build_option(b_generic_fast_relaxed_math),
                 get_ocl_build_option_length(b_generic_fast_relaxed_math) );

        char_added += get_ocl_build_option_length(b_generic_fast_relaxed_math);
        build_options_string[char_added++] = ' ';
    }

    if ((build_device_vendor_id == OCL_VENDOR_NVIDIA) && getenv("GMX_OCL_VERBOSE"))
    {
        strncpy(build_options_string + char_added,
                get_ocl_build_option(b_nvidia_verbose),
                get_ocl_build_option_length(b_nvidia_verbose));

        char_added += get_ocl_build_option_length(b_nvidia_verbose);
        build_options_string[char_added++] = ' ';
    }

    if ((build_device_vendor_id == OCL_VENDOR_AMD) && getenv("GMX_OCL_DUMP_INTERM_FILES"))
    {
        /* To dump OpenCL build intermediate files, caching must be off */
        if (NULL != getenv("GMX_OCL_NOGENCACHE"))
        {
            strncpy(build_options_string + char_added,
                    get_ocl_build_option(b_amd_dump_temp_files),
                    get_ocl_build_option_length(b_amd_dump_temp_files));

            char_added += get_ocl_build_option_length(b_amd_dump_temp_files);
            build_options_string[char_added++] = ' ';
        }
    }

    if ( ( build_device_vendor_id == OCL_VENDOR_AMD ) && getenv("GMX_OCL_DEBUG") && getenv("GMX_OCL_FORCE_CPU"))
    {
        strncpy( build_options_string+char_added,
                 get_ocl_build_option(b_generic_debug_symbols),
                 get_ocl_build_option_length(b_generic_debug_symbols) );

        char_added += get_ocl_build_option_length(b_generic_debug_symbols);
        build_options_string[char_added++] = ' ';
    }

    if (custom_build_options_append)
    {
        strncpy( build_options_string+char_added,
                 custom_build_options_append,
                 strlen(custom_build_options_append) );

        char_added += strlen(custom_build_options_append);
        build_options_string[char_added++] = ' ';
    }

    build_options_string[char_added++] = '\0';

    assert(char_added == build_options_length);

    return build_options_string;
}

/*! \brief Get the path to the main folder storing OpenCL kernels.
 *
 * By default, this function constructs the full path to the OpenCL from
 * the known location of the binary that is running, so that we handle
 * both in-source and installed builds. The user can override this
 * behavior by defining GMX_OCL_FILE_PATH environment variable.
 *
 * \return OS-normalized path string to the main folder storing OpenCL kernels
 *
 * \throws std::bad_alloc if out of memory.
 */
static std::string
get_ocl_root_path()
{
    const char *gmx_ocl_file_path;
    std::string ocl_root_path;

    /* Use GMX_OCL_FILE_PATH if the user has defined it */
    gmx_ocl_file_path = getenv("GMX_OCL_FILE_PATH");

    if (!gmx_ocl_file_path)
    {
        /* Normal way of getting ocl_root_dir. First get the right
           root path from the path to the binary that is running. */
        gmx::InstallationPrefixInfo info           = gmx::getProgramContext().installationPrefix();
        std::string                 dataPathSuffix = (info.bSourceLayout ?
                                                      "src/gromacs/mdlib/nbnxn_ocl" :
                                                      OCL_INSTALL_DIR);
        ocl_root_path = gmx::Path::join(info.path, dataPathSuffix);
    }
    else
    {
        ocl_root_path = gmx_ocl_file_path;
    }

    // Make sure we return an OS-correct path format
    return gmx::Path::normalize(ocl_root_path);
}

/*! \brief Get the size of the full kernel source file path and name
 *
 * The following full path size is computed:
 * strlen(ocl_root_path) + strlen(kernel_id.cl) + separator + null term
 *
 * \param kernel_src_id Id of the kernel source (auto,nvidia,amd,nowarp)
 * \return Size in bytes of the full kernel source file path and name including
 *          separators and null termination
 *
 * \throws std::bad_alloc if out of memory */
static size_t
get_ocl_kernel_source_file_info(kernel_source_index_t kernel_src_id)
{
    std::string ocl_root_path = get_ocl_root_path();

    if (ocl_root_path.empty())
    {
        return 0;
    }

    return (ocl_root_path.length() +                    /* Path to the main OpenCL folder*/
            1 +                                         /* Separator */
            strlen(kernel_filenames[kernel_src_id]) +   /* Kernel source file name */
            1                                           /* null char */
            );
}

/*! \brief Compose and the full path and name of the kernel src to be used
 *
 * \param ocl_kernel_filename   String where the full path and name will be saved
 * \param kernel_src_id         Id of the kernel source (default)
 * \param kernel_filename_len   Size of the full path and name string, as computed by get_ocl_kernel_source_file_info()
 * \return The ocl_kernel_filename complete with the full path and name; NULL if error.
 *
 * \throws std::bad_alloc if out of memory */
static char *
get_ocl_kernel_source_path(
        char *                  ocl_kernel_filename,
        kernel_source_index_t   kernel_src_id,
        size_t gmx_unused       kernel_filename_len)
{
    std::string ocl_root_path;

    assert(kernel_filename_len != 0);
    assert(ocl_kernel_filename != NULL);

    ocl_root_path = get_ocl_root_path();
    if (ocl_root_path.empty())
    {
        return NULL;
    }

    size_t chars_copied = 0;
    strncpy(ocl_kernel_filename, ocl_root_path.c_str(), ocl_root_path.length());
    chars_copied += ocl_root_path.length();

    ocl_kernel_filename[chars_copied++] = SEPARATOR;

    strncpy(&ocl_kernel_filename[chars_copied],
            kernel_filenames[kernel_src_id],
            strlen(kernel_filenames[kernel_src_id]) );
    chars_copied += strlen(kernel_filenames[kernel_src_id]);

    ocl_kernel_filename[chars_copied++] = '\0';

    assert(chars_copied == kernel_filename_len);

    return ocl_kernel_filename;
}

/* Undefine the separators */
#undef SEPARATOR

/*! \brief Loads the src inside the file filename onto a string in memory
 *
 * \param filename The name of the file to be read
 * \param p_source_length Pointer to the size of the source in bytes
 *                          (without null termination)
 * \return A string with the contents of the file with name filename,
 *  or NULL if there was a problem opening/reading the file
 */
static char*
load_ocl_source(const char* filename, size_t* p_source_length)
{
    FILE * filestream = NULL;
    char * ocl_source;
    size_t source_length;

    source_length = 0;

    if (!filename)
    {
        return NULL;
    }

    filestream    = fopen(filename, "rb");
    if (!filestream)
    {
        return NULL;
    }

    fseek(filestream, 0, SEEK_END);
    source_length = ftell(filestream);
    fseek(filestream, 0, SEEK_SET);

    ocl_source = (char*)malloc(source_length + 1);
    if (fread(ocl_source, source_length, 1, filestream) != 1)
    {
        fclose(filestream);
        free(ocl_source);
        return 0;
    }

    fclose(filestream);
    ocl_source[source_length] = '\0';

    *p_source_length = source_length;
    return ocl_source;
}

/*! \brief Handles the dumping of the OpenCL JIT compilation log
 *
 * In a debug build:
 *  -Success: Save to file kernel_id.SUCCEEDED in the run folder.
 *  -Fail   : Save to file kernel_id.FAILED in the run folder.
 *            Dump to stderr
 * In a release build:
 *  -Success: Nothing is logged.
 *  -Fail   : Save to a file kernel_id.FAILED in the run folder.
 * If GMX_OCL_DUMP_LOG is set, log is always dumped to file
 * If OCL_JIT_DUMP_STDERR is set, log is always dumped to stderr
 *
 * \param build_log String containing the OpenCL JIT compilation log
 * \param build_options_string String containing the options used for the build
 * \param build_status The OpenCL type status of the build (CL_SUCCESS etc)
 * \param kernel_src_id The id of the kernel src used for the build (default)
 *
 * \throws std::bad_alloc if out of memory */
static void
handle_ocl_build_log(
        const char        *   build_log,
        const char        *   build_options_string,
        cl_int                build_status,
        kernel_source_index_t kernel_src_id)
{
    bool dumpStdErr = false;
    bool dumpFile;
#ifdef NDEBUG
    dumpFile   = (build_status != CL_SUCCESS);
#else
    dumpFile   = true;
    if (build_status != CL_SUCCESS)
    {
        dumpStdErr = true;
    }
#endif

    /* Override default handling */
    if (getenv("GMX_OCL_DUMP_LOG") != NULL)
    {
        dumpFile = true;
    }
    if (getenv("OCL_JIT_DUMP_STDERR") != NULL)
    {
        dumpStdErr = true;
    }

    if (dumpFile || dumpStdErr)
    {
        FILE       *build_log_file       = NULL;
        const char *fail_header          = "Compilation of source file failed! \n";
        const char *success_header       = "Compilation of source file was successful! \n";
        const char *log_header           = "--------------LOG START---------------\n";
        const char *log_footer           = "---------------LOG END----------------\n";
        char       *build_info;
        std::string log_fname;

        build_info = (char*)malloc(32 + strlen(build_options_string) );
        sprintf(build_info, "-- Used build options: %s\n", build_options_string);

        if (dumpFile)
        {
            log_fname = gmx::formatString("%s.%s", kernel_filenames[kernel_src_id],
                                          (build_status == CL_SUCCESS) ? "SUCCEEDED" : "FAILED");
            build_log_file = fopen(log_fname.c_str(), "w");
        }

        size_t complete_message_size = 0;
        char * complete_message;


        complete_message_size  =  (build_status == CL_SUCCESS) ? strlen(success_header) : strlen(fail_header);
        complete_message_size += strlen(build_info) + strlen(log_header) + strlen(log_footer);
        complete_message_size += strlen(build_log);
        complete_message_size += 1; //null termination
        complete_message       = (char*)malloc(complete_message_size);

        sprintf(complete_message, "%s%s%s%s%s",
                (build_status == CL_SUCCESS) ? success_header : fail_header,
                build_info,
                log_header,
                build_log,
                log_footer);

        if (dumpFile)
        {
            if (build_log_file)
            {
                fprintf(build_log_file, "%s", complete_message);
            }

            printf("The OpenCL compilation log has been saved in \"%s\"\n", log_fname.c_str());
        }
        if (dumpStdErr)
        {
            if (build_status != CL_SUCCESS)
            {
                fprintf(stderr, "%s", complete_message);
            }
        }
        if (build_log_file)
        {
            fclose(build_log_file);
        }

        free(complete_message);
        free(build_info);
    }
}

/*!  \brief Get the warp size reported by device
 *
 *  This is platform implementation dependant and seems to only work on the Nvidia and Amd platforms!
 *  Nvidia reports 32, Amd for GPU 64. Ignore the rest
 *
 *  \param  context   Current OpenCL context
 *  \param  device_id OpenCL device with the context
 *  \return cl_int value of the warp size
 */
static cl_int
ocl_get_warp_size(cl_context context, cl_device_id device_id)
{
    cl_int      cl_error     = CL_SUCCESS;
    size_t      warp_size    = 0;
    const char *dummy_kernel = "__kernel void test(__global int* test){test[get_local_id(0)] = 0;}";

    cl_program  program =
        clCreateProgramWithSource(context, 1, (const char**)&dummy_kernel, NULL, &cl_error);

    cl_error =
        clBuildProgram(program, 0, NULL, NULL, NULL, NULL);

    cl_kernel kernel = clCreateKernel(program, "test", &cl_error);

    cl_error = clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                        sizeof(size_t), &warp_size, NULL);

    clReleaseKernel(kernel);
    clReleaseProgram(program);

    assert(warp_size != 0);
    assert(cl_error == CL_SUCCESS);
    return warp_size;

}

/*! \brief Automatically select vendor-specific kernel from vendor id
 *
 * \param vendor_id Vendor id enumerator (amd,nvidia,intel,unknown)
 * \return Vendor-specific kernel version
 */
static kernel_vendor_spec_t
ocl_autoselect_kernel_from_vendor(ocl_vendor_id_t vendor_id)
{
    kernel_vendor_spec_t kernel_vendor;
#ifndef NDEBUG
    printf("Selecting kernel source automatically\n");
#endif
    switch (vendor_id)
    {
        case OCL_VENDOR_AMD:
            kernel_vendor = amd_vendor_kernels;
            printf("Selecting kernel for AMD\n");
            break;
        case OCL_VENDOR_NVIDIA:
            kernel_vendor = nvidia_vendor_kernels;
            printf("Selecting kernel for NVIDIA\n");
            break;
        default:
            kernel_vendor = generic_vendor_kernels;
            printf("Selecting generic kernel\n");
            break;
    }
    return kernel_vendor;
}

/*! \brief Returns the compiler define string needed to activate vendor-specific kernels
 *
 * \param kernel_spec Kernel vendor specification
 * \return String with the define for the spec
 */
static const char *
ocl_get_vendor_specific_define(kernel_vendor_spec_t kernel_spec)
{
    assert(kernel_spec < auto_vendor_kernels );
#ifndef NDEBUG
    printf("Setting up kernel vendor spec definitions:  %s \n", kernel_vendor_spec_definitions[kernel_spec]);
#endif
    return kernel_vendor_spec_definitions[kernel_spec];
}

/*! \brief Check if there's a valid cache available, and return it if so
 *
 * \param[in]  ocl_binary_filename   Name of file containing the binary cache
 * \param[in]  build_options_string  Compiler command-line options to use (currently unused)
 * \param[in]  ocl_source            NULL-terminated string of OpenCL source code (currently unused)
 * \param[out] ocl_binary_size       Size of the binary file once loaded in memory
 * \param[out] ocl_binary            Pointer to the binary file bytes (valid only if return is true)
 * \return                           Whether the file reading was successful
 *
 * \todo Compare current build options and code against the build
 * options and the code corresponding to the cache. If any change is
 * detected this function must return false.
 */
bool
check_ocl_cache(char            *ocl_binary_filename,
                char gmx_unused *build_options_string,
                char gmx_unused *ocl_source,
                size_t          *ocl_binary_size,
                unsigned char  **ocl_binary)
{
    FILE  *f;
    size_t read_count;

    f = fopen(ocl_binary_filename, "rb");
    if (!f)
    {
        return false;
    }

    fseek(f, 0, SEEK_END);
    *ocl_binary_size = ftell(f);
    *ocl_binary      = (unsigned char*)malloc(*ocl_binary_size);
    fseek(f, 0, SEEK_SET);
    read_count = fread(*ocl_binary, 1, *ocl_binary_size, f);
    fclose(f);

    if (read_count != (*ocl_binary_size))
    {
        return false;
    }

    return true;
}

/*! \brief Builds a string with build options for the OpenCL kernels
 *
 * \throws std::bad_alloc if out of memory */
char*
ocl_get_build_options_string(cl_context           context,
                             cl_device_id         device_id,
                             kernel_vendor_spec_t kernel_vendor_spec,
                             ocl_vendor_id_t      ocl_device_vendor,
                             const char *         defines_for_kernel_types,
                             const char *         runtime_consts)
{
    char * build_options_string               = NULL;
    char   custom_build_options_prepend[1024] = { 0 };
    char  *custom_build_options_append        = NULL;
    cl_int warp_size = 0;

    /* Get the reported warp size. Compile a small dummy kernel to do so */
    warp_size = ocl_get_warp_size(context, device_id);

    /* Select vendor specific kernels automatically */
    if (kernel_vendor_spec == auto_vendor_kernels)
    {
        kernel_vendor_spec = ocl_autoselect_kernel_from_vendor(ocl_device_vendor);
    }

    /* Create include paths for kernel sources.
       All OpenCL kernel files are expected to be stored in one single folder. */
    {
        /* Apple does not seem to accept the quoted include paths other
         * OpenCL implementations are happy with. Since the standard still says
         * it should be quoted, we handle Apple as a special case.
         */
#ifdef __APPLE__
        std::string unescaped_ocl_root_path = get_ocl_root_path();
        std::string ocl_root_path;

        char        incl_opt_start[] = "-I";
        char        incl_opt_end[]   = "";

        for (std::string::size_type i = 0; i < unescaped_ocl_root_path.length(); i++)
        {
            if (unescaped_ocl_root_path[i] == ' ')
            {
                ocl_root_path.push_back('\\');
            }
            ocl_root_path.push_back(unescaped_ocl_root_path[i]);
        }
        // Here the Apple ocl_root_path has all spaces prepended with a backslash
#else
        std::string ocl_root_path = get_ocl_root_path();

        char        incl_opt_start[] = "-I\"";
        char        incl_opt_end[]   = "\"";

#endif
        size_t      chars            = 0;

        custom_build_options_append =
            (char*)calloc((ocl_root_path.length()   /* Path to the OpenCL folder */
                           + strlen(incl_opt_start) /* -I" */
                           + strlen(incl_opt_end)   /* " */
                           + 1                      /* null char */
                           ), 1);

        strncpy(&custom_build_options_append[chars], incl_opt_start, strlen(incl_opt_start));
        chars += strlen(incl_opt_start);

        strncpy(&custom_build_options_append[chars], ocl_root_path.c_str(), ocl_root_path.length());
        chars += ocl_root_path.length();

        strncpy(&custom_build_options_append[chars], incl_opt_end, strlen(incl_opt_end));
    }

    /* Get vendor specific define (amd,nvidia,nowarp) */
    const char * kernel_vendor_spec_define =
        ocl_get_vendor_specific_define(kernel_vendor_spec);

    /* Compose the build options to be prepended. */
    sprintf(custom_build_options_prepend,
            "-DWARP_SIZE_TEST=%d %s %s %s",
            warp_size,
            kernel_vendor_spec_define,
            defines_for_kernel_types,
            runtime_consts ? runtime_consts : ""
            );

    /* Get the size of the complete build options string */
    size_t build_options_length =
        create_ocl_build_options_length(
                ocl_device_vendor,
                custom_build_options_prepend,
                custom_build_options_append
                );

    build_options_string = (char *)malloc(build_options_length);

    /* Compose the complete build options */
    create_ocl_build_options(
            build_options_string,
            build_options_length,
            ocl_device_vendor,
            custom_build_options_prepend,
            custom_build_options_append
            );

    if (custom_build_options_append)
    {
        free(custom_build_options_append);
    }

    return build_options_string;
}

/*! \brief Implement caching of OpenCL binaries
 *
 * \param[in] program     Index of program to cache
 * \param[in] file_name  Name of file to use for the cache
 */
void
print_ocl_binaries_to_file(cl_program program, char* file_name)
{
    size_t         ocl_binary_size = 0;
    unsigned char *ocl_binary      = NULL;

    clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &ocl_binary_size, NULL);

    ocl_binary = (unsigned char*)malloc(ocl_binary_size);

    clGetProgramInfo(program, CL_PROGRAM_BINARIES, sizeof(unsigned char *), &ocl_binary, NULL);

    FILE *f = fopen(file_name, "wb");
    fwrite(ocl_binary, 1, ocl_binary_size, f);
    fclose(f);

    free(ocl_binary);
}

/*! \brief Compile the kernels as described by kernel src id and vendor spec
 *
 * \param[in]  kernel_source_file        Index of the kernel src to be used (default)
 * \param[in]  kernel_vendor_spec        Vendor-specific compilation (auto,nvidia,amd,nowarp)
 * \param[in]  defines_for_kernel_types  Preprocessor defines that trigger the compilation of the kernels
 * \param[out] result_str                Gromacs error string
 * \param[in]  context                   Current context on the device to compile for
 * \param[in]  device_id                 OpenCL device id of the device to compile for
 * \param[in]  ocl_device_vendor         Enumerator of the device vendor to compile for
 * \param[out] p_program                 Pointer to the cl_program where the compiled
 *                                       cl_program will be stored
 * \param[in]  runtime_consts            Optional string with runtime constants.
 *                                       Each constant is given according to the following
 *                                       format: "-Dname=value".
 *                                       Multiple defines are separated by blanks.
 *
 * \return cl_int with the build status AND any other OpenCL error appended to it
 *
 * \todo Consider whether we can parallelize the compilation of all
 * the kernels by compiling them in separate programs - but since the
 * resulting programs can't refer to each other, that might lead to
 * bloat of util code?
 *
 * \throws std::bad_alloc if out of memory
 */
cl_int
ocl_compile_program(
        kernel_source_index_t kernel_source_file,
        kernel_vendor_spec_t  kernel_vendor_spec,
        const char *          defines_for_kernel_types,
        char *                result_str,
        cl_context            context,
        cl_device_id          device_id,
        ocl_vendor_id_t       ocl_device_vendor,
        cl_program *          p_program,
        const char *          runtime_consts
        )
{
    char         * build_options_string   = NULL;
    cl_int         cl_error               = CL_SUCCESS;

    char         * ocl_source              = NULL;
    size_t         ocl_source_length       = 0;
    size_t         kernel_filename_len     = 0;

    bool           bCacheOclBuild           = false;
    bool           bOclCacheValid           = false;

    char           ocl_binary_filename[256] = { 0 };
    size_t         ocl_binary_size          = 0;
    unsigned char *ocl_binary               = NULL;

    /* Load OpenCL source files */
    {
        char* kernel_filename = NULL;

        /* Get the size of the kernel source filename */
        kernel_filename_len = get_ocl_kernel_source_file_info(kernel_source_file);
        if (kernel_filename_len)
        {
            kernel_filename = (char*)malloc(kernel_filename_len);
        }

        /* Get the actual full path and name of the source file with the kernels */
        get_ocl_kernel_source_path(kernel_filename, kernel_source_file, kernel_filename_len);

        /* Load the above source file and store its contents in ocl_source */
        ocl_source = load_ocl_source(kernel_filename, &ocl_source_length);

        if (!ocl_source)
        {
            sprintf(result_str, "Error loading OpenCL code %s", kernel_filename);
            return CL_BUILD_PROGRAM_FAILURE;
        }

        /* The sources are loaded so the filename is not needed anymore */
        free(kernel_filename);
    }

    /* Allocate and initialize the string with build options */
    build_options_string =
        ocl_get_build_options_string(context, device_id, kernel_vendor_spec,
                                     ocl_device_vendor,
                                     defines_for_kernel_types,
                                     runtime_consts);

    /* Check if OpenCL caching is ON - currently caching is disabled
       until we resolve concurrency issues. */
    /* bCacheOclBuild = (NULL == getenv("GMX_OCL_NOGENCACHE"));*/
    if (bCacheOclBuild)
    {
        clGetDeviceInfo(device_id, CL_DEVICE_NAME, sizeof(ocl_binary_filename), ocl_binary_filename, NULL);
        strcat(ocl_binary_filename, ".bin");

        /* Check if there's a valid cache available */
        bOclCacheValid = check_ocl_cache(ocl_binary_filename,
                                         build_options_string,
                                         ocl_source,
                                         &ocl_binary_size, &ocl_binary);
    }

    /* Create OpenCL program */
    if (bCacheOclBuild && bOclCacheValid)
    {
        /* Create program from pre-built binaries */
        *p_program =
            clCreateProgramWithBinary(
                    context,
                    1,
                    &device_id,
                    &ocl_binary_size,
                    (const unsigned char**)&ocl_binary,
                    NULL,
                    &cl_error);
    }
    else
    {
        /* Create program from source code */
        *p_program =
            clCreateProgramWithSource(
                    context,
                    1,
                    (const char**)(&ocl_source),
                    &ocl_source_length,
                    &cl_error
                    );
    }

    /* Build program */
    cl_int build_status         = CL_SUCCESS;
    {
        /* Now we are ready to launch the build */
        build_status =
            clBuildProgram(*p_program, 0, NULL, build_options_string, NULL, NULL);

        if (build_status == CL_SUCCESS)
        {
            if (bCacheOclBuild)
            {
                /* If OpenCL caching is ON, but the current cache is not
                   valid => update it */
                if (!bOclCacheValid)
                {
                    print_ocl_binaries_to_file(*p_program, ocl_binary_filename);
                }
            }
            else
            if ((OCL_VENDOR_NVIDIA == ocl_device_vendor) && getenv("GMX_OCL_DUMP_INTERM_FILES"))
            {
                /* If dumping intermediate files has been requested and this is an NVIDIA card
                   => write PTX to file */
                char ptx_filename[256];

                clGetDeviceInfo(device_id, CL_DEVICE_NAME, sizeof(ptx_filename), ptx_filename, NULL);
                strcat(ptx_filename, ".ptx");

                print_ocl_binaries_to_file(*p_program, ptx_filename);
            }
        }

        // Get log string size
        size_t build_log_size       = 0;
        cl_error =
            clGetProgramBuildInfo(
                    *p_program,
                    device_id,
                    CL_PROGRAM_BUILD_LOG,
                    0,
                    NULL,
                    &build_log_size
                    );

        /* Regardless of success or failure, if there is something in the log
         *  we might need to display it */
        if (build_log_size && (cl_error == CL_SUCCESS) )
        {
            char *build_log = NULL;

            /* Allocate memory to fit the build log,
                it can be very large in case of errors */
            build_log = (char*)malloc(build_log_size);

            if (build_log)
            {
                /* Get the actual compilation log */
                cl_error =
                    clGetProgramBuildInfo(
                            *p_program,
                            device_id,
                            CL_PROGRAM_BUILD_LOG,
                            build_log_size,
                            build_log,
                            NULL
                            );

                /* Save or display the log */
                if (!cl_error)
                {
                    handle_ocl_build_log(
                            build_log,
                            build_options_string,
                            build_status,
                            kernel_source_file
                            );
                }

                /* Build_log not needed anymore */
                free(build_log);
            }
        }
    }

    /*  Final clean up */
    if (ocl_binary)
    {
        free(ocl_binary);
    }

    if (build_options_string)
    {
        free(build_options_string);
    }

    if (ocl_source)
    {
        free(ocl_source);
    }

    /* Append any other error to the build_status */
    return build_status | cl_error;
}

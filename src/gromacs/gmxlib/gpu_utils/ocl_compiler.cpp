/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

#include "gromacs/utility/stringutil.h"

/* This path is defined by CMake and it depends on the install prefix option.
   The opencl kernels are installed in bin/opencl.*/
#if !defined(OCL_INSTALL_FULL_PATH)
#pragma error "OCL_INSTALL_FULL_PATH has not been defined"
#endif

#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__))
#define SEPARATOR '/'
#define SSEPARATOR "/"
#elif defined(_WIN32)
#define SEPARATOR '\\'
#define SSEPARATOR "\\"
#endif

/*! \brief Compiler options index
 */
typedef enum {
    _invalid_option_          = 0,
    _amd_cpp_,
    _nvidia_verbose_,
    _generic_cl11_,
    _generic_cl12_,
    _generic_fast_relaxed_math_,
    _generic_noopt_compilation_,
    _generic_debug_symbols_,
    _amd_dump_temp_files_,
    _include_install_opencl_dir_,
    _include_source_opencl_dirs_,
    _num_build_options_
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
    "-save-temps",                      /**< AMD option to dump intermediate temporary
                                             files such as IL or ISA code */
    "-I" "\""OCL_INSTALL_FULL_PATH "\"" /**< Include path to kernel sources */
    /*,
       "-I../../src/gromacs/gmxlib/ocl_tools           -I../../src/gromacs/mdlib/nbnxn_ocl            -I../../src/gromacs/pbcutil            -I../../src/gromacs/mdlib"
       -I../../../gromacs/src/gromacs/gmxlib/ocl_tools -I../../../gromacs/src/gromacs/mdlib/nbnxn_ocl -I../../../gromacs/src/gromacs/pbcutil -I../../../gromacs/src/gromacs/mdlib" */
};

/*! \brief Include paths when using the GMX_OCL_FILE_PATH to point to the gromacs source tree */
#define INCLUDE_PATH_COUNT 4

/*! \brief List of include paths relative to the path defined by \p GMX_OCL_FILE_PATH
 */
static const char* include_path_list[] =
{
    "gromacs" SSEPARATOR "mdlib" SSEPARATOR "nbnxn_ocl",
    "gromacs" SSEPARATOR "gmxlib" SSEPARATOR "ocl_tools",
    "gromacs" SSEPARATOR "pbcutil",
    "gromacs" SSEPARATOR "mdlib"
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

/*! \brief Array of the defines needed to generate a specific eel flavour
 */
static const char * kernel_electrostatic_family_definitions[] =
{
    "-DEL_CUTOFF -D_EELNAME=_ElecCut",
    "-DEL_RF -D_EELNAME=_ElecRF",
    "-DEL_EWALD_TAB -D_EELNAME=_ElecEwQSTab",
    "-DEL_EWALD_TAB -DLJ_CUTOFF_CHECK -D_EELNAME=_ElecEwQSTabTwinCut",
    "-DEL_EWALD_ANA -D_EELNAME=_ElecEw",
    "-DEL_EWALD_ANA -DLJ_CUTOFF_CHECK -D_EELNAME=_ElecEwTwinCut"
};

/*! \brief Array of the defines needed to generate a specific vdw flavour
 */
static const char * kernel_VdW_family_definitions[] =
{
    "-D_VDWNAME=_VdwLJ",
    "-DLJ_FORCE_SWITCH -D_VDWNAME=_VdwLJFsw",
    "-DLJ_POT_SWITCH -D_VDWNAME=_VdwLJPsw",
    "-DLJ_EWALD_COMB_GEOM -D_VDWNAME=_VdwLJEwCombGeom",
    "-DLJ_EWALD_COMB_LB -D_VDWNAME=_VdwLJEwCombLB"
};


/*! \brief Get the string of a build option of the specific id
 * \param  build_option_id  The option id as defines in the header
 * \return String containing the actual build option string for the compiler
 */
static const char* get_ocl_build_option(build_options_index_t build_option_id)
{
    if (build_option_id < _num_build_options_)
    {
        return build_options_list[build_option_id];
    }
    else
    {
        return build_options_list[_invalid_option_];
    }
}

/*! \brief Get the size of the string (without null termination) required
 *  for the build option of the specific id
 * \param  build_option_id  The option id as defines in the header
 * \return size_t containing the size in bytes of the build option string
 */
static size_t get_ocl_build_option_length(build_options_index_t build_option_id)
{

    if (build_option_id < _num_build_options_)
    {
        return strlen(build_options_list[build_option_id]);
    }
    else
    {
        return strlen(build_options_list[_invalid_option_]);
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

    assert(build_device_vendor_id <= _OCL_VENDOR_UNKNOWN_);

    if (custom_build_options_prepend)
    {
        build_options_length +=
            strlen(custom_build_options_prepend)+whitespace;
    }

    if ( (build_device_vendor_id == _OCL_VENDOR_AMD_) && getenv("GMX_OCL_DEBUG") && getenv("GMX_OCL_FORCE_CPU") )
    {
        build_options_length += get_ocl_build_option_length(_generic_debug_symbols_)+whitespace;
    }

    if (getenv("GMX_OCL_NOOPT"))
    {
        build_options_length +=
            get_ocl_build_option_length(_generic_noopt_compilation_)+whitespace;
    }

    if (getenv("GMX_OCL_FASTMATH"))
    {
        build_options_length +=
            get_ocl_build_option_length(_generic_fast_relaxed_math_)+whitespace;
    }

    if ((build_device_vendor_id == _OCL_VENDOR_NVIDIA_) && getenv("GMX_OCL_VERBOSE"))
    {
        build_options_length +=
            get_ocl_build_option_length(_nvidia_verbose_) + whitespace;
    }

    if ((build_device_vendor_id == _OCL_VENDOR_AMD_) && getenv("GMX_OCL_DUMP_INTERM_FILES"))
    {
        /* To dump OpenCL build intermediate files, caching must be off */
        if (NULL != getenv("GMX_OCL_NOGENCACHE"))
        {
            build_options_length +=
                get_ocl_build_option_length(_amd_dump_temp_files_) + whitespace;
        }
    }

    build_options_length +=
        get_ocl_build_option_length(_include_install_opencl_dir_)+whitespace;

    /*build_options_length +=
        get_ocl_build_option_length(_include_source_opencl_dirs_)+whitespace; */

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
                 get_ocl_build_option(_generic_noopt_compilation_),
                 get_ocl_build_option_length(_generic_noopt_compilation_) );

        char_added += get_ocl_build_option_length(_generic_noopt_compilation_);
        build_options_string[char_added++] = ' ';

    }

    if (getenv("GMX_OCL_FASTMATH") )
    {
        strncpy( build_options_string+char_added,
                 get_ocl_build_option(_generic_fast_relaxed_math_),
                 get_ocl_build_option_length(_generic_fast_relaxed_math_) );

        char_added += get_ocl_build_option_length(_generic_fast_relaxed_math_);
        build_options_string[char_added++] = ' ';
    }

    if ((build_device_vendor_id == _OCL_VENDOR_NVIDIA_) && getenv("GMX_OCL_VERBOSE"))
    {
        strncpy(build_options_string + char_added,
                get_ocl_build_option(_nvidia_verbose_),
                get_ocl_build_option_length(_nvidia_verbose_));

        char_added += get_ocl_build_option_length(_nvidia_verbose_);
        build_options_string[char_added++] = ' ';
    }

    if ((build_device_vendor_id == _OCL_VENDOR_AMD_) && getenv("GMX_OCL_DUMP_INTERM_FILES"))
    {
        /* To dump OpenCL build intermediate files, caching must be off */
        if (NULL != getenv("GMX_OCL_NOGENCACHE"))
        {
            strncpy(build_options_string + char_added,
                    get_ocl_build_option(_amd_dump_temp_files_),
                    get_ocl_build_option_length(_amd_dump_temp_files_));

            char_added += get_ocl_build_option_length(_amd_dump_temp_files_);
            build_options_string[char_added++] = ' ';
        }
    }

    if ( ( build_device_vendor_id == _OCL_VENDOR_AMD_ ) && getenv("GMX_OCL_DEBUG") && getenv("GMX_OCL_FORCE_CPU"))
    {
        strncpy( build_options_string+char_added,
                 get_ocl_build_option(_generic_debug_symbols_),
                 get_ocl_build_option_length(_generic_debug_symbols_) );

        char_added += get_ocl_build_option_length(_generic_debug_symbols_);
        build_options_string[char_added++] = ' ';
    }

    strncpy( build_options_string+char_added,
             get_ocl_build_option(_include_install_opencl_dir_),
             get_ocl_build_option_length(_include_install_opencl_dir_)
             );
    char_added += get_ocl_build_option_length(_include_install_opencl_dir_);
    build_options_string[char_added++] = ' ';

    /* strncpy( build_options_string+char_added,
              get_ocl_build_option(_include_source_opencl_dirs_),
              get_ocl_build_option_length(_include_source_opencl_dirs_)
       );
       char_added += get_ocl_build_option_length(_include_source_opencl_dirs_);
       build_options_string[char_added++]=' '; */

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

/*! \brief Get the size of the full kernel source file path and name
 *
 * If GMX_OCL_FILE_PATH is defined in the environment the following full path size is returned:
 *  strlen($GMX_OCL_FILE_PATH) + strlen(kernel_id.cl) + separator + null term
 * Otherwise the following full path size is returned (OCL_INSTALL_FULL_PATH is provided by CMAKE
 *  installation prefix path) :
 *  strlen( OCL_INSTALL_FULL_PATH ) + strlen(kernel_id.cl) + separator + null term
 *
 * \param kernel_src_id Id of the kernel source (auto,nvidia,amd,nowarp)
 * \return Size in bytes of the full kernel source file path and name including
 *          separators and null termination
 */
static size_t
get_ocl_kernel_source_file_info(kernel_source_index_t kernel_src_id)
{
    char * gmx_pathname = NULL;
    if ( (gmx_pathname = getenv("GMX_OCL_FILE_PATH")) != NULL)
    {
        return (strlen(gmx_pathname)                        /* Path to gromacs source    */
                + strlen(kernel_filenames[kernel_src_id])   /* Kernel source filename   */
                + strlen(include_path_list[0])              /* Path from gmx src to kernel src */
                + 2                                         /* separator and null char  */
                );
    }
    else
    {
        /* Kernel source are located in the installation folder of Gromacs */
        return( strlen(kernel_filenames[kernel_src_id])     /* Kernel source filename     */
                + strlen(OCL_INSTALL_FULL_PATH)             /* Full path to kernel source */
                + 2);                                       /* separator and null char    */
    }
}

/*! \brief Compose and the full path and name of the kernel src to be used
 *
 * If GMX_OCL_FILE_PATH is defined in the environment the following full path size is composed:
 *  $GMX_OCL_FILE_PATH/kernel_id.cl
 * Otherwise the following full path is composed (OCL_INSTALL_FULL_PATH is provided by CMAKE
 *  installation prefix path):
 *  OCL_INSTALL_FULL_PATH/kernel_id.cl
 *
 * \param ocl_kernel_filename   String where the full path and name will be saved
 * \param kernel_src_id         Id of the kernel source (default)
 * \param kernel_filename_len   Size of the full path and name string
 * \return The ocl_kernel_filename complete with the full path and name
 */
static char *
get_ocl_kernel_source_path(
        char *                  ocl_kernel_filename,
        kernel_source_index_t   kernel_src_id,
        size_t gmx_unused       kernel_filename_len)
{
    char *filepath = NULL;

    assert(kernel_filename_len != 0);
    assert(ocl_kernel_filename != NULL);

    if ( (filepath = getenv("GMX_OCL_FILE_PATH")) != NULL)
    {
        FILE  *file_ok = NULL;

        size_t chars_copied = 0;
        strncpy(ocl_kernel_filename, filepath, strlen(filepath));
        chars_copied += strlen(filepath);

        strncpy(&ocl_kernel_filename[chars_copied],
                include_path_list[0],
                strlen(include_path_list[0]));
        chars_copied += strlen(include_path_list[0]);

        ocl_kernel_filename[chars_copied++] = SEPARATOR;

        strncpy(&ocl_kernel_filename[chars_copied],
                kernel_filenames[kernel_src_id],
                strlen(kernel_filenames[kernel_src_id]));
        chars_copied += strlen(kernel_filenames[kernel_src_id]);

        ocl_kernel_filename[chars_copied++] = '\0';

        assert(chars_copied == kernel_filename_len);


        //Try to open the file to check that it exists
        file_ok = fopen(ocl_kernel_filename, "rb");
        if (file_ok)
        {
            fclose(file_ok);
        }
        else
        {
            printf("Warning, you seem to have misconfigured the GMX_OCL_FILE_PATH environent variable: %s\n",
                   filepath);
        }
    }
    else
    {
        size_t chars_copied = 0;
        strncpy(ocl_kernel_filename, OCL_INSTALL_FULL_PATH, strlen(OCL_INSTALL_FULL_PATH));
        chars_copied += strlen(OCL_INSTALL_FULL_PATH);

        ocl_kernel_filename[chars_copied++] = SEPARATOR;

        strncpy(&ocl_kernel_filename[chars_copied],
                kernel_filenames[kernel_src_id],
                strlen(kernel_filenames[kernel_src_id]) );
        chars_copied += strlen(kernel_filenames[kernel_src_id]);

        ocl_kernel_filename[chars_copied++] = '\0';

        assert(chars_copied == kernel_filename_len);

    }
    return ocl_kernel_filename;
}

/* Undefine the separators */
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__))
#undef SEPARATOR
#elif defined(_WIN32)
#undef SEPARATOR
#endif

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
 */
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
        FILE      * build_log_file       = NULL;
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
            fclose(build_log_file);
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
    printf("Selecting kernel source automatically\n");
    switch (vendor_id)
    {
        case _OCL_VENDOR_AMD_:
            kernel_vendor = _amd_vendor_kernels_;
            printf("Selecting kernel for AMD\n");
            break;
        case _OCL_VENDOR_NVIDIA_:
            kernel_vendor = _nvidia_vendor_kernels_;
            printf("Selecting kernel for Nvidia\n");
            break;
        default:
            kernel_vendor = _generic_vendor_kernels_;
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
    assert(kernel_spec < _auto_vendor_kernels_ );
    printf("Setting up kernel vendor spec definitions:  %s \n", kernel_vendor_spec_definitions[kernel_spec]);
    return kernel_vendor_spec_definitions[kernel_spec];
}

/*! \brief Populates algo_defines with the compiler defines required to avoid all flavour generation
 *
 * For example if flavour eelOclRF with evdwOclFSWITCH, the output will be such that the corresponding
 * kernel flavour is generated:
 * -D_OCL_FASTGEN_ (will replace flavour generator kernels.clh with the fastgen one kernels_fastgen.clh)
 * -DEL_RF             (The eelOclRF flavour)
 * -D_EELNAME=_ElecRF  (The first part of the generated kernel name )
 * -DLJ_EWALD_COMB_GEOM (The evdwOclFSWITCH flavour)
 * -D_VDWNAME=_VdwLJEwCombGeom (The second part of the generated kernel name )
 * prune/energy are still generated as originally. It is only the the flavour-level that has changed, so that
 * only the required flavour for the simulation is compiled.
 *
 * \param p_gmx_algo_family    Pointer to algo_family structure (eel,vdw)
 * \param p_algo_defines       String to populate with the defines
 */
static void
ocl_get_fastgen_define(
        gmx_algo_family_t * p_gmx_algo_family,
        char *              p_algo_defines)
{
    printf("Setting up kernel fastgen definitions: ");
    sprintf(p_algo_defines, "-D_OCL_FASTGEN_ %s %s ",
            kernel_electrostatic_family_definitions[p_gmx_algo_family->eeltype],
            kernel_VdW_family_definitions[p_gmx_algo_family->vdwtype]
            );
    printf(" %s \n", p_algo_defines);
}

/*! \brief TODO */
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

    // TODO: Compare current build options and code against the builds options
    // and the code corresponding to the cache. If any change is detected this
    // function must return false.
    return true;
}

/*! \brief TODO */
char*
ocl_get_build_options_string(cl_context           context,
                             cl_device_id         device_id,
                             kernel_vendor_spec_t kernel_vendor_spec,
                             ocl_vendor_id_t      ocl_device_vendor,
                             gmx_algo_family_t  * p_gmx_algo_family,
                             int                  DoFastGen)
{
    char * build_options_string              = NULL;
    char   custom_build_options_prepend[512] = { 0 };
    char  *custom_build_options_append       = NULL;
    char  *oclpath   = NULL;
    cl_int warp_size = 0;

    /* Get the reported warp size. Compile a small dummy kernel to do so */
    warp_size = ocl_get_warp_size(context, device_id);

    /* Select vendor specific kernels automatically */
    if (kernel_vendor_spec == _auto_vendor_kernels_)
    {
        kernel_vendor_spec = ocl_autoselect_kernel_from_vendor(ocl_device_vendor);
    }

    /* Create include paths for non-standard location of the kernel sources */
    if ((oclpath = getenv("GMX_OCL_FILE_PATH")) != NULL)
    {
        size_t chars = 0;
        custom_build_options_append =
            (char*)calloc((strlen(oclpath) + 32)*INCLUDE_PATH_COUNT, 1);

        for (int i = 0; i < INCLUDE_PATH_COUNT; i++)
        {
            strncpy(&custom_build_options_append[chars], "-I", strlen("-I"));
            chars += strlen("-I");

            strncpy(&custom_build_options_append[chars], oclpath, strlen(oclpath));
            chars += strlen(oclpath);

            strncpy(
                    &custom_build_options_append[chars],
                    include_path_list[i],
                    strlen(include_path_list[i])
                    );
            chars += strlen(include_path_list[i]);

            strncpy(&custom_build_options_append[chars], " ", 1);
            chars += 1;
        }
        printf("GMX_OCL_FILE_PATH includes: %s\n", custom_build_options_append);
    }

    /* Get vendor specific define (amd,nvidia,nowarp) */
    const char * kernel_vendor_spec_define =
        ocl_get_vendor_specific_define(kernel_vendor_spec);

    /* Use the fastgen flavour-level kernel generator instead of the original */
    char kernel_fastgen_define[128] = { 0 };
    if (DoFastGen)
    {
        ocl_get_fastgen_define(p_gmx_algo_family, kernel_fastgen_define);
    }

    /* Compose the build options to be prepended */
    sprintf(custom_build_options_prepend,
            "-DWARP_SIZE_TEST=%d %s %s",
            warp_size,
            kernel_vendor_spec_define,
            kernel_fastgen_define
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

/*! \brief TODO */
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
 * \param kernel_source_file Index of the kernel src to be used (default)
 * \param kernel_vendor_spec Vendor specific compilation (auto,nvidia,amd,nowarp)
 * \param p_gmx_algo_family  Flavour of kernels requested
 * \param DoFastGen          Use the info on the required flavour of kernels to reduce
 *                            JIT compilation time
 * \param result_str         Gromacs error string
 * \param context            Current context on the device to compile for
 * \param device_id          OpenCL device id of the device to compile for
 * \param ocl_device_vendor  Enumerator of the device vendor to compile for
 * \param p_program          Pointer to the cl_program where the compiled
 *                            cl_program will be stored
 * \return cl_int with the build status AND any other OpenCL error appended to it
 */
cl_int
ocl_compile_program(
        kernel_source_index_t       kernel_source_file,
        kernel_vendor_spec_t        kernel_vendor_spec,
        gmx_algo_family_t *         p_gmx_algo_family,
        int                         DoFastGen,
        char *                      result_str,
        cl_context                  context,
        cl_device_id                device_id,
        ocl_vendor_id_t             ocl_device_vendor,
        cl_program *                p_program
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
                                     ocl_device_vendor, p_gmx_algo_family,
                                     DoFastGen);

    /* Check if OpenCL caching is ON */
    bCacheOclBuild = (NULL == getenv("GMX_OCL_NOGENCACHE"));
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
            if ((_OCL_VENDOR_NVIDIA_ == ocl_device_vendor) && getenv("GMX_OCL_DUMP_INTERM_FILES"))
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

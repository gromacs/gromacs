/* ************************************************************************
 * Copyright 2014 Advanced Micro Devices, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * ************************************************************************/

#ifndef __CLFFT_CLFFT_BINARY_LOOKUP__
#define __CLFFT_CLFFT_BINARY_LOOKUP__

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <string>
#include <vector>

#include "generator.h"
#include "plan.h"

//
// FFTBinaryLookup defines an API to manage the kernel cache on the disk
//
// The FFTBinaryLookup object provides methods to:
//  * check if a cache file exists on the disk or not
//  * fill-up the signature to characterize the program beeing built on the disk
//  * build a cl_program from a string kernel or from a binary
// 
// A cache entry is a file stored on the disk which contains 3 sections:
//  * A header section (providing information about file structure)
//  * The binary contained in the cl_program
//  * A signature which provides additionnal informations about the kernel
//    and allows to characterize the kernel in the disk cache
//
// The environment variable CLFFT_CACHE_PATH defines the location of the
// cache on the disk. If the variable CLFFT_CACHE_PATH is not defined, no
// cache file is written on the disk, but the cl_program can be built and
// remains on memory
//
// Concerning multithreading, the policy is that every thread build the
// cl_program from the source, but only the first one writes it on the
// disk. Other threads continue with the cl_program in memory.
//
// A typical cache query shall be composed of the following steps:
//
//  (1) Create a local instance of FFTBinaryLookup 
//
//  (2) Specify the additional characteristics (i.e. variants) of the
//      requested program. Those information combined with the program
//      name and the OpenCL context and device shall form a unique
//      signature for the binary program.
// 
//  (3) Perform the effective search by calling the 'found' method
// 
//  (4) if the search was successfull then cl_program can be retreived 
//      by a call to the 'getProgram' method
// 
//  (5) if the search was not successfull then a cl_program 
//      must be created  and populated in the cache by a call 
//      to the 'setProgram' method.
// 
//  (6) Destroy the FFTBinaryLookup local instance.
// 
// For instance, that could be 
//
//     cl_program program  ;
//   
//     The program name is part of the signature and shall be unique 
//     const char * program_name = "... my unique program name ... " ;
//   
//     FFTBinaryLookup bl(context, device, program_name);
//   
//     //  Specify additionnal information used to build a
//     //  signature signature for that cache entry 
//                  
//     bl.variantInt( vectorSize );
//     bl.variantInt( hasBorder );
//     ... 
//   
//     // Perform the query 
//     if ( bl.found() ) 
//     {
//        //   Success! use the cl_program retreived from the cache
//        program = bl.getProgram();
//     }
//     else 
//     {
//        //   Failure! we need to build the program ourself
//        program = build_the_program(context,device,vectorSize,...) ; 
//        //   Inform the lookup object of the program
//        bl.setProgram(program);  
//        //   And populate the cache
//        bl.populateCache() 
//     }
// 
// Remark: The members buildFromSource, buildFromBinary etc are utility 
//         functions that can be used to build the cl_program from either 
//         sources or binary (e.g. SPIR). Their use is optionnal. 
//
//
class FFTBinaryLookup
 {
public:
    // Constructor
    // \param ctxt the context for which the cl_program should be built
    // \param device the device for which the cl_program should be built
    // \param kernel_name the kernel identifier
    FFTBinaryLookup(const clfftGenerators gen, const clfftPlanHandle plHandle, cl_context ctxt, cl_device_id device);
    ~FFTBinaryLookup();

    // Methods to fill up the signature of the cache entry
    void variantInt(int num);
    void variantDouble(double num);
    void variantCompileOptions(const std::string & opts);
    void variantRaw(const void * data, size_t bytes);

    // Indicates whether or not the cache entry was found on the disk
    // If the cache entry was found and is complete on the disk, its content
    // is loaded
    // \return true if a cache entry was found, false else
    bool found();

    // Build a cl_program from the source code and init attributes
    // of the current structure
    // so that the program can be accessed with the getProgram method
    // Write the file to the cache
    cl_int buildFromSource(const char * source);

    // Build a cl_program from the source code and init attributes
    // so that the program can be accessed with the getProgram method
    // Write the file to the cache
    cl_int buildFromBinary(const void * data, 
                           size_t len);

    // Returns the cl_program built from binary or loaded from disk
    cl_program getProgram();

    // Set the current m_program to the given program
    void setProgram(cl_program program, const char * source);

    // Build a cl_program from a text
    static cl_program buildProgramFromSource(const char * filename,
                                             cl_context context,
                                             cl_device_id device,
                                             cl_int & err,
                                             const char * options = 0);

    // Build a cl_program from binary
    static cl_program buildProgramFromBinary(const char * data,
                                             size_t data_size,
                                             cl_context context,
                                             cl_device_id device,
                                             cl_int & err,
                                             const char * options = 0);

    // Initialize the whole cache file information (magic_key, header and program)
    // and dump on the disk
    cl_int populateCache();

private:

    // Serialize variants and compute the checksum to load the file from cache
    void finalizeVariant();

    // Build a cl_program from the source code and init attributes
    // so that the program can be accessed with the getProgram method
    // Do not write the file to the cache
    cl_int buildFromLoadedBinary(const void * data, 
                                 size_t len);

    // Try to retrieve the header of the cache file
    // Returns: ok if the header sections was successfully loaded, false else
    bool loadHeader(std::ifstream &file, size_t length);

    // Try to retrieve the cl_program and its signature in file
    // Returns: ok if the binary and signature sections were successfully loaded, false else
    bool loadBinaryAndSignature(std::ifstream &file);

    // Try to create a file associated to the current program/variant in the cache folder
    // Returns true if the file was successfully opened and loaded, false else
    bool tryLoadCacheFile();

    // Dump the file on the disk with the name stored in this->m_cache_entry_name
    cl_int writeCacheFile(std::vector<unsigned char*> &data);

    // Retrieve device name, device vendor and driver number by calling
    // clGetDeviceInfo
    cl_int retrieveDeviceAndDriverInfo();

    // Cache entry name 
    std::string m_cache_entry_name;

    // Path for the cache entry name
    std::string m_path;

    // Header structure of a cache entry
    typedef struct Header_
    {
        char magic_key[4]; // = |C|L|F|\0, useful to know that we are loading a clfft cache entry
        size_t whole_file_size; // the whole file of the size to know if the current file is complete or not
        size_t header_size; // = sizeof(Header)
        size_t binary_size; // kernel binary size
        size_t signature_size; // variant information
    } Header;

    Header m_header;

    cl_context   m_context;
    cl_device_id m_device;

    cl_program   m_program;

    std::string  m_source;

    unsigned char * m_binary;
    char *          m_signature;

    enum VariantKind
    {
        INT,
        DOUBLE,
        STRING,
        DATA
    };

    struct Variant
    {
        Variant();
        Variant(VariantKind kind, char * data, size_t size);
		Variant(const Variant &obj);
		Variant &operator=(const Variant &obj);


        ~Variant();

        VariantKind m_kind;
        size_t      m_size;
        char *      m_data;
        
    };

    // Cache entry, useful to abstract Windows and linux 
    // cache entry file descriptor
    struct CacheEntry
    {
        CacheEntry(const std::string & filename);
        bool exclusive_create();
        void close();
        bool successful_creation();

    private:
        std::string m_filename;
        bool        m_successful_creation;
        void *      m_handle;
    };

    // Variants
    std::vector<Variant> m_variants;

    // Indicates whether the cache should be used or not
    bool m_cache_enabled;
};

#undef SIZE

#endif

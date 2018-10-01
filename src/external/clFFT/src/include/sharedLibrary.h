/* ************************************************************************
 * Copyright 2013 Advanced Micro Devices, Inc.
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


#pragma once
#ifndef _SHAREDLIBRARY_H_
#define _SHAREDLIBRARY_H_
#include <string>

//	_WIN32 is defined for both 32 & 64 bit environments
#if defined( _WIN32 )
	#define WIN32_LEAN_AND_MEAN			// Exclude rarely-used stuff from Windows headers
	// Windows Header Files:
	#include <windows.h>
#else
	#include <dlfcn.h>
#endif

inline void* LoadSharedLibrary( std::string unixPrefix, std::string libraryName, bool quiet )
{
#if defined( _WIN32 )
	libraryName += ".dll";

	//	HMODULE is actually the load address; function returns NULL if it cannot find the shared library
	HMODULE fileHandle	= ::LoadLibraryExA( libraryName.c_str( ), NULL, NULL );
#elif defined(__linux__) || defined(__GNU__) || (defined(__FreeBSD_kernel__) && defined(__GLIBC__))
        tstring linuxName = unixPrefix;
	linuxName += libraryName += ".so";
	void* fileHandle = ::dlopen( linuxName.c_str( ), RTLD_NOW );
	if( !quiet && !fileHandle )
	{
		std::cerr << ::dlerror( ) << std::endl;
	}
#elif defined(__APPLE__)
  tstring appleName = unixPrefix;
  appleName += libraryName += ".dylib";
  void* fileHandle = ::dlopen( appleName.c_str( ), RTLD_NOW );
  if( !quiet && !fileHandle )
  {
          std::cerr << ::dlerror( ) << std::endl;
  }
#elif defined(__FreeBSD__)
        tstring freebsdName = unixPrefix;
        freebsdName += libraryName += ".so";
        void* fileHandle = ::dlopen( freebsdName.c_str( ), RTLD_NOW );
        if( !quiet && !fileHandle )
        {
                std::cerr << ::dlerror( ) << std::endl;
        }
#else
        #error "unsupported platform"
#endif

	return fileHandle;
}

//	If the function succeeds, the return value is nonzero.
//	If the function fails, the return value is zero.
inline int FreeSharedLibrary( void*& libHandle )
{
	int result	= 0;

#if defined( _WIN32 )
	if( libHandle != 0 )
		result = ::FreeLibrary( reinterpret_cast< HMODULE >( libHandle ) );
#else
	if( libHandle != 0 )
		result = ( ::dlclose( libHandle ) == 0 );
#endif

	libHandle	= NULL;

	return result;
}

//	This takes a shared module handle returned from LoadSharedLibrary, and a text string of a symbol
//	to load from the module, and returns a pointer to that symbol.  If the symbol is not found, NULL
//	is returned.  If the module handle is NULL, NULL is returned.
inline void* LoadFunctionAddr( void* libHandle, std::string funcName )
{
	if( libHandle == NULL )
		return NULL;

#if defined( _WIN32 )
	HMODULE fileHandle = reinterpret_cast< HMODULE >( libHandle );

	void* pFunc	= reinterpret_cast< void* >( ::GetProcAddress( fileHandle, funcName.c_str( ) ) );
#else
	void* pFunc = ::dlsym( libHandle, funcName.c_str( ) );
#endif

	return pFunc;
}

#endif // _SHAREDLIBRARY_H_

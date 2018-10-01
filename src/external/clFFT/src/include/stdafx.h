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


// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <memory>
#include <vector>
#include <valarray>
#include <cstring>
#include <stdarg.h>
#include <assert.h>
#include <complex>

//	_WIN32 is defined for both 32 & 64 bit environments
#if defined( _WIN32 )
	#include <tchar.h>
	#include "targetver.h"

#if !defined( NOMINMAX )
	#define NOMINMAX
#endif

    #define WIN32_LEAN_AND_MEAN			// Exclude rarely-used stuff from Windows headers
	// Windows Header Files:
	#include <windows.h>
#endif

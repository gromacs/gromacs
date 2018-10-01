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


// StatTimer.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "statisticalTimer.extern.h"
#include "statisticalTimer.CPU.h"
#include "statisticalTimer.GPU.h"

//	Even though the individual getInstance functions of the timer classes return references,
//	we convert those to pointers before returning from here so that the clients can initialize
//	their local variables to NULL, which refernces do not allow.
baseStatTimer* getStatTimer( const clfftTimerType type )
{
	if( type == CLFFT_CPU )
		return	&CpuStatTimer::getInstance( );

	return	&GpuStatTimer::getInstance( );
}

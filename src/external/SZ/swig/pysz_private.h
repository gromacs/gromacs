/**
 *
 * A set of python bindings for SZ 
 * 
 * Developed by Robert Underwood while he was at Clemson University
 * This material is based upon work supported by the National Science 
 * Foundation under Grant No. 1633608.
 * 
 * Copyright © 2019 Robert Underwood
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY Robert Underwood ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL Robert Underwood BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation
 * are those of the authors and should not be interpreted as representing
 * official policies, either expressed or implied, of Robert Underwood.
 */

#pragma once
#include "sz.h"

template <class T>
class SZTypeToTypeID
{};

#define MAKE_TYPE_TO_ID(type, id)                                              \
  template <>                                                                  \
  class SZTypeToTypeID<type>                                                   \
  {                                                                            \
  public:                                                                      \
    static const int value = id;                                               \
  };

MAKE_TYPE_TO_ID(float, SZ_FLOAT);
MAKE_TYPE_TO_ID(double, SZ_DOUBLE);
MAKE_TYPE_TO_ID(uint8_t, SZ_UINT8);
MAKE_TYPE_TO_ID(int8_t, SZ_INT8);
MAKE_TYPE_TO_ID(uint16_t, SZ_UINT16);
MAKE_TYPE_TO_ID(int16_t, SZ_INT16);
MAKE_TYPE_TO_ID(uint32_t, SZ_UINT32);
MAKE_TYPE_TO_ID(int32_t, SZ_INT32);
MAKE_TYPE_TO_ID(uint64_t, SZ_UINT64);
MAKE_TYPE_TO_ID(int64_t, SZ_INT64);

//TODO add the other types as they are supported by the dispatch method

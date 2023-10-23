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

%module pysz
%feature("autodoc", 1);

%{
#define SWIG_FILE_WITH_INIT
#define SWIG_PYTHON_STRICT_BYTE_CHAR
#include "pysz.h"
#include "defines.h"
#include <vector>
#include <cstdint>
%}

%init %{
import_array();
%}

%include <std_vector.i>
%include <std_string.i>
%include "numpy.i"

namespace std {
  %template(vectori8) vector<int8_t>;
  %template(vectori16) vector<int16_t>;
  %template(vectori32) vector<int>;
  %template(vectori64) vector<int64_t>;
  %template(vectorui8) vector<uint8_t>;
  %template(vectorui16) vector<uint16_t>;
  %template(vectorui32) vector<uint32_t>;
  %template(vectorui64) vector<uint64_t>;
  %template(vectorf) vector<float>;
  %template(vectord) vector<double>;
};

%apply (float* INPLACE_ARRAY1, int DIM1 ) {(float* data, size_t r1)}
%apply (float* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(float* data, size_t r1, size_t r2)}
%apply (float* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(float* data, size_t r1, size_t r2, size_t r3)}
%apply (float* INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4 ) {(float* data, size_t r1, size_t r2, size_t r3, size_t r4)}
%apply (double* INPLACE_ARRAY1, int DIM1 ) {(double* data, size_t r1)}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(double* data, size_t r1, size_t r2)}
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(double* data, size_t r1, size_t r2, size_t r3)}
%apply (double* INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4 ) {(double* data, size_t r1, size_t r2, size_t r3, size_t r4)}
%apply (uint8_t* INPLACE_ARRAY1, int DIM1 ) {(uint8_t* data, size_t r1)}
%apply (uint8_t* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(uint8_t* data, size_t r1, size_t r2)}
%apply (uint8_t* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(uint8_t* data, size_t r1, size_t r2, size_t r3)}
%apply (uint8_t* INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4 ) {(uint8_t* data, size_t r1, size_t r2, size_t r3, size_t r4)}
%apply (int8_t* INPLACE_ARRAY1, int DIM1 ) {(int8_t* data, size_t r1)}
%apply (int8_t* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(int8_t* data, size_t r1, size_t r2)}
%apply (int8_t* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(int8_t* data, size_t r1, size_t r2, size_t r3)}
%apply (int8_t* INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4 ) {(int8_t* data, size_t r1, size_t r2, size_t r3, size_t r4)}
%apply (uint16_t* INPLACE_ARRAY1, int DIM1 ) {(uint16_t* data, size_t r1)}
%apply (uint16_t* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(uint16_t* data, size_t r1, size_t r2)}
%apply (uint16_t* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(uint16_t* data, size_t r1, size_t r2, size_t r3)}
%apply (uint16_t* INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4 ) {(uint16_t* data, size_t r1, size_t r2, size_t r3, size_t r4)}
%apply (int16_t* INPLACE_ARRAY1, int DIM1 ) {(int16_t* data, size_t r1)}
%apply (int16_t* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(int16_t* data, size_t r1, size_t r2)}
%apply (int16_t* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(int16_t* data, size_t r1, size_t r2, size_t r3)}
%apply (int16_t* INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4 ) {(int16_t* data, size_t r1, size_t r2, size_t r3, size_t r4)}
%apply (uint32_t* INPLACE_ARRAY1, int DIM1 ) {(uint32_t* data, size_t r1)}
%apply (uint32_t* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(uint32_t* data, size_t r1, size_t r2)}
%apply (uint32_t* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(uint32_t* data, size_t r1, size_t r2, size_t r3)}
%apply (uint32_t* INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4 ) {(uint32_t* data, size_t r1, size_t r2, size_t r3, size_t r4)}
%apply (int32_t* INPLACE_ARRAY1, int DIM1 ) {(int32_t* data, size_t r1)}
%apply (int32_t* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(int32_t* data, size_t r1, size_t r2)}
%apply (int32_t* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(int32_t* data, size_t r1, size_t r2, size_t r3)}
%apply (int32_t* INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4 ) {(int32_t* data, size_t r1, size_t r2, size_t r3, size_t r4)}
%apply (uint64_t* INPLACE_ARRAY1, int DIM1 ) {(uint64_t* data, size_t r1)}
%apply (uint64_t* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(uint64_t* data, size_t r1, size_t r2)}
%apply (uint64_t* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(uint64_t* data, size_t r1, size_t r2, size_t r3)}
%apply (uint64_t* INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4 ) {(uint64_t* data, size_t r1, size_t r2, size_t r3, size_t r4)}
%apply (int64_t* INPLACE_ARRAY1, int DIM1 ) {(int64_t* data, size_t r1)}
%apply (int64_t* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(int64_t* data, size_t r1, size_t r2)}
%apply (int64_t* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(int64_t* data, size_t r1, size_t r2, size_t r3)}
%apply (int64_t* INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4 ) {(int64_t* data, size_t r1, size_t r2, size_t r3, size_t r4)}


%include "pysz.h"
%include "defines.h"
%ignore ExaFELConfigBuilder::peaks;
%ignore ExaFELConfigBuilder::calibPanel;
%newobject ExaFELConfigBuilder::build;

%extend Compressor {
  %template(CompressFloat1) Compress1<float>;
  %template(CompressFloat2) Compress2<float>;
  %template(CompressFloat3) Compress3<float>;
  %template(CompressFloat4) Compress4<float>;
  %template(CompressDouble1) Compress1<double>;
  %template(CompressDouble2) Compress2<double>;
  %template(CompressDouble3) Compress3<double>;
  %template(CompressDouble4) Compress4<double>;
  %template(CompressUInt8_1) Compress1<uint8_t>;
  %template(CompressUInt8_2) Compress2<uint8_t>;
  %template(CompressUInt8_3) Compress3<uint8_t>;
  %template(CompressUInt8_4) Compress4<uint8_t>;
  %template(CompressInt8_1) Compress1<int8_t>;
  %template(CompressInt8_2) Compress2<int8_t>;
  %template(CompressInt8_3) Compress3<int8_t>;
  %template(CompressInt8_4) Compress4<int8_t>;
  %template(CompressUInt16_1) Compress1<uint16_t>;
  %template(CompressUInt16_2) Compress2<uint16_t>;
  %template(CompressUInt16_3) Compress3<uint16_t>;
  %template(CompressUInt16_4) Compress4<uint16_t>;
  %template(CompressInt16_1) Compress1<int16_t>;
  %template(CompressInt16_2) Compress2<int16_t>;
  %template(CompressInt16_3) Compress3<int16_t>;
  %template(CompressInt16_4) Compress4<int16_t>;
  %template(CompressUInt32_1) Compress1<uint32_t>;
  %template(CompressUInt32_2) Compress2<uint32_t>;
  %template(CompressUInt32_3) Compress3<uint32_t>;
  %template(CompressUInt32_4) Compress4<uint32_t>;
  %template(CompressInt32_1) Compress1<int32_t>;
  %template(CompressInt32_2) Compress2<int32_t>;
  %template(CompressInt32_3) Compress3<int32_t>;
  %template(CompressInt32_4) Compress4<int32_t>;
  %template(CompressUInt64_1) Compress1<uint64_t>;
  %template(CompressUInt64_2) Compress2<uint64_t>;
  %template(CompressUInt64_3) Compress3<uint64_t>;
  %template(CompressUInt64_4) Compress4<uint64_t>;
  %template(CompressInt64_1) Compress1<int64_t>;
  %template(CompressInt64_2) Compress2<int64_t>;
  %template(CompressInt64_3) Compress3<int64_t>;
  %template(CompressInt64_4) Compress4<int64_t>;






  %template(DecompressFloat) Decompress<float>;
  %template(DecompressDouble) Decompress<double>;
  %template(DecompressUInt8) Decompress<uint8_t>;
  %template(DecompressInt8) Decompress<int8_t>;
  %template(DecompressUInt16) Decompress<uint16_t>;
  %template(DecompressInt16) Decompress<int16_t>;
  %template(DecompressUInt32) Decompress<uint32_t>;
  %template(DecompressInt32) Decompress<int32_t>;
  %template(DecompressUInt64) Decompress<uint64_t>;
  %template(DecompressInt64) Decompress<int64_t>;

  %pythoncode %{
    import numpy

    __Compress = {
      (1, numpy.dtype('float64')): CompressDouble1,
      (2, numpy.dtype('float64')): CompressDouble2,
      (3, numpy.dtype('float64')): CompressDouble3,
      (4, numpy.dtype('float64')): CompressDouble4,
      (1, numpy.dtype('float32')): CompressFloat1,
      (2, numpy.dtype('float32')): CompressFloat2,
      (3, numpy.dtype('float32')): CompressFloat3,
      (4, numpy.dtype('float32')): CompressFloat4,
      (1, numpy.dtype('uint8')): CompressUInt8_1,
      (2, numpy.dtype('uint8')): CompressUInt8_2,
      (3, numpy.dtype('uint8')): CompressUInt8_3,
      (4, numpy.dtype('uint8')): CompressUInt8_4,
      (1, numpy.dtype('int8')): CompressInt8_1,
      (2, numpy.dtype('int8')): CompressInt8_2,
      (3, numpy.dtype('int8')): CompressInt8_3,
      (4, numpy.dtype('int8')): CompressInt8_4,
      (1, numpy.dtype('uint16')): CompressUInt16_1,
      (2, numpy.dtype('uint16')): CompressUInt16_2,
      (3, numpy.dtype('uint16')): CompressUInt16_3,
      (4, numpy.dtype('uint16')): CompressUInt16_4,
      (1, numpy.dtype('int16')): CompressInt16_1,
      (2, numpy.dtype('int16')): CompressInt16_2,
      (3, numpy.dtype('int16')): CompressInt16_3,
      (4, numpy.dtype('int16')): CompressInt16_4,
      (1, numpy.dtype('uint32')): CompressUInt32_1,
      (2, numpy.dtype('uint32')): CompressUInt32_2,
      (3, numpy.dtype('uint32')): CompressUInt32_3,
      (4, numpy.dtype('uint32')): CompressUInt32_4,
      (1, numpy.dtype('int32')): CompressInt32_1,
      (2, numpy.dtype('int32')): CompressInt32_2,
      (3, numpy.dtype('int32')): CompressInt32_3,
      (4, numpy.dtype('int32')): CompressInt32_4,
      (1, numpy.dtype('uint64')): CompressUInt64_1,
      (2, numpy.dtype('uint64')): CompressUInt64_2,
      (3, numpy.dtype('uint64')): CompressUInt64_3,
      (4, numpy.dtype('uint64')): CompressUInt64_4,
      (1, numpy.dtype('int64')): CompressInt64_1,
      (2, numpy.dtype('int64')): CompressInt64_2,
      (3, numpy.dtype('int64')): CompressInt64_3,
      (4, numpy.dtype('int64')): CompressInt64_4,
    }

    __Decompress = {
      numpy.dtype('float32'): DecompressFloat,
      numpy.float32: DecompressFloat,

      numpy.dtype('float64'): DecompressDouble,
      numpy.float64: DecompressDouble,

      numpy.dtype('uint8'): DecompressUInt8,
      numpy.uint8: DecompressUInt8,
      numpy.dtype('int8'): DecompressInt8,
      numpy.int8: DecompressInt8,

      numpy.dtype('uint16'): DecompressUInt16,
      numpy.uint16: DecompressUInt16,
      numpy.dtype('int16'): DecompressInt16,
      numpy.int16: DecompressInt16,

      numpy.dtype('uint32'): DecompressUInt32,
      numpy.uint32: DecompressUInt32,
      numpy.dtype('int32'): DecompressInt32,
      numpy.int32: DecompressInt32,

      numpy.dtype('uint64'): DecompressUInt64,
      numpy.uint64: DecompressUInt64,
      numpy.dtype('int64'): DecompressInt64,
      numpy.int64: DecompressInt64,
    }

    def Compress(self, array, userparams=None):
      length = len(array.shape)
      dtype = array.dtype
      return self.__Compress[length, dtype](self, array, userparams)

    def Decompress(self, bytes, dims, dtype, userparams=None):
      try:
        values = self.__Decompress[dtype](self, bytes, list(dims), userparams)
        return self.numpy.reshape(values, dims)
      except KeyError as e:
        raise TypeError("type {} not supported".format(i.args[0]))


      

  %}

}

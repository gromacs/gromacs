/* ************************************************************************
* Copyright 2016 Advanced Micro Devices, Inc.
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
#if !defined( AMD_CLFFT_GENERATOR_TRANSPOSE_HEADER )
#define AMD_CLFFT_GENERATOR_TRANSPOSE_HEADER
#include <iomanip>
#include "private.h"
#include "repo.h"
#include "plan.h"
#include "generator.stockham.h"
#include "action.h"

#define AVAIL_MEM_SIZE 32768 

inline std::stringstream& clKernWrite(std::stringstream& rhs, const size_t tabIndex)
{
	rhs << std::setw(tabIndex) << "";
	return rhs;
}

namespace clfft_transpose_generator
{
//generate transepose kernel with sqaure 2d matrix of row major with arbitrary batch size
/*
Below is a matrix(row major) containing three sqaure sub matrix along column
The transpose will be done within each sub matrix.
[M0
M1
M2]
*/
clfftStatus genTransposeKernelBatched(const FFTGeneratedTransposeSquareAction::Signature & params, std::string& strKernel, const size_t& lwSize, const size_t reShapeFactor);

//generate transpose kernel with square 2d matrix of row major with blocks along the leading dimension
//aka leading dimension batched
/*
Below is a matrix(row major) contaning three square sub matrix along row
[M0 M2 M2]
*/
clfftStatus genTransposeKernelLeadingDimensionBatched(const FFTGeneratedTransposeNonSquareAction::Signature & params, std::string& strKernel, const size_t& lwSize, const size_t reShapeFactor);

//swap lines. This kind of kernels are using with combination of square transpose kernels to perform nonsqaure transpose 1:2 ratio
clfftStatus genSwapKernel(const FFTGeneratedTransposeNonSquareAction::Signature & params, std::string& strKernel, std::string& KernelFuncName, const size_t& lwSize, const size_t reShapeFactor);

clfftStatus genSwapKernelGeneral(const FFTGeneratedTransposeNonSquareAction::Signature & params, std::string& strKernel, std::string& KernelFuncName, const size_t& lwSize, const size_t reShapeFactor);

void get_cycles(size_t *cycle_map, size_t num_reduced_row, size_t num_reduced_col);

void permutation_calculation(size_t m, size_t n, std::vector<std::vector<size_t> > &permutationVec);
}//end of namespace clfft_transpose_generator

#endif
#include <algorithm>

#include <vector>
#include <iterator>
#include "sz.h"
#include "zlib.h"

namespace {
	const int r = 64;
	const int size = r*r*r;
}

sz_params
sz_default_config()
{
  sz_params params;
  params.dataType = SZ_FLOAT;
  params.max_quant_intervals = 65536;
  params.quantization_intervals = 0;
  params.predThreshold = 0.99;
  params.szMode = SZ_BEST_SPEED;
  params.gzipMode = Z_BEST_SPEED;
  params.errorBoundMode = ABS;
  params.absErrBound = 1e-4;
  params.relBoundRatio = 1e-4;
  params.psnr = 80;
  params.pw_relBoundRatio = 1e-2;
  params.segment_size = 25;
  params.pwr_type = SZ_PWR_MIN_TYPE;
  params.sampleDistance = 100;
  params.sol_ID = SZ;
  params.predictionMode = SZ_PREVIOUS_VALUE_ESTIMATE;
	params.randomAccess = true;
	params.losslessCompressor = GZIP_COMPRESSOR;
  return params;
}

template <class OutputIt>
unsigned int count_non_equal(unsigned char* it1, unsigned char* it2, OutputIt out, size_t size)
{
	unsigned int not_equal_elms = 0; 
	for (size_t i = 0; i < size; ++i) {
		bool is_equal;
		if((is_equal = it1[i] != it2[i])) {
			not_equal_elms++;
		}
		*out = is_equal;
		++out;
	}
	return not_equal_elms;
}

template <class Function>
std::vector<decltype(std::declval<Function>()(std::declval<int>()))>
evaluate(int size, Function func)
{
	int i = 0;
	std::vector<decltype(func(i))> values (size);
	std::generate(std::begin(values), std::end(values), [&i, &func]() { return func(i++); });
	return values;
}

template <class Func, class Func2>
void run_test(std::vector<float> &test_array, Func sz_compress1, Func2 sz_compress2) {
	std::vector<float> test_array2(test_array);
	std::vector<bool> not_equal_entries;

	size_t size1 = 0;
	sz_params conf_params1{sz_default_config()};
	SZ_Init_Params(&conf_params1);
	unsigned char* data1 = sz_compress1(test_array.data(), r, r, r, .99, &size1);
	SZ_Finalize();

	size_t size2 = 0;
	sz_params conf_params2{sz_default_config()};
	SZ_Init_Params(&conf_params2);
	unsigned char* data2 = sz_compress2(test_array2.data(), r, r, r, .99, &size2);
	SZ_Finalize();
	

	//testing using:
	// CU_ASSERT_EQUAL_ARRAY_BYTE(data1, data2, std::min(size1, size2));
	//outputs too many failures to be useful prefer the following instead
	
	CU_ASSERT_EQUAL(size1, size2);
	auto num_not_equal = count_non_equal(data1, data2, std::back_inserter(not_equal_entries),  std::min(size1, size2));
	CU_ASSERT_EQUAL(num_not_equal, 0);

	free(data1);
	free(data2);
}

template <class Compressor, class Compressor2>
void test_identical_output_random(int num_random_test_cases, Compressor func, Compressor2 func2)
{
	std::seed_seq seed{0};
	std::mt19937 gen{seed};
	std::uniform_real_distribution<float> dist(0.0f,1.0f);
	auto rand = [&gen,&dist](){return dist(gen);};

	for (int test_case = 0; test_case < num_random_test_cases; ++test_case) {
		std::vector<float> test_array(size, 0.0f);
		std::generate(std::begin(test_array), std::end(test_array), rand);
		run_test(test_array, func, func2);
	}

}

template <class Compressor, class Compressor2>
void test_identical_output_deterministic(Compressor func, Compressor2 func2)
{
	std::vector<std::vector<float>> deterministic_test_cases{
		evaluate(size, [](int ) { return 0.0f; }),
		evaluate(size, [](int i){ return 3.14f * i - 2.1f; }),
		evaluate(size, [](int i){ return 7.2f * (i*i) + 3.14f*i - 2.1f; }),
		evaluate(size, [](int i){ return 1.1f * (i*i*i) + 7.2f*(i*i) + 3.14f*i - 2.1f;}),
		evaluate(size, [](int i){ return -2.1f * (i*i*i*i) + 1.1f*(i*i*i) + 7.2f*(i*i) + 3.14f*i - 2.1f;}),
	};

	for(auto& test_case: deterministic_test_cases)
	{
		run_test(test_case, func, func2);
	}

}

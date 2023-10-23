#include <vector>
#include <random>
#include <algorithm>
#include <functional>

#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include "CUnit_Array.h"
#include "RegressionTest.hpp"

extern "C" {
int init_suite()
{
	return 0;
}

int clean_suite()
{
	return 0;
}

void test_all_functions()
{
	auto functions = {
		SZ_compress_float_3D_MDQ_decompression_random_access_with_blocked_regression,
	};
	auto num_random_test_cases = 4;

	for (auto& func : functions) {
		test_identical_output_random(num_random_test_cases, func, func);
		test_identical_output_deterministic(func, func);
	}
}

int main(int argc, char *argv[])
{
	unsigned int num_failures = 0;
	if (CUE_SUCCESS != CU_initialize_registry())
	{
		return CU_get_error();
	}

	CU_pSuite suite = CU_add_suite("test_opencl_suite", init_suite, clean_suite);
	if(suite == nullptr) {
		goto error;
	}

	if(CU_add_test(suite, "test_all_functions", test_all_functions) == nullptr) {
		goto error;
	}

	CU_basic_set_mode(CU_BRM_VERBOSE);
	CU_basic_run_tests();
	CU_basic_show_failures(CU_get_failure_list());
	num_failures = CU_get_number_of_failures();

error:
	CU_cleanup_registry();
	return  num_failures ||  CU_get_error();
}

}

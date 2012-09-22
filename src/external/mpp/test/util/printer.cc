#include "printer.h"

// Since most methods are very similar, use macros to reduce boilerplate.
// This defines a member that forwards the call only on Master
#define GTEST_FORWARD_ONROOT_METHOD_(Name, Type)    \
void MPIEventForward::Name(const Type& parameter) { \
  if (comm_.rank()==0) {              \
    listener_->Name(parameter); \
  } \
}

//copied declarations
namespace testing {
namespace internal {
enum GTestColor {
    COLOR_DEFAULT,
    COLOR_RED,
    COLOR_GREEN,
    COLOR_YELLOW
};

void ColoredPrintf(GTestColor color, const char* fmt, ...);
void PrintTestName(const char * test_case, const char * test) {
    printf("%s.%s", test_case, test);
}
void PrintFullTestCommentIfPresent(const TestInfo& test_info);
}
}

using namespace testing;
using namespace testing::internal;

GTEST_FORWARD_ONROOT_METHOD_(OnTestProgramStart, UnitTest)
GTEST_FORWARD_ONROOT_METHOD_(OnEnvironmentsSetUpStart, UnitTest)
GTEST_FORWARD_ONROOT_METHOD_(OnTestStart, TestInfo)
GTEST_FORWARD_ONROOT_METHOD_(OnEnvironmentsTearDownStart, UnitTest)
GTEST_FORWARD_ONROOT_METHOD_(OnEnvironmentsSetUpEnd, UnitTest)
GTEST_FORWARD_ONROOT_METHOD_(OnEnvironmentsTearDownEnd, UnitTest)
GTEST_FORWARD_ONROOT_METHOD_(OnTestCaseEnd, TestCase)
GTEST_FORWARD_ONROOT_METHOD_(OnTestProgramEnd, UnitTest)

//copied static functions

// Formats a countable noun.  Depending on its quantity, either the
// singular form or the plural form is used. e.g.
//
// FormatCountableNoun(1, "formula", "formuli") returns "1 formula".
// FormatCountableNoun(5, "book", "books") returns "5 books".
static internal::String FormatCountableNoun(int count,
                                            const char * singular_form,
                                            const char * plural_form) {
  return internal::String::Format("%d %s", count,
                                  count == 1 ? singular_form : plural_form);
}

// Formats the count of tests.
static internal::String FormatTestCount(int test_count) {
  return FormatCountableNoun(test_count, "test", "tests");
}

// Formats the count of test cases.
static internal::String FormatTestCaseCount(int test_case_count) {
  return FormatCountableNoun(test_case_count, "test case", "test cases");
}

MPIEventForward* MPIEventForward::instance = NULL;

void MPIEventForward::ComputeFailedTests(const UnitTest& unit_test, std::vector<int> &results) {
  for (int i = 0; i < unit_test.total_test_case_count(); ++i) {
    const TestCase& test_case = *unit_test.GetTestCase(i);
    if (!test_case.should_run()) {
        continue;
    }
    for (int j = 0; j < test_case.total_test_count(); ++j) {
      const TestInfo& test_info = *test_case.GetTestInfo(j);
      if (!test_info.should_run()) {
          continue;
      }
      results.push_back(!test_info.result()->Passed());
    }
  }
  comm_.reduce(results.begin(), results.end(),
                          MPI_BOR, 0);
}

static void PrintFailedTests(const UnitTest& unit_test, std::vector<int> &results) {
  std::vector<int>::iterator result_it = results.begin();
  //has to loop exactly like ComputePassedTests to match results vector
  for (int i = 0; i < unit_test.total_test_case_count(); ++i) {
    const TestCase& test_case = *unit_test.GetTestCase(i);
    if (!test_case.should_run()) {
        continue;
    }
    for (int j = 0; j < test_case.total_test_count(); ++j) {
      const TestInfo& test_info = *test_case.GetTestInfo(j);
      if (!test_info.should_run() || !*(result_it++)) {
          continue;
      }
      ColoredPrintf(COLOR_RED, "[  FAILED  ] ");
      printf("%s.%s", test_case.name(), test_info.name());
      PrintFullTestCommentIfPresent(test_info);
      printf("\n");
    }
  }
}

void MPIEventForward::OnTestCaseStart(const TestCase& test_case) {
    test_case_name_ = test_case.name();
    if (comm_.rank()==0) {
        listener_->OnTestCaseStart(test_case);
    }
}

void MPIEventForward::OnTestIterationStart(const UnitTest& unit_test,
                                             int iteration) {
  if (comm_.rank()==0) {
    listener_->OnTestIterationStart(unit_test, iteration);
  }
}

void MPIEventForward::OnTestIterationEnd(const UnitTest& unit_test,
                                                     int /*iteration*/) {
  std::vector<int> failed_tests;
  ComputeFailedTests(unit_test, failed_tests);
  int num_failures = std::accumulate(failed_tests.begin(),failed_tests.end(),0);
  if (comm_.rank()!=0) return;
  ColoredPrintf(COLOR_GREEN,  "[==========] ");
  printf("%s from %s ran.",
         FormatTestCount(unit_test.test_to_run_count()).c_str(),
         FormatTestCaseCount(unit_test.test_case_to_run_count()).c_str());
  if (GTEST_FLAG(print_time)) {
    printf(" (%s ms total)",
           internal::StreamableToString(unit_test.elapsed_time()).c_str());
  }
  printf("\n");
  ColoredPrintf(COLOR_GREEN,  "[  PASSED  ] ");
  printf("%s.\n", FormatTestCount(unit_test.test_to_run_count()-num_failures).c_str());

  if (num_failures!=0) {
    ColoredPrintf(COLOR_RED,  "[  FAILED  ] ");
    printf("%s, listed below:\n", FormatTestCount(num_failures).c_str());
    PrintFailedTests(unit_test,failed_tests);
    printf("\n%2d FAILED %s\n", num_failures,
                        num_failures == 1 ? "TEST" : "TESTS");
  }

  int num_disabled = unit_test.disabled_test_count();
  if (num_disabled && !GTEST_FLAG(also_run_disabled_tests)) {
    if (!num_failures) {
      printf("\n");  // Add a spacer if no FAILURE banner is displayed.
    }
    ColoredPrintf(COLOR_YELLOW,
                  "  YOU HAVE %d DISABLED %s\n\n",
                  num_disabled,
                  num_disabled == 1 ? "TEST" : "TESTS");
  }
  // Ensure that Google Test output is printed before, e.g., heapchecker output.
  fflush(stdout);
}

void MPIEventForward::OnTestPartResult(const TestPartResult& result) {
    //nothing - reporting in OnTestEnd
}

void MPIEventForward::OnTestEnd(const TestInfo& test_info) {
    int all_passed = true;
    //use computation of all_passed to make the printing sequentially
    int rank = comm_.rank(), size = comm_.size();
    if (rank>0 && size>1) {
        comm_(rank-1) >> all_passed;
    }
    const TestResult* result = test_info.result();
    if (result->Failed()) {
        printf("rank %d:\n", rank);
        for (int i=0;i<result->total_part_count();i++) {
            listener_->OnTestPartResult(result->GetTestPartResult(i));
        }
        all_passed = false;
    }
    if (size>1) {
        comm_((rank+1)%comm_.size()) << all_passed;
        if (rank==0) {
           comm_(comm_.size()-1) >> all_passed;
        }
    }

    if (comm_.rank()!=0) return;
    if (all_passed) {
        ColoredPrintf(COLOR_GREEN, "[       OK ] ");
    } else {
        ColoredPrintf(COLOR_RED, "[  FAILED  ] ");
    }
    PrintTestName(test_case_name_.c_str(), test_info.name());
    if (!all_passed)
        PrintFullTestCommentIfPresent(test_info);

    //Should we print max/avg?
    if (GTEST_FLAG(print_time)) {
        printf(" (%s ms)\n", internal::StreamableToString(
                             test_info.result()->elapsed_time()).c_str());
    } else {
        printf("\n");
    }
    fflush(stdout);
}

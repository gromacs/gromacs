#include "gtest/gtest.h"
#include "mpp.h"

using ::testing::UnitTest;

//Forwards as much as possible to the default Pretty Printer
//Reimplements anything which needs to be changed for MPI
//Code mostly copied from PrettyUnitTestResultPrinter
//Tried to copy as little as possible but the pretty printer
//isn't very extensible.
class MPIEventForward : public ::testing::TestEventListener {
public:
  explicit MPIEventForward(TestEventListener* listener) : listener_ (listener),
      comm_(mpi::comm::world) {}
  ~MPIEventForward() { if (listener_) delete listener_; }

  virtual void OnTestProgramStart(const UnitTest& unit_test);
  virtual void OnTestIterationStart(const UnitTest& unit_test, int iteration);
  virtual void OnEnvironmentsSetUpStart(const UnitTest& unit_test);
  virtual void OnEnvironmentsSetUpEnd(const UnitTest& unit_test);
  virtual void OnTestCaseStart(const ::testing::TestCase& test_case);
  virtual void OnTestStart(const ::testing::TestInfo& test_info);
  virtual void OnTestPartResult(const ::testing::TestPartResult& result);
  virtual void OnTestEnd(const ::testing::TestInfo& test_info);
  virtual void OnTestCaseEnd(const ::testing::TestCase& test_case);
  virtual void OnEnvironmentsTearDownStart(const UnitTest& unit_test);
  virtual void OnEnvironmentsTearDownEnd(const UnitTest& unit_test);
  virtual void OnTestIterationEnd(const UnitTest& unit_test, int iteration);
  virtual void OnTestProgramEnd(const UnitTest& unit_test);

  static void Insert() {
    UnitTest& unit_test = *UnitTest::GetInstance();
    ::testing::TestEventListeners& listeners = unit_test.listeners();
    ::testing::TestEventListener* defprinter = 
          listeners.Release(listeners.default_result_printer());
    assert(instance==NULL);
    instance = new MPIEventForward(defprinter);
    listeners.Append(instance);
  }
  static void Remove() {
    UnitTest& unit_test = *UnitTest::GetInstance();
    ::testing::TestEventListeners& listeners = unit_test.listeners();
    listeners.Release(instance);
    listeners.Append(instance->listener_);
    instance->listener_ = NULL;
    delete instance;
  }

private:
  TestEventListener* listener_;
  GTEST_DISALLOW_COPY_AND_ASSIGN_(MPIEventForward);
  testing::internal::String test_case_name_;
  mpi::comm comm_;
  static MPIEventForward* instance;

  void ComputeFailedTests(const UnitTest&, std::vector<int> &);

};  // class TersePrinter

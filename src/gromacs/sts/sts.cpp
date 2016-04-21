#include "sts.h"

std::deque<Thread> STS::threads_ = {};
std::atomic<int> STS::stepCounter_(0);
STS* STS::defaultInstance_ = nullptr;
std::map<std::string, STS *> STS::stsInstances_ = {};
STS* STS::instance_ = nullptr;

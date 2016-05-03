#include "sts.h"

#include <memory>

std::unique_ptr<STS> STS::instance_(new STS());

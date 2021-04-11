#ifndef __INTEL_LLVM_COMPILER
#error
#endif
#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)
#pragma message(VALUE(__INTEL_LLVM_COMPILER))
int main() {}

//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#ifndef RCSBPLATFORM_H
#define RCSBPLATFORM_H


#include "rcsb_types.h"


class RcsbPlatform
{
  public:
    RcsbPlatform();
    ~RcsbPlatform();

    static bool IsLittleEndian();

  private:
    static const UInt16 _ENDIANNESS_TEST_INT;

};


#endif // RCSBPLATFORM_H not defined


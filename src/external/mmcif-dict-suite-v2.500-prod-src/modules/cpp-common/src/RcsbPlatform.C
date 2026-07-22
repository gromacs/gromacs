//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#include "RcsbPlatform.h"


const UInt16 RcsbPlatform::_ENDIANNESS_TEST_INT = 0x0001;


RcsbPlatform::RcsbPlatform()
{

}


RcsbPlatform::~RcsbPlatform()
{

}


bool RcsbPlatform::IsLittleEndian()
{
    char* leastSigByte = (char*) &_ENDIANNESS_TEST_INT;

    return (leastSigByte[0] ? true : false);
}


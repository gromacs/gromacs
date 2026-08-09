//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$

#ifndef RCSB_TYPES_H
#define RCSB_TYPES_H

// Basic types
typedef char SInt8;
typedef unsigned char UInt8;
typedef short int SInt16;
typedef unsigned short int UInt16;
typedef int SInt32;
typedef unsigned int UInt32;

// Basic types sizes in octets
const UInt8 UINT32_SIZE = sizeof(UInt32);

enum eTypeCode
{
    eTYPE_CODE_NONE = 0,
    eTYPE_CODE_INT,        // 1
    eTYPE_CODE_FLOAT,      // 2
    eTYPE_CODE_STRING,     // 3
    eTYPE_CODE_TEXT,       // 4
    eTYPE_CODE_DATETIME,   // 5
    eTYPE_CODE_BIGINT,     // 6
    eTYPE_CODE_FLOAT_SCI   // 7 - float representing scientific notation - used for XML types
};

#endif // RCSB_TYPES_H not defined


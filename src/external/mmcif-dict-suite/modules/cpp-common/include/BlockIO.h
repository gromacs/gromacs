//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#ifndef BLOCKIO_H
#define BLOCKIO_H

#include "rcsb_types.h"

const unsigned int WORDSPERBLOCK = 2048;
const unsigned int BLKSIZE = 8192;

class BlockIO
{
public:
    BlockIO();
    ~BlockIO();

    void AssociateBuffer(char** newBuffer);

    unsigned int ReadBlock(const int fd, const UInt32 blockNum);
    unsigned int WriteBlock(const int fd, const UInt32 blockNum);

private:
    UInt32 _buffer[WORDSPERBLOCK]; // A buffer for reading/writing blocks

};

#endif


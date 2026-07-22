//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*
 ** Class for blockwise IO operations
 */

#ifdef VLAD_DECIDE
#include <strings.h>
#endif

#include <string.h>
#include <unistd.h>
#include <fcntl.h>

#include <string>

#include "rcsb_types.h"
#include "Exceptions.h"
#include "GenString.h"
#include "BlockIO.h"

using std::string;

BlockIO::BlockIO()
{
    memset(_buffer, 0, BLKSIZE);
}

BlockIO::~BlockIO()
{

}

unsigned int BlockIO::ReadBlock(const int fd, const UInt32 blockNum)
{
    if (lseek(fd, BLKSIZE * blockNum, SEEK_SET) == -1)
    {
        throw FileException(string("Could not seek in the file with fd: ") + \
          String::IntToString(fd), "BlockIO::ReadBlock");
    }

    // VLAD - CHECK FOR -1 AND THROW EXCEPTION
    return read(fd, _buffer, BLKSIZE);
}

unsigned int BlockIO::WriteBlock(const int fd, const UInt32 blockNum)
{
    if (lseek(fd, blockNum * BLKSIZE, SEEK_SET) == -1)
    {
        throw FileException(string("Could not seek in the file with fd: ") + \
          String::IntToString(fd), "BlockIO::WriteBlock");
    }

    // VLAD - CHECK FOR -1 AND THROW EXCEPTION
    return write(fd, _buffer, BLKSIZE);
}

void BlockIO::AssociateBuffer(char** newBuffer)
{
    *newBuffer = (char*)_buffer;
}


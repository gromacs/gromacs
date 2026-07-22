//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#ifndef SERIALIZER_H
#define SERIALIZER_H


#include <string>
#include <vector>
#include <fstream>

#include "rcsb_types.h"
#include "BlockIO.h"


const int NO_TYPE = 0; // This is reserved
const unsigned int STRINGS_TYPE = 1;
const unsigned int STRING_TYPE = 2;
const int INT_TYPE = 3;
const int LONG_TYPE = 4;
const int FLOAT_TYPE = 5;
const int DOUBLE_TYPE = 6;
const unsigned int WORD_TYPE = 7;
const unsigned int WORDS_TYPE = 8;
const unsigned int UWORD_TYPE = 9;
const unsigned int UWORDS_TYPE = 10;

const int INDEX_INCREMENT = 1024;

enum eFileMode
{
    NO_MODE = 0,
    READ_MODE,
    CREATE_MODE,
    UPDATE_MODE,
    VIRTUAL_MODE
};


class Serializer
{
  public:
    // Constructors and destructor
    Serializer(const std::string& fileName, const eFileMode fileMode,
      const bool verbose = false);
    ~Serializer();

    inline unsigned int GetNumDataIndices();

    // Read methods
    UInt32 ReadUInt32(const UInt32 index);
    void ReadUInt32s(std::vector<UInt32>& UInt32s, const UInt32 index);
    void ReadString(std::string& retString, const UInt32 index);
    void ReadStrings(std::vector<std::string>& theStrings, const UInt32 index);

    // Write methods 
    UInt32 WriteUInt32(const UInt32 theWord);
    UInt32 WriteUInt32s(const std::vector<UInt32>& theWords);
    UInt32 WriteString(const std::string& theString);
    UInt32 WriteStrings(const std::vector<std::string>& theStrings);

    // Update methods
    UInt32 UpdateUInt32(const UInt32 theWord, const UInt32 oldIndex);
    UInt32 UpdateUInt32s(const std::vector<UInt32>& theWords,
      const UInt32 oldIndex);
    UInt32 UpdateString(const std::string& theString, const UInt32 oldIndex);
    UInt32 UpdateStrings(const std::vector<std::string>& theStrings,
      const UInt32 oldIndex);

  private:
    typedef struct
    {
        // Block number of the start of the indices info
        UInt32 fileIndexBlock;

        // Number of blocks that hold the indices
        UInt32 fileIndexNumBlocks;

        // Total size in bytes of all indices
        UInt32 fileIndexLength;

        // Number of indices
        UInt32 numIndices;

        // Reserved information
        UInt32 reserved[3];

        // File version
        UInt32 version;
    } tFileHeader;

    // Represents an entry index. Entry is a value of some type that has
    // been stored in to the file. Index shows entry location (in which
    // block and offset from the start of the block), its length, dataType.
    typedef struct
    {
        UInt32 blockNumber; // Block in which data is located
        UInt32 offset;      // Data offset in the block
        UInt32 length;      // The length of the data 
        UInt32 dataType;    // Type of data
        UInt32 vLength;     // Virtual length (length adjusted for word size)
        UInt32 reserved[3];
    } EntryIndex;

    static bool _littleEndian;

    static const UInt32 _version = 1;
    static const UInt32 _indicesPerBlock = BLKSIZE / sizeof(EntryIndex);

    // An array of index entries (i.e., these are indices)
    std::string _fileName;

    // Stored in block 0, this holds the info about the index
    tFileHeader _fileHeader;

    std::vector<EntryIndex> _indices;

    std::ofstream _log;

    bool _verbose;

    UInt32 _currentBlock;  // The current block number of the current buffer
    UInt32 _currentOffset; // The offset into the current buffer

    char* _buffer;

    eFileMode _mode;

    BlockIO _theBlock; // A block for doing read/write a block at a time

    void Init();

    void WriteUInt32AtIndex(const UInt32 theWord, const UInt32 index);
    void WriteUInt32sAtIndex(const std::vector<UInt32>& Words,
      const UInt32 index);
    void WriteStringAtIndex(const std::string& theString, const UInt32 index);
    void WriteStringsAtIndex(const std::vector<std::string>& theStrings,
      const UInt32 index);

    void Delete(const UInt32 index);

    void GetLastDataBuffer(void);
    void GetDataBufferAtIndex(const UInt32 index);

    void _GetHeader(const char* where);
    void _PutHeader(char* where);

    void _GetIndex(EntryIndex& outIndex, const char* where);
    void _PutIndex(const EntryIndex& inIndex, char* where);

    UInt32 _GetUInt32(const char* where);
    void _PutUInt32(const UInt32 inWord, char* where);

    void SwapHeader(tFileHeader& out, const tFileHeader& in);
    void SwapIndex(EntryIndex& out, const EntryIndex& in);
    UInt32 SwapUInt32(const UInt32 theWord);

    void _ReadFileHeader();
    void _WriteFileHeader();

    void AllocateIndices(const UInt32 index);

    UInt32 ReadBlock(const UInt32 blockNum);
    UInt32 WriteBlock(const UInt32 blockNum);

    char* GetWritingPoint(const UInt32 index);
    void WriteLast(const char* const where);

    void SetVirtualLength(const UInt32 index);

    int _fd; // The file descriptor of the file that is opened, -1 if unopened

    UInt32 _numBlocksIO; // The number of blocks in the currently opened file
    UInt32 _currentBlockIO; // The block that is currently read into buffer

    void OpenFileIO(const std::string& filename, const eFileMode fileMode);
    void CloseFileIO();

    inline UInt32 GetCurrentBlockNumberIO() const;
    inline UInt32 GetNumBlocksIO() const;

    void PrintIndex();
    void PrintIndexPosition(const UInt32 position);
    void DumpFile();
    void PrintBuffer();
};


inline UInt32 Serializer::GetNumDataIndices()
{
    return (_indices.size());
}


inline UInt32 Serializer::GetCurrentBlockNumberIO() const
{
    return (_currentBlockIO);
}


inline UInt32 Serializer::GetNumBlocksIO() const
{
    return (_numBlocksIO);
}

#endif


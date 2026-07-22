//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#include <stdio.h>
#include <string>
#include <ios>
#include <fstream>

#include "Exceptions.h"
#include "RcsbFile.h"


using std::string;
using std::ios;
using std::ofstream;


const string RcsbFile::DIR_SEPARATOR = "/";


RcsbFile::RcsbFile()
{

}


RcsbFile::~RcsbFile()
{

}


bool RcsbFile::IsEmpty(ofstream& fileStream)
{
    fileStream.seekp(0, ios::end);

    unsigned int fileSize = fileStream.tellp();

    if (fileSize != 0)
        return (false);

    return (true); 
}


void RcsbFile::Delete(const string& fileName)
{
    int ret = remove(fileName.c_str());

    if (ret != 0)
    {
        throw NotFoundException("File \"" + fileName + "\" not found",
          "RcsbFile::Delete");
    }
}


void RcsbFile::RelativeFileName(string& relName, const string& absName)
{
    relName.clear();

    if (absName.empty())
    {
        return;
    }

    // Find the directory separator doing a reverse search
    string::size_type dirSepIndex = absName.rfind(DIR_SEPARATOR,
      absName.size() - 1);

    if (dirSepIndex != string::npos)
    {
        // Found it.
        relName = absName.substr(dirSepIndex + 1);
    }
    else
    {
        relName = absName;
    }
}


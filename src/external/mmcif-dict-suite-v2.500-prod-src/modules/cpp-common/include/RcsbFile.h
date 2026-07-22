//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#ifndef RCSBFILE_H
#define RCSBFILE_H


#include <string>
#include <fstream>


class RcsbFile
{
  public:
    static const std::string DIR_SEPARATOR;

    static bool IsEmpty(std::ofstream& fileStream);
    static void Delete(const std::string& fileName);

    static void RelativeFileName(std::string& relName,
      const std::string& absName);

  private:
    RcsbFile();

    ~RcsbFile();
};


#endif // RCSBFILE_H not defined


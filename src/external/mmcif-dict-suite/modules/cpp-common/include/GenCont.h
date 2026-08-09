//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#ifndef GENCONT_H
#define GENCONT_H


#include <string>
#include <vector>

#include "GenString.h"


std::ostream& operator<<(std::ostream& out,
  const std::vector<std::string>& contVector);


class GenCont
{
  public:
    static bool IsInVector(const std::string& element,
      const std::vector<std::string>& contVector, const Char::eCompareType
      compareType = Char::eCASE_SENSITIVE);
    static bool IsInVectorCi(const std::string& element,
      const std::vector<std::string>& contVector);

  private:
    GenCont();

    ~GenCont();
};


#endif // GENCONT_H not defined


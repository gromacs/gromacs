#include <string>
#include <sstream>
#include <math.h>
#include "stringutil.hpp"

std::vector<std::string> &split(const std::string &s, 
                                char delim, 
                                std::vector<std::string> &elems) 
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s,
                               char delim) 
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

std::string gmx_ftoa(double f)
{
    char buf[32];
    
    if (fabs(f) < 100)
    {
        sprintf(buf, "%.3f", f);
    }
    else
    {
        sprintf(buf, "%g", f);
    }
    return std::string(buf);
}

std::string gmx_itoa(int f)
{
    char a[32];

    sprintf(a, "%d", f);

    return std::string(a);
}

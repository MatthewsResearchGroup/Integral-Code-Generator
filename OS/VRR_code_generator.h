#ifndef _VRR_CODE_GENERATOR_H
#define _VRR_CODE_GENERATOR_H

#include <array>
// #include <string>


std::string NewNameScheme(const std::array<int, 3>& ang_mom);

int dirchoose(const std::array<int, 3> a);

std::string namemap(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c,const std::array<int, 3> d, int m);

char center_decrease(std::array<int, 4> angmom);


#endif


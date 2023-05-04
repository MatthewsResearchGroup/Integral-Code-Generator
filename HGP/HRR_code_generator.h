#ifndef _HRR_CODE_GENERATOR_H
#define _HRR_CODE_GENERATOR_H

#include "name.h"
#include "VRR_code_generator.h"

std::string sent_gen(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c,const std::array<int, 3> d, int m, char center, int xyz);

void hgp(const std::array<int, 3>& a, const std::array<int, 3>& b,const std::array<int, 3>& c,const std::array<int, 3>& d, int m, std::map<std::array<int, 13>, std::string>& osmap);

#endif

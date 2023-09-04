#ifndef _VRR_CODE_GENERATOR_H
#define _VRR_CODE_GENERATOR_H

#include <array>
#include <map>

void vrr_code_print(int la , int lb, int lc, int ld, std::map<std::array<int, 13>, std::string>& osmap);

void eri(int la, int lb, int lc, int ld, std::map<std::array<int, 13>, std::string>& osmap);

void obara_saika(const std::array<int, 3>& a, const std::array<int, 3>& b,const std::array<int, 3>& c,const std::array<int, 3>& d, int m, std::map<std::array<int, 13>, std::string>& osmap);

#endif


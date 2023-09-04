#ifndef _VRR_3C_2E_H
#define _VRR_3C_2E_H

#include <array>
#include <map>
void vrr_code_print(int la , int lb, int lc, std::map<std::array<int, 10>, std::string>& osmap);
void eri(int la, int lb, int lc, std::map<std::array<int, 10>, std::string>& osmap);
void obara_saika(const std::array<int, 3>& a, const std::array<int, 3>& b,const std::array<int, 3>& c, int m, std::map<std::array<int, 10>, std::string>& osmap);

#endif

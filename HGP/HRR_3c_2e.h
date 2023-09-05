#ifndef _HRR_3C_2E_H
#define _HRR_3C_2E_H

#include <name.h>
#include <VRR_code_generator.h>
#include <VRR_3c_2e.h>
#include <array>
#include <string>
#include <iostream>
#include <map>


std::string sent_gen(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c, int m, int xyz);

void hgp(const std::array<int, 3>& a, const std::array<int, 3>& b,const std::array<int, 3>& c, int m, std::map<std::array<int, 10>, std::string>& osmap);

void hgp_eri(int la, int lb, int lc, std::map<std::array<int, 10>, std::string>& osmap);
#endif

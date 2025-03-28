#ifndef _NAME_H
#define _NAME_H
#include <array>
#include <string>
#include <iostream>
#include "name.h"
#include <cmath>
#include <iomanip>
#include <climits>

std::string NewNameScheme(const std::array<int, 3>& ang_mom);
int dirchoose(const std::array<int, 3> a);
// four centers
std::string namemap(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c,const std::array<int, 3> d, int m);
// three centers
std::string namemap(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c, int m);
// center choose of obara saika algorithm
char center_choose(std::array<int, 4> angmom);
// center choose of os algorithm for 3 centers
char center_choose(std::array<int, 3> angmom);
// cneter choose for HGP algorithm
char center_choose(int a, int b);
bool intex(const std::array<int,3>& a, const std::array<int,3>& b, const std::array<int,3>& c, const std::array<int,3>& d);
bool intex(const std::array<int,3>& a, const std::array<int,3>& b, const std::array<int,3>& c);
int addrsear(std::array<int, 3> a);
std::string orbname(const std::array<int, 3>& ang_mom);
void save_int(int la, int lb, int lc, int ld);
#endif

#ifndef _NAME_H
#define _NAME_H
#include <array>

std::string NewNameScheme(const std::array<int, 3>& ang_mom);

int dirchoose(const std::array<int, 3> a);

std::string namemap(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c,const std::array<int, 3> d, int m);

char center_choose(std::array<int, 4> angmom);

char center_choose(int a, int b);

bool intex(const std::array<int,3>& a, const std::array<int,3>& b, const std::array<int,3>& c, const std::array<int,3>& d);

int addrsear(std::array<int, 3> a);

std::string orbname(const std::array<int, 3>& ang_mom);

void save_int(int la, int lb, int lc, int ld);
#endif

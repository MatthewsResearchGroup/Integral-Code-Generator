#include <string>
#include <iostream>
#include "name.h"

std::string NewNameScheme(const std::array<int, 3>& ang_mom)
{
    // return orbital name, like for d(2,0,0), return d200
    std::string s;
    int tot_ang = ang_mom[0] + ang_mom[1] + ang_mom[2];
    if (tot_ang == 0)
        s.append("s");
    else if (tot_ang == 1)
        s.append("p");
    else if (tot_ang == 2)
        s.append("d");
    else if (tot_ang == 3)
        s.append("f");
    else if (tot_ang == 4)
        s.append("g");
    else if (tot_ang == 5)
        s.append("h");
    else if (tot_ang == 6)
        s.append("i");
    s.append(std::to_string(ang_mom[0]));
    s.append(std::to_string(ang_mom[1]));
    s.append(std::to_string(ang_mom[2]));
    s.append("_");
    return s;
}

int dirchoose(const std::array<int, 3> a)
{
    int xyz = 0;
    int minx = INT_MAX;
    for (auto i=0; i < 3; i++)
    {
        if ((a[i] != 0) and (a[i] < minx))
        {
            xyz = i;
            minx = a[i];
        }
    }
    return xyz;
}


std::string namemap(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c,const std::array<int, 3> d, int m)
{
    // auto stra = xyz_tran(a);
    // auto strb = xyz_tran(b);
    // auto strc = xyz_tran(c);
    // auto strd = xyz_tran(d);
    auto stra = NewNameScheme(a);
    auto strb = NewNameScheme(b);
    auto strc = NewNameScheme(c);
    auto strd = NewNameScheme(d);
    if (stra == "Error" || strb == "Error" || strc == "Error" || strd == "Error")
        // return "TERM_DOES_NOT_EXIST";
        return "0.f";
    else
    {   
        stra.append(strb);
        stra.append(strc);
        stra.append(strd);
        stra.append(std::to_string(m));
        return stra;
    }   
}


char center_choose(std::array<int, 4> angmom)
{// Array, the total momentum for four centers.
    auto num_centers = 4;
    auto total_ang = 0;
    auto min_ang = INT_MAX;
    auto center = -1;
    for (auto i = 0; i < num_centers; i++)
    {
        total_ang += angmom[i];
        if (min_ang > angmom[i] && angmom[i] != 0){
            min_ang = angmom[i];
            center = i;
            // std::cout<< i << min_ang <<"\n";
        }
    }
    if (total_ang == 0)
        return 'z' ; // means all orbitals are s orbital.

    else if (center == 0){
        return 'a';
    }
    else if (center == 1){
        return 'b';
    }
    else if (center == 2){
        return 'c';
    }
    else if (center == 3){
        return 'd';
    } // return the center having the lowest anglar momentum
}


char center_choose(int a, int b)
{
    if (b == 0)
        return 'b';
    else if ( (a == 0) and (b != 0))
        return 'd';
    else if (a < b)
        return 'b';
    else if (a >= b)
        return 'd';
}

bool intex(const std::array<int,3>& a, const std::array<int,3>& b, const std::array<int,3>& c, const std::array<int,3>& d) // integral exists or not
{
    if (a[0] < 0 || a[1] < 0 || a[2] < 0 || \
        b[0] < 0 || b[1] < 0 || b[2] < 0 || \
        c[0] < 0 || c[1] < 0 || c[2] < 0 || \
        d[0] < 0 || d[1] < 0 || d[2] < 0
        )
        return false;
    else
        return true;
}

int addrsear(std::array<int, 3> a)
{
    int la = a[0] + a[1] + a[2];
    int ia = 0;
    for (auto ax = la; ax >= 0; ax--)
    for (auto ay = la - ax; ay >= 0; ay--)
    {
        auto az = la - ax - ay;
        if (a == std::array<int, 3> {ax, ay, az})
            return ia;
        ia++;
    }
}

std::string orbname(const std::array<int, 3>& ang_mom)
{
    // return orbital name, like for (2,0,0) or (1,1,0), return d
    std::string s;
    int tot_ang = ang_mom[0] + ang_mom[1] + ang_mom[2];
    if (tot_ang == 0)
        s.append("s");
    else if (tot_ang == 1)
        s.append("p");
    else if (tot_ang == 2)
        s.append("d");
    else if (tot_ang == 3)
        s.append("f");
    else if (tot_ang == 4)
        s.append("g");
    else if (tot_ang == 5)
        s.append("h");
    else if (tot_ang == 6)
        s.append("i");
    return s;
}


void save_int(int la, int lb, int lc, int ld)
{
    auto lla = (la+1)*(la+2)/2;
    auto llb = (lb+1)*(lb+2)/2;
    auto llc = (lc+1)*(lc+2)/2;
    auto lld = (ld+1)*(ld+2)/2;
    auto llabcd = lla*llb*llc*lld;
    auto id = 0;
    for (auto dx = ld; dx >= 0; dx--)
    for (auto dy = ld-dx; dy >= 0; dy--)
    {
        auto dz = ld - dx - dy;
        auto ic = 0 ;
        for (auto cx = lc; cx >=0; cx-- )
        for (auto cy = lc-cx; cy >=0; cy-- )
        {
            auto cz = lc - cx -cy;
            auto ib = 0;
            for (auto bx = lb; bx >=0; bx--)
            for (auto by = lb-bx; by >=0; by--)
            {
                auto bz = lb - bx - by;
                auto ia = 0;
                for (auto ax = la; ax >=0; ax--)
                for (auto ay = la-ax; ay >=0; ay--)
                {
                    // auto ia = 0;
                    auto az = la - ax - ay;
                    std::array<int, 3> axyz = {ax,ay,az};
                    std::array<int, 3> bxyz = {bx,by,bz};
                    std::array<int, 3> cxyz = {cx,cy,cz};
                    std::array<int, 3> dxyz = {dx,dy,dz};
                    // cab + (ia+ib*lla+ic*lla*llb+id*lla*llb*llc)*nab + ccd*lla*llb*llc*lld*nab
                    printf("    I_[idx + (%d + %d * %d + %d * %d * %d + %d * %d * %d * %d) * nab + idy * nab * %d]", ia, ib, lla, ic, lla, llb, id, lla, llb, llc, llabcd);
                    std::cout << " += " << namemap(axyz,bxyz,cxyz,dxyz,0) << " ;"<< std::endl;
                    ia++;
                }
                ib++;
            }
            ic++;
        }
        id++;
    }
}


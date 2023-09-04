#include <iostream>
#include <set>
#include <limits.h>
#include <map>
#include "VRR_3c_2e.h"
#include "name.h"


std::string os_sent_gen(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c, int m, char center, int xyz)
{
    auto amins1 = a;
    auto amins2 = a;  
    auto bmins1 = b;  
    auto bmins2 = b;  
    auto cmins1 = c;  
    auto cmins2 = c;  

    std::string s;
    s.append("        auto ");
    std::array<int, 3> s_orb = {0, 0, 0};
    // base case for three center two electrons integrals
    if (a == s_orb and b == s_orb and c == s_orb)
    {
        s += namemap(a, b, c, m) + " = Fm[" + std::to_string(m) + "];\n";
        return s;
    }

    if (center == 'a')
    {
        amins1[xyz] = a[xyz] - 1;
        amins2[xyz] = a[xyz] - 2;
        bmins1[xyz] = b[xyz] - 1;
        cmins1[xyz] = c[xyz] - 1;

        s += namemap(a, b, c, m) + " = PA[" + std::to_string(xyz) + "] * " + namemap(amins1, b, c, m) + " + WP[" + std::to_string(xyz) + "] * " + namemap(amins1, b, c, m+1);

        if (intex(amins2, b, c))
        {
            s += " + 0.5 * zab_inv * " + std::to_string(amins1[xyz]) + " * " + namemap(amins2, b, c, m) + " + (- rho) * zab_inv * 0.5 * zab_inv * " + std::to_string(amins1[xyz]) + " * " + namemap(amins2, b, c, m+1);
        }

        if (intex(amins1, bmins1, c))
        {
            s += " + 0.5 * zab_inv * " + std::to_string(b[xyz]) + " * " + namemap(amins1, bmins1, c, m) + " + (- rho) * zab_inv * 0.5 * zab_inv * " + std::to_string(b[xyz]) + " * " + namemap(amins1, bmins1, c,  m+1); 
        }

        if (intex(amins1, b, cmins1))
        {
            s += " + 0.5 * zabcd_inv * " + std::to_string(c[xyz]) + " * " + namemap(amins1, b, cmins1, m+1);
        }
        s.append(";\n");
    }

    if (center == 'b')
    {
        amins1[xyz] = a[xyz] - 1;
        bmins1[xyz] = b[xyz] - 1;
        bmins2[xyz] = b[xyz] - 2;
        cmins1[xyz] = c[xyz] - 1;

        s += namemap(a, b, c, m) + " = PB[" + std::to_string(xyz) + "] * " + namemap(a, bmins1, c, m) + " + WP[" + std::to_string(xyz) + "] * " + namemap(a, bmins1, c, m+1);

        if (intex(a, bmins2, c))
        {
            s += " + 0.5 * zab_inv * " + std::to_string(bmins1[xyz]) + " * " + namemap(a, bmins2, c, m) + " + (- rho) * zab_inv * 0.5 * zab_inv * " + std::to_string(bmins1[xyz]) + " * " + namemap(a, bmins2, c, m+1);
        }

        if (intex(amins1, bmins1, c))
        {
            s += " + 0.5 * zab_inv * " + std::to_string(a[xyz]) + " * " + namemap(amins1, bmins1, c, m) + " + (- rho) * zab_inv * 0.5 * zab_inv * " + std::to_string(a[xyz]) + " * " + namemap(amins1, bmins1, c, m+1);
        }

        if (intex(a, bmins1, cmins1))
        {
            s += " + 0.5 * zabcd_inv * " + std::to_string(c[xyz]) + " * " + namemap(a, bmins1, cmins1, m+1);
        }
        s.append(";\n");
    }

    if (center == 'c')
    {
        amins1[xyz] = a[xyz] - 1;
        bmins1[xyz] = b[xyz] - 1;
        cmins1[xyz] = c[xyz] - 1;
        cmins2[xyz] = c[xyz] - 2;

        s += namemap(a,b,c,m) + " = QC[" + std::to_string(xyz) + "] * " + namemap(a,b,cmins1,m) + " + WQ[" + std::to_string(xyz) + "] * " + namemap(a,b,cmins1,m+1);

        if (intex(a, b, cmins2))
        {
            s += " + 0.5 * zcd_inv * " + std::to_string(cmins1[xyz]) + " * " + namemap(a, b, cmins2, m) + " + (- rho) * zcd_inv * 0.5 * zcd_inv * " + std::to_string(cmins1[xyz]) + " * " + namemap(a, b, cmins2, m+1);
        }

        if (intex(amins1, b, cmins1))
        {
            s += " + 0.5 * zabcd_inv * " + std::to_string(a[xyz]) + " * " + namemap(amins1, b, cmins1, m+1);
        }

        if (intex(a, bmins1, cmins1))
        {
            s += " + 0.5 * zabcd_inv * " + std::to_string(b[xyz]) + " * " + namemap(a, bmins1, cmins1, m+1);
        }
        s.append(";\n");
    }
    return s;
}



void obara_saika(const std::array<int, 3>& a, const std::array<int, 3>& b,const std::array<int, 3>& c, int m, std::map<std::array<int, 10>, std::string>& osmap)
{
    
    int xyz;
    std::array<int, 10> os_id = {a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], m};
    
    if (osmap.count(os_id) != 0)
        return ;
    
    if ((a == std::array<int, 3> {0,0,0}) and
        (b == std::array<int, 3> {0,0,0}) and 
        (c == std::array<int, 3> {0,0,0})) 
    {
        auto hello = os_sent_gen(a,b,c,m,'a',0);
        // std::cout << "auto " << hello ;
        if (osmap.count(os_id) == 0)
        {
            osmap.insert({os_id, hello});
        }
        return ;
    }
    else if (!intex(a,b,c)) // if a,b,c,d don't exist
    {
        return ;
    }
    std::array<int, 3> lambda_abc = {(a[0] + a[1] + a[2]),
                                      (b[0] + b[1] + b[2]),
                                      (c[0] + c[1] + c[2])};

    char center_start = center_choose(lambda_abc);
    
    auto amins1 = a;
    auto amins2 = a; 
    auto bmins1 = b; 
    auto bmins2 = b; 
    auto cmins1 = c; 
    auto cmins2 = c; 
    
    // A center
    if (center_start == 'a')
    {
        // determine x,y,z direction to decrease
        auto xyz = dirchoose(a);
        // for (auto i = 0; i < 3; i++)
        // {
        //     if (a[i] != 0)
        //     {
        //         xyz = i;
        //         break;
        //     }
        // }
        amins1[xyz] = a[xyz] - 1;
        amins2[xyz] = a[xyz] - 2;
        bmins1[xyz] = b[xyz] - 1;
        cmins1[xyz] = c[xyz] - 1;

        // s.append(namemap(a,b,c,d,m)); 
        auto hello = os_sent_gen(a,b,c,m,center_start, xyz);
        if (osmap.count(os_id) == 0)
        {
            osmap.insert({os_id, hello});
        }

        //std::cout << "auto " << hello << std::endl;

        obara_saika(amins1,      b,      c,       m  , osmap);
        obara_saika(amins1,      b,      c,       m+1, osmap);
        obara_saika(amins2,      b,      c,       m  , osmap);
        obara_saika(amins2,      b,      c,       m+1, osmap);
        obara_saika(amins1,      bmins1, c,       m+1, osmap);
        obara_saika(amins1,      b,      cmins1,  m+1, osmap);
    }
    // B center
    else if (center_start == 'b')
    {
        // determine x,y,z direction to decrease
        auto xyz = dirchoose(b);

        amins1[xyz] = a[xyz] - 1;
        bmins1[xyz] = b[xyz] - 1;
        bmins2[xyz] = b[xyz] - 2;
        cmins1[xyz] = c[xyz] - 1;

        // s.append(namemap(a,b,c,d,m)); 
        auto hello = os_sent_gen(a,b,c,m,center_start, xyz);
        if (osmap.count(os_id) == 0)
        {
            osmap.insert({os_id, hello});
        }

        //std::cout << "auto " << hello << std::endl;

        obara_saika(a     ,      bmins1, c,       m  , osmap);
        obara_saika(a     ,      bmins1, c,       m+1, osmap);
        obara_saika(a     ,      bmins2, c,       m  , osmap);
        obara_saika(a     ,      bmins2, c,       m+1, osmap);
        obara_saika(amins1,      bmins1, c,       m+1, osmap);
        obara_saika(a     ,      bmins1, cmins1,  m+1, osmap);
    }
    // C center
    else if (center_start == 'c')
    {
        // determine x,y,z direction to decrease
        auto xyz = dirchoose(c);

        amins1[xyz] = a[xyz] - 1;
        bmins1[xyz] = b[xyz] - 1;
        cmins1[xyz] = c[xyz] - 1;
        cmins2[xyz] = c[xyz] - 2;

        // s.append(namemap(a,b,c,d,m)); 
        auto hello = os_sent_gen(a,b,c,m,center_start, xyz);
        if (osmap.count(os_id) == 0)
        {
            osmap.insert({os_id, hello});
        }

        //std::cout << "auto " << hello << std::endl;

        obara_saika(a     ,      b, cmins1,        m  , osmap);
        obara_saika(a     ,      b, cmins1,        m+1, osmap);
        obara_saika(a     ,      b, cmins2,        m  , osmap);
        obara_saika(a     ,      b, cmins2,        m+1, osmap);
        obara_saika(amins1,      b, cmins1,        m+1, osmap);
        obara_saika(a     , bmins1, cmins1,        m+1, osmap);
    }
}


static void code_print(std::array<int, 3>& a , std::array<int, 3>& b, std::array<int, 3>& c, int m, std::map<std::array<int, 10>, std::string>& osmap)
{
    int labcm = a[0] + a[1] + a[2] + b[0] + b[1] + b[2] + c[0] + c[1] + c[2] + m; 
    for (auto ax = 0; ax <= a[0]; ax++)
    for (auto ay = 0; ay <= a[1]; ay++)
    for (auto az = 0; az <= a[2]; az++)
    for (auto bx = 0; bx <= b[0]; bx++)
    for (auto by = 0; by <= b[1]; by++)
    for (auto bz = 0; bz <= b[2]; bz++)
    for (auto cx = 0; cx <= c[0]; cx++)
    for (auto cy = 0; cy <= c[1]; cy++)
    for (auto cz = 0; cz <= c[2]; cz++)
    for (auto mi = 0; mi <= labcm; mi++)
    {
        std::array<int, 10> os_id = {ax, ay, az, bx, by, bz, cx, cy, cz, mi};
        if (osmap.count(os_id) != 0)
            std::cout << osmap.at(os_id); // << "hello" << std::endl;
    }
}


void vrr_code_print(int la , int lb, int lc, std::map<std::array<int, 10>, std::string>& osmap)
{
    // int labcdm = a[0] + a[1] + a[2] + b[0] + b[1] + b[2] + c[0] + c[1] + c[2] + d[0] + d[1] + d[2]; 
    int labcdm = la + lb + lc; 
    // for (auto llabcdm = 0; llabcdm <= labcdm; llabcdm++)
    for (auto lla = 0; lla <= la; lla++) 
    for (auto llb = 0; llb <= lb; llb++) 
    for (auto llc = 0; llc <= lc; llc++) 
    {
        // printf("%d,%d,%d,%d\n", lla, llb,llc, lld);
        for (auto ax = lla; ax >= 0; ax--)
        for (auto ay = lla - ax; ay >= 0; ay--)
        // for (auto az = 0; az <= a[2]; az++)
        for (auto bx = llb; bx >= 0; bx--)
        for (auto by = llb - bx; by >= 0; by--)
        // for (auto bz = 0; bz <= b[2]; bz++)
        for (auto cx = llc; cx >= 0; cx--)
        for (auto cy = llc - cx; cy >= 0; cy--)
        // for (auto cz = 0; cz <= c[2]; cz++)
        for (auto mi = 0; mi <= labcdm; mi++)
        {
            std::array<int, 10> os_id = {ax, ay, lla-ax-ay, bx, by, llb-bx-by, cx, cy, llc-cx-cy, mi};
            if (osmap.count(os_id) != 0)
                std::cout << osmap.at(os_id); // << "hello" << std::endl;
        }
    }
}


void eri(int la, int lb, int lc, std::map<std::array<int, 10>, std::string>& osmap)
{
    // auto id = 0;
    // auto ib = 0 ;
    for (auto cx = lc; cx >=0; cx-- )
    for (auto cy = lc-cx; cy >=0; cy-- )
    {
        auto cz = lc - cx -cy;
        // auto ic = 0;
        for (auto bx = lb; bx >=0; bx--)
        for (auto by = lb-bx; by >=0; by--)
        {
            auto bz = lb - bx - by;
            // auto ib = 0;
            for (auto ax = la; ax >=0; ax--)
            for (auto ay = la-ax; ay >=0; ay--)
            {
                // auto ia = 0;
                auto az = la - ax - ay;
                std::array<int, 3> axyz = {ax,ay,az};
                std::array<int, 3> bxyz = {bx,by,bz};
                std::array<int, 3> cxyz = {cx,cy,cz};
                obara_saika(axyz,bxyz,cxyz,0, osmap);
            }
        }
    }
}



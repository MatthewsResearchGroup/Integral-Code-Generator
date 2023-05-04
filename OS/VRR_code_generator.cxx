#include <iostream>
#include <set>
// #include <array>
#include <limits.h>
#include <map>
#include "VRR_code_generator.h"

void f(int num);

#if 0
std::string xyz_tran(const std::array<int, 3>& ang_mom)
{
    // s orbital
    std::array<int, 3> s = {0,0,0};
    // p orbital
    std::array<int, 3> x = {1,0,0};
    std::array<int, 3> y = {0,1,0};
    std::array<int, 3> z = {0,0,1};
    // d orbital
    std::array<int, 3> xx = {2,0,0};
    std::array<int, 3> yy = {0,2,0};
    std::array<int, 3> zz = {0,0,2};
    std::array<int, 3> xy = {1,1,0};
    std::array<int, 3> xz = {1,0,1};
    std::array<int, 3> yz = {0,1,1};
    // f orbital
    std::array<int, 3> xxx = {3,0,0};
    std::array<int, 3> yyy = {0,3,0};
    std::array<int, 3> zzz = {0,0,3};
    std::array<int, 3> xxy = {2,1,0};
    std::array<int, 3> xxz = {2,0,1};
    std::array<int, 3> xyy = {1,2,0};
    std::array<int, 3> xzz = {1,0,2};
    std::array<int, 3> xyz = {1,1,1};
    std::array<int, 3> yyz = {0,2,1};
    std::array<int, 3> yzz = {0,1,2};
    // g orbital
    std::array<int, 3> xxxx = {4,0,0};
    std::array<int, 3> yyyy = {0,4,0};
    std::array<int, 3> zzzz = {0,0,4};
    std::array<int, 3> xxxy = {3,1,0};
    std::array<int, 3> xxxz = {3,0,1};
    std::array<int, 3> yyyz = {0,3,1};
    std::array<int, 3> xyyy = {1,3,0};
    std::array<int, 3> xzzz = {1,0,3};
    std::array<int, 3> yzzz = {0,1,3};
    std::array<int, 3> xxyy = {2,2,0};
    std::array<int, 3> yyzz = {0,2,2};
    std::array<int, 3> xxzz = {2,0,2};
    std::array<int, 3> xxyz = {2,1,1};
    std::array<int, 3> xyyz = {1,2,1};
    std::array<int, 3> xyzz = {1,1,2};

    if (ang_mom == s)
        return "s_";
    // p orbital
    else if (ang_mom == x)
        return "x_";
    else if (ang_mom == y)
        return "y_";
    else if (ang_mom == z)
        return "z_";
    // d orbital
    else if (ang_mom == xx)
        return "xx_";
    else if (ang_mom == yy)
        return "yy_";
    else if (ang_mom == zz)
        return "zz_";
    else if (ang_mom == xy)
        return "xy_";
    else if (ang_mom == xz)
        return "xz_";
    else if (ang_mom == yz)
        return "yz_";
    // f orbital
    else if (ang_mom == xxx)
        return "xxx_";
    else if (ang_mom == yyy)
        return "yyy_";
    else if (ang_mom == zzz)
        return "zzz_";
    else if (ang_mom == xxy)
        return "xxy_";
    else if (ang_mom == xxz)
        return "xxz_";
    else if (ang_mom == xyy)
        return "xyy_";
    else if (ang_mom == xzz)
        return "xzz_";
    else if (ang_mom == yyz)
        return "yyz_";
    else if (ang_mom == yzz)
        return "yzz_";
    else if (ang_mom == xyz)
        return "xyz_";
    // g orbital
    else if (ang_mom == xxxx)
        return "xxxx_";
    else if (ang_mom == yyyy)
        return "yyyy_";
    else if (ang_mom == zzzz)
        return "zzzz_";
    else if (ang_mom == xxxy)
        return "xxxy_";
    else if (ang_mom == xxxz)
        return "xxxz_";
    else if (ang_mom == xxyy)
        return "xxyy_";
    else if (ang_mom == xxyz)
        return "xxyz_";
    else if (ang_mom == xxzz)
        return "xxzz_";
    else if (ang_mom == xyyy)
        return "xyyy_";
    else if (ang_mom == xyyz)
        return "xyyz_";
    else if (ang_mom == xyzz)
        return "xyzz_";
    else if (ang_mom == xzzz)
        return "xzzz_";
    else if (ang_mom == yyyz)
        return "yyyz_";
    else if (ang_mom == yyzz)
        return "yyzz_";
    else if (ang_mom == yzzz)
        return "yzzz_";
    else
        return "Error";
}
#endif

static std::string NewNameScheme(const std::array<int, 3>& ang_mom)
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

static int dirchoose(const std::array<int, 3> a)
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

static std::string namemap(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c,const std::array<int, 3> d, int m)
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

static char center_decrease(std::array<int, 4> angmom)
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
 
static bool intex(const std::array<int,3>& a, const std::array<int,3>& b, const std::array<int,3>& c, const std::array<int,3>& d) // integral exists or not
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

static std::string sent_gen(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c,const std::array<int, 3> d, int m, char center, int xyz)
{
    auto amins1 = a;
    auto amins2 = a;  
    auto bmins1 = b;  
    auto bmins2 = b;  
    auto cmins1 = c;  
    auto cmins2 = c;  
    auto dmins1 = d;  
    auto dmins2 = d;
    
    std::string s;
    s.append("        auto ");
    std::array<int,3> s_orbital {0,0,0};
    if (a == s_orbital and b == s_orbital and c == s_orbital and d == s_orbital)
    {
        s.append(namemap(a,b,c,d,m));
        s.append(" = ");
        s.append("fm[");
        s.append(std::to_string(m));
        s.append("] ;\n");
        return s;
    }

    // A center
    if (center == 'a') 
    {
        amins1[xyz] = a[xyz] - 1;
        amins2[xyz] = a[xyz] - 2;
        bmins1[xyz] = b[xyz] - 1;
        cmins1[xyz] = c[xyz] - 1;
        dmins1[xyz] = d[xyz] - 1;
    
        s.append(namemap(a,b,c,d,m));
        s.append(" = ");
        s.append("PA[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(amins1,b,c,d,m));
        s.append(" + ");
        s.append("WP[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(amins1,b,c,d,m+1));
        if (intex(amins2, b, c, d))
        {
            s.append(" + ");
            s.append("0.5 * zab_inv * ");
            s.append(std::to_string(amins1[xyz]));
            s.append(" * ");
            s.append(namemap(amins2, b, c, d, m));
            s.append(" + ");
            s.append("(- rho) * zab_inv * 0.5 * zab_inv * ");
            s.append(std::to_string(amins1[xyz]));
            s.append(" * ");
            s.append(namemap(amins2, b, c, d, m+1));
        }
        if (intex(amins1, bmins1, c, d))
        {
            s.append(" + ");
            s.append("0.5 * zab_inv * ");
            s.append(std::to_string(b[xyz]));
            s.append(" * ");
            s.append(namemap(amins1, bmins1, c, d, m));
            s.append(" + ");
            s.append("(- rho) * zab_inv * 0.5 * zab_inv * ");
            s.append(std::to_string(b[xyz]));
            s.append(" * ");
            s.append(namemap(amins1, bmins1, c, d, m+1));
        }
        if (intex(amins1, b, cmins1, d))
        {
            // s.append(" + ");
            // s.append(" 0.0 * ");
            // s.append(namemap(amins1, b, cmins1, d, m));
            s.append(" + ");
            s.append("0.5 * zabcd_inv * ");
            s.append(std::to_string(c[xyz]));
            s.append(" * ");
            s.append(namemap(amins1, b, cmins1, d, m+1));
        }
        if (intex(amins1, b, c, dmins1))
        {
            // s.append(" + ");
            // s.append(" 0.0 * ");
            // s.append(namemap(amins1, b, c, dmins1, m));
            s.append(" + ");
            s.append("0.5 * zabcd_inv * ");
            s.append(std::to_string(d[xyz]));
            s.append(" * ");
            s.append(namemap(amins1, b, c, dmins1, m+1));
        }
        s.append(";\n");
    }
    
    // B center
    
    if (center == 'b') 
    {
        amins1[xyz] = a[xyz] - 1; 
        bmins1[xyz] = b[xyz] - 1; 
        bmins2[xyz] = b[xyz] - 2; 
        cmins1[xyz] = c[xyz] - 1; 
        dmins1[xyz] = d[xyz] - 1; 
    
        s.append(namemap(a,b,c,d,m));
        s.append(" = ");
        s.append("PB[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(a,bmins1,c,d,m));
        s.append(" + ");
        s.append("WP[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(a,bmins1,c,d,m+1));
        if (intex(a, bmins2, c, d))
        {
            s.append(" + ");
            s.append("0.5 * zab_inv * ");
            s.append(std::to_string(bmins1[xyz]));
            s.append(" * ");
            s.append(namemap(a, bmins2, c, d, m));
            s.append(" + ");
            s.append("(- rho) * zab_inv * 0.5 * zab_inv * ");
            s.append(std::to_string(bmins1[xyz]));
            s.append(" * ");
            s.append(namemap(a, bmins2, c, d, m+1));
        }
        if (intex(amins1, bmins1, c, d))
        {
            s.append(" + ");
            s.append("0.5 * zab_inv * ");
            s.append(std::to_string(a[xyz]));
            s.append(" * ");
            s.append(namemap(amins1, bmins1, c, d, m));
            s.append(" + ");
            s.append("(- rho) * zab_inv * 0.5 * zab_inv * ");
            s.append(std::to_string(a[xyz]));
            s.append(" * ");
            s.append(namemap(amins1, bmins1, c, d, m+1));
        }
        if (intex(a, bmins1, cmins1, d))
        {
            s.append(" + ");
            s.append(" 0.0 * ");
            s.append(namemap(a, b, cmins1, d, m));
            s.append(" + ");
            s.append("0.5 * zabcd_inv * ");
            s.append(std::to_string(c[xyz]));
            s.append(" * ");
            s.append(namemap(amins1, b, cmins1, d, m+1));
        }
        if (intex(a, bmins1, c, dmins1))
        {
            s.append(" + ");
            s.append(" 0.0 * ");
            s.append(namemap(amins1, b, c, dmins1, m));
            s.append(" + ");
            s.append("0.5 * zabcd_inv * ");
            s.append(std::to_string(d[xyz]));
            s.append(" * ");
            s.append(namemap(amins1, b, c, dmins1, m+1));
        }
        s.append(";\n");
    }
    // C center
    if (center == 'c') 
    {
        amins1[xyz] = a[xyz] - 1;
        bmins1[xyz] = b[xyz] - 1;
        cmins1[xyz] = c[xyz] - 1;
        cmins2[xyz] = c[xyz] - 2;
        dmins1[xyz] = d[xyz] - 1;
    
        s.append(namemap(a,b,c,d,m));
        s.append(" = ");
        s.append("QC[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(a,b,cmins1,d,m));
        s.append(" + ");
        s.append("WQ[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(a,b,cmins1,d,m+1));
        if (intex(a, b, cmins2, d))
        {
            s.append(" + ");
            s.append("0.5 * zcd_inv * ");
            s.append(std::to_string(cmins1[xyz]));
            s.append(" * ");
            s.append(namemap(a, b, cmins2, d, m));
            s.append(" + ");
            s.append("(- rho) * zcd_inv * 0.5 * zcd_inv * ");
            s.append(std::to_string(cmins1[xyz]));
            s.append(" * ");
            s.append(namemap(a, b, cmins2, d, m+1));
        }
        if (intex(a, b, cmins1, dmins1))
        {
            s.append(" + ");
            s.append("0.5 * zcd_inv * ");
            s.append(std::to_string(d[xyz]));
            s.append(" * ");
            s.append(namemap(a, b, cmins1, dmins1, m));
            s.append(" + ");
            s.append("(- rho) * zcd_inv * 0.5 * zcd_inv * ");
            s.append(std::to_string(d[xyz]));
            s.append(" * ");
            s.append(namemap(a, b, cmins1, dmins1, m+1));
        }
        if (intex(amins1, b, cmins1, d))
        {
            s.append(" + ");
            s.append(" 0.0 * ");
            s.append(namemap(amins1, b, cmins1, d, m));
            s.append(" + ");
            s.append("0.5 * zabcd_inv * ");
            s.append(std::to_string(a[xyz]));
            s.append(" * ");
            s.append(namemap(amins1, b, cmins1, d, m+1));
        }
        if (intex(a, bmins1, cmins1, d))
        {
            s.append(" + ");
            s.append(" 0.0 * ");
            s.append(namemap(a, bmins1, cmins1, d, m));
            s.append(" + ");
            s.append("0.5 * zabcd_inv * ");
            s.append(std::to_string(b[xyz]));
            s.append(" * ");
            s.append(namemap(a, bmins1, cmins1, d, m+1));
        }
        s.append(";\n");
    }
    // D center 
    if (center == 'd') 
    {
        amins1[xyz] = a[xyz] - 1;
        bmins1[xyz] = b[xyz] - 1;
        cmins1[xyz] = c[xyz] - 1;
        dmins1[xyz] = d[xyz] - 1;
        dmins2[xyz] = d[xyz] - 2;
    
        s.append(namemap(a,b,c,d,m));
        s.append(" = ");
        s.append("QD[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(a,b,c,dmins1,m));
        s.append(" + ");
        s.append("WQ[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(a,b,c,dmins1,m+1));
        if (intex(a, b, c, dmins2))
        {
            s.append(" + ");
            s.append("0.5 * zcd_inv * ");
            s.append(std::to_string(dmins1[xyz]));
            s.append(" * ");
            s.append(namemap(a, b, c, dmins2, m));
            s.append(" + ");
            s.append("(- rho) * zcd_inv * 0.5 * zcd_inv * ");
            s.append(std::to_string(dmins1[xyz]));
            s.append(" * ");
            s.append(namemap(a, b, c, dmins2, m+1));
        }
        if (intex(a, b, cmins1, dmins1))
        {
            s.append(" + ");
            s.append("0.5 * zcd_inv * ");
            s.append(std::to_string(c[xyz]));
            s.append(" * ");
            s.append(namemap(a, b, cmins1, dmins1, m));
            s.append(" + ");
            s.append("(- rho) * zcd_inv * 0.5 * zcd_inv * ");
            s.append(std::to_string(c[xyz]));
            s.append(" * ");
            s.append(namemap(a, b, cmins1, dmins1, m+1));
        }
        if (intex(amins1, b, c, dmins1))
        {
            s.append(" + ");
            s.append(" 0.0 * ");
            s.append(namemap(amins1, b, c, dmins1, m));
            s.append(" + ");
            s.append("0.5 * zabcd_inv * ");
            s.append(std::to_string(a[xyz]));
            s.append(" * ");
            s.append(namemap(amins1, b, c, dmins1, m+1));
        }
        if (intex(a, bmins1, c, dmins1))
        {
            s.append(" + ");
            s.append(" 0.0 * ");
            s.append(namemap(a, bmins1, c, dmins1, m));
            s.append(" + ");
            s.append("0.5 * zabcd_inv * ");
            s.append(std::to_string(b[xyz]));
            s.append(" * ");
            s.append(namemap(a, bmins1, c, dmins1, m+1));
        }
        s.append(";\n");
    }
    
    return s;
       
}

void obara_saika(const std::array<int, 3>& a, const std::array<int, 3>& b,const std::array<int, 3>& c,const std::array<int, 3>& d, int m, std::map<std::array<int, 13>, std::string>& osmap)
{
    
    int xyz;
    std::array<int, 13> os_id = {a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], d[0], d[1], d[2], m};
    
    if (osmap.count(os_id) != 0)
        return ;
    
    if ((a == std::array<int, 3> {0,0,0}) and
        (b == std::array<int, 3> {0,0,0}) and 
        (c == std::array<int, 3> {0,0,0}) and 
        (d == std::array<int, 3> {0,0,0}))
    {
        auto hello = sent_gen(a,b,c,d,m,'a',0);
        // std::cout << "auto " << hello ;
        if (osmap.count(os_id) == 0)
        {
            osmap.insert({os_id, hello});
        }
        return ;
    }
    else if (!intex(a,b,c,d)) // if a,b,c,d don't exist
    {
        return ;
    }
    std::array<int, 4> lambda_abcd = {(a[0] + a[1] + a[2]),
                                      (b[0] + b[1] + b[2]),
                                      (c[0] + c[1] + c[2]),
                                      (d[0] + d[1] + d[2])};

    char center_start = center_decrease(lambda_abcd);
    
    auto amins1 = a;
    auto amins2 = a; 
    auto bmins1 = b; 
    auto bmins2 = b; 
    auto cmins1 = c; 
    auto cmins2 = c; 
    auto dmins1 = d; 
    auto dmins2 = d;
    
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
        dmins1[xyz] = d[xyz] - 1;

        // s.append(namemap(a,b,c,d,m)); 
        auto hello = sent_gen(a,b,c,d,m,center_start, xyz);
        if (osmap.count(os_id) == 0)
        {
            osmap.insert({os_id, hello});
        }

        //std::cout << "auto " << hello << std::endl;

        obara_saika(amins1,      b,      c,      d,        m  , osmap);
        obara_saika(amins1,      b,      c,      d,        m+1, osmap);
        obara_saika(amins2,      b,      c,      d,        m  , osmap);
        obara_saika(amins2,      b,      c,      d,        m+1, osmap);
        obara_saika(amins1,      bmins1, c,      d,        m  , osmap);
        obara_saika(amins1,      bmins1, c,      d,        m+1, osmap);
        obara_saika(amins1,      b,      cmins1, d,        m  , osmap);
        obara_saika(amins1,      b,      cmins1, d,        m+1, osmap);
        obara_saika(amins1,      b,      c,      dmins1,   m  , osmap);
        obara_saika(amins1,      b,      c,      dmins1,   m+1, osmap);
    }
    // B center
    else if (center_start == 'b')
    {
        // determine x,y,z direction to decrease
        auto xyz = dirchoose(b);
        // for (auto i = 0; i < 3; i++)
        // {
        //     if (b[i] != 0)
        //     {
        //         xyz = i;
        //         break;
        //     }
        // }
        amins1[xyz] = a[xyz] - 1;
        bmins1[xyz] = b[xyz] - 1;
        bmins2[xyz] = b[xyz] - 2;
        cmins1[xyz] = c[xyz] - 1;
        dmins1[xyz] = d[xyz] - 1;

        // s.append(namemap(a,b,c,d,m)); 
        auto hello = sent_gen(a,b,c,d,m,center_start, xyz);
        if (osmap.count(os_id) == 0)
        {
            osmap.insert({os_id, hello});
        }

        //std::cout << "auto " << hello << std::endl;

        obara_saika(a     ,      bmins1, c,      d,        m  , osmap);
        obara_saika(a     ,      bmins1, c,      d,        m+1, osmap);
        obara_saika(a     ,      bmins2, c,      d,        m  , osmap);
        obara_saika(a     ,      bmins2, c,      d,        m+1, osmap);
        obara_saika(amins1,      bmins1, c,      d,        m  , osmap);
        obara_saika(amins1,      bmins1, c,      d,        m+1, osmap);
        obara_saika(a     ,      bmins1, cmins1, d,        m  , osmap);
        obara_saika(a     ,      bmins1, cmins1, d,        m+1, osmap);
        obara_saika(a     ,      bmins1, c,      dmins1,   m  , osmap);
        obara_saika(a     ,      bmins1, c,      dmins1,   m+1, osmap);
    }
    // C center
    else if (center_start == 'c')
    {
        // determine x,y,z direction to decrease
        auto xyz = dirchoose(c);
        // for (auto i = 0; i < 3; i++)
        // {
        //     if (c[i] != 0)
        //     {
        //         xyz = i;
        //         break;
        //     }
        // }
        amins1[xyz] = a[xyz] - 1;
        bmins1[xyz] = b[xyz] - 1;
        cmins1[xyz] = c[xyz] - 1;
        cmins2[xyz] = c[xyz] - 2;
        dmins1[xyz] = d[xyz] - 1;

        // s.append(namemap(a,b,c,d,m)); 
        auto hello = sent_gen(a,b,c,d,m,center_start, xyz);
        if (osmap.count(os_id) == 0)
        {
            osmap.insert({os_id, hello});
        }

        //std::cout << "auto " << hello << std::endl;

        obara_saika(a     ,      b, cmins1,      d,        m  , osmap);
        obara_saika(a     ,      b, cmins1,      d,        m+1, osmap);
        obara_saika(a     ,      b, cmins2,      d,        m  , osmap);
        obara_saika(a     ,      b, cmins2,      d,        m+1, osmap);
        obara_saika(a     ,      b, cmins1, dmins1,        m  , osmap);
        obara_saika(a     ,      b, cmins1, dmins1,        m+1, osmap);
        obara_saika(amins1,      b, cmins1,      d,        m  , osmap);
        obara_saika(amins1,      b, cmins1,      d,        m+1, osmap);
        obara_saika(a     , bmins1, cmins1,      d,        m  , osmap);
        obara_saika(a     , bmins1, cmins1,      d,        m+1, osmap);
    }
    // D center
    else if (center_start == 'd')
    {
        // determine x,y,z direction to decrease
        auto xyz = dirchoose(d);
        // for (auto i = 0; i < 3; i++)
        // {
        //     if (d[i] != 0)
        //     {
        //         xyz = i;
        //         break;
        //     }
        // }
        amins1[xyz] = a[xyz] - 1;
        bmins1[xyz] = b[xyz] - 1;
        cmins1[xyz] = c[xyz] - 1;
        dmins1[xyz] = d[xyz] - 1;
        dmins2[xyz] = d[xyz] - 2;

        // s.append(namemap(a,b,c,d,m)); 
        auto hello = sent_gen(a,b,c,d,m,center_start, xyz);
        if (osmap.count(os_id) == 0)
        {
            osmap.insert({os_id, hello});
        }

        obara_saika(a     ,      b,      c, dmins1,        m  , osmap);
        obara_saika(a     ,      b,      c, dmins1,        m+1, osmap);
        obara_saika(a     ,      b,      c, dmins2,        m  , osmap);
        obara_saika(a     ,      b,      c, dmins2,        m+1, osmap);
        obara_saika(a     ,      b, cmins1, dmins1,        m  , osmap);
        obara_saika(a     ,      b, cmins1, dmins1,        m+1, osmap);
        obara_saika(amins1,      b,      c, dmins1,        m  , osmap);
        obara_saika(amins1,      b,      c, dmins1,        m+1, osmap);
        obara_saika(a     , bmins1,      c, dmins1,        m  , osmap);
        obara_saika(a     , bmins1,      c, dmins1,        m+1, osmap);
    }
             
}

static void code_print(std::array<int, 3>& a , std::array<int, 3>& b, std::array<int, 3>& c, std::array<int, 3>& d, int m, std::map<std::array<int, 13>, std::string>& osmap)
{
    int labcdm = a[0] + a[1] + a[2] + b[0] + b[1] + b[2] + c[0] + c[1] + c[2] + d[0] + d[1] + d[2] + m; 
    for (auto ax = 0; ax <= a[0]; ax++)
    for (auto ay = 0; ay <= a[1]; ay++)
    for (auto az = 0; az <= a[2]; az++)
    for (auto bx = 0; bx <= b[0]; bx++)
    for (auto by = 0; by <= b[1]; by++)
    for (auto bz = 0; bz <= b[2]; bz++)
    for (auto cx = 0; cx <= c[0]; cx++)
    for (auto cy = 0; cy <= c[1]; cy++)
    for (auto cz = 0; cz <= c[2]; cz++)
    for (auto dx = 0; dx <= d[0]; dx++)
    for (auto dy = 0; dy <= d[1]; dy++)
    for (auto dz = 0; dz <= d[2]; dz++)
    for (auto mi = 0; mi <= labcdm; mi++)
    {
        std::array<int, 13> os_id = {ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, mi};
        if (osmap.count(os_id) != 0)
            std::cout << osmap.at(os_id); // << "hello" << std::endl;
    }
}


void vrr_code_print(int la , int lb, int lc, int ld, std::map<std::array<int, 13>, std::string>& osmap)
{
    // int labcdm = a[0] + a[1] + a[2] + b[0] + b[1] + b[2] + c[0] + c[1] + c[2] + d[0] + d[1] + d[2]; 
    int labcdm = la + lb + lc + ld;
    // for (auto llabcdm = 0; llabcdm <= labcdm; llabcdm++)
    for (auto lla = 0; lla <= la; lla++) 
    for (auto llb = 0; llb <= lb; llb++) 
    for (auto llc = 0; llc <= lc; llc++) 
    for (auto lld = 0; lld <= ld; lld++) 
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
        for (auto dx = lld; dx >= 0; dx--)
        for (auto dy = lld - dx; dy >= 0; dy--)
        // for (auto dz = 0; dz <= d[2]; dz++)
        for (auto mi = 0; mi <= labcdm; mi++)
        {
            // printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", ax,ay,lla-ax-ay, bx, by, llb-bx-by, cx, cy, llc-cx-cy, dx, dy, lld-dx-dy, mi);
            std::array<int, 13> os_id = {ax, ay, lla-ax-ay, bx, by, llb-bx-by, cx, cy, llc-cx-cy, dx, dy, lld-dx-dy, mi};
            if (osmap.count(os_id) != 0)
                std::cout << osmap.at(os_id); // << "hello" << std::endl;
        }
    }
}

static void save_int(int la, int lb, int lc, int ld)
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


void eri(int la, int lb, int lc, int ld, std::map<std::array<int, 13>, std::string>& osmap)
{
    // auto id = 0;
    for (auto dx = ld; dx >= 0; dx--)
    for (auto dy = ld-dx; dy >= 0; dy--)
    {
        auto dz = ld - dx - dy;
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
                    std::array<int, 3> dxyz = {dx,dy,dz};
                    obara_saika(axyz,bxyz,cxyz,dxyz,0, osmap);
                }
            }
        }
    }
}

#if 0
int main()
{
    std::array<int, 3> a {0,1,1};
    std::array<int, 3> b {0,0,0};
    std::array<int, 3> c {0,0,0};
    std::array<int, 3> d {0,0,0};
    // std::cout << dirchoose(b) << std::endl;
    // auto a_name = NewNameScheme(a);
    // std::cout << "a_name = " << a_name << std::endl;
    // int m = 0;
    // auto name_str = namemap(a,b,c,d,1);
    // std::cout << name_str << std::endl;

    int la = 2 ;
    int lb = 1 ;
    int lc = 0 ;
    int ld = 0 ;
    std::map<std::array<int, 13>, std::string> osmap;
    eri(la,lb,lc,ld,osmap);
    code_print(la , lb, lc, ld, osmap);
    save_int(la,lb,lc,ld);
    // obara_saika(a,b,c,d,m, osmap);
    // code_print(a, b, c, d, m, osmap);
    return 0;
}
#endif

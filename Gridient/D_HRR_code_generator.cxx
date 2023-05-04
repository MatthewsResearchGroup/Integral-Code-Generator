#include <map>
#include <array>
#include <string>
#include "name.h"
#include <iostream>
#include "D_HRR_code_generator.h"

std::string sent_gen(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c,const std::array<int, 3> d, int m, char center, int xyz, int deri_center)
// deri_center: the center where the derivation applies on, alpha means center A, beta means center B, gamma means center C
{
    auto aplus1 = a;
    auto bmins1 = b;
    auto cplus1 = c;
    auto dmins1 = d;

    std::string string_deri_center;
    if (deri_center == 0)
        string_deri_center = "alpha";
    else if (deri_center == 1)
        string_deri_center = "beta";
    else if (deri_center == 2)
        string_deri_center = "gamma";

    std::string s;
    s.append("    auto ");
    std::array<int,3> s_orbital {0,0,0};

    if (b == s_orbital and d == s_orbital)
    {
        int lla = (a[0]+a[1]+a[2]+1)*(a[0]+a[1]+a[2]+2)/2;
        int llb = (b[0]+b[1]+b[2]+1)*(b[0]+b[1]+b[2]+2)/2;
        int llc = (c[0]+c[1]+c[2]+1)*(c[0]+c[1]+c[2]+2)/2;
        int lld = (d[0]+d[1]+d[2]+1)*(d[0]+d[1]+d[2]+2)/2;
        s.append(namemap(a,b,c,d,m));
        s.append("_");
        s.append(string_deri_center);
        s.append("_con = 0.0;\n");
        return s;
    }
    // B center
    if (center == 'b')
    {
        aplus1[xyz] = a[xyz] + 1;
        bmins1[xyz] = b[xyz] - 1;

        s.append(namemap(a,b,c,d,m));
        s.append("_");
        s.append(string_deri_center);
        s.append("_con = ");
        s.append(namemap(aplus1,bmins1,c,d,m));
        s.append("_");
        s.append(string_deri_center);
        s.append("_con + ");
        s.append("AB[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(a,bmins1,c,d,m));
        s.append("_");
        s.append(string_deri_center);
        s.append("_con ; \n");
    }

    else if (center == 'd')
    {
        cplus1[xyz] = c[xyz] + 1;
        dmins1[xyz] = d[xyz] - 1;

        s.append(namemap(a,b,c,d,m));
        s.append("_");
        s.append(string_deri_center);
        s.append("_con = ");
        s.append(namemap(a,b,cplus1,dmins1,m));
        s.append("_");
        s.append(string_deri_center);
        s.append("_con + ");
        s.append("CD[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(a,b,c,dmins1,m));
        s.append("_");
        s.append(string_deri_center);
        s.append("_con ; \n");
    }
    return s;
}


void hgp(const std::array<int, 3>& a, const std::array<int, 3>& b,const std::array<int, 3>& c,const std::array<int, 3>& d, int m, std::map<std::array<int, 14>, std::string>& osmap, int deri_center)
// deri_center: the center where the derivation applies on, 0 means center A, 1 means center B, 2 means center C
{
    int xyz;
    std::array<int, 14> os_id = {a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], d[0], d[1], d[2], m, deri_center};
    if (osmap.count(os_id) != 0)
        return ;
    if ((b == std::array<int, 3> {0,0,0}) and
        (d == std::array<int, 3> {0,0,0}))
    {
        auto sen = sent_gen(a,b,c,d,m,'a',0, deri_center);
        // std::cout  << sen ;
        if (osmap.count(os_id) == 0)
        {
           osmap.insert({os_id, sen});
        }
        return ;
    }
    else if (!intex(a,b,c,d)) // if a,b,c,d don't exist
    {
        return ;
    }
    
    int lambda_b = b[0] + b[1] + b[2];
    int lambda_d = d[0] + d[1] + d[2];

    auto center = center_choose(lambda_b, lambda_d);
    // std::cout << center <<std::endl;
    auto aplus1 = a;
    auto bmins1 = b;
    auto cplus1 = c;
    auto dmins1 = d;
    //B center
    if (center == 'b')
    {
        auto xyz = dirchoose(b);
        aplus1[xyz] = a[xyz] + 1;
        bmins1[xyz] = b[xyz] - 1;

        auto sen = sent_gen(a,b,c,d,m,center,xyz, deri_center);
        std::cout << sen ;

        if (osmap.count(os_id) == 0)
        {
            osmap.insert({os_id, sen});
        }
        hgp(aplus1, bmins1, c, d, m, osmap, deri_center);
        hgp(a     , bmins1, c, d, m, osmap, deri_center);
    }
    else if (center == 'd')
    {
        auto xyz = dirchoose(d);
        cplus1[xyz] = c[xyz] + 1;
        dmins1[xyz] = d[xyz] - 1;

        auto sen = sent_gen(a,b,c,d,m,center,xyz,deri_center);
        std::cout << sen ;

        if (osmap.count(os_id) == 0)
        {
            osmap.insert({os_id, sen});
        }
        hgp(a, b, cplus1, dmins1, m, osmap, deri_center);
        hgp(a, b, c     , dmins1, m, osmap, deri_center);
    }
}


std::string deri_sen(std::array<int, 3 > a, std::array<int, 3 > b, std::array<int, 3 > c, std::array<int, 3 > d, char center, int xyz)
{
    std::string s;
    if (center == 'A')
    {
        auto aplus1 = a;
        auto amins1 = a;
        aplus1[xyz] = a[xyz] + 1;
        amins1[xyz] = a[xyz] - 1;
        s.append("    auto d_A");
        s.append(std::to_string(xyz));
        s.append("_");
        s.append(NewNameScheme(a));
        s.append(NewNameScheme(b));
        s.append(NewNameScheme(c));
        s.append(NewNameScheme(d));
        s.append("0_");
        s.append("con");
        s.append(" = ");
        s.append(NewNameScheme(aplus1));
        s.append(NewNameScheme(b));
        s.append(NewNameScheme(c));
        s.append(NewNameScheme(d));
        s.append("0_");
        s.append("alpha_con");
        if (intex(amins1, b, c, d))
        {
            s.append(" - ");
            s.append(std::to_string(a[xyz]));
            s.append(" * ");
            s.append(NewNameScheme(amins1));
            s.append(NewNameScheme(b));
            s.append(NewNameScheme(c));
            s.append(NewNameScheme(d));
            s.append("0_");
            s.append("con");
        }
        s.append(" ;\n");
    }
    // std::cout << s << std::endl;
    return s;
}

void deri_rr(std::array<int, 3 > a, std::array<int, 3 > b, std::array<int, 3 > c, std::array<int, 3 > d, int deri_center, int xyz, std::map<std::array<int, 14>, std::string>& scaled_map, std::map<std::array<int, 13>, std::string>& normal_map, std::map<std::array<int, 13>, std::string>& deri_map)
{

    auto aplus1 = a;
    auto amins1 = a;
    auto bplus1 = b;
    auto bmins1 = b;
    auto cplus1 = c;
    auto cmins1 = c;
    auto dplus1 = d;
    auto dmins1 = d;
        
    if (deri_center == 0)
    {
        aplus1[xyz] = a[xyz] + 1;
        amins1[xyz] = a[xyz] - 1;
        auto sen = deri_sen(a,  b,  c,  d, 'A', xyz);
        hgp(aplus1, b, c, d, 0, scaled_map, 0);
        if (intex(amins1, b, c, d))
        {
            hgp(amins1, b, c, d, 0, normal_map);
        }
        std::array<int, 13> deri_id = {a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], d[0], d[1], d[2], xyz};
        if (deri_map.count(deri_id) == 0 )
        {
            deri_map.insert({deri_id, sen});
        }
    }
}

void hrr_code_print(int la, int lb, int lc, int ld,std::map<std::array<int, 14>, std::string>& map)
{
    int deri_center = 0; // 0 means A center. 
    for (auto na = la + lb; na >= la; na--)
    for (auto nb = 0; nb <= lb; nb++)
    for (auto nc = lc + ld; nc >= lc; nc--)
    for (auto nd = 0; nd <= ld; nd++)
    {   
        for (auto ax = na; ax >= 0; ax--)
        for (auto ay = na - ax; ay >= 0; ay--)
        for (auto bx = nb; bx >= 0; bx--)
        for (auto by = nb - bx; by >= 0; by--)
        for (auto cx = nc; cx >= 0; cx--)
        for (auto cy = nc - cx; cy >= 0; cy--)
        for (auto dx = nd; dx >= 0; dx--)
        for (auto dy = nd - dx; dy >= 0; dy--)
        {   
            if (bx != 0 or by != 0 or nb-bx-by != 0 or dx != 0 or dy != 0 or nd-dx-dy != 0 ) 
            {   
            std::array<int, 14> os_id = {ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, dx, dy, nd-dx-dy, 0, 0};
            if (map.count(os_id) != 0)
                std::cout << map.at(os_id);
            }   
        }   

    }   
}


void hrr_code_print(int la, int lb, int lc, int ld,std::map<std::array<int, 13>, std::string>& map)
{
    for (auto na = la + lb; na >= la; na--)
    for (auto nb = 0; nb <= lb; nb++)
    for (auto nc = lc + ld; nc >= lc; nc--)
    for (auto nd = 0; nd <= ld; nd++)
    {   
        for (auto ax = na; ax >= 0; ax--)
        for (auto ay = na - ax; ay >= 0; ay--)
        for (auto bx = nb; bx >= 0; bx--)
        for (auto by = nb - bx; by >= 0; by--)
        for (auto cx = nc; cx >= 0; cx--)
        for (auto cy = nc - cx; cy >= 0; cy--)
        for (auto dx = nd; dx >= 0; dx--)
        for (auto dy = nd - dx; dy >= 0; dy--)
        {   
            if (bx != 0 or by != 0 or nb-bx-by != 0 or dx != 0 or dy != 0 or nd-dx-dy != 0 ) 
            {   
            std::array<int, 13> os_id = {ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, dx, dy, nd-dx-dy, 0};
            // printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, dx, dy, nd-dx-dy);
            if (map.count(os_id) != 0)
                std::cout << map.at(os_id);
            }   
        }   

    }   
}


static void vrr(std::map<std::array<int, 14>, std::string>& hrr_scaled_map, 
         std::map<std::array<int, 13>, std::string>& hrr_normal_map,
         std::map<std::array<int, 13>, std::string>& os_map)
{
    for(auto it = hrr_scaled_map.cbegin(); it != hrr_scaled_map.cend(); ++it)
    {
        auto k = it->first;
        std::array<int, 3> a = {k[0], k[1], k[2]};
        std::array<int, 3> b = {k[3], k[4], k[5]};
        std::array<int, 3> c = {k[6], k[7], k[8]};
        std::array<int, 3> d = {k[9], k[10], k[11]};
        std::array<int, 3> s = {0, 0, 0};
        if (b == s and d == s)
        {
            obara_saika(a, b, c, d, 0, os_map);   
        }
    }
    for(auto it = hrr_normal_map.cbegin(); it != hrr_normal_map.cend(); ++it)
    {
        auto k = it->first;
        std::array<int, 3> a = {k[0], k[1], k[2]};
        std::array<int, 3> b = {k[3], k[4], k[5]};
        std::array<int, 3> c = {k[6], k[7], k[8]};
        std::array<int, 3> d = {k[9], k[10], k[11]};
        std::array<int, 3> s = {0, 0, 0};
        if (b == s and d == s)
        {
            obara_saika(a, b, c, d, 0, os_map);   
        }
    }
}


void di_hrr(int la, int lb, int lc, int ld, int deri_center,
            std::map<std::array<int, 14>, std::string>& scaled_map,
            std::map<std::array<int, 13>, std::string>& normal_map,
            std::map<std::array<int, 13>, std::string>& os_map,
            std::map<std::array<int, 13>, std::string>& deri_map)
{
    if (deri_center == 0)
    {
        for (auto ax = la; ax >= 0; ax--)
        for (auto ay = la - ax; ay >= 0; ay--)
        for (auto bx = lb; bx >= 0; bx--)
        for (auto by = lb - bx; by >= 0; by--)
        for (auto cx = lc; cx >= 0; cx--)
        for (auto cy = lc - cx; cy >= 0; cy--)
        for (auto dx = ld; dx >= 0; dx--)
        for (auto dy = ld - dx; dy >= 0; dy--)
        {
            std::array<int, 3> a = {ax, ay, la - ax - ay};
            std::array<int, 3> b = {bx, by, lb - bx - by};
            std::array<int, 3> c = {cx, cy, lc - cx - cy};
            std::array<int, 3> d = {dx, dy, ld - dx - dy};
            for (auto xyz = 0; xyz <= 2; xyz++)
            {
            deri_rr(a,b,c,d, 0, xyz, scaled_map, normal_map, deri_map);
            }
            hgp(a,b,c,d,0,normal_map);
        }
        
    }
    // std::cout << "scaled_map size: " << scaled_map.size() << std::endl;
    // std::cout << "normal_map size: " << normal_map.size() << std::endl;
//    std::cout << "os_map size: " << os_map.size() << std::endl;
//get the integral first, including the normal map (hrr) and os map (vrr)
    for (auto l = la +lb; l >= la; l--)
    for (auto r = lc +ld; r >= lc; r--)
    {   
        eri(l,0,r,0,os_map);
    }   
    // std::cout << "os_map size: " << os_map.size() << std::endl;

    // generate the vrr for scaled map and normal_map 
    vrr(scaled_map, normal_map, os_map);
    // std::cout << "scaled_map size: " << scaled_map.size() << std::endl;
    // std::cout << "normal_map size: " << normal_map.size() << std::endl;
    // std::cout << "os_map size: " << os_map.size() << std::endl;
}

void save_deri(int la, int lb, int lc, int ld, int deri_center)
{
    char center;
    if (deri_center == 0)
        center = 'A';
    else if (deri_center == 1)
        center = 'B';
    else if (deri_center == 2)
        center = 'C';
    auto lla = (la+1)*(la+2)/2;
    auto llb = (lb+1)*(lb+2)/2;
    auto llc = (lc+1)*(lc+2)/2;
    auto lld = (ld+1)*(ld+2)/2;
    auto llabcd = lla*llb*llc*lld;
    auto xyz = 3;
    for (auto ixyz = 0; ixyz < xyz; ixyz++)
    {
        auto id = 0;
        // printf("ld = %d\n", ld);
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
                        //  cab + nab * (ia + lla * (ib + llb * (ic + llc * (id + lld * (ccd + ncd * (xyz + 3 * center)))))) 
                        //
                        printf("    dI_[idx + nab * (%d + %d * (%d + %d * (%d + %d * (%d + %d * (idy + ncd * (%d + %d * %d))))))]", ia, lla, ib, llb, ic, llc, id, lld, ixyz, xyz, deri_center);
                        std::cout << " = "; 
                        std::cout << "d_" << center << ixyz << "_"<< namemap(axyz,bxyz,cxyz,dxyz,0) << "_con ;"<< std::endl;
                        ia++;
                    }
                    ib++;
                }
                ic++;
            }
            id++;
            // printf("id = %d\n ", id);
            // printf("dx, dy, dz = %d, %d, %d\n ", dx, dy, dz);
        }
    }
    std::cout << "    }" << std::endl;
}

void d_code_print(int la, int lb, int lc, int ld,
                  std::map<std::array<int, 14>, std::string>& scaled_map,
                  std::map<std::array<int, 13>, std::string>& normal_map,
                  std::map<std::array<int, 13>, std::string>& os_map,
                  std::map<std::array<int, 13>, std::string>& deri_map)
{
    // std::cout << HrrFuncName(la,lb,lc,ld) << std::endl;
    std::cout << "KERNEL void di_ssss(shell_pairs_gpu& ab, shell_pairs_gpu&cd, double* I_, double* dI_)\n{\n";
    std::cout << "    auto idx = threadIdx.x + blockDim.x * blockIdx.x; // idx for shell_pairs ab \n\
    auto idy = threadIdx.y + blockDim.y * blockIdx.y; // idy for shell_pairs cd \n\n";
    
    std::cout << "    auto nab = ab.na * ab.nb; //  sizeof(ab.S)/sizeof(ab.S[0]);\n";
    std::cout << "    auto ncd = cd.na * cd.nb; //  sizeof(cd.S)/sizeof(cd.S[0]);\n\n";
    std::cout << "    auto mab = ab.ma * ab.mb; \n";
    std::cout << "    auto mcd = cd.ma * cd.mb; \n\n";
    std::cout << "    if (idx < nab  and idy < ncd)\n    {\n";

    std::cout << "    double AB[3] = \n    {\n        ab.AB[0][idx],\n        ab.AB[1][idx],\n        ab.AB[2][idx]\n    };\n";
    std::cout << "    double CD[3] = \n    {\n        cd.AB[0][idx],\n        cd.AB[1][idx],\n        cd.AB[2][idx]\n    };\n\n";
    for (auto na = la + lb + 1; na >= la; na--)
    for (auto nb = 0; nb <= lb; nb++)
    for (auto nc = lc + ld; nc >= lc; nc--)
    for (auto nd = 0; nd <= ld; nd++)
    {   
        for (auto ax = na; ax >= 0; ax--)
        for (auto ay = na - ax; ay >= 0; ay--)
        for (auto bx = nb; bx >= 0; bx--)
        for (auto by = nb - bx; by >= 0; by--)
        for (auto cx = nc; cx >= 0; cx--)
        for (auto cy = nc - cx; cy >= 0; cy--)
        for (auto dx = nd; dx >= 0; dx--)
        for (auto dy = nd - dx; dy >= 0; dy--)
        {   
            // std::array<int, 13> os_id = {ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, dx, dy, nd-dx-dy, 0};
            if (bx == 0 and by == 0 and nb-bx-by == 0 and dx == 0 and dy == 0 and nd-dx-dy == 0)
            {   
            std::array<int, 14> scaled_id = {ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, dx, dy, nd-dx-dy, 0, 0}; 
            if (scaled_map.count(scaled_id) != 0)
                std::cout << scaled_map.at(scaled_id);
            }   
        }   
    }


    for (auto na = la + lb; na >= la-1; na--)
    for (auto nb = 0; nb <= lb; nb++)
    for (auto nc = lc + ld; nc >= lc-1; nc--)
    for (auto nd = 0; nd <= ld; nd++)
    {   
        for (auto ax = na; ax >= 0; ax--)
        for (auto ay = na - ax; ay >= 0; ay--)
        for (auto bx = nb; bx >= 0; bx--)
        for (auto by = nb - bx; by >= 0; by--)
        for (auto cx = nc; cx >= 0; cx--)
        for (auto cy = nc - cx; cy >= 0; cy--)
        for (auto dx = nd; dx >= 0; dx--)
        for (auto dy = nd - dx; dy >= 0; dy--)
        {   
            // std::array<int, 13> os_id = {ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, dx, dy, nd-dx-dy, 0};
            if (bx == 0 and by == 0 and nb-bx-by == 0 and dx == 0 and dy == 0 and nd-dx-dy == 0)
            {   
            std::array<int, 13> normal_id = {ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, dx, dy, nd-dx-dy, 0}; 
            if (normal_map.count(normal_id) != 0)
                std::cout << normal_map.at(normal_id);
            }   
        }   
    }

    std::cout << "\n    for (int pab = 0; pab < mab; pab++)" << std::endl;
    std::cout << "    for (int pcd = 0; pcd < mcd; pcd++)" << std::endl;
    std::cout << "\
    {\n\
        auto iab = idx + pab * nab; \n\
        auto icd = idy + pcd * ncd; \n\
                                    \n\
        auto zab = ab.Z[iab];\n\
        auto zcd = cd.Z[icd];\n\
                              \n\
        auto za = ab.ZA[iab];\n\
        auto zb = zab - za;\n\
        auto zc = cd.ZA[icd];\n\
                             \n\ 
        auto zab_inv = 1 / zab;\n\
        auto zcd_inv = 1 / zcd;\n\ 
        auto zabcd_inv_sqrt = rsqrt(zab + zcd);\n\ 
        auto zabcd_inv = zabcd_inv_sqrt * zabcd_inv_sqrt;\n\ 
        auto rho = zab * zcd * zabcd_inv;\n\
                                         \n\
        double P[3] = \n\
        {\n\
            ab.P[0][iab],\n\
            ab.P[1][iab],\n\
            ab.P[2][iab],\n\
        };\n\
        \n\
        double Q[3] = \n\
        {\n\
            cd.P[0][icd],\n\
            cd.P[1][icd],\n\
            cd.P[2][icd],\n\
        };\n\
        \n\
        double PA[3] =\n\
        {\n\
            ab.PA[0][iab],\n\
            ab.PA[1][iab],\n\
            ab.PA[2][iab],\n\
        };\n\
        \n\
        double PB[3] =\n\
        {\n\
            ab.PA[0][iab] + ab.AB[0][iab],\n\
            ab.PA[1][iab] + ab.AB[1][iab],\n\
            ab.PA[2][iab] + ab.AB[2][iab],\n\
        };\n\
        \n\
        double QC[3] = \n\
        {\n\
            cd.PA[0][icd],\n\
            cd.PA[1][icd],\n\
            cd.PA[2][icd],\n\
        };\n\
        \n\
        double QD[3] =\n\
        {\n\
            cd.PA[0][icd] + cd.AB[0][icd],\n\
            cd.PA[1][icd] + cd.AB[1][icd],\n\
            cd.PA[2][icd] + cd.AB[2][icd],\n\
        };\n\
        \n\
        double W[3] =\n\
        {\n\
            (zab * P[0] + zcd * Q[0]) * zabcd_inv,\n\
            (zab * P[1] + zcd * Q[1]) * zabcd_inv,\n\
            (zab * P[2] + zcd * Q[2]) * zabcd_inv,\n\
        };\n\
        double WP[3] =\n\
        {\n\
            W[0] - P[0],\n\
            W[1] - P[1],\n\
            W[2] - P[2]\n\
        };\n\
        \n\
        double WQ[3] =\n\
        {\n\
            W[0] - Q[0],\n\
            W[1] - Q[1],\n\
            W[2] - Q[2]\n\
        };\n\
        \n\
        auto T = rho * ((P[0]-Q[0])*(P[0]-Q[0]) + (P[1]-Q[1])*(P[1]-Q[1]) + (P[2]-Q[2])*(P[2]-Q[2]));\n\
        auto K = zabcd_inv_sqrt * (ab.K[iab]) * cd.K[icd]; \n\
        \n\    
        double fm[" << la+lb+lc+ld+2 << "]; \n\
        vgamma<double>(" << la+lb+lc+ld+1 << ", T, fm);\n\
        \n\
        for (auto i = 0; i < " << la+lb+lc+ld+2 << 
        " ; i++) \n\            fm[i] *= K;\n" << std::endl;
        vrr_code_print(la+lb+1,0,lc+ld,0,os_map);
    // print vrr (os map )
     // vrr_code_print(la+lb+1,0,lc+ld,0,os_map);

    // sumup the promitive functions integrals
    // normal part
    for (auto l = la +lb; l >= la-1; l--)
    for (auto r = lc +ld; r >= lc; r--)
    {   
        for (auto ax = l; ax >= 0; ax--)
        for (auto ay = l - ax; ay >= 0; ay--)
        for (auto cx = r; cx >= 0; cx--)
        for (auto cy = r - cx; cy >= 0; cy--)
        {   
            std::array<int, 3> a = {ax, ay, l - ax - ay};
            std::array<int, 3> b = {0, 0, 0}; 
            std::array<int, 3> c = {cx, cy, r - cx - cy};
            std::array<int, 3> d = {0, 0, 0}; 
            std::array<int, 13> normal_id = {ax, ay, l-ax-ay, 0, 0, 0, cx, cy, r-cx-cy, 0, 0, 0, 0};
            if (normal_map.count(normal_id) != 0)
            {
                std::cout << "        " << namemap(a,b,c,d,0) << "_con += " << namemap(a,b,c,d,0) << ";" << std::endl;
            }
        }   
    }

    // scaled part
    for (auto l = la +lb+1; l >= la; l--)
    for (auto r = lc +ld; r >= lc; r--)
    {   
        for (auto ax = l; ax >= 0; ax--)
        for (auto ay = l - ax; ay >= 0; ay--)
        for (auto cx = r; cx >= 0; cx--)
        for (auto cy = r - cx; cy >= 0; cy--)
        {   
            std::array<int, 3> a = {ax, ay, l - ax - ay};
            std::array<int, 3> b = {0, 0, 0}; 
            std::array<int, 3> c = {cx, cy, r - cx - cy};
            std::array<int, 3> d = {0, 0, 0}; 
            std::array<int, 14> scaled_id = {ax, ay, l-ax-ay, 0, 0, 0, cx, cy, r-cx-cy, 0, 0, 0, 0, 0};
            if (scaled_map.count(scaled_id) != 0)
            {
                std::cout << "        " <<  namemap(a,b,c,d,0) << "_alpha_con += " << "2 * za * " << namemap(a,b,c,d,0) << ";" << std::endl;
            }
        }   
    }    

        std::cout << "        }" << std::endl;

        std::cout << "         //****************//" << std::endl;
    // print hrr part (normal map)

    for (auto na = la + lb; na >= la; na--)
    for (auto nb = 0; nb <= lb; nb++)
    for (auto nc = lc + ld; nc >= lc; nc--)
    for (auto nd = 0; nd <= ld; nd++)
    {   
        for (auto ax = na; ax >= 0; ax--)
        for (auto ay = na - ax; ay >= 0; ay--)
        for (auto bx = nb; bx >= 0; bx--)
        for (auto by = nb - bx; by >= 0; by--)
        for (auto cx = nc; cx >= 0; cx--)
        for (auto cy = nc - cx; cy >= 0; cy--)
        for (auto dx = nd; dx >= 0; dx--)
        for (auto dy = nd - dx; dy >= 0; dy--)
        {   
            if (bx != 0 or by != 0 or nb-bx-by != 0 or dx != 0 or dy != 0 or nd-dx-dy != 0 ) 
            {   
            std::array<int, 13> normal_id = {ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, dx, dy, nd-dx-dy, 0}; 
            if (normal_map.count(normal_id) != 0)
                std::cout << normal_map.at(normal_id);
            }   
        }   

    }

    // print scaled hrr part (scaled map)
    for (auto na = la + lb + 1; na >= la; na--)
    for (auto nb = 0; nb <= lb; nb++)
    for (auto nc = lc + ld; nc >= lc; nc--)
    for (auto nd = 0; nd <= ld; nd++)
    {   
        for (auto ax = na; ax >= 0; ax--)
        for (auto ay = na - ax; ay >= 0; ay--)
        for (auto bx = nb; bx >= 0; bx--)
        for (auto by = nb - bx; by >= 0; by--)
        for (auto cx = nc; cx >= 0; cx--)
        for (auto cy = nc - cx; cy >= 0; cy--)
        for (auto dx = nd; dx >= 0; dx--)
        for (auto dy = nd - dx; dy >= 0; dy--)
        {   
            if (bx != 0 or by != 0 or nb-bx-by != 0 or dx != 0 or dy != 0 or nd-dx-dy != 0 ) 
            {   
            std::array<int, 14> scaled_id = {ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, dx, dy, nd-dx-dy, 0, 0}; 
            if (scaled_map.count(scaled_id) != 0)
                std::cout << scaled_map.at(scaled_id);
            }   
        }   
    }

    std::cout << "    // here is the contracted integral part" << std::endl;

    // print the integral
    for (auto ax = la; ax >= 0; ax--)
    for (auto ay = la - ax; ay >= 0; ay--)
    for (auto bx = lb; bx >= 0; bx--)
    for (auto by = lb - bx; by >= 0; by--)
    for (auto cx = lc; cx >= 0; cx--)
    for (auto cy = lc - cx; cy >= 0; cy--)
    for (auto dx = ld; dx >= 0; dx--)
    for (auto dy = ld - dx; dy >= 0; dy--)
    {
        if (bx != 0 or by != 0 or lb-bx-by != 0 or dx != 0 or dy != 0 or ld-dx-dy != 0 )
        {
        std::array<int, 13> normal_id = {ax, ay, la - ax - ay, bx, by, lb - bx - by, cx, cy, lc - cx - cy, dx, dy, ld - dx - dy, 0};
        if (normal_map.count(normal_id) != 0 )
            std::cout << normal_map.at(normal_id);
        }
    }

    // print the derivative equation
    for (auto xyz = 0; xyz <= 2; xyz++)
    for (auto ax = la; ax >= 0; ax--)
    for (auto ay = la - ax; ay >= 0; ay--)
    for (auto bx = lb; bx >= 0; bx--)
    for (auto by = lb - bx; by >= 0; by--)
    for (auto cx = lc; cx >= 0; cx--)
    for (auto cy = lc - cx; cy >= 0; cy--)
    for (auto dx = ld; dx >= 0; dx--)
    for (auto dy = ld - dx; dy >= 0; dy--)
    {
        std::array<int, 13> deri_id = {ax, ay, la - ax - ay, bx, by, lb - bx - by, cx, cy, lc - cx - cy, dx, dy, ld - dx - dy, xyz};
        if (deri_map.count(deri_id) != 0 )
            std::cout << deri_map.at(deri_id);
    }


    // save the integral
    save_int(la,lb,lc,ld);

    // save the derivative of integral
    save_deri(la,lb,lc,ld,0);

    std::cout << "}" << std::endl;
}



int main()
{
    std::array<int, 3> a = {1, 0, 0};
    std::array<int, 3> b = {0, 0, 0};
    std::array<int, 3> c = {0, 0, 0};
    std::array<int, 3> d = {0, 0, 0};
    // auto sen = deri_sen(a,b,c,d,'A',0);
    // std::cout << sen << std::endl;
    std::map<std::array<int, 14>, std::string> scaled_map;
    std::map<std::array<int, 13>, std::string> normal_map;
    std::map<std::array<int, 13>, std::string> os_map;
    std::map<std::array<int, 13>, std::string> deri_map;
    // deri_rr(a,b,c,d, 0, 0, scaled_map, normal_map);
    auto la = a[0] + a[1] + a[2]; 
    auto lb = b[0] + b[1] + b[2]; 
    auto lc = c[0] + c[1] + c[2]; 
    auto ld = d[0] + d[1] + d[2]; 
    // std::cout << la << "," << lb << ", " << lc << ", " << ld << std::endl;
    // hrr_code_print(la+1, lb, lc, ld, scaled_map);
    // hrr_code_print(la-1, lb, lc, ld, normal_map);
    // vrr(scaled_map, normal_map, os_map);
//    std::cout << "scaled_map size: " << scaled_map.size() << std::endl;
//    std::cout << "normal_map size: " << normal_map.size() << std::endl;
//    std::cout << "os_map size: " << os_map.size() << std::endl;
    di_hrr(la, lb, lc, ld, 0, scaled_map, normal_map, os_map, deri_map);
    // std::cout << "deri_map size: " << deri_map.size() << std::endl;
    d_code_print(la,lb,lc,ld,scaled_map, normal_map, os_map, deri_map);
// for(auto it =os_map.cbegin(); it != os_map.cend(); ++it)
// {
//     std::cout << it->first[0]  << "\n";
// }

}
    


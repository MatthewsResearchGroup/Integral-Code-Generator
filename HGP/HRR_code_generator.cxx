#include <iostream>
#include <map>
#include "HRR_code_generator.h"

std::string sent_gen(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c,const std::array<int, 3> d, int m, char center, int xyz)
{
    auto aplus1 = a;
    auto bmins1 = b;
    auto cplus1 = c;
    auto dmins1 = d;

    std::string s;
    s.append("    double ");
    std::array<int,3> s_orbital {0,0,0};

    if (b == s_orbital and d == s_orbital)
    {
        int lla = (a[0]+a[1]+a[2]+1)*(a[0]+a[1]+a[2]+2)/2; 
        int llb = (b[0]+b[1]+b[2]+1)*(b[0]+b[1]+b[2]+2)/2; 
        int llc = (c[0]+c[1]+c[2]+1)*(c[0]+c[1]+c[2]+2)/2; 
        int lld = (d[0]+d[1]+d[2]+1)*(d[0]+d[1]+d[2]+2)/2; 
        s.append(namemap(a,b,c,d,m));
        s.append(" = ");
        s.append("I_");
        s.append(orbname(a));
        s.append(orbname(b));
        s.append(orbname(c));
        s.append(orbname(d));
        s.append("[");
        // cab + (ia+ib*lla+ic*lla*llb+id*lla*llb*llc)*nab + ccd*lla*llb*llc*lld*nab
        // s.append(std::to_string(m));
        s.append("idx + (");
        s.append(std::to_string(addrsear(a)));
        s.append(" + ");
        s.append(std::to_string(addrsear(b)));
        s.append(" * ");
        s.append(std::to_string(lla));
        s.append(" + ");
        s.append(std::to_string(addrsear(c)));
        s.append(" * ");
        s.append(std::to_string(lla));
        s.append(" * ");
        s.append(std::to_string(llb));
        s.append(" + ");
        s.append(std::to_string(addrsear(d)));
        s.append(" * ");
        s.append(std::to_string(lla));
        s.append(" * ");
        s.append(std::to_string(llb));
        s.append(" * ");
        s.append(std::to_string(llc));
        s.append(") * nab");
        s.append(" + idy *");
        s.append(std::to_string(lla * llb * llc *lld));
        s.append(" * nab");
        s.append("] ;\n");
        return s;
    }
    // B center
    if (center == 'b')
    {
        aplus1[xyz] = a[xyz] + 1;
        bmins1[xyz] = b[xyz] - 1;
        
        s.append(namemap(a,b,c,d,m));
        s.append(" = ");
        s.append(namemap(aplus1,bmins1,c,d,m));
        s.append(" + ");
        s.append("AB[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(a,bmins1,c,d,m));
        s.append(" ; \n");
    }

    else if (center == 'd')
    {
        cplus1[xyz] = c[xyz] + 1;
        dmins1[xyz] = d[xyz] - 1;

        s.append(namemap(a,b,c,d,m));
        s.append(" = ");
        s.append(namemap(a,b,cplus1,dmins1,m));
        s.append(" + ");
        s.append("CD[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(a,b,c,dmins1,m));
        s.append(" ; \n");
    }
    return s;
}

void hgp(const std::array<int, 3>& a, const std::array<int, 3>& b,const std::array<int, 3>& c,const std::array<int, 3>& d, int m, std::map<std::array<int, 13>, std::string>& osmap)
{
    int xyz;
    std::array<int, 13> os_id = {a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], d[0], d[1], d[2], m};
    if (osmap.count(os_id) != 0)
        return ;
    if ((b == std::array<int, 3> {0,0,0}) and
        (d == std::array<int, 3> {0,0,0}))
    {
        auto sen = sent_gen(a,b,c,d,m,'a',0);
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

        auto sen = sent_gen(a,b,c,d,m,center,xyz);
        // std::cout << sen << std::endl;

        if (osmap.count(os_id) == 0)
        {
            osmap.insert({os_id, sen});
        }
        hgp(aplus1, bmins1, c, d, m, osmap);
        hgp(a     , bmins1, c, d, m, osmap);
    }
    else if (center == 'd')
    {
        auto xyz = dirchoose(d);
        cplus1[xyz] = c[xyz] + 1;
        dmins1[xyz] = d[xyz] - 1;

        auto sen = sent_gen(a,b,c,d,m,center,xyz);
        // std::cout << sen << std::endl;

        if (osmap.count(os_id) == 0)
        {
            osmap.insert({os_id, sen});
        }
        hgp(a, b, cplus1, dmins1, m, osmap);
        hgp(a, b, c     , dmins1, m, osmap);
    }
}

void hgp_eri(int la, int lb, int lc, int ld, std::map<std::array<int, 13>, std::string>& osmap)
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
                    hgp(axyz,bxyz,cxyz,dxyz,0, osmap);
                }
            }
        }
    }
}

std::string HrrFuncName(int la, int lb, int lc, int ld)
{
    std::string s;
    s.append("KERNEL void hrr_");
    s.append(orbname(std::array<int, 3> {la,0,0}));
    s.append(orbname(std::array<int, 3> {lb,0,0}));
    s.append(orbname(std::array<int, 3> {lc,0,0}));
    s.append(orbname(std::array<int, 3> {ld,0,0}));
    s.append("(shell_pairs_gpu& ab, shell_pairs_gpu& cd, int nab, double *I_, ");
    for (auto l = la +lb; l >= la; l--)
    for (auto r = lc +ld; r >= lc; r--)
    {
        s.append("double *I_");
        s.append(orbname(std::array<int, 3> {l,0,0}));
        s.append("s");
        s.append(orbname(std::array<int, 3> {r,0,0}));
        s.append("s");
        s.append(", ");
    }
    s.append(")\n{");

    return s;
}


void code_print(int la , int lb, int lc, int ld, std::map<std::array<int, 13>, std::string>& osmap)
{
    std::cout << HrrFuncName(la,lb,lc,ld) << std::endl;
    
    std::cout << "    auto idx = threadIdx.x + blockDim.x * blockIdx.x; // idx for shell_pairs ab \n\
    auto idy = threadIdx.y + blockDim.y * blockIdx.y; // idy for shell_pairs cd \n\n";
    
    std::cout << "    double AB[3] = \n    {\n        ab.AB[0][idx],\n        ab.AB[1][idx],\n        ab.AB[2][idx]\n    };\n";
    std::cout << "    double CD[3] = \n    {\n        cd.AB[0][idx],\n        cd.AB[1][idx],\n        cd.AB[2][idx]\n    };\n\n";

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
            std::array<int, 13> os_id = {ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, dx, dy, nd-dx-dy, 0};
            if (osmap.count(os_id) != 0)
                std::cout << osmap.at(os_id);
        }
    }
    save_int(la, lb, lc, ld);
    std::cout << "}" << std::endl;
}


int main()
{
    std::array<int, 3> c = {1,2,0};
    std::array<int, 3> d = {0,0,0};
    std::array<int, 3> a = {1,0,0};
    std::array<int, 3> b = {0,0,0};
    // auto h = addrsear(c);
    // printf("%d\n",h);
    // auto a_name = NewNameScheme(a);
    // std::cout << "a_name = " << a_name << std::endl;
    // auto hello = sent_gen(a,b,c,d,0,'d',1);
    // std::cout << hello << std::endl;

    int la = 1; 
    int lb = 1; 
    int lc = 1; 
    int ld = 1; 

    //auto funcname = HrrFuncName(la,lb,lc,ld);
    //std::cout << funcname << std::endl;
    std::map<std::array<int, 13>, std::string> osmap;
    hgp_eri(la, lb, lc, ld, osmap);
    code_print(la, lb, lc, ld, osmap);
    // save_int(la, lb, lc, ld);
    // std::cout << hello << std::endl;
    // printf("Hello\n");
}

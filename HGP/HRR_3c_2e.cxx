#include <HRR_3c_2e.h>

/* 
 * Since there is no center D in 3 center 2 electron integral
 * then we only need to apply HGP algorithm to center B
 */

std::string sent_gen(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c, int m, int xyz)
{
    auto aplus1 = a;
    auto bmins1 = b;
    auto cplus1 = c;

    std::string s;
    s.append("    auto ");
    std::array<int,3> s_orbital {0,0,0};

    if (b == s_orbital)
    {
        s.append(namemap(a,b,c,m));
        s.append("_con = 0.0;\n");
        return s;
    }
    // B center
    aplus1[xyz] = a[xyz] + 1;
    bmins1[xyz] = b[xyz] - 1;

    s += namemap(a,b,c,m) + "_con = " + namemap(aplus1,bmins1,c,m) + "_con + " + "AB[" + std::to_string(xyz) + "] * " + namemap(a,bmins1,c,m) + "_con ; \n";

    return s;
}


void hgp(const std::array<int, 3>& a, const std::array<int, 3>& b,const std::array<int, 3>& c, int m, std::map<std::array<int, 10>, std::string>& osmap)
{
    // int xyz;
    std::array<int, 10> os_id = {a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], m};
    if (osmap.count(os_id) != 0)
        return ;
    if (b == std::array<int, 3> {0,0,0})
    {
        auto sen = sent_gen(a,b,c,m,0);
        if (osmap.count(os_id) == 0)
        {
           osmap.insert({os_id, sen});
        }
        return ;
    }
    else if (!intex(a,b,c)) // if a,b,c,d don't exist
    {
        return ;
    }
    
    auto aplus1 = a;
    auto bmins1 = b;
    auto xyz = dirchoose(b);
    aplus1[xyz] = a[xyz] + 1;
    bmins1[xyz] = b[xyz] - 1;

    auto sen = sent_gen(a,b,c,m,xyz);
    // std::cout << sen << std::endl;

    if (osmap.count(os_id) == 0)
    {
        osmap.insert({os_id, sen});
    }
    hgp(aplus1, bmins1, c, m, osmap);
    hgp(a     , bmins1, c, m, osmap);
}


void hgp_eri(int la, int lb, int lc, std::map<std::array<int, 10>, std::string>& osmap)
{
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
                hgp(axyz, bxyz, cxyz, 0, osmap);
            }
        }
    }
}


void code_print(int la , int lb, int lc, int ld, std::map<std::array<int, 10>, std::string>& osmap)
{
    std::string s = "KERNEL void hrr_" + orbname(std::array<int, 3> {la,0,0}) + orbname(std::array<int, 3> {lb,0,0}) + orbname(std::array<int, 3> {lc,0,0}) + "(shell_pairs_gpu ab, shell_pairs_gpu cd,  double* I_)\n{";
    std::cout << s << std::endl;
    
    std::cout << "    auto idx = threadIdx.x + blockDim.x * blockIdx.x; // idx for shell_pairs ab \n\
    auto idy = threadIdx.y + blockDim.y * blockIdx.y; // idy for shell_pairs cd \n\n";
    
    std::cout << "    auto nab = ab.na * ab.nb; //  sizeof(ab.S)/sizeof(ab.S[0]);\n";
    std::cout << "    auto ncd = cd.na * cd.nb; //  sizeof(cd.S)/sizeof(cd.S[0]);\n\n";
    std::cout << "    auto mab = ab.ma * ab.mb; \n";
    std::cout << "    auto mcd = cd.ma * cd.mb; \n\n";
    std::cout << "    if (idx < nab  and idy < ncd)\n    {\n";

    std::cout << "    double AB[3] = \n    {\n        ab.AB[0][idx],\n        ab.AB[1][idx],\n        ab.AB[2][idx]\n    };\n";
    std::cout << "    double CD[3] = \n    {\n        cd.AB[0][idy],\n        cd.AB[1][idy],\n        cd.AB[2][idy]\n    };\n\n";
    // initial the contracted integral
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
            // std::array<int, 13> os_id = {ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, dx, dy, nd-dx-dy, 0};
            if (bx == 0 and by == 0 and nb-bx-by == 0)
            {
            std::array<int, 10> os_id = {ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, 0};
            if (osmap.count(os_id) != 0)
                std::cout << osmap.at(os_id);
            }
        }
    }

    // vrr code here (maybe a function name)
    
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
        double fm[" << la+lb+lc+1 << "]; \n\
        vgamma<double>(" << la+lb+lc << ", T, fm);\n\
        \n\
        for (auto i = 0; i <" << la+lb+lc+1 << 
        " ; i++) \n\            fm[i] *= K;\n" << std::endl;
        // merge_vrr_code( la,  lb,  lc,  ld);
        std::cout << "        }" << std::endl;

    std::cout << "         //****************//" << std::endl;


    // cotracted integrals equation

    for (auto na = la + lb; na >= la; na--)
    for (auto nb = 0; nb <= lb; nb++)
    for (auto nc = lc; nc >= lc; nc--)
    {
        for (auto ax = na; ax >= 0; ax--)
        for (auto ay = na - ax; ay >= 0; ay--)
        for (auto bx = nb; bx >= 0; bx--)
        for (auto by = nb - bx; by >= 0; by--)
        for (auto cx = nc; cx >= 0; cx--)
        for (auto cy = nc - cx; cy >= 0; cy--)
        {
            if (bx != 0 or by != 0 or nb-bx-by != 0)
            {
            std::array<int, 10> os_id = {ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, 0};
            if (osmap.count(os_id) != 0)
                std::cout << osmap.at(os_id);
            }
        }   

    }

    // save_int(la, lb, lc, ld);
    std::cout << "    }\n}" << std::endl;
}

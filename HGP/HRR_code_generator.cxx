#include  <array>
#include <iostream>
#include "HRR_code_generator.h"

void merge_vrr_code(int la, int lb, int lc, int  ld);

std::string sent_gen(const std::array<int, 3> a, const std::array<int, 3> b,const std::array<int, 3> c,const std::array<int, 3> d, int m, char center, int xyz)
{
    auto aplus1 = a;
    auto bmins1 = b;
    auto cplus1 = c;
    auto dmins1 = d;

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
        s.append("_con = 0.0;\n");
        // s.append("I_");
        // s.append(orbname(a));
        // s.append(orbname(b));
        // s.append(orbname(c));
        // s.append(orbname(d));
        // s.append("[");
        // // cab + (ia+ib*lla+ic*lla*llb+id*lla*llb*llc)*nab + ccd*lla*llb*llc*lld*nab
        // // s.append(std::to_string(m));
        // s.append("idx + (");
        // s.append(std::to_string(addrsear(a)));
        // s.append(" + ");
        // s.append(std::to_string(addrsear(b)));
        // s.append(" * ");
        // s.append(std::to_string(lla));
        // s.append(" + ");
        // s.append(std::to_string(addrsear(c)));
        // s.append(" * ");
        // s.append(std::to_string(lla));
        // s.append(" * ");
        // s.append(std::to_string(llb));
        // s.append(" + ");
        // s.append(std::to_string(addrsear(d)));
        // s.append(" * ");
        // s.append(std::to_string(lla));
        // s.append(" * ");
        // s.append(std::to_string(llb));
        // s.append(" * ");
        // s.append(std::to_string(llc));
        // s.append(") * nab");
        // s.append(" + idy *");
        // s.append(std::to_string(lla * llb * llc *lld));
        // s.append(" * nab");
        // s.append("] ;\n");
        return s;
    }
    // B center
    if (center == 'b')
    {
        aplus1[xyz] = a[xyz] + 1;
        bmins1[xyz] = b[xyz] - 1;
        
        s.append(namemap(a,b,c,d,m));
        s.append("_con = ");
        s.append(namemap(aplus1,bmins1,c,d,m));
        s.append("_con + ");
        s.append("AB[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(a,bmins1,c,d,m));
        s.append("_con ; \n");
    }

    else if (center == 'd')
    {
        cplus1[xyz] = c[xyz] + 1;
        dmins1[xyz] = d[xyz] - 1;

        s.append(namemap(a,b,c,d,m));
        s.append("_con = ");
        s.append(namemap(a,b,cplus1,dmins1,m));
        s.append("_con + ");
        s.append("CD[");
        s.append(std::to_string(xyz));
        s.append("] * ");
        s.append(namemap(a,b,c,dmins1,m));
        s.append("_con ; \n");
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
    s.append("(shell_pairs_gpu ab, shell_pairs_gpu cd,  double* I_ ");
    // for (auto l = la +lb; l >= la; l--)
    // for (auto r = lc +ld; r >= lc; r--)
    // {
    //     s.append(", double* I_");
    //     s.append(orbname(std::array<int, 3> {l,0,0}));
    //     s.append("s");
    //     s.append(orbname(std::array<int, 3> {r,0,0}));
    //     s.append("s");
    //     // s.append(", ");
    // }
    s.append(")\n{");

    return s;
}


void code_print(int la , int lb, int lc, int ld, std::map<std::array<int, 13>, std::string>& osmap)
{
    std::cout << HrrFuncName(la,lb,lc,ld) << std::endl;
    
    std::cout << "    auto idx = threadIdx.x + blockDim.x * blockIdx.x; // idx for shell_pairs ab \n\
    auto idy = threadIdx.y + blockDim.y * blockIdx.y; // idy for shell_pairs cd \n\n";
    
    std::cout << "    auto nab = ab.na * ab.nb; //  sizeof(ab.S)/sizeof(ab.S[0]);\n";
    std::cout << "    auto ncd = cd.na * cd.nb; //  sizeof(cd.S)/sizeof(cd.S[0]);\n\n";
    std::cout << "    auto mab = ab.ma * ab.mb; \n";
    std::cout << "    auto mcd = cd.ma * cd.mb; \n\n";
    std::cout << "    if (idx < nab  and idy < ncd)\n    {\n";

    std::cout << "    double AB[3] = \n    {\n        ab.AB[0][idx],\n        ab.AB[1][idx],\n        ab.AB[2][idx]\n    };\n";
    std::cout << "    double CD[3] = \n    {\n        cd.AB[0][idx],\n        cd.AB[1][idx],\n        cd.AB[2][idx]\n    };\n\n";
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
            if (bx == 0 and by == 0 and nb-bx-by == 0 and dx == 0 and dy == 0 and nd-dx-dy == 0)
            {
            std::array<int, 13> os_id = {ax, ay, na-ax-ay, bx, by, nb-bx-by, cx, cy, nc-cx-cy, dx, dy, nd-dx-dy, 0};
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
        double fm[" << la+lb+lc+ld+1 << "]; \n\
        vgamma<double>(" << la+lb+lc+ld << ", T, fm);\n\
        \n\
        for (auto i = 0; i <" << la+lb+lc+ld+1 << 
        " ; i++) \n\            fm[i] *= K;\n" << std::endl;
        merge_vrr_code( la,  lb,  lc,  ld);
        std::cout << "        }" << std::endl;

    std::cout << "         //****************//" << std::endl;


    // cotracted integrals equation

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
            if (osmap.count(os_id) != 0)
                std::cout << osmap.at(os_id);
            }
        }   

    }

    save_int(la, lb, lc, ld);
    std::cout << "    }\n}" << std::endl;
}


void merge_vrr_code(int la, int lb, int lc, int ld)
{
    std::map<std::array<int, 13>, std::string> osmap;
    for (auto l = la +lb; l >= la; l--)
    for (auto r = lc +ld; r >= lc; r--)
    {
        eri(l,0,r,0,osmap);
    }
    // printf("osmap size is %d\n", osmap.size());
    vrr_code_print(la+lb,0,lc+ld,0,osmap);

    for (auto l = la +lb; l >= la; l--)
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
            std::cout << "        " << namemap(a,b,c,d,0) << "_con += " << namemap(a,b,c,d,0) << ";" << std::endl;
        }   
    }


}

void hrr_host(int la, int lb, int lc, int ld)
{
    // print the host function for hrr
    std::string intename;
    intename.append(orbname(std::array<int, 3> {la,0,0}));
    intename.append(orbname(std::array<int, 3> {lb,0,0}));
    intename.append(orbname(std::array<int, 3> {lc,0,0}));
    intename.append(orbname(std::array<int, 3> {ld,0,0}));
    std::cout << "__host__ void " << intename << "(shell_pairs ab, shell_pairs cd)" << std::endl;
    std::cout << "{" << std::endl;
    std::cout << "    int lgab = ab.S.size(); // length of ab" << std::endl;
    std::cout << "    int lgcd = cd.S.size(); // length of cd\n" << std::endl;

    std::cout << "    shell_pairs_gpu ab_d (ab);" << std::endl;
    std::cout << "    shell_pairs_gpu cd_d (cd)\n;" << std::endl;

    std::cout << "    int block_size = 16;" << std::endl;
    std::cout << "    dim3 block_dim(block_size, block_size);" << std::endl;
    std::cout << "    dim3 grid_dim((((lgab-1)/block_size)+1), (((lgcd-1)/block_size)+1));\n" << std::endl;
    
    // for (auto l = la +lb; l >= la; l--)
    // for (auto r = lc +ld; r >= lc; r--)
    // {
    //     std::string s;
    //     std::string name;
    //     s.append("\n    // ");
    //     name.append(orbname(std::array<int, 3> {l,0,0}));
    //     name.append("s");
    //     name.append(orbname(std::array<int, 3> {r,0,0}));
    //     name.append("s");
    //     s.append(name);
    //     s.append(" integral\n");
    //     s.append("    auto llabcd_");
    //     s.append(name);
    //     s.append(" = ");
    //     s.append(std::to_string((l+1)*(l+2)/2*1*1*(r+1)*(r+2)/2));
    //     s.append(" ;\n");
    //     s.append("    auto lintegral_");
    //     s.append(name);
    //     s.append(" = llabcd_");
    //     s.append(name);
    //     s.append(" * lgab * lgcd;\n\n");
    //     s.append("    // allocate memory for gpu");
    //     s.append("    double *I_");
    //     s.append(name);
    //     s.append("_d\n");
    //     s.append("    cudaMalloc((void**) &I_");
    //     s.append(name);
    //     s.append("_d, lintegral_");
    //     s.append(name);
    //     s.append(" * sizeof(double));\n");
    //     s.append("    vrr_con_");
    //     s.append(name);
    //     s.append("<<<grid_dim,block_dim>>>(ab_d, cd_d, I_");
    //     s.append(name);
    //     s.append("_d);");
    //     std::cout << s << std::endl;
    // }
    std::cout << "    // HRR for " << intename << std::endl;
    std::cout << "    auto llabcd_" << intename << " = " << (la+1)*(la+2)*(lb+1)*(lb+2)*(lc+1)*(lc+2)*(ld+1)*(ld+2)/16 << ";"<< std::endl;
    std::cout << "    auto lintegral_" << intename << " = llabcd_" << intename << " * lgab * lgcd;\n" << std::endl;
    std::cout << "    // allocate memory for gpu" << std::endl;
    std::cout << "    double *I_" << intename << "_d;"<< std::endl;
    std::cout << "    cudaMalloc((void**) &I_" << intename << "_d, lintegral_" << intename << " * sizeof(double));\n" << std::endl;
    std::cout << "    // pinned memory" << std::endl;
    std::cout << "    double *I_" << intename << "_h;" << std::endl;
    std::cout << "    cudaHostAlloc((void**) &I_" << intename << "_h, lintegral_" << intename << " * sizeof(double), cudaHostAllocDefault);\n" << std::endl;
    std::cout << "    hrr_" << intename << "<<<grid_dim,block_dim>>>(ab_d, cd_d, I_" << intename << "_d);";
    // for (auto l = la +lb; l >= la; l--)
    // for (auto r = lc +ld; r >= lc; r--)
    // {
    //     std::string intermediate_name;
    //     intermediate_name.append(orbname(std::array<int, 3> {l,0,0}));
    //     intermediate_name.append("s");
    //     intermediate_name.append(orbname(std::array<int, 3> {r,0,0}));
    //     intermediate_name.append("s");
    //     std::cout << " I_" << intermediate_name << "_d,";
    // }
    std::cout << "\n\n    cudaDeviceSynchronize();\n" << std::endl;
    std::cout << "    // copy data back" << std::endl;
    std::cout << "    cudaMemcpy(I_" << intename << "_h, I_" << intename << "_d, lintegral_" << intename << " * sizeof(double),  cudaMemcpyDeviceToHost);" << std::endl;
    std::cout << "    // free gpu memory" << std::endl;
    // for (auto l = la +lb; l >= la; l--)
    // for (auto r = lc +ld; r >= lc; r--)
    // {
    //     std::string intermediate_name;
    //     intermediate_name.append(orbname(std::array<int, 3> {l,0,0}));
    //     intermediate_name.append("s");
    //     intermediate_name.append(orbname(std::array<int, 3> {r,0,0}));
    //     intermediate_name.append("s");
    //     std::cout << "    cudaFree(I_" << intermediate_name << "_d);" << std::endl;
    // }
    std::cout << "    cudaFree(I_" << intename << "_d);" << std::endl;
    std::cout << "    cudaFreeHost(I_" << intename << "_h);" << std::endl;
    std::cout << "}" << std::endl;
}

#if 0
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

    int la = 2; 
    int lb = 2; 
    int lc = 2; 
    int ld = 2;

    for(auto a = 0 ; a <= la; a++)
    for(auto b = 0 ; b <= a; b++)
    for(auto c = 0 ; c <= a; c++)
    for(auto d = 0 ; d <= c and d <= b; d++)
    {
        hrr_host(a, b, c, d);
        // std::map<std::array<int, 13>, std::string> osmap;
        // hgp_eri(a, b, c, d, osmap);
        // code_print(a, b, c, d, osmap);
        // printf("\n\n");
    }

    // for (auto l = la +lb; l >= la; l--)
    // for (auto r = lc +ld; r >= lc; r--)
    // {
    //     std::map<std::array<int, 13>, std::string> osmap;
    //     eri(l,0,r,0,osmap);
    //     vrr_code_print(l,0,r,0,osmap);
    //     // printf("osmap size is %d\n", osmap.size());
    // }
    // merge_vrr_code( la,  lb,  lc,  ld);
    //auto funcname = HrrFuncName(la,lb,lc,ld);
    //std::cout << funcname << std::endl;
    //  std::map<std::array<int, 13>, std::string> osmap;
    //  // eri(la,lb,lc,ld,osmap);
    //  // vrr_code_print(la , lb, lc, ld, osmap);
    //  hgp_eri(la, lb, lc, ld, osmap);
    //  code_print(la, lb, lc, ld, osmap);
    //  hrr_host(la,  lb, lc, ld);
    // save_int(la, lb, lc, ld);
    // std::cout << hello << std::endl;
    // printf("Hello\n");
}
#endif

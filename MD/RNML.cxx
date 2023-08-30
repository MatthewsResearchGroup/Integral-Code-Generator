/*
 * This code is desiged to generate the code for evluate RNLM intermidiates in MD algorithm,
 * which will be transfered to a matrix form to accrelate the following calculation.
 */

#include <string>
#include <iostream>
#include <array>
#include <map>

typedef std::string str;
// using tostr = std::to_string;

// return the name of RNLM, whose form will be R_1_1_1_0 if N = 1, M = 1, L = 1 and j = 0.
str RNLM_naming(int N, int L, int M, int j)
{
    str RNLMj;
    
    RNLMj = "R_" + std::to_string(N) + "_" + std::to_string(L) + "_" + std::to_string(M) + "_" + std::to_string(j);

    return RNLMj;
}

str RNLM_equ(int N, int L, int M, int j)
{
    str sen;

    // if N, L, and M are 0, return the gamma function
    if (N == 0 and L == 0 and M == 0)
    {
        sen = RNLM_naming(N, L, M, j) + " = Fm[" + std::to_string(j) + "];"; 
    }

    // if N = 0 and L = 0, R,0,0,M+1,j = c*R,0,0,M,j+1 + M R,0,0,M-1,j+1
    else if (N == 0 and L == 0)
    {
        sen = RNLM_naming(N, L, M, j) + " = " + "c*" + RNLM_naming(N, L, M-1, j+1);
        // the last term will be 0 if M < 2, which is not needed to be included
        if ( M >= 2 )
            sen += " + " + std::to_string(M) + " * " + RNLM_naming(N, L, M-2, j+1); 
        sen += ";" ;
    }

    // if N = 0, R,0,L+1,M,j = c*R,0,L,M,j+1 + M R,0,L-1,M,j+1
    else if (N == 0)
    {
        sen = RNLM_naming(N, L, M, j) + " = " + "b*" + RNLM_naming(N, L-1, M, j+1);
        // the last term will be 0 if M < 2, which is not needed to be included
        if ( L >= 2 )
            sen += " + " + std::to_string(L) + " * " + RNLM_naming(N, L-2, M, j+1); 
        sen += ";" ;
    }
    else
    {
        sen = RNLM_naming(N, L, M, j) + " = " + "a*" + RNLM_naming(N-1, L, M, j+1);
        // the last term will be 0 if M < 2, which is not needed to be included
        if ( N >= 2 )
            sen += " + " + std::to_string(N) + " * " + RNLM_naming(N-2, L, M, j+1); 
        sen += ";" ;
    }
    
    return sen;
}

// build the recursive relations for RNLM
// template <class T>
void RNLM_rr(int N, int L, int M, int j, std::map<std::array<int, 4>, std::string >& RNLMmap)
{
    std::array<int, 4> RNLM_id = {N, L, M, j};
    str RNLM_sen = RNLM_equ(N, L, M, j);
    if (RNLMmap.count(RNLM_id) > 0)
    {
        return ;
    }
    RNLMmap.insert(std::pair<std::array<int, 4>, std::string>(RNLM_id, RNLM_equ(N, L, M, j)));

    if (N == 0 and L ==0 and M == 0)
    {   
        // The key for RNLM map using to find the value in the map
        // auto RNLM_SEN = 
        return;
    }

    else if (N == 0 and L == 0)
    {
        RNLM_rr(N, L, M-1, j+1, RNLMmap);
        if (M >= 2)
            RNLM_rr(N, L, M-2, j+1, RNLMmap);
    }

    else if (N == 0)
    {
        RNLM_rr(N, L-1, M, j+1, RNLMmap);
        if (L >= 2)
            RNLM_rr(N, L-2, M, j+1, RNLMmap);
    } 

    else 
    {
        RNLM_rr(N-1, L, M, j+1, RNLMmap);
        if (N >= 2)
            RNLM_rr(N-2, L, M, j+1, RNLMmap);
    } 
}

void RNLM_matrix(int la, int lb, int lc, int ld,  std::map<std::array<int, 4>, std::string >& RNLMmap)
{
    // get ax, ay ,az
    for (auto ax = la; ax >= 0; ax--)
    for (auto ay = la - ax; ay >= 0; ay--)
    {
        auto az = la - ax - ay;

        // get all combinations of lb
        for (auto bx = lb; bx >= 0; bx--)
        for (auto by = lb - bx; by >= 0; by--)
        {
            auto bz = lb - bx - by;

            // get all combinations of lc
            for (auto cx = lc; cx >= 0; cx--)
            for (auto cy = lc - cx; cy >= 0; cy--)
            {
                auto cz = lc - cx - cy;

                // get all combinations of ld
                for (auto dx = ld; dx >= 0; dx--)
                for (auto dy = ld - dx; dy >= 0; dy--)
                {
                    auto dz = ld - dx - dy;
                    
                    auto abx = ax + bx;
                    auto aby = ay + by;
                    auto abz = az + bz;
                    auto cdx = cx + dx;
                    auto cdy = cy + dy;
                    auto cdz = cz + dz;
                    // get all combinations of NLM and NPLPMP
                    for (auto N = 0; N <= abx; N++)
                    for (auto L = 0; L <= aby; L++)
                    for (auto M = 0; M <= abz; M++)
                    for (auto NP = 0; NP <= cdx; NP++)
                    for (auto LP = 0; LP <= cdy; LP++)
                    for (auto MP = 0; MP <= cdz; MP++)
                    {
                        auto NNP = N + NP;
                        auto LLP = L + LP;
                        auto MMP = M + MP;
                        RNLM_rr(NNP, LLP, MMP, 0, RNLMmap);
                    }
                }
            }
        }

    }
}

void RNML_print(int la, int lb, int lc, int ld, std::map<std::array<int, 4>, std::string >& RNLMmap)
{
    int count = 0;
    auto llabcd = la + lb + lc + ld;
    for (auto NNP = 0; NNP <= llabcd; NNP++)
    for (auto LLP = 0; LLP <= llabcd; LLP++)
    for (auto MMP = 0; MMP <= llabcd; MMP++)
    for (auto j = 0; j <= llabcd; j++)
    {
        std::array<int, 4> RNML_id = {NNP, LLP, MMP, j};
        if (RNLMmap.count(RNML_id) > 0)
        {
            std::cout << RNLMmap.at(RNML_id) << std::endl;
            count += 1;
        }
    }
    printf("The size of RR is %d\n", count);
}

int main()
{
    int la = 4;
    int lb = 4;
    int lc = 4;
    int ld = 4;
    int N = 6;
    int M = 6;
    int L = 6;
    int j = 0;

    std::map<std::array<int, 4>, std::string> RNLMmap;
    // std::cout << RNLM_naming(N, L, M, j) << std::endl;
    // std::cout << RNLM_equ(2, 1, 2, 2) << std::endl;
    // RNLM_rr(N, M, L, j, RNLMmap);   
    RNLM_matrix(la, lb, lc, ld, RNLMmap);
    RNML_print(la, lb, lc, ld, RNLMmap);
    std::cout << "The size of map: " << RNLMmap.size() << std::endl;
    return 0;
}

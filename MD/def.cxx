/* This code is designed to generate the recursive code for def coefficients
 * in MD algorithms ()
 */

#include <iostream>
#include <string>

typedef std::string str;


/********************************************************
 *
 * The Revurssion Relation for def 
 *
 *  d_{N}^{n+1, n'} = 1/(2\zeta_{ab}) * d_{N-1}^{n, n'} +
 *                   (P_x - A_x)      * d_{N  }^{n, n'} +
 *                   (N + 1)          * d_{N+1}^{n, n'};
 *
 *
 *******************************************************/


str def_naming(int N, int n, int np)
{
    str s;
    if (n < 0 || np < 0 || N < 0)
        return "0.0 ";
    return "def_" + std::to_string(N) + "_" + std::to_string(n) + "_"  + std::to_string(np);
}

std::string def_equ(int n, int np, int N)
{
    str s;
    if (N == 0 && n == 0 && np == 0)
    {
        s = def_naming(n , np, N)  +  " = 1.0;";
    }


    // n > np, we could decrement np first. 
    else if (n < np || np == 0)
    {
        // if N + 2 == n + np; the last term will be 0
        if (N + 1 == n + np)
        {
            s = def_naming(n, np, N) + " = 1/(2 * zeta_ab) * " + def_naming(n - 1, np, N - 1) + " + (Px - Ax) * " + def_naming(n - 1, np, N) + ";"; 
        }
        else if ( N == n + np)
        {
            s = def_naming(n, np, N) + " = 1/(2 * zeta_ab) * " + def_naming(n - 1, np, N - 1) + ";" ; 
        }
        else if (N > n + np)
        {
            s = def_naming(n, np, N) + " = 0.0;";
        }
        else
        {
            s = def_naming(n, np, N) + " = 1/(2 * zeta_ab) * " + def_naming(n - 1, np, N - 1) + " + (Px - Ax) * " + def_naming(n - 1, np, N) +  " + " + std::to_string(N - 1) + " * " + def_naming(n - 1, np, N + 1) + ";"; 
        }
    }

    else if (n > np || n == 0)
    {
        // if N + 2 == n + np; the last term will be 0
        if (N + 1 == n + np)
        {
            s = def_naming(n, np, N) + " = 1/(2 * zeta_ab) * " + def_naming(n, np - 1, N - 1) + " + (Px - Ax) * " + def_naming(n, np - 1, N) + ";"; 
        }
        else if ( N == n + np)
        {
            s = def_naming(n, np, N) + " = 1/(2 * zeta_ab) * " + def_naming(n, np - 1, N - 1) + ";" ; 
        }
        else if (N > n + np)
        {
            s = def_naming(n, np, N) + " = 0.0;";
        }
        else
        {
            s = def_naming(n, np, N) + " = 1/(2 * zeta_ab) * " + def_naming(n, np - 1, N - 1) + " + (Px - Ax) * " + def_naming(n, np - 1, N) +  " + " + std::to_string(N - 1) + " * " + def_naming(n, np - 1, N + 1) + ";"; 
        }

    }
    return s;
}


double def_rr(int N, int n, int np)
{
    if (N == 0 && n == 0 && np == 0)
    {
        std::cout << "Here 1." << std::endl;
        return 1.0;
    }

    else if (N > n + np || n < 0 || np < 0 || N < 0)
    {
        std::cout << "Here 2." << std::endl;
        return 0.0;
    }

    std::cout << N  << "," << n << ", " << np << std::endl;
    return def_rr(N - 1, n - 1, np) + def_rr(N, n - 1, np) + def_rr(N + 1, n - 1, np);      
}



int main()
{
    int N = 1;
    int n = 2;
    int np = 0;

    std::cout << def_rr(N, n, np) << std::endl;

    str s = def_equ(N, n, np);

    std::cout << s << std::endl;

    return 0;

}

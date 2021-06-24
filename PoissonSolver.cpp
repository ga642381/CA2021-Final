#include "classes.h"
#include <stdio.h>
#include <string>

using namespace std;

int main(int argc, char const *argv[])
{
    string method;
    int M, alg_num, func, steps = 0;

    /* === Input Parameters === */
    // algorithms
    printf("\n");
    printf("1. Choose algorithms:\n\
    [1] SOR\n\
    [2] Two Grid\n\
    [3] V-Cycle\n\
    [4] W-Cycle\n");
    printf("    algoritm : ");
    scanf("%d", &alg_num);

    // Matrix Size
    // h : interval (delta)
    // N : grid number
    printf("\n2. Choose M. h=1/M (interval) ; N=M-1 (grid number in x and y direction):\n");
    printf("    M : ");
    scanf("%d", &M);

    // Equation
    printf("\n3. Choose equation:\n\
    [1] -nabla u(x,y) = -4, and u(x,y) = x^2 + y^2.\n\
    [2] -nabla u(x,y) = 0, and u(x,y) = 1.\n");
    printf("    equation : ");
    scanf("%d", &func);
    M = M - 1;

    /* === Initialize === */
    PoissonMatrix A(M);
    LowerMatrix L(M);
    UpperMatrix U(M);
    Boundary B(M, func);
    Startvector X(M, 0.0);
    Algorithms Algs(M);

    double timer, start = 0.0, end = 0.0;
    start = clock();

    /* === Algorithms === */
    if (alg_num == 1)
    {
        method = "SOR";
        steps = Algs.SORMethod(X.x, B.b, B.solved);
    }

    if (alg_num == 2)
    {
        method = "Two-Grid";
        steps = Algs.MultiGridMethod(A,X.x, B.b, B.solved, method);
    }
    if (alg_num == 3)
    {
        method = "V-Cycle";
        steps = Algs.MultiGridMethod(A,X.x, B.b, B.solved, method); 
    }
    if (alg_num == 4)
    {
        method = "W-Cycle";
        steps = Algs.MultiGridMethod(A,X.x, B.b, B.solved, method);
    }

    /* === Information ===*/
    end = clock();
    timer = (end - start) / CLOCKS_PER_SEC;
    cout << "==========\n";
    cout << "Method : " << method << endl;
    cout << "Steps : " << steps << endl;
    cout << "Time : " << timer << " (s)" << endl;

    /* === For Plotting ===*/
    X.WriteToFile();

    return EXIT_SUCCESS;
}

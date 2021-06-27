#include "classes.h"

using namespace std;

int main(int argc, char const *argv[])
{
    int num_thread;

    printf("\n");
    printf("0. Set number of thread for OpenMP\n");
    printf("    Number of Thread : ");
    scanf("%d", &num_thread);

    printf("\n");
#pragma omp parallel num_threads(num_thread)
    {
        int id = omp_get_thread_num();
        printf("Hi I am thread number : %d\n", id);
    }

    //==============================================================================================================

    string method;
    int M, alg_num, func, steps, fixed_steps = 0;

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
    [2] -nabla u(x,y) = 0, and u(x,y) = 1.\n\
    [3] -nabla u(x,y) = 2[(6x^2-6x+1)y^2(y-1)^2+(6y^2-6y+1)x^2(x-1)^2], and u(x,y) = x^2(x-1)^2y^2(y-1)^2.\n\
    [4] -nabla u(x,y) = 2*x*(x-1) + 2*y*(y-1), and u(x,y) = x*(x-1)*y*(y-1).\n\
    [5] -nabla u(x,y) = 0 and u(x,y) = exp(-2pi * x) * sin(2pi * y)\n");
    printf("    equation : ");
    scanf("%d", &func);

    // If fixed steps
    printf("\n4. Update fixed steps (if run until converge to analytical sol, please enter -1)\n");
    printf("    fixed steps : ");
    scanf("%d", &fixed_steps);
    printf("\n");

    M = M - 1;

    /* === Initialize === */
    PoissonMatrix A(M);
    Boundary B(M, func);
    Startvector X(M, 0.0);
    Algorithms Algs(M, num_thread);

    double start_time, end_time;
    start_time = omp_get_wtime();

    /* === Algorithms === */
    if (alg_num == 1)
    {
        method = "SOR";
        steps = Algs.SORMethod(X.x, B.b, B.solved, fixed_steps);
    }

    if (alg_num == 2)
    {
        method = "Two-Grid";
        steps = Algs.MultiGridMethod(A, X.x, B.b, B.solved, method, fixed_steps);
    }
    if (alg_num == 3)
    {
        method = "V-Cycle";
        steps = Algs.MultiGridMethod(A, X.x, B.b, B.solved, method, fixed_steps);
    }
    if (alg_num == 4)
    {
        method = "W-Cycle";
        steps = Algs.MultiGridMethod(A, X.x, B.b, B.solved, method, fixed_steps);
    }

    /* === Information ===*/
    end_time = omp_get_wtime();

    cout << "==========\n";
    cout << "Thread Number : " << num_thread << endl;
    cout << "Method : " << method << endl;
    cout << "Grid Size : " << M << endl;
    cout << "Steps : " << steps << endl;
    cout << "Update Steps : " << Algs.update_step << endl;
    cout << "Time : " << (end_time - start_time) << " (s)" << endl;

    /* === For Plotting ===*/
    //X.WriteToFile();

    return EXIT_SUCCESS;
}

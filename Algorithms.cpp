#include "classes.h"

using namespace std;

double cal_error_abs(vector<double> &x, const vector<double> &y)
{
    vector<double> r = x - y;
    int N = (int)(r.size());
    double error = 0;
    for (int i = 0; i < N; i++)
    {
        error += abs(r[i]);
    }

    return error / (double)(r.size());
}

double cal_error(vector<double> &x, const vector<double> &y)
{
    vector<double> r = x - y;
    return pow((r | r), 2) / double(r.size());
}

Algorithms::Algorithms(int n, int num_threads)
{
    this->n = n;
    this->num_threads = num_threads;
    this->dim = n * n;
    this->h = 1.0 / (double)(n + 1);
    this->Vcounter = 0;
    this->update_step = 0;
}

Algorithms::~Algorithms() {}

int Algorithms::SORMethod(vector<double> &x, const vector<double> &b, const vector<double> &solved, const int fixed_step)
{
    //double Pi = 3.141592654;
    //double omega = 2 / (1 + sqrt(1 - pow(cos(Pi * h), 2)));
    double omega = 1.5;
    vector<double> tmp = x;
    vector<double> sum_rec = x;
    int n = sqrt(x.size());
    int steps = 0;
    double h = 1.0 / (double)(n + 1);
    // residual
    vector<double> r(x.size());
    r = x - solved;

    //  === iterate until converge === //
    if (fixed_step < 0)
    {

        double TOL = (r | r) * pow(10, -3);
        while (TOL <= (r | r))
        {
#pragma omp parallel num_threads(this->num_threads)
            {
#pragma omp for
                for (int i = 0; i < n * n; i++) // even cell
                {
                    double sum = 0.0;
                    if ((i / n + i % n) % 2 == 0)
                    {
                        if (i >= n && i < dim - n)
                            sum += x[i - n] + x[i + n];
                        if (i < n)
                            sum += x[i + n];
                        if (i >= dim - n)
                            sum += x[i - n];
                        if (i % n != 0)
                            sum += x[i - 1];
                        if (i % n != n - 1)
                            sum += x[i + 1];

                        tmp[i] = x[i] * (1 - omega);
                        x[i] = tmp[i] + omega * 1.0 / (4.0 * pow(1.0 / h, 2)) * (b[i] + pow(1.0 / h, 2) * sum);
                        r[i] = x[i] - solved[i];
                    }
                } // end even cell
#pragma omp for
                for (int i = 0; i < n * n; i++) // odd cell
                {
                    double sum = 0.0;
                    if ((i / n + i % n) % 2 == 1)
                    {

                        if (i >= n && i < dim - n)
                            sum += x[i - n] + x[i + n];
                        if (i < n)
                            sum += x[i + n];
                        if (i >= dim - n)
                            sum += x[i - n];
                        if (i % n != 0)
                            sum += x[i - 1];
                        if (i % n != n - 1)
                            sum += x[i + 1];

                        tmp[i] = x[i] * (1 - omega);
                        x[i] = tmp[i] + omega * 1.0 / (4.0 * pow(1.0 / h, 2)) * (b[i] + pow(1.0 / h, 2) * sum);
                        r[i] = x[i] - solved[i];
                    }
                } // end odd cell

            } // end parallel
            this->update_step++;

        } // end while
    }     // end if fixed_step <0

    // =================== Given Fixed Steps =========================//
    else
    {
        while (this->update_step < fixed_step)
        {
#pragma omp parallel num_threads(this->num_threads)
            {

#pragma omp for
                for (int i = 0; i < n * n; i++) // even cell
                {
                    double sum = 0.0;
                    if ((i / n + i % n) % 2 == 0)
                    {
                        if (i >= n && i < dim - n)
                            sum += x[i - n] + x[i + n];
                        if (i < n)
                            sum += x[i + n];
                        if (i >= dim - n)
                            sum += x[i - n];
                        if (i % n != 0)
                            sum += x[i - 1];
                        if (i % n != n - 1)
                            sum += x[i + 1];

                        tmp[i] = x[i] * (1 - omega);
                        x[i] = tmp[i] + omega * 1.0 / (4.0 * pow(1.0 / h, 2)) * (b[i] + pow(1.0 / h, 2) * sum);
                        r[i] = x[i] - solved[i];
                    }
                } // end even cell

#pragma omp for
                for (int i = 0; i < n * n; i++) // odd cell
                {
                    double sum = 0.0;
                    if ((i / n + i % n) % 2 == 1)
                    {

                        if (i >= n && i < dim - n)
                            sum += x[i - n] + x[i + n];
                        if (i < n)
                            sum += x[i + n];
                        if (i >= dim - n)
                            sum += x[i - n];
                        if (i % n != 0)
                            sum += x[i - 1];
                        if (i % n != n - 1)
                            sum += x[i + 1];

                        tmp[i] = x[i] * (1 - omega);
                        x[i] = tmp[i] + omega * 1.0 / (4.0 * pow(1.0 / h, 2)) * (b[i] + pow(1.0 / h, 2) * sum);
                        r[i] = x[i] - solved[i];
                    }
                } // end odd cell
            }     // end parallel
            steps++;
            this->update_step++;
        } // end while
    }
    //======== end SOR, print and return ========//
    cout << "==========\n";
    double error = cal_error(x, solved);
    cout << "Error (L2): " << error << endl;
    return steps;
    // end parallel decalre
}

int Algorithms::MultiGridMethod(Matrix &A, vector<double> &x, const vector<double> &b, const vector<double> &solved, string method, const int fixed_step)
{
    // MultiGrid Method
    // To be implemented
    // input should be Matrix
    int steps = 0, dim = x.size();
    vector<double> r(dim);
    r = x - solved;
    int numberOfGrids = 0;
    int VW = 0;

    if (method == "Two-Grid")
    {
        numberOfGrids = 1;
        VW = 0;
    }

    else if (method == "V-Cycle")
    {
        numberOfGrids = 2;
        VW = 0;
    }

    else if (method == "W-Cycle")
    {
        numberOfGrids = 2;
        VW = 1;
    }

    // iterate until converge
    if (fixed_step < 0)
    {
        // step 0
        //double error = cal_error_abs(x, solved);
        //this->error_array.push_back(error);

        double TOL = pow(10, -3) * (r | r);
        while (TOL <= (r | r))
        {
            x = Cycle(A, x, b, numberOfGrids, VW);
            r = x - solved;

            // double error = cal_error_abs(x, solved);
            // this->error_array.push_back(error);
            steps++;
        }
    }

    // run fixed step (might be more, we run until the whole cycle is finished)
    else
    {
        while (this->update_step < fixed_step)
        {
            x = Cycle(A, x, b, numberOfGrids, VW);
            r = x - solved;
            steps++;
        }
    }
    cout << "==========\n";
    double error = cal_error(x, solved);
    cout << "error : " << error << endl;
    return steps;
}

vector<double> Algorithms::Cycle(Matrix &A, vector<double> &x, const vector<double> &b, int level, int w_cycle)
{
    // level : number of grids
    // w_cycle : VW

    // Multigrid : (level, w_cycle)
    // Two-Grid  : (1, 0)
    // V-Cycle   : (2, 0)
    // W-Cycle   : (2, 1)

    // using "recursion" to achieve V-Cycle and W-Cycle
    // each Cycle contains :
    // 1. Restriction
    // 2. Cycle
    // 3. Cycle (if W-Cycle)
    // 4. Prolongation

    // when at the bottom level, do SORRelaxation

    int dim = x.size(), n = sqrt(dim), N2h = (n + 1) / 2 - 1, dim2h = pow(N2h, 2);
    vector<double> r(dim, 0), E(dim, 0), sol(dim, 0), r2h(dim2h, 0), E2h(dim2h, 0), sol2h(dim2h, 0);

    int smoothing_step = 3;
    if (this->Vcounter == level)
    {
        SORRelaxation(A, x, b, -1);
        return x;
    }
    else
    {

        this->Vcounter++;
        // pre-smooth
        JacobiRelaxation(A, x, b, smoothing_step);

        r = b - A * x; // residual
        Restriction(r, r2h, n);

        E2h = Cycle(A, E2h, r2h, level, w_cycle);

        // for W-Cycle
        if (w_cycle == 1)
            E2h = Cycle(A, E2h, r2h, level, w_cycle);

        Prolongation(E2h, E, n);
        x += E; // correction

        // post smooth
        JacobiRelaxation(A, x, b, smoothing_step);

        this->Vcounter--;
        return x;
    }
}

void Algorithms::JacobiMethod(Matrix &A, vector<double> &x, const vector<double> &b, const vector<double> &solved)
{
    vector<double> tmp;
    double sum;
    int n = sqrt(x.size()), dim = n * n, steps = 0;
    double h = 1.0 / (double)(n + 1);
    vector<double> r(x.size());
    r = x - solved;
    double TOL = (r | r) * pow(10, -3);
    while (TOL <= (r | r))
    {
        cout << (r | r) / TOL << endl;
        tmp = x;
        for (int i = 0; i < n * n; i++)
        {
            sum = 0.0;
            if (i >= n && i < dim - n)
                sum += tmp[i - n] + tmp[i + n];
            if (i < n)
                sum += tmp[i + n];
            if (i >= dim - n)
                sum += tmp[i - n];
            if (i % n != 0)
                sum += tmp[i - 1];
            if (i % n != n - 1)
                sum += tmp[i + 1];
            x[i] = 1.0 / (4.0 * pow(1.0 / h, 2)) * (b[i] + pow(1.0 / h, 2) * sum);
            r[i] = x[i] - solved[i];
        }
        steps++;
    }
    //return steps;
}

void Algorithms::JacobiRelaxation(Matrix &A, vector<double> &x, const vector<double> &b, int steps)
{ //a,x,b,3
    vector<double> tmp;
    tmp = x;
    double sum, omega = 4.0 / 5.0;
    int n = sqrt(x.size()), dim = n * n;
    for (int k = 0; k < steps; k++)
    {
        for (int i = 0; i < n * n; i++)
        {
            sum = 0.0;
            if (i >= n && i < dim - n)
                sum += tmp[i - n] + tmp[i + n];
            if (i < n)
                sum += tmp[i + n];
            if (i >= dim - n)
                sum += tmp[i - n];
            if (i % n != 0)
                sum += tmp[i - 1];
            if (i % n != n - 1)
                sum += tmp[i + 1];
            x[i] = omega * 1.0 / (4.0 * pow(1.0 / h, 2)) * (b[i] + pow(1.0 / h, 2) * sum) + tmp[i] * (1 - omega);
            tmp[i] = x[i];
        }
        this->update_step++;
    }
}

void Algorithms::SORRelaxation(Matrix &A, vector<double> &x, const vector<double> &b, int steps)
{
    vector<double> x_prev;
    double Pi = 3.141592654, omega = 2 / (1 + sqrt(1 - pow(cos(Pi * h), 2)));
    int n = sqrt(x.size()), dim = n * n;

    if (steps < 0)
    {
        double residual = 999;
        while (residual > 1e-8)
        {
            x_prev = x;
#pragma omp parallel num_threads(this->num_threads)
            {
#pragma omp for
                for (int i = 0; i < n * n; i++) // even cell
                {
                    double sum = 0.0;
                    if ((i / n + i % n) % 2 == 0)
                    {
                        if (i >= n && i < dim - n)
                            sum += x[i - n] + x[i + n];
                        if (i < n)
                            sum += x[i + n];
                        if (i >= dim - n)
                            sum += x[i - n];
                        if (i % n != 0)
                            sum += x[i - 1];
                        if (i % n != n - 1)
                            sum += x[i + 1];
                        x[i] = omega * 1.0 / (4.0 * pow(1.0 / h, 2)) * (b[i] + pow(1.0 / h, 2) * sum) + x[i] * (1 - omega);
                    }
                } // end even cell

#pragma omp for
                for (int i = 0; i < n * n; i++) // odd cell
                {
                    double sum = 0.0;
                    if ((i / n + i % n) % 2 == 1)
                    {
                        if (i >= n && i < dim - n)
                            sum += x[i - n] + x[i + n];
                        if (i < n)
                            sum += x[i + n];
                        if (i >= dim - n)
                            sum += x[i - n];
                        if (i % n != 0)
                            sum += x[i - 1];
                        if (i % n != n - 1)
                            sum += x[i + 1];
                        x[i] = omega * 1.0 / (4.0 * pow(1.0 / h, 2)) * (b[i] + pow(1.0 / h, 2) * sum) + x[i] * (1 - omega);
                    }
                } // end odd cell
            }     // end parallel

            //====== Cal residual ===========//
            vector<double> r = x - x_prev;
            double r_tmp = 0;
            for (int i = 0; i < (int)(r.size()); i++)
            {
                r_tmp += abs(r[i]);
            }
            residual = r_tmp / (double)(r.size());
            this->update_step++;
        } // end while
    }     // end if step < 0

    // Relax given steps
    else
    {
        for (int k = 0; k < steps; k++) // step
        {
#pragma omp parallel num_threads(this->num_threads)
            {
#pragma omp for
                for (int i = 0; i < n * n; i++) // even cell
                {
                    double sum = 0.0;
                    if ((i / n + i % n) % 2 == 0)
                    {
                        if (i >= n && i < dim - n)
                            sum += x[i - n] + x[i + n];
                        if (i < n)
                            sum += x[i + n];
                        if (i >= dim - n)
                            sum += x[i - n];
                        if (i % n != 0)
                            sum += x[i - 1];
                        if (i % n != n - 1)
                            sum += x[i + 1];
                        x[i] = omega * 1.0 / (4.0 * pow(1.0 / h, 2)) * (b[i] + pow(1.0 / h, 2) * sum) + x[i] * (1 - omega);
                    }
                } // end even cell

#pragma omp for
                for (int i = 0; i < n * n; i++) // odd cell
                {
                    double sum = 0.0;
                    if ((i / n + i % n) % 2 == 1)
                    {
                        if (i >= n && i < dim - n)
                            sum += x[i - n] + x[i + n];
                        if (i < n)
                            sum += x[i + n];
                        if (i >= dim - n)
                            sum += x[i - n];
                        if (i % n != 0)
                            sum += x[i - 1];
                        if (i % n != n - 1)
                            sum += x[i + 1];
                        x[i] = omega * 1.0 / (4.0 * pow(1.0 / h, 2)) * (b[i] + pow(1.0 / h, 2) * sum) + x[i] * (1 - omega);
                    }
                } // end odd cell
                this->update_step++;
            }
        } // end for step
    }     // end step >0
}

void Algorithms::Restriction(const vector<double> &r, vector<double> &r2h, int n)
{
    for (int i = 1, l = 0, k = 0; i <= n; i++)
    {
        for (int j = 1; j <= n; j++, k++)
        {
            if (i % 2 == 0 && j % 2 == 0)
            {
                r2h[l] = 1.0 / 16.0 * (4.0 * r[k] + 2.0 * (r[k - 1] + r[k + 1] + r[k - n] + r[k + n]) + r[k + n - 1] + r[k + n + 1] + r[k - n - 1] + r[k - n + 1]);
                l++;
            }
        }
    }
}

void Algorithms::Prolongation(const vector<double> &E2h, vector<double> &E, int n)
{
    for (int i = 1, k = 0, l = 0; i <= n; i++)
    {
        for (int j = 1; j <= n; j++, k++)
        {
            if (i % 2 == 0 && j % 2 == 0)
            {
                E[k] = E2h[l];
                l++;
            }
        }
    }
    for (int i = 1, k = 0; i <= n; i++)
    {
        for (int j = 1; j <= n; j++, k++)
        {
            if (i % 2 == 0 && j % 2 != 0)
            {
                if (j != 1 && j != n)
                    E[k] = 1.0 / 2.0 * (E[k - 1] + E[k + 1]);
                if (j == 1)
                    E[k] = 1.0 / 1.0 * (E[k + 1]);
                if (j == n)
                    E[k] = 1.0 / 1.0 * (E[k - 1]);
            }
            if (i % 2 != 0 && j % 2 == 0)
            {
                if (i != 1 && i != n)
                    E[k] = 1.0 / 2.0 * (E[k - n] + E[k + n]);
                if (i == 1)
                    E[k] = 1.0 / 1.0 * (E[k + n]);
                if (i == n)
                    E[k] = 1.0 / 1.0 * (E[k - n]);
            }
            if (i % 2 != 0 && j % 2 != 0)
            {
                if (i != 1 && j != 1 && i != n && j != n)
                    E[k] = 1.0 / 4.0 * (E[k + n - 1] + E[k + n + 1] + E[k - n + 1] + E[k - n - 1]);
                if (i == 1 && j == 1)
                    E[k] = 1.0 / 1.0 * (E[k + n + 1]);
                if (i == n && j == n)
                    E[k] = 1.0 / 1.0 * (E[k - n - 1]);
                if (i == 1 && j == n)
                    E[k] = 1.0 / 1.0 * (E[k + n - 1]);
                if (i == n && j == 1)
                    E[k] = 1.0 / 1.0 * (E[k - n + 1]);
                if (i == 1 && j != 1 && j != n)
                    E[k] = 1.0 / 2.0 * (E[k + n + 1] + E[k + n - 1]);
                if (i == n && j != 1 && j != n)
                    E[k] = 1.0 / 2.0 * (E[k - n - 1] + E[k - n + 1]);
                if (j == 1 && i != 1 && i != n)
                    E[k] = 1.0 / 2.0 * (E[k + n + 1] + E[k - n + 1]);
                if (j == n && i != 1 && i != n)
                    E[k] = 1.0 / 2.0 * (E[k + n - 1] + E[k - n - 1]);
            }
        }
    }
}

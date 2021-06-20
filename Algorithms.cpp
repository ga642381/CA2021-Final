#include "classes.h"

Algorithms::Algorithms(int n)
{
    this->n = n;
    this->dim = n * n;
    this->h = 1.0 / (double)(n + 1);
    this->Vcounter = 0;
}

Algorithms::~Algorithms() {}

int Algorithms::MultiGridMethod(vector<double> &x, const vector<double> &b, const vector<double> &solved, string method)
{
    // MultiGrid Method
    // To be implemented
    // input should be Matrix
    if (method == "Two-Grid")
        return 0;
    else if (method == "V-Cycle")
        return 0;
    else if (method == "W-Cycle")
        return 0;
    else
        return 0;
}

int Algorithms::SORMethod(vector<double> &x, const vector<double> &b, const vector<double> &solved)
{

    double sum;
    double Pi = 3.141592654;
    double omega = 2 / (1 + sqrt(1 - pow(cos(Pi * h), 2)));

    int n = sqrt(x.size());
    int steps = 0;
    double h = 1.0 / (double)(n + 1);

    // residual
    vector<double> r(x.size());
    r = x - solved;
    double TOL = (r | r) * pow(10, -3);

    // update loop
    while (TOL <= (r | r))
    {
        // iterate each cell
        for (int i = 0; i < n * n; i++)
        {
            sum = 0.0;
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
            r[i] = x[i] - solved[i];
        }
        steps++;
    }
    return steps;
}